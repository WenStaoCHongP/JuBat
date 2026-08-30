mutable struct CZMResult
    displacement::Vector{Float64}
    damage::Vector{Float64}
    traction_n::Vector{Float64}
    traction_t::Vector{Float64}
    separation_n::Vector{Float64}
    separation_t::Vector{Float64}
    converged::Bool
    iterations::Int64
    residual_norm::Float64
    
    CZMResult(ndof::Int, n_coh::Int) = new(
        zeros(ndof), zeros(n_coh), zeros(n_coh), zeros(n_coh),
        zeros(n_coh), zeros(n_coh), false, 0, Inf)
end

function clone_damage_states(damage_states::AbstractVector{<:AbstractDamageState})
    return [begin
        new_state = DamageState()
        new_state.D = s.D
        new_state.D_visc = s.D_visc
        new_state.δ_max_n = s.δ_max_n
        new_state.δ_max_t = s.δ_max_t
        new_state.δ_max_eff = s.δ_max_eff
        new_state.fractured = s.fractured
        new_state.accumulated_damage = s.accumulated_damage
        new_state
    end for s in damage_states]
end

"""
    update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta=1.0)

按 cohesive 单元的 interface_type 分组，分别调用 `update_damage`
（界面参数宿主：:PE_PCC→param.PCC、:NE_NCC→param.NCC）。
当所有单元属于同一界面时，退化为单次 `update_damage` 调用。
"""
function update_damage_per_interface(czm_mesh::CohesiveMesh, damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, param::Params, czm_model::String; visc_beta::Float64=1.0)
    n_coh = czm_mesh.n_cohesive
    @assert length(damage_states) == n_coh "damage_states length mismatch"
    @assert length(separations) == n_coh "separations length mismatch"

    # 按 interface_type 分批（保持原始顺序）
    new_states = Vector{DamageState}(undef, n_coh)
    for iface in (:PE_PCC, :NE_NCC)
        idx = findall(i -> czm_mesh.cohesive_elements[i].interface_type == iface, 1:n_coh)
        isempty(idx) && continue
        ds_sub = damage_states[idx]
        sep_sub = separations[idx]
        updated = update_damage(ds_sub, sep_sub, collector_params(param, iface), czm_model; visc_beta=visc_beta)
        for (k, i) in enumerate(idx)
            new_states[i] = updated[k]
        end
    end
    return new_states
end

"""
    extract_bc_dofs(czm_mesh, param; fix_inner=true)

从 czm_mesh 提取 Dirichlet BC 的自由度列表和对应值（每次求解入口现算，
不缓存——identify_bc_nodes_czm 为 O(nnode)，成本可忽略）。
"""
function extract_bc_dofs(czm_mesh::CohesiveMesh, param; fix_inner::Bool=true)
    bc_nodes, _, _ = identify_bc_nodes_czm(czm_mesh, param; fix_inner=fix_inner)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        end
    end
    return bc_dofs, bc_vals
end

"""
    backtrack_line_search!(u, Δu, czm_mesh, param, damage_states, F_ext, F_thermo_chem, R_norm_current, bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws; max_halvings=8)

回溯线搜索（零化式 BC 残差）。仅用于 solve_czm_basic_step。
返回 (u_new, R_new_norm, accepted, α_used)。
accepted 时返回的 u_new 已含 BC 赋值，外部无需再执行 u = u + α*Δu。
未 accepted 时返回原始 u（未修改），外部应 break。
"""
function backtrack_line_search!(u::Vector{Float64}, Δu::Vector{Float64},czm_mesh::CohesiveMesh, param::Params, damage_states,F_ext::Vector{Float64}, F_thermo_chem::Vector{Float64},R_norm_current::Float64,bc_dofs::Vector{Int64}, bc_vals::Vector{Float64},K_bulk_cached, geom_cache, ws;max_halvings::Int=8, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    α = 1.0
    for _ in 1:max_halvings
        u_trial = u + α * Δu
        apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

        _, f_int_trial, _, _ = assemble_coupled_system(czm_mesh, u_trial, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain,
        plasticity=plasticity, mech_state=mech_state, prestress=prestress)

        R_trial = F_ext + F_thermo_chem - f_int_trial
        for (dof, val) in zip(bc_dofs, bc_vals)
            R_trial[dof] = val - u_trial[dof]
        end

        R_trial_norm = norm(R_trial)
        if !isnan(R_trial_norm) && R_trial_norm < R_norm_current
            return u_trial, R_trial_norm, true, α
        end

        α *= 0.5
    end
    return u, R_norm_current, false, 0.0
end

function apply_czm_dirichlet!(u::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64}, bc_vals::AbstractVector{Float64})
    for (dof, val) in zip(bc_dofs, bc_vals)
        u[dof] = val
    end
    return u
end

function zero_czm_bc_entries!(v::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64})
    for dof in bc_dofs
        v[dof] = 0.0
    end
    return v
end

function fill_czm_result!(result::CZMResult, u::Vector{Float64}, damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, tractions::Vector{Tuple{Float64, Float64}})
    result.displacement = u
    for i in eachindex(damage_states)
        result.damage[i] = damage_states[i].D
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end
    return result
end

function build_arc_length_augmented_matrix(K_bc::SparseMatrixCSC{Float64, Int64}, load_vector::Vector{Float64}, delta_u::Vector{Float64}, delta_lambda::Float64, arc_length_alpha::Float64)
    ndof = length(load_vector)
    A = spzeros(Float64, ndof + 1, ndof + 1)
    A[1:ndof, 1:ndof] = K_bc
    for i in 1:ndof
        A[i, ndof + 1] = -load_vector[i]
        A[ndof + 1, i] = 2.0 * delta_u[i]
    end
    A[ndof + 1, ndof + 1] = 2.0 * arc_length_alpha^2 * delta_lambda
    return A
end

function spherical_arc_length_correction(delta_u_bar::AbstractVector{<:Real},
        delta_lambda_bar::Real, delta_u_R::AbstractVector{<:Real},
        delta_u_F::AbstractVector{<:Real}, arc_length_alpha::Real, arc_radius::Real)
    length(delta_u_bar) == length(delta_u_R) == length(delta_u_F) ||
        throw(DimensionMismatch("spherical arc correction vectors must have equal length"))
    isfinite(arc_length_alpha) && arc_length_alpha > 0.0 || throw(ArgumentError(
        "spherical arc correction requires finite positive alpha, got $arc_length_alpha"))
    isfinite(arc_radius) && arc_radius > 0.0 || throw(ArgumentError(
        "spherical arc correction requires finite positive radius, got $arc_radius"))
    g = dot(delta_u_bar, delta_u_bar) +
        arc_length_alpha^2 * delta_lambda_bar^2 - arc_radius^2
    denominator = 2.0 * dot(delta_u_bar, delta_u_F) +
                  2.0 * arc_length_alpha^2 * delta_lambda_bar
    denominator_scale = 2.0 * (norm(delta_u_bar) * norm(delta_u_F) +
                                arc_length_alpha^2 * abs(delta_lambda_bar))
    abs(denominator) > eps(Float64) * max(1.0, denominator_scale) || error(
        "spherical arc correction is singular (denominator=$denominator, g=$g)")
    delta_lambda = (-g - 2.0 * dot(delta_u_bar, delta_u_R)) / denominator
    delta_u = delta_u_R .+ delta_lambda .* delta_u_F
    all(isfinite, delta_u) && isfinite(delta_lambda) || error(
        "spherical arc correction produced a non-finite update")
    return delta_u, delta_lambda, g
end

function solve_czm_basic_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param, ms::MechState; dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(ms.u_prev)
        u_start = copy(u)
        damage_states = clone_damage_states(ms.damage_states)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)

        # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
        if geo_nl
            F_thermo_chem = zeros(Float64, ndof)
            K_bulk_cached = nothing
        else
            F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
            K_bulk_cached = bulk_stiffness(czm_mesh, param)
        end
        eig_kwargs = geo_nl ? (geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress) : ()
        geom_cache = cohesive_geometry(czm_mesh)
        ws_basic = assembly_workspace(czm_mesh)

        total_iter = 0
        R_norm_0 = 1.0
        R_norm = Inf
        converged = false
        converged_R_norm = Inf
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

            R = F_ext + F_thermo_chem - f_int_total

            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            if iter == 1
                R_norm_0 = max(R_norm, 1e-10)
            end
            rel_norm = R_norm / R_norm_0

            if R_norm < tol || rel_norm < tol
                converged = true
                converged_R_norm = R_norm
                damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta=visc_beta)
                break
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            Δu = try
                K_bc \ R_bc
            catch
                break
            end

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            u, R_norm, ls_accepted, α_used = backtrack_line_search!(u, Δu, czm_mesh, param,damage_states, F_ext, F_thermo_chem, R_norm,bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws_basic;visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress)

            if !ls_accepted
                break
            end
        end

        if !converged
            u = u_start
            damage_states = clone_damage_states(ms.damage_states)   # 未收敛不触碰 ms（试探态丢弃）
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

        R = F_ext + F_thermo_chem - f_int_total

        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        R_norm = norm(R)

        # 收敛时报告收敛时的残差（损伤更新前），
        # 未收敛时报告重新组装后的残差
        final_R_norm = converged ? converged_R_norm : R_norm

        result.converged = converged
        result.iterations = total_iter
        result.residual_norm = final_R_norm
        result.displacement = u
        fill_czm_result!(result, u, damage_states, separations, tractions)
        if converged
            ms.damage_states = damage_states   # 收敛提交（D-提交语义）
            ms.u_prev = copy(u)
        end
        return result
    end

    function solve_czm_arc_length_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param, ms::MechState; dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(ms.u_prev)
        damage_states = clone_damage_states(ms.damage_states)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)

        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = bulk_stiffness(czm_mesh, param)
        geom_cache = cohesive_geometry(czm_mesh)
        ws = assembly_workspace(czm_mesh)

        # 增量载荷参考
        _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model)
        F_target = F_ext + F_thermo_chem_total
        F_delta = F_target - f_int_ref

        total_iter = 0
        load_progress = 0.0
        load_step = 0
        step_size = 1.0 / max(1, n_load_steps)
        step_size_min = step_size / 128.0
        step_size_max = step_size
        last_residual = Inf
        converged_substep = false
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        while load_progress < 1.0 - 1e-12
            load_step += 1
            load_start = load_progress
            target_progress = min(1.0, load_start + step_size)

            u_start = copy(u)
            damage_start = clone_damage_states(damage_states)   # 子步试探回滚（ms 不受影响）
            converged_substep = false
            last_residual = Inf

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model)

            F_load_bc = copy(F_delta)
            zero_czm_bc_entries!(F_load_bc, bc_dofs)

            F_applied = f_int_ref + load_start * F_delta
            R = F_applied - f_int_total
            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end
            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            tangent = try
                K_bc \ F_load_bc
            catch
                nothing
            end

            if tangent === nothing || any(isnan, tangent) || any(isinf, tangent)
                break
            end

            delta_lambda_pred = target_progress - load_start
            if delta_lambda_pred <= 0.0
                delta_lambda_pred = step_size
            end

            delta_u_pred = tangent * delta_lambda_pred
            arc_target = sqrt(sum(abs2, delta_u_pred))
            if !isfinite(arc_target) || arc_target <= 0.0
                break
            end

            u = u_start + delta_u_pred
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            load_progress = target_progress

            for iter in 1:max_iter
                total_iter += 1

                F_applied = f_int_ref + load_progress * F_delta
                K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model)

                R = F_applied - f_int_total
                for (dof, val) in zip(bc_dofs, bc_vals)
                    R[dof] = val - u[dof]
                end

                delta_u = u - u_start
                delta_lambda = load_progress - load_start
                arc_constraint = dot(delta_u, delta_u) - arc_target^2
                residual_norm = sqrt(norm(R)^2 + arc_constraint^2)
                last_residual = residual_norm

                substep_tol = tol * 10.0
                if norm(R) < substep_tol && abs(arc_constraint) < substep_tol
                    converged_substep = true
                    load_progress = min(load_progress, 1.0)
                    step_size = min(step_size * 1.25, step_size_max)
                    break
                end

                K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

                # Crisfield cylindrical arc-length: solve two linear systems
                delta_u_R = try
                    K_bc \ R_bc
                catch
                    nothing
                end
                if delta_u_R === nothing || any(isnan, delta_u_R) || any(isinf, delta_u_R)
                    break
                end

                delta_u_F = try
                    K_bc \ F_load_bc
                catch
                    nothing
                end
                if delta_u_F === nothing || any(isnan, delta_u_F) || any(isinf, delta_u_F)
                    break
                end

                # Quadratic coefficients for ||delta_u + delta_u_R + dl * delta_u_F||^2 = arc_target^2
                du_bar = delta_u + delta_u_R
                a_q = dot(delta_u_F, delta_u_F)
                b_q = 2.0 * dot(du_bar, delta_u_F)
                c_q = dot(du_bar, du_bar) - arc_target^2

                discriminant = b_q^2 - 4.0 * a_q * c_q
                if discriminant < 0.0
                    break
                end

                sqrt_disc = sqrt(discriminant)
                dl1 = (-b_q + sqrt_disc) / (2.0 * a_q)
                dl2 = (-b_q - sqrt_disc) / (2.0 * a_q)

                # Pick root closest to predictor direction
                du_new_1 = du_bar + dl1 * delta_u_F
                du_new_2 = du_bar + dl2 * delta_u_F
                dot1 = dot(du_new_1, delta_u_pred)
                dot2 = dot(du_new_2, delta_u_pred)
                delta_lambda_corr = dot1 >= dot2 ? dl1 : dl2

                delta_u_corr = delta_u_R + delta_lambda_corr * delta_u_F

                if any(isnan, delta_u_corr) || any(isinf, delta_u_corr)
                    break
                end

                u = u + delta_u_corr
                load_progress = load_progress + delta_lambda_corr
                apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            end

            if !converged_substep
                u = u_start
                damage_states = damage_start
                load_progress = load_start
                step_size *= 0.5

                if step_size < step_size_min
                    @warn "CZM arc-length stepping stalled" load_progress=load_progress target_progress=target_progress residual=last_residual step_size=step_size
                    break
                end

                @debug "Arc-length substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_residual step_size=step_size
                continue
            end
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(
            czm_mesh, u, param;
            damage_states=damage_states, K_bulk_cached=K_bulk_cached,
            geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model)

        F_applied_final = f_int_ref + load_progress * F_delta
        R = F_applied_final - f_int_total
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        R_norm = norm(R)

        final_tol = tol * 100.0
        result.converged = load_progress >= 1.0 - 1e-12 && converged_substep
        result.iterations = total_iter
        result.residual_norm = R_norm
        result.displacement = u

        # 所有子步完成后，统一更新损伤（与 basic 方法一致）
        if result.converged
            damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta=visc_beta)
        end

        fill_czm_result!(result, u, damage_states, separations, tractions)
        if result.converged
            ms.damage_states = damage_states   # 收敛提交
            ms.u_prev = copy(u)
        end
        return result
    end

# ========================================================================
# 7. Newton-Raphson solver
# ========================================================================

"""
    newton_raphson_czm(czm_mesh, F_ext, param; dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,max_iter=50, tol=1e-8, u0=nothing, n_load_steps=10)

Newton-Raphson nonlinear solver with load substeps.

# Returns
- `result`: CZMResult
- `new_czm_mesh`: updated CZM mesh with damage states
"""
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param, ms::MechState; dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)

    u = copy(ms.u_prev)
    damage_states = clone_damage_states(ms.damage_states)

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)

    # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
    if geo_nl
        F_thermo_chem_total = zeros(Float64, ndof)
        K_bulk_cached = nothing
    else
        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = bulk_stiffness(czm_mesh, param)
    end
    eig_kwargs = geo_nl ? (geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress) : ()
    geom_cache = cohesive_geometry(czm_mesh)
    ws = assembly_workspace(czm_mesh)

    # 增量载荷参考：u_prev 近似在上一时间步的平衡态
    # f_int(u_prev) ≈ 上一步外力，F_delta = 目标载荷 - 平衡内力（增量，通常很小）
    _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)
    F_target = F_ext + F_thermo_chem_total
    F_delta = F_target - f_int_ref

    total_iter = 0
    load_progress = 0.0
    load_step = 0
    step_size = 1.0 / max(1, n_load_steps)
    step_size_min = step_size / 128.0
    step_size_max = step_size

    last_R_norm = Inf
    converged_substep = false

    while load_progress < 1.0 - 1e-12
        load_step += 1
        target_progress = min(1.0, load_progress + step_size)
        # 从平衡态 f_int_ref 逐步推进到目标态 F_target
        F_applied = f_int_ref + target_progress * F_delta

        u_start = copy(u)
        damage_start = clone_damage_states(damage_states)
        converged_substep = false
        last_R_norm = Inf

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, param;
                damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

            R = F_applied - f_int_total

            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            last_R_norm = R_norm
            substep_tol = tol * 10.0

            if R_norm < substep_tol
                converged_substep = true
                load_progress = target_progress
                step_size = min(step_size * 1.25, step_size_max)
                break
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
            Δu = K_bc \ R_bc

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            α = 1.0
            ls_accepted = false
            for _ in 1:8
                u_trial = u + α * Δu
                for (dof, val) in zip(bc_dofs, bc_vals)
                    u_trial[dof] = val
                end

                _, f_int_trial, _, _ = assemble_coupled_system(
                    czm_mesh, u_trial, param;
                    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                    geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

                R_trial = F_applied - f_int_trial
                for (dof, val) in zip(bc_dofs, bc_vals)
                    R_trial[dof] = val - u_trial[dof]
                end

                R_trial_norm = norm(R_trial)
                if !isnan(R_trial_norm) && R_trial_norm < R_norm
                    ls_accepted = true
                    break
                end

                α *= 0.5
            end

            if !ls_accepted
                break
            end

            u = u + α * Δu

            for (dof, val) in zip(bc_dofs, bc_vals)
                u[dof] = val
            end
        end

        if !converged_substep
            u = u_start
            damage_states = damage_start
            step_size *= 0.5

            if step_size < step_size_min
                @warn "CZM adaptive load stepping stalled" load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
                break
            end

            @debug "Load substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
            continue
        end
    end

    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

    F_applied_final = f_int_ref + load_progress * F_delta
    R = F_applied_final - f_int_total
    for (dof, val) in zip(bc_dofs, bc_vals)
        R[dof] = val - u[dof]
    end
    R_norm = norm(R)

    final_tol = tol * 100.0

    result.converged = load_progress >= 1.0 - 1e-12 && R_norm < final_tol
    result.iterations = total_iter
    result.residual_norm = R_norm
    result.displacement = u

    # 所有子步完成后，统一更新损伤（与 basic 方法一致：冻结损伤求解位移，收敛后更新）
    if result.converged
        damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta=visc_beta)
    end

    for i in 1:n_coh
        result.damage[i] = damage_states[i].D
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end

    if result.converged
        ms.damage_states = damage_states   # 收敛提交
        ms.u_prev = copy(u)
    end
    return result
end

"""
    solve_czm_arc_geo_step(czm_mesh, F_ext, param, ms; ...) -> CZMResult

Crisfield 球面弧长 geo 路径（Theory §6.10；λ 缩放本征应变增量）。
约束为 `‖Δu‖² + α²Δλ² = Δl²`。平衡残差采用代码约定
`R = F_ext - f_int`，故修正分解为 `δu = K⁻¹R + δλ K⁻¹f̂`，
其中 `f̂ = ∂R/∂λ = -∂f_int/∂λ` 在每个当前迭代状态重新差分。
失败子步回滚位移/λ；损伤和塑性状态只在最终 `λ=1` 平衡验收后提交。
"""
function solve_czm_arc_geo_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64},
        param, ms::MechState;
        dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,
        max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10,
        arc_length_alpha::Float64=1.0,
        visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true,
        eigenstrain=nothing, eigenstrain_ref=nothing,
        plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    eigenstrain === nothing && error("solve_czm_arc_geo_step: geo 弧长需要 eigenstrain（λ 的缩放对象）")
    isfinite(arc_length_alpha) && arc_length_alpha > 0.0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: arc_length_alpha must be finite and positive, got $arc_length_alpha"))
    max_iter > 0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: max_iter must be positive, got $max_iter"))
    n_load_steps > 0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: n_load_steps must be positive, got $n_load_steps"))
    isfinite(tol) && tol > 0.0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: tol must be finite and positive, got $tol"))
    ndof = 2 * czm_mesh.nnode
    n_coh = czm_mesh.n_cohesive
    result = CZMResult(ndof, n_coh)
    u = copy(ms.u_prev)
    damage_states = clone_damage_states(ms.damage_states)
    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)
    zero_bc_vals = zeros(Float64, length(bc_vals))
    apply_czm_dirichlet!(u, bc_dofs, bc_vals)
    ws = assembly_workspace(czm_mesh)
    geom_cache = cohesive_geometry(czm_mesh)

    eig_kw = (geo_nl=true, plasticity=plasticity, mech_state=mech_state, prestress=prestress)
    mix(lam) = (dT=eigenstrain_ref === nothing ? lam .* eigenstrain.dT :
                    eigenstrain_ref.dT .+ lam .* (eigenstrain.dT .- eigenstrain_ref.dT),
                Δsn=eigenstrain_ref === nothing ? lam .* eigenstrain.Δsn :
                    eigenstrain_ref.Δsn .+ lam .* (eigenstrain.Δsn .- eigenstrain_ref.Δsn),
                Δsp=eigenstrain_ref === nothing ? lam .* eigenstrain.Δsp :
                    eigenstrain_ref.Δsp .+ lam .* (eigenstrain.Δsp .- eigenstrain_ref.Δsp))
    function assemble_at(ul, lam)
        return assemble_coupled_system(czm_mesh, ul, param;
            damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model,
            eig_kw..., eigenstrain=mix(lam))
    end
    function residual_at(f_int, ul)
        R = F_ext .- f_int
        for dof in bc_dofs
            R[dof] = 0.0
        end
        return R
    end
    function load_direction_at(ul, lam)
        hλ = cbrt(eps(Float64)) * max(1.0, abs(lam))
        f_plus = copy(assemble_at(ul, lam + hλ)[2])
        f_minus = copy(assemble_at(ul, lam - hλ)[2])
        f_hat = -(f_plus .- f_minus) ./ (2.0 * hλ)
        f_hat[bc_dofs] .= 0.0
        all(isfinite, f_hat) || error(
            "solve_czm_arc_geo_step: non-finite load direction at λ=$lam")
        return f_hat
    end

    # λ=0 参考态可能含卷绕预应力；自由芯部下通常需要多次非线性 Newton，
    # 不能只做一次修正后就判失败。损伤/塑性在此仍为 trial，不提交历史状态。
    K0 = spzeros(Float64, ndof, ndof)
    f0 = zeros(Float64, ndof)
    R0 = fill(Inf, ndof)
    reference_iter = 0
    reference_converged = false
    zero_external = zeros(Float64, ndof)
    for _ in 1:max_iter
        K0, f0, _, _ = assemble_at(u, 0.0)
        R0 = residual_at(f0, u)
        R0_norm = norm(R0)
        isfinite(R0_norm) || error(
            "solve_czm_arc_geo_step: non-finite reference residual")
        if R0_norm <= tol
            reference_converged = true
            break
        end

        K0_bc, R0_bc = apply_bc_czm(
            K0, R0; bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        delta_u = K0_bc \ R0_bc
        all(isfinite, delta_u) || error(
            "solve_czm_arc_geo_step: non-finite reference-state Newton correction")
        u_trial, _, accepted, _ = backtrack_line_search!(
            u, delta_u, czm_mesh, param, damage_states,
            F_ext, zero_external, R0_norm, bc_dofs, bc_vals,
            nothing, geom_cache, ws;
            max_halvings=12, visc_beta=visc_beta, czm_model=czm_model, geo_nl=true,
            eigenstrain=mix(0.0), plasticity=plasticity,
            mech_state=mech_state, prestress=prestress)
        accepted || error(
            "solve_czm_arc_geo_step: reference-state line search failed " *
            "(iteration=$(reference_iter + 1), residual=$R0_norm)")
        u .= u_trial
        reference_iter += 1
    end
    if !reference_converged
        K0, f0, _, _ = assemble_at(u, 0.0)
        R0 = residual_at(f0, u)
        reference_converged = norm(R0) <= tol
    end
    reference_converged || error(
        "solve_czm_arc_geo_step: reference state is not in equilibrium " *
        "after $max_iter iterations (residual=$(norm(R0)))")
    f_hat0 = load_direction_at(u, 0.0)
    K0_bc, _ = apply_bc_czm(K0, zeros(Float64, ndof);
                            bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
    tangent0 = K0_bc \ f_hat0
    tangent_norm0 = sqrt(dot(tangent0, tangent0) + arc_length_alpha^2)
    isfinite(tangent_norm0) && tangent_norm0 > 0.0 || error(
        "solve_czm_arc_geo_step: invalid initial augmented tangent norm $tangent_norm0")

    λ = 0.0
    Δl0 = tangent_norm0 / n_load_steps
    Δl = Δl0
    Δl_min = Δl0 / 128.0
    previous_tangent = Vector{Float64}()
    total_iter = reference_iter
    step_count = 0
    residual_history = Float64[]
    lambda_history = Float64[λ]
    step_history = Float64[]
    arc_step_tol = 10.0 * tol  # 中间路径点容差；最终 λ=1 平衡仍严格使用 tol

    while λ < 1.0 - 1e-10
        step_count += 1
        step_count <= 10000 || error(
            "solve_czm_arc_geo_step: exceeded 10000 arc steps " *
            "(λ_history=$lambda_history, residual_history=$residual_history)")
        u_start = copy(u)
        λ_start = λ

        K_start, _, _, _ = assemble_at(u_start, λ_start)
        f_hat_start = load_direction_at(u_start, λ_start)
        K_start_bc, _ = apply_bc_czm(K_start, zeros(Float64, ndof);
                                     bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        tangent = K_start_bc \ f_hat_start
        augmented_norm = sqrt(dot(tangent, tangent) + arc_length_alpha^2)
        isfinite(augmented_norm) && augmented_norm > 0.0 || error(
            "solve_czm_arc_geo_step: invalid augmented tangent at λ=$λ_start")

        direction = 1.0
        if !isempty(previous_tangent)
            augmented_tangent = vcat(tangent, arc_length_alpha)
            dot(augmented_tangent, previous_tangent) < 0.0 && (direction = -1.0)
        end
        dλ_pred = direction * Δl / augmented_norm
        if direction > 0.0 && dλ_pred > 1.0 - λ_start
            dλ_pred = 1.0 - λ_start
        end
        Δl_step = abs(dλ_pred) * augmented_norm
        Δl_step > 0.0 || error(
            "solve_czm_arc_geo_step: zero predictor step at λ=$λ_start")
        u .= u_start .+ dλ_pred .* tangent
        λ = λ_start + dλ_pred
        apply_czm_dirichlet!(u, bc_dofs, bc_vals)

        step_ok = false
        last_residual = Inf
        last_step_error = nothing
        for _ in 1:max_iter
            total_iter += 1
            try
                K, f_int, _, _ = assemble_at(u, λ)
                R = residual_at(f_int, u)
                Δu_bar = u .- u_start
                Δλ_bar = λ - λ_start
                g = dot(Δu_bar, Δu_bar) +
                    arc_length_alpha^2 * Δλ_bar^2 - Δl_step^2
                last_residual = norm(R)
                push!(residual_history, last_residual)
                if last_residual <= arc_step_tol &&
                   abs(g) <= arc_step_tol * max(Δl_step^2, eps(Float64))
                    step_ok = true
                    break
                end

                K_bc, R_bc = apply_bc_czm(
                    K, R; bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
                f_hat = load_direction_at(u, λ)
                delta_u_R = K_bc \ R_bc
                delta_u_F = K_bc \ f_hat
                Δu, dλ, _ = spherical_arc_length_correction(
                    Δu_bar, Δλ_bar, delta_u_R, delta_u_F,
                    arc_length_alpha, Δl_step)
                u .+= Δu
                λ += dλ
                apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            catch err
                err isa ErrorException || rethrow()
                last_step_error = sprint(showerror, err)
                break
            end
        end

        if step_ok
            Δu_step = u .- u_start
            Δλ_step = λ - λ_start
            previous_tangent = vcat(Δu_step, arc_length_alpha * Δλ_step)
            previous_norm = norm(previous_tangent)
            previous_norm > 0.0 || error(
                "solve_czm_arc_geo_step: accepted a zero arc step at λ=$λ")
            previous_tangent ./= previous_norm
            push!(lambda_history, λ)
            push!(step_history, Δl_step)
            continue
        end

        u .= u_start
        λ = λ_start
        Δl *= 0.5
        if Δl < Δl_min
            error("solve_czm_arc_geo_step: arc stepping stalled at λ=$λ " *
                  "(residual=$last_residual, Δl=$Δl, λ_history=$lambda_history, " *
                  "step_error=$(repr(last_step_error)), step_history=$step_history, " *
                  "residual_history=$residual_history)")
        end
    end

    # 精确落到 λ=1，并以固定 λ Newton 复算最终平衡；不能用弧长预测残差代替。
    λ = 1.0
    final_residual = Inf
    separations = Vector{Tuple{Float64,Float64}}(undef, n_coh)
    tractions = Vector{Tuple{Float64,Float64}}(undef, n_coh)
    for _ in 1:max_iter
        total_iter += 1
        K, f_int, separations, tractions = assemble_at(u, λ)
        R = residual_at(f_int, u)
        final_residual = norm(R)
        push!(residual_history, final_residual)
        final_residual <= tol && break
        K_bc, R_bc = apply_bc_czm(K, R;
            bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        Δu = K_bc \ R_bc
        all(isfinite, Δu) || error(
            "solve_czm_arc_geo_step: non-finite final λ=1 correction")
        u .+= Δu
        apply_czm_dirichlet!(u, bc_dofs, bc_vals)
    end
    final_residual <= tol || error(
        "solve_czm_arc_geo_step: final λ=1 equilibrium did not converge " *
        "(residual=$final_residual, λ_history=$lambda_history, " *
        "step_history=$step_history, residual_history=$residual_history)")

    _, f_final, separations, tractions = assemble_at(u, 1.0)
    final_residual = norm(residual_at(f_final, u))
    final_residual <= tol || error(
        "solve_czm_arc_geo_step: final residual recomputation failed (residual=$final_residual)")
    damage_states = update_damage_per_interface(
        czm_mesh, damage_states, separations, param, czm_model; visc_beta=visc_beta)
    result.converged = true
    result.iterations = total_iter
    result.residual_norm = final_residual
    fill_czm_result!(result, u, damage_states, separations, tractions)
    ms.damage_states = damage_states   # λ=1 平衡验收后提交
    ms.u_prev = copy(u)
    return result
end

"""
    solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; dT/Δsn/Δsp...) -> CZMResult

CZM 单步统一入口（2026-08-30 终态签名）：求解配置从 `czm_opt::CzmOptions` 展开，
演化状态在 `ms::MechState` 上收敛提交（失败/试探不触碰 ms）。
"""
function solve_czm_step(czm_mesh::CohesiveMesh, ms::MechState, param, F_ext::Vector{Float64},
        czm_opt::CzmOptions;
        dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,
        eigenstrain=nothing, mech_state=nothing, prestress=nothing)
    method = lowercase(czm_opt.iter_method)
    max_iter = czm_opt.max_iter
    tol = czm_opt.tol
    n_load_steps = czm_opt.load_steps
    arc_length_alpha = czm_opt.arc_length_alpha
    visc_beta = 1.0
    if czm_opt.viscous_enabled && czm_opt.viscous_tau > 0.0
        tau_visc = czm_opt.viscous_tau / param.scale.t0
        delta_s = method == "basic" ? 1.0 : 1.0 / max(1, n_load_steps)
        visc_beta = delta_s / (tau_visc + delta_s)
    end
    czm_model = czm_opt.model
    fix_inner = czm_opt.fix_inner
    geo_nl = czm_opt.geo_nonlinear

    if method == "load_substep"
        return newton_raphson_czm(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, mech_state=mech_state, prestress=prestress)
    elseif method == "basic"
        return solve_czm_basic_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, mech_state=mech_state, prestress=prestress)
    elseif (method == "arc_length" || method == "arclength" || method == "arc-length") && geo_nl
        return solve_czm_arc_geo_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
            arc_length_alpha=arc_length_alpha,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, mech_state=mech_state, prestress=prestress)
    elseif method == "arc_length" || method == "arclength" || method == "arc-length"
        return solve_czm_arc_length_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner)
    else
        error("Unknown CZM iteration method: $(czm_opt.iter_method). Use 'basic', 'load_substep', or 'arc_length'.")
    end
end
