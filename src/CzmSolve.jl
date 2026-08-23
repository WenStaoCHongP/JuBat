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

function clone_czm_mesh_with_damage(czm_mesh::CohesiveMesh, damage_states::AbstractVector{<:AbstractDamageState})
    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = damage_states
    return new_czm_mesh
end

"""
    update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta=1.0)

按 cohesive 单元的 interface_type 分组，分别调用 `update_damage`。
当所有单元属于同一界面时，退化为单次 `update_damage` 调用。
"""
function update_damage_per_interface(czm_mesh::CohesiveMesh, damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, param_cache::CzmParamCache; visc_beta::Float64=1.0)
    n_coh = czm_mesh.n_cohesive
    @assert length(damage_states) == n_coh "damage_states length mismatch"
    @assert length(separations) == n_coh "separations length mismatch"

    # 按 interface_type 分批（保持原始顺序）
    new_states = Vector{DamageState}(undef, n_coh)
    for iface in keys(param_cache.by_interface)
        params = param_cache.by_interface[iface]
        idx = findall(i -> czm_mesh.cohesive_elements[i].interface_type == iface, 1:n_coh)
        isempty(idx) && continue
        ds_sub = damage_states[idx]
        sep_sub = separations[idx]
        updated = update_damage(ds_sub, sep_sub, params; visc_beta=visc_beta)
        for (k, i) in enumerate(idx)
            new_states[i] = updated[k]
        end
    end
    return new_states
end

"""
    extract_bc_dofs(czm_mesh, param; cache=nothing)

从 czm_mesh 提取 Dirichlet BC 的自由度列表和对应值。
优先使用缓存中的 bc_dofs/bc_vals，否则从 identify_bc_nodes_czm 重新计算。
"""
function extract_bc_dofs(czm_mesh::CohesiveMesh, param; cache::Union{Nothing, CZMAssemblyCache}=nothing, fix_inner::Bool=true)
    if cache !== nothing
        return cache.bc_dofs, cache.bc_vals
    end
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
    backtrack_line_search!(u, Δu, czm_mesh, param_cache, damage_states, F_ext, F_thermo_chem, R_norm_current, bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws; max_halvings=8)

回溯线搜索（零化式 BC 残差）。仅用于 solve_czm_basic_step。
返回 (u_new, R_new_norm, accepted, α_used)。
accepted 时返回的 u_new 已含 BC 赋值，外部无需再执行 u = u + α*Δu。
未 accepted 时返回原始 u（未修改），外部应 break。
"""
function backtrack_line_search!(u::Vector{Float64}, Δu::Vector{Float64},czm_mesh::CohesiveMesh, param_cache::CzmParamCache, damage_states,F_ext::Vector{Float64}, F_thermo_chem::Vector{Float64},R_norm_current::Float64,bc_dofs::Vector{Int64}, bc_vals::Vector{Float64},K_bulk_cached, geom_cache, ws;max_halvings::Int=8, visc_beta::Float64=1.0, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    α = 1.0
    for _ in 1:max_halvings
        u_trial = u + α * Δu
        apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

        _, f_int_trial, _, _ = assemble_coupled_system(czm_mesh, u_trial, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, geo_nl=geo_nl, eigenstrain=eigenstrain,
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

function solve_czm_basic_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param_cache::CzmParamCache, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, cache::Union{Nothing, CZMAssemblyCache}=nothing, visc_beta::Float64=1.0, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(u_prev)
        u_start = copy(u)
        damage_states = czm_mesh.damage_states
        damage_start = clone_damage_states(damage_states)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

        # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
        if geo_nl
            F_thermo_chem = zeros(Float64, ndof)
            K_bulk_cached = nothing
        else
            F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
            K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
        end
        eig_kwargs = geo_nl ? (geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress) : ()
        geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
        ws_basic = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(ndof, n_coh)

        total_iter = 0
        R_norm_0 = 1.0
        R_norm = Inf
        converged = false
        converged_R_norm = Inf
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta, eig_kwargs...)

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
                damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta=visc_beta)
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

            u, R_norm, ls_accepted, α_used = backtrack_line_search!(u, Δu, czm_mesh, param_cache,damage_states, F_ext, F_thermo_chem, R_norm,bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws_basic;visc_beta=visc_beta, geo_nl=geo_nl, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress)

            if !ls_accepted
                break
            end
        end

        if !converged
            u = u_start
            damage_states = damage_start
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta, eig_kwargs...)

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

        new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
        return result, new_czm_mesh
    end

    function solve_czm_arc_length_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param_cache::CzmParamCache, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, cache::Union{Nothing, CZMAssemblyCache}=nothing, visc_beta::Float64=1.0)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(u_prev)
        damage_states = czm_mesh.damage_states

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
        geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
        ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(ndof, n_coh)

        # 增量载荷参考
        _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)
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
            damage_start = clone_damage_states(damage_states)
            converged_substep = false
            last_residual = Inf

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

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
                K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

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
            czm_mesh, u, param_cache;
            damage_states=damage_states, K_bulk_cached=K_bulk_cached,
            geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

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
            damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta=visc_beta)
        end

        fill_czm_result!(result, u, damage_states, separations, tractions)

        new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
        return result, new_czm_mesh
    end

# ========================================================================
# 7. Newton-Raphson solver
# ========================================================================

"""
    newton_raphson_czm(czm_mesh, F_ext, param_cache, param;α_eff=0.0, β_n=0.0, β_p=0.0,dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,max_iter=50, tol=1e-8, u0=nothing, n_load_steps=10)

Newton-Raphson nonlinear solver with load substeps.

# Returns
- `result`: CZMResult
- `new_czm_mesh`: updated CZM mesh with damage states
"""
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param_cache::CzmParamCache, param; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, u0::Union{Vector{Float64},Nothing}=nothing, n_load_steps::Int=10, cache::Union{Nothing, CZMAssemblyCache}=nothing, visc_beta::Float64=1.0, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)

    u = u0 === nothing ? zeros(Float64, ndof) : copy(u0)
    damage_states = czm_mesh.damage_states

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

    # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
    if geo_nl
        F_thermo_chem_total = zeros(Float64, ndof)
        K_bulk_cached = nothing
    else
        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
    end
    eig_kwargs = geo_nl ? (geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity, mech_state=mech_state, prestress=prestress) : ()
    geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
    ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(2 * nnode, czm_mesh.n_cohesive)

    # 增量载荷参考：u_prev 近似在上一时间步的平衡态
    # f_int(u_prev) ≈ 上一步外力，F_delta = 目标载荷 - 平衡内力（增量，通常很小）
    _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, eig_kwargs...)
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
                czm_mesh, u, param_cache;
                damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, eig_kwargs...)

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
                    czm_mesh, u_trial, param_cache;
                    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                    geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, eig_kwargs...)

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

    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, eig_kwargs...)

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
        damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta=visc_beta)
    end

    for i in 1:n_coh
        result.damage[i] = damage_states[i].D
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end

    new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)

    return result, new_czm_mesh
end

"""
    solve_czm_arc_geo_step(czm_mesh, F_ext, param_cache, param, u_prev; ...) -> (CZMResult, CohesiveMesh)

Crisfield 柱面弧长 geo 路径（Batch 5，Theory §6.10；λ 缩放本征应变增量）。
ε*(λ) = ε_ref + λ·Δε*（ε_ref 缺省零向量）；f̂ = ∂f_int/∂λ（数值差分一步，f_int 对 λ 线性）。
增广 Newton：Δu = Δu_R + Δλ·Δu_F，柱面约束 (6.89) → Δλ 二次方程 (6.90–6.91)，
取与上步切向内积较大者；无实根 Δl/2 重试；λ 推进至 ≥1 且残差 < tol 提交；
步长下限 Δl/128 失败即报错终止（不伪造收敛）。
"""
function solve_czm_arc_geo_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64},
        param_cache::CzmParamCache, param, u_prev::Vector{Float64};
        α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0,
        dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,
        max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10,
        cache::Union{Nothing, CZMAssemblyCache}=nothing, visc_beta::Float64=1.0,
        eigenstrain=nothing, eigenstrain_ref=nothing,
        plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    eigenstrain === nothing && error("solve_czm_arc_geo_step: geo 弧长需要 eigenstrain（λ 的缩放对象）")
    ndof = 2 * czm_mesh.nnode
    n_coh = czm_mesh.n_cohesive
    result = CZMResult(ndof, n_coh)
    u = copy(u_prev)
    damage_states = czm_mesh.damage_states
    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)
    ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(ndof, n_coh)
    geom_cache = cache !== nothing ? cache.cohesive_geom : nothing

    eig_kw = (geo_nl=true, plasticity=plasticity, mech_state=mech_state, prestress=prestress)
    mix(lam) = (α_eff=eigenstrain.α_eff, β_n=eigenstrain.β_n, β_p=eigenstrain.β_p,
                dT=eigenstrain_ref === nothing ? lam .* eigenstrain.dT :
                    eigenstrain_ref.dT .+ lam .* (eigenstrain.dT .- eigenstrain_ref.dT),
                Δsn=eigenstrain_ref === nothing ? lam .* eigenstrain.Δsn :
                    eigenstrain_ref.Δsn .+ lam .* (eigenstrain.Δsn .- eigenstrain_ref.Δsn),
                Δsp=eigenstrain_ref === nothing ? lam .* eigenstrain.Δsp :
                    eigenstrain_ref.Δsp .+ lam .* (eigenstrain.Δsp .- eigenstrain_ref.Δsp))
    f_int_at(ul, lam) = assemble_coupled_system(czm_mesh, ul, param_cache;
        damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
        eig_kw..., eigenstrain=mix(lam))[2]
    hλ = 1e-8
    f_hat = -(f_int_at(u, hλ) .- f_int_at(u, 0.0)) ./ hλ   # f̂ = −∂R/∂λ（6.89；正号则平衡点在 −λΔu_F 侧，轨迹震荡）
    f_hat[bc_dofs] .= 0.0               # BC 自由度不参与 λ 方向（否则约束步漂移固定位移）
    # Δl 按解尺度初始化：‖K⁻¹·f̂‖/n_load_steps（满载切向位移长度；解位移 ~1e-3 时
    # 固定 Δl=1/n 会远超解范数，二次方程无合理根致长时间不收敛）
    K0, _ = assemble_coupled_system(czm_mesh, u, param_cache;
        damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
        eig_kw..., eigenstrain=mix(0.0))
    K0_bc, _ = apply_bc_czm(K0, zeros(Float64, ndof); bc_dofs=bc_dofs, bc_vals=bc_vals)
    Δl0 = norm(K0_bc \ f_hat) / max(1, n_load_steps)
    Δl0 > 0 || (Δl0 = 1e-6)

    λ = 0.0
    Δl = Δl0
    Δl_min = Δl / 128
    t_prev = zeros(Float64, ndof)
    total_iter = 0
    converged = false
    while λ < 1.0 - 1e-12
        ū = zeros(Float64, ndof)
        λ_step0 = λ
        step_ok = false
        for _ in 1:max_iter
            total_iter += 1
            K, f_int = assemble_coupled_system(czm_mesh, u, param_cache;
                damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
                eig_kw..., eigenstrain=mix(λ))
            R = F_ext .- f_int
            for dof in bc_dofs
                R[dof] = 0.0
            end
            # λ=0 处残差恒零（u=0、ε*=0 平凡平衡）——须本子步 λ 已推进才可判收敛，
            # 否则外层子步死循环（每步秒收敛而 λ 不动）
            norm(R) < 10 * tol && λ > λ_step0 + 1e-14 && (step_ok = true; break)   # 子步判据 10×tol（load_substep 惯例；残差略超 tol 时 K 病态方向放大 Δu_R 破坏约束根）
            K_bc, R_bc = apply_bc_czm(K, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
            Δu_R = -(K_bc \ R_bc)
            Δu_F = K_bc \ f_hat
            b = ū .+ Δu_R
            a2 = dot(Δu_F, Δu_F); a1 = 2 * dot(b, Δu_F); a0 = dot(b, b) - Δl^2
            disc = a2 == 0 ? 0.0 : a1 * a1 - 4 * a2 * a0
            disc < 0 && break
            sq = sqrt(disc)
            r1 = (-a1 + sq) / (2 * a2); r2 = (-a1 - sq) / (2 * a2)
            dλ = norm(t_prev) < 1e-30 ? max(r1, r2) :
                 (dot(t_prev, Δu_R .+ r1 .* Δu_F) ≥ dot(t_prev, Δu_R .+ r2 .* Δu_F) ? r1 : r2)
            # 根合理性守卫：期望步长 Δλ_e = Δl/‖Δu_F‖；大幅回跳（<−0.05Δλ_e）或
            # 停滞（|dλ|<1e-6Δλ_e）说明病态 Δu_R 污染二次方程——按子步失败处理（回退减半）
            dλ_e = Δl / max(norm(Δu_F), 1e-30)
            (dλ > -0.05 * dλ_e && abs(dλ) > 1e-6 * dλ_e) || break
            Δu = Δu_R .+ dλ .* Δu_F
            u .+= Δu
            ū .+= Δu
            λ += dλ
            t_prev = Δu
            λ ≥ 1.0 - 1e-12 && (step_ok = true; break)
            f_int2 = assemble_coupled_system(czm_mesh, u, param_cache;
                damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
                eig_kw..., eigenstrain=mix(λ))[2]
            R2 = F_ext .- f_int2
            for dof in bc_dofs
                R2[dof] = 0.0
            end
            norm(R2) < 10 * tol && (step_ok = true; break)
        end
        if step_ok && λ ≥ 1.0 - 1e-12
            converged = true
            break
        elseif step_ok
            norm(ū) > 1e-30 && (t_prev = ū ./ norm(ū))
            continue
        end
        # 本步失败：回退本步位移与 λ，弧长减半重试
        u .-= ū
        λ = λ_step0
        Δl /= 2
        if Δl < Δl_min
            @warn "CZM geo arc-length stepping stalled" λ=λ Δl=Δl
            break
        end
    end
    if converged
        _, _, separations, tractions = assemble_coupled_system(czm_mesh, u, param_cache;
            damage_states=damage_states, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
            eig_kw..., eigenstrain=mix(1.0))
        damage_states = update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta=visc_beta)
        fill_czm_result!(result, u, damage_states, separations, tractions)
        result.iterations = total_iter
        result.residual_norm = 0.0
    else
        fill_czm_result!(result, u, damage_states,
                         [(0.0, 0.0) for _ in 1:n_coh], [(0.0, 0.0) for _ in 1:n_coh])
        result.iterations = total_iter
    end
    new_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
    return result, new_mesh
end

"""
    solve_czm_step(czm_mesh, F_ext, param_cache, param, u_prev; ...)

Solve a single CZM step with selectable iteration method.
"""
function solve_czm_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param_cache::CzmParamCache, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, iter_method::String="load_substep", cache::Union{Nothing, CZMAssemblyCache}=nothing, visc_beta::Float64=1.0, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, mech_state=nothing, prestress=nothing)
    method = lowercase(iter_method)

    if method == "load_substep"
        return newton_raphson_czm(czm_mesh, F_ext, param_cache, param;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, u0=u_prev, n_load_steps=n_load_steps, cache=cache,
            visc_beta=visc_beta, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=plasticity, mech_state=mech_state, prestress=prestress)
    elseif method == "basic"
        return solve_czm_basic_step(czm_mesh, F_ext, param_cache, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, cache=cache,
            visc_beta=visc_beta, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=plasticity, mech_state=mech_state, prestress=prestress)
    elseif (method == "arc_length" || method == "arclength" || method == "arc-length") && geo_nl
        return solve_czm_arc_geo_step(czm_mesh, F_ext, param_cache, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, cache=cache,
            visc_beta=visc_beta, eigenstrain=eigenstrain,
            plasticity=plasticity, mech_state=mech_state, prestress=prestress)    elseif method == "arc_length" || method == "arclength" || method == "arc-length"
        return solve_czm_arc_length_step(czm_mesh, F_ext, param_cache, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, cache=cache,
            visc_beta=visc_beta)
    else
        error("Unknown CZM iteration method: $(iter_method). Use 'basic', 'load_substep', or 'arc_length'.")
    end
end
