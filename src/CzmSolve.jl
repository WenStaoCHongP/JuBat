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

function clone_damage_state(s::DamageState)
    new_state = DamageState()
    new_state.D = s.D
    new_state.D_visc = s.D_visc
    new_state.δ_max_n = s.δ_max_n
    new_state.δ_max_t = s.δ_max_t
    new_state.δ_max_eff = s.δ_max_eff
    new_state.fractured = s.fractured
    new_state.accumulated_damage = s.accumulated_damage
    return new_state
end

function clone_damage_states(damage_states::AbstractArray{DamageState})
    return map(clone_damage_state, damage_states)
end

"""
    update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta=1.0)

按 cohesive 单元的 interface_type 分组，分别调用 `update_damage`
（界面参数宿主：:PE_PCC→param.PCC、:NE_NCC→param.NCC）。
当所有单元属于同一界面时，退化为单次 `update_damage` 调用。
"""
function update_damage_per_interface(czm_mesh::CohesiveMesh, damage_states::AbstractVector{DamageState}, separations::Vector{Tuple{Float64, Float64}}, param::Params, czm_model::String; visc_beta::Float64=1.0)
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
    same_csc_contents(A, B) -> Bool

仅当两个 CSC 矩阵的尺寸、稀疏结构和存储数值完全一致时返回 `true`。
用于判定 BC 后刚度矩阵能否安全复用已有分解。
"""
function same_csc_contents(
    A::SparseMatrixCSC{Float64, Int64},
    B::SparseMatrixCSC{Float64, Int64},
)
    return size(A) == size(B) &&
           A.colptr == B.colptr &&
           A.rowval == B.rowval &&
           isequal(A.nzval, B.nzval)
end

"""
    solve_czm_linear_system_cached!(ws, K_bc, R_bc) -> Δu

非几何 basic/arc_length 路径的线性求解缓存。只有 BC 后矩阵内容完全一致时复用
`factorize(K_bc)`；否则重新分解，并在分解与回代均成功后提交新缓存。
"""
function solve_czm_linear_system_cached!(
    ws::CZMAssemblyWorkspace,
    K_bc::SparseMatrixCSC{Float64, Int64},
    R_bc::Vector{Float64},
)
    cached_matrix = ws.K_bc_factor_matrix
    cached_factorization = ws.K_bc_factorization
    if cached_matrix !== nothing && cached_factorization !== nothing &&
       same_csc_contents(cached_matrix, K_bc)
        return cached_factorization \ R_bc
    end

    new_factorization = factorize(K_bc)
    Δu = new_factorization \ R_bc
    ws.K_bc_factor_matrix = copy(K_bc)
    ws.K_bc_factorization = new_factorization
    return Δu
end

"""
    solve_equilibrated(K_bc, R_bc)

geo 路径线性求解的对角均衡化：D = diag(1/√colmax)，求解 (D·K·D)·y = D·R 后回代
Δu = D·y。装配切线在归一化体系下量级跨约 22 个数量级，直接分解时 UMFPACK 数值
主元退化（填充 ~69×，单次 ~2.5 s）；均衡化后填充 ~10×、~0.2 s（实测 12–14×）。
与 `K_bc \\ R_bc` 数学等价，仅舍入路径不同（解偏差 ~1e-9 相对）——geo 路径基线
自 2026-09-09 批次起按此数值口径重冻结。
"""
function solve_equilibrated(K_bc::SparseMatrixCSC{Float64, Int64}, R_bc::Vector{Float64})
    n = length(R_bc)
    d = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        vmax = 0.0
        for p in K_bc.colptr[j]:(K_bc.colptr[j+1]-1)
            a = abs(K_bc.nzval[p])
            a > vmax && (vmax = a)
        end
        d[j] = 1.0 / sqrt(vmax)
    end
    nzval = similar(K_bc.nzval)
    @inbounds for j in 1:n
        dj = d[j]
        for p in K_bc.colptr[j]:(K_bc.colptr[j+1]-1)
            nzval[p] = d[K_bc.rowval[p]] * K_bc.nzval[p] * dj
        end
    end
    Ks = SparseMatrixCSC(K_bc.m, K_bc.n, K_bc.colptr, K_bc.rowval, nzval)
    return d .* (Ks \ (d .* R_bc))
end

"""
    backtrack_line_search!(u, Δu, czm_mesh, param, damage_states, F_ext, F_thermo_chem, R_norm_current, bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws; max_halvings=8)

回溯线搜索（零化式 BC 残差）。仅用于 solve_czm_basic_step。
返回 (u_new, R_new_norm, accepted, α_used, assembly)。
accepted 时返回的 u_new 已含 BC 赋值，外部无需再执行 u = u + α*Δu；
未 accepted 时返回原始 u（未修改），外部应 break。
accepted 时 assembly = 被接受试探点的 (K_total, f_int, separations, tractions)。
接受试探点是本函数执行的最后一次装配，其输入（u_trial 与冻结的损伤/塑性试探态/
本征应变/预应力/粘性系数/BC）与下一轮 Newton 入口装配完全一致，故结果可被入口
逐位复用；未 accepted 时为 nothing。trial 塑性缓冲停留在接受点，同样与入口重写
逐位一致。
"""
function backtrack_line_search!(u::Vector{Float64}, Δu::Vector{Float64},czm_mesh::CohesiveMesh, param::Params, damage_states,F_ext::Vector{Float64}, F_thermo_chem::Vector{Float64},R_norm_current::Float64,bc_dofs::Vector{Int64}, bc_vals::Vector{Float64},K_bulk_cached, geom_cache, ws;max_halvings::Int=8, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, committed_plastic_states=nothing, trial_plastic_states=nothing, prestress=nothing, gp_damage_states=nothing, gp_trial_states=nothing)
    α = 1.0
    for _ in 1:max_halvings
        u_trial = u + α * Δu
        apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

        K_trial, f_int_trial, sep_trial, trac_trial = assemble_coupled_system(czm_mesh, u_trial, param;damage_states=damage_states, gp_damage_states=gp_damage_states, gp_trial_states=gp_trial_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain,
        plasticity=plasticity, committed_plastic_states=committed_plastic_states,
        trial_plastic_states=trial_plastic_states, prestress=prestress)

        R_trial = F_ext + F_thermo_chem - f_int_trial
        for (dof, val) in zip(bc_dofs, bc_vals)
            R_trial[dof] = val - u_trial[dof]
        end

        R_trial_norm = norm(R_trial)
        if !isnan(R_trial_norm) && R_trial_norm < R_norm_current
            return u_trial, R_trial_norm, true, α,
                (K_trial, f_int_trial, sep_trial, trac_trial)
        end

        α *= 0.5
    end
    return u, R_norm_current, false, 0.0, nothing
end

function gp_history_for_step(ms::MechState, czm_mesh::CohesiveMesh)
    weights, _ = NCweight(czm_mesh.bulk_mesh.gs.order)
    shape = (czm_mesh.n_cohesive, length(weights))
    if ms.gp_damage_states === nothing
        any(s -> s.D != 0.0 || s.D_visc != 0.0 || s.δ_max_eff != 0.0,
            ms.damage_states) && throw(ArgumentError(
            "cannot reconstruct GP damage history from a nonzero element-average history"))
        return [clone_damage_state(ms.damage_states[i]) for i in 1:shape[1], _ in 1:shape[2]]
    end
    size(ms.gp_damage_states) == shape || throw(DimensionMismatch(
        "GP damage history size $(size(ms.gp_damage_states)) must equal $shape"))
    return clone_damage_states(ms.gp_damage_states)
end

function aggregate_gp_damage_states(gp_states::Matrix{DamageState},
                                    weights::AbstractVector{<:Real})
    size(gp_states, 2) == length(weights) || throw(DimensionMismatch(
        "GP history and quadrature weights differ"))
    weight_sum = sum(weights)
    weight_sum > 0.0 || throw(ArgumentError("cohesive quadrature weights must sum to positive value"))
    element_states = Vector{DamageState}(undef, size(gp_states, 1))
    for i in axes(gp_states, 1)
        state = DamageState()
        state.D = sum(weights[g] * gp_states[i, g].D for g in eachindex(weights)) / weight_sum
        state.D_visc = sum(weights[g] * gp_states[i, g].D_visc for g in eachindex(weights)) / weight_sum
        state.δ_max_n = maximum(gp_states[i, g].δ_max_n for g in eachindex(weights))
        state.δ_max_t = maximum(gp_states[i, g].δ_max_t for g in eachindex(weights))
        state.δ_max_eff = maximum(gp_states[i, g].δ_max_eff for g in eachindex(weights))
        state.fractured = all(gp_states[i, g].fractured for g in eachindex(weights))
        state.accumulated_damage = sum(weights[g] * gp_states[i, g].accumulated_damage
                                       for g in eachindex(weights)) / weight_sum
        element_states[i] = state
    end
    return element_states
end

function prepare_plastic_state_buffers(ms::MechState, czm_mesh::CohesiveMesh,
                                       plasticity::Bool)
    plasticity || return nothing, nothing
    ne = size(czm_mesh.bulk_element, 1)
    committed = ms.plastic_states === nothing ?
        [PlasticState() for _ in 1:ne, _ in 1:4] : ms.plastic_states
    size(committed) == (ne, 4) || throw(DimensionMismatch(
        "plastic state size $(size(committed)) must equal ($ne, 4)"))
    return committed, clone_plastic_states(committed)
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

function fill_czm_result!(result::CZMResult, u::Vector{Float64}, damage_states::AbstractVector{DamageState}, separations::Vector{Tuple{Float64, Float64}}, tractions::Vector{Tuple{Float64, Float64}})
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

function solve_czm_basic_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param, ms::MechState; dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, prestress=nothing, gp_history::Bool=false)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(ms.u_prev)
        u_start = copy(u)
        damage_states = clone_damage_states(ms.damage_states)
        gp_states = gp_history ? gp_history_for_step(ms, czm_mesh) : nothing
        gp_trial_states = gp_history ? similar(gp_states) : nothing
        gp_weights = gp_history ? NCweight(czm_mesh.bulk_mesh.gs.order)[1] : nothing
        committed_plastic_states, trial_plastic_states =
            prepare_plastic_state_buffers(ms, czm_mesh, plasticity)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)

        # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
        if geo_nl
            F_thermo_chem = zeros(Float64, ndof)
            K_bulk_cached = nothing
        else
            F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
            K_bulk_cached = bulk_stiffness(czm_mesh, param)
        end
        eig_kwargs = geo_nl ? (
            geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity,
            committed_plastic_states=committed_plastic_states,
            trial_plastic_states=trial_plastic_states,
            prestress=prestress) : ()
        geom_cache = cohesive_geometry(czm_mesh)
        ws_basic = assembly_workspace(czm_mesh)

        total_iter = 0
        R_norm = Inf
        failure_R_norm = Inf
        converged = false
        converged_R_norm = Inf
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        cached_K = nothing
        cached_f_int = nothing
        cached_separations = nothing
        cached_tractions = nothing

        # One final assembly after the last Newton correction is needed to
        # judge that correction; it does not grant an extra Newton step.
        for iter in 1:(max_iter + 1)
            iter <= max_iter && (total_iter += 1)

            # 线搜索接受的试探装配在相同输入下与入口装配逐位一致，可直接复用；
            # 缓存生命周期内唯一的 ws 写者是线性求解的因子缓存字段（与装配
            # 输出别名 disjoint），任何 break 路径都不再消费缓存
            if cached_K === nothing
                K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, gp_damage_states=gp_states, gp_trial_states=gp_trial_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)
            else
                K_total = cached_K
                f_int_total = cached_f_int
                separations = cached_separations
                tractions = cached_tractions
                cached_K = nothing
                cached_f_int = nothing
                cached_separations = nothing
                cached_tractions = nothing
            end

            R = F_ext + F_thermo_chem - f_int_total

            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            failure_R_norm = R_norm
            if R_norm < tol
                trial_gp_states = gp_history ? clone_damage_states(gp_trial_states) : nothing
                trial_damage_states = gp_history ?
                    aggregate_gp_damage_states(trial_gp_states, gp_weights) :
                    update_damage_per_interface(czm_mesh, damage_states, separations,
                                                param, czm_model; visc_beta=visc_beta)
                _, f_int_committed, sep_committed, trac_committed =
                    assemble_coupled_system(czm_mesh, u, param;
                        damage_states=trial_damage_states,
                        gp_damage_states=trial_gp_states,
                        K_bulk_cached=K_bulk_cached, geom_cache=geom_cache,
                        ws=ws_basic, visc_beta=gp_history ? 0.0 : visc_beta,
                        czm_model=czm_model,
                        eig_kwargs...)
                R_committed = F_ext + F_thermo_chem - f_int_committed
                for (dof, val) in zip(bc_dofs, bc_vals)
                    R_committed[dof] = val - u[dof]
                end
                committed_R_norm = norm(R_committed)
                failure_R_norm = committed_R_norm
                if committed_R_norm < tol
                    damage_states = trial_damage_states
                    gp_history && (gp_states = trial_gp_states)
                    converged = true
                    converged_R_norm = committed_R_norm
                    separations, tractions = sep_committed, trac_committed
                    break
                end
                gp_history && break
                # Damage changed the residual. Re-equilibrate at the same load
                # with the trial history; ms is still untouched until success.
                damage_states = trial_damage_states
                continue
            end

            iter > max_iter && break
            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            Δu = try
                if geo_nl
                    solve_equilibrated(K_bc, R_bc)
                else
                    solve_czm_linear_system_cached!(ws_basic, K_bc, R_bc)
                end
            catch
                break
            end

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            u, R_norm, ls_accepted, α_used, ls_assembly = backtrack_line_search!(u, Δu, czm_mesh, param,damage_states, F_ext, F_thermo_chem, R_norm,bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws_basic;visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain, plasticity=plasticity, committed_plastic_states=committed_plastic_states, trial_plastic_states=trial_plastic_states, prestress=prestress, gp_damage_states=gp_states, gp_trial_states=gp_trial_states)
            failure_R_norm = R_norm

            if !ls_accepted
                break
            end
            cached_K, cached_f_int, cached_separations, cached_tractions = ls_assembly
        end

        if !converged
            u = u_start
            damage_states = clone_damage_states(ms.damage_states)   # 未收敛不触碰 ms（试探态丢弃）
        end

        if converged
            _, _, separations, tractions = assemble_czm_system(
                czm_mesh, u, param; damage_states=damage_states,
                gp_damage_states=gp_states,
                geom_cache=geom_cache, ws=ws_basic,
                visc_beta=gp_history ? 0.0 : visc_beta,
                czm_model=czm_model)
        else
            _, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, param; damage_states=damage_states,
                gp_damage_states=gp_states,
                K_bulk_cached=K_bulk_cached, geom_cache=geom_cache, ws=ws_basic,
                visc_beta=gp_history ? 0.0 : visc_beta,
                czm_model=czm_model, eig_kwargs...)
            R = F_ext + F_thermo_chem - f_int_total
            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end
            R_norm = norm(R)
        end

        # 失败时位移/历史回滚，但报告最后试算残差；回滚态可能恰好
        # 在旧损伤下平衡，其小残差不能解释本次求解为何失败。
        final_R_norm = converged ? converged_R_norm : failure_R_norm

        result.converged = converged
        result.iterations = total_iter
        result.residual_norm = final_R_norm
        result.displacement = u
        fill_czm_result!(result, u, damage_states, separations, tractions)
        if converged
            ms.damage_states = damage_states   # 收敛提交（D-提交语义）
            ms.gp_damage_states = gp_history ? gp_states : nothing
            ms.u_prev = copy(u)
            plasticity && (ms.plastic_states = trial_plastic_states)
        end
        return result
    end

    function cylindrical_arc_predictor(tangent::Vector{Float64}, load_start::Float64,
            step_size::Float64, previous_u_tangent::Union{Nothing,Vector{Float64}},
            previous_lambda_tangent::Float64, arc_length_alpha::Float64)
        delta_lambda = step_size
        if previous_u_tangent !== nothing &&
           dot(tangent, previous_u_tangent) +
           arc_length_alpha^2 * previous_lambda_tangent < 0.0
            delta_lambda = -delta_lambda
        end
        delta_u = tangent * delta_lambda
        return (delta_lambda=delta_lambda, delta_u=delta_u,
                lambda=load_start + delta_lambda, radius_sq=dot(delta_u, delta_u))
    end

    function cylindrical_arc_converged(R_norm::Float64, delta_u_sq::Float64,
            radius_sq::Float64, tol::Float64)
        isfinite(radius_sq) && radius_sq > 0.0 || return false
        return isfinite(R_norm) && R_norm < 10.0 * tol &&
               abs(delta_u_sq - radius_sq) / radius_sq < tol
    end

    function arc_damage_commit_check(czm_mesh::CohesiveMesh, u::Vector{Float64},
            param::Params, damage_start::Vector{DamageState},
            separations::Vector{Tuple{Float64,Float64}}, F_applied::Vector{Float64},
            bc_dofs::Vector{Int64}, bc_vals::Vector{Float64};
            K_bulk_cached=nothing, geom_cache=nothing, ws=nothing,
            visc_beta::Float64=1.0, czm_model::String="model1")
        trial_damage = update_damage_per_interface(czm_mesh, damage_start,
            separations, param, czm_model; visc_beta=visc_beta)
        K_total, f_int, sep_after, trac_after = assemble_coupled_system(
            czm_mesh, u, param; damage_states=trial_damage,
            K_bulk_cached=K_bulk_cached, geom_cache=geom_cache, ws=ws,
            visc_beta=visc_beta, czm_model=czm_model)
        R = F_applied - f_int
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        return (damage_states=trial_damage, K_total=K_total, R=R,
            residual_norm=norm(R), separations=sep_after, tractions=trac_after)
    end

    # Smooth-branch derivative of the element-average trial state. At a branch
    # boundary the local semismooth choice is zero; the residual remains the
    # exact constitutive evaluation and the arc substep can be reduced.
    function arc_trial_damage_gradient(δn::Float64, δt::Float64,
            previous::DamageState, trial::DamageState, ip::CurrentCollector,
            czm_model::String, visc_beta::Float64)
        zero_grad = (0.0, 0.0)
        (previous.fractured || trial.fractured) && return zero_grad, zero_grad
        δn_pos = max(δn, 0.0)
        if czm_model == "model1"
            e = δn_pos
            a, c = ip.δ_0, ip.δ_c
            e_n = if δn > 0.0
                1.0
            else
                0.0
            end
            e_t = 0.0
            a_n = a_t = c_n = c_t = 0.0
        else
            e = hypot(δn_pos, δt)
            e > 1e-15 || return zero_grad, zero_grad
            β = abs(δt) / e
            βη = β^ip.eta
            a = sqrt(ip.δ_0^2 + (ip.δ_0_t^2 - ip.δ_0^2) * βη)
            c = sqrt(ip.δ_c^2 + (ip.δ_c_t^2 - ip.δ_c^2) * βη)
            e_n = if δn > 0.0
                δn_pos / e
            else
                0.0
            end
            e_t = δt / e
            a_n = a_t = c_n = c_t = 0.0
            if δt != 0.0
                β_n = -β * e_n / e
                β_t = sign(δt) / e - β * e_t / e
                a_β = (ip.δ_0_t^2 - ip.δ_0^2) * ip.eta * β^(ip.eta - 1.0) / (2.0a)
                c_β = (ip.δ_c_t^2 - ip.δ_c^2) * ip.eta * β^(ip.eta - 1.0) / (2.0c)
                a_n, a_t = a_β * β_n, a_β * β_t
                c_n, c_t = c_β * β_n, c_β * β_t
            end
        end
        (e > previous.δ_max_eff && a < e < c) || return zero_grad, zero_grad

        den = c - a
        D_e = c * a / (e^2 * den)
        D_a = c * (e - c) / (e * den^2)
        D_c = -a * (e - a) / (e * den^2)
        q_eq = (D_e * e_n + D_a * a_n + D_c * c_n,
                D_e * e_t + D_a * a_t + D_c * c_t)
        raw_visc = previous.D_visc + visc_beta * (trial.D - previous.D_visc)
        q_visc = if raw_visc > previous.D_visc
            (visc_beta * q_eq[1], visc_beta * q_eq[2])
        else
            zero_grad
        end
        return q_eq, q_visc
    end

    function arc_local_separation_operator(R::Matrix{Float64}, ξ::Float64)
        N1, N2 = 0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)
        B = zeros(Float64, 2, 8)
        B[1, 1] = -N1; B[2, 2] = -N1
        B[1, 3] = -N2; B[2, 4] = -N2
        B[1, 5] = N2;  B[2, 6] = N2
        B[1, 7] = N1;  B[2, 8] = N1
        return R * B
    end

    # The first arc trial assembly was used only for these kinematic averages.
    # Compute them directly, then assemble force and tangent once with the
    # resulting element-average trial history.
    function arc_element_separations(czm_mesh::CohesiveMesh,
            u::Vector{Float64}, param::Params,
            geom_cache::Vector{CohesiveElementGeom})
        Λ = param.scale.L / param.scale.δ_czm
        separations = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
        for i in eachindex(separations)
            geom = geom_cache[i]
            dofs = geom.dofs
            R = geom.R
            δn_sum = 0.0
            δt_sum = 0.0
            w_sum = 0.0
            for (ξ, w) in zip(geom.gauss_pts, geom.gauss_wts)
                N1, N2 = 0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)
                dx = -N1 * u[dofs[1]] - N2 * u[dofs[3]] +
                     N2 * u[dofs[5]] + N1 * u[dofs[7]]
                dy = -N1 * u[dofs[2]] - N2 * u[dofs[4]] +
                     N2 * u[dofs[6]] + N1 * u[dofs[8]]
                δn_sum += w * Λ * (R[1, 1] * dx + R[1, 2] * dy)
                δt_sum += w * Λ * (R[2, 1] * dx + R[2, 2] * dy)
                w_sum += w
            end
            separations[i] = (δn_sum / w_sum, δt_sum / w_sum)
        end
        return separations
    end

    # The cohesive GP tangent holds the element-average trial history fixed.
    # Add (∂T_gp/∂D_trial)(dD_trial/dδ_avg) Λ B_avg to that local tangent.
    function arc_damage_chain_tangent!(K_total::SparseMatrixCSC{Float64,Int64},
            czm_mesh::CohesiveMesh, u::Vector{Float64}, param::Params,
            damage_start::Vector{DamageState}, trial_damage::Vector{DamageState},
            separations::Vector{Tuple{Float64,Float64}},
            geom_cache::Vector{CohesiveElementGeom}, visc_beta::Float64,
            czm_model::String)
        Λ = param.scale.L / param.scale.δ_czm
        for i in eachindex(trial_damage)
            geom = geom_cache[i]
            ip = collector_params(param, czm_mesh.cohesive_elements[i].interface_type)
            δavg_n, δavg_t = separations[i]
            q_eq, q_visc = arc_trial_damage_gradient(δavg_n, δavg_t,
                damage_start[i], trial_damage[i], ip, czm_model, visc_beta)
            all(iszero, q_eq) && all(iszero, q_visc) && continue

            w_sum = sum(geom.gauss_wts)
            ξ_avg = sum(w * ξ for (ξ, w) in zip(geom.gauss_pts, geom.gauss_wts)) / w_sum
            B_avg = arc_local_separation_operator(geom.R, ξ_avg)
            u_e = u[geom.dofs]
            for (ξ, w) in zip(geom.gauss_pts, geom.gauss_wts)
                B_gp = arc_local_separation_operator(geom.R, ξ)
                δ_gp = Λ * B_gp * u_e
                δn, δt = δ_gp
                gp_state = last(bilinear_traction_state(δn, δt,
                    trial_damage[i], ip, czm_model; visc_beta=visc_beta))
                e_gp = if czm_model == "model1"
                    max(δn, 0.0)
                else
                    hypot(max(δn, 0.0), δt)
                end
                q_gp_eq = if e_gp > trial_damage[i].δ_max_eff
                    (0.0, 0.0)
                else
                    q_eq
                end
                raw_gp_visc = trial_damage[i].D_visc +
                    visc_beta * (gp_state.D - trial_damage[i].D_visc)
                q_gp_visc = if raw_gp_visc > trial_damage[i].D_visc
                    ((1.0 - visc_beta) * q_visc[1] + visc_beta * q_gp_eq[1],
                     (1.0 - visc_beta) * q_visc[2] + visc_beta * q_gp_eq[2])
                else
                    q_visc
                end
                normal_scale = if δn >= 0.0
                    -ip.K_n * δn
                else
                    0.0
                end
                shear_scale = if czm_model == "model1"
                    0.0
                else
                    -ip.K_t * δt
                end
                wJΛ = w * geom.length * 0.5 * Λ
                for a in 1:8, b in 1:8
                    force_n = (B_gp[1, a] * normal_scale + B_gp[2, a] * shear_scale) * q_gp_visc[1]
                    force_t = (B_gp[1, a] * normal_scale + B_gp[2, a] * shear_scale) * q_gp_visc[2]
                    K_total[geom.dofs[a], geom.dofs[b]] +=
                        wJΛ * (force_n * B_avg[1, b] + force_t * B_avg[2, b])
                end
            end
        end
        return K_total
    end

    # Each Newton candidate derives its damage from the last accepted substep.
    # A rejected candidate must not become the history for the next candidate.
    function arc_equilibrium_trial(czm_mesh::CohesiveMesh, u::Vector{Float64},
            param::Params, damage_start::Vector{DamageState},
            F_applied::Vector{Float64}, bc_dofs::Vector{Int64},
            bc_vals::Vector{Float64}; K_bulk_cached=nothing,
            geom_cache=nothing, ws=nothing, visc_beta::Float64=1.0,
            czm_model::String="model1")
        geometry = geom_cache === nothing ? cohesive_geometry(czm_mesh) : geom_cache
        separations = arc_element_separations(czm_mesh, u, param, geometry)
        trial = arc_damage_commit_check(czm_mesh, u, param, damage_start,
            separations, F_applied, bc_dofs, bc_vals;
            K_bulk_cached=K_bulk_cached, geom_cache=geometry, ws=ws,
            visc_beta=visc_beta, czm_model=czm_model)
        arc_damage_chain_tangent!(trial.K_total, czm_mesh, u, param,
            damage_start, trial.damage_states, separations, geometry,
            visc_beta, czm_model)
        return trial
    end

    # Fixed-load endpoint of the nongeometric arc path. Every line-search
    # candidate uses the same accepted history, including after a rejected
    # candidate; only a balanced candidate advances MechState.
    function solve_czm_arc_target_step(czm_mesh::CohesiveMesh,
            F_ext::Vector{Float64}, param::Params, ms::MechState;
            dT_elem::Union{Vector{Float64},Nothing}=nothing,
            Δsoc_n_elem::Union{Vector{Float64},Nothing}=nothing,
            Δsoc_p_elem::Union{Vector{Float64},Nothing}=nothing,
            max_iter::Int=50, tol::Float64=1e-8, visc_beta::Float64=1.0,
            czm_model::String="model1", fix_inner::Bool=true)
        u_start = copy(ms.u_prev)
        u = copy(u_start)
        damage_start = clone_damage_states(ms.damage_states)
        result = CZMResult(length(u), czm_mesh.n_cohesive)
        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)
        F_target = F_ext + assemble_thermal_chemical_load(czm_mesh, param,
            dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = bulk_stiffness(czm_mesh, param)
        geom_cache = cohesive_geometry(czm_mesh)
        ws = assembly_workspace(czm_mesh)
        last_residual = Inf
        corrections = 0

        for iter in 0:max_iter
            trial = arc_equilibrium_trial(czm_mesh, u, param, damage_start,
                F_target, bc_dofs, bc_vals; K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws, visc_beta=visc_beta,
                czm_model=czm_model)
            last_residual = trial.residual_norm
            if last_residual < tol
                ms.u_prev = copy(u)
                ms.damage_states = trial.damage_states
                result.converged = true
                result.iterations = corrections
                result.residual_norm = last_residual
                fill_czm_result!(result, u, trial.damage_states,
                    trial.separations, trial.tractions)
                return result
            end
            iter == max_iter && break

            K_bc, R_bc = apply_bc_czm(trial.K_total, trial.R;
                bc_dofs=bc_dofs, bc_vals=bc_vals)
            Δu = try
                solve_czm_linear_system_cached!(ws, K_bc, R_bc)
            catch
                break
            end
            all(isfinite, Δu) || break
            corrections += 1

            accepted = false
            α = 1.0
            for _ in 1:11
                candidate = u + α * Δu
                apply_czm_dirichlet!(candidate, bc_dofs, bc_vals)
                candidate_trial = arc_equilibrium_trial(czm_mesh, candidate,
                    param, damage_start, F_target, bc_dofs, bc_vals;
                    K_bulk_cached=K_bulk_cached, geom_cache=geom_cache,
                    ws=ws, visc_beta=visc_beta, czm_model=czm_model)
                if isfinite(candidate_trial.residual_norm) &&
                   candidate_trial.residual_norm < (1.0 - 1e-4 * α) * last_residual
                    u = candidate
                    accepted = true
                    break
                end
                α *= 0.5
            end
            accepted || break
        end

        _, f_int_start, separations, tractions = assemble_coupled_system(czm_mesh,
            u_start, param; damage_states=damage_start,
            K_bulk_cached=K_bulk_cached, geom_cache=geom_cache, ws=ws,
            visc_beta=visc_beta, czm_model=czm_model)
        R_start = F_target - f_int_start
        for (dof, val) in zip(bc_dofs, bc_vals)
            R_start[dof] = val - u_start[dof]
        end
        result.converged = false
        result.iterations = corrections
        result.residual_norm = norm(R_start)
        fill_czm_result!(result, u_start, damage_start, separations, tractions)
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
        F_load_bc = copy(F_delta)
        zero_czm_bc_entries!(F_load_bc, bc_dofs)
        if all(iszero, F_load_bc)
            return solve_czm_arc_target_step(czm_mesh, F_ext, param, ms;
                dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem,
                Δsoc_p_elem=Δsoc_p_elem, max_iter=max_iter, tol=tol,
                visc_beta=visc_beta, czm_model=czm_model,
                fix_inner=fix_inner)
        end

        total_iter = 0
        load_progress = 0.0
        load_step = 0
        max_substep_attempts = max(100, 20 * n_load_steps)
        # A failed physical step must not consume an unbounded sequence of
        # expensive sparse solves while repeatedly shrinking arc substeps.
        max_total_iter = max(2 * max_iter, 8 * n_load_steps)
        step_size = 1.0 / max(1, n_load_steps)
        step_size_min = step_size / 128.0
        step_size_max = step_size
        last_residual = Inf
        converged_substep = false
        previous_u_tangent = nothing
        previous_lambda_tangent = 0.0
        last_target_attempt_start = NaN
        stop_reason = :none
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        while load_progress < 1.0 - 1e-12
            if total_iter >= max_total_iter
                stop_reason = :correction_budget
                break
            end
            load_step += 1
            if load_step > max_substep_attempts
                stop_reason = :substep_budget
                break
            end
            load_start = load_progress
            target_progress = min(1.0, load_start + step_size)

            u_start = copy(u)
            damage_start = clone_damage_states(damage_states)   # 子步试探回滚（ms 不受影响）
            converged_substep = false
            last_residual = Inf

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param;damage_states=damage_states, K_bulk_cached=K_bulk_cached,geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, czm_model=czm_model)

            F_applied = f_int_ref + load_start * F_delta
            R = F_applied - f_int_total
            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end
            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            tangent = try
                solve_czm_linear_system_cached!(ws, K_bc, F_load_bc)
            catch
                nothing
            end

            if tangent === nothing || any(isnan, tangent) || any(isinf, tangent)
                stop_reason = :invalid_predictor_solve
                break
            end

            predictor = cylindrical_arc_predictor(tangent, load_start,
                target_progress - load_start, previous_u_tangent,
                previous_lambda_tangent, arc_length_alpha)
            delta_lambda_pred = predictor.delta_lambda
            delta_u_pred = predictor.delta_u
            arc_target_sq = predictor.radius_sq
            if !isfinite(arc_target_sq) || arc_target_sq <= 0.0
                stop_reason = :invalid_arc_radius
                break
            end

            u = u_start + delta_u_pred
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            load_progress = predictor.lambda

            for iter in 1:max_iter
                if total_iter >= max_total_iter
                    stop_reason = :correction_budget
                    break
                end
                total_iter += 1

                F_applied = f_int_ref + load_progress * F_delta
                trial = arc_equilibrium_trial(czm_mesh, u, param,
                    damage_start, F_applied, bc_dofs, bc_vals;
                    K_bulk_cached=K_bulk_cached, geom_cache=geom_cache,
                    ws=ws, visc_beta=visc_beta, czm_model=czm_model)
                K_total, R = trial.K_total, trial.R
                separations, tractions = trial.separations, trial.tractions

                delta_u = u - u_start
                delta_lambda = load_progress - load_start
                delta_u_sq = dot(delta_u, delta_u)
                arc_constraint = delta_u_sq - arc_target_sq
                residual_norm = sqrt(trial.residual_norm^2 + arc_constraint^2)
                last_residual = residual_norm

                if cylindrical_arc_converged(trial.residual_norm,
                        delta_u_sq, arc_target_sq, tol)
                    last_residual = trial.residual_norm
                    if trial.residual_norm < tol && load_progress <= 1.0 + 1e-12
                        damage_states = trial.damage_states
                        converged_substep = true
                        previous_u_tangent = u - u_start
                        previous_lambda_tangent = load_progress - load_start
                        step_size = min(step_size * 1.25, step_size_max)
                        break
                    end
                end

                K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

                # Crisfield cylindrical arc-length: solve two linear systems
                delta_u_R = try
                    solve_czm_linear_system_cached!(ws, K_bc, R_bc)
                catch
                    nothing
                end
                if delta_u_R === nothing || any(isnan, delta_u_R) || any(isinf, delta_u_R)
                    break
                end

                delta_u_F = try
                    solve_czm_linear_system_cached!(ws, K_bc, F_load_bc)
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
                c_q = dot(du_bar, du_bar) - arc_target_sq

                discriminant = b_q^2 - 4.0 * a_q * c_q
                if discriminant < 0.0
                    break
                end

                sqrt_disc = sqrt(discriminant)
                dl1 = (-b_q + sqrt_disc) / (2.0 * a_q)
                dl2 = (-b_q - sqrt_disc) / (2.0 * a_q)

                # Continue along the last accepted path tangent. At the first
                # substep the predictor supplies the direction instead.
                du_new_1 = du_bar + dl1 * delta_u_F
                du_new_2 = du_bar + dl2 * delta_u_F
                direction_u = previous_u_tangent === nothing ? delta_u_pred : previous_u_tangent
                direction_lambda = previous_u_tangent === nothing ?
                    delta_lambda_pred : previous_lambda_tangent
                dot1 = dot(du_new_1, direction_u) +
                       arc_length_alpha^2 * (load_progress + dl1 - load_start) * direction_lambda
                dot2 = dot(du_new_2, direction_u) +
                       arc_length_alpha^2 * (load_progress + dl2 - load_start) * direction_lambda
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
                if target_progress == 1.0 &&
                   load_start >= 1.0 - step_size_max &&
                   load_start != last_target_attempt_start
                    # The arc constraint could not reach the prescribed load
                    # from this accepted state. Try its strict fixed-load
                    # equilibrium once before retrying smaller arc radii.
                    last_target_attempt_start = load_start
                    trial_ms = MechState(czm_mesh)
                    trial_ms.u_prev = copy(u_start)
                    trial_ms.damage_states = clone_damage_states(damage_start)
                    corrected = solve_czm_arc_target_step(czm_mesh, F_ext,
                        param, trial_ms; dT_elem=dT_elem,
                        Δsoc_n_elem=Δsoc_n_elem,
                        Δsoc_p_elem=Δsoc_p_elem,
                        max_iter=min(max_iter, max_total_iter - total_iter),
                        tol=tol, visc_beta=visc_beta,
                        czm_model=czm_model, fix_inner=fix_inner)
                    total_iter += corrected.iterations
                    if corrected.converged
                        ms.damage_states = trial_ms.damage_states
                        ms.u_prev = trial_ms.u_prev
                        corrected.iterations = total_iter
                        return corrected
                    end
                end
                if total_iter >= max_total_iter
                    stop_reason = :correction_budget
                    break
                end
                step_size *= 0.5

                if step_size < step_size_min
                    stop_reason = :minimum_step_size
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

        # On failure report the residual against the requested target load;
        # the accepted continuation point may itself be perfectly balanced.
        R = F_target - f_int_total
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        R_norm = norm(R)

        result.converged = false
        result.iterations = total_iter
        result.residual_norm = R_norm
        result.displacement = u

        fill_czm_result!(result, u, damage_states, separations, tractions)
        if load_progress >= 1.0 - step_size_max &&
           load_progress <= 1.0 + 1e-12 &&
           load_progress != last_target_attempt_start
            # The arc constraint no longer applies at the prescribed final
            # load. A stalled final arc substep may still have reached the
            # target's neighborhood; only the strict fixed-load correction
            # may accept and commit that state.
            trial_ms = MechState(czm_mesh)
            trial_ms.u_prev = copy(u)
            trial_ms.damage_states = damage_states
            corrected = solve_czm_arc_target_step(czm_mesh, F_ext, param, trial_ms;
                dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem,
                Δsoc_p_elem=Δsoc_p_elem,
                max_iter=min(max_iter, max_total_iter - total_iter), tol=tol,
                visc_beta=visc_beta, czm_model=czm_model,
                fix_inner=fix_inner)
            corrected.iterations += total_iter
            if corrected.converged
                ms.damage_states = trial_ms.damage_states
                ms.u_prev = trial_ms.u_prev
            end
            return corrected
        end
        @warn "CZM arc-length path did not reach a balanced target load" reason=stop_reason load_progress=load_progress attempts=load_step corrections=total_iter residual=R_norm
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
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, param, ms::MechState; dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, prestress=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)

    u = copy(ms.u_prev)
    damage_states = clone_damage_states(ms.damage_states)
    committed_plastic_states, trial_plastic_states =
        prepare_plastic_state_buffers(ms, czm_mesh, plasticity)

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)

    # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
    if geo_nl
        F_thermo_chem_total = zeros(Float64, ndof)
        K_bulk_cached = nothing
    else
        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = bulk_stiffness(czm_mesh, param)
    end
    eig_kwargs = geo_nl ? (
        geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity,
        committed_plastic_states=committed_plastic_states,
        trial_plastic_states=trial_plastic_states,
        prestress=prestress) : ()
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
            # R1-d：geo 路径与 basic 统一走对角均衡化（装配切线量级失衡，未均衡化
            # 分解填充 ~69×；非 geo 保持原路径）。失败语义不变：奇异即抛出
            Δu = geo_nl ? solve_equilibrated(K_bc, R_bc) : (K_bc \ R_bc)

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
        plasticity && (ms.plastic_states = trial_plastic_states)
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
        plasticity::Bool=false, prestress=nothing)
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
    committed_plastic_states, trial_plastic_states =
        prepare_plastic_state_buffers(ms, czm_mesh, plasticity)
    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner)
    zero_bc_vals = zeros(Float64, length(bc_vals))
    apply_czm_dirichlet!(u, bc_dofs, bc_vals)
    ws = assembly_workspace(czm_mesh)
    geom_cache = cohesive_geometry(czm_mesh)

    eig_kw = (
        geo_nl=true, plasticity=plasticity,
        committed_plastic_states=committed_plastic_states,
        trial_plastic_states=trial_plastic_states,
        prestress=prestress)
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
        u_trial, _, accepted, _, _ = backtrack_line_search!(
            u, delta_u, czm_mesh, param, damage_states,
            F_ext, zero_external, R0_norm, bc_dofs, bc_vals,
            nothing, geom_cache, ws;
            max_halvings=12, visc_beta=visc_beta, czm_model=czm_model, geo_nl=true,
            eigenstrain=mix(0.0), plasticity=plasticity,
            committed_plastic_states=committed_plastic_states,
            trial_plastic_states=trial_plastic_states, prestress=prestress)
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
    plasticity && (ms.plastic_states = trial_plastic_states)
    return result
end

"""
    czm_viscous_beta(czm_opt, dt_seconds)

按两次 CZM 更新间实际经过的物理秒数计算后向 Euler 粘性松弛系数。
"""
function czm_viscous_beta(czm_opt::CzmOptions, dt_seconds::Union{Nothing, Real})
    czm_opt.viscous_enabled || return 1.0
    τ = czm_opt.viscous_tau
    isfinite(τ) && τ >= 0.0 || throw(ArgumentError(
        "CZM viscous_tau must be finite and nonnegative seconds, got $τ"))
    τ == 0.0 && return 1.0
    dt_seconds === nothing && throw(ArgumentError(
        "CZM viscosity requires the elapsed physical time in seconds"))
    Δt = Float64(dt_seconds)
    isfinite(Δt) && Δt > 0.0 || throw(ArgumentError(
        "CZM viscosity requires finite positive elapsed seconds, got $Δt"))
    return Δt / (τ + Δt)
end

"""
    solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; dT/Δsn/Δsp...) -> CZMResult

CZM 单步统一入口（2026-08-30 终态签名）：求解配置从 `czm_opt::CzmOptions` 展开，
演化状态在 `ms::MechState` 上收敛提交（失败/试探不触碰 ms）。
"""
function solve_czm_step(czm_mesh::CohesiveMesh, ms::MechState, param, F_ext::Vector{Float64},
        czm_opt::CzmOptions;
        dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,
        eigenstrain=nothing, prestress=nothing, dt_seconds::Union{Nothing, Real}=nothing)
    method = lowercase(czm_opt.iter_method)
    if method != "gp_basic" && ms.gp_damage_states !== nothing
        throw(ArgumentError("CZM GP damage history cannot be represented by the element-average solver; keep gp_basic or reset damage history explicitly"))
    end
    if czm_opt.viscous_enabled && czm_opt.viscous_tau > 0.0 && method != "gp_basic"
        throw(ArgumentError("Physical-time CZM viscosity requires gp_basic so damage relaxes once per accepted time step"))
    end
    max_iter = czm_opt.max_iter
    tol = czm_opt.tol
    n_load_steps = czm_opt.load_steps
    arc_length_alpha = czm_opt.arc_length_alpha
    visc_beta = czm_viscous_beta(czm_opt, dt_seconds)
    czm_model = czm_opt.model
    fix_inner = czm_opt.fix_inner
    geo_nl = czm_opt.geo_nonlinear

    if method == "load_substep"
        return newton_raphson_czm(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, prestress=prestress)
    elseif method == "basic" || method == "gp_basic"
        return solve_czm_basic_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, geo_nl=geo_nl, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, prestress=prestress,
            gp_history=method == "gp_basic")
    elseif (method == "arc_length" || method == "arclength" || method == "arc-length") && geo_nl
        return solve_czm_arc_geo_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
            arc_length_alpha=arc_length_alpha,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, eigenstrain=eigenstrain,
            plasticity=czm_opt.j2_plasticity, prestress=prestress)
    elseif method == "arc_length" || method == "arclength" || method == "arc-length"
        return solve_czm_arc_length_step(czm_mesh, F_ext, param, ms;
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha,
            visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner)
    else
        error("Unknown CZM iteration method: $(czm_opt.iter_method). Use 'basic', 'gp_basic', 'load_substep', or 'arc_length'.")
    end
end
