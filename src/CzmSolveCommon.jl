mutable struct CZMResult
    displacement::Vector{Float64}
    damage::Vector{Float64}
    damage_visc::Vector{Float64}
    fractured::BitVector
    traction_n::Vector{Float64}
    traction_t::Vector{Float64}
    separation_n::Vector{Float64}
    separation_t::Vector{Float64}
    converged::Bool
    iterations::Int64
    residual_norm::Float64
    
    CZMResult(ndof::Int, n_coh::Int) = new(
        zeros(ndof), zeros(n_coh), zeros(n_coh), falses(n_coh), zeros(n_coh), zeros(n_coh),
        zeros(n_coh), zeros(n_coh), false, 0, Inf)
end

mechanical_eigenstrain(load::MechanicalLoad) =
    (dT=load.dT, Δsn=load.Δsoc_n, Δsp=load.Δsoc_p)

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
    extract_bc_dofs(czm_mesh, param; fix_inner=true)

从 czm_mesh 提取 Dirichlet BC 的自由度列表和对应值（每次求解入口现算，
不缓存——identify_bc_nodes_czm 为 O(nnode)，成本可忽略）。
"""
function extract_bc_dofs(czm_mesh::CohesiveMesh, param; fix_inner::Bool=true, boundary=nothing)
    boundary === nothing || return mechanical_boundary_dofs(boundary)
    bc_nodes, _, _, _ = identify_bc_nodes_czm(czm_mesh, param; fix_inner=fix_inner)
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
    backtrack_line_search!(u, Δu, czm_mesh, param, F_ext, F_thermo_chem, R_norm_current, bc_dofs, bc_vals; czm_workspace, max_halvings=8)

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
function backtrack_line_search!(u::Vector{Float64}, Δu::Vector{Float64},czm_mesh::CohesiveMesh, param::Params, F_ext::Vector{Float64}, F_thermo_chem::Vector{Float64},R_norm_current::Float64,bc_dofs::AbstractVector{Int64}, bc_vals::Vector{Float64}; czm_workspace::CZMWorkspace,max_halvings::Int=8, visc_beta::Float64=1.0, czm_model::String="model1", geo_nl::Bool=false, eigenstrain=nothing, plasticity::Bool=false, committed_plastic_states=nothing, trial_plastic_states=nothing, prestress=nothing, gp_damage_states::Matrix{DamageState}, gp_trial_states=nothing)
    α = 1.0
    for _ in 1:max_halvings
        u_trial = u + α * Δu
        apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

        K_trial, f_int_trial, sep_trial, trac_trial = assemble_coupled_system(czm_mesh, u_trial, param; czm_workspace=czm_workspace,gp_damage_states=gp_damage_states, gp_trial_states=gp_trial_states, visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain,
        plasticity=plasticity, committed_plastic_states=committed_plastic_states,
        trial_plastic_states=trial_plastic_states, prestress=prestress)

        R_trial = F_ext + F_thermo_chem - f_int_trial
        set_czm_residual!(R_trial, u_trial, bc_dofs, bc_vals)

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
    require_mechanical_load_history(ms)
    length(ms.damage_states) == length(czm_mesh.cohesive_elements) || throw(DimensionMismatch(
        "element damage history length must equal $(length(czm_mesh.cohesive_elements))"))
    weights, _ = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)
    shape = (length(czm_mesh.cohesive_elements), length(weights))
    if ms.gp_damage_states === nothing
        ms.mechanics_initialized && throw(ArgumentError("initialized mechanics lacks required GP history"))
        any(s -> s.D != 0.0 || s.D_visc != 0.0 || s.δ_max_eff != 0.0 ||
                 s.δ_max_n != 0.0 || s.δ_max_t != 0.0 ||
                 s.accumulated_damage != 0.0 || s.fractured,
            ms.damage_states) && throw(ArgumentError(
            "cannot reconstruct GP damage history from a nonzero element-average history"))
        return [DamageState() for _ in 1:shape[1], _ in 1:shape[2]]
    end
    size(ms.gp_damage_states) == shape || throw(DimensionMismatch(
        "GP damage history size $(size(ms.gp_damage_states)) must equal $shape"))
    all(s->all(f->isfinite(getfield(s,f)),fieldnames(typeof(s))) && 0<=s.D<=1 && 0<=s.D_visc<=1,ms.gp_damage_states) || throw(ArgumentError("invalid GP damage history"))
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
    committed = if ms.plastic_states === nothing
        ms.mechanics_initialized && throw(ArgumentError("initialized mechanics lacks required plastic history"))
        [PlasticState() for _ in 1:ne, _ in 1:4]
    else
        ms.plastic_states
    end
    size(committed) == (ne, 4) || throw(DimensionMismatch(
        "plastic state size $(size(committed)) must equal ($ne, 4)"))
    all(s->all(isfinite,s.eps_p) && isfinite(s.kappa),committed) || throw(ArgumentError("invalid committed plastic history"))
    return committed, clone_plastic_states(committed)
end

function apply_czm_dirichlet!(u::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64}, bc_vals::AbstractVector{Float64})
    for (dof,val) in zip(bc_dofs,bc_vals)
        u[dof] = val
    end
    return u
end

# Optional internal observer for bounded diagnostic capture. It is inert in
# production and restored by the diagnostic entry point in a finally block.
const mechanical_step_observer = Ref{Any}(nothing)

# Prepare the actual successful load before publishing any mechanical field.
# The caller supplies the strain fields actually used by the accepted assembly.
function commit_czm_state!(ms, mesh, u, damage, gp, strain::NamedTuple;
        plastic_states=ms.plastic_states)
    ne=size(mesh.bulk_element,1)
    load = (dT=copy(strain.dT),Δsn=copy(strain.Δsn),Δsp=copy(strain.Δsp))
    all(a->length(a)==ne,(load.dT,load.Δsn,load.Δsp)) || throw(DimensionMismatch("successful mechanical load shape"))
    all(a->all(isfinite,a),(load.dT,load.Δsn,load.Δsp)) || throw(ArgumentError("non-finite successful mechanical load"))
    displacement=copy(u)
    ms.u_prev=displacement
    ms.damage_states=damage
    ms.gp_damage_states=gp
    ms.plastic_states=plastic_states
    ms.last_eigenstrain=load
    ms.mechanics_initialized=true
    return ms
end

function zero_czm_bc_entries!(v::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64})
    for dof in bc_dofs
        v[dof] = 0.0
    end
    return v
end

function fill_czm_result!(result::CZMResult, u::Vector{Float64}, damage_states::AbstractVector{DamageState}, separations::Vector{Tuple{Float64, Float64}}, tractions::Vector{Tuple{Float64, Float64}})
    length(result.damage) == length(damage_states) == length(separations) == length(tractions) ||
        throw(DimensionMismatch("CZM result fields must have the same element count"))
    result.displacement = copy(u)
    for i in eachindex(damage_states)
        result.damage[i] = damage_states[i].D
        result.damage_visc[i] = damage_states[i].D_visc
        result.fractured[i] = damage_states[i].fractured
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end
    return result
end

