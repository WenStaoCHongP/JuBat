# CouplingState.jl — 类型安全的状态布局与网格几何定义
# 替代 Dict{String,Any} 的 multi_spme_layout

"""
    MultiSPMeLayout

多SPMe状态向量的布局索引。初始化后不可变。
"""
struct MultiSPMeLayout
    ne::Int                        # 热单元数
    n_chem::Int                    # 每单元电化学 DOF 数 = Nrn + Nrp + Nel
    nT::Int                        # 热节点 DOF 数
    n_total::Int                   # 全局状态向量总长 = ne*n_chem + nT
    chem_range::UnitRange{Int}     # 1:(ne*n_chem)
    thermal_range::UnitRange{Int}  # (ne*n_chem+1):(ne*n_chem+nT)
    areas::Vector{Float64}         # 预计算的单元面积（网格不变量）
end

"""便捷构造器：自动计算 range 和 n_total，areas 延迟填充"""
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)
    MultiSPMeLayout(
        ne, n_chem, nT,
        ne * n_chem + nT,
        1:(ne * n_chem),
        (ne * n_chem + 1):(ne * n_chem + nT),
        zeros(Float64, ne)
    )
end

"""便捷构造器：接收 mesh 计算单元面积"""
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int, mesh_th)
    thermal_range = (ne * n_chem + 1):(ne * n_chem + nT)
    areas = zeros(Float64, ne)
    ngs = length(mesh_th.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh_th.gs.ele[g]
        areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    return MultiSPMeLayout(ne, n_chem, nT, ne * n_chem + nT,
                           1:(ne * n_chem), thermal_range, areas)
end

"""
    BoundaryEdgeCache

预计算的外边界边列表（网格不变量），用于对流边界条件装配。
"""
struct BoundaryEdgeCache
    edges::Vector{Tuple{Int,Int}}   # (node_a, node_b) 对，a < b
    L_edge::Vector{Float64}         # 边长（无量纲）
end

"""
    compute_boundary_edge_cache(mesh, is_outer)

从网格和外部节点标记中提取去重的外边界边列表。
返回 BoundaryEdgeCache（边对 + 边长），仅计算一次。
"""
function compute_boundary_edge_cache(mesh, is_outer)
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    edges = Tuple{Int,Int}[]
    L_edge = Float64[]

    for e in 1:ne
        nodes = mesh.element[e, :]
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]),
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            (is_outer[a] && is_outer[b]) || continue
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            push!(edges, key)
            push!(L_edge, hypot(x[b] - x[a], y[b] - y[a]))
        end
    end
    return BoundaryEdgeCache(edges, L_edge)
end

"""
    MeshGeometry

Jellyroll 网格的几何拓扑信息。构建后不可变。
"""
struct MeshGeometry
    element_layer::Vector{Int}                      # 每个单元的层类型 (1=NE, 2=SP, 3=PE, 4=NCC, 5=PCC)
    is_inner_layer::Vector{Bool}                    # 是否为内层
    layer_weights::Matrix{Float64}                  # ne × 5 层面积权重 [NE, SP, PE, PCC, NCC]
    interface_pairs::Vector{Tuple{Int,Int}}         # CZM 界面配对 (top_elem, bot_elem)
    czm_element_map::Dict{Int,Vector{Int}}          # 热单元号 → CZM 单元索引向量（一对多映射）
    inner_nodes::Vector{Int}                        # 内边界节点索引
    outer_nodes::Vector{Int}                        # 外边界节点索引
    boundary_edges::Union{Nothing, BoundaryEdgeCache}  # 预计算的边界边
end

# ========================================================================
# CZM Assembly Cache & Workspace
# ========================================================================

"""
    CohesiveElementGeom

预计算的单个 cohesive 单元几何信息。构建后不变。
"""
struct CohesiveElementGeom
    length::Float64                    # 单元长度
    n_vec::Vector{Float64}             # 法向量 [2]
    t_vec::Vector{Float64}             # 切向量 [2]
    R::Matrix{Float64}                 # 旋转矩阵 [2×2]
    dofs::Vector{Int64}                # 全局 DOF 编号 [8]: [2n1-1,2n1,2n2-1,2n2,2n3-1,2n3,2n4-1,2n4]
    nodes_bottom::Vector{Int64}        # 底面节点 [n1, n2]
    nodes_top::Vector{Int64}           # 顶面节点 [n4, n3]
    gauss_wts::Vector{Float64}         # Gauss 权重
    gauss_pts::Vector{Float64}         # Gauss 点坐标
end

"""
    CZMAssemblyWorkspace

CZM 每轮 Newton 迭代复用的工作区。避免单元级临时分配。
所有中间矩阵/向量运算使用 mul! 复用预分配数组。
"""
mutable struct CZMAssemblyWorkspace
    # 单元级工作区
    u_e::Vector{Float64}              # 单元位移 [8]
    K_e::Matrix{Float64}              # 单元刚度矩阵 [8×8]
    f_int_e::Vector{Float64}          # 单元内力 [8]
    B_global::Matrix{Float64}         # 全局 B 矩阵 [2×8]
    B_local::Matrix{Float64}          # 局部 B 矩阵 [2×8]（mul! 输出）
    δ_local::Vector{Float64}          # 局部分离 [2]（mul! 输出）
    BL_dT::Matrix{Float64}            # B_local' * dT_dδ [8×2]（mul! 输出）
    BL_dT_B::Matrix{Float64}          # BL_dT * B_local [8×8]（mul! 输出）
    T_vec::Vector{Float64}            # [T_n, T_t] 牵引力向量
    BLtT::Vector{Float64}             # B_local' * T_vec [8]（mul! 输出）
    # 全局级工作区
    f_int_coh::Vector{Float64}        # 全局 cohesive 内力
    separations::Vector{Tuple{Float64, Float64}}  # 单元分离
    tractions::Vector{Tuple{Float64, Float64}}     # 单元牵引
    # 预分配稀疏矩阵（避免每轮 sparse() 重建）
    K_coh_buf::Vector{Float64}        # nonzero 值缓冲区
    K_coh::SparseMatrixCSC{Float64, Int64}  # 预分配结构的稀疏矩阵

    function CZMAssemblyWorkspace(ndof::Int, n_coh::Int)
        nnz_est = max(n_coh * 64, 1)
        # 预分配稀疏矩阵结构（空值，后续 fill! 复用）
        I_tmp = sizehint!(Int64[], nnz_est)
        J_tmp = sizehint!(Int64[], nnz_est)
        V_tmp = sizehint!(Float64[], nnz_est)
        K_coh_sp = sparse(I_tmp, J_tmp, V_tmp, max(ndof, 1), max(ndof, 1))

        new(zeros(8), zeros(8, 8), zeros(8), zeros(2, 8),
            zeros(2, 8), zeros(2), zeros(8, 2), zeros(8, 8),
            zeros(2), zeros(8),
            zeros(max(ndof, 1)),
            Vector{Tuple{Float64, Float64}}(undef, max(n_coh, 1)),
            Vector{Tuple{Float64, Float64}}(undef, max(n_coh, 1)),
            zeros(nnz_est),
            K_coh_sp)
    end
end

"""
    CZMAssemblyCache

CZM 求解器的静态/准静态缓存。当 E/ν 变化或 mesh 变更时需重建。

挂载在 `Case.czm_cache` 上，跨时间步复用。
"""
mutable struct CZMAssemblyCache
    K_bulk::SparseMatrixCSC{Float64, Int64}     # bulk 刚度矩阵
    bulk_dofs::Vector{Vector{Int64}}            # 每个 bulk 元素的 DOF 映射
    cohesive_geom::Vector{CohesiveElementGeom}  # 预计算的 cohesive 单元几何
    bc_dofs::Vector{Int64}                      # 边界条件 DOF
    bc_vals::Vector{Float64}                    # 边界条件值
    ws::CZMAssemblyWorkspace                    # 可复用工作区（跨时间步）
    E_eff::Float64                              # 缓存对应的 E_eff
    ν_eff::Float64                              # 缓存对应的 ν_eff
    valid::Bool                                 # 缓存是否有效

    function CZMAssemblyCache()
        empty_ws = CZMAssemblyWorkspace(0, 0)
        new(spzeros(0, 0), Vector{Vector{Int64}}(), CohesiveElementGeom[],
            Int64[], Float64[], empty_ws, 0.0, 0.0, false)
    end
end

"""
    CzmLayout

CZM 求解的布局信息和跨时间步状态。对标电化学的 MultiSPMeLayout。
"""
mutable struct CzmLayout
    n_coh::Int                    # cohesive 单元数
    ndof::Int                     # 总位移 DOF 数 (2 * nnode)
    u_prev::Vector{Float64}       # 上一步位移场（跨时间步持有）
end

"""便捷构造器：从 czm_mesh 初始化"""
function CzmLayout(czm_mesh::CohesiveMesh)
    n_coh = czm_mesh.n_cohesive
    ndof = 2 * czm_mesh.nnode
    CzmLayout(n_coh, ndof, zeros(Float64, ndof))
end

# ========================================================================
# CZM Coupling Helpers
# ========================================================================
# 参数计算、应变输入、损伤更新入口。
# 从 CzmSolve.jl 迁移至此，保持 CouplingState.jl 作为"状态+耦合"的统一归属点。
# 函数调用 solve_czm_step / ensure_czm_cache 等在运行时解析，无需调整 include 顺序。
# ========================================================================

"""
    compute_czm_effective_params(case)

计算 CZM 求解所需的有效材料参数，全部基于 `SetParams.NormaliseParam`
后的归一化参数。

# 返回
- `E_eff`: 有效弹性模量 [-]
- `ν_eff`: 有效泊松比 [-]
- `α_eff`: 有效热膨胀系数 [-]
- `β_n`: 负极扩散应变系数 [-]
- `β_p`: 正极扩散应变系数 [-]
"""
function compute_czm_effective_params(case)
    param = case.param

    # 有效弹性模量（厚度加权平均，均为归一化量）
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 有效泊松比（厚度加权平均）
    ν_eff = (param.NE.nu * param.NE.thickness + param.PE.nu * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 有效热膨胀系数（厚度加权平均，已按 T_ref 归一化）
    α_eff = (param.NE.alphaT * param.NE.thickness + param.PE.alphaT * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 扩散应变系数 β = (Ω * c_s,max) / 3 已在 SetParams 中完成归一化
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0

    return E_eff, ν_eff, α_eff, β_n, β_p
end

"""
    compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

计算 CZM 损伤计算所需的单元级应变输入。

# 返回
- `dT_elem`: 每个单元的温度变化 [K]
- `Δsoc_n_elem`: 每个单元的负极 SOC 变化 [-]
- `Δsoc_p_elem`: 每个单元的正极 SOC 变化 [-]
"""
function compute_czm_strain_inputs(case, variables::Dict, czm_mesh, T_nodes_carry)
    ne = size(czm_mesh.bulk_element, 1)
    param = case.param

    # 参考 SOC（归一化值）
    soc_ref_n = param.NE.cs0
    soc_ref_p = param.PE.cs0

    # 初始化输出数组
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)

    # 提取温度场（无量纲温度 T* = T / T_ref）
    if length(T_nodes_carry) >= czm_mesh.nnode
        for e in 1:ne
            nodes = czm_mesh.bulk_element[e, :]
            T_elem_nd = 0.0
            valid_nodes = 0
            for n in nodes
                if n <= length(T_nodes_carry)
                    T_elem_nd += T_nodes_carry[n]
                    valid_nodes += 1
                end
            end
            if valid_nodes > 0
                T_elem_nd /= valid_nodes
                dT_elem[e] = T_elem_nd - param.cell.T0
            end
        end
    end

    # 提取 SOC 分布（如果 variables 中有）
    soc_n_elem = variables["thermal2D element soc_n"]
    soc_p_elem = variables["thermal2D element soc_p"]

    # 处理数组维度（可能是 ne×1 或 ne×num）
    if isa(soc_n_elem, AbstractMatrix)
        soc_n_elem = soc_n_elem[:, end]
        soc_p_elem = soc_p_elem[:, end]
    end

    for e in 1:min(ne, length(soc_n_elem))
        Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
        Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
    end

    return dT_elem, Δsoc_n_elem, Δsoc_p_elem
end

"""
    update_czm_damage!(case, variables, T_nodes_carry)

更新 CZM 网格的损伤状态。

使用牛顿-拉弗森迭代求解力学平衡方程，通过载荷子步法处理软化收敛问题。
从 case.czm_mesh、case.param.cohesive、case.czm_layout 获取内部状态。

# 参数
- `case`: Case 对象（需已设置 czm_mesh 和 czm_layout）
- `variables`: 当前时间步的变量字典
- `T_nodes_carry`: 当前温度场

# 返回
- `u_czm`: 更新后的 CZM 位移场
- `converged`: 是否收敛
"""
function update_czm_damage!(case, variables, T_nodes_carry)
    czm_mesh = case.czm_mesh
    czm_params = case.param.cohesive
    param_dim = case.param_dim
    param = case.param

    # 同步CZM模型选项（model1 or mix）
    czm_params.czm_model = case.opt.czm_model

    # 计算有效材料参数
    E_eff, ν_eff, α_eff, β_n, β_p = compute_czm_effective_params(case)

    # 构建或复用 CZM 缓存
    cache = ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)

    # 计算应变输入
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

    # 外力向量（一般为零）
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)

    # 诊断：检测输入异常（帮助定位 NaN 来源）
    has_nan_T = any(isnan, T_nodes_carry)
    has_nan_soc_n = any(isnan, variables["thermal2D element soc_n"])
    has_nan_soc_p = any(isnan, variables["thermal2D element soc_p"])
    if has_nan_T || has_nan_soc_n || has_nan_soc_p || any(isnan, dT_elem) || any(isnan, Δsoc_n_elem) || any(isnan, Δsoc_p_elem)
        @warn "CZM inputs contain NaN" has_nan_T=has_nan_T has_nan_soc_n=has_nan_soc_n has_nan_soc_p=has_nan_soc_p n_nan_dT=count(isnan, dT_elem) n_nan_soc_n=count(isnan, Δsoc_n_elem) n_nan_soc_p=count(isnan, Δsoc_p_elem)
    end

    # 初始化位移（从 czm_layout 获取上一步值）
    u_czm_prev = case.czm_layout !== nothing ? case.czm_layout.u_prev : zeros(Float64, ndof)
    if length(u_czm_prev) != ndof
        u_czm_prev = zeros(Float64, ndof)
    elseif any(isnan, u_czm_prev)
        @warn "CZM u_czm_prev contains NaN, resetting to zeros"
        u_czm_prev = zeros(Float64, ndof)
    end

    # 调用 CZM 求解器（可选迭代方式）
    iter_method = case.opt.czm_iter_method
    max_iter = case.opt.czm_max_iter
    tol = case.opt.czm_tol
    n_load_steps = case.opt.czm_load_steps
    arc_length_alpha = case.opt.czm_arc_length_alpha

    result, updated_czm_mesh = solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, czm_params, param, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, iter_method=iter_method,
        cache=cache
    )

    # 诊断：检查求解结果是否异常
    has_nan_disp = any(isnan, result.displacement)
    has_nan_damage = any(ds -> isnan(ds.D), updated_czm_mesh.damage_states)
    if has_nan_disp || has_nan_damage || !result.converged
        @warn "CZM solve issue" converged=result.converged iterations=result.iterations residual=round(result.residual_norm; digits=4) has_nan_disp=has_nan_disp has_nan_damage=has_nan_damage
    end

    # Only commit damage states when the nonlinear solve converged.
    # This avoids propagating a partially converged or diverged state into the next time step.
    if result.converged
        czm_mesh.damage_states = updated_czm_mesh.damage_states
        if case.czm_layout !== nothing
            case.czm_layout.u_prev = result.displacement
        end
    end

    return result.displacement, result.converged
end

"""
    update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)

6 参数兼容入口：自动构建 CzmLayout 并委托给 3 参数版本。
"""
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
    # 确保 czm_layout 存在
    if case.czm_layout === nothing
        case.czm_layout = CzmLayout(czm_mesh)
    end
    # 同步外部传入的 u_prev
    if u_czm_prev !== nothing && length(u_czm_prev) == 2 * czm_mesh.nnode
        case.czm_layout.u_prev = u_czm_prev
    end
    return update_czm_damage!(case, variables, T_nodes_carry)
end
