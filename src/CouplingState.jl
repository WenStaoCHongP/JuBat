# CouplingState.jl — 类型安全的状态布局与网格几何定义
# 替代 Dict{String,Any} 的 multi_spme_layout

using Parameters: @with_kw

# ========================================================================
# CzmInterfaceParams & CzmParamCache (spec §3.5.1 v3 / §3.5.2)
# 按界面类型分组的 CZM 本构参数缓存，供 Chunk 2/4/6 使用。
# ========================================================================

"""
    CzmInterfaceParams

单一界面类型（如 :PE_PCC 或 :NE_NCC）的 CZM 本构参数（归一化后）。
按 spec §3.5.1 v3 定义 20 个字段，覆盖 Materialmatrix.jl 实际读取的全部字段。

所有字段都已通过 NormaliseParam 归一化：
- E_eff, σ_max, τ_max: / scale.σ_czm
- δ_0_n, δ_c_n, δ_0_t, δ_c_t: / scale.δ_czm（δ_czm = 2·G_c_pe_pcc/σ_max_pe_pcc，锚定界面 δ_c* ≡ 1）
- G_c, G_c_t: / scale.G_czm
- K_n, K_t: / scale.K_czm
- η, czm_model, h_c0, k_air, lambda_m, beta, threshold: 沿用原 Cohesive 归一化（无因次或已有尺度）

派生量（宏观力学无量纲化重设计 v2，见 docs/planning-with-files/力学模块修改/宏观力学模块无量纲化重设计.md）：
- Λ = scale.L / scale.δ_czm：位移空间（L 归一）→ 分离空间（δ_czm 归一）换算因子。
  装配时 δ̃ = Λ·B·ũ，切线刚度乘一次 Λ，内力不乘（虚功一致性，设计文档 §5）。
- E_star：界面双材料等效模量（调和平均），/ scale.σ_czm，用于 L_ch 与数值判据
- L_ch：界面内禀长度 E*·G_c/σ_max²，/ scale.L，用于网格分辨率判据
"""
@with_kw struct CzmInterfaceParams
    # ---- 体模量与热化学载荷（assemble_bulk_stiffness、assemble_thermal_chemical_load 用）----
    E_eff::Float64 = 0.0       # 涂层模量（PE.E_coat 或 NE.E_coat），非全栈均一化
    ν::Float64 = 0.0           # 涂层泊松比
    α::Float64 = 0.0           # 涂层热膨胀系数（归一化）

    # ---- 无量纲化派生量（重设计 v2）----
    Λ::Float64 = 1.0           # 位移→分离换算因子 scale.L/scale.δ_czm（旧方案 δ_czm=L 时为 1）
    E_star::Float64 = 0.0      # 双材料等效模量（调和平均，/ scale.σ_czm），诊断用
    L_ch::Float64 = 0.0        # 内禀长度 E*·G_c/σ_max²（/ scale.L），诊断用

    # ---- Mode I（法向）---- bilinear_* 与 compute_gap_conductance 用
    σ_max::Float64 = 0.0       # 最大法向牵引
    K_n::Float64 = 0.0         # 法向初始刚度
    δ_0_n::Float64 = 0.0       # 法向损伤起始位移
    δ_c_n::Float64 = 0.0       # 法向临界位移
    G_c::Float64 = 0.0         # Mode I 断裂能

    # ---- Mode II（切向）---- bilinear_* 用
    τ_max::Float64 = 0.0       # 最大切向牵引
    K_t::Float64 = 0.0         # 切向初始刚度
    δ_0_t::Float64 = 0.0       # 切向损伤起始位移
    δ_c_t::Float64 = 0.0       # 切向临界位移
    G_c_t::Float64 = 0.0       # Mode II 断裂能

    # ---- BK 混合模式 + 本构选择 ----
    η::Float64 = 1.45          # BK 准则指数（来自 Cohesive.eta）
    czm_model::String = "model1"   # 本构模型标识

    # ---- 界面热阻（compute_gap_conductance 用）----
    h_c0::Float64 = 1e7        # 完全接触界面传热系数
    k_air::Float64 = 0.026     # 空气导热系数
    lambda_m::Float64 = 70e-9  # 界面微观粗糙度尺度
    beta::Float64 = 1.0        # 粗糙度指数
    threshold::Float64 = 70e-9 # 间隙阈值
end

"""
    CzmParamCache

按界面类型分组的 CZM 参数缓存（spec §3.5.2）。
- param_ref: 保留 param 引用，供 assemble_bulk_stiffness 读 PE/NE.E_coat 等
- id: 内容哈希 hash((hash(pe_pcc), hash(ne_ncc)))，用于 ensure_czm_cache 快速失效判定。
  Task 4.4 fix：原为 objectid(param)，但原位修改 param 字段不改变 objectid，导致漏检。
  改为内容哈希后，任何 CzmInterfaceParams 字段值变化都能触发失效。
"""
struct CzmParamCache
    by_interface::Dict{Symbol, CzmInterfaceParams}
    param_ref::Params
    id::UInt64
end

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
    return MultiSPMeLayout(ne, n_chem, nT, ne * n_chem + nT,1:(ne * n_chem), thermal_range, areas)
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
    boundary_edges::BoundaryEdgeCache                  # 预计算的边界边
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

CZM 求解器的静态/准静态缓存。失效判据基于 `czm_mesh_id`（=`objectid(czm_mesh)`）
与 `param_cache_id`（=`param_cache.id`）——任一变化或 `fix_inner` 切换即重建。

挂载在 `Case.czm_cache` 上，跨时间步复用。
"""
mutable struct CZMAssemblyCache
    K_bulk::SparseMatrixCSC{Float64, Int64}     # bulk 刚度矩阵
    bulk_dofs::Vector{Vector{Int64}}            # 每个 bulk 元素的 DOF 映射
    cohesive_geom::Vector{CohesiveElementGeom}  # 预计算的 cohesive 单元几何
    bc_dofs::Vector{Int64}                      # 边界条件 DOF
    bc_vals::Vector{Float64}                    # 边界条件值
    ws::CZMAssemblyWorkspace                    # 可复用工作区（跨时间步）
    fix_inner::Bool                             # 是否固定内圈节点（影响 BC 构造）
    valid::Bool                                 # 缓存是否有效
    czm_mesh_id::UInt64                         # objectid(czm_mesh)，mesh 失效判据
    param_cache_id::UInt64                      # param_cache.id，参数失效判据

    function CZMAssemblyCache()
        empty_ws = CZMAssemblyWorkspace(0, 0)
        new(spzeros(0, 0), Vector{Vector{Int64}}(), CohesiveElementGeom[],
            Int64[], Float64[], empty_ws, true, false,
            UInt64(0), UInt64(0))
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
    plastic_states::Union{Nothing, Matrix{PlasticState}}  # PCC/NCC 高斯点塑性状态（Batch 3，[ne,4]；收敛才提交 D-B3-2）
    winding_prestress::Union{Nothing, Vector{NTuple{3, Float64}}}  # 卷绕预应力 σ₀ 场（Batch 2'，逐单元全局三分量；几何固定一次计算持久持有）
end

"""便捷构造器：从 czm_mesh 初始化"""
function CzmLayout(czm_mesh::CohesiveMesh)
    n_coh = czm_mesh.n_cohesive
    ndof = 2 * czm_mesh.nnode
    CzmLayout(n_coh, ndof, zeros(Float64, ndof), nothing, nothing)
end

# ========================================================================
# CZM Coupling Helpers
# ========================================================================
# 参数计算、应变输入、损伤更新入口。
# 从 CzmSolve.jl 迁移至此，保持 CouplingState.jl 作为"状态+耦合"的统一归属点。
# 函数调用 solve_czm_step / ensure_czm_cache 等在运行时解析，无需调整 include 顺序。
# ========================================================================

"""
    compute_czm_params_per_interface(case) -> CzmParamCache

按界面类型计算 CZM 参数。E_eff 用涂层模量（PE.E_coat / NE.E_coat），
不再做全栈均一化。

返回 `CzmParamCache`，包含 `:PE_PCC` 与 `:NE_NCC` 两个条目。
"""
function compute_czm_params_per_interface(case)
    param = case.param
    scale = case.param_dim.scale

    # 入口断言（统一使用归一化 param；positivity 与尺度无关）
    @assert param.PE.E_coat > 0 && param.NE.E_coat > 0 "CZM 应变驱动需要 PE/NE.E_coat > 0"
    @assert scale.σ_czm > 0 "scale.σ_czm = 0; must populate cohesive.σ_max_pe_pcc (Task 2.1) before enabling CZM"
    @assert scale.E_coat > 0 "scale.E_coat = 0; ChooseCell 检测到 PE/NE.E_coat 缺失，宏观力学分析不可用"
    @assert param.cohesive.σ_max_pe_pcc > 0 "cohesive.σ_max_pe_pcc 必须为正"
    @assert param.cohesive.σ_max_ne_ncc > 0 "cohesive.σ_max_ne_ncc 必须为正"
    @assert param.cohesive.G_c_pe_pcc > 0 && param.cohesive.G_c_ne_ncc > 0 "G_c_pe_pcc / G_c_ne_ncc 必须为正"
    @assert param.cohesive.K_n_pe_pcc > 0 && param.cohesive.K_n_ne_ncc > 0 "K_n_pe_pcc / K_n_ne_ncc 必须为正"

    # E_eff: 用涂层模量（非全栈均一化）
    E_eff_pe = param.PE.E_coat * scale.E_coat / scale.σ_czm
    E_eff_ne = param.NE.E_coat * scale.E_coat / scale.σ_czm

    # ---- 无量纲化派生量（重设计 v2 §2.2）----
    # Λ：位移空间（L 归一）→ 分离空间（δ_czm 归一）换算因子
    @assert scale.δ_czm > 0 "scale.δ_czm = 0; ChooseCell 未能锚定 δ_czm（检查 cohesive.σ_max_pe_pcc / G_c_pe_pcc）"
    Λ = scale.L / scale.δ_czm

    # E*：界面双材料等效模量（调和平均，量纲 Pa → / σ_czm）；L_ch：内禀长度（→ / L）
    coh_dim = case.param_dim.cohesive
    E_pe_dim  = case.param_dim.PE.E_coat
    E_ne_dim  = case.param_dim.NE.E_coat
    E_pcc_dim = case.param_dim.PCC.E
    E_ncc_dim = case.param_dim.NCC.E
    E_star_pe_dim = 2 * E_pe_dim * E_pcc_dim / (E_pe_dim + E_pcc_dim)
    E_star_ne_dim = 2 * E_ne_dim * E_ncc_dim / (E_ne_dim + E_ncc_dim)
    L_ch_pe = E_star_pe_dim * coh_dim.G_c_pe_pcc / coh_dim.σ_max_pe_pcc^2 / scale.L
    L_ch_ne = E_star_ne_dim * coh_dim.G_c_ne_ncc / coh_dim.σ_max_ne_ncc^2 / scale.L

    coh = param.cohesive

    # K_0 下界判据（重设计 v2 §6）：δ_0* ≤ 0.1（即 K_0* ≥ 10），否则损伤起始点
    # 贴近临界分离，本构在归一化空间不可分辨。
    if coh.δ_0_pe_pcc > 0.1 || coh.δ_0_ne_ncc > 0.1
        @warn "CZM 初始刚度过软：δ_0* > 0.1（PE-PCC=$(round(coh.δ_0_pe_pcc; sigdigits=3)), NE-NCC=$(round(coh.δ_0_ne_ncc; sigdigits=3))）。建议 K_0 ≥ 10·σ_max/δ_c。" maxlog=1
    end

    pe_pcc = CzmInterfaceParams(
        E_eff = E_eff_pe,
        Λ = Λ,
        E_star = E_star_pe_dim / scale.σ_czm,
        L_ch = L_ch_pe,
        ν = param.PE.nu_coat,
        α = param.PE.alphaT,
        σ_max = coh.σ_max_pe_pcc,
        K_n = coh.K_n_pe_pcc,
        δ_0_n = coh.δ_0_pe_pcc,
        δ_c_n = coh.δ_c_pe_pcc,
        G_c = coh.G_c_pe_pcc,
        τ_max = coh.τ_max_pe_pcc,
        K_t = coh.K_t_pe_pcc,
        δ_0_t = coh.δ_0_pe_pcc_t,
        δ_c_t = coh.δ_c_pe_pcc_t,
        G_c_t = coh.G_c_pe_pcc_t,
        η = coh.eta,
        czm_model = coh.czm_model,
        h_c0 = coh.h_c0,
        k_air = coh.k_air,
        lambda_m = coh.lambda_m,
        beta = coh.beta,
        threshold = coh.threshold,
    )

    ne_ncc = CzmInterfaceParams(
        E_eff = E_eff_ne,
        Λ = Λ,
        E_star = E_star_ne_dim / scale.σ_czm,
        L_ch = L_ch_ne,
        ν = param.NE.nu_coat,
        α = param.NE.alphaT,
        σ_max = coh.σ_max_ne_ncc,
        K_n = coh.K_n_ne_ncc,
        δ_0_n = coh.δ_0_ne_ncc,
        δ_c_n = coh.δ_c_ne_ncc,
        G_c = coh.G_c_ne_ncc,
        τ_max = coh.τ_max_ne_ncc,
        K_t = coh.K_t_ne_ncc,
        δ_0_t = coh.δ_0_ne_ncc_t,
        δ_c_t = coh.δ_c_ne_ncc_t,
        G_c_t = coh.G_c_ne_ncc_t,
        η = coh.eta,
        czm_model = coh.czm_model,
        h_c0 = coh.h_c0,
        k_air = coh.k_air,
        lambda_m = coh.lambda_m,
        beta = coh.beta,
        threshold = coh.threshold,
    )

    # spec §3.5.2 + Task 4.4 reviewer fix：id 用内容哈希（hash(pe_pcc), hash(ne_ncc)），
    # 而非 objectid(param)。原位修改 param 字段（mutating case.param.cohesive.σ_max_*）
    # 不改变 objectid，会导致缓存失效漏检。CzmInterfaceParams 是 @with_kw 不可变 struct，
    # 字段均为 Float64 / String，Julia 内置 hash 可直接处理。
    content_hash = hash((hash(pe_pcc), hash(ne_ncc)))
    return CzmParamCache(Dict(:PE_PCC => pe_pcc, :NE_NCC => ne_ncc), param, content_hash)
end

"""
    compute_czm_strain_inputs(case, variables, T_nodes) -> NamedTuple

按 spec v2 §5.1 计算 CZM 体单元粒度的 dT、Δsoc_p、Δsoc_n，以及 CZM 节点粒度的 T_czm_nodes。

# 数据流（与 spec §5.1 表格一致）
- T_czm_nodes：T_nodes (粗热节点 K) → thermal_to_czm 矩阵 → CZM 节点温度
- dT_czm：粗热单元 dT (= avg(T_nodes[thermal_elem[e]]) - T0)
          → thermal_elem_map 直接取值（element-to-element）
- Δsoc_p/n：variables["thermal2D element soc_p/n"] (粗热单元 mol/m³)
             → thermal_elem_map 直接取值 - cs0，按 material_type 分发

# 单位契约
- T_nodes, T0, dT_czm, T_czm_nodes: Kelvin
- soc_p/n, cs0, Δsoc: mol/m³
"""
function compute_czm_strain_inputs(case, variables, T_nodes)
    czm_mesh = case.czm_mesh
    submesh = czm_mesh.czm_submesh
    param = case.param
    ne_czm = size(submesh.mesh.element, 1)

    thermal_mesh = case.mesh["thermal2D"]
    thermal_elem = thermal_mesh.element
    n_thermal_node = thermal_mesh.nlen
    n_thermal_elem = size(thermal_elem, 1)

    # ---- 在任何乘法或索引前验证粗热场 → CZM 的完整接口契约 ----
    M = czm_mesh.thermal_to_czm
    M === nothing && error("compute_czm_strain_inputs: thermal_to_czm 矩阵未构造（请在 create_czm_mesh 中调用 build_thermal_to_czm_interp）")
    T_nodes isa AbstractVector || throw(DimensionMismatch(
        "compute_czm_strain_inputs: T_nodes 必须是一维活动热节点向量，实际 size=$(size(T_nodes))"))
    if length(T_nodes) != n_thermal_node || size(M, 1) != submesh.mesh.nlen || size(M, 2) != n_thermal_node
        throw(DimensionMismatch(
            "compute_czm_strain_inputs: 热节点接口尺寸不一致：thermal_to_czm size=$(size(M))，" *
            "CZM 节点数=$(submesh.mesh.nlen)，活动热网格节点数=$(n_thermal_node)，T_nodes 长度=$(length(T_nodes))"))
    end
    all(isfinite, T_nodes) || throw(ArgumentError(
        "compute_czm_strain_inputs: 活动热节点温度包含非有限值"))

    thermal_elem_map = submesh.thermal_elem_map
    length(thermal_elem_map) == ne_czm || throw(DimensionMismatch(
        "compute_czm_strain_inputs: thermal_elem_map 长度=$(length(thermal_elem_map))，CZM 单元数=$ne_czm"))
    invalid_parent = findfirst(e_th -> !(1 <= e_th <= n_thermal_elem), thermal_elem_map)
    invalid_parent === nothing || throw(ArgumentError(
        "compute_czm_strain_inputs: thermal_elem_map[$invalid_parent]=$(thermal_elem_map[invalid_parent]) " *
        "超出活动热单元范围 1:$n_thermal_elem"))

    length(submesh.material_type) == ne_czm || throw(DimensionMismatch(
        "compute_czm_strain_inputs: material_type 长度=$(length(submesh.material_type))，CZM 单元数=$ne_czm"))
    size(thermal_elem, 2) == 4 || throw(DimensionMismatch(
        "compute_czm_strain_inputs: 活动热网格必须使用 Q4 单元，实际每单元节点数=$(size(thermal_elem, 2))"))
    invalid_thermal_node = findfirst(node -> !(1 <= node <= n_thermal_node), thermal_elem)
    invalid_thermal_node === nothing || throw(ArgumentError(
        "compute_czm_strain_inputs: 活动热网格连接包含越界节点 $(thermal_elem[invalid_thermal_node])，" *
        "有效范围为 1:$n_thermal_node"))

    function current_soc_field(key::String, label::String)
        haskey(variables, key) || throw(KeyError(key))
        field = variables[key]
        current = if field isa AbstractVector
            field
        elseif field isa AbstractMatrix
            size(field, 2) > 0 || throw(DimensionMismatch(
                "compute_czm_strain_inputs: $label SOC 历史矩阵没有当前列，size=$(size(field))"))
            @view field[:, end]
        else
            throw(ArgumentError(
                "compute_czm_strain_inputs: $label SOC 必须是向量或历史矩阵，实际类型=$(typeof(field))"))
        end
        length(current) == n_thermal_elem || throw(DimensionMismatch(
            "compute_czm_strain_inputs: $label SOC 当前场长度=$(length(current))，活动热单元数=$n_thermal_elem"))
        all(isfinite, current) || throw(ArgumentError(
            "compute_czm_strain_inputs: $label SOC 当前场包含非有限值"))
        return current
    end

    soc_p_thermal = current_soc_field("thermal2D element soc_p", "正极")
    soc_n_thermal = current_soc_field("thermal2D element soc_n", "负极")

    # ---- T_czm_nodes 通过 thermal_to_czm 矩阵（nodal interpolation） ----
    T_czm_nodes = M * T_nodes

    # ---- dT_czm 通过 thermal_elem_map（element direct lookup） ----
    dT_thermal = zeros(n_thermal_elem)
    for e_th in 1:n_thermal_elem
        ns = thermal_elem[e_th, :]
        dT_thermal[e_th] = sum(T_nodes[ns]) / 4 - param.cell.T0
    end

    dT_czm = zeros(ne_czm)
    for e in 1:ne_czm
        dT_czm[e] = dT_thermal[thermal_elem_map[e]]
    end

    # ---- Δsoc 通过 thermal_elem_map（element direct lookup） ----
    Δsoc_p_czm = zeros(ne_czm)
    Δsoc_n_czm = zeros(ne_czm)

    for e in 1:ne_czm
        e_th = thermal_elem_map[e]
        mt = submesh.material_type[e]
        if mt == :PE
            Δsoc_p_czm[e] = soc_p_thermal[e_th] - param.PE.cs0
        elseif mt == :NE
            Δsoc_n_czm[e] = soc_n_thermal[e_th] - param.NE.cs0
        end
        # PCC/NCC/SP：保持 0
    end

    return (dT_czm = dT_czm, Δsoc_p_czm = Δsoc_p_czm, Δsoc_n_czm = Δsoc_n_czm,
            T_czm_nodes = T_czm_nodes)
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
- `result`: 完整的 `CZMResult`，包含位移、损伤、分离、牵引力和求解诊断

输入或求解结果包含非有限值、或非线性求解未收敛时抛出异常。
"""
function update_czm_damage!(case, variables, T_nodes_carry)
    czm_mesh = case.czm_mesh
    czm_params = case.param.cohesive
    param_dim = case.param_dim
    param = case.param

    # 同步CZM模型选项（model1 or mix）
    czm_params.czm_model = case.opt.czm_model

    # 计算有效材料参数（per-interface，见 spec §3.5.2）
    czm_param_cache = compute_czm_params_per_interface(case)
    # 缓存到 case 供后续 solve_czm_basic_step 等调用点复用
    case.czm_param_cache === nothing && (case.czm_param_cache = czm_param_cache)
    α_eff = czm_param_cache.by_interface[:PE_PCC].α  # cross-interface uniform per spec §7.1
    β_n = case.param.NE.Omega / 3.0
    β_p = case.param.PE.Omega / 3.0

    # 构建或复用 CZM 缓存（失效判据：objectid(czm_mesh) + param_cache.id + fix_inner）
    cache = ensure_czm_cache(case, czm_mesh, czm_param_cache; fix_inner=case.opt.czm_fix_inner)

    # 计算应变输入
    strain_in = compute_czm_strain_inputs(case, variables, T_nodes_carry)
    dT_elem = strain_in.dT_czm
    Δsoc_n_elem = strain_in.Δsoc_n_czm
    Δsoc_p_elem = strain_in.Δsoc_p_czm

    # 外力向量（一般为零）
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)

    all(isfinite, T_nodes_carry) || throw(ArgumentError("CZM temperature input contains non-finite values"))
    all(isfinite, variables["thermal2D element soc_n"]) || throw(ArgumentError("CZM negative-electrode SOC input contains non-finite values"))
    all(isfinite, variables["thermal2D element soc_p"]) || throw(ArgumentError("CZM positive-electrode SOC input contains non-finite values"))
    all(isfinite, dT_elem) || throw(ArgumentError("CZM temperature increment contains non-finite values"))
    all(isfinite, Δsoc_n_elem) || throw(ArgumentError("CZM negative-electrode SOC increment contains non-finite values"))
    all(isfinite, Δsoc_p_elem) || throw(ArgumentError("CZM positive-electrode SOC increment contains non-finite values"))

    # 初始化位移（从 czm_layout 获取上一步值）
    u_czm_prev = case.czm_layout.u_prev

    # 调用 CZM 求解器（可选迭代方式）
    iter_method = case.opt.czm_iter_method
    max_iter = case.opt.czm_max_iter
    tol = case.opt.czm_tol
    n_load_steps = case.opt.czm_load_steps
    arc_length_alpha = case.opt.czm_arc_length_alpha

    # Viscous regularization: compute visc_beta
    # β = Δs / (τ_v* + Δs), where Δs = 1/n_load_steps for load_substep, 1.0 for basic
    visc_beta = 1.0  # default: no regularization
    if case.opt.czm_viscous_enabled && czm_params.tau_visc > 0.0
        delta_s = if lowercase(iter_method) == "basic"
            1.0
        elseif lowercase(iter_method) in ("arc_length", "arclength", "arc-length")
            1.0 / max(1, n_load_steps)  # approximate; actual delta_lambda varies per substep
        else  # load_substep
            1.0 / max(1, n_load_steps)
        end
        visc_beta = delta_s / (czm_params.tau_visc + delta_s)
    end

    geo_nl = case.opt.czm_geo_nonlinear
    eig = geo_nl ? (α_eff=α_eff, β_n=β_n, β_p=β_p,
                    dT=dT_elem, Δsn=Δsoc_n_elem, Δsp=Δsoc_p_elem) : nothing
    prestress = nothing
    if case.opt.czm_winding_prestress
        geo_nl || error("update_czm_damage!: czm_winding_prestress=true 需要 czm_geo_nonlinear=true（D-B2'-2）。")
        (param.cell.winding_T_ne > 0.0 || param.cell.winding_T_pe > 0.0) || error(
            "update_czm_damage!: czm_winding_prestress=true 但 cell.winding_T_ne/pe 均为 0（未设置）。缺参即拦截（AGENTS 9.4/9.7）。")
        if case.czm_layout.winding_prestress === nothing
            case.czm_layout.winding_prestress = winding_prestress_field(czm_mesh, param)
        end
        prestress = case.czm_layout.winding_prestress
    end
    plastic_on = case.opt.czm_j2_plasticity
    mech_state = nothing
    if plastic_on
        geo_nl || error("update_czm_damage!: czm_j2_plasticity=true 需要 czm_geo_nonlinear=true（D-B3-1）。")
        (param.PCC.sigma_y > 0.0 && param.NCC.sigma_y > 0.0) || error(
            "update_czm_damage!: czm_j2_plasticity=true 但 PCC/NCC 的 sigma_y ≤ 0（未设置）。缺参即拦截，不默认、不置零（AGENTS 9.4/9.7）。")
        if case.czm_layout.plastic_states === nothing
            case.czm_layout.plastic_states = [PlasticState() for _ in 1:size(czm_mesh.bulk_element, 1), _ in 1:4]
        end
        mech_state = case.czm_layout.plastic_states
    end

    result, updated_czm_mesh = solve_czm_step(
        czm_mesh, F_ext, czm_param_cache, param, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, iter_method=iter_method,
        cache=cache, visc_beta=visc_beta,
        geo_nl=geo_nl, eigenstrain=eig,
        plasticity=plastic_on, mech_state=mech_state, prestress=prestress
    )

    if case.opt.debug_coupling
        pre_stats = get_damage_statistics(czm_mesh)
        trial_stats = get_damage_statistics(updated_czm_mesh)
        max_delta_n = isempty(result.separation_n) ? 0.0 : maximum(result.separation_n)
        max_delta_eff = 0.0
        for (δ_n, δ_t) in zip(result.separation_n, result.separation_t)
            if czm_params.czm_model == "model1"
                max_delta_eff = max(max_delta_eff, max(δ_n, 0.0))
            else
                max_delta_eff = max(max_delta_eff, sqrt(max(δ_n, 0.0)^2 + δ_t^2))
            end
        end
        D_visc_max_trial = isempty(updated_czm_mesh.damage_states) ? 0.0 : maximum(ds.D_visc for ds in updated_czm_mesh.damage_states)
        println("[CZM-Debug] max(δ_n)=$(round(max_delta_n; digits=6)), max(δ_eff)=$(round(max_delta_eff; digits=6)), converged=$(result.converged), β=$(round(visc_beta; digits=4)), D_max(before)=$(round(pre_stats.max_D; digits=6)), D_max(trial)=$(round(trial_stats.max_D; digits=6)), D_visc_max(trrial)=$(round(D_visc_max_trial; digits=6)), D_max(commit)=$(result.converged ? round(trial_stats.max_D; digits=6) : round(pre_stats.max_D; digits=6))")
    end

    all(isfinite, result.displacement) || error("CZM solve returned non-finite displacement")
    all(isfinite, result.damage) || error("CZM solve returned non-finite damage")
    all(isfinite, result.separation_n) || error("CZM solve returned non-finite normal separation")
    all(isfinite, result.separation_t) || error("CZM solve returned non-finite tangential separation")
    all(isfinite, result.traction_n) || error("CZM solve returned non-finite normal traction")
    all(isfinite, result.traction_t) || error("CZM solve returned non-finite tangential traction")
    isfinite(result.residual_norm) || error("CZM solve returned a non-finite residual norm")
    result.converged || error("CZM solve did not converge after $(result.iterations) iterations (residual=$(result.residual_norm))")

    czm_mesh.damage_states = updated_czm_mesh.damage_states
    if plastic_on
        assemble_coupled_system(czm_mesh, result.displacement, czm_param_cache;
            geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=mech_state,
            commit_plastic=true)
    end
    case.czm_layout.u_prev = result.displacement

    return result
end

# ========================================================================
# CZMSnapshot — per-step CZM state for CSV export
# ========================================================================

"""
    CZMSnapshot

Stores per-step CZM solver state for CSV export.
All physical values are stored in NORMALIZED (dimensionless) form.
Denormalization happens at CSV write time using `case.param.scale`.
"""
mutable struct CZMSnapshot
    time_s::Float64                     # physical time (already denormalized)
    cycle::Int                          # cycle number
    phase::String                       # phase name
    displacement::Vector{Float64}       # ndof-length, normalized
    damage::Vector{Float64}             # n_coh-length, [0,1]
    separation_n::Vector{Float64}       # n_coh-length, normalized
    separation_t::Vector{Float64}       # n_coh-length, normalized
    traction_n::Vector{Float64}         # n_coh-length, normalized
    traction_t::Vector{Float64}         # n_coh-length, normalized
    converged::Bool
    iterations::Int
    residual_norm::Float64
    method::String                      # "basic", "load_substep", or "arc_length"
end

# ========================================================================
# build_thermal_to_czm_interp — 粗热节点 → CZM 节点双线性插值矩阵
# (spec v2 §5.1)
# ========================================================================

"""
    build_thermal_to_czm_interp(thermal_mesh::Mesh, czm_submesh::CzmSubmesh) -> SparseMatrixCSC

构造粗热节点 → CZM 节点双线性插值矩阵（n_czm_node × n_thermal_node）。
每行 ≤4 个非零元（粗热 Q4 单元的 4 节点），权重 = 双线性形函数值，行和 = 1。

# 算法
力学网格与热网格共享周向节点。按结构化节点列确定唯一父热单元，
仅在该父单元内求等参坐标；定位失败立即报错。
"""
function build_thermal_to_czm_interp(thermal_mesh::Mesh, czm_submesh::CzmSubmesh)
    czm_node = czm_submesh.mesh.node
    n_czm_node = czm_submesh.mesh.nlen
    n_thermal_node = thermal_mesh.nlen
    n_thermal_elem = size(thermal_mesh.element, 1)

    n_theta_nodes = n_thermal_elem + 1
    n_czm_node % n_theta_nodes == 0 || throw(DimensionMismatch(
        "build_thermal_to_czm_interp: CZM node count $n_czm_node is not divisible by thermal angular node count $n_theta_nodes"))
    n_spirals = div(n_czm_node, n_theta_nodes)
    n_spirals >= 9 || throw(DimensionMismatch(
        "build_thermal_to_czm_interp: expected >= 9 mechanical spirals (8 base + thin subdiv), got $n_spirals"))

    # 双线性形函数 N_i(ξ, η) = 0.25 * (1 ± ξ)(1 ± η)
    function shape_funcs(ξ::Float64, η::Float64)
        return [
            0.25 * (1 - ξ) * (1 - η),
            0.25 * (1 + ξ) * (1 - η),
            0.25 * (1 + ξ) * (1 + η),
            0.25 * (1 - ξ) * (1 + η),
        ]
    end

    # Newton 迭代解等参坐标
    function solve_isoparametric(x_nodes, y_nodes, px, py; max_iter=20, tol=1e-13)
        ξ, η = 0.0, 0.0
        for _ in 1:max_iter
            N = shape_funcs(ξ, η)
            x_pred = sum(N .* x_nodes)
            y_pred = sum(N .* y_nodes)
            rx = px - x_pred
            ry = py - y_pred
            if abs(rx) < tol && abs(ry) < tol
                return ξ, η, true
            end
            dN_dξ = 0.25 * [-(1-η), (1-η), (1+η), -(1+η)]
            dN_dη = 0.25 * [-(1-ξ), -(1+ξ), (1+ξ), (1-ξ)]
            Jxξ = sum(dN_dξ .* x_nodes)
            Jxη = sum(dN_dη .* x_nodes)
            Jyξ = sum(dN_dξ .* y_nodes)
            Jyη = sum(dN_dη .* y_nodes)
            detJ = Jxξ * Jyη - Jxη * Jyξ
            abs(detJ) < 1e-20 && break
            ξ += (Jyη * rx - Jxη * ry) / detJ
            η += (-Jyξ * rx + Jxξ * ry) / detJ
        end
        return ξ, η, abs(ξ) <= 1.0 + 1e-6 && abs(η) <= 1.0 + 1e-6
    end

    I_rows = Int[]
    J_cols = Int[]
    V_vals = Float64[]

    for i in 1:n_czm_node
        px, py = czm_node[i, 1], czm_node[i, 2]
        angular_node = mod(i - 1, n_theta_nodes) + 1
        parent = min(angular_node, n_thermal_elem)
        ns = thermal_mesh.element[parent, :]
        x_nodes = thermal_mesh.node[ns, 1]
        y_nodes = thermal_mesh.node[ns, 2]
        ξ, η, ok = solve_isoparametric(x_nodes, y_nodes, px, py)
        if !ok || abs(ξ) > 1.0 + 1e-6 || abs(η) > 1.0 + 1e-6
            error("build_thermal_to_czm_interp: CZM node $i at ($px, $py) is outside its topological parent thermal element $parent (ξ=$ξ, η=$η)")
        end
        N = shape_funcs(clamp(ξ, -1, 1), clamp(η, -1, 1))
        for k in 1:4
            push!(I_rows, i)
            push!(J_cols, ns[k])
            push!(V_vals, N[k])
        end
    end

    return sparse(I_rows, J_cols, V_vals, n_czm_node, n_thermal_node)
end
