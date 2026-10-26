"""
    MechanicalLoad(czm_mesh; dT, Δsoc_n, Δsoc_p, F_ext, dt_seconds, sample_time_seconds)

One sampled mechanical load. dT = ΔT/T_ref; Δsoc_n and Δsoc_p are differences
of each electrode's c/c_max from its reference value. F_ext keeps the solver's
existing dimensionless nodal-force units. sample_time_seconds is the source
electrochemical/thermal state's global physical time; dt_seconds separately
accumulates elapsed mechanical time. Arrays are borrowed read-only; static
loads explicitly use dt_seconds=0.
"""
struct MechanicalLoad
    dT::Vector{Float64}
    Δsoc_n::Vector{Float64}
    Δsoc_p::Vector{Float64}
    F_ext::Vector{Float64}
    dt_seconds::Float64
    sample_time_seconds::Float64

    function MechanicalLoad(czm_mesh::CohesiveMesh;
            dT::Vector{Float64}, Δsoc_n::Vector{Float64}, Δsoc_p::Vector{Float64},
            F_ext::Vector{Float64}, dt_seconds::Real, sample_time_seconds::Real)
        ne = size(czm_mesh.bulk_element, 1)
        all(a -> length(a) == ne, (dT, Δsoc_n, Δsoc_p)) ||
            throw(DimensionMismatch("mechanical load fields must match $ne bulk elements"))
        length(F_ext) == 2 * size(czm_mesh.node, 1) ||
            throw(DimensionMismatch("mechanical external force must match the mesh DOFs"))
        all(isfinite, dT) || throw(ArgumentError("CZM temperature increment contains non-finite values"))
        all(isfinite, Δsoc_n) || throw(ArgumentError("CZM negative-electrode SOC increment contains non-finite values"))
        all(isfinite, Δsoc_p) || throw(ArgumentError("CZM positive-electrode SOC increment contains non-finite values"))
        all(isfinite, F_ext) || throw(ArgumentError("mechanical external force contains non-finite values"))
        dt = Float64(dt_seconds)
        sample_time = Float64(sample_time_seconds)
        isfinite(dt) && dt >= 0.0 ||
            throw(ArgumentError("mechanical elapsed time must be finite and nonnegative seconds"))
        isfinite(sample_time) ||
            throw(ArgumentError("mechanical sampling time must be finite global seconds"))
        return new(dT, Δsoc_n, Δsoc_p, F_ext, dt, sample_time)
    end
end

"""Session-owned CZM caches; never part of the accepted physical state."""
mutable struct CZMWorkspace
    mesh::CohesiveMesh
    param::Params
    K_bulk::Union{Nothing, SparseMatrixCSC}
    cohesive_geom::Union{Nothing, Vector}
    bulk_gp_geom::Union{Nothing, Vector}
    bulk_asm::Any
    ws::Any
end

CZMWorkspace(mesh::CohesiveMesh, param::Params) =
    CZMWorkspace(mesh, param, nothing, nothing, nothing, nothing, nothing)

function check_czm_workspace(mesh::CohesiveMesh, workspace::CZMWorkspace)
    workspace.mesh === mesh || throw(ArgumentError("CZM workspace belongs to a different mesh; create a new session workspace"))
    return workspace
end

using Parameters: @with_kw

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
            key = minmax(a, b)
            key in seen && continue
            push!(seen, key)
            push!(edges, key)
            push!(L_edge, hypot(x[b] - x[a], y[b] - y[a]))
        end
    end
    return BoundaryEdgeCache(edges, L_edge)
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
    # basic 非几何路径的 BC 后矩阵快照与原生 factorize 结果
    K_bc_factor_matrix::Union{Nothing, SparseMatrixCSC{Float64, Int64}}
    K_bc_factorization::Any

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
            K_coh_sp,
            nothing,
            nothing)
    end
end

"""
    DamageState - 内聚力单元的损伤状态

存储每个内聚力单元（或高斯点）的损伤历史。
（2026-08-30 重构自 czm.jl 前移至此：MechState 字段需要它在 include 顺序中
早于本文件可用。）

# 字段
- `D`: 不可逆包络损伤 D_eq [0, 1]，历次当前方向双线性损伤的最大值
- `D_visc`: 粘性有效损伤，用于牵引力和切线（粘性正则化关闭时 D_visc = D）
- `δ_max_n`: 历史最大正法向分离位移（几何统计，不控制损伤门）
- `δ_max_t`: 历史最大切向分离绝对值（几何统计）
- `δ_max_eff`: 历史最大等效分离位移（几何统计）
- `fractured`: 退化方向已饱和，规范状态须 D = D_visc = 1；model1 切向及压缩刚度保留
- `accumulated_damage`: 本构更新取历史最大值；外部循环接口另见73，不代表能量或独立疲劳模型
"""
mutable struct DamageState
    D::Float64                     # 等效损伤 D_eq
    D_visc::Float64                # 粘性有效损伤（用于牵引和切线）
    δ_max_n::Float64              # 历史最大法向分离
    δ_max_t::Float64              # 历史最大切向分离
    δ_max_eff::Float64            # 历史最大等效分离
    fractured::Bool               # 是否断裂
    accumulated_damage::Float64   # 累积损伤（循环）

    # 默认构造函数
    DamageState() = new(0.0, 0.0, 0.0, 0.0, 0.0, false, 0.0)
end

"""
    MechState

CZM 求解的演化状态聚合（2026-08-30 重构：吸收 CzmLayout；损伤态从 czm_mesh 迁入）。
缓存归会话 CZMWorkspace，演化状态随求解（本结构）。
求解器在收敛后原位提交；失败/试探不触碰本结构。
"""
mutable struct MechState
    u_prev::Vector{Float64}               # 上一步收敛位移（跨时间步持有）
    damage_states::Vector{DamageState}    # 逐 cohesive 单元输出聚合；不得反推GP恢复历史
    gp_damage_states::Union{Nothing, Matrix{DamageState}} # 真实逐GP接受史；仅首次推进前可缺席
    plastic_states::Union{Nothing, Matrix{PlasticState}}  # PCC/NCC 高斯点塑性状态（Batch 3，[ne,4]；收敛才提交 D-B3-2）
    winding_prestress::Union{Nothing, Vector{NTuple{3, Float64}}}  # 卷绕预应力 σ₀（Batch 2'；几何固定一次计算持久持有）
    node_ref::Union{Nothing, Matrix{Float64}}   # 初始螺旋节点快照（Δ_core 基准，永不重置）
    mechanics_initialized::Bool            # true after an actual successful mechanical advance
    last_eigenstrain::Union{Nothing,NamedTuple{(:dT,:Δsn,:Δsp),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}} # actual last successful bulk load
end

"""便捷构造器：从 czm_mesh 初始化"""
function MechState(czm_mesh::CohesiveMesh)
    ndof = 2 * size(czm_mesh.node, 1)
    MechState(zeros(Float64, ndof),
              [DamageState() for _ in 1:length(czm_mesh.cohesive_elements)],
              nothing, nothing, nothing, nothing, false, nothing)
end

function require_mechanical_load_history(ms::MechState)
    ms.mechanics_initialized && ms.last_eigenstrain===nothing &&
        throw(ArgumentError("initialized mechanics lacks successful eigenstrain history"))
    return ms.last_eigenstrain
end

"""
    ensure_node_ref!(case, mech) -> Matrix{Float64}

初始螺旋节点快照（惰性建立、永不重置）：Δ_core 始终相对初始完美螺旋计算（spec §3.5）。
"""
function ensure_node_ref!(case, mech::MechState)
    ref = mech.node_ref === nothing ? copy(case.mesh["czm"].node) : mech.node_ref
    size(ref)==size(case.mesh["czm"].node) || throw(DimensionMismatch("reference-node shape"))
    all(isfinite,ref) || throw(ArgumentError("non-finite reference nodes"))
    mech.node_ref=ref
    return ref
end

"""
    core_ovalization(czm_mesh, u, ref_node) -> (w_core, Δ_core)

spec §3.5（D8 冻结定义）位移基 Δ_core：Γ_in,free = 内螺旋第一匝窗口节点，
u_n = u·r̂（r̂ 取初始螺旋节点向径，始终相对初始螺旋）；窗口内离散 Fourier 去除
n=0（均值）与 n=1（刚体平移）后取 max|ũ_n|；Δ_core = w_core / cell.Rin。
"""
function core_ovalization(czm_mesh, u, ref_node::AbstractMatrix{<:Real})
    bonded = czm_mesh.czm_submesh.mesh_bonded
    # 单匝窗口 = 每匝角节点数（第 1 匝分段数，由 winding_turn 反推；element[1,2] 是整条
    # 螺旋角节点总数跨全部匝，不得使用）；节点 1..nθn 天然按角序（θ>π 处 atan 断跳，禁排序）
    wt = czm_mesh.czm_submesh.winding_turn
    nθn = count(==(1), wt[1:count(==(wt[1]), wt)])
    r_ref = hypot(ref_node[1, 1], ref_node[1, 2])
    u_n = zeros(nθn)
    for (i, k) in enumerate(1:nθn)        x, y = ref_node[k, 1], ref_node[k, 2]
        c, s = x / hypot(x, y), y / hypot(x, y)
        u_n[i] = c * u[2*k-1] + s * u[2*k]
    end
    u_n .-= sum(u_n) / nθn                    # 去 n=0
    # 去 n=1（最小二乘一阶谐波）
    a1 = 2 / nθn * sum(u_n[i] * cos(2π * (i - 1) / nθn) for i in 1:nθn)
    b1 = 2 / nθn * sum(u_n[i] * sin(2π * (i - 1) / nθn) for i in 1:nθn)
    u_n .-= [a1 * cos(2π * (i - 1) / nθn) + b1 * sin(2π * (i - 1) / nθn) for i in 1:nθn]
    w_core = maximum(abs.(u_n))
    return w_core, w_core / r_ref
end

# ========================================================================
# CZM Coupling Helpers
# ========================================================================
# 参数计算、应变输入、损伤更新入口。
# 从 CzmSolve.jl 迁移至此，保持 CouplingState.jl 作为"状态+耦合"的统一归属点。
# 函数调用 solve_czm_step 及 mesh 缓存访问器在运行时解析，无需调整 include 顺序。
# ========================================================================


"""
    compute_czm_strain_inputs(case, variables, T_nodes) -> NamedTuple

计算 CZM 体单元粒度的 dT、Δsoc_p、Δsoc_n。

# 数据流
- dT_czm：粗热单元 dT (= avg(T_nodes[thermal_elem[e]]) - T0)
          → thermal_elem_map 直接取值（element-to-element）
- Δsoc_p/n：variables["thermal2D element soc_p/n"] (粗热单元无量纲 c/c_max)
             → thermal_elem_map 直接取值 - cs0，按 material_type 分发

# 单位契约
- T_nodes、param.cell.T0 为 T/T_ref；dT_czm 为 ΔT/T_ref，例如 12.3 K 输入 12.3/T_ref。
- soc_p/n, cs0, Δsoc: 无量纲 c/c_max；Omega 已乘 cs_max，ΩΔsoc/3 的尺度契约不变。
"""
function compute_czm_strain_inputs(case, variables, T_nodes)
    czm_mesh = case.mesh["czm"]
    submesh = czm_mesh.czm_submesh
    param = case.param
    ne_czm = size(submesh.mesh.element, 1)

    thermal_mesh = case.mesh["thermal2D"]
    thermal_elem = thermal_mesh.element
    n_thermal_node = thermal_mesh.nlen
    n_thermal_elem = size(thermal_elem, 1)

    # ---- 在任何索引前验证父热单元场 → CZM 的完整接口契约 ----
    T_nodes isa AbstractVector || throw(DimensionMismatch(
        "compute_czm_strain_inputs: T_nodes 必须是一维活动热节点向量，实际 size=$(size(T_nodes))"))
    length(T_nodes) == n_thermal_node || throw(DimensionMismatch(
        "compute_czm_strain_inputs: 热节点接口尺寸不一致：活动热网格节点数=$(n_thermal_node)，" *
        "T_nodes 长度=$(length(T_nodes))"))
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

    return (dT_czm = dT_czm, Δsoc_p_czm = Δsoc_p_czm, Δsoc_n_czm = Δsoc_n_czm)
end

"""
    update_czm_damage!(case, mech, variables, T_nodes_carry; dt_seconds, sample_time_seconds, czm_workspace)

更新 CZM 网格的损伤状态。

使用牛顿-拉弗森迭代求解力学平衡方程，通过载荷子步法处理软化收敛问题。
从 case.mesh["czm"]、case.param、mech 与 case.opt.czm 获取拓扑、参数、演化状态和求解选项。

# 参数
- `case`: Case 对象（需已设置 czm_mesh 和 mech）
- `variables`: 当前时间步的变量字典
- `T_nodes_carry`: 已接受左端温度场
- `dt_seconds`: 自上次机械推进累计的实际秒数
- `sample_time_seconds`: 载荷采样的全局物理秒数

# 返回
- `result`: 完整的 `CZMResult`，包含位移、损伤、分离、牵引力和求解诊断

输入或求解结果包含非有限值、或非线性求解未收敛时抛出异常。
"""
function update_czm_damage!(case, mech::MechState, variables, T_nodes_carry; czm_workspace::CZMWorkspace, dt_seconds::Real, sample_time_seconds::Real)
    czm_mesh = case.mesh["czm"]
    param = case.param
    czm_opt = case.opt.czm

    # 计算应变输入
    strain_in = compute_czm_strain_inputs(case, variables, T_nodes_carry)
    dT_elem = strain_in.dT_czm
    Δsoc_n_elem = strain_in.Δsoc_n_czm
    Δsoc_p_elem = strain_in.Δsoc_p_czm

    # 外力向量（一般为零）
    ndof = 2 * size(czm_mesh.node, 1)
    F_ext = zeros(Float64, ndof)

    all(isfinite, T_nodes_carry) || throw(ArgumentError("CZM temperature input contains non-finite values"))
    all(isfinite, variables["thermal2D element soc_n"]) || throw(ArgumentError("CZM negative-electrode SOC input contains non-finite values"))
    all(isfinite, variables["thermal2D element soc_p"]) || throw(ArgumentError("CZM positive-electrode SOC input contains non-finite values"))
    load = MechanicalLoad(czm_mesh; dT=dT_elem, Δsoc_n=Δsoc_n_elem,
        Δsoc_p=Δsoc_p_elem, F_ext=F_ext, dt_seconds=dt_seconds,
        sample_time_seconds=sample_time_seconds)

    # 初始化位移与配置由 ms/czm_opt 携带；粘性系数按实际经过的物理秒数计算。
    geo_nl = czm_opt.geo_nonlinear
    prestress = nothing
    if czm_opt.winding_prestress
        geo_nl || error("update_czm_damage!: winding_prestress=true 需要 geo_nonlinear=true（D-B2'-2）。")
        (param.cell.winding_T_ne > 0.0 || param.cell.winding_T_pe > 0.0) || error(
            "update_czm_damage!: winding_prestress=true 但 cell.winding_T_ne/pe 均为 0（未设置）。缺参即拦截（AGENTS 9.4/9.7）。")
        if mech.winding_prestress === nothing
            mech.mechanics_initialized && throw(ArgumentError("initialized mechanics lacks winding prestress"))
            mech.winding_prestress = winding_prestress_field(czm_mesh, param)
        end
        prestress = mech.winding_prestress
    end
    if czm_opt.j2_plasticity
        geo_nl || error("update_czm_damage!: j2_plasticity=true 需要 geo_nonlinear=true（D-B3-1）。")
        (param.PCC.sigma_y > 0.0 && param.NCC.sigma_y > 0.0) || error(
            "update_czm_damage!: j2_plasticity=true 但 PCC/NCC 的 sigma_y ≤ 0（未设置）。缺参即拦截，不默认、不置零（AGENTS 9.4/9.7）。")
    end

    # 求解器收敛后在 mech 上原位提交损伤/位移（失败不触碰）
    result = solve_czm_step(
        czm_mesh, mech, param, load, czm_opt; czm_workspace=czm_workspace,
        prestress=prestress
    )

    all(isfinite, result.displacement) || error("CZM solve returned non-finite displacement")
    all(isfinite, result.damage) || error("CZM solve returned non-finite damage")
    all(isfinite, result.separation_n) || error("CZM solve returned non-finite normal separation")
    all(isfinite, result.separation_t) || error("CZM solve returned non-finite tangential separation")
    all(isfinite, result.traction_n) || error("CZM solve returned non-finite normal traction")
    all(isfinite, result.traction_t) || error("CZM solve returned non-finite tangential traction")
    isfinite(result.residual_norm) || error("CZM solve returned a non-finite residual norm")
    result.converged || error("CZM solve did not converge after $(result.iterations) iterations (residual=$(result.residual_norm))")

    if case.opt.debug_coupling
        stats = get_damage_statistics(mech.damage_states)
        max_delta_n = isempty(result.separation_n) ? 0.0 : maximum(result.separation_n)
        max_delta_eff = 0.0
        for (δ_n, δ_t) in zip(result.separation_n, result.separation_t)
            if czm_opt.model == "model1"
                max_delta_eff = max(max_delta_eff, max(δ_n, 0.0))
            else
                max_delta_eff = max(max_delta_eff, sqrt(max(δ_n, 0.0)^2 + δ_t^2))
            end
        end
        D_visc_max = isempty(mech.damage_states) ? 0.0 : maximum(ds.D_visc for ds in mech.damage_states)
        println("[CZM-Debug] max(δ_n)=$(round(max_delta_n; digits=6)), max(δ_eff)=$(round(max_delta_eff; digits=6)), converged=$(result.converged), D_max(commit)=$(round(stats.max_D; digits=6)), D_visc_max=$(round(D_visc_max; digits=6))")
    end

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
    time_s::Float64                     # phase-local physical time [s]; add PhaseResult.t_start for global CSV time
    cycle::Int                          # cycle number
    phase::String                       # phase name
    displacement::Vector{Float64}       # ndof-length, normalized
    damage::Vector{Float64}             # n_coh-length, [0,1]
    damage_visc::Vector{Float64}         # real accepted viscous damage
    fractured::BitVector               # all actual GP flags, never a D threshold
    separation_n::Vector{Float64}       # n_coh-length, normalized
    separation_t::Vector{Float64}       # n_coh-length, normalized
    traction_n::Vector{Float64}         # n_coh-length, normalized
    traction_t::Vector{Float64}         # n_coh-length, normalized
    converged::Bool
    iterations::Int
    residual_norm::Float64
    method::String                      # "basic", "load_substep", or "arc_length"
    mechanical_dt_s::Float64            # actual accumulated mechanical interval [s]
    czm_model::String                   # model used by this accepted step
    gp_damage_states::Matrix{DamageState} # copied accepted history, indexed (COH element, GP)
    gp_points::Vector{Float64}          # actual NCweight point order
    gp_separation_n::Matrix{Float64}    # current separations, not history peaks
    gp_separation_t::Matrix{Float64}
    gp_traction_n::Matrix{Float64}      # fixed accepted history, visc_beta=0
    gp_traction_t::Matrix{Float64}
    plasticity::Bool
    bulk_kappa::Union{Nothing,Matrix{Float64}} # separate (bulk element, bulk GP); absent with J2 off
end
