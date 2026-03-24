# MultiSPMe 状态结构重构设计

## 概述

重构 `Initialisation.jl` 中的初始化逻辑，将电化学状态从一维向量改为矩阵形式，简化代码并消除重复。

**范围限定**：本重构仅针对 SPMe + distributed2D 热模型的多SPMe模式。SPM、单SPMe、P2D 模式保持原有向量形式不变。

## 动机

### 当前问题

1. **代码重复**：`ModelInitialisation` 和 `ModelInitialisation_MultiSPMe` 有大量重复逻辑
2. **状态结构复杂**：一维向量 `[yt_chem[1]; ...; yt_chem[ne]; T_nodes]` 需要复杂的索引管理
3. **布局信息冗余**：`multi_spme_layout` 存储的可计算信息（`ne`, `n_chem`, `nT`, `chem_range`, `thermal_range`）
4. **辅助函数过多**：`MultiSPMe_extract_element_state`, `MultiSPMe_get_thermal_dofs`, `MultiSPMe_update_state` 仅用于状态切片

### 目标

1. 合并初始化函数为一个 `ModelInitialisation`
2. 电化学状态改为矩阵 `Y_chem::Matrix{Float64}` (n_chem × ne)
3. 热场状态独立为 `T_nodes::Vector{Float64}` (nT)
4. 删除 `multi_spme_layout` 中的可计算信息
5. 时间积分器支持分块状态结构

## 设计

### 1. 状态结构

新增 `MultiSPMeState` 结构体（定义在 `Initialisation.jl` 顶部）：

```julia
"""
多SPMe状态结构：电化学矩阵 + 热场向量

仅用于 SPMe + distributed2D 热模型模式。
"""
struct MultiSPMeState
    Y_chem::Matrix{Float64}   # (n_chem × ne) 每列是一个单元的电化学状态
    T_nodes::Vector{Float64}  # (nT,) 热场节点温度
end

# 辅助函数
Base.length(s::MultiSPMeState) = length(s.Y_chem) + length(s.T_nodes)
n_elements(s::MultiSPMeState) = size(s.Y_chem, 2)
n_chem_dofs(s::MultiSPMeState) = size(s.Y_chem, 1)

# 转换函数（用于与现有代码兼容）
function to_vector(state::MultiSPMeState)
    return [vec(state.Y_chem); state.T_nodes]
end

function from_vector(y::Vector{Float64}, n_chem::Int, ne::Int, nT::Int)
    Y_chem = reshape(y[1:n_chem*ne], n_chem, ne)
    T_nodes = y[n_chem*ne+1:end]
    return MultiSPMeState(Y_chem, T_nodes)
end
```

**设计决策**：
- 矩阵每列对应一个热单元的电化学状态
- Julia 按列存储，访问 `Y_chem[:, e]` 高效
- 热场保持向量形式（全局耦合）
- 使用不可变 `struct`（函数式风格，每次返回新状态）

### 2. 初始化函数（合并后）

```julia
function ModelInitialisation(case::Case; initial_soc_distribution::Union{Nothing, Vector{Float64}}=nothing)
    # 获取基础电化学状态（单单元）
    y0_chem = init_chem_state(case)

    if should_use_multi_spme(case)
        # 多SPMe模式：矩阵形式
        ne = get_ne(case)
        n_chem = length(y0_chem)
        nT = get_nT(case)

        if initial_soc_distribution === nothing
            # 均匀初始化
            Y_chem = repeat(y0_chem, 1, ne)
        else
            # 非均匀初始化：根据 SOC 分布调整各单元浓度
            Y_chem = zeros(n_chem, ne)
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            for e in 1:ne
                soc_e = initial_soc_distribution[e]
                cn_surf_e = case.param.NE.cs0 * soc_e
                cp_surf_e = case.param.PE.cs0 * (1.0 - soc_e)
                Y_chem[1:Nrn, e] .= cn_surf_e
                Y_chem[Nrn+1:Nrn+Nrp, e] .= cp_surf_e
                Y_chem[Nrn+Nrp+1:end, e] .= case.param.EL.ce0
            end
        end

        T_nodes = fill(case.param.cell.T0, nT)
        return MultiSPMeState(Y_chem, T_nodes)
    else
        # 单SPMe模式：保持原有向量形式
        y0 = y0_chem
        if case.opt.thermalmodel == "lumped"
            y0 = [y0; case.param.cell.T0]
        elseif case.opt.thermalmodel == "distributed2D"
            y0 = [y0; fill(case.param.cell.T0, get_nT(case))]
        end
        # P2D 额外自由度处理...
        return y0
    end
end

# 判断是否使用多SPMe
function should_use_multi_spme(case::Case)
    return case.opt.model == "SPMe" && case.opt.thermalmodel == "distributed2D"
end

# 子函数：初始化单单元电化学状态
function init_chem_state(case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nel = case.mesh["electrolyte"].nlen

    csn0 = ones(Float64, Nrn) * case.param.NE.cs0
    csp0 = ones(Float64, Nrp) * case.param.PE.cs0

    if case.opt.model in ["SPMe", "P2D"]
        ce0 = ones(Float64, Nel) * case.param.EL.ce0
        return [csn0; csp0; ce0]
    else  # SPM
        return [csn0; csp0]
    end
end
```

### 3. 辅助函数（替代 multi_spme_layout）

```julia
# 实时计算布局信息
function get_ne(case::Case)
    haskey(case.mesh, "thermal2D") ? size(case.mesh["thermal2D"].element, 1) : 1
end

function get_nT(case::Case)
    haskey(case.mesh, "thermal2D") ? case.mesh["thermal2D"].nlen : 0
end

function get_n_chem(case::Case)
    case.mesh["negative particle"].nlen +
    case.mesh["positive particle"].nlen +
    case.mesh["electrolyte"].nlen
end

# 计算单元均温（从节点温度）
function compute_element_temperatures(mesh_th, T_nodes::Vector{Float64})
    ne = size(mesh_th.element, 1)
    Te = zeros(Float64, ne)
    for e in 1:ne
        nds = mesh_th.element[e, :]
        Te[e] = sum(T_nodes[nds]) / length(nds)
    end
    return Te
end
```

### 4. 时间积分器

```julia
function solve_step_multi_spme(state::MultiSPMeState, case, t, dt, θ)
    ne = n_elements(state)
    n_chem = n_chem_dofs(state)

    # 1. 调用 CallModel_MultiSPMe 获取分块矩阵
    M_blocks, K_blocks, F_matrix, MT, KT, FT, vars = CallModel_MultiSPMe_blocked(state, case, t)

    # 2. 并行求解电化学（单元间无耦合）
    Y_chem_new = similar(state.Y_chem)
    Threads.@threads for e in 1:ne
        Me, Ke = M_blocks[e], K_blocks[e]
        fe = F_matrix[:, e]
        LHS = Me + θ * dt * Ke
        RHS = (Me - (1-θ) * dt * Ke) * state.Y_chem[:, e] + dt * fe
        Y_chem_new[:, e] = LHS \ RHS
    end

    # 3. 求解热场（全局耦合）
    LHS_T = MT + θ * dt * KT
    RHS_T = (MT - (1-θ) * dt * KT) * state.T_nodes + dt * FT
    T_nodes_new = LHS_T \ RHS_T

    return MultiSPMeState(Y_chem_new, T_nodes_new), vars
end
```

**设计决策**：
- 电化学单元间无耦合，可独立求解（并行）
- 热场全局耦合，需要全局求解
- 保持强耦合策略（同一时间步内求解）

### 5. CallModel_MultiSPMe 修改

```julia
function CallModel_MultiSPMe_blocked(state::MultiSPMeState, case::Case, t::Float64; jacobi::String="update")
    ne = n_elements(state)
    n_chem = n_chem_dofs(state)
    T_nodes = state.T_nodes
    mesh_th = case.mesh["thermal2D"]

    # 初始化 variables
    variables = StandardVariables(case, 1)

    # 1. 计算单元面积和均温
    A = compute_element_areas(mesh_th)
    variables["thermal2D element area"] = A
    Te_prev = compute_element_temperatures(mesh_th, T_nodes)

    # 2. 分流求解
    I_total = case.opt.Current(t * case.param.scale.t0) / case.param.scale.I_typ
    variables["cell current"] = I_total
    variables["T_nodes"] = T_nodes

    yt_representative = mean(state.Y_chem, dims=2)[:, 1]
    T_rep = mean(Te_prev)
    vars_rep = SPMe_variables(case, yt_representative, t; I_app=I_total, T_e=T_rep)
    for (k, v) in vars_rep
        variables[k] = v
    end
    variables["T_nodes"] = T_nodes
    variables["thermal2D element area"] = A

    # 获取 CZM 失效单元列表
    deactivated_elements = Int64[]
    if case.opt.czm_enabled && haskey(case.multi_spme_layout, "czm_mesh")
        deactivated_elements = get_deactivated_elements(case.multi_spme_layout["czm_mesh"])
    end

    variables, I_e, Vc = solve_branch_currents_newton(
        case, variables, yt_representative, t, I_total, A, Te_prev;
        deactivated_elements=deactivated_elements
    )

    # 3. 并行求解每个单元
    M_blocks = Vector{Matrix{Float64}}(undef, ne)
    K_blocks = Vector{Matrix{Float64}}(undef, ne)
    F_matrix = zeros(n_chem, ne)
    vars_elems = Vector{Dict}(undef, ne)

    Threads.@threads for e in 1:ne
        yt_e = state.Y_chem[:, e]  # 直接取列，无需切片计算
        Me, Ke, Fe, vars_e = SPMe_element(case, yt_e, t, e; I_e=I_e[e], T_e=Te_prev[e], jacobi=jacobi)
        M_blocks[e] = Me
        K_blocks[e] = Ke
        F_matrix[:, e] = vec(Fe)
        vars_elems[e] = vars_e
    end

    # 4. 计算热源（使用正确的函数签名）
    if case.opt.czm_enabled && haskey(case.multi_spme_layout, "czm_mesh")
        czm_mesh = case.multi_spme_layout["czm_mesh"]
        variables = compute_heat_sources_with_czm(case, variables, vars_elems, I_e, Te_prev, A, czm_mesh, mesh_th)
    else
        variables = compute_heat_sources(case, variables, vars_elems, I_e, Te_prev, A; per_element_spme=true)
    end

    # 5. 装配热学矩阵
    MT, KT, FT = ThermalDistributed2D(case, variables)
    MT, KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)

    # 6. 更新 variables
    variables["cell voltage"] = Vc
    variables["time"] = t
    variables["temperature"] = thermal2D_volume_average_temperature(mesh_th, T_nodes)
    variables["T_nodes"] = T_nodes
    variables["thermal2D element current"] = I_e
    variables["thermal2D element voltages"] = [vars_elems[e]["cell voltage"] for e in 1:ne]

    # 保存辅助变量
    for e in 1:ne
        vars_e = vars_elems[e]
        variables["thermal2D eta_n_e"][e] = vars_e["negative electrode overpotential"][1]
        variables["thermal2D eta_p_e"][e] = vars_e["positive electrode overpotential"][end]
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        variables["thermal2D dUdT_n_e"][e] = case.param.NE.dUdT(cn_surf_e)[1]
        variables["thermal2D dUdT_p_e"][e] = case.param.PE.dUdT(cp_surf_e)[1]
        csn_data = vars_e["negative particle lithium concentration"]
        csp_data = vars_e["positive particle lithium concentration"]
        variables["thermal2D element soc_n"][e] = mean(vec(csn_data))
        variables["thermal2D element soc_p"][e] = mean(vec(csp_data))
    end

    return M_blocks, K_blocks, F_matrix, MT, KT, FT, variables
end

# 辅助函数：计算单元面积
function compute_element_areas(mesh_th)
    ne = size(mesh_th.element, 1)
    A = zeros(Float64, ne)
    ngs = length(mesh_th.gs.detJ)
    for g in 1:ngs
        e = mesh_th.gs.ele[g]
        A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    return A
end
```

### 6. 主循环适配（Solve.jl）

```julia
# 在主求解循环中
if should_use_multi_spme(case)
    # 使用新的分块状态结构
    state = ModelInitialisation(case; initial_soc_distribution=initial_soc_distribution)

    while t < t_end
        state, vars = solve_step_multi_spme(state, case, t, dt, θ)

        # 检查截止电压（使用 vars 中的电压）
        V = vars["cell voltage"]
        if V <= v_lower || V >= v_upper
            break
        end

        # 记录历史
        push!(time_hist, t)
        push!(voltage_hist, V)
        # ...

        t += dt
    end

    # 返回结果时转换格式（保持与现有结果结构兼容）
    return convert_to_result_format(state, vars_hist)
else
    # 原有单SPMe求解逻辑不变
    y0 = ModelInitialisation(case)
    # ...
end
```

### 7. 删除的内容

| 文件 | 删除内容 |
|------|----------|
| `Initialisation.jl` | `ModelInitialisation_MultiSPMe` 函数 |
| `Initialisation.jl` | `MultiSPMe_extract_element_state` |
| `Initialisation.jl` | `MultiSPMe_get_thermal_dofs` |
| `Initialisation.jl` | `MultiSPMe_update_state` |
| `SetCase.jl` | `multi_spme_layout` 中的 `ne`, `n_chem`, `nT`, `n_total`, `chem_range`, `thermal_range` |

**保留** `multi_spme_layout` 中的几何/求解器信息：
- `czm_mesh`
- `czm_element_map`
- `interface_pairs`
- `element_layer`
- `is_inner_layer`
- `inner_nodes`, `outer_nodes`
- `thermal_variables`, `thermal_update_fn`, `thermal_record`
- `polar_mesh_data`

## 改动文件清单

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `Initialisation.jl` | 重构 | 新增 `MultiSPMeState`，合并函数，删除辅助函数 |
| `Solve.jl` | 重构 | 新增 `solve_step_multi_spme`，修改主循环分发逻辑 |
| `SPMe.jl` | 不变 | `SPMe_element` 接口不变（接受向量） |
| `CycleSolver.jl` | 适配 | 适配新状态结构 |
| `CycleData.jl` | 适配 | 适配新状态结构 |
| `PostProcessing.jl` | 适配 | 适配新状态结构 |
| `SetCase.jl` | 清理 | 删除 `multi_spme_layout` 中的可计算字段 |

## 迁移计划

### 阶段 1：添加新结构（不破坏现有代码）
1. 添加 `MultiSPMeState` 结构体
2. 添加辅助函数 `get_ne`, `get_nT`, `get_n_chem`, `compute_element_temperatures`
3. 添加 `CallModel_MultiSPMe_blocked` 函数
4. 添加 `solve_step_multi_spme` 函数

### 阶段 2：修改初始化
1. 修改 `ModelInitialisation` 支持返回 `MultiSPMeState`
2. 保留 `ModelInitialisation_MultiSPMe` 作为别名（临时兼容）

### 阶段 3：修改求解器
1. 修改 `Solve.jl` 主循环，根据 `should_use_multi_spme` 分发到新函数
2. 修改 `CycleSolver.jl` 适配新状态结构

### 阶段 4：清理
1. 删除 `ModelInitialisation_MultiSPMe`
2. 删除 `MultiSPMe_extract_element_state`, `MultiSPMe_get_thermal_dofs`, `MultiSPMe_update_state`
3. 清理 `multi_spme_layout` 中的可计算字段

## 向后兼容性

- 单SPMe模式（SPM, SPMe without distributed2D）保持原有向量形式
- API 变化：`ModelInitialisation` 返回类型可能改变（新增关键字参数 `initial_soc_distribution`）
- `multi_spme_layout` 中删除的字段可通过辅助函数实时计算

## 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| 时间积分器改动引入数值误差 | 保持相同的时间离散格式（Crank-Nicolson），只改变数据结构 |
| 并行求解的线程安全问题 | 确保每个线程操作独立的内存区域，使用 `similar` 创建新数组 |
| 现有代码依赖 `multi_spme_layout` | 逐文件检查并替换为辅助函数调用，添加 `to_vector` 转换函数 |
| CZM 集成问题 | 保持 `multi_spme_layout` 中的 CZM 相关字段不变 |
| 调试日志丢失 | 在新函数中保留相同的 NaN 检测和温度范围警告逻辑 |

## 验收标准

1. 所有现有测试通过
2. 多SPMe仿真的数值结果与重构前一致（相对误差 < 1e-10）
3. `Initialisation.jl` 代码行数减少 ≥ 80 行（从 ~240 行减少到 ~160 行）
4. 初始化函数数量从 2 个减少到 1 个
5. `multi_spme_layout` 不再存储 `ne`, `n_chem`, `nT`, `chem_range`, `thermal_range`
6. `MultiSPMeState` 结构体正确定义且辅助函数可用

### 测试文件

| 文件 | 说明 |
|------|------|
| `example/SPMe_Thermal_example.jl` | 电化学-热耦合基础测试 |
| `example/czm_cycle_example.jl` | CZM 循环测试 |
| `example/testexample.jl` | 全耦合测试 |
| `example/热模块验证/thermal_verify.jl` | 热模型验证 |
