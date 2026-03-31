# MultiSPMe 状态结构与求解器重构设计

## 概述

重构 `Initialisation.jl` 和 `Solve.jl`，简化代码并消除重复逻辑。

**范围限定**：本重构针对 SPMe + distributed2D 热模型的多SPMe模式。SPM、单SPMe、P2D 模式保持原有向量形式不变。

## 动机

### 当前问题

#### Initialisation.jl 问题
1. **代码重复**：`ModelInitialisation` 和 `ModelInitialisation_MultiSPMe` 有大量重复逻辑
2. **状态结构复杂**：一维向量 `[yt_chem[1]; ...; yt_chem[ne]; T_nodes]` 需要复杂的索引管理
3. **布局信息冗余**：`multi_spme_layout` 存储的可计算信息（`ne`, `n_chem`, `nT`, `chem_range`, `thermal_range`）
4. **辅助函数过多**：`MultiSPMe_extract_element_state`, `MultiSPMe_get_thermal_dofs`, `MultiSPMe_update_state` 仅用于状态切片

#### Solve.jl 问题
1. **纯热模型分支**：`if case.opt.model == "thermal"` (行38-106) 嵌入主函数，增加复杂度
2. **theta 系数计算重复**：在 thermal 分支 (行46-54) 和 electrochemical 分支 (行161-169) 各赋值一次
3. **CallModel 分发复杂**：`CallModel_MultiSPMe` 是独立函数但只有一个调用点
4. **distributed2D 热处理重复**：`CallModel` 中的 distributed2D 分支 (行766-816) 与 `CallModel_MultiSPMe` 逻辑重复
5. **元素面积计算重复**：在 `CallModel_MultiSPMe` (行576-582) 和 `CallModel` distributed2D 分支 (行774-780) 重复
6. **元素均温计算重复**：在 `CallModel_MultiSPMe` (行585-589) 和 `CallModel` distributed2D 分支 (行784-788) 重复
7. **调试代码分散**：NaN 检查、温度异常检查等分散在多处，增加主循环噪声

### 目标

1. 合并初始化函数为一个 `ModelInitialisation`（已完成）
2. 电化学状态改为矩阵 `Y_chem::Matrix{Float64}` (n_chem × ne)（已完成）
3. 热场状态独立为 `T_nodes::Vector{Float64}` (nT)（已完成）
4. 删除 `multi_spme_layout` 中的可计算信息（已完成）
5. 时间积分器支持分块状态结构
6. **提取纯热模型分支为独立函数 `solve_thermal_only()`**
7. **将 `CallModel_MultiSPMe` 内联到 `CallModel`**
8. **theta 赋值合并到 Solve 函数顶部一处**
9. **面积计算统一用 `jellyroll_element_properties`**
10. **调试代码提取为辅助函数**

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

    # 1. 调用 CallModel 获取分块矩阵
    M_blocks, K_blocks, F_matrix, MT, KT, FT, vars = CallModel_multi_spme_inner(state, case, t)

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

### 5. Solve.jl 重构

#### 5.1 theta 合并到 Solve 函数顶部

**当前问题**：theta 在 thermal 分支 (行46-54) 和 electrochemical 分支 (行161-169) 各赋值一次。

**解决方案**：将 theta 赋值移到 Solve 函数开头（thermal 分支判断之前），thermal-only 的 early return 和 electrochemical 循环共享同一个 theta。

```julia
function Solve(case::Case; initial_state=nothing, return_final_state=false)
    # ... 日志初始化代码不变 ...

    result = nothing
    try
        # theta 时间离散系数（thermal 和 electrochemical 共用）
        if case.opt.solveType == "Crank-Nicolson"
            theta = 0.5
        elseif case.opt.solveType == "forward"
            theta = 0.0
        elseif case.opt.solveType == "backward"
            theta = 1.0
        else
            error("Error: $(case.opt.solveType) difference scheme has not been implemented!")
        end

        if case.opt.model == "thermal"
            return solve_thermal_only(case; theta=theta)
        end

        # ... electrochemical 求解逻辑（不再需要重复赋值 theta）...
```

**不提取 `get_theta()` 函数**：theta 赋值逻辑足够简单（3 个分支），提取为函数反而增加不必要的调用层。合并到函数顶部一处赋值即可自然消除重复。

#### 5.2 纯热模型分支提取为 solve_thermal_only()

**当前问题**：行38-106 的 thermal-only 逻辑（约70行）内嵌在 Solve 函数中，增加函数复杂度。

**解决方案**：提取为独立函数。

```julia
"""
    solve_thermal_only(case::Case; theta::Float64)

纯热模型求解（case.opt.model == "thermal"）。
从 Solve 中提取，用于分离关注点。
"""
function solve_thermal_only(case::Case; theta::Float64)
    t0 = case.opt.time[1] / case.param.scale.t0
    t_end = case.opt.time[end] / case.param.scale.t0
    dt = case.opt.dt[1] / case.param.scale.t0
    vars = case.multi_spme_layout["thermal_variables"]
    update_fn = case.multi_spme_layout["thermal_update_fn"]
    record = case.multi_spme_layout["thermal_record"]

    mesh = case.mesh["thermal2D"]
    nnode = mesh.nlen
    T_nodes = vars["T_nodes"]

    times = collect(range(t0, step=dt, length=Int(cld(t_end - t0, dt)) + 1))
    T_hist = record ? zeros(Float64, nnode, length(times)) : zeros(Float64, 0, 0)
    record && (T_hist[:, 1] .= T_nodes)

    update_fn !== nothing && update_fn(t0, vars)
    if case.opt.thermalmodel == "ring2D_polar"
        mesh_data = case.multi_spme_layout["polar_mesh_data"]
        MT_old, KT_old, FT_old = ThermalPolar2D_Ring(case, vars, mesh_data)
    elseif case.opt.thermalmodel == "ring2D"
        MT_old, KT_old, FT_old = ThermalDistributed2D_Ring(case, vars)
        KT_old, FT_old = ThermalRing2D_BC(KT_old, FT_old, case, vars["thermal2D outer_nodes"], t0)
    else
        MT_old, KT_old, FT_old = ThermalDistributed2D(case, vars)
        KT_old, FT_old = ThermalDistributed2D_BC(KT_old, FT_old, case, t0)
    end

    for step in 1:(length(times) - 1)
        t = times[step + 1]
        update_fn !== nothing && update_fn(t, vars)

        if case.opt.thermalmodel == "ring2D_polar"
            mesh_data = case.multi_spme_layout["polar_mesh_data"]
            MT_new, KT_new, FT_new = ThermalPolar2D_Ring(case, vars, mesh_data)
        elseif case.opt.thermalmodel == "ring2D"
            MT_new, KT_new, FT_new = ThermalDistributed2D_Ring(case, vars)
            KT_new, FT_new = ThermalRing2D_BC(KT_new, FT_new, case, vars["thermal2D outer_nodes"], t)
        else
            MT_new, KT_new, FT_new = ThermalDistributed2D(case, vars)
            KT_new, FT_new = ThermalDistributed2D_BC(KT_new, FT_new, case, t)
        end

        A = MT_new - theta * dt * KT_new
        rhs = (MT_new + (1.0 - theta) * dt * KT_old) * T_nodes +
              dt * (theta * FT_new + (1.0 - theta) * FT_old)
        T_nodes = A \ rhs
        vars["T_nodes"] = T_nodes

        record && (T_hist[:, step + 1] .= T_nodes)

        MT_old = MT_new
        KT_old = KT_new
        FT_old = FT_new
    end

    return (time = times, T_nodes = T_nodes, T_hist = T_hist)
end
```

Solve 中简化为：
```julia
if case.opt.model == "thermal"
    return solve_thermal_only(case; theta=theta)
end
```

#### 5.3 CallModel 内联 CallModel_MultiSPMe

**当前问题**：`CallModel_MultiSPMe` (行530-710) 是独立函数，只有一个调用点在 `CallModel` 中（行744）。中间还有布局自动推断逻辑（行717-738）。

**解决方案**：将多SPMe逻辑直接内联到 `CallModel` 的 `multi_spme_enabled` 分支，删除独立函数和布局自动推断逻辑。

```julia
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    multi_spme_enabled = (case.opt.model == "SPMe" && case.opt.thermalmodel == "distributed2D" && !isempty(case.multi_spme_layout))

    if multi_spme_enabled
        # === 多SPMe 模式（原 CallModel_MultiSPMe 逻辑内联） ===
        # 使用 jellyroll_element_properties 计算面积
        mesh_th = case.mesh["thermal2D"]
        areas, layer_weights = jellyroll_element_properties(mesh_th, case.param)
        variables["thermal2D element area"] = areas

        # 提取热场
        T_nodes = MultiSPMe_get_thermal_dofs(yt, case)
        Te_prev = compute_element_temperatures(mesh_th, T_nodes)

        # 分流求解
        I_total = case.opt.Current(t * case.param.scale.t0) / case.param.scale.I_typ
        variables = StandardVariables(case, 1)
        variables["T_nodes"] = T_nodes
        variables["thermal2D element area"] = areas
        variables["cell current"] = I_total

        # ... 代表性状态 + 分流求解 ...
        # ... 并行 SPMe_element 求解 ...
        # ... 热源计算 + 热矩阵装配 ...

        return M, K, F, variables, y_phi
    end

    # 原有逻辑（单SPMe模式）
    if case.opt.model == "SPM"
        M, K, F, variables = SPM(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "SPMe"
        M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "P2D"
        M, K, F, variables, y_phi = P2D(case, yt, t, jacobi=jacobi)
    elseif case.opt.model == "sP2D"
        M, K, F, variables, y_phi = sP2D(case, yt, t, jacobi=jacobi)
    else
        error("Error: $(case.opt.model) model has not been implemented!")
    end

    if case.opt.thermalmodel == "lumped"
        MT, FT = ThermalLumped(case, variables)
        M = blockdiag(M, sparse(MT))
        K = blockdiag(K, sparse(zeros(1,1)))
        F = [F; FT]
    elseif case.opt.thermalmodel == "distributed2D"
        # 单 SPMe + distributed2D：使用 jellyroll_element_properties
        mesh_th = case.mesh["thermal2D"]
        areas, layer_weights = jellyroll_element_properties(mesh_th, case.param)
        variables["thermal2D element area"] = areas

        nT = mesh_th.nlen
        variables["T_nodes"] = yt[(end - nT + 1):end]
        Te_prev = compute_element_temperatures(mesh_th, variables["T_nodes"])

        # 分流求解 + 热源计算 + 热矩阵装配
        # ...
    end

    return M, K, F, variables, y_phi
end
```

#### 5.4 面积计算统一用 jellyroll_element_properties

**当前重复**：
- `CallModel_MultiSPMe` 行576-583：手动循环 `A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]`
- `CallModel` distributed2D 分支行774-779：完全相同的循环

**解决方案**：两处统一替换为已存在于 `Jellyrollmodel.jl` 中的函数：

```julia
# 替代手动循环
areas, layer_weights = jellyroll_element_properties(mesh_th, case.param)
variables["thermal2D element area"] = areas
```

**不新增 `compute_element_areas` 函数**：`jellyroll_element_properties` 已提供完全相同的功能，新增函数属于冗余。

#### 5.5 调试代码提取为辅助函数

**当前分散的调试代码**：
- Solve 行192-195：初始状态 NaN 检查
- Solve 行244-248：初始求解步骤异常检查
- CallModel_MultiSPMe 行559-567：温度场异常检查
- CallModel_MultiSPMe 行590-595：温度场 NaN 检查

**提取为两个辅助函数**：

```julia
"""
    check_state_validity(yt; context="") -> Int

检查状态向量中 NaN/Inf 的数量，打印警告。返回异常值数量。
"""
function check_state_validity(yt; context="")
    nan_count = sum(.!isfinite.(yt))
    if nan_count > 0
        prefix = isempty(context) ? "" : "[$context] "
        @warn "$(prefix)状态向量包含 $nan_count 个 NaN/Inf，长度 $(length(yt))"
    end
    return nan_count
end

"""
    check_temperature_field(T_nodes; context="")

检查温度场是否包含异常值（NaN/Inf、过大值、过大偏差），打印警告。
"""
function check_temperature_field(T_nodes; context="")
    nan_count = sum(.!isfinite.(T_nodes))
    abnormal = sum(abs.(T_nodes) .> 10.0)
    large_dev = sum(abs.(T_nodes .- 1.0) .> 5.0)
    if nan_count > 0 || abnormal > 0 || large_dev > 0
        T_min, T_max = extrema(T_nodes)
        prefix = isempty(context) ? "" : "[$context] "
        @warn "$(prefix)温度场异常" range=(T_min,T_max) nan=nan_count abnormal=abnormal large_dev=large_dev
    end
end
```

使用方式：
```julia
# 替代行192-195
check_state_validity(y0; context="初始化")

# 替代行559-567
check_temperature_field(T_nodes; context="CallModel")

# 替代行590-595
check_temperature_field(Te_prev; context="单元均温")
```

#### 5.6 简化后的代码结构

```
Solve.jl 重构后:
├── check_state_validity()              # 新增：状态向量检查
├── check_temperature_field()           # 新增：温度场检查
├── solve_thermal_only()                # 新增：纯热模型求解（从 Solve 提取）
├── Solve()                             # 简化：theta 统一赋值，thermal 分支委托
│   ├── theta 赋值（一处，thermal/electro 共用）
│   ├── thermal-only → solve_thermal_only()
│   └── electrochemical 主循环
├── CallModel()                         # 整合：多SPMe 逻辑内联，面积用 jellyroll_element_properties
│   ├── 多SPMe 分支（内联，原 CallModel_MultiSPMe）
│   ├── SPM/SPMe/P2D 分支
│   └── 热模型分支（面积统一）
├── RecordMatrix!()                     # 不变
└── ErrorEstimation()                   # 不变
```

### 6. 删除的内容

| 文件 | 删除内容 |
|------|----------|
| `Initialisation.jl` | `ModelInitialisation_MultiSPMe` 函数 |
| `Initialisation.jl` | `MultiSPMe_extract_element_state` |
| `Initialisation.jl` | `MultiSPMe_get_thermal_dofs` |
| `Initialisation.jl` | `MultiSPMe_update_state` |
| `Solve.jl` | `CallModel_MultiSPMe` 独立函数（内联到 CallModel） |
| `Solve.jl` | theta 重复赋值（行161-169） |
| `Solve.jl` | 两处手动面积计算循环（改用 jellyroll_element_properties） |
| `SetCase.jl` | `multi_spme_layout` 中的 `ne`, `n_chem`, `nT`, `n_total`, `chem_range`, `thermal_range` |

**不新增** `get_theta()` 和 `compute_element_areas()` 函数。

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
| `Solve.jl` | 重构 | theta 合并、纯热分支提取、CallModel 内联、面积统一、调试分离 |
| `SPMe.jl` | 不变 | `SPMe_element` 接口不变（接受向量） |
| `CycleSolver.jl` | 适配 | 适配新状态结构 |
| `CycleData.jl` | 适配 | 适配新状态结构 |
| `PostProcessing.jl` | 适配 | 适配新状态结构 |
| `SetCase.jl` | 清理 | 删除 `multi_spme_layout` 中的可计算字段 |

## 实施顺序

每项优化独立执行并测试，按依赖顺序排列：

### 阶段 1：Initialisation.jl（已完成）

1. ~~添加 `MultiSPMeState` 结构体~~
2. ~~重写 `ModelInitialisation` 函数（合并逻辑）~~
3. ~~删除旧函数~~

### 阶段 2：Solve.jl（每项独立）

| 步骤 | 任务 | 风险 | 测试 |
|------|------|------|------|
| 4 | theta 合并到 Solve 顶部 | 低 | 运行 SPMe_Thermal_example.jl |
| 5 | CallModel 内联 CallModel_MultiSPMe | 中 | 运行全耦合 testexample.jl |
| 6 | 面积计算统一用 jellyroll_element_properties | 低 | 对比热源计算结果 |
| 7 | 使用 Variables.jl 中的正确变量名 | 低 | 检查输出变量 |
| 8 | 纯热模型分支提取为 solve_thermal_only() | 低 | 运行 thermal_verify.jl |
| 8.5 | 调试代码提取为辅助函数 | 低 | 运行任意仿真确认日志输出 |

### 阶段 3：调用点更新

1. 更新 `CycleSolver.jl`
2. 更新 `PostProcessing.jl`
3. 更新 `ThermalDistributed.jl`（如有需要）

### 阶段 4：验证

1. 运行所有测试
2. 对比重构前后数值结果

## 向后兼容性

- 单SPMe模式（SPM, SPMe without distributed2D）保持原有向量形式
- API 变化：`ModelInitialisation` 返回类型可能改变（新增关键字参数 `initial_soc_distribution`）
- `multi_spme_layout` 中删除的字段可通过辅助函数实时计算

## 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| CallModel 内联可能遗漏逻辑 | 逐步内联，每步对比数值结果 |
| 时间积分器改动引入数值误差 | 保持相同的时间离散格式（Crank-Nicolson），只改变数据结构 |
| 并行求解的线程安全问题 | 确保每个线程操作独立的内存区域，使用 `similar` 创建新数组 |
| 现有代码依赖 `multi_spme_layout` | 逐文件检查并替换为辅助函数调用 |
| CZM 集成问题 | 保持 `multi_spme_layout` 中的 CZM 相关字段不变 |
| 调试日志丢失 | 辅助函数保留相同的 NaN 检测和温度范围警告逻辑 |

## 验收标准

1. 所有现有测试通过
2. 多SPMe仿真的数值结果与重构前一致（相对误差 < 1e-10）
3. `Initialisation.jl` 代码行数减少 >= 80 行
4. `Solve.jl` 代码行数减少 >= 150 行
5. 初始化函数数量从 2 个减少到 1 个
6. `multi_spme_layout` 不再存储 `ne`, `n_chem`, `nT`, `chem_range`, `thermal_range`
7. `MultiSPMeState` 结构体正确定义且辅助函数可用
8. 不新增 `get_theta()` 函数（theta 内联合并）
9. 不新增 `compute_element_areas()` 函数（使用 `jellyroll_element_properties`）
10. 纯热模型逻辑在独立函数 `solve_thermal_only()` 中
11. 调试检查通过 `check_state_validity()` 和 `check_temperature_field()` 执行

### 测试文件

| 文件 | 说明 |
|------|------|
| `example/SPMe_Thermal_example.jl` | 电化学-热耦合基础测试 |
| `example/czm_cycle_example.jl` | CZM 循环测试 |
| `example/testexample.jl` | 全耦合测试 |
| `example/热模块验证/thermal_verify.jl` | 热模型验证 |
