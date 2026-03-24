# MultiSPMe 状态结构重构设计

## 概述

重构 `Initialisation.jl` 中的初始化逻辑，将电化学状态从一维向量改为矩阵形式，简化代码并消除重复。

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

新增 `MultiSPMeState` 结构体：

```julia
"""
多SPMe状态结构：电化学矩阵 + 热场向量
"""
struct MultiSPMeState
    Y_chem::Matrix{Float64}   # (n_chem × ne) 每列是一个单元的电化学状态
    T_nodes::Vector{Float64}  # (nT,) 热场节点温度
end

# 辅助函数
Base.length(s::MultiSPMeState) = length(s.Y_chem) + length(s.T_nodes)
n_elements(s::MultiSPMeState) = size(s.Y_chem, 2)
n_chem_dofs(s::MultiSPMeState) = size(s.Y_chem, 1)
```

**设计决策**：
- 矩阵每列对应一个热单元的电化学状态
- Julia 按列存储，访问 `Y_chem[:, e]` 高效
- 热场保持向量形式（全局耦合）

### 2. 初始化函数（合并后）

```julia
function ModelInitialisation(case::Case)
    # 获取基础电化学状态（单单元）
    y0_chem = init_chem_state(case)  # 提取为子函数

    if should_use_multi_spme(case)
        # 多SPMe模式：矩阵形式
        ne = get_ne(case)
        n_chem = length(y0_chem)
        nT = get_nT(case)

        Y_chem = repeat(y0_chem, 1, ne)  # 复制为矩阵
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
        # P2D 额外自由度...
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
```

### 4. 时间积分器

```julia
function solve_step_multi_spme(state::MultiSPMeState, case, t, dt, θ)
    ne = n_elements(state)
    n_chem = n_chem_dofs(state)

    # 1. 调用 CallModel_MultiSPMe 获取分块矩阵
    M_blocks, K_blocks, F_matrix, MT, KT, FT, vars = CallModel_MultiSPMe_blocked(state, case, t)
    # M_blocks: Vector{Matrix} (ne个)
    # F_matrix: Matrix (n_chem × ne)

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
function CallModel_MultiSPMe_blocked(state::MultiSPMeState, case::Case, t::Float64)
    ne = n_elements(state)
    T_nodes = state.T_nodes

    # 计算单元均温
    Te_prev = compute_element_temperatures(case.mesh["thermal2D"], T_nodes)

    # 分流求解
    I_e = solve_branch_currents(case, state, t, Te_prev)

    # 并行求解每个单元
    M_blocks = Vector{Matrix{Float64}}(undef, ne)
    K_blocks = Vector{Matrix{Float64}}(undef, ne)
    F_matrix = zeros(n_chem_dofs(state), ne)
    vars_elems = Vector{Dict}(undef, ne)

    Threads.@threads for e in 1:ne
        yt_e = state.Y_chem[:, e]  # 直接取列，无需切片计算
        Me, Ke, Fe, vars_e = SPMe_element(case, yt_e, t, e; I_e=I_e[e], T_e=Te_prev[e])
        M_blocks[e] = Me
        K_blocks[e] = Ke
        F_matrix[:, e] = vec(Fe)
        vars_elems[e] = vars_e
    end

    # 计算热源和热矩阵
    vars = compute_heat_sources(case, vars_elems, I_e, Te_prev)
    MT, KT, FT = ThermalDistributed2D(case, vars)
    MT, KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)

    return M_blocks, K_blocks, F_matrix, MT, KT, FT, vars
end
```

### 6. 删除的内容

| 文件 | 删除内容 |
|------|----------|
| `Initialisation.jl` | `ModelInitialisation_MultiSPMe` 函数 |
| `Initialisation.jl` | `MultiSPMe_extract_element_state` |
| `Initialisation.jl` | `MultiSPMe_get_thermal_dofs` |
| `Initialisation.jl` | `MultiSPMe_update_state` |
| `SetCase.jl` | `multi_spme_layout` 中的 `ne`, `n_chem`, `nT`, `chem_range`, `thermal_range` |

**保留** `multi_spme_layout` 中的几何信息：
- `czm_mesh`
- `interface_pairs`
- `element_layer`
- `is_inner_layer`
- `inner_nodes`, `outer_nodes`
- `thermal_variables`, `thermal_update_fn`, `thermal_record`

## 改动文件清单

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `SetCase.jl` | 新增 | `MultiSPMeState` 结构体 |
| `Initialisation.jl` | 重构 | 合并函数，删除辅助函数 |
| `Solve.jl` | 重构 | 新增 `solve_step_multi_spme`，修改主循环 |
| `SPMe.jl` | 不变 | `SPMe_element` 接口不变（接受向量） |
| `CycleSolver.jl` | 适配 | 适配新状态结构 |
| `CycleData.jl` | 适配 | 适配新状态结构 |
| `PostProcessing.jl` | 适配 | 适配新状态结构 |

## 向后兼容性

- 单SPMe模式（SPM, SPMe without distributed2D）保持原有向量形式
- API 变化：`ModelInitialisation` 返回类型可能改变
- 需要更新测试脚本和示例代码

## 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| 时间积分器改动引入数值误差 | 保持相同的时间离散格式（Crank-Nicolson），只改变数据结构 |
| 并行求解的线程安全问题 | 确保每个线程操作独立的内存区域 |
| 现有代码依赖 `multi_spme_layout` | 逐文件检查并替换为辅助函数调用 |

## 验收标准

1. 所有现有测试通过
2. 多SPMe仿真的数值结果与重构前一致（误差 < 1e-10）
3. 代码行数减少（删除重复代码）
4. 初始化函数数量从2个减少到1个
5. `multi_spme_layout` 不再存储可计算信息
