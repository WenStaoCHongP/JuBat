# CallModel.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 188
- 路径: `src/CallModel.jl`
- 来源: 从 `Solve.jl` 中整体迁出（commit 3ff3027）

### 主要函数列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)` | L1 | 多SPMe并行架构 CallModel，每个热单元独立求解 |
| `CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)` | L162 | 原有调度入口，自动路由到 MultiSPMe 或单模型 |

## 功能描述

本文件从 Solve.jl 中提取，负责模型调度（根据模型类型调用对应的物理模型）。Solve.jl 保留纯步进器职责。

### CallModel_MultiSPMe
多 SPMe 并行架构的核心求解函数（约 161 行），工作流程：
1. 验证 `case.layout` 前提条件（fail-fast）
2. 解析全局状态向量：`extract_element_state` 提取每个单元的电化学状态，`get_thermal_dofs` 提取热场
3. 计算元素面积和均温
4. 调用分流求解器 `solve_branch_currents` 获取逐单元电流 I_e
5. 并行（`Threads.@threads`）调用 `SPMe_element` 求解每个单元的电化学响应
6. 计算逐单元热源（含 CZM 路径分支）
7. 装配热学矩阵（`ThermalDistributed2D` + BC）
8. 全局装配 `blockdiag(M_chem, MT)`，返回统一的 M/K/F

### CallModel
模型调度入口（约 27 行）：
- 当 `case.opt.per_element_spme` 为 true 时，委托给 `CallModel_MultiSPMe`
- 否则根据 `case.opt.model` 调用 SPM/SPMe/P2D
- 支持 lumped 热模型耦合

## 依赖关系

### 该文件依赖哪些其他文件
- `CouplingState.jl` — 使用 `MultiSPMeLayout`（通过 `case.layout`）
- `Initialisation.jl` — 调用 `extract_element_state`, `get_thermal_dofs`
- `SPMe.jl` — 调用 `SPMe()`, `SPMe_element()`, `SPMe_variables()`
- `SPM.jl`, `P2D.jl` — 被 `CallModel` 调度
- `Parallelsolution.jl` — 调用 `solve_branch_currents()`
- `ThermalDistributed.jl` — 调用 `ThermalDistributed2D()`, `ThermalDistributed2D_BC()`, `compute_heat_sources()`, `compute_heat_sources_with_czm()`
- `Tools.jl` — 调用 `thermal2D_volume_average_temperature()`

### 哪些文件依赖该文件
- `Solve.jl` — 调用 `CallModel()` 和 `CallModel_MultiSPMe()`
- `CycleData.jl` — 调用 `CallModel()`（步后变量获取）

## 耦合分析

此文件是整个耦合架构的**模型调度层**：
- **multi-SPMe 与 distributed2D 耦合**: `CallModel_MultiSPMe` 中先求解电化学得到逐单元热源，再装配热学矩阵
- **CZM 耦合**: 热源计算时区分 `compute_heat_sources_with_czm` 和 `compute_heat_sources` 两条路径
- **分流求解器耦合**: 通过 `solve_branch_currents` 获取逐单元电流分布

## 设计原则

- 从 Solve.jl **整体搬家**，不做逻辑拆分或函数重命名
- CallModel.jl 独立于 Solve.jl，可单独维护和测试
- Solve.jl 保留纯步进器职责（时间推进、自适应步长、截止检测）
