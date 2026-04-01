# Solve.jl

## 文件状态
修改

## main分支
- **行数**: 155
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `Solve(case::Case)` | 主求解器入口 |
  | 103 | `CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)` | 根据模型类型调用对应模型（SPM/SPMe/P2D/sP2D） |
  | 126 | `RecordMatrix!(case::Case, M::SparseMatrixCSC, K::SparseMatrixCSC)` | 稀疏矩阵记录 |
  | 136 | `ErrorEstimation(case::Case, y_old, y_new, coeff::Float64)` | 自适应时间步误差估计 |

## Parameters_Design分支
- **行数**: 850
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `Solve(case::Case; initial_state=nothing, return_final_state=false)` | 主求解器入口，增加状态传递与循环衔接支持 |
  | 530 | `CallModel_MultiSPMe(case::Case, yt, t; jacobi)` | 多SPMe并行架构 CallModel，每个热单元独立求解 |
  | 711 | `CallModel(case::Case, yt, t; jacobi)` | 原有调度入口，新增 multi-SPMe 自动路由 |
  | 821 | `RecordMatrix!(...)` | 稀疏矩阵记录（未变更） |
  | 831 | `ErrorEstimation(...)` | 误差估计（未变更） |

## 变更详情

### 新增函数

#### `CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)`
- **位置**: 第 530 行
- **功能**: 多 SPMe 并行架构的核心求解函数。工作流程：
  1. 解析全局状态向量，提取每个单元的电化学状态和热场
  2. 计算元素面积和均温
  3. 调用分流求解器 `solve_branch_currents_newton` 获取逐单元电流 I_e
  4. 并行（`Threads.@threads`）调用 `SPMe_element` 求解每个单元的电化学响应
  5. 计算逐单元热源（含 CZM 路径）
  6. 装配热学矩阵（`ThermalDistributed2D` + BC）
  7. 全局装配 `blockdiag(M_chem, MT)`，返回统一的 M/K/F

### 修改函数

#### `Solve(case::Case; ...)`
- **签名变更**: 新增关键字参数 `initial_state::Union{Dict{String,Any},Nothing}=nothing` 和 `return_final_state::Bool=false`
- **新增功能**:
  - 调试日志重定向到 `output/debug_TIMESTAMP.log` 文件
  - 纯热模式 `case.opt.model == "thermal"` 的独立求解路径
  - 多 SPMe 模式检测（`model == "SPMe" && thermalmodel == "distributed2D"`）
  - 外部状态恢复（用于循环仿真的相位衔接）
  - 多 SPMe 初始化路由（`ModelInitialisation_MultiSPMe`）
  - 热场初始化：从状态向量或外部传入恢复 `T_nodes`
  - Collector-seeded `layer_weights` 计算
  - 分布式热模型在主循环中的温度提取与回写
  - 单元级截止电压追踪（`first_cutoff_detected`, `cutoff_time`, `cutoff_element`）
  - 终止原因记录（`termination_reason`）
  - 结果中大量新增热学后处理字段（分层热源 W/m3、节点温度 K、单元温度 K）
  - `return_final_state` 支持：将最终状态打包为 Dict 返回，供循环求解器使用
  - 时间步最大步数限制（multi-SPMe: 50000, 单 SPMe: 100000）

#### `CallModel(case::Case, yt, t; jacobi)`
- **新增逻辑**:
  - 自动检测是否应启用多 SPMe 模式（`should_use_multi_spme`）
  - 当 `multi_spme_layout` 为空但状态向量长度匹配时自动初始化布局
  - 多 SPMe 模式下直接路由到 `CallModel_MultiSPMe`
  - 新增 `distributed2D` 热模型分支（非 multi-SPMe 的单 SPMe + distributed2D 路径）：
    - 从状态向量提取热自由度
    - 计算元素面积和均温
    - 调用非线性分流求解器
    - 计算热源并装配热学矩阵
    - 拼接到主系统矩阵

### 删除函数
无

## 依赖关系

### 该文件依赖哪些其他文件
- `Initialisation.jl` — 调用 `ModelInitialisation` 和 `ModelInitialisation_MultiSPMe`
- `SPMe.jl` — 通过 `CallModel` 间接调用 `SPMe()` 和 `SPMe_element()`
- `Parallelsolution.jl` — 调用 `solve_branch_currents_newton()`
- `ThermalDistributed.jl` — 调用 `ThermalDistributed2D()`, `ThermalDistributed2D_BC()`, `compute_heat_sources()`, `compute_heat_sources_with_czm()`, `thermal2D_volume_average_temperature()`, `element_nodal_mean()`
- `Jellyrollmodel.jl` — 调用 `jellyroll_element_properties()`（collector-seeded 模式）
- `Variables.jl` — 调用 `StandardVariables()`, `Variable_update!()`
- `Postprocess.jl` — 调用 `PostProcessing()`

### 哪些文件依赖该文件
- `JuBat.jl` — export `Solve`
- `CycleSolver.jl` / `CycleData.jl` — 调用 `Solve(case; initial_state=..., return_final_state=true)`
- `example/` 下的所有示例脚本

### 新增的外部依赖
- `Threads` — 多线程并行求解（`Threads.@threads`）
- 无其他新的外部包依赖

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系
Solve.jl 是整个耦合架构的**调度中枢**，所有耦合逻辑都在此文件中协调：
- **multi-SPMe 与 distributed2D 的耦合**: `CallModel_MultiSPMe` 中，先求解电化学得到逐单元热源，再装配热学矩阵，实现电-热双向耦合
- **CZM 耦合**: 在热源计算时区分 `compute_heat_sources_with_czm` 和 `compute_heat_sources` 两条路径
- **分流求解器耦合**: 通过 `solve_branch_currents_newton` 获取逐单元电流分布

### 哪些变更是耦合相关的
- `CallModel_MultiSPMe` 整个新增函数（第 530-694 行）— 这是 multi-SPMe + distributed2D 的核心耦合实现
- `Solve` 中的 `thermal-distributed` 初始化、温度提取与回写（第 210-304 行）
- `CallModel` 中的 `distributed2D` 分支（第 754-804 行）
- 截止电压追踪逻辑（第 342-380 行）
- 热源后处理和物理单位转换（第 403-457 行）

### 哪些变更是独立的
- 调试日志重定向机制（第 13-46 行）— 独立的基础设施改进
- 时间步安全限制（第 183-190 行）
- 初始状态 NaN 检查（第 201-207 行）
- `RecordMatrix!` 和 `ErrorEstimation` 未变更
