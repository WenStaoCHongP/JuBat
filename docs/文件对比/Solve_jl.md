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
- **行数**: 446
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `Solve(case::Case; initial_state, return_final_state, thermal_variables, thermal_update_fn, thermal_record, polar_mesh_data)` | 纯步进器入口（CallModel已迁出） |
  | 417 | `RecordMatrix!(...)` | 稀疏矩阵记录（未变更） |
  | 427 | `ErrorEstimation(...)` | 误差估计（未变更） |

## 变更详情

### 已迁出函数 (→ CallModel.jl)

#### `CallModel_MultiSPMe` 和 `CallModel`
- **原位置**: Solve.jl 第 530-850 行
- **新位置**: `CallModel.jl` 第 1-189 行（整体迁出，不做逻辑拆分或函数重命名）
- **迁出原因**: Solve.jl 保留纯步进器职责，模型调度逻辑独立维护

### 修改函数

#### `Solve(case::Case; ...)`
- **签名变更**: 新增关键字参数 `initial_state`, `return_final_state`, `thermal_variables`, `thermal_update_fn`, `thermal_record`, `polar_mesh_data`
- **新增功能**:
  - 纯热模式 `case.opt.model == "thermal"` 的独立求解路径
  - 多 SPMe 模式检测（`case.opt.per_element_spme`）
  - 外部状态恢复（用于循环仿真的相位衔接）
  - 多 SPMe 初始化路由（`ModelInitialisation_MultiSPMe`）
  - 热场初始化：从状态向量或外部传入恢复 `T_nodes`
  - 分布式热模型在主循环中的温度提取与回写
  - 单元级截止电压追踪（`first_cutoff_detected`, `cutoff_time`, `cutoff_element`）
  - 终止原因记录（`termination_reason`）
  - 结果中大量新增热学后处理字段（分层热源 W/m3、节点温度 K、单元温度 K）
  - `return_final_state` 支持：将最终状态打包为 Dict 返回，供循环求解器使用
  - 时间步最大步数限制（multi-SPMe: 50000, 单 SPMe: 100000）
  - 耗时统计（`timing_totals`），用于识别性能瓶颈
  - 使用 `case.layout`（`MultiSPMeLayout`）替代 `multi_spme_layout` Dict

### 删除函数
无（`CallModel` 和 `CallModel_MultiSPMe` 迁出到 `CallModel.jl`，非删除）

## 依赖关系

### 该文件依赖哪些其他文件
- `CallModel.jl` — 调用 `CallModel()` 和 `CallModel_MultiSPMe()`
- `Initialisation.jl` — 调用 `ModelInitialisation` 和 `ModelInitialisation_MultiSPMe`
- `CouplingState.jl` — 使用 `MultiSPMeLayout`（通过 `case.layout`）
- `Variables.jl` — 调用 `StandardVariables()`, `Variable_update!()`
- `PostProcessing.jl` — 调用 `PostProcessing()`

### 哪些文件依赖该文件
- `JuBat.jl` — export `Solve`
- `CycleSolver.jl` / `CycleData.jl` — 调用 `Solve(case; initial_state=..., return_final_state=true)`
- `example/` 下的所有示例脚本

### 新增的外部依赖
- 无新的外部包依赖

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系
Solve.jl 现在是耦合架构的**纯步进器**（模型调度已迁出到 CallModel.jl）：
- **时间推进**: 自适应步长控制、误差估计、截止电压检测
- **状态管理**: 外部状态恢复、热场携带、最终状态导出
- **后处理**: 热源物理单位转换、节点/单元温度输出

### 哪些变更是耦合相关的
- 多 SPMe 模式检测和初始化路由
- 热场温度提取与回写
- 截止电压追踪逻辑
- 热源后处理和物理单位转换
- 使用 `case.layout` 替代 Dict 访问

### 哪些变更是独立的
- 时间步安全限制
- `RecordMatrix!` 和 `ErrorEstimation` 未变更
- 耗时统计（通用基础设施）

## 后续变更 (2026-04-01)

- **移除调试日志重定向逻辑**：删除了文件开头的 `debug_to_file` 相关代码（约 35 行），包括 `output/debug_TIMESTAMP.log` 文件创建、`stdout`/`stderr` 重定向逻辑。`Solve` 函数不再自动创建调试日志文件。
- **移除初始状态 NaN/Inf 检查**：删除了约 5 行对 `initial_state` 中 NaN/Inf 值的检查代码。
- **移除单元温度跟踪调试代码**：删除了约 20 行用于调试目的的逐单元温度跟踪打印逻辑。
- **行数变化**：从约 851 行减少至约 710 行。
- **核心功能不变**：`Solve` 主函数仍处理所有模型类型（SPM/SPMe/P2D），`CallModel_MultiSPMe` 和 `CallModel` 调度逻辑未受影响。

## 后续变更 (2026-04-07)

- **P1 CallModel 拆分**: `CallModel_MultiSPMe` (~170行) + `CallModel` (~80行) 整体迁出到新文件 `CallModel.jl`
- **消除 thermal_extras Dict**: 移除 `case.thermal_extras::Dict{String,Any}` 字段及其所有访问
- **使用 typed 字段**: `case.layout` 替代 Dict 键值访问，`case.layout === nothing` 替代 `isempty` 检查
- **行数变化**: 从约 710 行减少至 446 行
- **职责变化**: Solve.jl 现在是纯步进器，模型调度逻辑在 CallModel.jl 中
