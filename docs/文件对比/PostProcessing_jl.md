# PostProcessing.jl

## 文件状态: 修改

## main分支
- 行数: 48
- 主要函数列表:
  - `PostProcessing(case, variables, v)` — 后处理函数，将无量纲变量转换为有量纲结果字典

## Parameters_Design分支
- 行数: 311
- 主要函数列表:
  - `PostProcessing(case, variables, v)` — 原后处理函数（扩展了distributed2D热输出）
  - `_phase_termination_symbol(phase_type, reason)` — **新增**：确定阶段终止类型符号
  - `_state_concentration_variance(case, y_state)` — **新增**：计算多SPMe浓度分布方差
  - `_postprocess_phase_result(case, phase_type, solve_result, initial_state, I_current, t_start, D_max_init, D_mean_init, czm_mesh)` — **新增**：阶段结果后处理
  - `_postprocess_cycle_result!(cycle_result, charge_result, discharge_result, rest1_result, rest2_result, czm_mesh)` — **新增**：循环结果汇总
  - `_append_cycle_result!(result, cycle, cycle_result; save_detailed)` — **新增**：追加循环结果
  - `_update_soh_and_capacity!(result, cycle, cycle_result, initial_capacity)` — **新增**：SOH更新
  - `_print_cycle_summary(cycle, cycle_result, current_soh)` — **新增**：打印循环摘要
  - `_check_cycle_termination(cycle, cycle_result, czm_mesh, current_soh, soh_threshold; verbose)` — **新增**：检查循环终止条件
  - `_print_cycling_summary(result, initial_capacity, soh_terminated)` — **新增**：打印最终循环摘要
  - `plot_cycling_results(result; save_path)` — **新增**：绘制循环结果图（4子图）

## 变更详情

### 修改函数

#### `PostProcessing(case, variables, v)` (L1-L70)
新增distributed2D热模型的后处理输出：

**新增输出项（distributed2D模式）**：
| 键 | 说明 | 量纲转换 |
|----|------|----------|
| `"thermal2D Q_rxn_NE [W/m3]"` | 负极反应热 | `* scale.q` |
| `"thermal2D Q_rev_NE [W/m3]"` | 负极可逆热 | `* scale.q` |
| `"thermal2D Q_ohm_s_NE [W/m3]"` | 负极固相欧姆热 | `* scale.q` |
| `"thermal2D Q_ohm_e_NE [W/m3]"` | 负极液相欧姆热 | `* scale.q` |
| `"thermal2D Q_SP [W/m3]"` | 隔膜欧姆热 | `* scale.q` |
| `"thermal2D Q_rxn_PE [W/m3]"` | 正极反应热 | `* scale.q` |
| `"thermal2D Q_rev_PE [W/m3]"` | 正极可逆热 | `* scale.q` |
| `"thermal2D Q_ohm_s_PE [W/m3]"` | 正极固相欧姆热 | `* scale.q` |
| `"thermal2D Q_ohm_e_PE [W/m3]"` | 正极液相欧姆热 | `* scale.q` |
| `"thermal2D Q_PCC [W/m3]"` | 正极集流体欧姆热 | `* scale.q` |
| `"thermal2D Q_NCC [W/m3]"` | 负极集流体欧姆热 | `* scale.q` |
| `"thermal2D temperature at nodes [K]"` | 节点温度场 | `* scale.T_ref` |

**新增输出项（lumped模式）**：
| 键 | 说明 |
|----|------|
| `"thermal lumped internal heat [W/m3]"` | 集总模型内热源 |

**注释掉的CZM输出**：CZM相关后处理代码存在但被注释掉（用三引号`"""`包裹），包括位移、损伤、牵引力、分离位移。这可能是因为CZM结果的量纲转换尚未完全验证。

### 新增函数

#### `_phase_termination_symbol(phase_type, reason)` (L72)
- 确定阶段终止类型符号：静置阶段返回 `:time`，电压截止返回 `:voltage`

#### `_state_concentration_variance(case, y_state)` (L79)
- 计算多SPMe场景下各单元的浓度分布方差
- 用于静置阶段评估锂扩散松弛程度
- 支持 multi_spme 和单 SPMe 两种模式

#### `_postprocess_phase_result(...)` (L104)
- 从solve_result提取阶段关键指标：持续时间、电压范围、容量、最高温度、损伤变化
- 计算容量 = abs(I_current) * duration / 3600
- 静置阶段额外计算浓度松弛百分比（`cs_relaxation_n/p`）

#### `_postprocess_cycle_result!(...)` (L165)
- 计算库伦效率 = 100 * discharge_capacity / charge_capacity
- 更新循环末损伤统计（max_D, mean_D, n_fractured）
- 取四个阶段的最高温度作为循环最高温度

#### `_append_cycle_result!(...)` (L178)
- 将CycleResult追加到CyclingResult的向量中
- 可选保存详细结果（`save_detailed`控制）

#### `_update_soh_and_capacity!(...)` (L193)
- 第一个循环的放电容量作为初始容量
- SOH = current_discharge_capacity / initial_capacity

#### `_print_cycle_summary(...)` (L202)
- 格式化输出单循环摘要：充电/放电容量、库伦效率、损伤、SOH

#### `_check_cycle_termination(...)` (L208)
- 终止条件1：SOH <= soh_threshold（且cycle > 1）
- 终止条件2：断裂CZM单元 > 50%总数
- 返回 (should_stop, soh_terminated) 元组

#### `_print_cycling_summary(...)` (L225)
- 格式化输出全部循环的最终摘要

#### `plot_cycling_results(result; save_path)` (L241)
- 绘制4个子图：容量衰减、损伤演化、库伦效率、温度历史
- 使用Plots.jl保存为PNG文件
- 返回组合图对象

### 删除函数

无删除。

## 依赖关系

### main分支依赖
- `src/Variables.jl` — Case类型、variables字典结构
- `src/SetParams.jl` — scale参数（量纲转换）

### Parameters_Design分支新增依赖
- `src/CzmSolve.jl` — `get_damage_statistics`（用于阶段/循环后处理）
- `src/Jellyrollmodel.jl`（间接） — 通过case.multi_spme_layout使用布局信息
- `Plots.jl` — `plot_cycling_results`绘图功能（外部依赖）

### 调用该文件的文件
- `src/JuBat.jl` — `include("PostProcessing.jl")`（L21）
- `src/CycleSolver.jl` — 调用所有 `_postprocess_*`、`_print_*`、`_check_*`、`_update_*`、`_append_*` 辅助函数
- `src/Solve.jl` — 调用 `PostProcessing` 主函数
- 用户脚本 — 调用 `plot_cycling_results` 绘制循环结果

## 耦合分析

### main分支
- PostProcessing仅作为纯量纲转换工具，将无量纲内部变量转换为物理单位输出
- 不参与耦合计算，是求解流程的最后一步

### Parameters_Design分支扩展
- **从后处理工具扩展为循环仿真基础设施**：新增的大量函数是CycleSolver的后端支撑
  - `_postprocess_phase_result`、`_postprocess_cycle_result!` 是循环求解器调用的核心后处理
  - `_check_cycle_termination` 实现了SOH和CZM断裂的终止判断逻辑
  - `_print_*` 系列函数提供循环过程的实时监控输出

- **与distributed2D热模型耦合**：新增了完整的热源后处理（12个热源分量），支持分层热源分析和验证

- **与CZM耦合**：通过 `get_damage_statistics` 获取损伤统计信息。CZM后处理代码已存在但被注释，表明CZM输出功能处于开发中。

- **与multi-SPMe耦合**：`_state_concentration_variance` 支持多SPMe布局下的浓度方差计算，用于评估静置阶段的扩散松弛效果。

- **设计模式**：辅助函数以 `_` 前缀命名（Julia惯例表示内部函数），由CycleSolver.jl调用。这避免了CycleSolver过度膨胀，将后处理逻辑解耦到PostProcessing.jl中。

## 后续变更 (2026-04-07)

- **Dict 访问替换**: 3 处 `multi_spme_layout` Dict 访问替换为 `case.layout` 和 `case.geometry` 字段直接访问
- 核心后处理逻辑不变
- 行数保持约 310 行
