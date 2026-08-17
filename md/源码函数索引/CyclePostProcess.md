# CyclePostProcess.jl

- **源文件**: `src/CyclePostProcess.jl`
- **行数**: 235 行
- **函数/struct 计数**: 10 个函数
- **职责**: 循环阶段与周期结果汇总、SOH、终止判据、用户可见摘要和循环结果绘图。
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/13_耦合验证方案.md`

## 数据结构

本文件无独立 struct；消费 `PhaseResult`、`CycleResult`、`CyclingResult` 兼容对象。

## 函数清单

- `phase_termination_symbol` — L4-L12：把求解终止原因映射为 `:time` / `:voltage`。
- `state_concentration_variance` — L14-L42：计算 REST 阶段前后颗粒浓度方差。
- `postprocess_phase_result` — L44-L102：从阶段 `solve_result` 汇总容量、电压、温度、损伤和最终状态。
- `postprocess_cycle_result!` — L104-L119：汇总循环容量、效率、损伤和最高温度。
- `append_cycle_result!` — L121-L137：追加循环序列并按配置保存详细结果。
- `update_soh_and_capacity!` — L139-L148：锁定初始容量并更新 SOH。
- `print_cycle_summary` — L150-L154：打印单循环摘要。
- `check_cycle_termination` — L156-L175：检查 SOH 与 CZM 断裂终止条件。
- `print_cycling_summary` — L177-L193：打印全部循环总结。
- `plot_cycling_results` — L195-L235：输出容量、损伤、效率、温度及组合 PNG。

## 职责边界

- 不执行时间积分；求解流程仍在 `CycleSolver.jl`。
- 不采集原始时间步快照；该职责在 `CycleData.jl`。
- 不写 CSV；该职责统一在 `CsvExport.jl`。
- CZM 统计使用 `CzmPostProcess.get_damage_statistics`，不重复实现损伤统计。

## 省略项

无。

### [DEBUG]

| 行号 | 内容 | 用途 |
|---|---|---|
| L161/L169 | SOH/CZM 终止 `@warn` | 用户可见终止原因 |
| L179-L191 | 循环完成摘要打印 | 用户可见运行总结 |
| L232 | 绘图保存提示 | 用户可见输出位置 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|---|---|---|
| L49 | 空电压历史时 `V_start` 可回退为 `NaN` | 失败阶段可能传播 NaN |
| L59/L61 | 无温度数据时回退初始温度 | 纯电化学场景的合理默认 |
| L107 | 零充电容量时效率回退 0 | 不完整循环的显示语义 |
| L145 | 初始容量为零时 SOH 回退 1 | 第一循环失败时可能掩盖无容量 |
| L178 | 无有效循环时最终 SOH 回退 1 | 仅影响打印摘要 |
| L220 | 库仑效率图范围固定 95–105% | 异常效率可能被裁剪 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L54-L61 | 温度最大值有三路数据回退 | 可抽取温度摘要 helper |
| L107 | 库仑效率三元除零保护 | 可读性可接受 |
| L167 | CZM 断裂终止含两个条件 | 可读性可接受 |
| L178 | 最终 SOH 有双条件三元回退 | 可读性可接受 |
