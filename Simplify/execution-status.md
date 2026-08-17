# src 简化计划执行状态

日期：2026-08-05

## 结果摘要

- 37 个 plan 已完成或审计后保留。
- 5 个 plan 因 D1/D2 设计缺陷暂不执行。
- 修改 4 个生产文件，新增 2 个聚焦测试文件。
- 36 个顶层 `src/*.jl`：10,027 → 9,889 行，净减 138 行。
- 22/22 个测试文件通过；冻结 `example/testexample.jl` 的全部严格指标与 PNG SHA-256 一致。

## 已执行源码批次

| Concern | 文件 | 结果 |
|---|---|---|
| D3 热边界双胞胎 | `src/ThermalDistributed.jl` | 保留两个原函数名，直接采用原位实现，删除两个 `!` 变体 |
| CZM 旧兼容入口 | `src/CouplingState.jl` | 删除无调用、未导出的六参数适配器 |
| Silent catch | `src/Solve.jl` | 删除最终热节点结果写入的静默吞错 |
| CSV 写入守卫 | `src/CsvExport.jl` | 七处同构 try/catch 收敛为一个私有 helper |

## 阻塞计划

| 重复簇 | Plan | 原因 |
|---|---|---|
| D1 | `CycleData.md`、`CycleSolver.md` | 当前 callback 方案改变求解器行为且缺少可靠逐步 characterization；不能直接删约 400 行复制实现 |
| D2 | `VariableKeys.md`、`Variables.md`、`CallModel.md` | 当前常量化范围、键数量与测试设计不一致；新增全量常量层可能增加概念和代码量 |

这些阻塞项需要先重写为以现有行为测试为起点的单一 concern，再进入下一轮源码修改。
