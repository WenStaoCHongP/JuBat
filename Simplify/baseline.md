# 代码简化行为基线

代码简化必须保持现有行为。每个修改批次完成后，先执行对应基线，再进入下一批次。

| 基线 | 入口 | 状态 | 关键判定 | 档案 |
|---|---|---|---|---|
| `testexample-20260815T011730-0600` | `example/testexample.jl` | PASS | exit code、网格/步数、科学结果、CZM 收敛与 PNG SHA-256 完全一致；耗时除外 | [详细记录](baseline/testexample/README.md) |

## 固定运行环境

```powershell
$env:GKSwstype = '100'
$env:JULIA_NUM_THREADS = '1'
& 'D:\Julia-1.11.2\bin\julia.exe' --startup-file=no example\testexample.jl
```

结构化期望值见 [`baseline/testexample/metrics.toml`](baseline/testexample/metrics.toml)，完整控制台输出见 [`baseline/testexample/run.log`](baseline/testexample/run.log)。

## 基线重建记录

| 日期 | 原因 | 旧基线 | 新基线 | 关键变化 |
|---|---|---|---|---|
| 2026-08-15 | 力学周向离散改为直接继承热网格角段，并纳入已审阅的严格循环状态契约 | `testexample-20260806T031217-0600` | `testexample-20260815T011730-0600` | 力学周向段数与 1682 个父热单元严格对应；最大法向分离更新为 `1.2557e-14 m`；PNG 与 46 文件源码清单同步重建 |
| 2026-08-06 | 修复 CZM 插值矩阵绑定未合并候选网格、而温度状态来自活动合并网格的尺寸错误 | `testexample-20260805T031305-0600` | `testexample-20260806T031217-0600` | CZM 19 次更新实际执行并全部收敛；最大法向分离由错误路径的 0 更新为 `1.3527e-14 m`；PNG 基线同步更新 |

## 简化批次记录

| 日期 | 批次 | 生产代码变化 | 局部验证 | 强制基线 | 结果 |
|---|---|---:|---|---|---|
| 2026-08-05 | ThermalDistributed D3：保留原函数名并直接采用原位变体、删除 `!` 变体 | `+6/-133`，净删 127 行 | `smoke_thermal_bc.jl` 18/18 | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | CouplingState 六参数 CZM 兼容入口清理 | `-17` 行 | CZM 定向测试通过（4+9 assertions，1 项既有 broken；smoke OK） | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | Solve 最终热数据 silent catch 清理 | `-4` 行 | 热边界 smoke 18/18（最终结构） | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | CsvExport 七处写入容错统一 | `+22/-37`，净删 15 行 | CSV guard/public path 10/10 | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-21 | 堆芯塌陷 Batch 1：bulk 残差/切线统一入口 `assemble_bulk_residual_tangent` + Option 六默认关子选项 | `+62/-5`（Option +7，czm +55/-5） | `test_czm_mech_core.jl` 8/8 testset（含接线前后逐位等价）、`test_czm_option_defaults.jl` 3/3、`test_assemble_coupled_system.jl` 通过；全套 24/24；`verify_czm_standalone.jl` 快照逐位一致 | 所有科学指标及 PNG SHA-256 完全一致（`4ba6207c…`；另经 stash A/B 复核接线前后同 SHA） | PASS |

本轮 36 个顶层 `src/*.jl` 的 PowerShell 物理行统计由 10,027 降至 9,889，净减 138 行。

## 审计决策记录

| 日期 | 文件/入口 | 决策 | 证据 |
|---|---|---|---|
| 2026-08-05 | `SetCase.jl` 五参数构造器 | 保留 | `SetCase` 主路径内部直接调用 |
| 2026-08-05 | `parameters/LGM50.jl` | 保留 | 多个活跃示例调用 |
| 2026-08-05 | `parameters/Ring.jl` + `ring.jl` | 保留 | 多个热验证与工具脚本调用 |
| 2026-08-05 | `parameters/Enertech.jl`、`parameters/Northrop.jl` | 保留 | `ChooseCell` 已公开并文档化；仓库 grep 不能证明外部无调用 |
| 2026-08-05 | `install.jl` | 保留 | 独立、明确的依赖安装入口，无重复职责 |
