# CZM 分离输出修复发现

## Requirements

- `mix` 模式打印最大有效分离位移 `sqrt(max(δ_n,0)^2 + δ_t^2)`。
- `model1` 模式打印最大法向分离位移。
- 两种模式均必须打印各自对应的全时空真实峰值。
- 保留现有结果键，不破坏现有调用者。

## Research Findings

- `DamageState` 已有 `δ_max_eff`，`bilinear_traction_state` 在创建新有效分离包络时更新它。
- `Variables.jl` 当前只初始化 `czm δ_max_n`，`PostProcessing.jl` 也只还原 `czm δ_max_n [m]`。
- `CallModel.jl` 在力学求解前从已提交状态读取 `δ_max_n`；`Solve.jl` 后续更新 CZM，但在 `Variable_update!` 前未刷新摘要。
- `czm_output_to_variables` 能写回当前位移、损伤、牵引与法/切向分离，但当前 `Solve` 没有调用它。
- 60 s `mix` 诊断已证明：历史法向峰值少报一步，逐单元分离历史为零。

## Test Contract

- 真实对象、无 mock。
- 期望值使用手工指定的法/切向分离和物理尺度，不用生产峰值 helper 构造。
- 能捕获的生产缺陷：缺少有效值键、漏写本步分离、打印选错分量。

## RED Evidence

- `test/test_czm_separation_output.jl` 首次运行退出码 1。
- 受控 `CZMResult` 的手算有效分离为 `[5.0, 3.0]`，期望空间峰值 `[5.0]`；当前两个键均返回 `nothing`。
- 失败为 2 个明确 assertion failure，0 error，证明测试捕获的是缺失功能。

## Implementation Findings

- `czm_output_to_variables` 现按单元计算 `hypot(max(δ_n,0), δ_t)`，写入 `czm separation effective` 及本步空间峰值。
- `Variables.jl`/`PostProcessing.jl` 已对称新增 `czm δ_max_eff` 和逐单元有效分离，物理还原仍只乘一次 `scale.δ_czm`。
- `Solve.jl` 在 `update_czm_damage!` 收敛返回后立即写回当前 CZM 场，再从已提交 `DamageState` 刷新法向/有效历史峰值，最后才 `Variable_update!`。
- `czm_max_separation_key` 将 `mix` 映射到 `czm δ_max_eff [m]`、`model1` 映射到保留的 `czm δ_max_n [m]`，未知模式显式失败。
- `testexample.jl` 与 `couple_example.jl` 现按模式选择结果键与文字标签。

## GREEN Evidence

- 受控有效分离 2/2 pass。
- 模式选择 4/4 pass，包括未知模式 `ArgumentError`。
- nθ=8、60 s 真实 `Solve` 5/5 pass：末步法向场不再为零，法向/有效历史峰值与全部 `CZMSnapshot` 手算一致。
- 正式 60 s `testexample.jl` (`mix`) 退出 0，19 个记录层，打印 `CZM 最大有效分离位移: 1.1207e-12 m`。
- 正式 60 s `couple_example.jl` (`model1`) 退出 0，打印 `CZM 最大法向分离位移: 1.1011e-12 m`，并成功输出三张最终场图。
- nθ=80 全尺寸快照交叉验证：法向的结果历史/快照峰值/已提交状态均为 `1.101106737768e-12 m`；有效分离三路均为 `1.120658115988e-12 m`。
- 有效峰值位于 `t=60 s`、NE–NCC 单元 5048；其法向 `1.092584023242e-12 m`、切向 `-2.492684598731e-13 m`，独立端点平均投影复算的有效值为 `1.120658115988e-12 m`。
- 全尺寸历史对齐误差：法向严格 `0.0 m`，有效值 `1.262177448354e-29 m`（浮点舍入）。
- 回归：CZM 尺度 28/28、basic 因子分解 17/17、分界面本构 71/71 通过。
- 本次 `couple_example` 图像 SHA-256：温度 `540fe42f...f978c`，环向应力 `165793bb...f1aa`，切向剪应力 `fadcb7b8...b285`。冻结档案已显式标记原应力 PNG 自 v9 起过期，本任务不擅自重冻结。
