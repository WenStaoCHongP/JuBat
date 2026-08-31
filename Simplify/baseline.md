# 代码简化行为基线

代码简化必须保持现有行为。每个修改批次完成后，先执行对应基线，再进入下一批次。

| 基线 | 入口 | 状态 | 关键判定 | 档案 |
|---|---|---|---|---|
| `testexample-20260831T212819+0800` | `example/testexample.jl`（纯文字快速门，60 s） | PASS | 结构层热膨胀系数取物理值（SP 30e-6 / PCC 23e-6 / NCC 17e-6 1/K）；电/热逐位不变，力学输出移动；三张 final 云图哈希自 v9 起失效，须重跑 `example/couple_example.jl` 后方可作为门禁 | [详细记录](baseline/testexample/README.md) |

## 固定运行环境

```powershell
$env:GKSwstype = '100'
$env:JULIA_NUM_THREADS = '1'
& 'D:\Julia-1.11.2\bin\julia.exe' --startup-file=no --project=. example\testexample.jl
```

> v9（2026-08-31）在 `D:` 盘不可用的机器上重跑，改用同版本 Julia 1.11.2
> （`C:\Users\19303\AppData\Local\Programs\Julia-1.11.2`）；线程数、`GKSwstype`、
> `--startup-file=no` 与 `--project=.` 均未变。Julia 版本一致即满足门禁环境要求，
> 路径本身不是判定项。

结构化期望值见 [`baseline/testexample/metrics.toml`](baseline/testexample/metrics.toml)，完整控制台输出见 [`baseline/testexample/run.log`](baseline/testexample/run.log)。

## 基线重建记录

| 日期 | 原因 | 旧基线 | 新基线 | 关键变化 |
|---|---|---|---|---|
| 2026-08-31 | 用户指定：结构层热膨胀系数由显式置零改为物理值（SP.alphaT→30e-6、PCC.alphaT→23e-6、NCC.alphaT→17e-6 1/K） | `testexample-20260830T172856+0800` | `testexample-20260831T212819+0800` | 电/热指标、网格（1682/1763）、19 步、零损伤与 19/19 收敛**逐位不变**（alphaT 不进热残差）；分离 9.6486e-13→7.4202e-13 m；环向应力 −1.3954~+3.8603→−1.8433~+4.1344 MPa；切向剪 −0.3743~+0.3955→−0.2618~+0.2167 MPa。A/B 隔离：alphaT=0 一侧逐位复现 v8 环向范围，同时证明 `mix` 脚本漂移对宏观应力无影响。全量测试 32/32。三张应力 PNG 哈希失效（本批未重跑绘图），已标记非门禁。源码清单 18 行变化，其中 8 行为遗留 LF 归一化哈希的约定统一 |
| 2026-08-30 | 用户授权：修复 CZM/应力历史时间错位与非更新步伪零，删除额外径向散点图；alphaT 保持 0 | `testexample-20260830T005629+0800` | `testexample-20260830T172856+0800` | 电/热、网格、19 步、零损伤不变；分离 9.4407e-13→9.6486e-13 m；环向/切向应力小幅更新；绘图收敛为三张 Q4 云图 |
| 2026-08-30 | 用户授权：α/β 同批分层化——`eigenstrain_of(param, mt)` 取代跨层均匀 α_eff/β_n/β_p，电极膨胀只作用于本层涂层，SP/PCC/NCC.alphaT 显式置零；α_eff/β 死参链删除 | `testexample-20260829T231308+0800` | `testexample-20260830T005629+0800` | 电/热指标、网格、步数、零损伤、19/19 收敛**逐位不变**；separation 7.0037e-13→9.4407e-13 m；环向应力 −1.7766~+3.9971→−1.3952~+3.8570 MPa；温度 PNG 哈希不变、三张应力 PNG 更新（`ecdc9f58/2c29e35b/5124dec3`）；47 文件清单重建；j2 测试按分层物理重标定（箔屈服驱动 0.3→1.0） |
| 2026-08-29 | 验证入口拆分（减负）：testexample.jl 改纯文字 60 s 快速门，全部代码移至新建 couple_example.jl；常规修改仅跑受影响验证 + 快速门 | `testexample-20260829T224907+0800` | `testexample-20260829T231308+0800` | 求解指标与 v5 全同；新增应力范围文字指标；四张 PNG 由 couple_example 产出且哈希与 v5 一致；源码清单 47 文件（+couple_example.jl） |
| 2026-08-29 | 用户授权：PE.Omega 由 −7.28e-7 更正为 +7.88e-7（嵌锂膨胀，原值为 placeholder）；宏观应力层分辨化并耦合在线导出 | `testexample-20260824T043411-0600` | `testexample-20260829T224907+0800` | 电/热指标、网格、步数、零损伤、19/19 收敛全部不变；`maximum normal separation` 更新为 `7.0037e-13 m`；绘图产物改为四张 final PNG（层分辨环向应力 −1.78~+4.00 MPa，NE 拉/PE 压交替）；46 文件源码清单同步重建；J2 积分测试载荷按新 Ω 符号重标定 |
| 2026-08-15 | 力学周向离散改为直接继承热网格角段，并纳入已审阅的严格循环状态契约 | `testexample-20260806T031217-0600` | `testexample-20260815T011730-0600` | 力学周向段数与 1682 个父热单元严格对应；最大法向分离更新为 `1.2557e-14 m`；PNG 与 46 文件源码清单同步重建 |
| 2026-08-24 | 接受 cohesive 法向 host-inner→host-outer 定向，并收敛热—力温度路径 | `testexample-20260815T011730-0600`（含 v2/v3 修订） | `testexample-20260824T043411-0600` | 用户要求不重跑，直接复用已完成运行；分离更新为 `1.5174e-12 m`，PNG SHA 更新为 `0946646a...`，46 文件源码清单同步刷新 |
| 2026-08-06 | 修复 CZM 插值矩阵绑定未合并候选网格、而温度状态来自活动合并网格的尺寸错误 | `testexample-20260805T031305-0600` | `testexample-20260806T031217-0600` | CZM 19 次更新实际执行并全部收敛；最大法向分离由错误路径的 0 更新为 `1.3527e-14 m`；PNG 基线同步更新 |

## 简化批次记录

| 日期 | 批次 | 生产代码变化 | 局部验证 | 强制基线 | 结果 |
|---|---|---:|---|---|---|
| 2026-08-30 | 四层力学结构体重构（spec `2026-08-30-mechanics-struct-refactor`）：界面参数挂 CurrentCollector（Cohesive/CzmInterfaceParams/CzmParamCache/内容哈希全删）、`opt.czm::CzmOptions` 嵌套收敛 20 字段、装配缓存惰性挂 CohesiveMesh（K_bulk/标架/ws）+ BC 现算、MechState 聚合演化状态（damage_states 迁出网格、克隆链删除、收敛原位提交）、Λ 内联 scale.L/δ_czm、`solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt)` 终态签名 | 14 个 src 文件重构 + 16 测试适配 + 3 示例适配 | 全套 32/32（8m04s）；单步 A/B 探针与父提交 66c718e 17 位逐位一致；中途门×2 绿；顺带修复 bilinear_tangent 混合模式 4 处裸 K_n/K_t（model2 崩溃级） | **相对父提交 66c718e 全指标逐位一致**（分离 9.6486e-13 m、环向 −1.3954~+3.8603 MPa、三张 PNG 哈希不变 `540fe42f/9726a7c8/b3b43d7c`，用户核验）；基线数值零漂移，仅刷新源码清单（46 文件）与脚本哈希 | PASS（v8 基线沿用） |
| 2026-08-05 | ThermalDistributed D3：保留原函数名并直接采用原位变体、删除 `!` 变体 | `+6/-133`，净删 127 行 | `smoke_thermal_bc.jl` 18/18 | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | CouplingState 六参数 CZM 兼容入口清理 | `-17` 行 | CZM 定向测试通过（4+9 assertions，1 项既有 broken；smoke OK） | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | Solve 最终热数据 silent catch 清理 | `-4` 行 | 热边界 smoke 18/18（最终结构） | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-05 | CsvExport 七处写入容错统一 | `+22/-37`，净删 15 行 | CSV guard/public path 10/10 | 所有科学指标及 PNG SHA-256 完全一致 | PASS |
| 2026-08-21 | 堆芯塌陷 Batch 1：bulk 残差/切线统一入口 `assemble_bulk_residual_tangent` + Option 六默认关子选项 | `+62/-5`（Option +7，czm +55/-5） | `test_czm_mech_core.jl` 8/8 testset（含接线前后逐位等价）、`test_czm_option_defaults.jl` 3/3、`test_assemble_coupled_system.jl` 通过；全套 24/24；`verify_czm_standalone.jl` 快照逐位一致 | 所有科学指标及 PNG SHA-256 完全一致（`4ba6207c…`；另经 stash A/B 复核接线前后同 SHA） | PASS |
| 2026-08-21 | 堆芯塌陷 Batch 2：完全 GL/TL geo_nl 槽位（SVK + 初应力 K_G，ε* 内嵌 D-B2-1）+ basic/load_substep 求解器接线 | `+182/-38`（czm +159、CzmSolve +55/-38、CouplingState +6/-2） | `test_czm_geometric_stiffness.jl` 5/5（FD 切线/刚体转动/线性退化/自由膨胀/K_G 方向）、`test_czm_geo_c1.jl` 2/2（C1）、mech_core 8/8；全套 26/26；`verify_czm_standalone.jl` 快照逐位一致（含 b551dac arc 误伤被门禁捕获后回滚 `1bb0026`） | 所有科学指标及 PNG SHA-256 完全一致（geo_nl 默认关） | PASS |
| 2026-08-21 | 堆芯塌陷 Batch 3：PCC/NCC 平面应力一致 J2 塑性（物理箔本构 D-B3-0、mech_state 消费、收敛提交 D-B3-2） | `+246/-39`（CzmPlasticity 新建 +101） | `unit_czm_j2.jl` 5/5（含一致切线 FD）、`test_czm_j2_integration.jl` 4/4（C2-lite 屈服/回归锚/不可逆/报错）；全套 28/28；`verify_czm_standalone.jl` 快照逐位一致 | 所有科学指标及 PNG SHA-256 完全一致（塑性默认关） | PASS |
| 2026-08-22 | 堆芯塌陷 Batch 2'：卷绕预应力 σ₀（等应变卷入张力 + 对数累积压力，D-B2'-1）+ 用户宏观参数修正（PCC 70GPa/NCC 110GPa/E_coat 1GPa/SP 750MPa）→ **基线重冻结 v2** | `+165/-52`（链路 5 文件 + Jellyroll/SetParams 参数） | `test_czm_winding_prestress.jl` 6/6；全套 29/29 | **参数刻意变更，基线 v2 重冻结**：testexample 电/热指标不变、分离 1.2557e-14→1.2572e-13 m、PNG SHA→`b31ffb49…`；探针景观移动（basic 8/8 D=0.55、load_sub/arc 6/8）——v1/v2 档案同存 | PASS（v2 基线） |
| 2026-08-22 | 堆芯塌陷 Batch 2''：D13 网格探针（厚涂层细分 + split_KG + μ_crit 二分/主模态/单匝 DFT，八组双参考态） | `+302`（mesh/czm/probe/test） | `test_czm_thin_subdiv.jl` 4/4；全套 30/30；探针 Summary 与 v2 逐位 | 所有科学指标及 PNG SHA-256 与 v2 一致（细分/分离默认关） | PASS |
| 2026-08-24 | cohesive 法向定向 + 删除细力学节点温度插值场 | 当前脏工作区窄批次 | 温度/CZM 专项通过；全套 33/33；`testexample` exit 0、19/19 CZM 收敛 | 用户接受当前结果并直接重冻结；未在重冻结阶段再次运行 | PASS（v4 基线） |

本轮 36 个顶层 `src/*.jl` 的 PowerShell 物理行统计由 10,027 降至 9,889，净减 138 行。

## 审计决策记录

| 日期 | 文件/入口 | 决策 | 证据 |
|---|---|---|---|
| 2026-08-05 | `SetCase.jl` 五参数构造器 | 保留 | `SetCase` 主路径内部直接调用 |
| 2026-08-05 | `parameters/LGM50.jl` | 保留 | 多个活跃示例调用 |
| 2026-08-05 | `parameters/Ring.jl` + `ring.jl` | 保留 | 多个热验证与工具脚本调用 |
| 2026-08-05 | `parameters/Enertech.jl`、`parameters/Northrop.jl` | 保留 | `ChooseCell` 已公开并文档化；仓库 grep 不能证明外部无调用 |
| 2026-08-05 | `install.jl` | 保留 | 独立、明确的依赖安装入口，无重复职责 |
