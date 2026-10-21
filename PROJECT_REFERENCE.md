# JuBat 按需技术参考

这是项目入口的辅助导航，不是新的科学验收报告。规则以 [AGENTS.md](AGENTS.md) 为准；配置默认值、历史计数和实现状态以当前源码、对应任务证据核实。本次精简不修改模型或重跑实验。

## 源码与文档入口

`KB = D:/OneDrive/Desktop/lib知识库/lib-knowledge-base`；手册底稿在 `KB/raw/repos/jubat/md-raw/`，旧 docs 在 `KB/raw/repos/jubat/docs-archive/`。

| 主题 | 主要源码 | KB 底稿 |
|---|---|---|
| 参数与归一化 | src/Option.jl、src/parameters/Jellyroll.jl | 01_参数定义与归一化.md、15_颗粒与极片模量区分.md |
| 几何、边界 | src/Jellyrollmodel.jl、src/CzmMesh.jl、src/CzmBC.jl | 02_几何与网格.md、03_边界条件.md |
| 电化学 | src/SPMe.jl、src/CallModel.jl | 04_电化学模型_SPMe.md |
| 热学 | src/ThermalDistributed.jl、src/ThermalPolar2D.jl | 05_热模型_二维分布式.md |
| 力学、CZM | src/Mechanical.jl、src/czm.jl、src/CzmSolve.jl | 06_内聚力模型_CZM.md |
| 热阻 | src/czm.jl 及热耦合调用路径 | 07_界面热阻模型.md |
| 多 SPMe、分流、循环 | src/Parallelsolution.jl、src/Solve.jl、src/CycleSolver.jl | 08_逐单元算法.md、09_分流求解器.md、10_参数传递与模块架构.md |
| 验证方案 | 相关测试、example 与任务记录 | 11_电化学验证方案.md、12_热模型验证方案.md、13_耦合验证方案.md |

手册和历史状态可能早于源码；先核对调用路径，不能以手册替代实现证据。理论推导从 [Theory/README.md](Theory/README.md) 按章读取。

## 模型数据流

每热单元一个 SPMe 时，状态按各单元化学向量加节点温度组合：
`yt_global=[yt_chem[1]; …; yt_chem[ne]; T_nodes]`，化学块包含负/正颗粒浓度及电解液浓度。

电化学产热传给热模型，单元温度反馈电化学；化学/热本征应变驱动力学与 CZM。损伤对热阻、有效面积等支路的实际作用须按当前开关和装配路径核实，不能从示意图推断所有支路已开启。

CZM 接口：
`solve_czm_step(czm_mesh, case.mech, param, F_ext, opt.czm; dT_elem=…, Δsoc_n_elem=…, Δsoc_p_elem=…) -> CZMResult`。

宏观应力与颗粒应力分开；历史设计中 `export_macro_stress` 输出 xx/yy/xy/vonMises，czm-off 的 `thermal_diffusion_stress_2D` 使用 bonded 域。具体字段与有效性检查以当前实现为准。

## 热学归一化与模量

参考尺度与禁止混用的约束见 AGENTS。常用还原关系：

| 量 | 无量纲表达 |
|---|---|
| 密度 | ρ* = ρ/ρ_ref |
| 比热容 | c* = c·ρ_ref·L³·T_ref/(t₀·P_ref) |
| 体积热容 | (ρc)* = ρc·L³·T_ref/(t₀·P_ref) |
| 导热率 | k* = k·L·T_ref/P_ref |
| 热源 | q* = q·L³/P_ref = q/scale.q |
| 面积 | A_phys = A_mesh·scale.L² |

颗粒应力尺度为 `scale.E_p/E_n`；宏观参考为 `scale.E_coat`。CZM 锚点在量纲阶段定义：`scale.σ_czm=PCC.σ_max`，`scale.δ_czm=2*PCC.G_c/PCC.σ_max`。实际量值与材料参数从当前参数集读取，不用旧表“典型值”替代。

## 配置检索提示

Jellyroll 耦合示例通常使用 `model="SPMe"`、`thermal_enabled=true`、`thermalmodel="distributed2D"`、`per_element_spme=true`、`mechanicalmodel="full"`、`czm.enabled=true`。电流、时间、步长与边界须显式匹配工况；这些不是全局默认值。

`Option` 与 `CzmOptions` 的字段、默认值和校验均以 [src/Option.jl](src/Option.jl) 为入口。常查项：`model`、`update_interval`、`iter_method`、`max_iter/tol`、`load_steps`、`arc_length_alpha`、`viscous_enabled`、`viscous_tau`、`geo_nonlinear`、`winding_prestress`、`j2_plasticity`、`area_loss_enabled`、`continuous_feedback`。

边界设计的易错语义（修改前核对 `resolve_mechanical_bc` 及测试）：

- `fix_inner=true` 的旧规则固定 O∪I∪S∪E；false 固定 O∪(E\{b})。bonded 路径须同步消费相应圈属性。
- `outer_bc=:fixed_xy|:radial_slide`；径向小位移滑移与 geo=true 的组合按入口契约拒绝。
- `fix_start/fix_end` 为显式布尔或 nothing；nothing 沿用旧规则。
- `endpoint_variant=:keep|:omit_a|:omit_b|:omit_ab`：omit 去除对应端点的人工圈属性补充与 S/E 固定贡献，保留真实几何支承；非法值入口报错。
- a 是卷绕起点的人工外圈端点；b 是末端的人工内圈端点，仅在 fix_inner=true 激活。历史 nθ80 计数或 omit_b 等价结果不是任意网格的常数。
- `area_loss_enabled` 的任务59契约是无阈值 factor=1−D；D=1 支路面积失活、I_e=0，全失活带载报错、静置 status=4 诊断。修改时连同分流/循环调用者与失败路径验证，不能把后续面积反馈方案直接当成已实现。
- `friction_mu` 是预留参数；字段存在不代表接触/摩擦已经装配。

## 输出与运行导航

常用结果有 `time [s]`、`cell voltage [V]`、`temperature`、`thermal2D temperature [K]`、`thermal2D element current/area`、`heat_source_fields`。分层热源包括 NE/PE 反应热和可逆热、SP 欧姆热、PCC/NCC 集流体热（单位 W/m³）。CZM 包括 `czm D_max/D_mean`、`czm n_fractured` 和法向分离统计。导出前检查实际键名、单位和有限性；不在此复制全部字段表。

- 快速行为对照：[example/testexample.jl](example/testexample.jl) 与 [Simplify/baseline.md](Simplify/baseline.md)。
- 独立绘图：[example/couple_example.jl](example/couple_example.jl)，按用户要求运行。
- 其他示例按任务检索 `example/`；热学验证在 `example/热模块验证/`，测试/工具归档及恢复流程见 AGENTS。
- 调试字段：`debug_coupling`、`debug_log_path`。调试日志/图像/CSV 遵循脚本输出目录约定。
- 历史提速配置 `-t8` 与 BLAS=8 来自特定 14700K/源码配置，不能直接认定当前最优；行为基线始终按规定的单线程配置。

原项目/全局/Theory 入口的完整原件、哈希与本次整理记录位于用户指定备份目录：
`D:/Graduate student/codex账户/记忆备份/20261021T150210+0800/`。
