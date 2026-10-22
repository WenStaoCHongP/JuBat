# Case2 nθ20 放电–静置独立损伤基线

冻结日期：2026-10-22。用户明确要求使用 `czm_cycling_n20` 循环脚本的基准工况，只做一次放电–静置半循环。**实际推进一次，两个阶段均完成，退出0**；不进入充电/rest2，不覆盖原60s、v11或v12指标。

## 身份与原件

- HEAD：`57755d4681b3719ad36f5d3672f36c95e01f6b2f`；分支 `codex/src-physics-modularization`。运行源码包含72本构修复和 `czm_cycling_n20.jl` 的共享setup提取；未实施73/74。
- 原件目录：[output/check_czm_n20_halfcycle/20261022_190021_ca21e2](../../../output/check_czm_n20_halfcycle/20261022_190021_ca21e2/)。[source.diff](../../../output/check_czm_n20_halfcycle/20261022_190021_ca21e2/source.diff) 保留父提交到候选差异；[parent_cycling_script.jl](../../../output/check_czm_n20_halfcycle/20261022_190021_ca21e2/parent_cycling_script.jl) 是原循环入口的父提交blob。
- [source_manifest.csv](source_manifest.csv)：运行前73个源码/依赖/脚本完整相对路径→SHA256；[source_snapshot](../../../output/check_czm_n20_halfcycle/20261022_190021_ca21e2/source_snapshot/) 保存相应原件，包含所有src、Project/Manifest、新薄runner和循环脚本。复制后73项及运行结束后当前73项均与运行前manifest完全一致。
- [artifact_manifest.csv](artifact_manifest.csv)：完整日志、真实初态/相位结果/GP端点、快照及源码包等原件的路径、SHA256和字节数。
- [metrics.toml](metrics.toml) 冻结全部phase_summary值、初态、网格、步数、按界面损伤与退出码；[phase_summary.csv](phase_summary.csv) 保持运行原字节。全部数值是本次实测，不借旧校准注释作通过门。

## 固定工况

复用 `example/czm_cycling_n20.jl` 的Case2（Case4/8的重复基准）：单正极耳位于累计角范围50%，W=4mm、L=7.425mm；负极耳沿原默认44π；总极耳冷却面积59.4mm²，`tab.h=802.9425303710864 W/(m²K)`。极耳电导/厚度、几何和所有其余材料值沿原脚本/参数源。

SOC0.65；初始温度及环境298.15K；SPMe的Nn/Ns/Np=10/5/10，Nrn/Nrp=10/10，分布式二维热、多SPMe、极耳冷却。电流：放电+5A，静置0A；电压边界2.5/4.2V；放电1800s、静置600s；阶段步长范围[1,10]s，Crank–Nicolson/自动步长。原opt.dt=[0.5,10]在相位入口按生产接口暂换为dt_cycle，退出恢复引用。

CZM mix/gp_basic，强度×0.085、刚度×0.1、Gc×2，PCC/NCC各自材料及I/II等参关系不变；黏性ON、τ=5s；tol=1e-3、max_iter=100、load_steps=10、update_interval=1。geo、J2、area_loss、continuous_feedback、winding_prestress均OFF。不调参、不加预算。

边界沿原脚本：fix_inner=false，outer_bc=fixed_xy，fix_start/end=nothing，endpoint_variant=keep；旧默认实际固定外圈O及E\{b}，补a外圈属性，b自由；不替换为task63的omit_a边界。

薄runner枚举并清除外部JUBAT环境覆盖，`run_metadata.toml` 的外部JUBAT表为空；固定Julia1.11.2、-t1、--startup-file=no、OPENBLAS_NUM_THREADS=1、GKSwstype=100。完整实际选项/参数见configuration.txt及configuration.bin。

## 初态、衔接与实测结果

初态使用生产 `apply_initial_soc!(...,0.65)` 和 `ModelInitialisation_MultiSPMe`，真实新 `MechState`：u/D/Dv=0、mechanics_initialized=false、GP历史尚未初始化（记录缺席，不伪填GP零值）。放电final_state直接成为静置initial_state，温度、损伤、位移和GP历史连续继承。

| 实测项 | 放电 | 静置 |
|---|---|---|
| 全局范围 | 0–1800s | 1800–2400s |
| 真实结束原因 | time_limit | time_limit |
| 电压起/末 | 3.7674739081440656 / 3.22617285171737 V | 3.3422373454389236 / 3.4486740696871747 V |
| 容量 | 2.499999999999995Ah | 0Ah |
| 温度最高 / 末态节点均值 | 306.93091925753527 / 306.80655798854093K | 306.93091925753527 / 301.062717273598K |
| 单元Dmax / Dmean | 0 / 0 | 0.3175081217681377 / 0.0007650634726519105 |
| GP Dmax / Dvmax | 0 / 0 | 0.6350162435362754 / 0.6267639348948426 |
| 记录时刻数（含相位初态） | 210 | 96 |
| 真实推进时间区间 / 已接受机械快照 | 209 / 209 | 95 / 95 |

网格：nθ20、420热元/441热节点、3360力学bulk元/5072力学节点、1680cohesive元/3360GP。末温范围300.8201158364389–301.10575915553966K。损伤按严格D>0统计：静置末PE_PCC为10/840单元、10/1680GP，NE_NCC为0；全部单元及GP fractured均0。真实初/相位端点见accepted.bin与GP CSV；相位原件PhaseResult包含温度、电压、应力完整历史。

首个**记录到非零单元损伤的快照时刻**为全局2011.0s（静置局部210.99999999999991s，快照265），其Dmax=0.0005901158647954551；这不是精确连续起裂时刻。应力指标来自geo/J2关闭路径的**单元中心恢复应力**；vonMises全历史最大值不是J2 GP最大值。全部应力冻结值在phase_summary/metrics中。

这些结果证明本次配置确实经历损伤；静置增损是观测，不据此独断根因，也不等于完整循环、长期受损循环、生产精度或73/74验收。

## 运行与读取

实际物理命令（上述固定环境下，cwd为JuBat根）：

```powershell
& 'C:\Users\user\AppData\Local\Programs\Julia-1.11.2\bin\julia.exe' -t1 --startup-file=no --project=. example/check_czm_n20_halfcycle.jl 20261022_190021_ca21e2
```

stdout/stderr原件为run.log，退出码原件为exit_code.txt（0）。再次使用该入口应传新的RUN_ID；运行结果与冻结数值比较时须先匹配源码、依赖、环境和配置，不能以哈希相同替代物理判断。

仅分析已保存原件的可复用脚本是 [analyze_saved_results.jl](../../../output/check_czm_n20_halfcycle/analyze_saved_results.jl)，它不调用Solve。Julia参数在ChooseCell运行时include生成闭包；configuration.bin读取前按原probe→build两次ChooseCell注册相同真实定义。首个未注册闭包的读回错误日志保留；随后只构建配置（不推进物理）即成功读回全部原件，saved_evidence_analysis.exit=0。

首次启动目录 [20261022_185940_ca40e2](../../../output/check_czm_n20_halfcycle/20261022_185940_ca40e2/) 保留首部docstring解析失败及exit1；该次尚未加载JuBat或推进任何物理。修正入口并在启动前纠正J2关闭路径的应力键后，才进行上述唯一实际半循环。未为文书、评审或读取错误重复仿真。
