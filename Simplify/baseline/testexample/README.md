# example/testexample.jl 代码简化基线

- **Baseline ID**: `testexample-20260831T212819+0800`（v9）
- **状态**: PASS（exit code 0）
- **重冻结方式**: 2026-08-31 用户授权的结构层热膨胀系数变更后完整重跑（testexample.jl 纯文字快速门；本批未改绘图代码，未跑 couple_example.jl）
- **入口**: `example/testexample.jl`
- **命令**: `julia --startup-file=no --project=. example/testexample.jl`
- **环境**: Julia 1.11.2，Plots 1.40.9，1 thread，`GKSwstype=100`
- **Git HEAD**: `696f2aa2b4df4a22c01eb1cf7085e1cf5e53f3a2`（运行时 alphaT 变更尚未提交）
- **脚本 SHA-256**: `63edef136f21d1df99117d66dacddc845f6ba7587f264dec6967c24dfc3d61eb`
- **标准 TSV 源码聚合 SHA-256**: `0f015cdff6304715bd9391e3be911e76332b5871bf6a2a017c023ac1c3549de6`（TAB 分隔、LF 拼接、末尾无换行；尾行 `# aggregate_sha256` 不计入聚合）
- **说明**: 基线建立时工作树已有用户修改；因此以脚本哈希和 46 个 Julia 文件的内容清单为准，而不是仅以 HEAD 为准。

本基线取代 `testexample-20260830T172856+0800`（v8）及更早。本征应变仍按材料分层计算，但 **SP/PCC/NCC.alphaT 不再为零**（30e-6 / 23e-6 / 17e-6 1/K，用户指定）；CZM 与宏观应力使用历史列同一时间层的温度/SOC，非更新输出列保持最近一次有效力学解。`testexample` 的 19 次 CZM 更新全部实际执行并收敛。

## 冻结结果

| 指标 | 基线值 |
|---|---:|
| thermal elements | 1682 |
| thermal nodes | 1763 |
| result time steps | 19 |
| initial voltage | 4.0367 V |
| final voltage | 3.9438 V |
| voltage drop | 0.0929 V |
| final capacity | 0.0833 Ah |
| minimum temperature | 298.15 K |
| maximum temperature | 299.00 K |
| final CZM D_max | 0.0000% |
| final CZM D_mean | 0.0000% |
| maximum normal separation | 7.4202e-13 m |
| fractured elements | 0 |
| CZM converged updates | 19 / 19 |
| hoop/tangential 应力范围 | −1.8433~+4.1344 / −0.2618~+0.2167 MPa |
| result PNGs | 已移至 couple_example.jl（本门不含图；哈希自 v9 起失效，见下） |

## 后续比较规则

每个代码简化批次完成后，必须用同一 Julia 版本、线程数、入口和参数重新运行。

- 必须 exit code 0。
- 网格规模、时间步数和上表全部科学结果必须与基线在脚本打印精度下完全一致。
- 图片门移至 `example/couple_example.jl`（输出 `output/couple_example/`，三张 final PNG 哈希记录于 metrics.toml [artifact] 注释）；仅当修改涉及绘图/后处理代码时运行并核对。
- wall-clock、模块耗时、耗时占比不作严格相等要求，只记录趋势。
- 任一科学指标不一致时，该简化批次不得继续，先定位差异或回滚。
- 本次未生成 `output/simple_coupling_debug.log`，因此该文件不属于基线。

## 文件

- `metrics.toml`: 机器可读的关键指标与比较策略。
- `source_manifest.tsv`: 运行时所有 Julia 源文件和入口脚本的 SHA-256。
- `preflight.log`: 本次直接重冻结的源码与输出身份记录。
- `run.log`: 2026-08-31 v9 完整重跑的控制台记录与重冻结来源。

## 基线 v2（2026-08-22 重冻结，用户宏观参数修正）

- 触发：用户修正 `src/parameters/Jellyroll.jl` 宏观力学参数（PCC.E→70 GPa/0.33、NCC.E→110 GPa/0.34、PE/NE.E_coat→1 GPa、SP.E→750 MPa/0.35；提交 `1a74411`）。电化学/热全部指标与 v1 一致；仅两项力学输出移动：`maximum normal separation 1.2557e-14 → 1.2572e-13 m`（刚性箔下微小变形放大一个量级，仍为 ~零损伤工况）与 PNG SHA（上表已更新为 v2 值）。
- v1 SHA：`4ba6207c…e932`（历史比对记录见 `Simplify/baseline.md` 批次表）。

## 基线 v3（2026-08-22 重冻结，Φ 缝默认完美粘结 spec v1.5）

- 触发：Φ 配对节点默认合并（`8ab8863`，spec v1.5 §3.4）——网格拓扑变化。电/热全部指标与 v2 一致；仅 `separation 1.2572e-13 → 6.6820e-15 m`（Φ 粘结约束了匝间自由度）与 PNG SHA（上表已更新为 v3）。v1/v2 历史见上文与 baseline.md 批次表。

## 基线 v4（2026-08-24，cohesive 法向定向 + 父热单元温度路径）

- 触发：cohesive 局部法向统一按 `host_inner_elem → host_outer_elem` 定向，避免把物理张开记为受压；同时删除不进入残差的细力学节点温度插值场，热应变继续采用父热单元平均温差。
- 电/热、网格、步数、零损伤和 19/19 CZM 收敛保持不变；`maximum normal separation` 更新为 `1.5174e-12 m`，PNG 更新为 92,736 B / `0946646a...`。
- 用户明确要求不重跑；本次基线直接复用 2026-08-24 04:34:11 -06:00 已完成运行的输出，并同步刷新当时源码清单。

## 基线 v5（2026-08-29，PE.Omega 符号更正 + 层分辨宏观应力在线导出）

- 触发：用户授权两项变更——① `PE.Omega` 由 −7.28e-7 更正为 +7.88e-7（嵌锂膨胀；原值为未审 placeholder）；② 宏观应力彻底层分辨化：耦合流程求解中经 `export_macro_stress` 在线收割导出 `diffusion stress xx/yy/xy/vonMises [Pa]`，czm-off 走按需固体工具函数 `thermal_diffusion_stress_2D`（mesh_bonded 域）。
- 电/热全部指标、网格（1682/1763）、步数（19）、零损伤与 19/19 CZM 收敛保持不变；`maximum normal separation` 更新为 `7.0037e-13 m`（Ω_PE 载荷变化）。
- 绘图产物由单张 `testexample_results.png` 改为四张 final-field PNG（温度、环向应力、切向剪应力、环向应力径向剖面）；层分辨后环向应力范围 −1.78~+4.00 MPa（NE 涂层拉 / PE 涂层压交替）。
- 全量 test/ 套件与 test_electrode_coat_modulus.jl 同步通过（J2 积分测试载荷已按新 Ω 符号重标定）。

## 基线 v6（2026-08-29，验证入口拆分：减负两级流程）

- 触发：用户要求基线验证减负——`example/testexample.jl` 拆为纯文字结果输出版（60 秒快速门），全部代码（含绘图）移至新建 `example/couple_example.jl`（输出 `output/couple_example/`）。
- 常规修改的验证流程（AGENTS.md §9.6 已同步更新）：① 只运行受影响部分的验证（相关 test/ 文件与示例）；② 跑一次 60 秒 `testexample.jl`，文字指标须与本基线按记录精度一致。
- 全量绘图验证（`couple_example.jl`）仅在修改涉及绘图/后处理时运行；本次拆分后运行确认四张 PNG 哈希与 v5 完全一致（`540fe42f/6c970edf/9cd0ebff/87711dce`），行为无损。
- 新增文字指标：最终环向应力范围 −1.7766~+3.9971 MPa、切向剪应力范围 −0.2789~+0.2461 MPa（层分辨在线导出）。

## 基线 v7（2026-08-30，α/β 同批分层化物理变更）

- 触发：用户授权本征应变分层分辨率化——`eigenstrain_of(param, mt)` 取代跨层均匀 `α_eff/β_n/β_p`：α 按层取 `param.X.alphaT`（PE 1.5e-5、NE 8e-6，SP/PCC/NCC 在 Jellyroll.jl 显式置零），电极膨胀 β=Ω/3 只作用于本层涂层（NE→Δsn、PE→Δsp），箔/隔膜零本征应变；α_eff/β_n/β_p 死参链从 eigenstrain 元组与全部求解器签名中删除。
- 电/热全部指标、网格（1682/1763）、步数（19）、零损伤与 19/19 CZM 收敛**逐位不变**；couple_example 温度场 PNG 哈希逐位不变（`540fe42f`，热路径未动的硬证据）。
- 力学指标按预测移动：`maximum normal separation 7.0037e-13 → 9.4407e-13 m`；环向应力 −1.7766~+3.9971 → −1.3952~+3.8570 MPa；切向剪应力 −0.2789~+0.2461 → −0.3739~+0.3952 MPa；三张应力 PNG 哈希更新（`ecdc9f58/2c29e35b/5124dec3`）。
- 全量 `test/runtests.jl` 34/34 通过（8m00s）；`unit_czm_eigenstrain` 60/60（PE 层解析比 1.000）；j2 积分测试按分层物理重标定（箔屈服驱动 Δsoc 0.3→1.0，箔不再随电极膨胀）。

## 基线 v8（2026-08-30，应力时间对齐 + 非更新步有效保持 + 三图输出）

- 触发：用户要求修复应力历史错位和非更新步伪零；CZM 更新移到时间推进前，使用当前 `CallModel` 同时间层温度/SOC，`latest_macro_stress` 在非更新输出列保持最近一次实际力学解。
- SP/PCC/NCC.alphaT 按用户决定继续保持 0；电化学、热学、网格、19 步、零损伤和断裂数均与 v7 一致。
- 力学直接影响项更新：最大法向分离 `9.4407e-13 → 9.6486e-13 m`，环向应力 `−1.3952~+3.8570 → −1.3954~+3.8603 MPa`，切向剪应力 `−0.3739~+0.3952 → −0.3743~+0.3955 MPa`。
- 删除额外径向散点图；`couple_example.jl` 只生成温度、环向应力、切向剪应力三张 Q4 云图。

## 基线 v9（2026-08-31，结构层热膨胀系数由零改为物理值）

- 触发：用户指定 `src/parameters/Jellyroll.jl` 三个结构层的热膨胀系数——`SP.alphaT 0 → 30e-6`、`PCC.alphaT 0 → 23e-6`、`NCC.alphaT 0 → 17e-6`（1/K）。此前 v7/v8 一直按"显式置零、留待敏感性分析"处理，本批解除该置零。`eigenstrain_of` 因此对隔膜与两个集流体返回非零热应变项；求解器代码路径未改（早已是逐层查表），仅同步了 `src/czm.jl` 中已失效的文档串。
- 电化学、热学、网格（1682/1763）、步数（19）、零损伤（D_max=D_mean=0、断裂 0）与 19/19 CZM 收敛**逐位不变**——alphaT 不进入热残差，无热膨胀反馈回热模型。
- 力学指标移动：`maximum normal separation 9.6486e-13 → 7.4202e-13 m`；环向应力 `−1.3954~+3.8603 → −1.8433~+4.1344 MPa`；切向剪应力 `−0.3743~+0.3955 → −0.2618~+0.2167 MPa`。
- **A/B 隔离验证**：同进程内只切换这三个 alphaT 跑两次，`alphaT=0` 一侧**逐位复现 v8 冻结的环向范围 −1.3954~+3.8603 MPa**。这同时证明 `696f2aa` 给脚本加的 `opt.czm.model = "mix"` 对宏观应力无影响（D≡0，混合模式分支未激活），v8 数值在该脚本漂移下仍然有效。
- 分层实测（集流体屈服评估，一次性诊断）：峰值 von Mises 在箔上而非涂层——PCC `3.8439 → 3.9917 MPa`、NCC `3.6180 → 4.1124 MPa`（峰值 +4~14%）；但平均值近乎翻倍——PCC `0.4459 → 0.8155 MPa`、NCC `0.5623 → 0.9355 MPa`。本工况屈服利用率仍低：PCC 6.65%（σ_y 60 MPa）、NCC 2.06%（σ_y 200 MPa），**不屈服**。注意本工况 ΔT 仅 0.85 K；在 `|ΔT| ≤ 20 °C` 包络下五层 α 失配达 22e-6，热项将成为主导贡献之一。
- 门禁：`test/runtests.jl` **32/32** 通过（106 testsets，零 Fail/Error/Broken，20m01.7s，exit 0）；`testexample.jl` exit 0。
- **PNG 哈希自本版起失效**：本批未改绘图/后处理代码，按 AGENTS §9.6 未跑 `couple_example.jl`；但应力场已移动，metrics.toml `[artifact]` 记录的三张云图哈希不再对应当前代码，已标记 `png_hashes_stale_since = "v9"`，暂不作为门禁。下次涉及绘图的批次须重跑并刷新。
- **运行环境差异**：基线原记录的 `D:/Julia-1.11.2/bin/julia.exe` 在当前机器上不存在；v9 使用同版本 Julia 1.11.2（`C:/Users/19303/AppData/Local/Programs/Julia-1.11.2`），单线程、`GKSwstype=100`、`--startup-file=no` 不变；项目依赖本次 instantiate 到 `C:/Users/19303/.julia`。
- **源码清单修正（附带）**：`source_manifest.tsv` 共 18 行变化——2 行是本批改动（`src/parameters/Jellyroll.jl`、`src/czm.jl`）；7 行是 v8 冻结后累积但从未重新登记的源码漂移（`example/testexample.jl`、`src/CouplingState.jl`、`src/CsvExport.jl`、`src/Initialisation.jl`、`src/SetParams.jl`、`src/ThermalDistributed.jl`、`src/Tools.jl`）；**8 行是遗留的 LF 归一化哈希**，与其余 29 行的"工作树原始字节"约定不一致，本次统一为后者。清单此前处于两种约定混用状态，现已内部自洽。聚合规则未变并已重新验证（去 CR、TAB 分隔、LF 拼接、末尾无换行、尾行不计入）。
