# example/testexample.jl 代码简化基线

- **Baseline ID**: `testexample-20261007T152634+0800`（v11，任务59 P1+P2）
- **状态**: PASS（exit code 0）
- **重冻结方式**: 2026-10-07 用户在 JuBat 会话中明确指令"接受差异"（来源：当日 zcode 会话用户裁决，progress 已按时间序记录）后重冻结；复用技术评审通过的候选运行 `output/testexample/p2r2_gate/run.log`（其后源码未变，仅文档/KB；wall 23.638 s）；另批最终源码实跑 `couple_example.jl` 登记三图
- **入口**: `example/testexample.jl`
- **命令**: `OPENBLAS_NUM_THREADS=1 GKSwstype=100 julia -t1 --startup-file=no --project=. example/testexample.jl`（退出码 0 为执行方记录，评审未复跑）
- **环境**: Julia 1.11.2，1 thread，`GKSwstype=100`，`OPENBLAS_NUM_THREADS=1`
- **Git HEAD**: `ea70f341de37f5329e1c6a6ad236160afb9ccdaf`（任务59 P1+P2 为未提交工作树；源码身份以 46 文件内容清单为准——含 5 个参数文件；聚合口径见 preflight AGGREGATION_NOTE）
- **脚本 SHA-256**: `709f45ea2aa15747cd0997c435bb940dc57fc78a10646f0757875e63bbf0dffd`
- **标准 TSV 源码聚合 SHA-256**: 见 metrics.toml `source_manifest_tsv_sha256`（TAB 分隔、LF 拼接、末尾无换行；尾行 `# aggregate_sha256` 不计入聚合）；整文件 SHA 另记 `source_manifest_file_sha256`。未来运行须与本表记录精度严格一致（v10 迁移差异另列 metrics `scientific_metrics_changed`）
- **关键配置**: `fix_inner=false`（BC 92 节点/184 DOF）、`area_loss_enabled=false`、`geo_nonlinear=true`、`j2_plasticity=true`

本基线取代 `testexample-20260902T171549+0800`（v10，已归档 `archive/testexample-20260902T171549+0800/`）。

任务59 变更摘要：
1. **P1 无阈值面积权重**：`factor=1−D` 消费 max 映射损伤；独立 fractured/D≥0.99 停流与热源清零删除；全失活带载抛 `ErrorException`（含电流/时刻）、静置 `status=4` OCV 均值诊断；非法 D/映射显式失败；分流截止集与 Solve/相位同界限源（`param.cell`）；末步/相位/混合耗尽按方向电压截止终止。
2. **P2 机械边界**：CZM 域固定模式 `F=O∪I∪S∪E`（补回 a/b）；自由模式 `F=O∪(E\{b})`（释放 12 个层起点、a 经外圈属性固定、b 自由）。bonded 域同步 a_B/b_B 圈属性。BC 计数 nθ=80：CZM true 81/81/184、false 0/81/92。

科学变化（用户 2026-10-07 接受）：最大有效分离 1.1209e-12→1.0469e-12 m；环向 −1.6724/3.8807→−1.6795/3.9212 MPa；切向剪 −0.60032/0.99966→−0.60328/0.60634 MPa；GP Mises 4.9749→4.7052 MPa。其余指标（网格/步数/电压/容量/温度/损伤/断裂/塑性）按记录精度一致。

couple_example 三图为**独立配置门**（czm.model=model1、geo_nonlinear=false、j2_plasticity=false、mechanicalmodel=none）——与 testexample 的 mix/geo/J2-on 文字门不同物理；其自身日志：分离 1.0454e-12 m、环向 −1.6795/3.9203 MPa、剪 −0.60337/0.60593 MPa、wall 19.978 s。温度图与 v10 逐位一致；两张应力图为该独立配置在新边界下的自身结果。
