# 代码简化文字行为参照

2026-10-21：当前源码按用户决定回到3df40fe并仅保留独立P1诊断开关。**原60s简化入口尚未建立适用于新目标的正式冻结基线。** 不把旧v11/v12直接当作当前通过门，不自动重冻结。

- [v12历史说明](baseline/testexample/HISTORICAL-v12.md)：testexample-20261020T161116+0800，基于3df40fe加混合未提交候选；五份原件保持原字节，其README内“当前/后续”表述属于当时历史语境。
- [v12原始指标](baseline/testexample/metrics.toml)、[实际日志](baseline/testexample/run.log)、[源码身份](baseline/testexample/source_manifest.tsv)只供历史追溯，不用于要求3df独立迁移后复现混合候选科学结果。
- [v11原件归档](baseline/archive/testexample-20261007T152634+0800/README.md)保留；旧历史差异不被新结果补证。
- 本次P1验收：Julia1.11.2，单线程，nθ20、60s纯文字，同一3df40fe配置迁移前后17结果+3网格+退出码精确21/21；独立专项86项和短程对照9项通过。证据在output/validate_p1_migration/及[迁移任务](<D:/OneDrive/Desktop/lib知识库/lib-knowledge-base/raw/repos/jubat/docs-archive/planning-with-files/接口减层与诊断开关/旧候选清理_20261021/task_plan.md>)。这是本次实验对照，不写入新的metrics.toml，不等于重冻结。

入口仍为example/testexample.jl，nθ20、60s。运行环境为Julia1.11.2、-t1、--startup-file=no、OPENBLAS_NUM_THREADS=1、GKSwstype=100。计时不作科学指标；不调用couple_example或绘图/PNG门，独立绘图功能保留。

- 2026-10-22 用户授权新增[独立损伤半循环基线](baseline/czm_cycling_n20_halfcycle/README.md)：`example/czm_cycling_n20.jl` Case2 原配置，nθ20、SOC0.65，放电1800s→静置600s；入口 `example/check_czm_n20_halfcycle.jl`。源码为57755d4加72本构候选与共享setup提取，完整身份和原件在该基线记录；冻结两个阶段实际指标、网格/步数与退出码。仅用于该配置的后续回归，不改变旧60s/v11/v12指标，不等于73/74或全面物理验收。
