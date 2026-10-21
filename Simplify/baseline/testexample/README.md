# 纯文字行为基线（v12）

- Baseline ID：`testexample-20261020T161116+0800`，退出码0；Julia1.11.2，单线程，OPENBLAS_NUM_THREADS=1、GKSwstype=100、--startup-file=no。
- 入口：`example/testexample.jl`，物理时长60s，nθ=20；geo/J2开启，fix_inner=false，面积反馈关闭。
- 用户2026-10-20明确要求“重建基线，不生成图像，将θ设为20”“将生图代码移出基线”。基线只比较文字科学指标及网格/步数，计时不作强制一致门。
- 本次先运行当前源码，尚未实施接口减层。HEAD `3df40fee0b9bf9433736956f335f1cfdf9268b05` 加既有未提交修改；51份源码/配置逐文件SHA见[source_manifest.tsv](source_manifest.tsv)，冻结源码原件位于本次规划目录refreeze_source，不以HEAD单独复现脏工作树。
- 网格：420单元/441节点；20记录点。最终电压3.9386V，温度298.15–299.15K，容量0.0833Ah。
- 结构化指标：[metrics.toml](metrics.toml)；实际输出：[run.log](run.log)；环境：[preflight.log](preflight.log)。只冻结日志实际打印的精度，不把未采集的精确电压或收敛计数沿用为当前实测。
- 此次重冻结不证明现有物理修复整体正确。v11不一致是改前事实，旧原件保存在[历史v11](../archive/testexample-20261007T152634+0800/README.md)；新旧不同网格，不作单变量归因。
- `couple_example.jl`是独立绘图功能，已经移出本基线的运行及验收链；本次未调用或修改其生图代码，不生成图像。历史图像记录只存在旧档案中。

后续简化必须与本表科学指标按记录精度、网格和步数严格一致；不得再次重冻结来吸收本次接口减层差异。
