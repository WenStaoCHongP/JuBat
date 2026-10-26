# 接口收敛与单 SPMe 二维热旧路径清理

日期：2026-10-24；基于 `7dd8a9d69d1d3ceef67623842e3cd4d103698d72`，分支 `codex/src-physics-modularization`。本批包含 README、10 个 src 文件及三个测试更新。AGENTS.md 与历史快照 .DS_Store 保留，不纳入提交；没有推送。

## 源码差异审查

审查覆盖当前全部暂存、未暂存和未跟踪 src 差异；审查时暂存 src 和未跟踪 src 均为空。以下为修改后行号，完整差异和身份保存于 `output/interface_cleanup_20261024/`。

| 文件/位置 | 修改与数据合同判断 |
| --- | --- |
| src/CallModel.jl:94,165,183,201,232 | 代表量显式提供完整 StandardVariables，实际支路使用既有 workspace 池；热源直接调用唯一实现；删去 MT 乘 1。CZM 摘要和原物理模式分支保留，无缺失状态替代。 |
| src/CzmSolve.jl:405 | 仅删除 gp_history 开关及其 false 拒绝分支，真实 GP 历史读取、缺史失败、收敛和回滚不变。 |
| src/Initialisation.jl:2,33 | 在自定义 y0 判断前显式拒绝已删除单 SPMe+distributed2D 组合；非法模型消息读取真实 model，无默认替代。MultiSPMe 初始化不变。 |
| src/JuBat.jl:59 | 移除已删除薄包装的导出，无运行回退。 |
| src/Option.jl:30 | 删除无消费者的 inner_exit_only 字段及注释。 |
| src/SPMe.jl:2,40,123–130,149–157,209 | 单 SPMe 装配拒绝 distributed2D；SPMe_variables 合并重复公式，按实际热模式选择 workspace 或完整字典。distributed2D 直接使用传入 workspace，缺失不会重建；原生索引/赋值失败。none 使用有定义的初温，其他未注入温度使用 only 读取真实状态；I_app 的原 Current 来源不变。保留两个表面浓度赋值。非分布式额外电流键按同一物理模式判断。 |
| src/SetMesh.jl:57 | 删除把未合并网格冒充 bonded 的五参数兼容构造器；真实七参数构网调用不变。 |
| src/Solve.jl:147,191,471–474 | 错误消息读取真实 solveType；平均耗时除以实际调用次数，初始成功 CallModel 已计数，输出路径至少一次。比例零总耗时分支保留。 |
| src/ThermalDistributed.jl:295,304–307,334–347 | 删除单模型广播，逐单元变量和有效面积直接读取真实来源；保留原面积检查、公式、计算顺序、单位和结果键。删除仅转调包装，无零值/旧值/另一来源替代。 |
| src/Variables.jl:204,227 | 仅同步 workspace 模式注释；初始化键不变。 |

本批没有新增静默 `get`/`get!` 默认值、`haskey` 替代、三元补值、补零、沿用旧历史、吞异常或跳过物理计算。新增 `nothing` 参数表示显式工作区输入模式，分布式不会缺失后创建替代工作区；非分布式创建完整字典是明确模式。没有删除真实物理状态检查。结论来自实际差异及调用来源审查，不表示全库无错误。

## 唯一最终候选文字门

读取 `Simplify/baseline.md` 和 `docs/validation/solve_observer_cleanup.md` 后，核对保存父候选回执的全部 46 个 src/Project/Manifest/example 文件与 HEAD。45 项原始 blob 或 CRLF 哈希匹配；Solve 原件为混合换行，用 SHA `ec9ab33f63dd3f43235b659e5959c930f8658c6950c6f4555af872108bc8f116` 的已保存 `Solve.jl.before` 按已保存 task.diff 重建，重建 SHA 精确等于回执 `247a9ac2c6e37771460a447455fb040f538acb78e750c6c7c29cf1ad3d837f78`；仅统一 CRLF/LF 后与 HEAD 逐字节相同。因此保存父日志适用于当前 HEAD，无需重跑父版。

命令：`C:/Users/user/AppData/Local/Programs/Julia-1.11.2/bin/julia.exe -t1 --startup-file=no --project=. example/testexample.jl`。环境：`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。nθ20、60 s，无绘图、长循环或额外专项测试。

| 比较 | 结果 |
| --- | --- |
| 候选退出码 | 0 |
| 保存的当前父版文字指标 | 21/21 精确一致 |
| 正式 75 | 20/21，保留既知差异 |

正式 75 的唯一差异为 hoop_stress_min_mpa：参考 -1.9297e+00，候选 -1.9299e+00 MPa。没有重冻结，也不称正式 75 通过。候选为 420 热单元、441 热节点、20 记录时刻；本门仅证明此配置文字行为，未执行修改后的三个测试，不能扩展为所有入口或完整物理验收。

证据目录：`output/interface_cleanup_20261024/`，包含父身份核对、运行前后身份、实际命令、完整 diff、候选日志/退出码、父指标、两份比较和 receipt。运行前后全部跟踪文件及 AGENTS/.DS_Store 的 SHA256 一致。比较工具复用现有 `tools/verify_coupling_time75.py`。

三个测试按路径同步：`test_case_batch_a.jl`、`test_czm_arc_gp_history.jl`、`test_review_p1_voltage_consistency.jl`。不使用目录级同步覆盖其他 KB 文件；同步结果及 KB 提交身份另见同名交付回执。主仓提交须在审阅后进行。

静态格式检查补记：普通 `git diff --check` 报新增 CRLF 行的行尾诊断及 ThermalDistributed.jl 包装删除后留下的末尾空行；使用 `core.whitespace=cr-at-eol` 后只剩该末尾空行诊断。保留已验证源码字节，没有为格式再次改写源码，不将此检查描述为通过。
