# SPMe 装配合并与固定矩阵缓存

2026-10-24；父提交 `d3a41aa84aa44f7503162e8e66997e19f0a9bed0`，分支 `codex/src-physics-modularization`。范围为 10 个源码文件与三个既有测试，另附本记录和 KB 回执。AGENTS.md 与历史 .DS_Store 不提交。

## 差异及数据来源审查

提交前审查当前全部 src 暂存、未暂存及未跟踪实际差异；审查时暂存 src 与未跟踪 src 均为空。

| 文件/位置 | 审查结论 |
| --- | --- |
| src/SPMe.jl:1–58 | SPMe_element 合并到 SPMe，机械耦合、扩散、BC、blockdiag 共用。分布式 e/I_e/T_e/workspace 由 CallModel:165 及迁移的真实单元测试显式提供；普通模式不使用这些可选参数。nothing 不用来补回丢失状态。分布式质量及 none 刚度直接读 mesh 缓存；full 刚度和电解液刚度继续从当前状态重装配。普通 constant 分支沿用原 param 缓存，其质量已归一化而不重复缩放。 |
| src/SPMe.jl:62–77 | 化学缓存在初始化时从真实网格及冻结参数装配。三质量矩阵归一化一次；电解液质量按真实 NE/SP/PE 子网格高斯点数、原顺序拼接。none 固定 K 使用 Ds×Arrhenius(Eac_D,1)，不伪造浓度或温度状态；为已有 full 副本切换 none 对照预存，full 不消费固定 K。 |
| src/SetCase.jl:94–97 | Case 保持原五字段/五参数构造。仅 SPMe+per_element_spme 初始化专用化学矩阵，保存于现有 mesh；不是每步缺缓存后自动重建。 |
| src/ElectrodeDiffusion.jl:11–31 | 拆出质量与刚度核心，原联合接口保留。系数、积分及负号不变；没有修订原粒子温度未传入的模型现状。 |
| src/ElectrolyteDiffusion.jl:11–41 | 拆出质量与刚度核心，原联合接口保留；当前浓度/温度仍直接用于 De 和 Arrhenius。无替代来源或缺键默认值。 |
| src/ThermalDistributed.jl:1–62 | 初始化固定 MT 与体导热 KT；运行时直接读取，FT 仍逐次组装。缓存仅适用于冻结参数与最终网格，缺失键直接失败。ThermalDistributed2D_BC 保留 copy(KT)，SPMe blockdiag 分配新总阵；下游不原位改缓存。Ring/Polar 不变。 |
| src/Jellyrollmodel.jl:547 | 在最终 thermal2D 网格和 layer weights 设置后初始化热矩阵。手工 Cartesian fixture test/test_thermal_fem_rotation.jl:145 显式同步初始化；手工 Ring/Polar 示例不迁移。 |
| src/CallModel.jl:165 | 多单元调度保留，仅迁移统一 SPMe 调用，线程工作区与单元状态仍来自真实 caller。 |
| src/JuBat.jl:44 | 删除旧 SPMe_element 导出，不保留兼容包装。 |
| src/Variables.jl:191,240；src/Parallelsolution.jl:408 | 仅同步旧函数名注释，无物理逻辑改变。 |

未引入静默 haskey/get/get! 默认值、参数替代、历史补零、提前跳过或吞异常。none/full 是明确物理模式；固定矩阵是首次装配缓存，不是缺失运行历史的重建。缓存只读约定写于初始化函数，修改网格/参数须重新完成显式初始化。本审查不代表全部物理路径已验收。

## 唯一最终候选文字门

先读 Simplify/baseline.md 和上一批 interface_cleanup_20261024.md。上一批保存 identity_before 中所有 src 文件（含数据文件）及 Project.toml、Manifest.toml、example/testexample.jl 共 71 项，与父提交原始 blob、CRLF 表示或未改变工作区原件（规范化 CRLF/LF 后逐字节相同）全部匹配。保存父日志可用，不重跑父版。

唯一命令：`C:/Users/user/AppData/Local/Programs/Julia-1.11.2/bin/julia.exe -t1 --startup-file=no --project=. example/testexample.jl`；`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`；nθ20、60 s、纯文字。

运行退出 0；420 个热单元、441 个热节点、20 个记录时刻（包含初始列）。父版精确比较 21/21 通过；正式 75 比较 20/21，唯一差异仍是既知 hoop stress min：正式 -1.9297e+00 MPa，候选及父版 -1.9299e+00 MPa，未重冻结、未宣称正式门全通过。运行前后源码与保护项身份稳定。证据目录 `output/spme_matrix_cache_20261024/` 保存身份、命令、完整差异、日志、退出码及使用现有 verify_coupling_time75.py 的两份精确比较，不重冻结。

本批仅运行上述一次文字基线；三个修改测试未执行，未运行 full/constant 专项、绘图或长循环，不能把本配置通过扩大为这些路径已验收。测试按三路径同步 KB：test_case_batch_a.jl、test_review_p1_voltage_consistency.jl、test_thermal_fem_rotation.jl；不覆盖无关 KB dirty。
