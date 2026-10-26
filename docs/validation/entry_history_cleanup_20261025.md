# Solve 入口、闲置温度历史与 Cell 字段清理验证

日期：2026-10-25。分支：`codex/src-physics-modularization`。父提交：`5428de2f45249019bda02d116c80a7ec5331c7c6`。

## 本批范围

- `src/Solve.jl`：按用户裁决删除预期长度计算与长度断言，消除单模型为获得长度再次调用完整初始化。保留原初始化、外层 copy、自定义 opt.y0、恢复态复制及类型/有限性检查；`Initialisation.jl` 未修改。
- `src/Variables.jl`：仅删除 `thermal2D temperature history` 定义。节点温度历史、后处理单元平均温度及其余六个历史字段保留。
- `src/SetParams.jl`：删除无库内消费者的 `Cell.Nr_th`、`Nθ_th`、`n_windings` 及专属说明。
- 两个既有测试文件的三处短恢复向量断言由 DimensionMismatch 改为真实索引产生的 BoundsError。配置恢复断言保留；这不代表所有非法长度均会在入口被拦截。

## 验证及证据

命令：`python -X utf8 tools/verify_entry_history_cleanup.py`。

证据目录：`output/verify_entry_history_cleanup/`。运行器保存命令、日志、退出码、源码差异、前后身份、父版快照、逐项对照及 receipt.json。

固定环境：Julia 1.11.2，`-t1 --startup-file=no --project=.`；`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。

| 检查 | 结果 |
| --- | --- |
| 父版证据身份与 HEAD 对照 | 71 个输入匹配 |
| 新专项：单 SPMe+集总、多 SPMe，默认/自定义初态，恢复求解，节点/单元温度，Cell 其余字段 | 28/28 |
| 既有布局与状态测试 | 89/89 |
| 既有温度、电流、两次微循环、相位配置恢复测试 | 65/65 |
| 本批专项合计 | 182/182 |
| 60 秒 example/testexample.jl 与已验证父版文字结果 | 21/21 |
| 历史 formal75 文字对照 | 20/21，既有差异不变 |
| 五项 Julia 命令退出码 | 全部 0 |
| 验证期间输入及受保护文件身份 | 前后一致 |

新专项针对四组入口配置，父版和候选的完整保留变量及求解物理输出使用 isequal 精确比较，排除已批准删除的键和运行耗时（调用次数仍比较）。独立节点平均核对允许浮点运算次序造成的舍入差异（16eps）；没有放宽父版结果对照。

60 秒候选只运行一次，nθ=20，420 热单元、441 热节点，20 条时间记录（含初始记录）。父版复用 `output/verify_multispme_output/candidate_60s.log`，其 SHA256 为 `7e4d29e2a576adc66179d9e5120aedbcb8d7b9bc7cabd06b185309608bb9b021`。

formal75 的唯一差异仍为最小环向应力：历史 `-1.9297e+00 MPa`，父版与候选均为 `-1.9299e+00 MPa`。本批未新增偏差，没有重冻结基线，也不宣称历史 formal75 完全通过。本次通过表示相对已验证父版的清理回归通过，不代表其他配置或长期物理验收。

## 工具身份核对修正

首次执行在启动 Julia 前因父版 Solve.jl 混合 LF/CRLF 的字节身份未匹配而停止。原失败记录保留为 `initial_parent_identity_check.json`。运行器按 HEAD 内容与未变行的原始换行重建候选字节，只有完整 SHA256 与上次记录一致、且归一化内容与 HEAD 一致时才接受。恢复出的字节存于 `parent_original_bytes/src/Solve.jl`。未改变科学源码、基线指标或数值容差，未产生失败仿真或额外 60 秒运行。

## 提交审查与归档

提交前审查暂存、未暂存及未跟踪 src 范围：本批仅三个文件的删除，无新增兜底、默认替代、跳过计算或吞异常；没有额外防御分支。显式长度检查删除依据用户授权，短恢复向量仍由真实热自由度索引显式失败，既有回归已验证。差异空白检查通过。

测试/工具四个文件已逐项同步 KB，SHA256 一致，其他 136 个 test/tools 文件哈希不变。KB 提交：`ef80034`（`test: sync entry and history cleanup validation`）；清单及哈希见 `kb_sync_receipt.json`。

主仓库仅提交三个源码、三个测试、一个运行器及本记录。`AGENTS.md` 和历史快照中的 `.DS_Store` 保持原有未提交状态，不纳入本批。
