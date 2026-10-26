# M3：basic 成功出口复用已验收响应

2026-10-24；分支 `codex/src-physics-modularization`；父提交 `72757a7f28de0b545f7af97eb7a1f0677a4f2851`。

## 修改与提交前审查

- 仅修改 `src/CzmSolve.jl:538` 附近：删除成功出口重复的 `assemble_czm_system`，原失败装配分支改由 `if !converged` 控制，内容未变。
- `483–501` 的 committed-state 装配、原 tol 残差检查、trial/accepted 分离和响应保存保持不变；成功出口使用已验收的 `separations/tractions`。
- `533–535` 失败回滚和失败残差报告保持不变。成功路径从验收到 `fill_czm_result!` 没有新的装配写入；该函数逐项复制响应到结果数组，然后才调用 `commit_czm_state!`，保留共享 workspace 的消费时序。
- 已审查全部未提交 src：仅本文件有差异，暂存 src 和未跟踪 src 均为空。没有新增默认值、haskey/get 回退、nothing 替代、状态补零、异常吞掉或防呆。`!converged` 明确区分失败恢复与成功复用，不跳过平衡验收。
- 其他修改 `AGENTS.md` 及历史快照 `.DS_Store` 不纳入提交；测试源码未改。

## 验证

Julia 1.11.2，`-t1 --startup-file=no --project=.`；`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。

| 入口 | 结果 | 覆盖 |
| --- | --- | --- |
| `test/test_czm_arc_gp_history.jl` | 201/201，退出 0 | basic/gp_basic 等真实求解入口、非均匀 GP 历史、失败回滚/续算及原有状态合同 |
| `test/test_czm_output77.jl` | 151/151，退出 0 | CZM 结果发布、GP 响应、快照隔离和既有输出合同 |
| `example/testexample.jl` | 父版 21/21 精确一致，退出 0 | nθ20、60 s；420 热单元、441 热节点、20 个记录时刻（含初始列）；纯文字 |

只运行一次最终候选 60 s 基线。父版复用 `output/m2_trial_reuse_20261024/candidate_60s.log`；其保存的 src、Project/Manifest 与入口共 71 项身份已核对匹配父提交（允许 Git CRLF/LF 表示差异），没有重跑父版。

正式 75 仍为 20/21，唯一不同是 `hoop_stress_min_mpa`：正式 `-1.9297e+00`，父版及本候选 `-1.9299e+00`。差异列表与父批一致，没有新增偏差；未重冻结，不能称正式 75 全通过。测试日志中的预算耗尽/时间上限警告来自预期失败用例，相关断言通过。

验证前后源码、依赖、入口、测试及保护项 SHA-256 一致。未运行长循环、绘图或额外专项；本批通过不等于全部物理工况验收。

## 证据

目录 `output/m3_success_reuse_20261024/` 保存 runner、完整命令、源码差异、父版身份检查、运行前后身份、三份日志与退出码、比较 JSON、receipt 和 SHA-256 清单。

```text
python -X utf8 output/m3_success_reuse_20261024/verify_m3.py "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
```

- `candidate_60s.log` SHA-256：`464b984bb6c1b1db4477da7b7a3d9521230e4e827aafa1cc3efeacd1159e520d`
- `receipt.json` SHA-256：`db323a989efe89fea318d9940f18fb37871c7977eb57ff91a40822eb16f361a8`

状态：本批相关测试及父版文字回归通过，按用户授权提交；正式 75 的既知差异继续保留。
