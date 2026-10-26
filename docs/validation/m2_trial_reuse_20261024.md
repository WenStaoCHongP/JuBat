# M2：复用 cohesive GP 已计算的 trial 状态

2026-10-24；分支 `codex/src-physics-modularization`；父提交 `a9a02c7e288b05926cae497cb7dffbf1dead9cec`。

## 修改与提交前审查

- `src/Materialmatrix.jl:205–240`：公开 `bilinear_tangent` 签名及导出不变，仍自行调用 `bilinear_traction_state`；内部 `_bilinear_tangent` 接收 previous/trial，复用原切线计算主体。
- `src/czm.jl:254–264`：装配器第一次计算牵引后，将已有 `new_state` 传给内部核，不再重复映射。GP separation/traction 捕获、trial 保存及装配顺序不变。
- 加载/卸载判断仍由 `bilinear_damage_gradient` 比较 previous 与 trial；核只读二者，不原位覆盖历史。D/D_visc、压缩分支、model1/mix、饱和规范化和原非法历史检查保留。梯度内部重算包络仍保留，本批不改本构算术。
- 已审查当前全部未提交 src（暂存、未暂存、未跟踪）：仅上述两个文件有差异，暂存与未跟踪 src 均为空。未新增条件、防呆、haskey/get 默认值、nothing 回退、替代状态、跳过或吞异常；原 GP 可选捕获条件未改。
- 原有 `AGENTS.md` 修改及历史快照 `.DS_Store` 不纳入提交，验证前后 SHA-256 一致；测试源码未修改。

## 验证结果

环境：Julia 1.11.2，`-t1 --startup-file=no --project=.`；`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。

| 入口 | 结果 | 范围 |
| --- | --- | --- |
| `test/test_czm_constitutive_physics.jl` | 184/184，退出 0 | 本构切线有限差分、压缩/混合模式/黏性/饱和/非法史、真实 GP 装配及 accepted history 不被改写 |
| `test/test_czm_output77.jl` | 151/151，退出 0 | CZM 结果发布、GP 捕获与快照、历史隔离和既有输出合同 |
| `example/testexample.jl` | 与父版 21/21 精确一致，退出 0 | nθ20、60 s，420 热单元、441 热节点、20 个记录时刻（含初始列），纯文字 |

父版复用 `output/spme_matrix_cache_20261024/candidate_60s.log`；已核对该次保存的 src、Project/Manifest 与入口共 71 项身份匹配当前父提交（允许 Git CRLF/LF 表示差异）。本批只运行一次最终候选 60 s 仿真，没有重跑父版。

正式 75 基线为 20/21：唯一差异仍为 `hoop_stress_min_mpa`，正式值 `-1.9297e+00`、父版及本候选均为 `-1.9299e+00`。差异列表与父批完全一致，未新增偏差、未重冻结，不能称正式 75 全通过。本批证据不代表长循环或全物理验收。

初次汇总因 Python cp950 输出含中文的绝对路径而失败；仿真与两项测试均已退出 0。汇总改用 `-X utf8`，通过 `--compare-only` 复用日志完成比较，未重跑 Julia。相关源码、依赖、入口、测试及保护项验证前后身份不变。

## 可追溯证据

`output/m2_trial_reuse_20261024/` 保存 runner、命令、源码差异、父版身份检查、运行前后身份、测试及基线日志/退出码、两份精确比较、receipt 和 SHA-256 清单。运行命令：

```text
python -X utf8 output/m2_trial_reuse_20261024/verify_m2.py "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
# 仅重新汇总已保存结果：同命令加 --compare-only
```

- `candidate_60s.log` SHA-256：`6f9be1b9a7a7cca8a11e04ac353a23f5d77e9755c0a2939a6c383f62f63b67b6`
- `receipt.json` SHA-256：`32367122c84721d2328161d898c4c157e76c4cd3f9569e804106cbe12c37aca7`

状态：本批相关测试与父版文字回归通过，可按用户授权提交；保留正式 75 的既知差异。
