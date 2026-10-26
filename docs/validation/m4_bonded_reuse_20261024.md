# M4：bonded 装配复用索引与 eigenstrain 系数

2026-10-24；分支 `codex/src-physics-modularization`；父提交 `cc184b5f4e4e64007e6a8b2dd5233cbb36b05871`。

## 修改与提交前审查

- `src/Mechanical.jl:513–539`：删除与 Vi_u/Vi_v 完全相同的 Vj_u/Vj_v 数组，行列索引复用 u/v 两组 DOF；coeff_u/coeff_v 合为一组 coeff。
- 四块刚度、x/y 导数配对、四块相加及 Assemble1D 调用顺序未变；K_vu 仍独立计算，没有用转置替换。factor 表达式未改。
- 核实 `Assemble` / `Assemble1D` 只读索引和系数，写入的是函数内新建的装配缓冲/结果；复用输入不会导致原位覆盖。
- 提交前审查全部未提交 src：仅 Mechanical.jl 有差异，暂存 src 和未跟踪 src 均为空。没有新增条件、防呆、默认值回退、状态替代、跳过或吞异常。
- 保留原有 AGENTS.md 修改及历史快照 .DS_Store；不提交这两项，验证前后保护项哈希一致。现有测试文件未修改。

## 验证范围与结果

Julia 1.11.2，`-t1 --startup-file=no --project=.`；`OPENBLAS_NUM_THREADS=1`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。

标准 60 s 算例走 CZM，不能单独覆盖 bonded，因此增加本批有针对性的父版装配对照脚本，不扩大运行整套热 FEM 测试。

| 检查 | 结果 | 覆盖 |
| --- | --- | --- |
| `probe_bonded.jl` | 14/14，退出 0 | 复用 `test_thermal_fem_rotation.jl:consumer66` 的真实 nθ12 分层网格配置；均匀与非均匀温度/SOC 输入；父版和候选 K 的尺寸、CSC 指针/行索引、Float64 数值位模式及 F 位模式一致，两个方向载荷非零，输入字典不变 |
| `example/testexample.jl` | 父版 21/21 精确一致，退出 0 | nθ20、60 s、420 热单元、441 热节点、20 个记录时刻（含初始列），纯文字 |

bonded 对照矩阵为 4072×4072、72816 个存储项。父版函数从父提交的 Mechanical.jl 原始 blob 提取，仅更名后在同一 JuBat 模块中执行；依赖和算例输入一致，不改父版算术，也不替换工作区源码。提取身份及原函数 SHA-256 见 parent_extraction.json。

只运行一次最终候选 60 s 基线。父版复用 `output/m3_success_reuse_20261024/candidate_60s.log`，保存的 src、Project/Manifest 与入口共 71 项身份匹配父提交（允许 Git CRLF/LF 表示差异），不重跑父版。

正式 75 仍为 20/21：唯一差异为 hoop_stress_min_mpa，正式 `-1.9297e+00`、父版与本候选 `-1.9299e+00`。差异列表与父批一致，无新增偏差；未重冻结，不能称正式 75 全通过。本批没有验证长循环或全部物理配置，也没有测量性能收益。

## 可追溯证据

`output/m4_bonded_reuse_20261024/` 保存 runner、专项脚本、提取的父版函数及身份、命令、源码差异、运行前后身份、日志/退出码、两份比较、receipt 与 SHA-256 清单。源码、依赖、入口、专项输入及保护项运行前后身份一致。

```text
python -X utf8 output/m4_bonded_reuse_20261024/verify_m4.py "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
```

- candidate_60s.log：`b79aaf10f95e017db301ac59f8ca3d556622aec10ec0c6373fda4a599bdea20e`
- probe_bonded.log：`48333e5b1413efe6059c2f37d01e09e9e710a6b4fecb95e0b00bdc3b1cc88df8`
- receipt.json：`fe85f010200d272f256b4cfec40024ed1fa9498afa41d33cbed163975e8a0a19`

状态：本批 bonded 精确对照与父版文字回归通过，按用户授权提交；正式 75 的既知差异继续保留。
