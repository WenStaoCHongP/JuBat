# 层分辨应力求解：进度记录

## 会话：2026-08-29

### Phase 1：规划与事实核实

- **状态：** complete
- 计划经多轮评审收敛（系统导出 → 求解中导出 → 公共管线 → 双流程 + 保留原名 + 断言纪律），已获批准。
- 已核实尺度转换（v2 §3 双重再缩放 + §5 Λ）、执行顺序（update_czm_damage! 在 Solve.jl:272、u_prev 于 CouplingState.jl:696 更新）、历史写入机制（Variable_update! 仅搬运转存、收割需直写 variables_hist 列）。
- 已核实 Variables.jl:127 既有 opt.czm_enabled 外层门（原计划"收紧"不必要）。
- 已发现既有缺口：czm_output_to_variables 无调用点（czm displacement 等历史全零），不在本批范围。

### Phase 2-4：实现、接线与调用方迁移

- **状态：** complete
- src/Mechanical.jl 原址重写（325 行）：recover_bulk_stress / macro_eigenstrain / export_macro_stress / thermal_diffusion_stress_2D（保留原名，首行 E_coat 断言原样）。
- 接线：Variables.jl（czm_enabled+mesh+线性相容时分配 4 个应力历史）、Solve.jl CZM 更新块内收割（u_prev 与当步载荷配对）、PostProcessing 导出 [Pa] 键。
- 调用方迁移：testexample 直读结果键 + 99% 分位对称色标 + 新增环向应力径向剖面图；jellyroll 建力学网格（opt.czm_enabled 保持 false）+ T_nodes 归一化修正 + mesh_bonded 绘图；test_electrode_coat_modulus TEST 7 建网格直调（并修复其 TEST 1/3 的既有陈旧断言——HEAD 上已失效）。

### Ω_PE 参数更正（用户决策，会话中）

- 探针实测（1C 放电 60 s）：soc_n 0.9014→0.8934（脱锂）、soc_p 0.270→0.2752（嵌锂）；结合公式 ε₀=Ω·Δsoc/3 推导：旧 Ω_PE<0 时 PE/NE 特征应变永远同号（同缩同胀），"一缩一胀"需要 Ω_PE>0。
- 用户选定 **PE.Omega = +7.88e-7 m³/mol**（src/parameters/Jellyroll.jl 已改）。
- 验证：放电下 NE 涂层 +0.60 MPa（拉）/ PE 涂层 −0.19 MPa（压）——涂层一拉一压出现；耦合 nθ=80 算例 NE +0.23~+0.33 / PE −0.10~+0.05。
- 连带：test_czm_j2_integration.jl 两个用例标定于旧 β 反号（失配 ~1.5e-2），Ω 更正后失配减半且重解触发 Δγ=−2.3e-21 硬报错；已按新符号重标定（Δsn=−0.3/Δsp=+0.3，恢复 ~1.5e-2 失配），全绿。

### Phase 5-6：验证门

- **状态：** complete（jellyroll 重型示例按用户指示中止，未完成端到端运行；其余全过）
- testexample 完整运行 exit 0，求解指标与 v4 基线按记录精度一致（1682/1763/19 步/4.0367/3.9438/0.0833/298.15–299.00/D=0/断裂 0）；czm δ_max_n 7.0037e-13 m（v4 1.5174e-12，Ω_PE 授权变化）。
- 新层分辨环向应力范围 −1.78~+4.00 MPa（旧图 30–94 MPa 为 ~90× 单位放大）；径向剖面图确认 NE 带 +0.3 / PE ≈−0.05 / 集流体内圈 +2~+4 约束拉力过渡到外圈 −1.1 压缩。
- 全量 test 套件 34/34 通过（test_czm_j2_integration 载荷按新 Ω 符号重标定后）。
- jellyroll_stress_displacement.jl（nθ=360、3600 s 重型算例）后台运行 ~30 min 后按用户指示停止；脚本迁移已完成静态检查与解析验证，端到端运行留待后续需要时执行。

### 后续批次：验证入口拆分（用户指示，2026-08-29）

- **状态：** complete
- 新建 `example/couple_example.jl`：承接原 testexample.jl 全部代码（求解 + 四张绘图），输出按 §9.9 写入 `output/couple_example/`。
- `example/testexample.jl` 精简为纯文字结果输出（60 s 快速门）：保留参数/网格/求解/耗时统计/全部文字指标（含新增应力范围文字指标），移除 Plots 及全部绘图；`rotate_stress_to_polar` 为纯计算保留。
- 验证：60 s 快速门指标与 v5 基线完全一致；couple_example 出图四张 PNG 哈希与 v5 完全一致（拆分行为无损）。
- 基线 v6 重立（testexample-20260829T231308+0800）：快速门文字指标 + 图片门移至 couple_example 按需运行；源码清单 47 文件。
- AGENTS.md §9.6 更新为两级验证流程（受影响验证 + 60 s 快速门；涉绘图才跑全量）；§10 示例表补 couple_example.jl。

## 测试结果

| 检查 | 预期 | 实际 | 状态 |
|---|---|---|---|
| 语法解析（4 src + 4 脚本） | PARSE_OK | 全部 PARSE_OK | pass |
| 新测试 test_layer_resolved_stress.jl | 5 testset 全过 | 全过（Ω 更正后重跑） | pass |
| test_electrode_coat_modulus.jl | ALL TESTS PASSED | ALL TESTS PASSED | pass |
| testexample 基线比对 | 求解指标一致 | 一致；仅 czm δ 变化（授权） | pass |
| 视觉检查径向剖面 | 层间拉压交替可见 | NE +0.3 带 / PE ≈0 压 / 箔材内拉外压 | pass |
| test_czm_j2_integration.jl | 4 testset 全过 | 全过（载荷重标定后） | pass |
| 全量 test/runtests.jl | 34 文件全过 | 首跑 33/34（J2 标定问题已修）；复跑 34/34 全绿（8m56s） | pass |
| 60s 快速门 testexample.jl（拆分后） | 文字指标与 v5 一致 | 全部一致 exit 0 | pass |
| couple_example.jl（拆分后） | 四张 PNG 哈希与 v5 一致 | 完全一致 exit 0 | pass |
| jellyroll_stress_displacement.jl | exit 0 出图 | 用户指示中止（重型算例 ~30 min）；静态/解析检查通过 | stopped |
| 静态 grep 旧键/旧调用 | 零残留 | 零残留 | pass |
| git diff --check（本批涉及文件） | 干净 | 干净（仅既有 CRLF 提示） | pass |

## 错误日志

| 时间 | 错误 | 尝试 | 处理 |
|---|---|---:|---|
| 2026-08-29 | runtests.jl 需 --project=.（JuBat 包路径） | 1 | 加 --project=. 重跑 |
| 2026-08-29 | test_electrode_coat_modulus TEST 1/3 陈旧断言（HEAD 已失效） | 1 | 按当前参数值/全叠合公式修正 |
| 2026-08-29 | test_czm_j2_integration 标定于旧 Ω_PE 反号；Δγ=−2.3e-21 硬报错 | 2 | 载荷按新符号重标定（Δsp=+0.3） |
| 2026-08-29 | 视觉模型对环形云图色标读数不可靠（与代码 clims 矛盾） | 2 | 以数值统计+径向剖面分析为准 |
