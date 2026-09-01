# CZM 最大分离位移核查发现

## Requirements

- 核查最大法向分离位移是否计算正确。
- 运行单元级脚本并分析，而不是仅引用 `testexample.jl` 汇总输出。
- 区分单元分离公式、局部标架、时间历史归约和有量纲还原。
- 本轮只诊断；若发现缺陷，先报告根因，不擅自修复。

## Research Findings

- 当前完整示例在 `fix_inner=true` 与新增 22 个分层端点约束下输出 `CZM 最大法向分离位移 = 1.0667e-12 m`。
- 修改边界前同一示例输出 `7.4202e-13 m`；应先定位数据流和最大值来源，不能把数值变化直接解释为错误。
- 历史记忆指出当前权威入口应为 `compute_separation(czm_mesh, elem, u)`，法向由 host-inner 指向 host-outer；仍需以当前源码和实跑验证。
- 初步源码索引显示两条不同但关联的数据：`CZMResult.separation_n` 是当前步各 cohesive 单元的归一化分离；最终 `czm δ_max_n [m]` 来自 `DamageState.δ_max_n` 的历史最大值，而不是直接对最终步 `separation_n` 求最大。
- `CallModel.jl` 每步对全部 `DamageState.δ_max_n` 取空间最大，`PostProcessing.jl` 再乘一次 `scale.δ_czm` 还原为米；需继续验证状态更新与尺度是否一致。
- `result["czm separation normal [m]"]` 是逐单元当前分离历史，理论上可作为 `czm δ_max_n [m]` 的独立交叉来源，但二者语义分别为当前值与不可逆历史峰值。
- `compute_separation` 对 2 点线积分使用线性形函数，并返回单元平均跳跃；对于线性插值，该平均值应等于界面两端相对位移的算术平均投影，可独立手算。
- 局部旋转矩阵 `R=[nᵀ;tᵀ]`；法向先由 bottom 边切向左转得到，再依据 host-inner 到 host-outer 的质心向量必要时翻转，因此正法向语义明确为由内向外。
- 分离矩阵的节点排列为 bottom `(n1,n2)`、top `(n3,n4)`，其中 top 从 `elem.nodes_top=(n_lo_copy,n_hi_copy)` 解构为 `n4,n3` 后按 `[n3,n4]` 进入 `u_e`；需用受控端点非均匀位移测试确认没有次序互换。
- 当前未提交的 `CzmSolve.jl`/`CouplingState.jl` 改动只增加 basic 线性系统因子缓存及工作区字段，没有改变分离计算、状态更新或输出缩放链。
- `bilinear_traction_state` 仅在 `δ_eff > damage_state.δ_max_eff` 时更新 `δ_max_n=max(old, max(δ_n,0))`。对混合模式一般历史而言，这并不严格等价于逐时刻独立取 `max(δ_n,0)`：先出现较大切向有效分离后，较小有效值但更大法向值可能不会写入 `δ_max_n`。需通过单元序列测试判断这是潜在语义缺陷还是本算例未触发的边界情况。
- `CZMResult.separation_n` 与 `DamageState.δ_max_n` 均处于 `scale.δ_czm` 归一空间；后处理分别乘一次同一尺度，没有从静态链路看到重复缩放。
- 损伤状态只在求解收敛后提交；历史峰值摘要由 `CallModel` 读取已提交的 `case.mech.damage_states`。
- 静态调用链已确认两个输出时序问题：① `CallModel` 先用上一次已提交状态生成 `czm δ_max_n`，`Solve` 后调用 `update_czm_damage!`，却未在 `Variable_update!` 前刷新该摘要，所以记录历史落后一次 CZM 更新；② `update_czm_damage!` 只返回 `CZMResult`，`Solve` 没有调用 `czm_output_to_variables`，所以逐单元的位移/分离/牵引变量没有由本步求解结果写回。需用端到端快照确认实际数值影响。
- 现有三组单测全部通过：`test_create_czm_mesh.jl` 80,747 断言、`test_czm_scale_redesign.jl` 28 断言、`test_bilinear_per_interface.jl` 71 断言。
- 自定义受控单元探针在 PE–PCC、NE–NCC 和最后一个 cohesive 单元上均得到：非均匀端点跳跃的法向平均 `2.5e-9`、切向平均 `9.0e-10`，与独立算术平均投影一致；物理法向分离均为 `4.32e-13 m`。刚体平移在所有 cohesive 单元上给出严格 `0.0`。
- 单元探针同时确认 `compute_separation` 返回位移空间（以 `L` 归一）的跳跃；装配返回分离空间（以 `δ_czm` 归一）的值，两者满足 `raw·L = assembled·δ_czm`。
- 非比例混合模态序列实证一般性语义缺口：先施加归一切向分离 `4.429512516469e-4`，再施加法向分离 `2.657707509881e-4`时，后者小于旧 `δ_max_eff`，状态的 `δ_max_n` 仍为 `0.0`，而按字面的法向历史峰值应为 `2.657707509881e-4`。
- 用户更正边界契约：`fix_inner=false` 也必须叠加分层端点。已删除早返回，使 `fix_inner` 只控制内圈是否固定；两种模式共用同一端点集合和两个排除点。
- 实际 nθ=80 计数：`fix_inner=false` 原外圈 80 点，更新后共 103 点，净增 23；`fix_inner=true` 现有测试确认相对内外圈净增 22。
- 用户要求将 `example/testexample.jl` 从 3600 s 改回 60 s；诊断探针同步为 60 s。
- 60 s、nθ=80、`mix`、`fix_inner=false`+新增 23 端点、预应力关闭的端到端快照结果：真实全时空法向峰值为 `1.101106737768e-12 m`，最大值位于 `t=60 s`、cohesive 单元 5127、`:NE_NCC`，bottom 节点 `[10179,10180]`、top 节点 `[18674,18675]`。
- 对该最大单元，原始 `CZMResult` 快照为 `1.101106737768e-12 m`，不复用被测装配实现的节点端点平均投影手算为 `1.101106737769e-12 m`，`compute_separation·L` 也为 `1.101106737769e-12 m`；切向三路亦一致为约 `-8.5797670117e-14 m`。说明单元节点顺序、局部标架、平均积分和尺度换算正确。
- 但 `result["czm δ_max_n [m]"]` 的峰值/末值只有 `8.730984248720e-13 m`，比真实峰值少 `2.280083128963e-13 m`（约 20.7%）。整条记录与“累积峰值滞后一次 CZM 更新”完全一致：同步误差为 `2.280083128963e-13 m`，滞后一步的序列误差为严格 `0.0`。
- `result["czm separation normal [m]"]` 全历史最大绝对值为严格 `0.0 m`，证实 `Solve` 没有把返回的 `CZMResult` 写回逐单元分离历史。
- 正式 `example/testexample.jl` 退出码 0，19 个记录时间层，公开文字输出同样报告 `8.7310e-13 m`；该值不是本次 60 s 运行的真实全时空最大法向分离。

## Technical Decisions

| Decision | Rationale |
|---|---|
| 以独立公式复算，不复用被测 helper 构造期望 | 防止测试与实现同错 |
| 记录最大值的单元、时间、节点和位移 | 让端到端汇总值可追溯 |
| 分别检查归一化值与物理量值 | 防止 `scale.δ_czm` 重复或遗漏 |

## Issues Encountered

| Issue | Resolution |
|---|---|
| 暂无 | — |

## Resources

- `src/Tools.jl`
- `src/CzmSolve.jl`
- `src/CouplingState.jl`
- `src/Variables.jl`
- `src/PostProcessing.jl`
- `test/test_create_czm_mesh.jl`
- `test/unit_czm_bilinear.jl`
- `example/testexample.jl`

## Visual/Browser Findings

- 本任务不使用视觉或浏览器材料。
