# 层分辨应力求解：发现与决策

## Ω_PE 符号决策（2026-08-29，用户拍板）

- 实测探针（1C 放电 60 s）：soc_n 0.9014→0.8934（脱锂）、soc_p 0.270→0.2752（嵌锂）。
- 公式 ε₀ = (Ω/3)·Δsoc 下：充放电 Δsoc_n/Δsoc_p 必反号；Ω_PE<0 再翻一次 → PE/NE 特征应变永远**同号**（旧参数下同缩同胀）；"一缩一胀"需要 **Ω_PE>0**。
- 用户选定 `PE.Omega = +7.88e-7 m³/mol`（更正 src/parameters/Jellyroll.jl，原 −7.28e-7 为未审 placeholder）。
- 验证图案：放电 NE 涂层 +0.6 MPa（拉，脱锂收缩）/ PE 涂层 −0.2 MPa（压，嵌锂膨胀）；耦合 nθ=80：NE +0.23~+0.33 / PE −0.10~+0.05；集流体承担约束拉力（内圈集中 +2~+4 MPa，外圈转入 ~−1 压缩）。
- 连带影响与处理：test_czm_j2_integration 两个用例标定于旧 β 反号（设计失配 ~1.5e-2），新符号下同向 Δsoc 失配减半（~5e-3）且重解触发 return_mapping 的 Δγ=−2.3e-21 硬报错 → 用例改为放电向 Δsn=−0.3/Δsp=+0.3 恢复 ~1.5e-2 失配；C2-lite 仅改注释（5e-3 仍超屈服）。czm δ_max_n 基线由 1.5174e-12 → 7.0037e-13 m（v5 重立）。

## 用户要求与决策链

- 实现 per-interface（层分辨）应力求解；彻底删除原先粗网格的应力求解，直接完全替换。
- 尺度转换关切已核实解决（见下）；架构经多轮收敛：
  1. 不新建独立装配求解器（避免与在线系统重复、第二套应力定义）；
  2. 耦合流程从整体系统导出，且在求解过程中直接导出；
  3. 导出不放进 `update_czm_damage!`（内聚力独有步骤），放公共管线位置；
  4. czm_enabled=false 走仅固体力学流程，且作为工具函数按需调用，正常流程不经过；
  5. 不新建门控选项（固体力学单向、不影响热学，网格存在即门）；
  6. 固体工具函数保留原名 `thermal_diffusion_stress_2D`；新函数无 `_` 前缀、无 `!` 结尾；
  7. 不新增防御性断言（保留旧函数首行 E_coat 断言），避免代码臃肿。

## 尺度转换核实结论

- 耦合路径应力空间统一已解决：`moduli_of`（src/czm.jl:47-64）将 E_coat 归一的体模量乘 `scale.E_coat/scale.σ_czm` 转入 σ_czm 空间（重设计 v2 §3"双重再缩放"）。
- 位移/分离空间换算已解决：`Λ = scale.L/scale.δ_czm`（CouplingState.jl:363-365；czm.jl:178-249，切线刚度乘一次 Λ，重设计 v2 §5）。
- 旧路径装配端自洽（E_eff 同在 σ_czm 空间）；**遗留 bug 在还原端**：两个示例按 `scale.E_coat` 还原（jellyroll_stress_displacement.jl:206-213 注释声称 E_eff 按 E_coat 归一，过时错误），应力被放大 `E_coat/σ_czm ≈ 7.4e9/82e6 ≈ 90` 倍——旧环向应力图 30–94 MPa 的真实量级为 O(1) MPa。新实现输出在函数内部 ×scale.σ_czm → Pa、位移 ×scale.L → m。

## 关键代码事实

- `CzmSubmesh.material_type::Vector{Symbol}`（:PE/:PCC/:SP/:NE/:NCC 每 bulk 单元）与 `thermal_elem_map` 齐备（SetMesh.jl:42-50）；bulk 单元按 layer-outer/segment-inner 排序：e = (layer-1)*n_segments + seg。
- `czm_submesh.mesh_bonded` 为 Φ 合并完整 Mesh、自带 gs（Jellyrollmodel.jl:628-629），元素行序与 material_type 一致，无界面复制节点——固体分支求解域。
- `compute_czm_strain_inputs(case, variables, T_nodes)`（CouplingState.jl:461-549）逐层分发：PE 层只拿 Δsoc_p、NE 层只拿 Δsoc_n、PCC/NCC/SP 为 0；dT 按父热单元均值；接受向量或历史矩阵（取末列）。
- `moduli_of(param, mt)` 逐层平面应力模量（σ_czm 空间）。
- `SP/PCC/NCC.alphaT` 未在 Jellyroll.jl 设置（默认 0）；在线热载荷用跨界面 α_eff = :PE_PCC.α（CouplingState.jl:582）——新实现沿用，不伪造逐层 α。
- **执行顺序**：`update_czm_damage!` 在 Solve.jl 主循环（272 行）调用，位于 `Variable_update!` 写历史列（235 行）**之后**；`case.czm_layout.u_prev = result.displacement` 在 CouplingState.jl:696 更新为本步收敛值。⇒ 收割必须放 Solve.jl CZM 更新块内（放 CallModel 尾部会拿到滞后一步的位移），直接写 `variables_hist[...][:, v]` 列。与批准计划"CallModel 尾部"细节不同、意图一致（管线层、不进 update_czm_damage!）。
- `Variables.jl:127` 的 CZM 历史块已以 `opt.czm_enabled == true` 为外层门（138 行 mesh 条件嵌套在内）——计划中的"收紧"已天然存在，czm-off+建网格算例不会分配 CZM 死键。
- **既有缺口**：`czm_output_to_variables`（CzmPostProcess.jl:99）无任何调用点，"czm displacement/damage/traction/separation" 历史从未被写入（全零导出）。不在本批范围，仅记录。
- 旧函数调用方仅 3 处（testexample.jl:273、jellyroll_stress_displacement.jl:205、test_electrode_coat_modulus.jl:120/146），在线求解器与 test/ 零依赖；`q4_center_gradients`（Tools.jl:214-226）返回 (dNdx, dNdy, detJ)。
- Assemble/Assemble1D（Assemble.jl）：coeff 为逐高斯点向量（含 wJ 与材料因子）；mesh.gs.dNidx 列 1:4 = dN/dx、5:8 = dN/dy。
- 旧键 "thermal stress vonMises"、"diffusion stress vonMises only" 全仓零消费者。
- testexample 冻结基线（Simplify/baseline/testexample/metrics.toml）只含求解指标 + 旧 testexample_results.png；三张 final_* PNG 不在门内。

## 技术决策

| 决策 | 理由 |
|---|---|
| 耦合收割读 layout.u_prev 而非传 czm_result | 解耦数据来源；Solve.jl 调用点在 CZM 更新块内，u_prev 即本步收敛值 |
| 收割仅在 czm_update_interval 命中步执行 | 保证 (u, 载荷) 同步配对，避免间隔步载荷前进而位移滞后的错配 |
| 固体分支用 mesh_bonded 而非 czm_mesh.node | 无界面复制节点，Φ 合并连续体，避免无内聚力装配时的悬浮奇异 |
| BC 外圈固定 + 内圈按 opt.czm_fix_inner（默认 true） | 与在线 CZM 同一旋钮；默认值复现旧函数双固定语义 |
| 输出键 "diffusion stress xx/yy/xy/vonMises [Pa]"、"displacement x/y [m]" | 沿用旧词汇 + 单位后缀，两流程同名同单位，调用方零尺度接触 |
| ε₀ 用 α_eff 统一热应变 + 逐层 β·Δsoc | 与在线 CZM 载荷物理一致；SP/PCC/NCC alphaT 未定义，不伪造 |

## 资源

- src/Mechanical.jl、src/Variables.jl、src/Solve.jl、src/PostProcessing.jl
- src/CouplingState.jl（compute_czm_strain_inputs / CzmLayout / compute_czm_params_per_interface）
- src/czm.jl（moduli_of）、src/Assemble.jl、src/Tools.jl（q4_center_gradients、identify_boundary_nodes）
- example/testexample.jl、example/jellyroll_stress_displacement.jl、example/力学模块验证/test_electrode_coat_modulus.jl
- Simplify/baseline/testexample/
