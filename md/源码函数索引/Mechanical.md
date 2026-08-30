# Mechanical.jl

- **源文件**: `src/Mechanical.jl`
- **行数**: 325 行
- **函数/struct 计数**: 6 个独立函数
- **职责**: 颗粒扩散应力（粒子尺度）；层分辨宏观热-扩散应力（极片尺度，2D 平面应力）：共享恢复核、特征应变、耦合在线收割、固体按需求解
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`md/15_颗粒与极片模量区分.md`（CLAUDE.md §9.4）

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `Mechanicaloutput(case, variables) -> variables` — L1-L110

颗粒尺度扩散应力计算与电化学-力学耦合修正主入口。

- **SPM/SPMe 分支**（L3-L29）：
  - 调 `Calstressdisp(param.NE, mesh_n, c_n, T)` / `Calstressdisp(param.PE, mesh_p, c_p, T)`
  - 应力对过电位修正：`eta_p_new = eta_p - (2/3)·σ_θ_p·Omega_PE`（L14，PE 同理 NE）
  - 端电压耦合：`V_cell_new = V_cell - (2/3)·σ_θ_p·Omega_PE + (2/3)·σ_θ_n·Omega_NE`（L16）
- **P2D/sP2D 分支**（L30-L108）：
  - 每电极节点对应一个颗粒，逐节点调 `Calstressdisp`（L67-L76）
  - 高斯点应力由节点应力插值：`sum(gs.Ni · σ_θ_surf[element[gs.ele,:]])`（L81-L82）
  - 修正后重算界面电流密度：`j_n = j0_n · sinh(0.5·eta_n_new/T) · 2.0`（L79）
- 写入键（部分）："negative/positive particle center radial stress"、"negative/positive particle surface tangential stress"、"negative/positive particle surface displacement"、"negative/positive particle concentration at gauss point"、"negative/positive particle stress coupling diffusion coefficient"、"negative/positive electrode overpotential"、"cell voltage"
- 跨文件依赖：`Calstressdisp`、`PickElement`

### `Calstressdisp(electrode, mesh, cs, T) -> (stress_r_center, stress_theta_surf, disp_surf, theta_M, cs_gs)` — L112-L139

球形颗粒扩散应力解析解（颗粒尺度，CLAUDE.md §9.4：颗粒模量 `electrode.E`）。

- 输入：`electrode::Electrode`（含 `rs`、`nu`、`Omega`、`E`，L126-L129）、`mesh::Mesh`、`cs::Vector`（颗粒内锂浓度）、`T`
- 高斯点浓度：`cs_gs = sum(mesh.gs.Ni · cs[mesh.element[mesh.gs.ele,:]])`（L132）
- 平均浓度（球体积积分）：`cs_av = (3/(4π·rs³)) · IntV(cs_gs · 4π·r²)`（L133）
- 径向中心应力：`σ_r_center = (2·Ω·E·(cs_av - cs_center)) / (9·(1-ν))`（L134）
- 切向表面应力：`σ_θ_surf = (Ω·E·(cs_av - cs_surf)) / (3·(1-ν))`（L135）
- 表面位移：`u_surf = Ω·rs·cs_av/3`（L136）
- 应力-扩散耦合系数：`theta_M = 2·E·Ω² / (T·9·(1-ν))`（L137）
- 跨文件依赖：`IntV`

### `recover_bulk_stress(node, element, material_type, u, ε0, param)` — L151-L179

共享应力恢复核：逐 Q4 单元 `q4_center_gradients` 求 ε(B·u)，逐层平面应力
`σ = D·(ε − ε₀[1,1,0])`（D 由 `moduli_of(param, material_type[e])` 给出，σ_czm 归一空间），
返回 σ_xx/σ_yy/σ_xy/σ_vm。耦合收割与固体工具两流程共用，单一应力定义。

### `macro_eigenstrain(case, variables, T_nodes)` — L182-L192

逐力学体单元特征应变（与在线 CZM 热化学载荷同源，`eigenstrain_of(param, mt)` 分层
计算：α=该层 `alphaT`、β=该层 `Ω/3`；dT 按父热单元，Δsoc 经 `compute_czm_strain_inputs`，
电极膨胀只作用于本层涂层，集流体/隔膜仅热应变）。2026-08-29 α/β 分层化后不再依赖
`czm_param_cache`，旧 PE_PCC 占位与跨界面统一 α_eff 已消除。

### `export_macro_stress(case, variables, variables_hist, v, T_nodes)` — L200-L212

耦合流程（opt.czm_enabled=true）求解中收割：以 `case.czm_layout.u_prev`（本步 CZM
收敛位移）与当步载荷恢复层分辨应力，直写第 v 步历史列
`"diffusion stress xx/yy/xy/vonMises"`；`Solve.jl` 在 CZM 更新块内调用，
`PostProcessing` 导出 `result["diffusion stress xx/yy/xy/vonMises [Pa]"]`（×scale.σ_czm）。

### `thermal_diffusion_stress_2D(case, variables) -> new_variables` — L233-L325

仅固体力学流程（opt.czm_enabled=false）的层分辨宏观应力工具函数：在
`czm_submesh.mesh_bonded`（Φ 合并、无内聚力单元）上按逐层刚度与逐层特征应变
求解线性弹性静力平衡；仅显式调用时求解，不在正常求解流程中。

- **入口断言**（L235）：`PE.E_coat > 0 && NE.E_coat > 0`，缺失 `E_coat` 时阻断宏观力学分析
- **域**（L238-L240）：`czm_submesh.mesh_bonded`（自带 gs）；载荷经 `macro_eigenstrain`
- **本构/装配**（L251-L289）：逐高斯点材料 `moduli_of(material_type[ele_of_gp[g]])`
  → D11/D12/D33；`Assemble×4` 装配 K、`Assemble1D×2` 装配 F（系数含 ε₀·(1+ν)·wJ）
- **边界条件**（L292-L302）：外圈固定；内圈按 `opt.czm_fix_inner`（默认 true=内外均固定），
  相对罚 `penalty = 1e6·max(|diag(K)|)`
- **求解**（L305）：`U = K_mech \ F_mech`
- **应力恢复**（L307-L309）：共享核 `recover_bulk_stress`
- **输出**（L311-L321）：有量纲键 "diffusion stress xx/yy/xy/vonMises [Pa]"（×scale.σ_czm）、
  "displacement x/y [m]"（×scale.L）；不再产出 "thermal stress vonMises" 等死键
- 跨文件依赖：`compute_czm_strain_inputs`、`compute_czm_params_per_interface`、`moduli_of`、
  `Assemble`、`Assemble1D`、`identify_boundary_nodes`、`q4_center_gradients`

## 省略项

无。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info`；`@warn`（L297）为求解失败告警，属于运行时状态而非调试输出。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L273 | `penalty = dmax_bc > 0 ? 1e6 * dmax_bc : 1e12`（矩阵全零时回退 1e12） | 与 `czm.jl` L823 同模式：相对罚正常路径，回退固定值 1e12 主导条件数；正常路径不会触发 |
| L297 | `@warn "Mechanical solve failed, using zero displacement" e`（求解失败静默回退零位移） | 失败被静默吃掉：用户得到零位移场但 `variables` 仍写入，下游分析无法区分"零应力"与"求解失败"；建议返回 `converged` 标志或写入诊断键 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L67 | P2D 分支中两层嵌套循环（电极节点 × 颗粒内部节点，L67-L76），含 `Int64.(collect((i-1)*meshnum_perparticle_n .+ (1:meshnum_perparticle_n)))` 索引计算 | 抽出 `extract_particle_mesh(mesh_full, i, n_per_particle) -> Mesh` helper |
| L167 | `@assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0 "..."`（2 个 `&&` 条件 + 长 error message） | 简单入口断言，可保留；或拆为两条独立 `@assert` 分别给出更精确的错误信息 |
| L266 | `if is_inner[i] \|\| is_outer[i]` + 内部 `bc_nodes[i] = :fixed_xy`（嵌套 2 层，配合 2 个 `||` 谓词） | 简单循环，可保留 |
