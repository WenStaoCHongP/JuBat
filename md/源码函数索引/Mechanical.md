# Mechanical.jl

- **源文件**: `src/Mechanical.jl`
- **行数**: 360 行
- **函数/struct 计数**: 3 个独立函数
- **职责**: 颗粒扩散应力（粒子尺度）与宏观热-扩散应力（极片尺度，2D 平面应力 FEM）计算；输出应力/位移/耦合扩散系数到 `variables`
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

### `thermal_diffusion_stress_2D(case, variables) -> new_variables` — L165-L360

宏观 2D 平面应力 FEM 求解热 + 扩散应力（极片尺度，CLAUDE.md §9.4：涂层模量 `E_coat`）。

- **入口断言**（L167）：`PE.E_coat > 0 && NE.E_coat > 0`，缺失 `E_coat` 时阻断宏观力学分析
- **模量来源**（L184-L187）：`compute_czm_params_per_interface(case)` 取 `:PE_PCC` 界面的 `(E_eff, ν_eff, α_eff)`——非 CZM 路径暂用 PE_PCC 占位（注释 L183），per-interface 化由 CZM 路径负责
- **温度/SOC 提取**（L177-L203）：从 `variables["T_nodes"]`、`variables["thermal2D element soc_n/p"]` 计算每单元 `dT_elem` 与 `Δsoc_n/p_elem`
- **本构**（L215-L217）：平面应力 D 矩阵 `D11/12/33 = E_eff / (1-ν²) · [...]`
- **刚度装配**（L236-L241）：`K_uu + K_vv + K_uv + K_vu`，使用 `Assemble`
- **载荷装配**（L258-L260）：`F_u + F_v`，初始应变 `ε_0 = α·dT + β_n·Δsoc_n + β_p·Δsoc_p`（L246）
- **边界条件**（L263-L291）：内/外圈节点固定（`:fixed_xy`），相对罚方法 `penalty = 1e6·max(|diag(K)|)`（L273）
- **求解**（L294-L299）：`K_mech \ F_mech`，失败时静默回退零位移并 `@warn`（L297）
- **应力恢复**（L302-L345）：单元中心梯度 `q4_center_gradients`，按比例分热应力 / 扩散应力（L336-L340）
- **输出**（L348-L358）：转换为有量纲 `·L_ref`，写入键 "displacement x/y"、"diffusion stress xx/yy/xy/vonMises"、"thermal stress vonMises"、"diffusion stress vonMises only"
- 跨文件依赖：`compute_czm_params_per_interface`、`element_nodal_mean`、`Assemble`、`Assemble1D`、`identify_boundary_nodes`、`q4_center_gradients`

## 省略项

无。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info`；`@warn`（L297）为求解失败告警，属于运行时状态而非调试输出。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L183 | `# 获取材料参数：非 CZM 路径暂用 PE_PCC 占位（per-interface 化由 CZM 路径负责）`（注释明示 PE_PCC 占位） | 显式占位：非 CZM 路径下所有单元使用 PE_PCC 界面的 `(E_eff, ν_eff, α_eff)`，不区分 NE_NCC；当前简化但物理上可疑（NE 侧热膨胀与 PE 不同），恢复需逐单元查表 |
| L273 | `penalty = dmax_bc > 0 ? 1e6 * dmax_bc : 1e12`（矩阵全零时回退 1e12） | 与 `czm.jl` L823 同模式：相对罚正常路径，回退固定值 1e12 主导条件数；正常路径不会触发 |
| L297 | `@warn "Mechanical solve failed, using zero displacement" e`（求解失败静默回退零位移） | 失败被静默吃掉：用户得到零位移场但 `variables` 仍写入，下游分析无法区分"零应力"与"求解失败"；建议返回 `converged` 标志或写入诊断键 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L67-L76 | P2D 分支中两层嵌套循环（电极节点 × 颗粒内部节点），含 `Int64.(collect((i-1)*meshnum_perparticle_n .+ (1:meshnum_perparticle_n)))` 索引计算 | 抽出 `extract_particle_mesh(mesh_full, i, n_per_particle) -> Mesh` helper；当前两层嵌套跨 L67-L76 |
| L167 | `@assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0 "..."`（2 个 `&&` 条件 + 长 error message） | 简单入口断言，可保留；或拆为两条独立 `@assert` 分别给出更精确的错误信息 |
| L266 | `if is_inner[i] \|\| is_outer[i]` + 内部 `bc_nodes[i] = :fixed_xy`（嵌套 2 层，配合 2 个 `||` 谓词） | 简单循环，可保留 |
