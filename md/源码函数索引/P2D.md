# P2D.jl

- **源文件**: `src/P2D.jl`
- **行数**: 287 行
- **函数/struct 计数**: 5 个独立函数；0 个 struct
- **职责**: Pseudo-two-Dimensional (P2D / Doyle-Fuller-Newman) 电化学模型主体，求解颗粒扩散 + 电解液扩散 + 电极/电解液电势（电荷守恒）；`P2D_potentials` 用迭代法耦合电势方程
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `P2D(case::Case, yt, t; jacobi)` — L1-L47

全局 P2D 求解入口，返回 `(M, K, F, variables, phi_new)`。

- 先调 `P2D_variables`，若 `mechanicalmodel == "full"` 再调 `Mechanicaloutput` 取 `theta_Mn/theta_Mp`
- 颗粒扩散：`ElectrodeDiffusion`；电极电势：`ElectrodePotential`；电解液扩散/电势：`ElectrolyteDiffusion`/`ElectrolytePotential`
- 时间尺度归一化 `M .*= ts_*/t0`
- `jacobi == "constant"` 分支：复用 `param.{NE,PE}.{M_d,K_d,M_p,K_p}`
- 调 `P2D_potentials` 迭代求解电势；再调 `P2D_mass_BC` 构造源项
- 跨文件依赖：`P2D_variables`、`Mechanicaloutput`、`ElectrodeDiffusion`、`ElectrodePotential`、`ElectrolyteDiffusion`、`ElectrolytePotential`、`P2D_potentials`、`P2D_mass_BC`

### `P2D_mass_BC(case, variables)` — L49-L74

构造质量方程（颗粒扩散 + 电解液扩散）的源项 `F = [flux_np; flux_pp; flux_el]`。

- 颗粒表面通量：`flux_np[v_np_surf] = -j_n·rs^2`（用 `case.index` 定位表面 DOF）
- 电解液扩散：高斯点电流 `j_n_gs/j_p_gs`，`aj_el_gs = [j_n·as; 0; j_p·as]`，`coeff *= (1-t⁺)`，`Assemble1D` 装配
- 跨文件依赖：`Assemble1D`

### `P2D_charge_BC(case, variables)` — L76-L118

构造电荷守恒方程的三个源项 `(flux_ne, flux_pe, flux_elc)`。

- 负极电荷源：`coeff_ne = weight·detJ·as·j_n_gs`，`Assemble1D` 后 `flux_ne[1] += -I_app`（**注意**：该值随后在 `P2D_potentials` L160 被覆写为 0，注释 L92 标注 "wiped and not used"，见 [PLACEHOLDER]）
- 正极电荷源：同上，`flux_pe[end] += I_app`（同样被覆写）
- 电解液电荷源：含扩散-迁移耦合项 `kappa_D_eff`，涉及 `dlnf_dlnc(ce)`、`dcedx_gs`、Bruggeman 校正 `tau_el`
- 跨文件依赖：`Assemble1D`

### `P2D_potentials(case, yt, t, K_pot, variables)` — L120-L207

Newton 型迭代求解电极/电解液电势（耦合 Butler-Volmer）。返回 `(phi_new, variables)`。

- 直接强制边界条件：`K_pot[1,1]=-1`、`K_pot[ne+pe+1, ...]` 等多处置零（L130-L135）
- 迭代循环 `for i = 1:iter_max`（`iter_max=100`，`rel_tol=1e-9`，硬编码）
  - 解相对电势 `phi_new_rel = -K_pot \ F_pot`，提取 `phis_n_rel/phis_p_rel/phie_rel`
  - 用 B-V 显式积分算参考电势 `Ve`、`Vp`（L182-L183）
  - 绝对电势 `phi_new = [phis_n_rel; phis_p_rel.+Vp; phie_rel.+Ve]`
  - 收敛判据 `error_j = ‖j_new - j_old‖/‖j_old‖ < rel_tol`
- 达到 `iter_max` 时 `print(...)` 警告（见 [DEBUG]）
- 跨文件依赖：`IntV`、`P2D_variables`、`Mechanicaloutput`、`P2D_charge_BC`

### `P2D_variables(case, yt, t)` — L209-L287

计算 P2D 中间变量并填入 `StandardVariables`。

- 状态提取：从 `yt` 按 `case.index` 取浓度/电势
- 高斯点插值：`cs_gs/ce_gs/phis_gs/phie_gs`
- B-V 方程（高斯点与节点两套）：`j0 = k·Arrhenius·sqrt(cn·|1-cn|·ce)`，`j = 2·j0·sinh(η/2T)`
- 熵热项：`u = U(cs) + (T-T0)·dUdT(cs)`
- 端电压：`V = phis_p[end] - phis_n[1]`
- 跨文件依赖：`StandardVariables`、`Arrhenius`

## 省略项

无（本文件全部 5 个函数均含实质逻辑）。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L203 | `print("maximum iteration mumber has been reached for solving potential equations \n")` | `P2D_potentials` 迭代达上限时打印；属非结构化诊断输出，未升级为 `@warn`/`@error`，调用者无法感知收敛失败 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L92 | `flux_ne[1] += - I_app # this item is wiped and not used`；L99 `flux_pe[end] += I_app # this item is wiped and not used` | [UNCERTAIN] 该加法结果随后在 `P2D_potentials` L160/L161 被覆写为 `0.0`/`Vp0`，属已知死代码；注释自承"未使用"，可能遗留自旧版边界处理，应清理或明确语义 |
| L121 | `iter_max = 100;`；L122 `rel_tol = 1e-9`（硬编码收敛参数） | [UNCERTAIN] 魔数无注释说明取值依据；应提升为 `case.opt` 字段或函数参数以便调参 |

### [COMPLEX-CHECK]

无（迭代循环为单层 `for`，无 ≥3 层嵌套；`&&` 链均 ≤2 条件）。

---
