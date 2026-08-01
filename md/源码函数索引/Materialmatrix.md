# Materialmatrix.jl

- **源文件**: `src/Materialmatrix.jl`
- **行数**: 427 行
- **函数/struct 计数**: 12 个独立函数
- **职责**: 本构矩阵与间隙导热模型——热容量/各向异性导热分层装配、CZM 双线性牵引-分离律（含粘性正则化）、损伤批量更新、间隙导热（平行热路模型）、有效面积缩减因子
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/06_内聚力模型_CZM.md`、`md/07_界面热阻模型.md`、`md/14_粘性正则化.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ) -> Vector{Float64}` — L16-L23

每高斯点的体积热容权重（jellyroll 2D 热装配用）。

- 每单元按层权重 `fks[e,1..5]` 加权求和 NE/SP/PE/PCC/NCC 的 `rho·heat_Q`（L20）
- 网格已无量纲化，直接乘 `wJ`（L22）

### `thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy) -> (k_xx, k_xy, k_yy)` — L30-L57

每高斯点的各向异性导热分量（径向串联、切向并联）。

- 径向（L36-L39）：串联热阻 `1 / Σ(f_k/λ_k)`；`denom > 0` 否则 `error`（L38）
- 切向（L40）：并联 `Σ(f_k·λ_k)`
- 旋转矩阵（L48-L53）：`theta = atan(gy, gx)`，旋转到全局 `(k_xx, k_xy, k_yy)`

### `bilinear_traction_state(δ_n, δ_t, damage_state, params; visc_beta) -> (T_n, T_t, D_eq, new_state)` — L68-L161

双线性牵引-分离律核心计算，返回牵引力与更新后的 `DamageState`。

- **model1**（纯 Mode I，L102-L116）：`δ_eff = δ_n_pos`，切向不损伤（仅 K_t 弹性）
- **mixed-mode**（L107-L116）：`δ_eff = sqrt(δ_n² + δ_t²)`；`β = |δ_t|/δ_eff`；`δ_0/δ_c_eff = sqrt(δ_0/δ_c_n² + (δ_0/δ_c_t² - δ_0/δ_c_n²)·β^η)`
- **加载判定**（L121-L140）：`δ_eff > δ_max_hist` 时计算新 D_eq 并更新 `δ_max_*`
  - `δ_eff ≤ δ_0_eff` → `D_eq = 0`
  - `δ_eff ≥ δ_c_eff` → `D_eq = 1`（断裂）
  - 中间 → `D_eq = δ_c·(δ-δ_0)/(δ·(δ_c-δ_0))`
- **粘性正则化**（L142-L145，`md/14_粘性正则化.md`）：`D_visc = D_visc_old + visc_beta·(D_eq - D_visc_old)`，强制单调 `max(D_visc_old, ...)`
- **牵引力**（L148-L158）：使用 `D_visc`（而非 `D_eq`）；压缩 `δ_n < 0` 时无损伤衰减 `T_n = K_n·δ_n`
- 已断裂（L88-L100）：`D = D_visc = 1`，model1 保留切向 `K_t·δ_t`

### `bilinear_traction(δ_n, δ_t, damage_state, params; update, visc_beta) -> (T_n, T_t, D)` — L163-L175

`bilinear_traction_state` 的薄包装，`update=true` 时将 `new_state` 字段写回输入的 `damage_state`（原位修改）。

### `bilinear_tangent(δ_n, δ_t, damage_state, params; visc_beta) -> dT_dδ::Matrix{Float64,2}` — L182-L286

双线性切线刚度矩阵（2×2）。

- 内部复现 `bilinear_traction_state` 的 D_eq / D_visc 计算逻辑以保证一致性（L220-L233，注释明示）
- 加载判定（L222-L231）：`is_loading = δ_eff > δ_max_hist - 1e-15`（容差允许数值噪声）
- model1（L235-L256）：切向 `dT_dδ[2,2] = K_t`（不损伤）
  - 弹性段 / 卸载：`(1-D_visc)·K_n` 或 `K_n`（压缩）
  - 软化段（L247-L250）：`(1-D_visc)·K_n - K_n·δ_n·visc_beta·dD/dδn`，**关键**：dD/dδ 乘 visc_beta 一致线性化
  - 完全断裂：`1e-10·K_n`（避免奇异）
- mixed-mode（L257-L283）：含 `dδeff/dδn` / `dδeff/dδt` 链式求导
- 跨文件依赖：无（纯本构）

### `update_damage(damage_states, separations, params; visc_beta) -> new_states` — L293-L307

批量更新损伤状态，逐个调 `bilinear_traction_state` 取 `new_state`。

- 断言长度匹配（L295）；类型校验 `state isa DamageState`（L300）

### `compute_gap_conductance(D, δ_n, params) -> h_eff` — L329-L354

间隙导热系数（平行热路模型，`md/07_界面热阻模型.md`）。

- **单位契约**（重设计 v2，docstring L324-L327）：分离量 ÷Λ（×δ_czm/L）转换到 L 空间后再运算
- `inv_Λ = 1 / params.Λ`（L337）；`delta = max(δ_n, 0)·inv_Λ`
- `D_clamped = clamp(D, 0, 0.9999)`（L341，避免 D=1 时接触导热完全消失）
- 三段（L344-L350）：
  - `delta < delta0`：`h_c0 + k_air/(delta + 2βλ_m)`（无损伤衰减）
  - `delta0 ≤ delta < threshold`：`h_c0·(1-D) + k_air/(...)`（接触衰减）
  - `delta ≥ threshold`：`h_c0·(1-D) + k_air/(delta + delta0)`（gap 主导）
- `h_eff > 0` 否则 `error`（L352）

### `compute_element_gap_conductance(czm_mesh, elem_idx, params) -> h_eff` — L359-L364

`compute_gap_conductance` 的 cohesive mesh 包装：从 `damage_states[elem_idx]` 取 `(D, δ_max_n)`。

### `get_fractured_elements(czm_mesh) -> Vector{Int64}` — L369-L377

返回已断裂（`fractured || D >= 0.99`）的 cohesive 单元索引列表。

### `get_active_elements(czm_mesh, mesh_data) -> Vector{Int64}` — L382-L400

返回未被 CZM 失效隔离的热单元列表。

- 遍历内层热单元（`mesh_data.is_inner_layer`），检查其关联 cohesive 是否全断裂
- 跨文件依赖：`mesh_data.is_inner_layer`、`mesh_data.czm_element_map`

### `compute_all_gap_conductances(czm_mesh, params) -> Vector{Float64}` — L405-L412

所有 cohesive 单元的 `compute_element_gap_conductance` 批处理。

### `effective_area_factor(D, D_threshold) -> Float64` — L424-L427

热单元有效面积比例因子（CLAUDE.md §5.3 `czm_area_loss_*`）。

- `D ≤ D_threshold` → `1.0`（无缩减）
- `D > D_threshold` → `(1 - D) / (1 - D_threshold)`（线性缩减至 D=1 时为 0）

## 省略项

无。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info` / `@warn`；`error(...)`（L38, L245 of `moduli_of`）为参数验证，不计入 DEBUG。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L113 | `δ_0_eff = δ_0_n` / `δ_c_eff = δ_c_n`（mixed-mode 当 `δ_eff ≤ 1e-15` 时退化值，对应 `else` 分支 L113-L115） | 退化分支：零分离时无法计算 `β`，使用纯法向值合理；数值噪声下可能切换分支 |
| L194 | `dT_dδ[1, 1] = 1e-10 * K_n`（完全断裂时法向切线刚度回退，避免奇异） | 注释解释用途；1e-10 是 magic number，过小可能影响 Newton 矩阵条件数 |
| L198 | `dT_dδ[2, 2] = 1e-10 * K_t`（mixed-mode 完全断裂时切向切线刚度） | 同 L194 |
| L215 | `δ_0_eff = δ_0_n` / `δ_c_eff = δ_c_n`（mixed-mode 当 `δ_eff ≤ 1e-15` 时退化值，在 `bilinear_tangent` 内） | 与 L113 同模式，一致性要求；两处独立硬编码，可能不同步 |
| L245 | `dT_dδ[1, 1] = 1e-10 * K_n`（model1 完全断裂法向切线刚度） | 与 L194 同值（model1 在另一个分支） |
| L266 | `dT_dδ[1, 1] = 1e-10 * K_n` / `dT_dδ[2, 2] = 1e-10 * K_t`（mixed-mode 完全断裂时） | 与 L194/L198 重复，可抽出 `fracture_tangent(K, scale=1e-10)` helper |
| L341 | `D_clamped = clamp(D, 0.0, 0.9999)`（gap conductance 中避免 D=1 时接触完全消失） | magic number 0.9999：保留 0.01% 接触导热避免数值奇异；物理上断裂界面仍有少量接触 |
| L352 | `h_eff > 0 \|\| error(...)`（gap conductance 非正检查） | 参数验证，非占位 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L20 | `rho_c_e[e] = fks[e, 1] * (param.NE.rho * param.NE.heat_Q) + fks[e, 2] * (param.SP.rho * param.SP.heat_Q) + ...`（5 项求和，单行 ~250 字符） | 抽出 `layer_weighted_sum(fks[e,:], layers, field)` helper；同模式在 L37/L40 重复 |
| L80 | `iface = if (m1 == :PE && m2 == :PCC) \|\| (m1 == :PCC && m2 == :PE); :PE_PCC; elseif (m1 == :NE && m2 == :NCC) \|\| (m1 == :NCC && m2 == :NE); :NE_NCC; else; nothing; end`（4 个 `&&` + 2 个 `\|\|`，跨 L80-L86） | 与 `czm.jl` L80 同模式：抽出 `classify_interface(m1, m2)` 表驱动 helper，两文件共用 |
| L121-L140 | `if δ_eff > δ_max_hist; if δ_eff <= δ_0_eff; ...; elseif δ_eff >= δ_c_eff; ...; else; ...; end; ...; end`（嵌套 2 层 + 多个 `elseif`，跨 20 行） | 抽出 `compute_D_eq(δ_eff, δ_0_eff, δ_c_eff) -> Float64` 纯函数；该逻辑在 `bilinear_traction_state` (L121-L140) 与 `bilinear_tangent` (L225-L231) 重复 |
| L225-L231 | `if is_loading && δ_eff > δ_0_eff && δ_eff < δ_c_eff; ...; elseif is_loading && δ_eff >= δ_c_eff; ...; elseif is_loading && δ_eff <= δ_0_eff; ...; end`（多个 `&&` 链） | 与 L121-L140 同逻辑的独立实现，建议共用 helper |
| L258-L282 | mixed-mode 切线刚度的嵌套 `if δ_eff ≤ δ_0_eff \|\| !is_loading; if δ_n ≥ 0; ...; else; ...; end; elseif δ_eff ≥ δ_c_eff; ...; else; if δ_n ≥ 0 && δ_eff > 1e-15; ...; else; ...; end; end`（嵌套 3 层 + 多条件） | 抽出每段独立分支函数 `tangent_elastic / tangent_fractured / tangent_softening`，主流程分派 |
