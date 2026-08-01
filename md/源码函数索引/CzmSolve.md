# CzmSolve.jl

- **源文件**: `src/CzmSolve.jl`
- **行数**: 675 行
- **函数/struct 计数**: 1 个 struct + 12 个独立函数
- **职责**: CZM Newton-Raphson 求解器——`CZMResult` 结果对象、损伤状态克隆、按界面分批更新损伤、回溯线搜索、basic / load-substep / arc-length 三套求解策略、统一入口 `solve_czm_step`
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`md/14_粘性正则化.md`、`md/09_分流求解器.md`

## 数据结构

### `mutable struct CZMResult` — L1-L15

CZM 单步求解结果容器。

- 字段：`displacement`、`damage`、`traction_n/t`、`separation_n/t`、`converged::Bool`、`iterations::Int64`、`residual_norm::Float64`
- 默认构造 `CZMResult(ndof, n_coh)`：`residual_norm` 初始化为 `Inf`（L14），`converged=false`

## 函数清单

### `clone_damage_states(damage_states) -> Vector{DamageState}` — L17-L29

深拷贝损伤状态向量（用于 Newton 回滚、循环累积）。

### `clone_czm_mesh_with_damage(czm_mesh, damage_states) -> CohesiveMesh` — L31-L44

构建新的 `CohesiveMesh`（共享 bulk_mesh/node/etc. 引用，仅替换 `damage_states`）。

### `update_damage_per_interface(czm_mesh, damage_states, separations, param_cache; visc_beta) -> new_states` — L52-L71

按 `interface_type` 分组调 `update_damage`，保持原始顺序。

- 按 `keys(param_cache.by_interface)` 分批；空批跳过
- 跨文件依赖：`update_damage`（在 `Materialmatrix.jl`）

### `extract_bc_dofs(czm_mesh, param; cache, fix_inner) -> (bc_dofs, bc_vals)` — L79-L97

提取 Dirichlet BC 自由度列表与对应值。

- 优先返回 `cache.bc_dofs` / `cache.bc_vals`（L80-L82）
- 否则调 `identify_bc_nodes_czm` 现场计算

### `backtrack_line_search!(u, Δu, czm_mesh, param_cache, damage_states, F_ext, F_thermo_chem, R_norm_current, bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws; max_halvings, visc_beta) -> (u_new, R_new_norm, accepted, α_used)` — L107-L128

回溯线搜索（零化式 BC 残差），仅用于 `solve_czm_basic_step`。

- 迭代 `α ← α/2`，最多 `max_halvings=8` 次（L109）
- 接受判据：`R_trial_norm < R_norm_current && !isnan`（L121）
- accepted 时 `u_new` 已含 BC 赋值；未 accepted 时返回原始 `u`（未修改）
- 注：函数名带 `!` 但实际不原位修改 `u`，语义上"原位"指返回值会替代调用方 `u`

### `apply_czm_dirichlet!(u, bc_dofs, bc_vals) -> u` — L130-L135

原位强制 Dirichlet 值（`u[dof] = val`）。

### `zero_czm_bc_entries!(v, bc_dofs) -> v` — L137-L142

原位将 `bc_dofs` 位置置零（用于增量载荷向量）。

### `fill_czm_result!(result, u, damage_states, separations, tractions) -> result` — L144-L154

原位填充 `CZMResult` 各字段（位移、损伤、分离、牵引）。

### `build_arc_length_augmented_matrix(K_bc, load_vector, delta_u, delta_lambda, arc_length_alpha) -> A` — L156-L166

构造弧长法增广矩阵（spec §Crisfield cylindrical arc-length）。

- `A[1:ndof, 1:ndof] = K_bc`；`A[i, ndof+1] = -load_vector[i]`；`A[ndof+1, i] = 2·delta_u[i]`
- 末对角 `A[ndof+1, ndof+1] = 2·α²·delta_lambda`

### `solve_czm_basic_step(czm_mesh, F_ext, param_cache, param, u_prev; α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem, max_iter, tol, cache, visc_beta) -> (result, new_czm_mesh)` — L168-L263

基础 Newton-Raphson 单步求解（无载荷子步）。

- 每轮（L194-L235）：装配 `K_total, f_int_total` → 残差 `R = F_ext + F_thermo_chem - f_int_total` → BC 零化 → `K_bc \ R_bc` → 线搜索
- 收敛判据（L211）：`R_norm < tol` 或 `rel_norm < tol`（`rel_norm = R_norm / R_norm_0`）
- 失败回滚：未收敛时恢复 `u_start` / `damage_start`（L237-L240）
- 收敛后才更新损伤（L214，与 `newton_raphson_czm` 一致：冻结损伤求解位移，收敛后更新）
- 跨文件依赖：`assemble_coupled_system`、`assemble_thermal_chemical_load`、`apply_bc_czm`、`backtrack_line_search!`、`update_damage_per_interface`、`clone_damage_states`、`clone_czm_mesh_with_damage`

### `solve_czm_arc_length_step(czm_mesh, F_ext, param_cache, param, u_prev; α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem, max_iter, tol, n_load_steps, arc_length_alpha, cache, visc_beta) -> (result, new_czm_mesh)` — L265-L466

Crisfield 圆柱弧长法单步求解。

- 外层 `while load_progress < 1 - 1e-12`（L297）：逐步增加 `load_progress`
- 切线预测（L319-L333）：`tangent = K_bc \ F_load_bc` → `delta_u_pred = tangent · delta_lambda_pred`
- 内层 Newton 修正（L344-L421）：解两个线性系统 `K_bc \ R_bc` 与 `K_bc \ F_load_bc`，二次方程选根（L391-L410）
- 失败时 `step_size /= 2`，小于 `step_size_min` 触发 `@warn`（L430）
- 跨文件依赖：同 `solve_czm_basic_step` + `build_arc_length_augmented_matrix`

### `newton_raphson_czm(czm_mesh, F_ext, param_cache, param; α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem, max_iter, tol, u0, n_load_steps, cache, visc_beta) -> (result, new_czm_mesh)` — L481-L644

Newton-Raphson + 自适应载荷子步（iter_method="load_substep"）。

- 增量载荷参考（L502）：`F_delta = F_target - f_int_ref`，`f_int_ref = f_int(u_prev)`（u_prev 近似在上一步平衡态）
- 外层 `while load_progress < 1 - 1e-12`（L516）；子步容差 `substep_tol = tol * 10.0`（L543）
- 内层 Newton（L527-L595）：每轮含 8 次线搜索（L561-L584）
- 收敛后批量更新损伤（L629-L631）
- 收敛判据（L623）：`load_progress >= 1 - 1e-12 && R_norm < final_tol`（`final_tol = tol * 100.0`）
- 失败时 `step_size *= 0.5`，小于 `step_size_min` 触发 `@warn`（L602-L603）

### `solve_czm_step(czm_mesh, F_ext, param_cache, param, u_prev; iter_method, ...) -> (result, new_czm_mesh)` — L651-L675

CZM 单步求解统一入口，按 `iter_method` 分派。

- `"load_substep"`（默认）→ `newton_raphson_czm`
- `"basic"` → `solve_czm_basic_step`
- `"arc_length"` / `"arclength"` / `"arc-length"` → `solve_czm_arc_length_step`
- 未知方法 `error(...)`（L673）

## 省略项

无。所有函数与 struct 均独立列出。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试用途的 `@info`；`@debug`（L434, L607）是 Julia 结构化日志宏，按 S3 规则不计入；`@warn`（L430, L603）为求解器停滞告警，属于运行时状态而非调试输出。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L121 | `if !isnan(R_trial_norm) && R_trial_norm < R_norm_current`（NaN 防御性检查，回溯线搜索接受判据） | 静默丢弃 NaN 试验步，可能掩盖求解器数值不稳定；非占位而是数值防御 |
| L226 | `if any(isnan, Δu) \|\| any(isinf, Δu)`（basic step 失败检测，静默 break） | 静默退出循环，外部通过 `converged=false` 判定失败；可考虑加 `@warn` 报告原因 |
| L325 | `if tangent === nothing \|\| any(isnan, tangent) \|\| any(isinf, tangent)`（arc-length 切线失败，静默 break） | 同 L226：静默失败，外部 `converged=load_progress >= ...` 判定；建议加诊断信息 |
| L377 | `if delta_u_R === nothing \|\| any(isnan, delta_u_R) \|\| any(isinf, delta_u_R)`（arc-length 修正步失败） | 同 L325 |
| L386 | `if delta_u_F === nothing \|\| any(isnan, delta_u_F) \|\| any(isinf, delta_u_F)`（arc-length 修正步失败） | 同 L325 |
| L514 | `converged_substep = false`（子步初始化，配合 `last_R_norm = Inf`） | 初始值合理，非占位 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L107 | `backtrack_line_search!` 函数签名（14 个位置参数 + 2 个关键字参数，单行 ~150 字符） | 拆为 `(czm_mesh, param_cache, bc, loading, u_state)` 五个语义对象，或将 `bc_dofs/bc_vals/geom_cache/ws/K_bulk_cached` 包进 `cache::CZMAssemblyCache` 直接传 |
| L211 | `if R_norm < tol \|\| rel_norm < tol`（双判据，2 个 `\|\|`） | 简单判据，可保留 |
| L325 | `if tangent === nothing \|\| any(isnan, tangent) \|\| any(isinf, tangent)`（3 个 `\|\|`，配合 `any` 谓词） | 抽出 `is_finite_solution(x) = x !== nothing && all(isfinite, x)` 谓词函数，统一用于 L325/L377/L386 |
| L377 | `if delta_u_R === nothing \|\| any(isnan, delta_u_R) \|\| any(isinf, delta_u_R)` | 同 L325，可共用谓词 |
| L386 | `if delta_u_F === nothing \|\| any(isnan, delta_u_F) \|\| any(isinf, delta_u_F)` | 同 L325，可共用谓词 |
| L555 | `if any(isnan, Δu) \|\| any(isinf, Δu)`（load_substep 内层 Newton 失败检测） | 与 L325 同类，可共用 `is_finite_solution` 谓词 |
| L651 | `solve_czm_step` 函数签名（4 位置 + 15 关键字，单行 ~250 字符） | 拆为 `CzmSolveOptions` struct 收纳 α/β/tol/n_load_steps/iter_method 等 |
