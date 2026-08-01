# CsvExport.jl

- **源文件**: `src/CsvExport.jl`
- **行数**: 637 行
- **函数/struct 计数**: 16 个（1 struct + 2 便捷构造器 + 13 函数）
- **职责**: 循环仿真后处理 CSV 导出——按可配置采样策略（`:full`/`:phase_ends`/`:custom`）写出 cycle_summary、element_currents、node_temperature、cohesive_damage、node_displacement、cohesive_driving_force、czm_solver_diagnostics 七类 CSV 文件；安全访问 `solve_result` 字典中的 `variables["..."]` 字符串键并做越界保护
- **已知技术债**: 大量 `variables[...]` 字符串键访问（详见 `MEMORY.md`——全项目 17 文件 / 478 处），键名拼写错误在运行时才暴露；新增字段需同步多处
- **相关技术文档**: `md/10_参数传递与模块架构.md`（`variables` 字典结构）、`md/06_内聚力模型_CZM.md`（CZM 字段语义）、`CLAUDE.md` §6 输出变量表

## 数据结构

### `struct CsvExportOptions` — L20-L25

CSV 导出采样配置（不可变）。4 个字段：

- `mode::Symbol`：`:full`（所有步）/`:phase_ends`（每相位首末）/`:custom`（每 N 步）（L21）
- `save_every::Int`：`:custom` 模式下的步长间隔（L22）
- `full_output_cycles::Vector{Int}`：强制完整输出的循环号列表（L23）
- `skip_files::Vector{String}`：完全跳过的文件名（如 `["node_displacement.csv"]`）（L24）

### 便捷构造器 — L27-L28

- `CsvExportOptions()`：默认 `:full` 模式，`full_output_cycles=[1]`（L27）
- `CsvExportOptions(mode::Symbol)`：指定模式，其余默认（L28）

## 函数清单

### `_validate_csv_opt(csv_opt::CsvExportOptions)` — L30-L37

入口参数校验：`mode` 必须 ∈ `(:full, :phase_ends, :custom)`，`:custom` 模式下 `save_every ≥ 1`，否则 `error`。

### `_should_output_step(csv_opt, cycle, ti, n_steps)` — L40-L46

判断时间步 `ti` 是否应写入。短路链：`full_output_cycles` 命中 → `:full` → `:phase_ends` 首/末 → `:custom` 首/末/间隔。

### `_compute_last_snap_indices(czm_snapshots)` — L50-L65

预计算每个 (cycle, phase) 分组最后一条快照的索引集合。空数组返回空 `Set`。线性扫描 + key 比较。

### `_should_output_snapshot(csv_opt, cycle, snap_idx, last_indices)` — L68-L75

判断 CZM 快照是否应写入。`:phase_ends` 模式查 `last_indices` 集合。

### `export_cycling_csv(result, case, czm_mesh; output_dir, overwrite, csv_opt)` — L88-L208

主入口。逐文件 `try/catch` 包裹，单文件失败 `@warn` + 推入 `files_skipped`，不阻断其余文件。CZM 类文件需 `!isempty(result.czm_snapshots)`；`cohesive_driving_force.csv` 额外需 `case.geometry !== nothing`。末尾 `println` 总结 Written/Skipped/Directory。返回 `files_written::Vector{String}`。

### `_write_cycle_summary(result, output_dir, overwrite)` — L214-L241

写 `cycle_summary.csv`（11 列：cycle, phase, capacity_ah, soh, D_max, D_mean, n_fractured, T_max_K, T_mean_end_K, V_start, V_end）。`soh` 越界时填 `NaN`；`D_max/D_mean/T_mean_end/V_start/V_end` 的 NaN 改写为 0.0（L226-L235）。

### `_write_element_currents(result, case, output_dir, overwrite, csv_opt)` — L247-L314

写 `element_currents.csv`（14 列含电流、SOC、过电位、热源）。通过 `_safe_get` 二维数组越界保护读取 `solve_result` 字典中 10 个不同的字符串键（`"thermal2D element current"` 等，L279-L287）。物理面积由 `elem_areas[e] * scale.L^2` 从归一化还原（L305，详见 `CLAUDE.md` §9.1）。

### `_write_node_temperature(result, case, output_dir, overwrite, csv_opt)` — L320-L367

写 `node_temperature.csv`（7 列含 node_id, x, y, T_K）。节点坐标由 `mesh_th.node * scale.L` 还原（L331-L332）。仅当 `case.opt.thermal_enabled` 时调用。

### `_write_cohesive_damage(result, case, czm_mesh, output_dir, overwrite, csv_opt)` — L373-L421

写 `cohesive_damage.csv`（12 列含 coh_id, length, D, δ_n/δ_t, T_n/T_t, fractured, θ_deg）。预计算每个 cohesive 单元的物理长度与方位角（L386-L394）。分离/牵引量分别通过 `scale.δ_czm` / `scale.σ_czm` 还原（L411-L414）。`fractured = D >= 0.95`（L415）。

### `_write_node_displacement(result, case, czm_mesh, output_dir, overwrite, csv_opt)` — L427-L462

写 `node_displacement.csv`（8 列含 node_id, x, y, ux, uy）。位移由 `snap.displacement * scale.L` 还原（L455-L456）。

### `_write_driving_force(result, case, czm_mesh, output_dir, overwrite, csv_opt)` — L468-L567

写 `cohesive_driving_force.csv`（10 列含 dT, Δsoc, ε_thermal/diffusion/total）。**当前由于 `alpha_eff/beta_n/beta_p` 硬编码为 0.0（L484-L486），实际总是提前 return（L488-L491），不写出任何内容**。集成 `compute_czm_effective_params(case)` 后才能产出数据。

### `_write_czm_diagnostics(result, output_dir, overwrite)` — L573-L588

写 `czm_solver_diagnostics.csv`（7 列含 converged, iterations, residual_norm, method）。直接遍历 `result.czm_snapshots`，无采样过滤。

### `_safe_get(arr, row, col)` — L595-L601

二维数组越界保护，越界返回 `NaN`。两个方法：`AbstractArray{<:Real,2}` 走边界检查（L595-L600），其他类型 fallback 到 `NaN`（L601）。

### `_find_solve_result(result, cycle, phase)` — L604-L619

线性扫描 `cycle_results` 匹配 `cycle_idx` 与 phase 名（`discharge`/`rest1`/`charge`/`rest2`），未命中返回 `nothing`。

### `_compute_element_areas(mesh::Mesh)` — L622-L636

对 Q4 单元用 Shoelace 公式计算归一化面积（无量纲）。返回 `Vector{Float64}`，长度为单元数。

## 省略项

无。所有 struct 与 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L102 | `@warn "Failed to write cycle_summary.csv" exception=e` | 单文件写入失败的错误回执；结构化 `@warn` 但属运行时诊断（非参数验证） |
| L114 | `@warn "Failed to write element_currents.csv" exception=e` | 同上；element_currents 写入异常 |
| L131 | `@warn "Failed to write node_temperature.csv" exception=e` | 同上；node_temperature 写入异常 |
| L149 | `@warn "Failed to write cohesive_damage.csv" exception=e` | 同上；cohesive_damage 写入异常 |
| L160 | `@warn "Failed to write node_displacement.csv" exception=e` | 同上；node_displacement 写入异常 |
| L178 | `@warn "Failed to write cohesive_driving_force.csv" exception=e` | 同上；cohesive_driving_force 写入异常 |
| L194 | `@warn "Failed to write czm_solver_diagnostics.csv" exception=e` | 同上；czm_solver_diagnostics 写入异常 |
| L201 | `println("CSV export complete:")` | 用户可见的状态汇报，非 debug 打印 |
| L202 | `println("  Written: $(join(files_written, ", "))")` | 同上 |
| L204 | `println("  Skipped: $(join(files_skipped, ", "))")` | 同上 |
| L206 | `println("  Directory: $output_dir")` | 同上 |
| L217 | `println("  Skipping $filepath (already exists)")` | 文件已存在时跳过提示（重复模式） |
| L240 | `println("  Written: $filepath")` | 单文件完成提示（重复模式） |
| L252 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L313 | `println("  Written: $filepath")` | 同 L240 模式 |
| L325 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L366 | `println("  Written: $filepath")` | 同 L240 模式 |
| L378 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L420 | `println("  Written: $filepath")` | 同 L240 模式 |
| L432 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L461 | `println("  Written: $filepath")` | 同 L240 模式 |
| L474 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L489 | `println("  Skipping $filepath (no thermal/diffusion strain parameters)")` | 提示 driving_force 因参数为 0 提前退出 |
| L495 | `println("  Skipping $filepath (no geometry data)")` | 提示 driving_force 因缺几何数据退出 |
| L566 | `println("  Written: $filepath")` | 同 L240 模式 |
| L576 | `println("  Skipping $filepath (already exists)")` | 同 L217 模式 |
| L587 | `println("  Written: $filepath")` | 同 L240 模式 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L226 | `soh_val = i <= length(result.soh) ? result.soh[i] : NaN` | 越界 fallback 直接写 `NaN` 到 CSV——下游脚本若不做 NaN 处理会污染统计；建议改用空字段或显式 `"NA"` |
| L484 | `alpha_eff = 0.0` | 硬编码占位——`_write_driving_force` 在此前提下永远提前 return（L488-L491），CSV 文件不会生成。注释 L482-L483 标注 "default to 0.0 until compute_czm_effective_params(case) is integrated"，是显式 placeholder |
| L485 | `beta_n = 0.0` | 同 L484，硬编码占位 |
| L486 | `beta_p = 0.0` | 同 L484，硬编码占位 |
| L599 | `_safe_get` 越界返回 `NaN` | 越界 fallback——`variables[...]` 字符串键名拼写错误或 solve_result 中字段缺失时静默返回 NaN，写入 CSV 后下游才能发现；建议至少 `@warn` 一次 |
| L601 | `_safe_get(arr, row, col) = NaN` | 类型 fallback——非 `AbstractArray{<:Real,2}` 一律返回 NaN；若上游类型误传（如 `Matrix{Any}`），所有字段静默为 NaN |
| L415 | `frac = D >= 0.95` | magic number 0.95——`fractured` 阈值硬编码，未与 `CLAUDE.md` §5.3 `czm_area_loss_threshold=0.83` 对齐，可能造成 "D=0.9 标记未断裂但面积已开始缩减" 的不一致 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L107 | `has_data = !isempty(result.cycle_results) && any(cr -> cr.discharge.solve_result !== nothing, result.cycle_results)` | 2 个 `&&` + 嵌套 lambda——可读性可接受，无需重构 |
| L143 | 嵌套 3 层 `if` + 多个 `in csv_opt.skip_files` 判定（"4. cohesive_damage + 5. node_displacement" 块跨 L143-L169） | 抽出 `_try_write(filename, skip_files, write_fn)` helper，统一处理 skip 判定 + try/catch + push! 列表，可消除 ~7 处重复（参见 L99-L198 全部 7 个文件块结构同构） |
| L172 | `!isempty(result.czm_snapshots) && case.geometry !== nothing` | 2 个条件链，可读性可接受 |
| L296 | `I_e = _safe_get(I_e_all, e, ti); sn = _safe_get(soc_n_all, e, ti); sp = _safe_get(soc_p_all, e, ti); ...` 连续 9 个 `_safe_get` 调用（L296-L304） | 抽出 `extract_row_at(sr, keys, e, ti) = ...` helper 返回 NamedTuple，避免列名拼写错误与重复模板 |
| L408 | `for i in 1:min(n_coh, length(snap.damage)); dn = i <= length(snap.separation_n) ? ... : 0.0; dt_val = i <= length(snap.separation_t) ? ... : 0.0; tn = i <= length(snap.traction_n) ? ... : 0.0; tt = i <= length(snap.traction_t) ? ... : 0.0` | 4 个独立长度守卫——可抽 `_safe_index(v, i, scale) = i <= length(v) ? v[i] * scale : 0.0` helper |
| L548 | `if e_top <= size(T_elem, 1) && e_bot <= size(T_elem, 1)` 后跟 L551/L554 两个同构判定（driving_force 块跨 L548-L556） | 抽 `_mean_pair_safe(arr, e_top, e_bot, ti) = ...` helper |
