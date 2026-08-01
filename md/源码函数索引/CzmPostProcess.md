# CzmPostProcess.jl

- **源文件**: `src/CzmPostProcess.jl`
- **行数**: 117 行
- **函数/struct 计数**: 5 个独立函数
- **职责**: CZM 后处理、统计与损伤管理——从 `CzmSolve.jl` 拆分出来，无状态、纯计算，独立于求解器主循环
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `get_damage_statistics(czm_mesh) -> NamedTuple` — L12-L32

统计 cohesive 网格的损伤分布。

- 返回字段：`max_D`、`mean_D`、`min_D`、`n_fractured`、`fraction_damaged`、`total_accumulated`
- 空 mesh（`n_cohesive == 0`）早退返回全零 NamedTuple（L15-L17）
- `n_damaged = count(d -> d > 0.01, D_vals)`（L21，阈值 0.01 硬编码）
- 注释 L6：`Statistics` 已在 `JuBat.jl` 模块级 `using`，此处无需重复

### `check_fracture_criterion(czm_mesh; threshold=0.99) -> (is_fractured, fracture_info)` — L37-L52

按损伤阈值判定整体是否断裂。

- 双判据：平均损伤 `mean_D >= threshold` 或断裂比例 `(n_fractured / n) > 0.5`（L40-L41）
- `fracture_info.criterion` 返回 `:average_damage` / `:fractured_count` / `:none`（L47）

### `reset_damage_states(czm_mesh) -> new_czm_mesh` — L59-L63

返回新的 CZM mesh，损伤状态重置为默认 `DamageState()`。

- 调 `clone_czm_mesh_with_damage(czm_mesh, new_damage_states)`（L61）

### `accumulate_cycle_damage(czm_mesh, cycle_damage_increment) -> new_czm_mesh` — L70-L94

循环加载下的累积损伤更新。

- 对每个未断裂单元：`accumulated_damage += cycle_damage_increment`（L82）
- 累积量 `>= 1.0` 时标记断裂并设 `D = 1.0`（L83-L86）
- 已断裂单元不更新（L81）

### `czm_output_to_variables(czm_mesh, result, variables) -> new_variables` — L99-L117

将 CZM 结果写入 `variables` 字典（CLAUDE.md §6.4）。

- 位移拆分：`u_x = result.displacement[1:2:end]`、`u_y = result.displacement[2:2:end]`（L101-L102）
- 写入键：`"czm displacement x/y"`、`"czm damage"`、`"czm traction normal/tangent"`、`"czm separation normal/tangent"`
- 统计字段：`"czm D_max"`、`"czm D_mean"`、`"czm n_fractured"`
- 调 `get_damage_statistics`（L111）

## 省略项

无。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info`。

### [PLACEHOLDER]

无。本文件无 TODO/FIXME/占位注释，无静默 try-catch，无 NaN/Inf 防御性代码。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L47 | `criterion = is_fractured_avg ? :average_damage : (is_fractured_count ? :fractured_count : :none)`（嵌套三元） | 简单判据链，可保留；或显式 `if/elseif/else` 写法提高可读性 |
