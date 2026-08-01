# CycleData.jl

- **源文件**: `src/CycleData.jl`
- **行数**: 623 行
- **函数/struct 计数**: 2 struct + 4 函数 = 6 条目
- **职责**: 循环数据导出——单时间步数据容器（`TimeStepData`）、整循环导出容器（`CycleExportData`）、带导出的单阶段求解（`solve_phase_with_export`）、单循环导出（`solve_cycle_with_export`）、CSV 序列化（`export_cycle_data_to_csv`）、CSV 反序列化（`load_cycle_data_from_csv`）
- **相关技术文档**: `md/09_分流求解器.md`、`md/10_参数传递与模块架构.md`

## 数据结构

### `struct TimeStepData` — L7-L18

单个时间步的数据快照（不可变）。

- 时间：`time::Float64`（物理秒）、`phase::PhaseType`
- 电：`V`、`I`
- 热：`T_nodes::Vector{Float64}`（节点温度 K）、`T_max`、`T_mean`
- SOC：`soc_n::Vector{Float64}`（负极各单元）、`soc_p::Vector{Float64}`（正极各单元）、`soc_mean`

### `struct CycleExportData` — L23-L30

整循环的导出数据。

- `cycle_idx::Int`
- `timesteps::Vector{TimeStepData}`
- 网格：`node_coords::Matrix{Float64}`（nT × 2）、`element_connectivity::Matrix{Int}`（ne × 4）
- `ne`、`nT`

## 函数清单

### `solve_phase_with_export(case, phase_type, t_max, I_current, V_limit, initial_state; czm_mesh, czm_params, dt_range, export_interval) -> (PhaseResult, Vector{TimeStepData})` — L32-L258

带数据导出的单阶段求解器（独立于 `CycleSolver.solve_phase`，不调用 `Solve`，自行实现时间推进循环）。

- **初始化**（L36-L123）：
  - 时间缩放：`dt_min`、`dt_max`、`t_end_nd`、`T_ref`、`phi_scale`、`I_scale`（L37-L43）
  - PHASE_REST 强制 I=0（L46-L48）
  - 多 SPMe 检测（L50）
  - 从 initial_state 取 y0 / T_nodes_init（L53-L54）；按多 SPMe 分支处理（L56-L82）
  - theta 参数（L85）
  - **预计算单元面积**（L116-L122）：累加 Gauss 权重 × detJ；`area_sum = sum(elem_area)` 用于物理加权均温
- **主循环**（L126-L231）：
  - 更新 `case.param.cell.T0 = mean(T_nodes_carry)`（L130-L132）
  - 电化学步：`CallModel(...) → Mt, Kt, Ft → y_c → y_new`（L135-L140）
  - **以步后状态重算 variables**（L143-L144）：避免导出混用步前/步后状态导致锯齿
  - 温度提取（L147-L156）；`T_max_phase = max(T_max_phase, T_max_current)`
  - 容量累积：`capacity += abs(I_current) * dt_dim / 3600.0`（L161）
  - **导出**（每 `export_interval` 步，L165-L199）：
    - SOC：处理 matrix / vector 两种形式（L169-L170）
    - **温度中点滤波**（L173）：`step_count == 1 ? T_nodes_carry : 0.5 .* (T_nodes_prev_export .+ T_nodes_carry)`，抑制 Crank-Nicolson 奇偶步伪振荡
    - 单元均温：`dot(T_elem_K, elem_area) / area_sum`（L180），回退 `mean(T_nodes_K)`（L182）
    - 构造 `TimeStepData` 并 push（L186-L198）
  - 截止检测（L204-L210）：`PHASE_CHARGE && V >= V_limit` 或 `PHASE_DISCHARGE && V <= V_limit`
  - 自适应 dt（L213-L220）
  - 状态更新 + 末尾对齐（L223-L230）
- **结果装配**（L233-L257）：`PhaseResult` 各字段填充；`D_max / D_mean / ΔD_max = 0.0`（本函数不处理 CZM，见 [PLACEHOLDER]）；`final_state` 包含 y / T_nodes / V / t_global

### `solve_cycle_with_export(case, cycle_opt; verbose, export_interval) -> CycleExportData` — L268-L410

单循环导出，调用 `solve_phase_with_export` 4 次（放电 → 静置1 → 充电 → 静置2）。

- verbose banner（L273-L279）
- 应用初始 SOC（L282-L288）
- initial_state 默认（L291-L296；`"V" => 3.7` 占位）
- 4 个阶段顺序调用 + verbose 结果打印（L300-L391）
- 构造 `CycleExportData`（L394-L401）

### `export_cycle_data_to_csv(export_data, output_dir; prefix) -> (6 file paths)` — L425-L504

将 `CycleExportData` 序列化为 6 个 CSV 文件。

- `mkpath(output_dir)` 确保目录存在（L427）
- 文件：
  1. `{prefix}_timesteps.csv`（L434-L443）：step / time_s / phase / V / I / T_max / T_mean / soc_mean
  2. `{prefix}_T_nodes.csv`（L447-L456）：每行一个时间步，每列一个节点
  3. `{prefix}_soc_n.csv`（L459-L468）
  4. `{prefix}_soc_p.csv`（L471-L480）
  5. `{prefix}_mesh_nodes.csv`（L483-L490）
  6. `{prefix}_mesh_elements.csv`（L493-L501）
- 每个文件保存后 `println("  ✓ 保存: $file")` 反馈（L444 等）

### `load_cycle_data_from_csv(input_dir; prefix) -> Dict` — L514-L622

从 CSV 文件加载为 `Dict{String, Any}`（与 `CycleExportData` 类型不同，见 [PLACEHOLDER]）。

- 5 个文件分别加载，缺失文件跳过（`isfile` 检查）
- timesteps：解析为 times / phases / voltages / currents / T_max / T_mean / soc_mean（L520-L550）
- T_nodes：解析为 `n_steps × nT` 矩阵（L554-L569）
- soc_n / soc_p：解析为 `n_steps × ne` 矩阵（L572-L590）
- node_coords：`nT × 2`（L593-L604）
- element_connectivity：`ne × 4`（L607-L619）

## 省略项

无。所有 struct 与 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L273 | `println("="^60); println("电化学-热耦合单循环数据导出"); ...; @printf("  网格: %d 单元, %d 节点\n", ne, nT)`（跨 L273-L279，verbose 守卫下） | 进度提示：单循环导出 banner |
| L285 | `@printf("  初始SOC: %.1f%%\n", ...)`（跨 L285-L287，verbose 守卫下） | 进度提示：初始 SOC |
| L301 | `println("\n[放电阶段]")`（verbose 守卫下） | 进度提示：放电阶段开始 |
| L317 | `@printf("  完成: %.1fs, %d 个数据点, %.3fV→%.3fV\n", ...)`（跨 L317-L320，verbose 守卫下） | 进度提示：放电阶段结果 |
| L325 | `println("\n[静置1阶段]")`（跨 L325-L327，verbose + `t_rest1 > 0` 守卫下） | 进度提示：静置1 阶段开始 |
| L341 | `@printf("  完成: %.1fs, %d 个数据点\n", ...)`（跨 L341-L343） | 进度提示：静置1 结果 |
| L348 | `println("\n[充电阶段]")`（跨 L348-L350） | 进度提示：充电阶段开始 |
| L364 | `@printf("  完成: %.1fs, %d 个数据点, %.3fV→%.3fV\n", ...)`（跨 L364-L367） | 进度提示：充电阶段结果 |
| L372 | `println("\n[静置2阶段]")`（跨 L372-L374） | 进度提示：静置2 阶段开始 |
| L387 | `@printf("  完成: %.1fs, %d 个数据点\n", ...)`（跨 L387-L389） | 进度提示：静置2 结果 |
| L404 | `println("\n" * "="^60); @printf("导出完成: 共 %d 个时间步数据点\n", ...)`（跨 L404-L407） | 进度提示：导出完成总结 |
| L444 | `println("  ✓ 保存: $timesteps_file")`（每个 CSV 文件保存后） | 进度提示：文件保存反馈，共 6 处（L444 / L456 / L468 / L480 / L490 / L501） |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L101 | `# 直接使用一致的初始状态，避免导出首步出现非物理温度尖峰。` 注释（修复性说明） | 历史占位说明：注释提到曾出现"非物理温度尖峰"，但未说明根因；后续修改 y0 初始化逻辑时需保持此修复 |
| L245 | `result.D_max = 0.0; result.D_mean = 0.0; result.ΔD_max = 0.0`（跨 L245-L247，PhaseResult 损伤字段全部填 0） | 占位：本函数不集成 CZM，CZM 用户若误用此入口会得到虚假的损伤=0；与 `CycleSolver.solve_phase` 行为不一致 |
| L294 | `"V" => 3.7`（initial_state 默认电压占位） | 占位初值：与 CycleSolver.jl L256 同问题；3.7V 是典型 OCV 但非参数化 |
| L514 | `function load_cycle_data_from_csv(...) -> Dict`（返回 Dict 而非 `CycleExportData`） | 类型不对称：导出是 `CycleExportData`，加载是 Dict；调用方无法对称使用；`TimeStepData` 等结构体在加载路径丢失 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L56 | `if y0 === nothing; if multi_spme; y0 = ModelInitialisation_MultiSPMe(case); else; y0 = ModelInitialisation(case); end; else; y0 = vec(y0); if multi_spme; # layout is now stored as case.layout::Union{Nothing,MultiSPMeLayout}; end; end; if T_nodes_init !== nothing; T_nodes_carry = copy(T_nodes_init); if multi_spme; thermal_range = case.layout.thermal_range; y0[thermal_range] .= T_nodes_carry; end; elseif multi_spme; thermal_range = case.layout.thermal_range; T_nodes_carry = y0[thermal_range]; else; nT_mesh = case.mesh["thermal2D"].nlen; T_nodes_carry = y0[(end - nT_mesh + 1):end]; end`（嵌套 3 层 if + 多 SPMe 与 T_nodes 的 3 路组合，跨 L56-L82） | 抽出 `init_y0_with_thermal(case, initial_state, multi_spme) -> (y0, T_nodes_carry)` helper |
| L85 | `theta = case.opt.solveType == "Crank-Nicolson" ? 0.5 : (case.opt.solveType == "forward" ? 0.0 : 1.0)`（三元嵌套 + 默认 backward） | 接近阈值但未越界；可改为显式 if-elseif + error |
| L165 | 导出逻辑：`if step_count % export_interval == 0; soc_n_raw = ...; soc_n = isa(...) ? ... : ...; soc_p_raw = ...; soc_p = isa(...) ? ... : ...; T_nodes_out = step_count == 1 ? ... : ...; T_nodes_K = ...; T_max_K = ...; if ne > 0 && area_sum > 0.0; T_elem_K = [...]; T_mean_K = ...; else; T_mean_K = ...; end; ts_data = TimeStepData(...); push!(...); end`（嵌套 3 层 if + 三元运算 + 列表推导，跨 L165-L199） | 抽出 `build_timestep_data(...)` 函数；当前 SOC matrix/vector 判别 + 温度中点滤波 + 均温多路径全部内联 |
| L204 | `if phase_type == PHASE_CHARGE && V_current >= V_limit; terminated_by = :voltage; break; elseif phase_type == PHASE_DISCHARGE && V_current <= V_limit; terminated_by = :voltage; break; end`（跨 L204-L210，双 `&&` 链 × 2 分支） | 接近阈值但未越界（每分支 2 个条件） |
| L228 | `if t + dt > t_end_nd; dt = t_end_nd - t; end`（跨 L228-L230，边界对齐） | 简单条件，未越界 |
| L300 | 4 个阶段调用 `solve_phase_with_export` 的代码块结构高度相似（每个阶段约 20-25 行），仅在 phase_type / 电流符号 / verbose 消息上差异（跨 L300-L391） | 抽出 `run_exported_phase(case, phase_type, t, I, V, current_state; verbose, dt_range, export_interval) -> (result, data, new_state)` helper；当前 4 份近似代码 |
