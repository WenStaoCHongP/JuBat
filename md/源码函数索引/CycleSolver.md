# CycleSolver.jl

- **源文件**: `src/CycleSolver.jl`
- **行数**: 547 行
- **函数/struct 计数**: 3 struct + 2 默认构造函数 + 4 函数 = 9 条目
- **职责**: 循环仿真顶层——阶段结果与循环结果容器（`PhaseResult`/`CycleResult`/`CyclingResult`）、单阶段求解（`solve_phase`）、循环主循环（`solve_cycling`）、SOC 初始化辅助（`compute_cs0_from_soc`/`apply_initial_soc!`）
- **相关技术文档**: `md/09_分流求解器.md`、`md/10_参数传递与模块架构.md`

## 数据结构

### `mutable struct PhaseResult` — L6-L32

单阶段（充电/静置/放电）的求解结果。

- 电化学：`V_start`、`V_end`、`capacity`、`terminated_by::Symbol`（:time / :voltage / :none）
- 热：`T_max`、`T_mean_end`
- 损伤：`D_max`、`D_mean`、`ΔD_max`
- 链接：`final_state::Union{Dict, Nothing}`（传递给下一阶段）、`solve_result::Any`（CSV 导出用，仅 `save_detailed=true`）

### `mutable struct CycleResult` — L49-L70

单循环结果，包含四个 `PhaseResult`（charge / rest1 / discharge / rest2）+ 循环汇总（容量、库伦效率、损伤、温度）。

### `mutable struct CyclingResult` — L75-L98

全部循环汇总。每循环指标以 Vector 存储（`capacity_charge` / `capacity_discharge` / `coulombic_efficiency` / `D_max` / `D_mean` / `n_fractured` / `T_max` / `soh`）。还包含 `cycle_results::Vector{CycleResult}`（详细，可选）、`czm_snapshots::Vector{CZMSnapshot}`、`final_czm_mesh::Any`、`initial_capacity`、`soh_terminated`。

### `PhaseResult()` — L35-L44

默认构造函数，全 0 / `PHASE_REST` / `:none` / nothing 填充。

### `CyclingResult(n_cycles::Int)` — L100-L111

默认构造函数，向量初始化为空，`soh_terminated=false`、`initial_capacity=0.0`。

## 函数清单

### `solve_phase(case, phase_type, t_max, I_current, V_limit, initial_state; czm_mesh, czm_params, dt_range, czm_snapshots, czm_cycle) -> PhaseResult` — L118-L197

单阶段求解器。

- `PHASE_REST` 时强制 `I_current = 0.0`（L123-L125）
- 记录阶段开始损伤 `D_max_init` / `D_mean_init`（L128-L129）
- 临时覆盖 `case.opt.Current` / `time` / `dt`（L132-L134）
- **临时修改电压限制**（L136-L153）：保存 `old_v_l` / `old_v_h`；按阶段类型设定：放电设 `v_l=V_limit, v_h=1e6`；充电设 `v_l=-1e6, v_h=V_limit`；静置全放开；finally 块恢复（L191-L194）
- 调 `Solve(case; initial_state, return_final_state=true, czm_snapshots, czm_cycle, czm_phase=string(phase_type))`（L155-L157）
- 提取 `duration`、`final_state`，`final_state["t_global"]` 设为 `t_start + duration`（L159-L165）
- **CZM 损伤不再此处理**（L167-L170 注释）：旧版在此调用 `update_czm_damage!(..., u_czm_prev=nothing)` 会导致位移场从零重解，每阶段额外浪费约 12 s；现由 Solve.jl 主循环每步更新
- 调 `_postprocess_phase_result(...)`（L172-L176）计算 phase_data，写回 `result`（L178-L190）

### `solve_cycling(case, cycle_opt, czm_mesh=nothing; verbose, save_detailed) -> CyclingResult` — L217-L460

循环主循环。循环顺序：放电 → 静置1 → 充电 → 静置2。

- 创建 `CyclingResult(n_cycles)`（L220）、共享 `czm_snaps`（L223）
- **初始 SOC 应用**（L226-L234）：`soc_init ∈ [0,1]` 时调 `apply_initial_soc!`，verbose 时打印 cs0_NE / cs0_PE
- **verbose banner**（L237-L249）：打印循环次数 / 初始 SOC / 各阶段参数
- **initial_state 初始化**（L252-L257）：`"V" => 3.7`、`"t_global" => 0.0`、y / T_nodes 为 nothing
- **SOH 监控**（L262-L265）：阈值取 `case.opt.czm_soh_threshold`；`initial_capacity=0.0`、`soh_terminated=false`
- **循环主体**（L268-L447）：每个 cycle
  1. 创建 `CycleResult`（L275-L281）
  2. `reset_T_each_cycle && cycle > 1` 时清空 `current_state["T_nodes"]`（L284-L289）
  3. 放电：`solve_phase(..., PHASE_DISCHARGE, ...)`（L298-L310）
  4. 静置1：`t_rest1 > 0` 时执行（`solve_phase(..., PHASE_REST, I=0, ...)`）；否则跳过，`final_state = current_state`、`V_start = V_end = get(current_state, "V", NaN)`（L321-L356）
  5. 充电前可选 `reset_T_before_charge`（L360-L366）
  6. 充电：`solve_phase(..., PHASE_CHARGE, I=-I_charge, ...)`（L371-L390）
  7. 静置2：同静置1（L393-L426）
  8. 后处理：调 `_postprocess_cycle_result!` / `_append_cycle_result!` / `_update_soh_and_capacity!`（L429-L433）
  9. verbose 时 `_print_cycle_summary`（L435-L437）
  10. 终止检查 `_check_cycle_termination` → `should_stop / soh_hit`（L439-L446）
- 附加 `czm_snapshots`、`final_czm_mesh`，verbose 时 `_print_cycling_summary`（L449-L457）

### `compute_cs0_from_soc(param_dim, soc::Float64) -> (cs0_NE, cs0_PE)` — L495-L512

按目标 SOC 计算正负极初始锂浓度。

- 校验 `soc ∈ [0,1]`（L497-L499），否则 `error`
- 负极：`theta_n = NE.theta_0 + soc * (NE.theta_100 - NE.theta_0)`（L503）；`cs0_NE = theta_n * NE.cs_max`（L504）
- 正极：`theta_p = PE.theta_0 - soc * (PE.theta_0 - PE.theta_100)`（L508）；`cs0_PE = theta_p * PE.cs_max`（L509）

### `apply_initial_soc!(case, param_dim, soc::Float64) -> (cs0_NE, cs0_PE)` — L533-L546

将 SOC 设置写入 case。

- 调 `compute_cs0_from_soc`（L535）
- 更新维度参数 `param_dim.NE.cs0` / `param_dim.PE.cs0`（L538-L539）
- 重新归一化到 `case.param.NE.cs0 = param_dim.NE.cs0 / param_dim.NE.cs_max`（L542-L543）

## 省略项

以下 helper 函数定义在其他文件，本文件内仅引用，不单列：
- `_postprocess_phase_result` / `_postprocess_cycle_result!` / `_append_cycle_result!` / `_update_soh_and_capacity!` / `_print_cycle_summary` / `_check_cycle_termination` / `_print_cycling_summary`：循环求解器辅助函数（应在 CsvExport 或 PostProcessing 相关文件定义）

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L229 | `@printf("  初始SOC: %.1f%%\n", ...); @printf("    → 负极cs0: ...")` 等 3 行（跨 L229-L233，在 `verbose` 守卫下） | 进度提示：初始 SOC 详细信息 |
| L237 | `println("="^60); println("开始充放电循环仿真"); ...; @printf("  循环次数: %d\n", ...)` 等 8 行（跨 L237-L249，在 verbose 守卫下） | 进度提示：循环仿真 banner |
| L269 | `println("\n" * "-"^40); @printf("循环 %d/%d\n", ...)` 等 3 行（跨 L269-L273，在 verbose 守卫下） | 进度提示：每个循环开始 |
| L287 | `println("  温度场已重置")`（在 verbose + `reset_T_each_cycle` 守卫下） | 进度提示：温度场重置 |
| L294 | `print("  [放电] ")` 等 4 处阶段开始标签（跨 L294-L296，verbose 守卫下） | 进度提示：各阶段开始 |
| L312 | `@printf("%.1fs, %.3fV→%.3fV, %.3fAh (%s)\n", ...)`（跨 L312-L317，放电后） | 进度提示：放电阶段结果 |
| L324 | `print("  [静置1] 状态继承+锂扩散 ")`（verbose 守卫下） | 进度提示：静置阶段开始 |
| L342 | `@printf("%.1fs, T_max=%.2fK\n", ...)`（静置1 后） | 进度提示：静置1 结果 |
| L354 | `println("  [静置1] 跳过 (t=0，无扩散)")`（verbose + `t_rest1 == 0` 守卫下） | 进度提示：静置跳过说明 |
| L364 | `println("    (温度场已重置)")`（verbose + `reset_T_before_charge` 守卫下） | 进度提示：充电前温度重置 |
| L385 | `@printf("%.1fs, %.3fV→%.3fV, %.3fAh (%s)\n", ...)`（跨 L385-L390，充电后） | 进度提示：充电阶段结果 |
| L412 | `@printf("%.1fs, T_max=%.2fK\n", ...)`（跨 L412-L413，静置2 后） | 进度提示：静置2 结果 |

注：以上全部 `println / @printf` 均在 `verbose` 守卫下，属于用户可见的进度输出而非 debug 日志，但根据规则统一列入 [DEBUG] 表（非结构化 `println` / `print`）。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L148 | `case.param.cell.v_l = -1.0e6` / `v_h = 1.0e6`（充电/静置时放开电压限制） | magic number 占位：用 ±1e6 表示"无限制"；若真实电压意外接近 1e6 V（不可能但语义不清），更好的方式是用 `Inf` 或显式 `nothing` 标记 |
| L167 | `# CZM 损伤已由 Solve 主循环每步更新（Solve.jl:270-285），此处不再冗余调用。` 注释说明旧版逻辑已废弃（跨 L167-L170） | 历史代码注释：注释提到旧版每阶段调用 `update_czm_damage!(..., u_czm_prev=nothing)` 会浪费 12s/阶段，但未删除注释；维护时需注意跨文件依赖 |
| L256 | `"V" => 3.7`（initial_state 默认电压占位） | 占位初值：3.7V 是典型开路电压但并非来自参数；第一次 CallModel 调用会被实际计算值覆盖 |
| L264 | `initial_capacity = 0.0`（初始容量占位，第一个循环后由 `_update_soh_and_capacity!` 设置） | 占位初值：SOH 计算依赖 `initial_capacity`，若第一个循环失败会导致除零或 NaN SOH |
| L348 | `cycle_result.rest1.V_start = get(current_state, "V", NaN); cycle_result.rest1.V_end = cycle_result.rest1.V_start`（跨 L348-L349，跳过静置时用 NaN 兜底） | 占位 NaN：若 current_state 缺 "V" 键则记 NaN；下游 CSV / 分析可能未处理 NaN |
| L419 | 同上，rest2 跳过时的 NaN 兜底（跨 L419-L420） | 同上 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L138 | `try; if czm_mesh !== nothing; case.czm_mesh = czm_mesh; end; if phase_type == PHASE_DISCHARGE; case.param.cell.v_l = V_limit; case.param.cell.v_h = 1.0e6; elseif phase_type == PHASE_CHARGE; ...; else; ...; end; solve_result = Solve(...); ... finally; case.param.cell.v_l = old_v_l; case.param.cell.v_h = old_v_h; end`（嵌套 3 层：try + if + if-elseif-else，跨 L138-L194） | 抽出 `with_phase_overrides(case, phase_type, V_limit) do ... end` 高阶函数；当前 try-finally + 多分支混合 |
| L227 | `if soc_init >= 0.0 && soc_init <= 1.0`（双 `&&` 链） | 接近阈值但未越界（2 个条件） |
| L284 | 阶段1（放电）+ 阶段2（静置1）的整体控制流，含多个嵌套 if-else（verbose + t_rest1 > 0 + reset_T_each_cycle），跨 L284-L356 共 72 行 | 抽出 `run_phase_with_rest(case, phase_name, t, I, V, current_state; rest_enabled)` helper 减少重复 |
| L321 | `if cycle_opt.t_rest1 > 0; if verbose; print("  [静置1] ..."); end; rest1_result = solve_phase(...); ...; if verbose; @printf(...); end; else; cycle_result.rest1 = PhaseResult(); cycle_result.rest1.duration = 0.0; cycle_result.rest1.V_start = get(current_state, "V", NaN); cycle_result.rest1.V_end = cycle_result.rest1.V_start; cycle_result.rest1.final_state = current_state; if verbose; println("  [静置1] 跳过"); end; end`（嵌套 2 层 + 多赋值，跨 L321-L356） | 抽出 `skip_rest_phase("rest1", current_state)` 函数返回 PhaseResult；rest2 分支（L393-L426）几乎完全重复 |
| L393 | 同上模式的 rest2 分支，与 L321-L356 几乎完全重复（跨 L393-L426） | 共享同一个 helper（如 `skip_rest_phase`），避免双份代码 |
| L268 | 循环主体约 180 行，4 个阶段 + 后处理 + 终止检查全部内联（跨 L268-L447） | 抽出 `run_one_cycle(case, cycle_opt, cycle, current_state, ...)` 函数返回 `(cycle_result, current_state, should_stop, soh_hit)` |
