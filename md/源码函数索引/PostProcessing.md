# PostProcessing.jl

- **源文件**: `src/PostProcessing.jl`
- **行数**: 351 行
- **函数/struct 计数**: 11 个函数（无独立 struct）
- **职责**: 结果提取与去归一化还原、循环阶段/周期结果汇总、SOH/容量保持率计算、循环仿真终止判定、循环可视化绘图
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/13_耦合验证方案.md`

## 数据结构

本文件无独立 struct 定义。所有结果以 `Dict{String, Any}` 形式返回。

## 函数清单

### `PostProcessing(case, variables, v::Int64) -> Dict` — L1-L117

从无量纲 `variables` 字典提取前 `v` 个时间步，通过 `case.param.scale` / `case.param_dim.scale` 还原为物理单位。

- **通用还原**（L3-L12）：时间 `×t0`、电压 `×phi`、电流 `×I1C`、温度 `×T_ref`、应力 `×E_n/E_p`、位移 `×r0`
- **模型分发**（L13-L47）：
  - SPM（L13-L17）：粒子浓度 `×cn_max/cp_max`
  - SPMe（L18-L27）：同 SPM + 交换电流密度 `×j`、过电位 `×phi`、电解液浓度 `×ce`
  - P2D（L28-L47）：同 SPMe + 电位类变量 `×phi` + 高斯点 OCV/界面电流（注意 L43-L46 对这些键未切片 `[:, 1:v]`，直接全量乘 scale——若历史数组含未填充尾部会包含错误数据）
- **热还原**（L49-L98）：
  - lumped（L49-L50）：`×q`
  - distributed2D（L51-L98）：分层热源（`Q_rxn_NE` 等 11 个，`×q`）、节点温度（`×T_ref`）、单元温度（节点平均 `element_nodal_mean` 后 `×T_ref`）、热源场直传、单元级变量直传（9 个键）、截止/激活信息（5 个键）
- **CZM 还原**（L99-L115）：
  - 标量（L100-L105）：`D_max`、`D_mean`、`n_fractured` 无量纲直传
  - 分离量（L103-L104）：`×δ_czm`（重设计 v2，修正原误用 `scale.L`）
  - 位移（L107-L108）：`×L`
  - 牵引（L110-L111）：`×σ_czm`（重设计 v2，修正原误用 `E_n/E_p`）
  - 分离位移（L113-L114）：`×δ_czm`（重设计 v2，修正原误用 `scale.r0`）

### `_phase_termination_symbol(phase_type, reason) -> Symbol` — L119-L127

内部辅助：根据 phase_type 和终止原因返回 `:time` 或 `:voltage` Symbol。REST 阶段恒返回 `:time`；电压截止类原因返回 `:voltage`；其他返回 `:time`。

### `_state_concentration_variance(case, y_state) -> (var_n, var_p)` — L129-L157

内部辅助：从状态向量计算粒子浓度方差。用于 REST 阶段松弛度评估。

- `y_state === nothing` 时返回 `(0.0, 0.0)`（L130-L132）
- 多 SPMe 模式（L135）：逐单元提取 `cs_n_e`、`cs_p_e`，跨单元计算方差（L139-L152）
- 单 SPMe 模式（L154-L156）：直接从 `y` 切片计算方差

### `_postprocess_phase_result(case, phase_type, solve_result, initial_state, I_current, t_start, D_max_init, D_mean_init, czm_mesh) -> Dict` — L159-L217

内部辅助：从单阶段（charge/discharge/rest）的 `solve_result` 汇总出 `PhaseResult` 兼容的 Dict。

- 提取 `duration`（L162）、`V_start`/`V_end`（L164-L165）、`T_max`（L169-L175）、`T_mean_end`（L176）
- `final_state["t_global"] = t_start + duration`（L179）
- CZM 损伤：从 `czm_mesh` 取 `get_damage_statistics`（L183-L187）
- `capacity = abs(I_current) * duration / 3600`（L190，单位 Ah）
- REST 阶段扩散松弛度（L192-L202）：初始/终态方差对比，`cs_relaxation_n/p` 百分比

### `_postprocess_cycle_result!(cycle_result, charge_result, discharge_result, rest1_result, rest2_result, czm_mesh)` — L219-L234

内部辅助：汇总单循环结果到 `cycle_result`。计算 `coulombic_efficiency`（L222-L223）、从 `czm_mesh` 取损伤统计（L225-L230）、取四阶段 `T_max` 最大值（L232）。

### `_append_cycle_result!(result, cycle, cycle_result; save_detailed=false)` — L236-L252

内部辅助：将单循环结果追加到 `result`（CyclingResult 兼容对象）的各向量字段。`save_detailed=true` 时保留完整 `cycle_result` 对象。

### `_update_soh_and_capacity!(result, cycle, cycle_result, initial_capacity) -> (initial_capacity, current_soh)` — L254-L263

内部辅助：更新 SOH 与初始容量。

- 首次循环（L255-L258）：以首次放电容量锁定 `initial_capacity`
- `current_soh = capacity_discharge / initial_capacity`（L260），`initial_capacity = 0` 时回退为 `1.0`
- 返回更新后的 `(initial_capacity, current_soh)` 元组（注意：`initial_capacity` 在函数内被重新绑定，需通过返回值传回调用者）

### `_print_cycle_summary(cycle, cycle_result, current_soh)` — L265-L269

内部辅助：`@printf` 打印单循环摘要（充/放容量、CE、D_max、SOH）。

### `_check_cycle_termination(cycle, cycle_result, czm_mesh, current_soh, soh_threshold; verbose=true) -> (terminated, soh_terminated)` — L271-L290

内部辅助：检查循环终止条件。

- SOH ≤ 阈值且 `cycle > 1`（L274-L280）：`@warn` + 返回 `(true, true)`
- CZM 断裂比例 > 50%（L282-L287）：`@warn` + 返回 `(true, false)`
- 否则返回 `(false, false)`

### `_print_cycling_summary(result, initial_capacity, soh_terminated)` — L292-L308

内部辅助：`println` + `@printf` 打印循环仿真总体摘要。含完成循环数、初始/最终容量、SOH、损伤、终止原因。

### `plot_cycling_results(result; save_path="output/") -> Plots.Plot` — L310-L350

绘制循环仿真结果 4 子图：容量衰减、损伤演化、库仑效率、温度历史。组合为 2×2 布局，保存 PNG 到 `save_path`。

- `isdir(save_path) || mkdir(save_path)`（L311）
- 4 个子图分别 `savefig`（L321/L329/L336/L342）
- 组合图 `plot(p1, p2, p3, p4, layout=(2,2), size=(1200,900))`（L344）
- 库仑效率图 `ylims=(95, 105)`（L335，硬编码范围）

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L276 | `@warn "SOH降至$(round(current_soh*100, digits=1))%，低于阈值$(round(soh_threshold*100, digits=1))%，终止循环"` | SOH 终止告警；结构化但 `verbose=true` 默认开启，用户可见 |
| L284 | `@warn "超过50%的内聚力单元断裂，提前终止循环"` | CZM 断裂终止告警；同上，verbose 门控 |
| L294 | `println("\n" * "="^60)` / `println("循环仿真完成")` / `println("="^60)` / `@printf(...)` 系列（跨 L294-L306） | 循环仿真结束摘要打印；非调试用途，是用户可见的最终报告输出 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L164 | `V_start = get(initial_state, "V", isempty(voltage_hist) ? NaN : voltage_hist[1])`（NaN 兜底） | 空 voltage_hist 时 V_start = NaN；若下游消费者未处理 NaN 会传播。实践中 voltage_hist 为空意味着仿真未启动 |
| L174 | `case.param_dim.cell.T0`（T_max 退化分支：T_hist_K 与 T_nodes_K 均空时回退到初始温度） | 兜底：无温度数据时用 T0 作为合理估计；通常仅在纯电化学（无热耦合）仿真中出现 |
| L176 | `T_mean_end = !isempty(T_nodes_K) ? mean(T_nodes_K) : case.param_dim.cell.T0`（同 L174 模式） | 同 L174 的兜底模式 |
| L222 | `cycle_result.coulombic_efficiency = cycle_result.capacity_charge > 0 ? 100.0 * capacity_discharge / capacity_charge : 0.0`（capacity_charge = 0 时回退 0.0，跨 L222-L223） | 兜底：零充电容量时 CE = 0% 而非 NaN；通常仅在不完整循环中出现 |
| L260 | `current_soh = initial_capacity > 0 ? cycle_result.capacity_discharge / initial_capacity : 1.0`（initial_capacity = 0 时回退 1.0） | 兜底：initial_capacity = 0 时 SOH = 100%；与 L222 同模式 |
| L293 | `final_soh = initial_capacity > 0 && result.n_cycles > 0 ? result.capacity_discharge[end] / initial_capacity : 1.0`（多重兜底） | 兜底：无循环或零容量时 SOH = 100%；打印用，不影响数据 |
| L335 | `ylims=(95, 105)`（库仑效率图 y 轴范围硬编码） | magic number：假设 CE 在 95-105% 范围；异常 CE（如首循环不一致）会被裁剪不可见 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L169 | `T_max = if !isempty(T_hist_K); maximum(T_hist_K); elseif !isempty(T_nodes_K); maximum(T_nodes_K); else; case.param_dim.cell.T0; end`（嵌套 if-elseif-else + 多 `!isempty` 判定，跨 L169-L176） | 抽出 `compute_T_max(T_hist_K, T_nodes_K, T0_default) -> Float64` helper；同模式在 L176 T_mean_end 重复 |
| L222 | `cycle_result.coulombic_efficiency = cycle_result.capacity_charge > 0 ? 100.0 * cycle_result.capacity_discharge / cycle_result.capacity_charge : 0.0`（单行三元 + 除零保护，跨 L222-L223） | 阈值边缘（单条件），可读性可接受 |
| L282 | `if czm_mesh !== nothing && cycle_result.n_fractured > 0.5 * czm_mesh.n_cohesive`（2 个 `&&` + 浮点比较） | 阈值边缘（2 个条件），可读性可接受 |
| L293 | `final_soh = initial_capacity > 0 && result.n_cycles > 0 ? result.capacity_discharge[end] / initial_capacity : 1.0`（2 个 `&&` + 三元） | 阈值边缘（2 个条件），可读性可接受 |
| L13 | `if case.opt.model == "SPM"; ...8 个赋值...; elseif case.opt.model == "SPMe"; ...13 个赋值...; elseif case.opt.model == "P2D" \|\| case.opt.model == "sP2D"; ...20 个赋值...; end`（嵌套 1 层 + 单函数 ~35 行多模型分发，跨 L13-L47） | 抽出 `restore_chem_variables(case, variables, v) -> Dict` 独立函数；当前与热/CZM 还原混合在单函数中 |
| L85 | `for key in ["thermal2D element current", "thermal2D eta_n_e", "thermal2D eta_p_e", "thermal2D element soc_n", "thermal2D element soc_p", "thermal2D element voltages", "thermal2D element OCV", "thermal2D dUdT_n_e", "thermal2D dUdT_p_e"]; result[key] = variables[key][:, 1:v]; end`（9 元素字符串数组循环直传，跨 L85-L90） | 抽出常量 `DIMENSIONLESS_PASSTHROUGH_KEYS = [...]`；当前硬编码在函数体内，与 StandardVariables 键表不同步风险高 |
