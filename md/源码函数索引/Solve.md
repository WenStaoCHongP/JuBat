# Solve.jl

- **源文件**: `src/Solve.jl`
- **行数**: 448 行
- **函数/struct 计数**: 3 个函数（无独立 struct；含 1 个嵌套 closure）
- **职责**: 主求解器——`Solve`（含热模型独立路径与多 SPMe 主循环）、系统矩阵快照（`RecordMatrix!`）、误差估计（`ErrorEstimation`）
- **相关技术文档**: `md/08_逐单元算法.md`、`md/09_分流求解器.md`、`md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。`Solve` 返回 `Dict{String, Any}`（经 `PostProcessing` 包装），可附加 `final_state::Dict{String, Any}`。

## 函数清单

### `Solve(case; initial_state, return_final_state, thermal_variables, thermal_update_fn, thermal_record, polar_mesh_data, czm_snapshots, czm_cycle, czm_phase) -> Dict` — L1-L417

主求解器，含两条主要路径：

**路径 A：纯热模型**（`case.opt.model == "thermal"`，L2-L75）

- `thermal_variables` 为必填 kwarg，缺失即 `error`（L3-L5）
- 解析 `solveType` 为 `theta`（L13-L21），仅 Crank-Nicolson / forward / backward
- 时间序列 `times = range(t0, step=dt, ...)`（L27）
- 初始矩阵装配按 `thermalmodel` 三分支（L32-L42）：`ring2D_polar` → `ThermalPolar2D_Ring`；`ring2D` → `ThermalDistributed2D_Ring` + `ThermalRing2D_BC`；其他 → `ThermalDistributed2D` + `ThermalDistributed2D_BC`
- 主循环（L44-L71）：`(MT_new - θ·dt·KT_new) T = (MT_new + (1-θ)·dt·KT_old) T_old + dt·(θ·FT_new + (1-θ)·FT_old)`
- 返回命名元组 `(time, T_nodes, T_hist)`（L73）

**路径 B：电化学（含多 SPMe）主路径**（L76-L416）

- 时间步无量纲化（L76-L81）
- **初始化分支**（L89-L118，严格契约：外部状态非法即 throw，不回退）：
  - `y0_input === nothing`：按 `multi_spme_enabled` 选 `ModelInitialisation_MultiSPMe` / `ModelInitialisation`（L90-L95）
  - `isa(y0_input, Dict)`：缺 `"y"` 键或非数组直接 `throw(ArgumentError)`（L97-L99）；多 SPMe 模式下校验 `expected_multi_len = ne * n_chem + nT`，不匹配 `throw(DimensionMismatch)`（L101-L117）
  - `case.layout === nothing` 时惰性构造 `MultiSPMeLayout`（L110-L112）
- **solveType → theta 映射**（L119-L127），错误信息含拼写 `opt.solve_type`（应为 `case.opt.solveType`，L126）—— 见 [PLACEHOLDER]
- **预分配**（L129-L141）：`num_estimated = round(Int64, (t_end - t0)/dt * 1.5)`；`max_steps = multi_spme_enabled ? 50000 : 100000`；超出时 `@warn`（L138-L140）
- **嵌套 closure** `accumulate_callmodel_timing!`（L152-L158）：累加 4 个 timing 字段到 `timing_totals` Dict
- **初始 CallModel**（L160-L162）：取得 `M_old, K_old, F_old, variables, y_phi`，打印初始化信息（L164-L165）
- **热场持久化**（L167-L176）：多 SPMe 时从 variables 取 `thermal2D temperature at nodes` 存 `T_nodes_carry`；缺键/非数组/长度不符均 throw（严格契约）
- **初始步**（L178-L186）：`dt_init = 1e-8`；`y_c = (M_old - K_old*dt_init) \ (M_old*y0 + F_old*dt_init)`；`y_old = vcat(y_c, y_phi)`；若 `jacobi == "constant"` 调 `RecordMatrix!`
- **截止追踪变量**（L190-L196）：`first_cutoff_detected / time / element / ocv`、`total_cutoff_count`、`termination_reason`（默认 `"time_limit"`）
- **CZM 状态**（L198-L203）：`czm_active = czm_enabled && czm_mesh !== nothing`；缺 `czm_layout` 时调 `CzmLayout(...)` 构造
- **主循环**（L206-L341）：
  - 电化学步：`CallModel(...) → Mt, Kt, Ft`；`y_c = Mt \ (Kt*y_old + Ft)`；`y_new = vcat(y_c, y_phi)`（L208-L215）
  - 多 SPMe：长度校验后从 `y_c` 末尾 `nT` 个 DOF 提取 `T_nodes` 写回 variables 并更新 `T_nodes_carry`（L218-L224）
  - **误差估计 + 自适应 dt**（L225-L255）：`error_y > 2*dtThreshold` 且 `dt >= dt_min*4` 则回退 dt；否则记录变量、按 `error_y` 分层调整 dt（×2 / 重置 / ÷2）
  - **CZM 损伤更新**（L266-L296）：每 `czm_update_interval` 步触发；调 `update_czm_damage!(case, variables, T_nodes_carry)`；`czm_snapshots !== nothing` 时构造并 push `CZMSnapshot`；`try/finally` 仅保证计时累计，**异常正常传播**（无静默 catch）
  - **截止电压检测**（L299-L340）：多 SPMe 时读取 `n_cutoff_elements`，记录首个截止单元（L312-L319）；`n_cutoff >= ne_total` 则 `termination_reason = "all_elements_cutoff"` 并 break（L327-L330）；整体 V 检测在 L333-L340
- **后处理**（L343-L416）：
  - 调 `PostProcessing`（L347）
  - 写入 timing 字段（合计、平均、占比，L348-L373）
  - `case.opt.debug_coupling` 时打印各阶段耗时（L375-L386）
  - 写入 `termination_reason / first_cutoff_*` 字段（L388-L396）
  - 附加温度历史（L398-L403）：`per_element_spme && thermalmodel == "distributed2D"` 时直接写入（无 try/catch 兜底）
  - `return_final_state` 时构造 `final_state` Dict（L406-L415，含 `"y"/"T_nodes"/"V"/"t_global"`）

### `RecordMatrix!(case, M, K) -> case` — L419-L427

将单模型的全局 M、K 矩阵按电极切片缓存到 `case.param.NE.M_d / K_d` 与 `case.param.PE.M_d / K_d`。

- 切片范围：`M[1:l_np, 1:l_np]` → NE；`M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]` → PE（L422-L425）
- 仅在 `case.opt.jacobi == "constant"` 时由主求解器调用（L184-L186）

### `ErrorEstimation(case, y_old, y_new, coeff) -> Float64` — L429-L448

计算时间步误差估计 `error_y`。

- **SPM / SPMe 模式**（L431-L432）：`error_y = norm(y_new - y_old) / norm(y_old) * coeff`（整体相对误差）
- **其他模型**（L433-L446）：分 5 个字段独立计算（cn / cp / cel / phi_pp / phi_el），取最大值
- 用于自适应时间步控制（在 Solve 主循环与 `solve_phase_with_export` 中调用）

## 省略项

无。所有 function 均有独立条目。嵌套 closure `accumulate_callmodel_timing!` 因局部作用域不单列。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L165 | `println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")`（每次 Solve 调用都打印初始化电压与终止时间） | 进度提示：用户可见的初始化 banner，无 verbose 守卫 |
| L188 | `print("start to solve the problem \n")`（每次进入主循环前打印） | 进度提示：用户可见的开始消息，无 verbose 守卫 |
| L376-L385 | `println("\n[Solve-Timing] ...")` + 4 行 `@printf`（在 `case.opt.debug_coupling` 守卫下） | 调试输出：阶段累计耗时与占比，用于瓶颈定位 |
| L404 | `print("finish the simulation\n")`（仿真结束打印） | 进度提示：用户可见的完成消息，无 verbose 守卫 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L126 | `error("Error: $(opt.solve_type) difference scheme has not been implemented!\n")`（变量名拼写错误：`opt.solve_type` 应为 `case.opt.solveType`） | 引用未定义符号：若运行到此分支会抛 UndefVarError；正常路径不会触发，但代码维护时易混淆 |
| L135 | `max_steps = multi_spme_enabled ? 50000 : 100000`（magic number 阈值） | 硬编码上限：注释说明"避免内存溢出"，但阈值未参数化；超大仿真或未来机型可能再次触发 |
| L178 | `dt_init = 1e-8`（初始半步 dt） | magic number：注释无说明，数值经验性选取；不同物理尺度下可能需要调整 |
| L196 | `termination_reason = "time_limit"`（默认字符串占位） | 默认值：主循环因非时间原因退出但未更新该字段时可能误报；当前路径在 break/结尾均显式赋值 |
| L350 | `call_count_safe = max(timing_call_count, 1)`（防止除零） | 防御性兜底：正常路径 `timing_call_count >= 1`，但保留 max 防 0 除 |

> 2026-08-18 复核：早期的两处静默兜底已随严格契约移除——外部状态长度不匹配改为
> `throw(DimensionMismatch)`（L114-L116），CZM 更新改为 `try/finally`（异常传播，
> L271-L294），温度历史附加不再 try/catch（L398-L403）。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L89-L118 | 初始化分支：外部状态与多 SPMe 布局的组合判定（嵌套 if + 长度校验 throw） | 可抽出 `init_y0(case, initial_state, multi_spme_enabled)`；当前三处 `ModelInitialisation*` 调用模式重复度已低于旧版（回退分支删除后） |
| L110-L116 | `case.layout === nothing` 惰性构造 + `expected_multi_len` 校验 + throw | 接近 ≥2 阈值；layout 构造与状态校验可分离 |
| L227 | `if error_y > 2 * case.opt.dtThreshold && case.opt.dtType == "auto" && dt >= dt_min * 4`（3 个 `&&` 链） | 接近 ≥3 阈值；可抽出 `should_reject_step(error_y, dt, dt_min, opt)` helper |
| L254-L256 | `if t + dt > RunTime[vt] && t < RunTime[vt]`（双时间比较链） | 接近阈值但未越界；语义清晰 |
| L308-L331 | 单元级截止检测：嵌套 if + 多重 `&&`（`n_cutoff > 0 && !first_cutoff_detected`、`ne_total > 0 && n_cutoff >= ne_total`） | 抽出 `check_element_cutoff!(variables, ne_total, first_cutoff_state)` 函数返回 `(should_break, reason, updated_state)` |
| L399 | `if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D"`（2 个 `&&` 链） | 可抽出 `should_export_thermal_history(case)` helper；此类条件在多个文件重复出现 |
