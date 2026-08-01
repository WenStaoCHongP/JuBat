# Solve.jl

- **源文件**: `src/Solve.jl`
- **行数**: 472 行
- **函数/struct 计数**: 3 个函数（无独立 struct；含 1 个嵌套 closure）
- **职责**: 主求解器——`Solve`（含热模型独立路径与多 SPMe 主循环）、系统矩阵快照（`RecordMatrix!`）、误差估计（`ErrorEstimation`）
- **相关技术文档**: `md/09_分流求解器.md`、`md/08_逐单元算法.md`、`md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。`Solve` 返回 `Dict{String, Any}`（经 `PostProcessing` 包装），可附加 `final_state::Dict{String, Any}`。

## 函数清单

### `Solve(case; initial_state, return_final_state, thermal_variables, thermal_update_fn, thermal_record, polar_mesh_data, czm_snapshots, czm_cycle, czm_phase) -> Dict` — L1-L440

主求解器，含两条主要路径：

**路径 A：纯热模型**（`case.opt.model == "thermal"`，L2-L71）

- 解析 `solveType` 为 `theta`（L13-L21），仅 Crank-Nicolson / forward / backward
- 时间序列 `times = range(t0, step=dt, ...)`（L27）
- 初始矩阵装配按 `thermalmodel` 三分支（L32-L40）：`ring2D_polar` → `ThermalPolar2D_Ring`；`ring2D` → `ThermalDistributed2D_Ring` + `ThermalRing2D_BC`；其他 → `ThermalDistributed2D` + `ThermalDistributed2D_BC`
- 主循环（L42-L67）：`(MT_new - θ·dt·KT_new) T = (MT_new + (1-θ)·dt·KT_old) T_old + dt·(θ·FT_new + (1-θ)·FT_old)`
- 返回命名元组 `(time, T_nodes, T_hist)`（L69-L70）

**路径 B：电化学（含多 SPMe）主路径**（L73-L439）

- 时间步无量纲化（L73-L77）
- **初始化分支**（L86-L127）：
  - `y0_input === nothing`：按 `multi_spme_enabled` 选 `ModelInitialisation_MultiSPMe` / `ModelInitialisation`（L87-L91）
  - `isa(y0_input, Dict)` 且含 `"y"`：取出 vec；否则回退到模型初始化（L94-L108）
  - 多 SPMe 模式下额外校验 `expected_multi_len = ne * n_chem + nT`（L116），不匹配时 `@warn` 并回退（L122-L125）
- **solveType → theta 映射**（L128-L136），错误信息含拼写 `opt.solve_type`（应为 `case.opt.solveType`，L135）—— 见 [PLACEHOLDER]
- **预分配**（L142-L153）：`num_estimated = round(Int64, (t_end - t0)/dt * 1.5)`；`max_steps = multi_spme_enabled ? 50000 : 100000`；超出时 `@warn`（L147-L150）
- **嵌套 closure** `accumulate_callmodel_timing!`（L161-L167）：累加 4 个 timing 字段到 `timing_totals` Dict
- **初始 CallModel**（L169-L171）：取得 `M_old, K_old, F_old, variables, y_phi`，打印初始化信息（L173-L174）
- **初始步**（L178-L186）：`dt_init = 1e-8`；`y_c = (M_old - K_old*dt_init) \ (M_old*y0 + F_old*dt_init)`；`y_old = vcat(y_c, y_phi)`；若 `jacobi == "constant"` 调 `RecordMatrix!`
- **截止追踪变量**（L190-L196）：`first_cutoff_detected / time / element / ocv`、`total_cutoff_count`、`termination_reason`（默认 `"time_limit"`）
- **CZM 状态**（L198-L203）：`czm_active = czm_enabled && czm_mesh !== nothing`；缺 `czm_layout` 时调 `CzmLayout(...)` 构造
- **主循环**（L206-L360）：
  - 电化学步：`CallModel(...) → Mt, Kt, Ft`；`y_c = Mt \ (Kt*y_old + Ft)`；`y_new = vcat(y_c, y_phi)`（L208-L215）
  - 多 SPMe：从 `y_c` 末尾 `nT` 个 DOF 提取 `T_nodes` 写回 variables（L218-L225）
  - **误差估计 + 自适应 dt**（L226-L259）：`error_y > 2*dtThreshold` 且 `dt >= dt_min*4` 则回退 dt；否则记录变量、按 `error_y` 分层调整 dt（×2 / 重置 / ÷2）
  - **CZM 损伤更新**（L268-L315）：每 `czm_update_interval` 步触发；调 `update_czm_damage!(...)`；`czm_snapshots !== nothing` 时构造并 push `CZMSnapshot`；异常被 `try/catch` 捕获并 `@debug` 输出（L310-L312）
  - **截止电压检测**（L318-L350）：多 SPMe 时读取 `n_cutoff_elements`，记录首个截止单元；`n_cutoff >= ne_total` 则 `termination_reason = "all_elements_cutoff"` 并 break；整体 V 检测在 L353-L359
- **后处理**（L362-L439）：
  - 调 `PostProcessing`（L366）
  - 写入 timing 字段（合计、平均、占比，L371-L392）
  - `case.opt.debug_coupling` 时 `@printf` 输出各阶段耗时（L394-L405）
  - 写入 `termination_reason / first_cutoff_*` 字段（L408-L415）
  - 附加温度历史（L418-L426）—— 用 try/catch 静默失败（L424-L426，见 [PLACEHOLDER]）
  - `return_final_state` 时构造 `final_state` Dict（L429-L438）

### `RecordMatrix!(case, M, K) -> case` — L442-L450

将单模型的全局 M、K 矩阵按电极切片缓存到 `case.param.NE.M_d / K_d` 与 `case.param.PE.M_d / K_d`。

- 切片范围：`M[1:l_np, 1:l_np]` → NE；`M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]` → PE（L445-L448）
- 仅在 `case.opt.jacobi == "constant"` 时由主求解器调用（L184-L185）

### `ErrorEstimation(case, y_old, y_new, coeff) -> Float64` — L452-L471

计算时间步误差估计 `error_y`。

- **SPM / SPMe 模式**（L454-L455）：`error_y = norm(y_new - y_old) / norm(y_old) * coeff`（整体相对误差）
- **其他模型**（L456-L468）：分 5 个字段独立计算（cn / cp / cel / phi_pp / phi_el），取最大值
- 用于自适应时间步控制（在 Solve 主循环与 `solve_phase_with_export` 中调用）

## 省略项

无。所有 function 均有独立条目。嵌套 closure `accumulate_callmodel_timing!` 因局部作用域不单列。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L174 | `println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")`（每次 Solve 调用都打印初始化电压与终止时间） | 进度提示：用户可见的初始化 banner，无 verbose 守卫 |
| L188 | `print("start to solve the problem \n")`（每次进入主循环前打印） | 进度提示：用户可见的开始消息，无 verbose 守卫 |
| L311 | `@debug "CZM damage update failed at step $czm_step_count: $e"`（CZM 异常捕获时） | 调试日志：CZM 更新异常时低优先级输出，避免阻塞主循环 |
| L395 | `println("\n[Solve-Timing] CallModel 阶段累计耗时..."); @printf("  SPMe solve ...")` 等 4 行（跨 L395-L404，在 `case.opt.debug_coupling` 守卫下） | 调试输出：阶段累计耗时与占比，用于瓶颈定位 |
| L427 | `print("finish the simulation\n")`（仿真结束打印） | 进度提示：用户可见的完成消息，无 verbose 守卫 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L135 | `error("Error: $(opt.solve_type) difference scheme has not been implemented!\n")`（变量名拼写错误：`opt.solve_type` 应为 `case.opt.solveType`） | 引用未定义符号：若运行到此分支会抛 UndefVarError；正常路径不会触发，但代码维护时易混淆 |
| L144 | `max_steps = multi_spme_enabled ? 50000 : 100000`（magic number 阈值） | 硬编码上限：注释说明"避免内存溢出"，但阈值未参数化；超大仿真或未来机型可能再次触发 |
| L178 | `dt_init = 1e-8`（初始半步 dt） | magic number：注释无说明，数值经验性选取；不同物理尺度下可能需要调整 |
| L196 | `termination_reason = "time_limit"`（默认字符串占位） | 默认值：若主循环因非时间原因退出但未更新该字段，可能误报；正常路径会通过 break 时显式赋值 |
| L310 | `try; u_czm_new, czm_converged = update_czm_damage!(...); catch e; @debug ...; end`（跨 L310-L312，静默捕获 CZM 异常） | 静默 try-catch：CZM 失败被吞掉，仅 `@debug` 输出；运行中可能掩盖严重问题，损伤演化停滞 |
| L418 | `try; ... result["thermal2D final temperature at nodes [K]"] = ...; catch; # non-fatal; end`（跨 L418-L426，静默捕获温度历史附加失败） | 静默 try-catch：注释明示"non-fatal"；条件 `per_element_spme && thermalmodel == "distributed2D" && !isempty(T_nodes_carry)` 失败时静默跳过，运行中难诊断为何字段缺失 |
| L369 | `call_count_safe = max(timing_call_count, 1)`（防止除零） | 防御性兜底：正常路径 `timing_call_count >= 1`，但保留 max 防 0 除 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L86 | `if y0_input === nothing; if multi_spme_enabled; ...; else; ...; end; else; if isa(y0_input, Dict); y_from_state = get(...); if y_from_state !== nothing; y0 = vec(...); else; if multi_spme_enabled; ...; else; ...; end; end; else; y0 = vec(...); end; if multi_spme_enabled; ne = ...; ...; if length(y0) != expected_multi_len; @warn ...; ...; end; end; end`（嵌套 4 层 if + 重复 multi_spme 分支，跨 L86-L127） | 抽出 `init_y0(case, initial_state, multi_spme_enabled)` 函数；当前三处重复调用 `ModelInitialisation_MultiSPMe(case)` / `ModelInitialisation(case)` |
| L118 | `if case.layout === nothing; case.layout = MultiSPMeLayout(...); end; if length(y0) != expected_multi_len; @warn ...; y0 = ModelInitialisation_MultiSPMe(case); end`（跨 L118-L126，连续 `===nothing` 检查 + 长度校验 + 回退） | 接近 ≥2 阈值；layout 构造与状态校验可分离 |
| L228 | `if error_y > 2 * case.opt.dtThreshold && case.opt.dtType == "auto" && dt >= dt_min * 4`（3 个 `&&` 链） | 接近 ≥3 阈值；可抽出 `should_reject_step(error_y, dt, dt_min, opt)` helper |
| L255 | `if t + dt > RunTime[vt] && t < RunTime[vt]`（双时间比较链） | 接近阈值但未越界；语义清晰 |
| L327 | `if multi_spme_enabled; n_cutoff = Int(variables["thermal2D n_cutoff_elements"]); if n_cutoff > 0 && !first_cutoff_detected; first_cutoff_detected = true; ...; end; total_cutoff_count = n_cutoff; ne_total = size(...); if ne_total > 0 && n_cutoff >= ne_total; termination_reason = ...; break; end; end`（嵌套 3 层 + 多重 `&&`，跨 L327-L350） | 抽出 `check_element_cutoff!(variables, ne_total, first_cutoff_state)` 函数返回 `(should_break, reason, updated_state)` |
| L419 | `if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D" && !isempty(T_nodes_carry)`（3 个 `&&` 链） | 抽出 `should_export_thermal_history(case, T_nodes_carry)` helper；此类条件在多个文件重复出现 |
