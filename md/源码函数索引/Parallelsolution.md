# Parallelsolution.jl

- **源文件**: `src/Parallelsolution.jl`
- **行数**: 454 行
- **函数/struct 计数**: 1 struct + 1 类型辅助 + 10 函数 = 12 条目
- **职责**: 非线性分流求解器——单元级电流分配（`solve_branch_currents`）、Newton 主循环（`newton_iteration`）、线搜索（`line_search`）、电压截止检测（`detect_cutoff_elements`）、分支电压模型（`branch_voltage` / `branch_dVdI`）、电化学系数（`compute_prefactors` / `compute_element_coefficients` / `compute_all_coefficients`）、电流初始化与边界检查
- **相关技术文档**: `md/09_分流求解器.md`、`md/08_逐单元算法.md`

## 数据结构

### `struct CutoffInfo` — L6-L16

截止检测结果（替代 `Dict{String,Any}`，注释明示）。

- `active_mask::BitVector`：true 表示单元活跃
- `n_cutoff::Int`、`cutoff_elements::Vector{Int}`、`cutoff_ocv::Vector{Float64}`、`cutoff_type::Vector{Int}`（1=放电截止，2=充电截止）
- `all_ocv::Vector{Float64}`
- `nearest_element::Int`、`nearest_ocv::Float64`、`margin::Float64`（最近截止预警）

## 函数清单

### `scalarize(x) -> Float64` — L20

将数组取首元素转为 Float64，标量直接转。`isa(x, Number) ? Float64(x) : Float64(x[1])`。

### `compute_prefactors(variables, param, mesh_ne, mesh_pe) -> NamedTuple` — L23-L56

计算电化学预因子（与单元温度无关的部分）。

- 输入：cn_surf / cp_surf / ce_n_gs / ce_p_gs（从 variables 读取）
- `prefactor_n = IntV(abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5, mesh_ne) / NE.thickness`（L30）
- 开路电位 `u_n_ref = NE.U(cn_surf)`、温度导数 `du_n_dT = NE.dUdT(cn_surf)`，标量化（L42-L45）
- 固相电导 `c_sigma = (NE.thickness / NE.sig + PE.thickness / PE.sig) / 3.0`（L48）
- 返回 NamedTuple（13 字段）

### `compute_element_coefficients(e, T_e, param, prefactors, T_ref) -> NamedTuple` — L59-L84

单单元电化学系数（温度相关）。

- Arrhenius 因子（L61-L62）：`arr_n = Arrhenius(NE.Eac_k, T_e)`
- 交换电流密度（L63-L64）：`j0_n = NE.k * arr_n * prefactors.prefactor_n`
- 电解液电导（L67-L69）：有效电导率 = `EL.kappa(ce0, T_e) * eps^brugg`
- 电解液电阻 `R_EL = NE.t/(3κ_ne) + SP.t/κ_sp + PE.t/(3κ_pe)`（L70）
- 温度修正 OCP（L73-L74）
- 5 个系数：`C1 = (u_p - u_n) + ...`、`C2 = 2T_e`、`alpha_p`、`alpha_n`、`C5 = R_EL + c_sigma`（L77-L81）

### `compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref) -> Vector{NamedTuple}` — L87-L93

循环调用 `compute_element_coefficients` 得到所有单元的系数向量。`param` 的 `e` 实际未在 `compute_element_coefficients` 中使用（参数签名保留但内部不用，见 [PLACEHOLDER]）。

### `branch_voltage(coeff, I::Float64) -> Float64` — L97-L101

分支电压模型：`V = C1 + C2*(asinh(α_p*I) - asinh(α_n*I)) - C5*I`。

### `branch_dVdI(coeff, I::Float64) -> Float64` — L104-L110

分支电压对电流的导数：`dV/dI = C2*(α_p/√(1+(α_p*I)²) - α_n/√(1+(α_n*I)²)) - C5`。

### `initialize_currents(ne, w, I_total, x_prev) -> Vector{Float64}` — L114-L126

初始化单元电流猜测。

- 有 `x_prev` 且长度匹配：`copy(x_prev)`（L115）
- 否则：`w .* I_total`（L115）
- 归一化到满足总电流约束（L118-L123）：`sx = sum(w .* I_e)`；非零则 `I_e .*= (I_total / sx)`，否则退回 `w .* I_total`

### `check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="") -> (ok, cutoff_type, V_phys, V_threshold_phys)` — L129-L143

检查电压边界。`V < V_MIN` 返回 `(false, 1, V_phys, V_MIN_phys)`（放电截止）；`V > V_MAX` 返回 `(false, 2, V_phys, V_MAX_phys)`（充电截止）；否则 `(true, 0, V_phys, 0.0)`。

参数 `I_total / w / I_e / context` 在函数体内未使用（见 [PLACEHOLDER]）。

### `detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale) -> CutoffInfo` — L161-L193

检测达到截止电压的单元。

- `OCV[e] = coeffs[e].C1`（L163）
- 截止判定（L170）：`(I_total < 0 && OCV[e] >= V_MAX) || (I_total > 0 && OCV[e] <= V_MIN)`（充电时 OCV ≥ V_MAX 满，放电时 OCV ≤ V_MIN 空）
- 最近截止预警（L183-L190）：`I_total > 0` 时取 OCV 最小单元；`I_total < 0` 时取 OCV 最大单元；`I_total == 0` 时 `nearest_element=0, nearest_ocv=NaN, margin=NaN`

### `newton_iteration(I_e, V, ne, w, I_total, coeffs; tol_V, tol_I, max_iters, active_mask) -> (V, converged, last_iter)` — L200-L294

Newton 迭代主循环（支持部分单元截止）。

- 若 `active_mask === nothing`：所有单元活跃（L209-L211）
- `n_active == 0` 直接返回（L219-L221）
- **迭代**（L223-L291）：
  - 残差 `F[e] = V_e - V`、雅可比 `dFdI[e] = branch_dVdI(...)`（L228-L230）
  - **奇异雅可比兜底**（L233-L235）：`abs(dFdI) < 1e-12` 时强制设为 `sign * 1e-12` 或 `-C5`
  - 收敛判定：`res_V <= tol_V && abs(res_I) <= tol_I`（L246-L249）
  - **全活跃 vs 部分活跃**两条路径（L252-L273）：全活跃用向量操作，部分活跃只对 `active_idx` 计算
  - 线搜索（L276-L277）：`λ == 0.0 && break`
  - 更新 `I_e` / `V`（L280-L290）

### `line_search(I_e, V, ΔI, ΔV, I_trial, ne; max_attempts) -> (λ, V_trial)` — L297-L321

回溯线搜索。

- `λ = 1.0` 起始，每轮 `λ *= 0.5`（L300-L318）
- 失败条件：`!isfinite(V_trial)` 或 `!isfinite(val) || abs(val) > 1e12`（L303、L309）
- 失败返回 `(0.0, V)`（L320）

### `solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements, D_elem) -> (variables, I_e, V)` — L358-L453

分流求解器主入口。

- **渐进式面积损失**（L361-L372）：`czm_area_loss_enabled && D_elem !== nothing` 时 `A_eff = areas .* effective_area_factor.(D_elem, threshold)`；debug_coupling 时 println 损失统计（L365-L369）
- CZM 失效掩码合并（L377-L382）
- **预因子 + 系数**（L385-L386）：`compute_prefactors(...) → compute_all_coefficients(...)`
- **截止检测**（L389-L394）：`detect_cutoff_elements(...) → ci`；`active_mask[e] = deactivated_mask[e] ? false : ci.active_mask[e]`
- 初始化电流，非活跃单元置零（L397-L400）
- 初始电压（L403-L404）：无活跃单元用 OCV 均值，否则用活跃单元 branch_voltage 均值
- Newton 迭代（L407-L410）
- **归一化**（L413-L420）：`sx = sum(w[e] * I_e[e] for e in active_idx)`，非零则缩放，否则按权重分配
- 边界检查（L423）
- 写入 variables 约 20 个键（L426-L452）：电流、电压、收敛状态、active_mask、cutoff_details、inactive_reason 等

## 省略项

无。所有 struct 与 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L366 | `println("  [AreaLoss][Weight] 超阈值=$(length(loss_idx))单元 \| D=$(...) \| factor=$(...) \| w∈[...]")`（跨 L366-L369，在 `case.opt.debug_coupling` 守卫下） | 调试输出：渐进式面积损失的权重变化与超阈值单元数 |

注：本文件其他位置无 `println` / `@show` / `@info` / `@warn`。`check_voltage_bounds` 的 `context` 参数虽然可作 debug 用途，但函数体未引用，见 [PLACEHOLDER]。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L59 | `function compute_element_coefficients(e, T_e, param, prefactors, T_ref)` 中参数 `e` 未在函数体内使用 | 占位参数：签名保留单元索引 `e` 但实际不使用（系数与 e 无关，仅与 T_e 相关）；签名维护时易误以为按 e 分支 |
| L129 | `check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="")` 中 `I_total`、`w`、`I_e`、`context` 四个参数未在函数体内使用（跨 L129-L143） | 占位参数：4 个未使用参数，签名臃肿；调用者必须传值但实际上无效果，易误导 |
| L183 | `nearest_element, nearest_ocv, margin = 0, NaN, NaN`（`I_total == 0` 静置时的默认值） | 占位 NaN：静置阶段无方向，最近截止预警用 NaN 兜底；下游若未判 NaN 直接算术运算会传染 NaN |
| L220 | `if n_active == 0; return V, true, 0; end`（无活跃单元时直接返回 converged=true） | 占位语义：无活跃单元返回 converged=true 是"什么都不做即满足约束"的语义，但调用者可能误以为真的收敛到一个解 |
| L233 | `if abs(dFdI[e]) < 1e-12; dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5; end`（跨 L233-L235，奇异雅可比兜底） | 防止奇异：magic number `1e-12`；`sign(x) != 0.0` 的判断逻辑古怪（sign 返回 -1/0/1，与 0.0 比较是 Float 比较），当 dFdI 恰为 0 时用 `-C5` 兜底 |
| L255 | `abs(denom) < 1e-12 && break`（denom 接近 0 时静默退出 Newton 循环） | 静默 break：无 warning，converged 保持 false；调用方仅通过 `converged` 字段判断，但此路径下 last_iter 可能误指示迭代次数 |
| L309 | `if !isfinite(val) \|\| abs(val) > 1e12`（线搜索 val 范围检查） | magic number `1e12`：硬编码上限，注释无说明 |
| L362 | `A_eff = areas .* effective_area_factor.(D_elem, case.opt.czm_area_loss_threshold)`（渐进式面积损失，依赖外部函数 `effective_area_factor`） | 跨文件依赖：`effective_area_factor` 未在本文件定义；若未正确 export 会导致 UndefVarError |
| L403 | `V = isempty(active_idx) ? mean([coeffs[e].C1 for e in 1:ne]) : mean([branch_voltage(coeffs[e], I_e[e]) for e in active_idx])`（跨 L403-L404，无活跃单元用 OCV 均值占位） | 占位语义：无活跃单元的 V 仅为形式值，后续 check_voltage_bounds 会判定为截止；但变量名 V 在调用者视角可能被误用 |
| L414 | `sx = sum(w[e] * I_e[e] for e in active_idx); if abs(sx) > 1e-12; sf = I_total / sx; for e in active_idx; I_e[e] *= sf; end`（跨 L414-L416，归一化） | magic number `1e-12`：阈值未参数化 |
| L417 | `elseif abs(I_total) > 1e-12 && !isempty(active_idx); w_sum = sum(w[active_idx]); w_sum > 0 && (for e in active_idx; I_e[e] = I_total * w[e] / w_sum; end); end`（跨 L417-L420，sx 为 0 时的回退分配） | 兜底回退：当 sx 为 0 但 I_total 非零时按权重重新分配；语义合理但隐式 |
| L430 | `variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5`（magic number 3.0 / 3.5） | 占位编码：3.0=converged / 3.5=未收敛，但无注释说明，且用 Float 而非 Symbol/Bool |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L161 | `function detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)` 内部：`active_mask = trues(ne); OCV = [...]; cutoff_elements = Int[]; cutoff_ocv = Float64[]; cutoff_type = Int[]; for e in 1:ne; cutoff = (I_total < 0 && OCV[e] >= V_MAX) \|\| (I_total > 0 && OCV[e] <= V_MIN); if cutoff; ...; end; end; ...; if I_total > 0; idx = argmin(OCV); ...; elseif I_total < 0; idx = argmax(OCV); ...; end`（多层独立 if + 多 `&&`，跨 L161-L193） | 抽出 `find_active_units(coeffs, ne, V_MIN, V_MAX, I_total) -> (active_mask, cutoff_lists)` 与 `nearest_cutoff(OCV, I_total, V_MIN, V_MAX) -> (...)` 两个 helper |
| L170 | `cutoff = (I_total < 0 && OCV[e] >= V_MAX) \|\| (I_total > 0 && OCV[e] <= V_MIN)`（`||` + `&&` × 2 混合，单行 4 个逻辑运算符） | 抽出 `is_cutoff(I_total, OCV_e, V_MIN, V_MAX)` 函数；当前单行表达阅读成本高 |
| L200 | `function newton_iteration(...)`：单函数 ~95 行，含 if-active_mask 多分支、残差/J 计算循环、收敛判定、ΔV/ΔI 计算（全活跃 vs 部分活跃两份代码）、线搜索、更新（全活跃 vs 部分活跃两份代码）（跨 L200-L294） | 拆分为 `compute_residuals(I_e, V, coeffs, ne)` / `newton_step_all_active(...)` / `newton_step_partial_active(...)` 三个函数；当前全活跃/部分活跃分支几乎完全对称但代码重复 |
| L227 | `for e in 1:ne; V_e = branch_voltage(coeffs[e], I_e[e]); F[e] = V_e - V; dFdI[e] = branch_dVdI(coeffs[e], I_e[e]); if abs(dFdI[e]) < 1e-12; dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5; end; end`（跨 L227-L235，for 内含 if + 三元嵌套） | 抽出 `safe_dFdI(coeff, I_e) -> Float64` helper 处理奇异雅可比兜底 |
| L358 | `function solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements, D_elem)`（9 个位置参数 + 2 个关键字 = 11 参数） | 接近复杂度阈值；可考虑用 `BranchSolverInput` struct 封装 (areas, Te_prev, x_prev, deactivated_elements, D_elem) 减少参数数量 |
| L361 | `if case.opt.czm_area_loss_enabled && D_elem !== nothing; A_eff = ...; w = ...; loss_idx = findall(...); if !isempty(loss_idx) && case.opt.debug_coupling; factors = ...; println(...); end; else; w = areas ./ sum(areas); end`（跨 L361-L372，嵌套 2 层 if + 调试 println + 双赋值路径） | 抽出 `compute_branch_weights(areas, D_elem, opt) -> w` helper |
| L389 | `ci = detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale); active_mask = copy(ci.active_mask); for e in 1:ne; deactivated_mask[e] && (active_mask[e] = false); end; active_idx = findall(active_mask)`（跨 L389-L394，多步串联 + for 循环修改 mask） | 接近阈值但未越界；可改为 `active_mask = ci.active_mask .& .!deactivated_mask`（向量化） |
| L397 | `I_e = initialize_currents(ne, w, I_total, x_prev); for e in 1:ne; !active_mask[e] && (I_e[e] = 0.0); end`（跨 L397-L400，初始化 + 循环清零） | 接近阈值；可向量化 `I_e[.!active_mask] .= 0.0` |
| L413 | `sx = sum(w[e] * I_e[e] for e in active_idx); if abs(sx) > 1e-12; sf = I_total / sx; for e in active_idx; I_e[e] *= sf; end; elseif abs(I_total) > 1e-12 && !isempty(active_idx); w_sum = sum(w[active_idx]); w_sum > 0 && (for e in active_idx; I_e[e] = I_total * w[e] / w_sum; end); end`（跨 L413-L420，嵌套 if + 双 for + 生成器表达式） | 抽出 `normalize_currents!(I_e, active_idx, w, I_total)` helper |
| L438 | `inactive_reason = zeros(Float64, ne); for e in 1:ne; inactive_reason[e] = deactivated_mask[e] ? 2.0 : (!active_mask[e] ? 1.0 : 0.0); end; variables["thermal2D inactive_reason"] = inactive_reason; variables["thermal2D voltage_in_bounds"] = ...; ...`（跨 L438-L452，连续 ~15 行 variables 字典写入） | 接近阈值；可抽出 `write_branch_solver_outputs!(variables, I_e, V, ...)` 函数 |
