# Parallelsolution.jl

## 文件状态

新增

## main分支

- **状态**: 不存在（`git show main:src/Parallelsolution.jl` 返回 `fatal: path exists on disk, but not in 'main'`）
- **行数**: N/A

## Parameters_Design分支

- **行数**: 619
- **主要函数列表**:

  | 行号  | 函数签名                                                                                                          | 说明                                                                                        |
  | --- | ------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
  | 7   | `_debug_check_prefactors(...)`                                                                                | 调试：检查电化学预因子中的 NaN/Inf                                                                     |
  | 33  | `_debug_check_coefficients(...)`                                                                              | 调试：检查单元系数有效性                                                                              |
  | 55  | `_debug_check_initial_voltage(...)`                                                                           | 调试：检查初始电压分布合理性                                                                            |
  | 83  | `_compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)`                                     | 计算电化学预因子（交换电流密度、OCV 等）                                                                    |
  | 119 | `_compute_element_coefficients(e, T_e, param, prefactors, T_ref; debug_mode)`                                 | 计算单个单元的分流系数（C1, C2, alpha, C5）                                                            |
  | 152 | `_compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref; debug_mode)`                                | 批量计算所有单元的分流系数                                                                             |
  | 162 | `_branch_voltage(coeff, I::Float64)`                                                                          | 计算单元分支电压 V_e = coeff.C1 + coeff.C2*I + coeff.alpha_p*T + coeff.alpha_n*T + coeff.C5*ln(I) |
  | 169 | `_branch_dVdI(coeff, I::Float64)`                                                                             | 计算分支电压对电流的导数 dV_e/dI                                                                      |
  | 179 | `_initialize_currents(ne, w, I_total, x_prev)`                                                                | 初始化电流分布（均匀或使用上次结果）                                                                        |
  | 194 | `_check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e; context)`                                 | 检查电压是否越界，返回异常标志                                                                           |
  | 226 | `_detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)`                                       | 检测达到截止电压的单元，返回截止信息                                                                        |
  | 294 | `_newton_iteration!(I_e, V, ne, w, I_total, coeffs; tol_V, tol_I, max_iters, active_mask)`                    | Newton-Raphson 迭代求解分流电流，含线搜索                                                              |
  | 391 | `_line_search(I_e, V, dI, dV, I_trial, ne; max_attempts)`                                                     | 线搜索算法，防止 Newton 步过大导致发散                                                                   |
  | 452 | `solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements)` | 主入口：分流求解器，返回 (variables, I_e, V_common)                                                   |


## 变更详情

### 新增函数

#### `_debug_check_prefactors(prefactor_n, prefactor_p, csn_av, csp_av, u_n_ref_val, u_p_ref_val, du_n_dT_val, du_p_dT_val, c_sigma, cn_surf, cp_surf, ce_n_gs, ce_p_gs)`

- **位置**: 第 7 行
- **功能**: 检查电化学预因子是否包含 NaN/Inf，打印诊断信息

#### `_debug_check_coefficients(e, has_nan_prefactor, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n_val_T, u_p_val_T)`

- **位置**: 第 33 行
- **功能**: 检查单元分流系数的有效性，报告异常值

#### `_debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, coeffs, I_total, ne)`

- **位置**: 第 55 行
- **功能**: 验证初始电压分布的合理性，检测全零/全相同等异常

#### `_compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)`

- **位置**: 第 83 行
- **功能**: 从电化学变量中提取计算分流所需的预因子：
  - 负极/正极交换电流密度 `j0_n`, `j0_p`
  - 开路电压 `u_n`, `u_p`
  - 熵系数 `du_n_dT`, `du_p_dT`
  - 电导率 `c_sigma`
  - 高斯点浓度等

#### `_compute_element_coefficients(e, T_e, param, prefactors, T_ref; debug_mode=false)`

- **位置**: 第 119 行
- **功能**: 为单元 `e` 计算线性化的分流系数：
  - `C1 = u_p - u_n + ...`（OCV 差值相关）
  - `C2 = R_ohm + eta 相关项`（欧姆电阻相关）
  - `alpha_p`, `alpha_n`（温度系数）
  - `C5`（对数修正项）
- **支持**: 高/低温路径切换，调试模式输出

#### `_compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref; debug_mode=false)`

- **位置**: 第 152 行
- **功能**: 遍历所有单元调用 `_compute_element_coefficients`，返回系数向量

#### `_branch_voltage(coeff, I::Float64) -> Float64`

- **位置**: 第 162 行
- **功能**: 根据分流系数和电流计算分支电压

#### `_branch_dVdI(coeff, I::Float64) -> Float64`

- **位置**: 第 169 行
- **功能**: 分支电压对电流的解析导数（用于 Newton 迭代的 Jacobian）

#### `_initialize_currents(ne, w, I_total, x_prev) -> Vector{Float64}`

- **位置**: 第 179 行
- **功能**: 初始化电流分布。若有前一步结果则沿用，否则按面积权重均匀分配

#### `_check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e; context="")`

- **位置**: 第 194 行
- **功能**: 检查电压向量是否在合理范围内，打印越界诊断

#### `_detect_cutoff_elements(coeffs, ne::Int, V_MIN::Float64, V_MAX::Float64, I_total::Float64, phi_scale::Float64)`

- **位置**: 第 226 行
- **功能**: 检测已达到截止电压的单元。返回 NamedTuple：
  - `cutoff_elements`: 截止单元列表
  - `cutoff_ocv`: 截止单元的 OCV
  - `active_mask`: BitVector 活跃单元掩码
  - `n_cutoff`: 截止数量
  - `nearest_cutoff_element`: 最近截止的单元
  - `nearest_cutoff_ocv`: 最近截止的 OCV
  - `margin_to_cutoff`: 到截止的最小裕度

#### `_newton_iteration!(I_e, V, ne, w, I_total, coeffs; tol_V, tol_I, max_iters, active_mask)`

- **位置**: 第 294 行
- **功能**: Newton-Raphson 迭代求解分流方程：
  - 约束：所有分支电压相等 `V_e = V_c`
  - 约束：总电流守恒 `sum(w_e * I_e) = I_total`
  - 使用 `_branch_dVdI` 构造 Jacobian
  - 内置线搜索（调用 `_line_search`）
  - 支持活跃掩码（排除截止/失效单元）

#### `_line_search(I_e, V, dI, dV, I_trial, ne; max_attempts=12)`

- **位置**: 第 391 行
- **功能**: 回溯线搜索，确保 Newton 步不导致电压发散

#### `solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements=nothing)`

- **位置**: 第 452 行
- **功能**: 分流求解器的主入口函数
- **工作流程**:
  1. 从 `variables` 中提取物理量
  2. 调用 `_compute_electrochemical_prefactors` 获取预因子
  3. 调用 `_compute_all_coefficients` 获取所有单元的分流系数
  4. 检测截止/失效单元
  5. 调用 `_initialize_currents` 初始化电流
  6. 调用 `_newton_iteration!` 迭代求解
  7. 将结果（I_e, V_common, cutoff 信息）写入 `variables` 并返回
- **返回**: `(variables, I_e::Vector{Float64}, V_common::Float64)`

### 修改函数

不适用（新增文件）

### 删除函数

不适用（新增文件）

## 依赖关系

### 该文件依赖哪些其他文件

- `SetCase.jl` — 使用 `Case` 结构体（`case.param`, `case.opt`, `case.mesh`）
- `SPMe.jl` — 通过 `variables` 字典间接依赖 SPMe 变量（OCV、过电位等）
- 无显式 `using` 或 `include`（通过 `JuBat.jl` 统一加载）

### 哪些文件依赖该文件

- `Solve.jl` — 调用 `solve_branch_currents_newton()`（在 `CallModel_MultiSPMe` 和 `CallModel` 的 `distributed2D` 分支中）
- `JuBat.jl` — export `solve_branch_currents_newton`

### 新增的外部依赖

无（仅使用 Julia 标准库和 `LinearAlgebra`）

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系

Parallelsolution.jl 是**电-热耦合的关键桥梁**，它解决了多单元并联时电流如何分布的问题：

- **multi-SPMe**: 为每个热单元分配独立的电流 `I_e`，这是多 SPMe 并行架构的前提条件
- **distributed2D**: 逐单元温度 `Te_prev` 影响分流系数（通过温度对 OCV、交换电流密度、电导率的影响）
- **CZM**: `deactivated_elements` 参数支持将 CZM 损伤单元从分流计算中排除

### 哪些变更是耦合相关的

整个文件都是为 multi-SPMe + distributed2D 耦合而新增的。核心耦合点：

- 温度对分流系数的影响（`_compute_element_coefficients` 中的 `T_e` 参数）
- 截止电压检测（`_detect_cutoff_elements`），与 Solve.jl 中的截止逻辑联动
- CZM 失效单元排除（`deactivated_elements` 参数）

### 哪些变更是独立的

- `_debug_check_`* 系列调试函数 — 纯诊断工具，可独立移除
- `_line_search` — 通用数值算法，可独立使用

