# CzmSolve.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 515
- 路径: `src/CzmSolve.jl`

### 主要结构体

| 结构体 | 说明 |
|--------|------|
| `CZMResult` | CZM求解结果，包含位移场(displacement)、损伤(damage)、法向/切向牵引力、法向/切向分离、收敛标志、迭代次数、残差范数 |

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `newton_raphson_czm(czm_mesh, F_ext, E_eff, nu_eff, cohesive_params, param; alpha_eff, beta_n, beta_p, dT_elem, Dsoc_n_elem, Dsoc_p_elem, max_iter, tol, u0, n_load_steps)` | L30 | 带载荷子步的Newton-Raphson求解器 |
| `solve_czm_step(czm_mesh, F_ext, E_eff, nu_eff, cohesive_params, param, u_prev; ..., iter_method)` | L154 | 统一CZM求解入口，支持3种迭代方法 |
| `get_damage_statistics(czm_mesh)` | L382 | 统计损伤：max_D, mean_D, min_D, n_fractured, fraction_damaged, total_accumulated |
| `check_fracture_criterion(czm_mesh; threshold)` | L407 | 检查断裂准则：平均损伤 >= threshold 或断裂单元数 > 50% |
| `reset_damage_states(czm_mesh)` | L429 | 重置所有损伤状态为初始值 |
| `accumulate_cycle_damage(czm_mesh, cycle_damage_increment)` | L452 | 跨循环累积损伤增量 |
| `czm_output_to_variables(czm_mesh, result, variables)` | L496 | 将CZMResult转换为variables字典格式 |

## 功能描述

本文件实现了CZM（内聚力模型）的非线性求解器，提供三种迭代方法处理材料软化导致的收敛困难：

1. **载荷子步法**（`newton_raphson_czm`，iter_method="load_substep"）：
   - 将总载荷分为n_load_steps个子步逐步施加
   - 每个子步内进行Newton-Raphson迭代
   - 位移增量限制：`norm(du) > 1e-6` 时缩放
   - 收敛后更新损伤状态（`update_damage`）

2. **基本Newton-Raphson**（iter_method="basic"）：
   - 标准Newton-Raphson迭代，无载荷子步
   - 使用相对残差收敛准则（`R_norm / R_norm_0 < tol`）
   - 适合小载荷步或损伤较轻的情况

3. **弧长法**（iter_method="arc_length"）：
   - 通过约束方程控制载荷-位移路径追踪
   - 使用二次方程求解载荷增量 delta_lambda
   - 选择与预测值最接近的根确保路径连续性
   - 适合材料软化后载荷下降的路径追踪

4. **损伤管理**：
   - `get_damage_statistics`：统计损伤分布
   - `check_fracture_criterion`：基于平均损伤或断裂比例判断全局失效
   - `reset_damage_states`：重置损伤（用于测试）
   - `accumulate_cycle_damage`：跨循环累积疲劳损伤

5. **后处理**（`czm_output_to_variables`）：将位移、损伤、牵引力、分离位移写入variables字典，便于PostProcessing使用。

## 依赖关系

### 该文件依赖
- `src/czm.jl` — `assemble_coupled_system`、`assemble_thermal_chemical_load`、`apply_bc_czm`、`identify_bc_nodes_czm`
- `src/Materialmatrix.jl` — `update_damage`、`bilinear_traction_state`（间接）
- `src/SetMesh.jl` — `CohesiveMesh`、`Mesh`
- `src/SetParams.jl` — `Cohesive`参数类型
- `src/Variables.jl` — `Case`类型

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("CzmSolve.jl")`（L9）
- `src/CycleSolver.jl` — 调用 `solve_czm_step`（通过 `_update_czm_damage!`）、`get_damage_statistics`、`czm_output_to_variables`、`accumulate_cycle_damage`

## 耦合分析

本文件是CZM损伤演化的**求解引擎**，在电-热-CZM耦合中承担力学求解角色：

- **与CycleSolver耦合**：`solve_czm_step` 是循环求解器中CZM更新的核心入口。每个充电/放电/静置阶段结束时，CycleSolver调用 `_update_czm_damage!` -> `solve_czm_step` 更新损伤状态。

- **与热模型耦合**：`dT_elem` 从热模型的温度场提取，转化为热应力载荷驱动CZM。

- **与multi-SPMe耦合**：`Dsoc_n_elem`/`Dsoc_p_elem` 从各单元SPMe模型的SOC提取，扩散应变是CZM载荷来源之一。

- **与分流求解器耦合**：损伤状态通过 `check_fracture_criterion` 影响分流求解器——断裂的界面导致对应热单元不分配电流。

- **迭代方法选择**：通过 `case.opt.czm_iter_method` 控制迭代方法。"basic"适合快速计算，"load_substep"适合一般情况，"arc_length"适合严重软化和后峰路径追踪。
