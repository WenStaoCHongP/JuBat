# CycleSolver.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 730
- 路径: `src/CycleSolver.jl`

### 主要结构体

| 结构体 | 说明 |
|--------|------|
| `PhaseResult` | 单阶段（充电/静置/放电）求解结果，包含时间、电压、容量、温度、损伤、终止原因、最终状态 |
| `CycleResult` | 单个循环结果，包含4个PhaseResult（放电/静置1/充电/静置2）及循环汇总（容量、库伦效率、损伤、温度） |
| `CyclingResult` | 全部循环结果汇总，包含每循环指标向量、SOH轨迹、最终CZM网格 |

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `solve_phase(case, phase_type, t_max, I_current, V_limit, initial_state; czm_mesh, czm_params, dt_range)` | L110 | 求解单个阶段（放电/充电/静置），含CZM更新 |
| `solve_cycling(case, cycle_opt, czm_mesh; verbose, save_detailed)` | L209 | 主循环求解器：放电->静置1->充电->静置2 |
| `compute_cs0_from_soc(param_dim, soc)` | L486 | 根据SOC计算正负极初始锂浓度 |
| `apply_initial_soc!(case, param_dim, soc)` | L524 | 应用初始SOC并重新归一化参数 |
| `_compute_czm_effective_params(case, param_dim)` | L551 | 计算CZM有效材料参数（E_eff, nu_eff, alpha_eff, beta_n, beta_p） |
| `_compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)` | L581 | 计算CZM应变输入（温度变化、SOC变化） |
| `_update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)` | L658 | 调用CZM求解器更新损伤状态 |
| `_ensure_multi_spme_layout!(case)` | L708 | 确保multi_spme布局信息已初始化 |

## 功能描述

本文件实现了完整的充放电循环仿真框架，是JuBat项目循环老化仿真的顶层调度模块。主要功能包括：

1. **单阶段求解**（`solve_phase`）：
   - 根据阶段类型（放电/充电/静置）配置电压限制和电流方向
   - 调用 `Solve(case; initial_state, return_final_state)` 执行时间积分
   - 阶段末调用CZM求解器更新损伤状态
   - 通过 `_postprocess_phase_result` 后处理阶段结果

2. **循环求解**（`solve_cycling`）：
   - 循环顺序：放电 -> 静置1 -> 充电 -> 静置2
   - 状态在阶段间通过 `final_state` 字典传递（y, T_nodes, V, t_global）
   - 支持每循环温度重置（`reset_T_each_cycle`）和充电前温度重置（`reset_T_before_charge`）
   - 静置阶段：当 t_rest > 0 时执行零电流锂扩散求解；t_rest = 0 时跳过
   - SOH监控：基于放电容量衰减判断电池健康状态
   - 终止条件：SOH低于阈值（`czm_soh_threshold`）或超过50% CZM单元断裂

3. **SOC初始化**（`compute_cs0_from_soc`/`apply_initial_soc!`）：
   - 负极：`theta_n = theta_0_n + SOC * (theta_100_n - theta_0_n)`
   - 正极：`theta_p = theta_0_p - SOC * (theta_0_p - theta_100_p)`
   - 更新维度参数后重新归一化

4. **CZM接口**（`_update_czm_damage!`）：
   - 从当前variables提取温度场和SOC分布
   - 计算有效材料参数（厚度加权平均）
   - 调用 `solve_czm_step` 进行力学求解
   - 同步损伤状态到调用者的czm_mesh

5. **布局管理**（`_ensure_multi_spme_layout!`）：
   - 确保状态向量布局信息（ne, n_chem, nT, n_total）已初始化
   - 解决跨阶段状态传递时布局信息可能丢失的问题

## 依赖关系

### 该文件依赖
- `src/Solve.jl` — `Solve`主求解器、`CallModel`模型调用
- `src/CzmSolve.jl` — `solve_czm_step`、`get_damage_statistics`
- `src/czm.jl` — `assemble_thermal_chemical_load`（间接）
- `src/PostProcessing.jl` — `_postprocess_phase_result`、`_postprocess_cycle_result!`等辅助函数
- `src/Variables.jl` — `Case`、`CycleOption`、`PhaseType`等类型
- `src/SetParams.jl` — `Params`、`Cohesive`参数类型
- `src/Jellyrollmodel.jl` — 间接通过mesh数据使用

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("CycleSolver.jl")`（L31）
- `example/czm_cycle_example.jl` — 调用 `solve_cycling`
- 用户脚本 — 通过 `JuBat.solve_cycling` 接口调用

## 耦合分析

本文件是JuBat项目的**循环仿真调度中心**，协调所有物理场在时间循环中的交互：

- **与multi-SPMe耦合**：
  - 每个阶段调用 `Solve` 时，multi-SPMe求解器在每个时间步计算各单元的电化学状态
  - `_ensure_multi_spme_layout!` 保证跨阶段状态向量布局正确
  - 静置阶段（PHASE_REST）设置电流为0，SPMe仅执行颗粒扩散

- **与distributed2D热模型耦合**：
  - 温度场通过 `T_nodes` 在 `initial_state`/`final_state` 中跨阶段传递
  - 支持 `reset_T_each_cycle` 在每个循环开始时重置温度场
  - 支持 `reset_T_before_charge` 在充电前重置温度场

- **与CZM耦合**：
  - 每个阶段结束调用 `_update_czm_damage!` 更新损伤
  - 损伤状态通过 `czm_mesh.damage_states` 跨循环累积
  - SOH终止判断基于CZM断裂状态（`check_fracture_criterion`）
  - `_compute_czm_strain_inputs` 将热-电化学状态转化为力学应变输入

- **与分流求解器耦合**：
  - CZM损伤状态通过 `get_active_elements` 影响下一阶段的电流分配
  - 断裂界面导致对应热单元变为非活跃，不分配电流

## 后续变更 (2026-04-01)

- **移除 `solve_cycling` 中的调试打印语句**（约 13 行）：
  - 删除了放电阶段结束后的 `V_out` / `y_len` 打印语句
  - 删除了静置阶段结束后的 `V_out` / `y_len` 打印语句
  - 删除了充电阶段开始前的 `V_in` / `y_len` 打印语句
- **核心循环逻辑不变**：放电 -> 静置1 -> 充电 -> 静置2 的阶段调度、状态传递、CZM 更新等均未受影响。
