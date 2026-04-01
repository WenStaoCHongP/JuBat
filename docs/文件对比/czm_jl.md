# czm.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 504
- 路径: `src/czm.jl`

### 主要结构体

| 结构体 | 说明 |
|--------|------|
| `CohesiveElement <: AbstractCohesiveElement` | 内聚力单元，包含4个节点（底面/顶面各2个）、单元长度、层间界面索引 |
| `DamageState <: AbstractDamageState` | 损伤状态变量，包含D（损伤变量[0,1]）、历史最大法向/切向/等效分离位移、断裂标志、累积损伤 |

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `create_czm_mesh(thermal_mesh::Mesh, param_dim; tol)` | L54 | 基于热网格创建内聚力网格：检测界面节点对、生成CohesiveElement、初始化DamageState |
| `assemble_czm_system(czm_mesh, u, cohesive_params; damage_states)` | L154 | 组装内聚力单元的全局刚度矩阵和内力向量，返回K_coh, f_int_coh, separations, tractions |
| `assemble_bulk_stiffness(czm_mesh, E_eff, nu_eff)` | L283 | 组装固体Q4单元的平面应力刚度矩阵 |
| `assemble_thermal_chemical_load(czm_mesh, E_eff, nu_eff, alpha_eff, beta_n, beta_p, dT_elem, Dsoc_n_elem, Dsoc_p_elem)` | L349 | 组装热-化学载荷向量：epsilon_0 = alpha*deltaT + beta_n*deltaSoc_n + beta_p*deltaSoc_p |
| `assemble_coupled_system(czm_mesh, u, E_eff, nu_eff, cohesive_params; ...)` | L405 | 组装耦合系统：K_total = K_bulk + K_coh, f_int_total = f_int_bulk + f_int_coh |
| `assemble_coupled_system_full(czm_mesh, u, ...; ...)` | L427 | 完整耦合组装：包含热化学载荷，残差 R = F_ext + F_thermo_chem - f_int_total |
| `apply_bc_czm(K, F; bc_nodes, bc_dofs, bc_vals)` | L450 | 施加边界条件（罚函数法） |
| `identify_bc_nodes_czm(czm_mesh, param; opt)` | L484 | 识别CZM网格的边界节点（内外壁固定约束） |

## 功能描述

本文件实现了内聚力模型（CZM）的网格创建和有限元系统组装，是电-热-CZM全耦合仿真的力学核心模块。主要功能包括：

1. **内聚力网格创建**（`create_czm_mesh`）：从热分析网格出发，通过坐标重合检测识别层间界面节点对。在外螺旋节点（n_out）与内螺旋节点（n_in）重合处创建CohesiveElement（4节点2D单元）。每个CZM单元关联一个DamageState用于追踪损伤演化。

2. **内聚力系统组装**（`assemble_czm_system`）：对每个内聚力单元进行高斯积分，计算：
   - 分离位移（法向delta_n、切向delta_t）通过旋转变换矩阵R从全局坐标转换
   - 牵引力通过双线性本构模型（`bilinear_traction_state`）计算
   - 切线刚度通过`bilinear_tangent`计算
   - 使用Newton-Cotes积分（`NCweight`）沿单元长度积分

3. **固体单元组装**（`assemble_bulk_stiffness`）：标准Q4平面应力有限元，弹性矩阵 `D = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]`。

4. **热-化学载荷**（`assemble_thermal_chemical_load`）：将温度变化（alpha*deltaT）和锂浓度变化（beta*deltaSoc）转化为等效节点力向量，`F = integral B^T * D * epsilon_0 dOmega`。

5. **耦合系统**（`assemble_coupled_system_full`）：组装完整的力学系统，K_total = K_bulk + K_coh，残差 = 外力 + 热化学力 - 内力。

## 依赖关系

### 该文件依赖
- `src/SetMesh.jl` — `Mesh`结构体、`GetGS`、`IntQ4`、`NCweight`、`PickElement`
- `src/Variables.jl` — `Case`类型、`CohesiveMesh`结构体
- `src/SetParams.jl` — `Params`、`Cohesive`参数类型
- `src/Materialmatrix.jl` — `bilinear_traction_state`、`bilinear_tangent`本构模型
- `src/ring.jl`（间接） — `identify_boundary_nodes`边界识别函数

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("czm.jl")`（L8）
- `src/CzmSolve.jl` — 调用 `assemble_czm_system`、`assemble_bulk_stiffness`、`assemble_thermal_chemical_load`、`assemble_coupled_system`
- `src/CycleSolver.jl` — 通过 `create_czm_mesh` 创建CZM网格

## 耦合分析

本文件是 **CZM求解器** 的核心，连接了热模型和电化学模型：

- **与distributed2D热模型耦合**：
  - `create_czm_mesh` 接收热网格作为输入，在热单元的界面节点处插入内聚力单元
  - 热模型的温度场通过 `dT_elem` 传入 `assemble_thermal_chemical_load`，计算热应力载荷
  - 损伤状态通过界面热阻模型（`compute_gap_conductance`）反向影响热传导

- **与multi-SPMe耦合**：
  - SPMe计算的锂浓度分布通过 `Dsoc_n_elem`、`Dsoc_p_elem` 传入载荷计算
  - 扩散应变 `beta * deltaSoc` 是CZM载荷的重要组成部分

- **与分流求解器耦合**：
  - CZM断裂的单元（`fractured`）通过 `get_active_elements` 影响电流分布
  - 断裂单元对应的热单元被视为非活跃，分流求解器不向其分配电流

- **关键数据流**：温度场 + SOC分布 -> 热化学应变 -> CZM载荷 -> 损伤演化 -> 界面热阻变化 -> 温度场更新
