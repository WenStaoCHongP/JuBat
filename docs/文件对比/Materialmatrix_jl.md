# Materialmatrix.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 376
- 路径: `src/Materialmatrix.jl`

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)` | L16 | 计算高斯点热容权重：rho_c = sum(fk * rho * heat_Q) |
| `thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy)` | L30 | 计算高斯点各向异性导热系数（k_xx, k_xy, k_yy），径向串联+周向并联 |
| `bilinear_traction_state(delta_n, delta_t, damage_state, cohesive_params)` | L68 | 双线性牵引-分离本构：计算牵引力并返回更新损伤状态 |
| `bilinear_traction(delta_n, delta_t, damage_state, cohesive_params; update)` | L150 | 双线性牵引（可选择性更新损伤状态） |
| `bilinear_tangent(delta_n, delta_t, damage_state, cohesive_params)` | L168 | 双线性切线刚度矩阵 dT/d_delta (2x2) |
| `update_damage(damage_states, separations, cohesive_params)` | L262 | 批量更新损伤状态 |
| `compute_gap_conductance(D, delta_n, cohesive)` | L287 | 计算界面有效导热系数（含损伤耦合） |
| `compute_element_gap_conductance(czm_mesh, elem_idx, cohesive)` | L322 | 获取指定CZM单元的界面导热系数 |
| `get_fractured_elements(czm_mesh)` | L332 | 获取所有断裂（D >= 0.99）的CZM单元索引 |
| `get_active_elements(czm_mesh, mesh_data)` | L345 | 获取非活跃热单元列表（受CZM断裂影响） |
| `compute_all_gap_conductances(czm_mesh, cohesive)` | L368 | 计算所有CZM单元的界面导热系数 |

## 功能描述

本文件实现了JuBat项目的两大材料模型体系：热材料矩阵和CZM本构模型。主要功能包括：

### 一、热材料矩阵

1. **热容权重**（`thermal_capacity_weights_2d`）：
   - 基于5层权重矩阵 `fks[ne, 5]` 计算每个单元的等效体积热容
   - `rho_c = f_NE*rho_NE*Q_NE + f_SP*rho_SP*Q_SP + f_PE*rho_PE*Q_PE + f_PCC*rho_PCC*Q_PCC + f_NCC*rho_NCC*Q_NCC`
   - 乘以高斯点权重得到组装用系数

2. **各向异性导热系数**（`thermal_anisotropic_conductivity_2d`）：
   - **径向导热**（层间串联）：`lambda_r = 1 / sum(f_k / lambda_k)`，热阻串联模型
   - **周向导热**（层间并联）：`lambda_t = sum(f_k * lambda_k)`，面积加权平均
   - **坐标变换**：基于高斯点位置角度theta，将极坐标各向异性转换为笛卡尔坐标
   - `k_xx = lr*cos^2 + lt*sin^2`，`k_xy = (lt-lr)*sin*cos`，`k_yy = lr*sin^2 + lt*cos^2`

### 二、CZM双线性本构模型

3. **双线性牵引-分离律**（`bilinear_traction_state`）：
   - 支持两种模式：
     - `model1`：仅法向分离控制损伤（`delta_eff = max(0, delta_n)`），切向不损伤
     - `mix`（默认）：混合模式（BK准则），`delta_eff = sqrt(delta_n^2 + delta_t^2)`
   - 损伤变量：`D = delta_c * (delta_eff - delta_0) / (delta_eff * (delta_c - delta_0))`
   - 法向牵引力：`T_n = (1-D) * K_n * delta_n`（delta_n >= 0），受压时不损伤
   - 切向牵引力：model1中 `T_t = K_t * delta_t`（不损伤），mix中 `T_t = (1-D) * K_t * delta_t`

4. **切线刚度**（`bilinear_tangent`）：
   - 计算dT/d_delta的2x2矩阵
   - 加载阶段（delta_eff > delta_max_hist）：考虑损伤演化对刚度的影响
   - 卸载阶段：使用当前损伤值D降低刚度

5. **界面热阻模型**（`compute_gap_conductance`）：
   - 有效导热系数 h_eff 考虑损伤D和分离位移delta_n
   - `h_eff = 1 / (h_c0*(1-D) + k_air / (delta + 2*beta*lambda_m))`
   - 损伤增大 -> 热阻增大 -> 层间导热减弱

### 三、损伤-热耦合工具

6. **断裂管理**：
   - `get_fractured_elements`：返回D >= 0.99的CZM单元
   - `get_active_elements`：基于CZM断裂状态和热单元映射，确定活跃热单元
   - `compute_all_gap_conductances`：批量计算界面导热系数

## 依赖关系

### 该文件依赖
- `src/SetParams.jl` — `Params`、`Cohesive`、`DamageState`等类型
- `src/Variables.jl` — `CohesiveMesh`结构体
- `src/SetMesh.jl` — 高斯积分数据结构

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("Materialmatrix.jl")`（L22）
- `src/czm.jl` — 调用 `bilinear_traction_state`、`bilinear_tangent`（在`assemble_czm_system`中）
- `src/CzmSolve.jl` — 调用 `update_damage`（在Newton-Raphson迭代中）
- `src/ThermalDistributed.jl` — 调用 `thermal_capacity_weights_2d`、`thermal_anisotropic_conductivity_2d`
- `src/Parallelsolution.jl` — 调用 `get_active_elements`、`compute_gap_conductance`

## 后续变更 (2026-04-20)

- **`get_active_elements` 修复**: `ne = mesh_data.ne` 改为 `ne = length(mesh_data.is_inner_layer)`
  - 修复原因：`mesh_data` 可能没有 `ne` 字段，使用 `is_inner_layer` 向量长度获取单元数更可靠
  - 确保与 `MeshGeometry` 结构兼容（`is_inner_layer` 是 `MeshGeometry` 的必选字段）
