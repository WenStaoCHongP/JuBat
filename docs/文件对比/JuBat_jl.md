# JuBat.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 30
- 主要内容: 模块入口文件，定义 `module JuBat`，导入依赖包，include 所有源文件，export 公开符号

### include 文件列表（按顺序）
1. `Option.jl`
2. `SetMesh.jl`
3. `SetParams.jl`
4. `SetCase.jl`
5. `Assemble.jl`
6. `ElectrodeDiffusion.jl`
7. `ElectrolyteDiffusion.jl`
8. `ElectrodePotential.jl`
9. `ElectrolytePotential.jl`
10. `SPM.jl`
11. `SPMe.jl`
12. `P2D.jl`
13. `Solve.jl`
14. `PostProcessing.jl`
15. `Tools.jl`
16. `Thermal.jl`
17. `Variables.jl`
18. `Initialisation.jl`
19. `sP2D.jl`
20. `Citation.jl`
21. `Mechanical.jl`

### 导出列表（3组）
- `Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell`
- `Mesh1D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables`
- `SPM, Solve, SPMe`
- `Arrhenius, IntV`

### 依赖包
`LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator`

## Parameters_Design分支
- 行数: 74 (+44, +147%)
- 主要内容: 大幅扩展的模块入口文件，新增 11 个 include 文件和大量 export 声明

### 新增 include 文件（11个）
1. `czm.jl` - 内聚力区域模型
2. `CzmSolve.jl` - CZM 求解器
3. `Parallelsolution.jl` - 分流求解器（**从 Solve.jl 之前加载**）
4. `Tools.jl` - **移动到 Solve.jl 之前**（确保工具函数先定义）
5. `Materialmatrix.jl` - 材料矩阵
6. `ThermalDistributed.jl` - 二维分布式热模型
7. `ThermalPolar2D.jl` - 极坐标 FVM 求解器
8. `mechanical.jl` - 应力计算（**替代 Mechanical.jl**，注意大小写变化）
9. `Jellyrollmodel.jl` - Jellyroll 几何与网格
10. `ring.jl` - 环形网格
11. `CycleSolver.jl` - 循环求解器
12. `CycleData.jl` - 循环数据导出

### 删除的 include 文件（3个）
1. `sP2D.jl` - 已移除
2. `Citation.jl` - 已移除
3. `Mechanical.jl` - 被 `mechanical.jl`（小写）替代

### 新增依赖包
- `Statistics` - 用于温度场统计计算
- `Printf` - 用于格式化输出

### 新增导出（按类别分组）

**2D 网格与几何**
- `Mesh2D, IntQ4`
- `jellyroll_collector_seed_mesh, jellyroll_element_properties, jellyroll_tab_node_indices, edge_boundary, jellyroll_element_centers`
- `ring_mesh`

**多 SPMe 模型**
- `SPMe_element, ModelInitialisation`
- `ModelInitialisation_MultiSPMe, MultiSPMe_extract_element_state, MultiSPMe_get_thermal_dofs, MultiSPMe_update_state`

**热模型**
- `setup_thermal2D_mesh`
- `thermal2D_volume_average_temperature`
- `ThermalDistributed2D, ThermalDistributed2D_BC, ThermalDistributed2D_Ring, ThermalRing2D_BC, ThermalPolar2D_Ring`
- `identify_boundary_nodes, apply_convection_bc, apply_cool_method`
- `ThermalModel, ThermalLumpedModel, ThermalDistributed2DModel`

**分流求解**
- `compute_heat_sources, compute_heat_sources_with_czm, solve_branch_currents_newton`

**力学**
- `thermal_diffusion_stress_2D`

**CZM (内聚力模型) - 完整导出**
- `CohesiveElement, CohesiveMesh, DamageState, CZMResult`
- `create_czm_mesh, compute_separation`
- `bilinear_traction, bilinear_tangent, update_damage`
- `assemble_czm_system, assemble_coupled_system, assemble_bulk_stiffness`
- `assemble_thermal_chemical_load, assemble_coupled_system_full`
- `apply_bc_czm, identify_bc_nodes_czm`
- `newton_raphson_czm, solve_czm_step`
- `get_damage_statistics, check_fracture_criterion, reset_damage_states, accumulate_cycle_damage`
- `czm_output_to_variables`
- `compute_gap_conductance, compute_element_gap_conductance`
- `get_fractured_elements, get_active_elements, compute_all_gap_conductances`

**循环求解**
- `CycleOption, PhaseType, PHASE_CHARGE, PHASE_REST, PHASE_DISCHARGE`
- `PhaseResult, CycleResult, CyclingResult`
- `solve_phase, solve_cycling, plot_cycling_results`
- `compute_cs0_from_soc, apply_initial_soc!`

**数据导出**
- `TimeStepData, CycleExportData`
- `solve_phase_with_export, solve_cycle_with_export`
- `export_cycle_data_to_csv, load_cycle_data_from_csv`

## 变更详情

### 核心变更
1. **include 顺序重组**: `Tools.jl` 移到 `Solve.jl` 之前，`Parallelsolution.jl` 也移到 `Solve.jl` 之前，确保依赖函数先定义
2. **新增 11 个 include 文件**: 涵盖 CZM、热模型、Jellyroll、循环求解等完整功能链
3. **删除 3 个 include 文件**: `sP2D.jl`, `Citation.jl`, `Mechanical.jl`
4. **导出列表从 3 组扩展到约 70+ 个符号**: 几乎所有新增模块的核心函数都被导出
5. **新增依赖**: `Statistics`, `Printf`

## 依赖关系

### 被依赖关系
- 作为模块入口，所有其他文件都被此文件 include
- include 顺序决定了函数/类型的可见性

### 关键顺序约束
- `czm.jl` 和 `CzmSolve.jl` 必须在 `SetCase.jl` 之后（依赖 Case 类型）
- `Tools.jl` 必须在 `Solve.jl` 之前（提供 IntQ4 等工具函数）
- `Parallelsolution.jl` 必须在 `Solve.jl` 之前（分流求解器被主求解器调用）
- `CycleData.jl` 必须在 `CycleSolver.jl` 之后（依赖 PhaseType 等类型）

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 是（核心入口文件）

此文件是整个 Jellyroll 多物理场耦合框架的入口点。新增的所有 include 文件和 export 声明都直接服务于 multi-SPMe + distributed2D + CZM 耦合架构。关键耦合路径：
- 多 SPMe 并行 -> 分流求解器 -> 热模型 -> CZM 损伤 -> 界面热阻 -> 反馈到热模型
