# JuBat.jl

- **源文件**: `src/JuBat.jl`
- **行数**: 89 行
- **函数/struct 计数**: 0 个独立函数；0 个 struct 定义
- **职责**: JuBat 模块入口文件，负责按依赖顺序 `include` 全部子模块，并通过 `export` 对外暴露公共 API
- **相关技术文档**: `md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

本文件无独立函数定义。全部内容为 `include` 加载顺序与 `export` 列表，归类如下：

### include 顺序（按文件出现顺序）

| 行号 | 文件 | 用途分组 |
|------|------|----------|
| L4 | `Option.jl` | 配置 |
| L5 | `SetMesh.jl` | 网格 |
| L6 | `SetParams.jl` | 参数 |
| L7 | `CouplingState.jl` | 类型定义（MultiSPMeLayout + MeshGeometry） |
| L8 | `SetCase.jl` | 案例装配 |
| L9-L12 | `czm.jl`、`CzmUnitMesh.jl`、`CzmSolve.jl`、`CzmPostProcess.jl` | CZM 内聚力模型 |
| L13 | `Assemble.jl` | 装配 |
| L14-L17 | `ElectrodeDiffusion.jl`、`ElectrolyteDiffusion.jl`、`ElectrodePotential.jl`、`ElectrolytePotential.jl` | 电化学子系统 |
| L18-L20 | `SPM.jl`、`SPMe.jl`、`P2D.jl` | 模型实现 |
| L21-L22 | `Parallelsolution.jl`、`Tools.jl` | 分流求解/工具 |
| L23-L24 | `CallModel.jl`、`Solve.jl` | 模型调度与主求解器 |
| L25-L26 | `PostProcessing.jl`、`Materialmatrix.jl` | 后处理/材料矩阵 |
| L27-L29 | `Thermal.jl`、`ThermalDistributed.jl`、`ThermalPolar2D.jl` | 热模型 |
| L30-L31 | `Variables.jl`、`Initialisation.jl` | 状态/初始化 |
| L32-L34 | `mechanical.jl`、`Jellyrollmodel.jl`、`ring.jl` | 力学/Jellyroll 几何/圆环 |
| L35-L37 | `CycleSolver.jl`、`CycleData.jl`、`CsvExport.jl` | 循环求解/数据导出 |

### export 分组（按出现顺序）

| 行号 | 分组 | 内容 |
|------|------|------|
| L40 | 基础 | `Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell` |
| L41 | 网格/参数 | `Mesh1D, Mesh2D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables` |
| L42-L44 | 模型/初始化 | `SPM, Solve, SPMe, SPMe_element, ModelInitialisation, ModelInitialisation_MultiSPMe, extract_element_state, get_thermal_dofs, update_state` |
| L45-L46 | 耦合类型 | `MultiSPMeLayout, MeshGeometry, CohesiveElementGeom, CZMAssemblyCache, CZMAssemblyWorkspace, CzmInterfaceParams, CzmParamCache` |
| L47-L49 | Jellyroll 几何 | `Arrhenius, IntV, IntQ4, jellyroll_collector_seed_mesh, jellyroll_element_properties, jellyroll_tab_node_indices, edge_boundary, jellyroll_element_centers` |
| L50-L53 | 网格/温度工具 | `ring_mesh, setup_thermal2D_mesh, thermal2D_volume_average_temperature` |
| L54-L59 | 热模型 | `ThermalDistributed2D, ThermalDistributed2D_BC, ThermalDistributed2D_Ring, ThermalRing2D_BC, ThermalPolar2D_Ring, identify_boundary_nodes, apply_convection_bc, apply_cool_method, compute_heat_sources, compute_heat_sources_with_czm, solve_branch_currents, ThermalModel, ThermalLumpedModel, ThermalDistributed2DModel` |
| L60-L62 | 力学/CZM 装配 | `thermal_diffusion_stress_2D, compute_czm_params_per_interface, build_thermal_to_czm_interp` |
| L64-L78 | CZM 本构/装配/查询 | `CohesiveElement, CohesiveMesh, DamageState, CZMResult, CzmSubmesh, create_czm_mesh, create_unit_czm_strip, compute_separation, bilinear_traction, bilinear_tangent, update_damage, assemble_czm_system, assemble_coupled_system, assemble_bulk_stiffness, assemble_thermal_chemical_load, assemble_coupled_system_full, apply_bc_czm, identify_bc_nodes_czm, newton_raphson_czm, solve_czm_step, solve_czm_basic_step, solve_czm_arc_length_step, build_czm_cache, ensure_czm_cache, get_damage_statistics, check_fracture_criterion, reset_damage_states, accumulate_cycle_damage, czm_output_to_variables, compute_gap_conductance, compute_element_gap_conductance, get_fractured_elements, get_active_elements, compute_all_gap_conductances, effective_area_factor, map_czm_damage_to_thermal` |
| L80-L83 | 循环求解器 | `CycleOption, PhaseType, PHASE_CHARGE, PHASE_REST, PHASE_DISCHARGE, PhaseResult, CycleResult, CyclingResult, solve_phase, solve_cycling, plot_cycling_results, compute_cs0_from_soc, apply_initial_soc!` |
| L85-L87 | 数据导出 | `TimeStepData, CycleExportData, solve_phase_with_export, solve_cycle_with_export, export_cycle_data_to_csv, load_cycle_data_from_csv` |
| L89 | CSV 导出 | `CZMSnapshot, export_cycling_csv, CsvExportOptions` |

## 省略项

- `using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator, Statistics, Printf`（L2）：标准依赖加载
- `module JuBat ... end`（L1, L90）：模块包装

### [DEBUG]

无

### [PLACEHOLDER]

无

### [COMPLEX-CHECK]

无
