# JuBat.jl

- **源文件**: `src/JuBat.jl`
- **行数**: 94 行
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
| L8-L9 | `SetCase.jl`、`StateAccess.jl` | 案例装配与状态访问规则 |
| L10-L15 | `CzmBC.jl`、`Czm.jl`、`CzmMesh.jl`、`CzmUnitMesh.jl`、`CzmSolve.jl`、`CzmPostProcess.jl` | CZM 边界、核心装配、网格、求解与后处理 |
| L16-L20 | `Assemble.jl` 及四个电化学基础方程文件 | 装配与电化学子系统 |
| L21-L23 | `SPM.jl`、`SPMe.jl`、`P2D.jl` | 模型实现 |
| L24-L27 | `Parallelsolution.jl`、`Tools.jl`、`CallModel.jl`、`Solve.jl` | 分流、工具、模型调度与主求解器 |
| L28-L29 | `PostProcessing.jl`、`Materialmatrix.jl` | 后处理/材料矩阵 |
| L30-L32 | `Thermal.jl`、`ThermalDistributed.jl`、`ThermalPolar2D.jl` | 热模型 |
| L33-L34 | `Variables.jl`、`Initialisation.jl` | 状态/初始化 |
| L35-L37 | `mechanical.jl`、`Jellyrollmodel.jl`、`ring.jl` | 力学/Jellyroll 几何/圆环 |
| L38-L41 | `CycleSolver.jl`、`CyclePostProcess.jl`、`CycleData.jl`、`CsvExport.jl` | 循环求解、循环后处理、在线数据采集与 CSV 导出 |

### export 分组（按出现顺序）

| 行号 | 分组 | 内容 |
|------|------|------|
| L44 | 基础 | `Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell` |
| L45 | 网格/参数 | `Mesh1D, Mesh2D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables` |
| L46-L48 | 模型/初始化 | `SPM, Solve, SPMe, SPMe_element, ModelInitialisation, ModelInitialisation_MultiSPMe, extract_element_state, get_thermal_dofs, update_state` |
| L49-L50 | 耦合类型 | `MultiSPMeLayout, MeshGeometry, CohesiveElementGeom, CZMAssemblyCache, CZMAssemblyWorkspace, CzmInterfaceParams, CzmParamCache` |
| L51-L54 | Jellyroll 几何 | `Arrhenius, IntV, IntQ4, jellyroll_collector_seed_mesh, jellyroll_element_properties, jellyroll_tab_node_indices, edge_boundary, jellyroll_element_centers` |
| L55-L57 | 网格/温度工具 | `ring_mesh, setup_thermal2D_mesh, thermal2D_volume_average_temperature` |
| L58-L63 | 热模型 | `ThermalDistributed2D, ThermalDistributed2D_BC, ThermalPolar2D_Ring` 等热模型与边界/热源接口 |
| L64-L66 | 力学/CZM 装配 | `thermal_diffusion_stress_2D, compute_czm_params_per_interface, build_thermal_to_czm_interp` |
| L68-L82 | CZM 本构/装配/查询 | CZM 网格、本构、装配、求解、统计与间隙导热公开 API |
| L84-L87 | 循环求解器 | 循环配置、结果类型、求解、绘图和初始 SOC API |
| L89-L91 | 在线数据导出 | `TimeStepData, CycleExportData, solve_phase_with_export, solve_cycle_with_export, export_cycle_data_to_csv` |
| L93 | CSV 导出 | `CZMSnapshot, export_cycling_csv, CsvExportOptions` |

## 省略项

- `using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator, Statistics, Printf`（L2）：标准依赖加载
- `module JuBat ... end`（L1, L94）：模块包装

### [DEBUG]

无

### [PLACEHOLDER]

无

### [COMPLEX-CHECK]

无
