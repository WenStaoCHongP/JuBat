# CZM 子模块独立化设计

> 日期: 2026-05-19
> 状态: 设计评审 (v3 — 已修复 v2 审查问题)

## 1. 目标

将 CZM（内聚力模型）代码从 JuBat 耦合层中提取为独立的 Julia submodule `JuBat.CZM`，放到 `src/CZM/` 目录下。

**设计约束：**
- CZM 模块可独立运行内聚力仿真（不依赖电化学/热模型）
- CZM 模块可通过纯数据参数接口接入 JuBat 耦合模型
- CZM 模块可依赖 JuBat 基础类型（`Mesh`, `Params`, `Scale`, `Cohesive`），但不依赖耦合层类型（`Case`, `Option`, `variables Dict`, `MultiSPMeLayout`）

## 2. 当前状态

### 2.1 CZM 代码分布（4 个文件）

| 文件 | 行数 | 内容 |
|------|------|------|
| `src/czm.jl` | ~662 | 结构体 + 网格构建 + 系统组装 + 缓存 + BC |
| `src/CzmSolve.jl` | ~663 | CZMResult + 3 种 Newton 求解器 |
| `src/CzmPostProcess.jl` | ~118 | 损伤统计/管理/输出转换 |
| `src/Materialmatrix.jl` (CZM 部分) | ~344 | 双线性本构 + 切线刚度 + 间隙热导 |

### 2.2 CZM 与耦合层的交互点

**CZM → 耦合层（CZM 输出的数据）：**
- 损伤状态 D → 间隙热导系数 `compute_gap_conductance` → 热模型 BC
- 断裂单元列表 `get_fractured_elements` → 热源屏蔽 `compute_heat_sources_with_czm`
- 活跃单元列表 `get_active_elements` → 分流求解 `solve_branch_currents`
- 损伤统计 `get_damage_statistics` → 输出和诊断
- 位移场 `result.displacement` → 跨时间步保持 `case.czm_layout.u_prev`

**耦合层 → CZM（CZM 需要的输入）：**
- 有效材料参数 E_eff, ν_eff, α_eff, β_n, β_p（从 param 计算）
- 单元级应变输入 dT_elem, Δsoc_n_elem, Δsoc_p_elem（从温度场和 SOC 分布提取）
- CZM 求解选项 max_iter, tol, iter_method, visc_beta（从 opt 读取）
- 上一步位移 u_prev（从 czm_layout 获取）

## 3. 目标结构

### 3.1 文件组织

```
src/CZM/
├── CZM.jl              # module CZM ... end，include 子文件
├── Types.jl             # 所有 CZM 类型定义
├── Mesh.jl              # 网格构建
├── Constitutive.jl      # 本构模型 + 间隙热导
├── Assembly.jl          # 系统组装 + 缓存
├── Boundary.jl          # 边界条件
├── Solve.jl             # 求解器
└── PostProcess.jl       # 后处理 + 统计
```

### 3.2 各文件详细内容

#### `src/CZM/CZM.jl` — 模块入口

```julia
module CZM

using ..JuBat: Mesh, GaussPoint, NCweight, IntQ4
using ..JuBat: Params, Scale, Cohesive
using ..JuBat: CohesiveMesh, AbstractCohesiveElement, AbstractDamageState
using ..JuBat: identify_boundary_nodes, Assemble, Assemble1D
using LinearAlgebra, SparseArrays, Statistics

include("Types.jl")
include("Constitutive.jl")
include("Mesh.jl")
include("Assembly.jl")
include("Boundary.jl")
include("Solve.jl")
include("PostProcess.jl")

# 导出列表
export CohesiveElement, DamageState, CohesiveMesh, CZMResult
export CohesiveElementGeom, CZMAssemblyCache, CZMAssemblyWorkspace, CzmLayout
export create_czm_mesh
export bilinear_traction_state, bilinear_traction, bilinear_tangent
export update_damage, compute_gap_conductance
export assemble_czm_system, assemble_bulk_stiffness
export assemble_thermal_chemical_load, assemble_coupled_system, assemble_coupled_system_full
export build_czm_cache, ensure_czm_cache
export apply_bc_czm, identify_bc_nodes_czm
export solve_czm_step, solve_czm_basic_step, newton_raphson_czm, solve_czm_arc_length_step
export get_damage_statistics, check_fracture_criterion
export reset_damage_states, accumulate_cycle_damage
export czm_output_to_variables
export get_fractured_elements
export compute_all_gap_conductances, compute_element_gap_conductance
export clone_damage_states, clone_czm_mesh_with_damage

end # module CZM
```

#### `src/CZM/Types.jl` — 类型定义

从当前代码中提取以下类型（不修改逻辑）：

- `CohesiveElement` — 从 `czm.jl:1`
- `DamageState` — 从 `czm.jl:24`
- `CZMResult` — 从 `CzmSolve.jl:1`
- `CohesiveElementGeom` — 从 `CouplingState.jl:106`
- `CZMAssemblyWorkspace` — 从 `CouplingState.jl:124`
- `CZMAssemblyCache` — 从 `CouplingState.jl:170`
- `CzmLayout` — 从 `CouplingState.jl:193`

**保留在 `SetMesh.jl` 的类型（不移动）：**
- `AbstractCohesiveElement` — `SetMesh.jl:23`
- `AbstractDamageState` — `SetMesh.jl:24`
- `CohesiveMesh` (mutable struct + 默认构造器) — `SetMesh.jl:26-46`

**原因：** `Tools.jl:identify_boundary_nodes` 使用 `mesh isa CohesiveMesh` 进行类型分派。由于 `Tools.jl` 必须在 CZM 模块之前加载（CZM 需要 `IntQ4` 和 `identify_boundary_nodes`），而 `CohesiveMesh` 在 `SetMesh.jl` 中定义（`Tools.jl` 之前加载），因此这三个类型必须保留在 `SetMesh.jl`。CZM 模块通过 `using ..JuBat: CohesiveMesh, AbstractCohesiveElement, AbstractDamageState` 引用它们。

**`SetMesh.jl` 无变更：** 上述 3 个类型定义保留在原位，不删除。

#### `src/CZM/Constitutive.jl` — 本构模型

从 `Materialmatrix.jl` 中移入 CZM 相关函数：
- `bilinear_traction_state` (Materialmatrix.jl:68)
- `bilinear_traction` (Materialmatrix.jl:156)
- `bilinear_tangent` (Materialmatrix.jl:175)
- `update_damage` (Materialmatrix.jl:282)
- `compute_gap_conductance` (Materialmatrix.jl:313)
- `compute_element_gap_conductance` (Materialmatrix.jl:348)

**留在 Materialmatrix.jl 的函数：**
- `thermal_capacity_weights_2d` — 纯热学函数，不属于 CZM
- `thermal_anisotropic_conductivity_2d` — 纯热学函数，不属于 CZM

#### `src/CZM/Mesh.jl` — 网格构建

从 `czm.jl` 中移入：
- `create_czm_mesh(thermal_mesh::Mesh, param_dim; tol)` (czm.jl:56)

**依赖：** `Mesh`, `GetGS`, `NCweight`（从 JuBat 顶层 import）

#### `src/CZM/Assembly.jl` — 系统组装

从 `czm.jl` 中移入：
- `assemble_czm_system` (czm.jl:156)
- `assemble_bulk_stiffness` (czm.jl:331)
- `assemble_thermal_chemical_load` (czm.jl:397)
- `build_czm_cache` (czm.jl:463)
- `assemble_coupled_system` (czm.jl:562)
- `assemble_coupled_system_full` (czm.jl:585)

**注意：** `ensure_czm_cache` (czm.jl:551) 当前接受 `case::Case` 参数，违反"CZM 不依赖 Case"的约束。将重写为不接受 Case 的版本：

```julia
# 在 CZM 模块中（纯数据接口）
function ensure_czm_cache(cache, czm_mesh, E_eff, ν_eff, param)
    if cache === nothing || !cache.valid ||
       cache.E_eff != E_eff || cache.ν_eff != ν_eff ||
       length(cache.cohesive_geom) != czm_mesh.n_cohesive
        cache = build_czm_cache(czm_mesh, E_eff, ν_eff, param)
    end
    return cache
end
```

耦合层的适配器 `update_czm_damage!` 负责传递 `case.czm_cache` 和 `case.param`：
```julia
cache = CZM.ensure_czm_cache(case.czm_cache, case.czm_mesh, E_eff, ν_eff, case.param)
case.czm_cache = cache  # 写回（如果重建了）
```

**依赖：** `IntQ4`, `Assemble`, `Assemble1D`

#### `src/CZM/Boundary.jl` — 边界条件

从 `czm.jl` 中移入：
- `apply_bc_czm` (czm.jl:608)
- `identify_bc_nodes_czm` (czm.jl:642)

从 `CzmSolve.jl` 中移入：
- `apply_czm_dirichlet!` (CzmSolve.jl:116)
- `zero_czm_bc_entries!` (CzmSolve.jl:123)
- `extract_bc_dofs` (CzmSolve.jl:52)

**依赖：** `identify_boundary_nodes`（从 JuBat 顶层 import）

#### `src/CZM/Solve.jl` — 求解器

从 `CzmSolve.jl` 中移入：
- `clone_damage_states` (CzmSolve.jl:17)
- `clone_czm_mesh_with_damage` (CzmSolve.jl:31)
- `fill_czm_result!` (CzmSolve.jl:130)
- `build_arc_length_augmented_matrix` (CzmSolve.jl:142)
- `backtrack_line_search!` (CzmSolve.jl:83)
- `solve_czm_basic_step` (CzmSolve.jl:154)
- `solve_czm_arc_length_step` (CzmSolve.jl:261)
- `newton_raphson_czm` (CzmSolve.jl:473)
- `solve_czm_step` (CzmSolve.jl:639)

#### `src/CZM/PostProcess.jl` — 后处理

从 `CzmPostProcess.jl` 中移入全部：
- `get_damage_statistics` (CzmPostProcess.jl:12)
- `check_fracture_criterion` (CzmPostProcess.jl:37)
- `reset_damage_states` (CzmPostProcess.jl:59)
- `accumulate_cycle_damage` (CzmPostProcess.jl:70)
- `czm_output_to_variables` (CzmPostProcess.jl:99)

从 `Materialmatrix.jl` 中移入：
- `get_fractured_elements` (Materialmatrix.jl:358)
- `compute_all_gap_conductances` (Materialmatrix.jl:394)

**保留在 `Materialmatrix.jl` 的函数：**
- `get_active_elements` (Materialmatrix.jl:371) — 依赖 `MeshGeometry` 类型（定义在 `CouplingState.jl`），属于耦合层逻辑（根据断裂单元过滤热源/分流），不属于 CZM 内部算法。外部调用点直接通过 `JuBat.get_active_elements` 调用，无需 CZM 前缀。

### 3.3 CouplingState.jl 的变更

`CouplingState.jl` 中与 CZM 相关的代码需要重写为**薄适配层**：

#### 保留在 CouplingState.jl 中的函数（重写为适配器）

```julia
# compute_czm_effective_params — 不变，仍从 case.param 计算
function compute_czm_effective_params(case)
    # 同现有逻辑
end

# compute_czm_strain_inputs — 不变，仍从 variables 提取
function compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)
    # 同现有逻辑
end

# update_czm_damage! — 重写为薄适配层
function update_czm_damage!(case, variables, T_nodes_carry)
    E_eff, ν_eff, α_eff, β_n, β_p = compute_czm_effective_params(case)
    cache = CZM.ensure_czm_cache(case.czm_cache, case.czm_mesh, E_eff, ν_eff, case.param)
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = compute_czm_strain_inputs(...)

    # 计算粘性参数
    visc_beta = ...
    # 调用 CZM 模块
    result, updated_czm_mesh = CZM.solve_czm_step(
        case.czm_mesh, F_ext, E_eff, ν_eff, case.param.cohesive, case.param,
        case.czm_layout.u_prev;
        α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem,
        max_iter=case.opt.czm_max_iter, tol=case.opt.czm_tol,
        iter_method=case.opt.czm_iter_method,
        n_load_steps=case.opt.czm_load_steps,
        arc_length_alpha=case.opt.czm_arc_length_alpha,
        cache=cache, visc_beta=visc_beta
    )
    # 提交损伤状态（同现有逻辑）
    if result.converged
        case.czm_mesh.damage_states = updated_czm_mesh.damage_states
        case.czm_layout.u_prev = result.displacement
    end
    return result.displacement, result.converged
end
```

#### 从 CouplingState.jl 移除的类型定义

以下类型定义移入 `src/CZM/Types.jl`：
- `CohesiveElementGeom`
- `CZMAssemblyWorkspace`
- `CZMAssemblyCache`
- `CzmLayout`

`CouplingState.jl` 保留：
- `MultiSPMeLayout` 及其构造器
- `BoundaryEdgeCache` 及其 `compute_boundary_edge_cache`
- `MeshGeometry`
- `compute_czm_effective_params`
- `compute_czm_strain_inputs`
- `update_czm_damage!` (适配器)

#### Case 结构体中 CZM 字段的类型

```julia
mutable struct Case
    # ...
    czm_mesh::Union{Nothing, CZM.CohesiveMesh}      # 改为 CZM.CohesiveMesh
    czm_cache::Union{Nothing, CZM.CZMAssemblyCache}  # 改为 CZM.CZMAssemblyCache
    czm_layout::Union{Nothing, CZM.CzmLayout}        # 改为 CZM.CzmLayout
end
```

### 3.4 JuBat.jl 的变更

```julia
module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator, Statistics, Printf

include("Option.jl")
include("SetMesh.jl")
include("SetParams.jl")
include("Tools.jl")            # ← 提前：CZM 需要 IntQ4, identify_boundary_nodes
include("CouplingState.jl")   # MultiSPMeLayout + MeshGeometry + CZM适配器
include("SetCase.jl")
include("CZM/CZM.jl")         # ← 新增：CZM 子模块（依赖 Mesh, Params, Tools）
include("Assemble.jl")
include("ElectrodeDiffusion.jl")
include("ElectrolyteDiffusion.jl")
include("ElectrodePotential.jl")
include("ElectrolytePotential.jl")
include("SPM.jl")
include("SPMe.jl")
include("P2D.jl")
include("Parallelsolution.jl")
include("CallModel.jl")
include("Solve.jl")
include("PostProcessing.jl")
include("Materialmatrix.jl")   # 仅保留热学函数
include("Thermal.jl")
include("ThermalDistributed.jl")
include("ThermalPolar2D.jl")
include("Variables.jl")
include("Initialisation.jl")
include("mechanical.jl")
include("Jellyrollmodel.jl")
include("ring.jl")
include("CycleSolver.jl")
include("CycleData.jl")

# CZM 导出 — 通过顶层 re-export 保持现有接口不变
using .CZM
export CohesiveElement, CohesiveMesh, DamageState, CZMResult
export create_czm_mesh, compute_separation
export bilinear_traction, bilinear_tangent, update_damage
# ... (保持现有 export 列表)
end
```

**include 顺序说明：**
- `Tools.jl` 从原位置（Parallelsolution 之后）提前到 `SetParams.jl` 之后，因为 CZM 需要 `IntQ4` 和 `identify_boundary_nodes`
- `Tools.jl` 本身不依赖 CZM 之后的任何文件，所以提前是安全的
- `compute_separation` 留在 `Tools.jl` 中不变（仅做单元分离计算，是通用工具函数）

### 3.5 Materialmatrix.jl 的变更

移除 CZM 相关函数（bilinear_traction_state, bilinear_traction, bilinear_tangent, update_damage, compute_gap_conductance 等）。

保留纯热学函数：
- `thermal_capacity_weights_2d`
- `thermal_anisotropic_conductivity_2d`

### 3.6 外部调用点的更新

| 调用位置 | 当前调用 | 变更后调用 |
|----------|---------|-----------|
| `CallModel.jl:109-123` | `compute_heat_sources_with_czm` | 内部调用 `CZM.get_damage_statistics` |
| `CallModel.jl:111` | `compute_heat_sources_with_czm(case, ...)` | 内部调用 `get_active_elements`（留在 Materialmatrix.jl） |
| `Solve.jl:206` | `CzmLayout(case.czm_mesh)` | `CZM.CzmLayout(case.czm_mesh)` |
| `Solve.jl:278` | `update_czm_damage!(case, variables, ...)` | 不变（适配器） |
| `ThermalDistributed.jl:292-310` | `compute_gap_conductance(D, δ_n, ...)` | `CZM.compute_gap_conductance(...)` |
| `CallModel.jl:61-70` | `variables["deactivated_elements"]` | 内部调用 `CZM.get_fractured_elements` |

## 4. 不变的文件

以下文件不需要修改（或仅做类型前缀微调）：
- `Option.jl` — 不变
- `SetParams.jl` — 不变
- `SPMe.jl` — 不变
- `Parallelsolution.jl` — 不变
- `mechanical.jl` — 不变
- `Jellyrollmodel.jl` — 不变
- `CycleSolver.jl` — 仅更新 CZM 类型前缀
- `Variables.jl` — 不变
- `Initialisation.jl` — 仅更新 CZM 类型前缀
- `PostProcessing.jl` — 不变

## 5. 实施风险与缓解

### 5.1 include 顺序敏感性

CZM 模块依赖 `Mesh`, `Params`, `CouplingState` 中的 `MeshGeometry`。需确保在 `CouplingState.jl` 之后加载。

**缓解：** 在 `CZM/CZM.jl` 中使用 `using ..JuBat` 引用已加载的类型。

### 5.2 循环依赖

CZM 的 `get_active_elements` 接受 `MeshGeometry` 参数。`MeshGeometry` 定义在 `CouplingState.jl`。

**缓解：** `get_active_elements` 保留在 `Materialmatrix.jl`（耦合层），不进入 CZM 模块。它本质上是耦合逻辑（根据损伤状态过滤活跃单元），不是 CZM 内部算法。`MeshGeometry` 保持在 JuBat 顶层（`CouplingState.jl`）。

### 5.3 `identify_boundary_nodes` 与 `CohesiveMesh` 依赖

CZM 的 `identify_bc_nodes_czm` 调用 `identify_boundary_nodes`（定义在 `Tools.jl`）。`identify_boundary_nodes` 内部使用 `mesh isa CohesiveMesh` 进行类型分派。

**缓解：**
- `CohesiveMesh` 保留在 `SetMesh.jl`（在 `Tools.jl` 之前加载），不移动到 CZM
- `Tools.jl` 在 `SetParams.jl` 之后、`CZM/CZM.jl` 之前加载
- CZM 模块通过 `using ..JuBat` 引用 `CohesiveMesh` 和 `identify_boundary_nodes`
- 加载顺序：`SetMesh.jl` (定义 CohesiveMesh) → `Tools.jl` (使用 CohesiveMesh) → `CZM/CZM.jl` (使用两者)

## 6. 验证标准

1. **独立运行**：CZM 模块可独立创建网格、组装系统、求解、后处理，无需 Case/variables
2. **耦合回归**：现有 SPMe+热+CZM 全耦合仿真的输出数值不变
3. **导出兼容**：JuBat 顶层 re-export 所有 CZM 函数，外部代码无需修改
4. **文件清晰**：每个 CZM 文件 < 400 行，职责单一
