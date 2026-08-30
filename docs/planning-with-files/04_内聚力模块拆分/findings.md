# Findings & Decisions

## Requirements
- 用户希望制定内聚力模块的拆分和简化计划，而不是先做大规模代码改动。
- 计划需要覆盖模块边界、迁移顺序、兼容性和验证步骤。
- 目标是降低后续维护成本，减少单文件臃肿和重复逻辑。

## Quantitative Baseline

### 文件规模
| 文件 | 行数 | 函数数 | 类型数 |
|------|------|--------|--------|
| `src/CzmSolve.jl` | 1001 | 19 | 1 (`CZMResult`) |
| `src/czm.jl` | 659 | 11 | 2 (`CohesiveElement`, `DamageState`) |
| `src/Materialmatrix.jl` (CZM 部分) | ~120 | 3 | 0 |
| **合计** | **1660+** | **30+** | **3** |

### CzmSolve.jl 函数清单 (1001 行)
| # | 函数 | 行号 | 行数 | 职责 | 调用者 | 最大嵌套 |
|---|------|------|------|------|--------|----------|
| 1 | `CZMResult(ndof, n_coh)` | 12-15 | 4 | 结果结构体构造 | 内部 | 0 |
| 2 | `clone_damage_states` | 17-28 | 12 | 深拷贝损伤状态 | solve_czm_basic_step, solve_czm_arc_length_step, newton_raphson_czm | 1 |
| 3 | `clone_czm_mesh_with_damage` | 30-43 | 14 | 用新损伤状态重建 czm_mesh | solve_czm_basic_step, solve_czm_arc_length_step | 1 |
| 4 | `apply_czm_dirichlet!` | 45-50 | 6 | 对位移向量施加 Dirichlet BC | 所有求解器 | 1 |
| 5 | `zero_czm_bc_entries!` | 52-57 | 6 | 对向量零化 BC 自由度 | arc_length 求解器 | 1 |
| 6 | `fill_czm_result!` | 59-69 | 11 | 填充 CZMResult | 所有求解器 | 1 |
| 7 | `build_arc_length_augmented_matrix` | 71-81 | 11 | 构造弧长法增广矩阵 | arc_length 求解器 | 1 |
| 8 | `solve_czm_basic_step` | 83-219 | **137** | Basic Newton 求解（单步） | solve_czm_step | **4** |
| 9 | `solve_czm_arc_length_step` | 221-440 | **220** | 弧长法求解 | solve_czm_step | **5** |
| 10 | `newton_raphson_czm` | 455-647 | **193** | Newton-Raphson + 载荷子步 | solve_czm_step | **4** |
| 11 | `solve_czm_step` | 654-675 | 22 | 方法调度器 | update_czm_damage! | 1 |
| 12 | `get_damage_statistics` | 684-704 | 21 | 损伤统计 | PostProcessing, CallModel, czm_output_to_variables, 外部脚本 | 1 |
| 13 | `check_fracture_criterion` | 709-724 | 16 | 断裂准则判定 | 外部脚本 | 1 |
| 14 | `reset_damage_states` | 731-747 | 17 | 重置损伤 | 外部脚本 | 1 |
| 15 | `accumulate_cycle_damage` | 754-789 | 36 | 累积循环损伤 | 外部脚本 | 2 |
| 16 | `czm_output_to_variables` | 798-816 | 19 | 后处理：结果写入 variables | 外部 | 1 |
| 17 | `compute_czm_effective_params` | 835-855 | 21 | 计算有效材料参数 | update_czm_damage!, 外部脚本 | 1 |
| 18 | `compute_czm_strain_inputs` | 867-915 | 49 | 计算单元级应变输入 | update_czm_damage!, 外部脚本 | 2 |
| 19 | `update_czm_damage!` | 936-1001 | **66** | CZM 损伤更新入口（被 Solve.jl 调用） | Solve.jl:276 | 1 |

### czm.jl 函数清单 (659 行)
| # | 函数 | 行号 | 行数 | 职责 | 调用者 | 最大嵌套 |
|---|------|------|------|------|--------|----------|
| 1 | `DamageState()` | 32 | 1 | 默认构造 | 多处 | 0 |
| 2 | `create_czm_mesh` | 54-134 | 81 | 从热网格创建 CZM 网格 | 外部 (SetCase) | 2 |
| 3 | `assemble_czm_system` | 154-318 | **165** | 组装内聚力全局刚度和内力 | assemble_coupled_system | 4 |
| 4 | `assemble_bulk_stiffness` | 329-393 | 65 | 组装 Q4 固体刚度 | build_czm_cache, assemble_coupled_system | 3 |
| 5 | `assemble_thermal_chemical_load` | 395-449 | 55 | 热-化学载荷向量 | assemble_coupled_system_full, 求解器 | 2 |
| 6 | `build_czm_cache` | 461-542 | 82 | 构建装配缓存 | ensure_czm_cache | 2 |
| 7 | `ensure_czm_cache` | 549-558 | 10 | 缓存有效性检查 | update_czm_damage! | 1 |
| 8 | `assemble_coupled_system` | 560-581 | 22 | 组装耦合系统 (K_bulk + K_coh) | 求解器 | 1 |
| 9 | `assemble_coupled_system_full` | 583-599 | 17 | 含热化学载荷的完整系统 | (当前未使用) | 1 |
| 10 | `apply_bc_czm` | 606-638 | 33 | 惩罚法施加 BC | 求解器 | 2 |
| 11 | `identify_bc_nodes_czm` | 640-659 | 20 | 识别边界节点 | build_czm_cache, 求解器 | 2 |

### 依赖的类型（定义在其他文件）
| 类型 | 定义文件 | 行号 |
|------|----------|------|
| `AbstractCohesiveElement` | `src/SetMesh.jl` | 23 |
| `AbstractDamageState` | `src/SetMesh.jl` | 24 |
| `CohesiveMesh` | `src/SetMesh.jl` | 26 |
| `CohesiveElementGeom` | `src/CouplingState.jl` | 106 |
| `CZMAssemblyWorkspace` | `src/CouplingState.jl` | 124 |
| `CZMAssemblyCache` | `src/CouplingState.jl` | 170 |
| `bilinear_traction_state` | `src/Materialmatrix.jl` | 68 |
| `bilinear_tangent` | `src/Materialmatrix.jl` | 168 |
| `update_damage` | `src/Materialmatrix.jl` | 262 |

### include 顺序
```
SetMesh.jl → SetParams.jl → CouplingState.jl → SetCase.jl → czm.jl → CzmSolve.jl → ... → Materialmatrix.jl → ... → Solve.jl → ...
```

### 重复代码统计
| 重复模式 | 出现位置 | 估计行数 |
|----------|----------|----------|
| BC 提取逻辑 (`bc_nodes → bc_dofs/bc_vals`) | solve_czm_basic_step:94-115, solve_czm_arc_length_step:230-251, newton_raphson_czm:466-488 | ~60 行 (3×~20) |
| `clone_czm_mesh` 样板代码 | clone_czm_mesh_with_damage:31-43, newton_raphson_czm:634-645, reset_damage_states:734-744, accumulate_cycle_damage:776-787 | ~44 行 (4×~11) |
| Newton + 线搜索循环 | solve_czm_basic_step:129-193, solve_czm_arc_length_step:317-403, newton_raphson_czm:517-586 | ~180 行 (3×~60) |
| 最终系统重组 + 残差计算 | 三个求解器各自末尾 | ~30 行 (3×~10) |
| **总重复估计** | | **~314 行** |

### 外部调用点
| 调用方文件 | 调用的 CZM 函数 |
|-----------|----------------|
| `src/Solve.jl:276` | `update_czm_damage!` |
| `src/CallModel.jl:114` | `get_damage_statistics` |
| `src/PostProcessing.jl:147,189` | `get_damage_statistics` |
| `src/JuBat.jl:59,63,65,66` | exports |
| `tools/verify_czm_system.jl` | `solve_czm_step` |
| `example/循环验证/czm_from_precomputed_example.jl` | `solve_czm_step`, `get_damage_statistics` |
| `docs/planning-with-files/06_损伤异常/*.jl` | `compute_czm_effective_params`, `compute_czm_strain_inputs`, `solve_czm_step` |

## Research Findings
- `src/CzmSolve.jl` (1001 行) 目前同时承载了损伤状态克隆、边界条件处理、载荷步推进、Newton 求解、弧长法求解以及结果回填，职责明显过宽。
- `src/czm.jl` (659 行) 主要集中在内聚力网格、内聚力单元、损伤状态、本构牵引和组装相关逻辑，已经比求解器更接近"模型层"。
- `src/JuBat.jl` 目前通过先 include `czm.jl` 再 include `CzmSolve.jl` 的顺序把两者串起来，说明拆分时必须显式处理 include 顺序。
- 仓库记忆显示，当前 CZM 路径已经有归一化、后处理键名和收敛性方面的历史坑，拆分计划不能回退到旧的硬限幅或静默跳过模式。
- 仓库记忆还提示，`DamageState` 只有零参数构造器，直接按字段重建的路径容易出错，这类对象复制逻辑需要在拆分时特别留意。
- 三个求解器 (basic, arc_length, newton_raphson) 有约 180 行重复的 Newton + 线搜索逻辑，是最大的重复来源。
- reviewer 指出 `newton_raphson_czm` 的 BC 残差语义是惩罚式 `R[dof] = val - u[dof]`，而 `backtrack_line_search!` 采用零化式 BC 处理；因此该 helper 只能用于 `solve_czm_basic_step`，不能直接共享给 `newton_raphson_czm`。
- `clone_czm_mesh` 样板代码重复 4 次（11 行×4），是第二大的重复来源。
- BC 提取逻辑重复 3 次（20 行×3），有缓存时跳过但逻辑相同。
- `CZMResult` 和 `DamageState` 都只有零参数/双参数构造器，缺少通用构造器，导致复制逻辑必须逐字段赋值。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 保留现有公开接口，优先做兼容层 | 降低对现有例程和调用点的影响 |
| 将求解调度与本构/组装分开 | 让每个文件只负责一类变化原因 |
| 先抽重复 helper，再改文件结构 | 先降低复杂度，再搬迁更稳妥 |
| 将 `backtrack_line_search!` 限定在 `solve_czm_basic_step` | `newton_raphson_czm` 使用惩罚式 BC 残差，不能复用零化式 helper |
| 保持线搜索与失败回滚语义不变 | 这是当前稳定性的关键，不应在重构中丢失 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| 结构边界一开始不够清晰 | 通过阅读 `src/CzmSolve.jl`、`src/czm.jl` 和 `src/JuBat.jl` 先确定职责分层 |
| 可能的 include 顺序风险 | 将 include 顺序作为显式设计约束写入计划 |
| Task 6 初稿把 `newton_raphson_czm` 也纳入 helper 提取 | 按 reviewer 建议收缩范围，保留零化式 BC helper 只服务 basic 路径 |

## Resources
- `src/CzmSolve.jl` (1001 行, 19 函数)
- `src/czm.jl` (659 行, 11 函数)
- `src/Materialmatrix.jl` (375 行, CZM 部分约 120 行)
- `src/CouplingState.jl` (CZMAssemblyWorkspace, CZMAssemblyCache, CohesiveElementGeom)
- `src/SetMesh.jl` (CohesiveMesh, AbstractCohesiveElement, AbstractDamageState)
- `src/JuBat.jl` (include 顺序和 export 声明)
- `src/Solve.jl` (update_czm_damage! 唯一内部调用点)
- `src/CallModel.jl` (get_damage_statistics 调用点)
- `src/PostProcessing.jl` (get_damage_statistics 调用点)

## Call Dependency Graph

```
外部入口 (Solve.jl / CycleSolver / 外部脚本)
    │
    ▼
update_czm_damage!  ◄── Solve.jl:276 (唯一内部调用点)
    ├── compute_czm_effective_params
    ├── ensure_czm_cache ──► build_czm_cache
    │                           ├── assemble_bulk_stiffness
    │                           ├── identify_bc_nodes_czm
    │                           └── CZMAssemblyWorkspace (构造)
    ├── compute_czm_strain_inputs
    └── solve_czm_step (调度器)
        ├── newton_raphson_czm   ◄─── 默认方法 "load_substep"
        ├── solve_czm_basic_step ◄─── "basic"
        └── solve_czm_arc_length_step ◄─── "arc_length"

三个求解器共同调用:
    ├── assemble_coupled_system
    │       ├── assemble_czm_system
    │       │       ├── bilinear_traction_state  [Materialmatrix.jl]
    │       │       └── bilinear_tangent         [Materialmatrix.jl]
    │       └── assemble_bulk_stiffness (via cache)
    ├── assemble_thermal_chemical_load
    ├── apply_bc_czm / apply_czm_dirichlet!
    ├── clone_damage_states
    ├── clone_czm_mesh_with_damage  (basic + arc_length)
    ├── update_damage               [Materialmatrix.jl]
    └── fill_czm_result!

后处理/统计 (独立分支):
    get_damage_statistics ◄── CallModel.jl, PostProcessing.jl, czm_output_to_variables
    check_fracture_criterion
    reset_damage_states
    accumulate_cycle_damage
    czm_output_to_variables
```

### 按职责分层的依赖关系

```
┌───────────────────────────────────────────────────────────┐
│ 调度层 (CzmSolve.jl)                                      │
│   update_czm_damage! → solve_czm_step → {3个求解器}       │
├───────────────────────────────────────────────────────────┤
│ 求解层 (CzmSolve.jl)                                      │
│   solve_czm_basic_step / solve_czm_arc_length_step        │
│   / newton_raphson_czm                                    │
├───────────────────────────────────────────────────────────┤
│ 组装层 (czm.jl)                                           │
│   assemble_czm_system / assemble_bulk_stiffness           │
│   assemble_thermal_chemical_load / assemble_coupled_system│
├───────────────────────────────────────────────────────────┤
│ 本构层 (Materialmatrix.jl)                                │
│   bilinear_traction_state / bilinear_tangent              │
│   update_damage                                           │
├───────────────────────────────────────────────────────────┤
│ 网格/状态层 (czm.jl + SetMesh.jl)                         │
│   create_czm_mesh / CohesiveMesh / DamageState            │
│   CohesiveElement / CohesiveElementGeom                   │
├───────────────────────────────────────────────────────────┤
│ 缓存/工作区层 (czm.jl + CouplingState.jl)                 │
│   CZMAssemblyCache / CZMAssemblyWorkspace                 │
│   build_czm_cache / ensure_czm_cache                      │
├───────────────────────────────────────────────────────────┤
│ 后处理/统计层 (CzmSolve.jl)                               │
│   get_damage_statistics / check_fracture_criterion        │
│   reset_damage_states / accumulate_cycle_damage           │
│   czm_output_to_variables                                 │
└───────────────────────────────────────────────────────────┘
```

## Visual/Browser Findings
- None.
