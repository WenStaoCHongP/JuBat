# CouplingState.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 95
- 路径: `src/CouplingState.jl`
- 来源: 新增（commit ce0cf62），替代 `Dict{String,Any}` 的 `multi_spme_layout`

### 主要结构体列表

|| 结构体 | 行号 | 说明 |
|--------|------|------|
| `MultiSPMeLayout` | L9 | 多SPMe状态向量的布局索引（初始化后不可变） |
| `BoundaryEdgeCache` | L56 | 预计算的外边界边列表（网格不变量） |
| `MeshGeometry` | L81 | Jellyroll网格的几何拓扑信息（构建后不可变） |

### 主要函数列表

|| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)` | L19 | 便捷构造器：自动计算 range 和 n_total，areas 延迟填充 |
| `MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int, mesh_th)` | L26 | 便捷构造器：从网格计算单元面积 |
| `compute_boundary_edge_cache(mesh, is_outer)` | L69 | 预计算外边界边列表和边长（网格不变量） |

## 结构体字段详情

### MultiSPMeLayout
|| 字段 | 类型 | 说明 |
|------|------|------|
| `ne` | Int | 热单元数 |
| `n_chem` | Int | 每单元电化学 DOF 数 = Nrn + Nrp + Nel |
| `nT` | Int | 热节点 DOF 数 |
| `n_total` | Int | 全局状态向量总长 = ne*n_chem + nT |
| `chem_range` | UnitRange{Int} | 1:(ne*n_chem) |
| `thermal_range` | UnitRange{Int} | (ne*n_chem+1):(ne*n_chem+nT) |
| `areas` | Vector{Float64} | 预计算的单元面积（网格不变量） |

### BoundaryEdgeCache
|| 字段 | 类型 | 说明 |
|------|------|------|
| `edges` | Vector{Tuple{Int,Int}} | 外边界边对 (node_a, node_b)，a < b，已去重 |
| `L_edge` | Vector{Float64} | 边长（无量纲），与 edges 一一对应 |

### MeshGeometry
|| 字段 | 类型 | 说明 |
|------|------|------|
| `element_layer` | Vector{Int} | 每个单元的层类型 (1=NE, 2=SP, 3=PE, 4=NCC, 5=PCC) |
| `is_inner_layer` | Vector{Bool} | 是否为内层 |
| `layer_weights` | Matrix{Float64} | ne x 5 层面积权重 [NE, SP, PE, PCC, NCC] |
| `interface_pairs` | Vector{Tuple{Int,Int}} | CZM 界面配对 (top_elem, bot_elem) |
| `czm_element_map` | Dict{Int,Vector{Int}} | 热单元号 → CZM 单元索引向量 |
| `inner_nodes` | Vector{Int} | 内边界节点索引 |
| `outer_nodes` | Vector{Int} | 外边界节点索引 |
| `boundary_edges` | Union{Nothing, BoundaryEdgeCache} | 预计算的边界边缓存 |

## 设计原则

### Fail-Fast 替代 Dict 查询
- 禁止 `haskey`、`get(dict, key, default)`、`try/catch` 吞错误
- 直接访问字段，未初始化就让 Julia 抛异常
- 替代原先 `Dict{String,Any}` 的 `multi_spme_layout` 字典

### 不可变设计
- 使用 `struct`（非 `mutable struct`），字段初始化后不可变
- 避免运行时意外修改布局信息

## 使用位置

|| 文件 | 使用方式 |
|------|----------|
| `SetCase.jl` | Case struct 的 `layout` 和 `geometry` 字段类型 |
| `Initialisation.jl` | `ModelInitialisation_MultiSPMe` 中构造 `MultiSPMeLayout`（含 areas） |
| `CallModel.jl` | `CallModel_MultiSPMe` 中使用 `case.layout.ne/n_chem/nT/areas` |
| `Solve.jl` | 外部状态恢复时构造 `MultiSPMeLayout`（含 areas） |
| `Jellyrollmodel.jl` | `setup_thermal2D_mesh` 中构造 `MeshGeometry`（含 boundary_edges） |
| `ThermalDistributed.jl` | 通过 `case.geometry.layer_weights` 获取层权重，通过 `boundary_edges` 缓存加速对流 BC |
| `CycleData.jl` | 通过 `case.layout` 访问布局信息 |

## 依赖关系

### 该文件依赖哪些其他文件
- `SetMesh.jl` — `Mesh` 类型（`compute_boundary_edge_cache` 和面积构造器中使用 `mesh.gs.detJ`、`mesh.gs.weight`、`mesh.gs.ele`）
- `Tools.jl` — `identify_boundary_nodes`（在 `compute_boundary_edge_cache` 中间接使用）

### 哪些文件依赖该文件
- `SetCase.jl` — Case struct 引用 `MultiSPMeLayout`、`BoundaryEdgeCache` 和 `MeshGeometry` 类型
- `Initialisation.jl` — 构造 `MultiSPMeLayout` 实例
- `CallModel.jl` — 使用 `case.layout` 字段（含 `areas`）
- `Solve.jl` — 构造 `MultiSPMeLayout` 实例
- `Jellyrollmodel.jl` — 构造 `MeshGeometry` 实例（含 `boundary_edges`）
- `ThermalDistributed.jl` — 使用 `edge_cache` 参数加速对流边界条件
- `JuBat.jl` — export `MultiSPMeLayout`, `MeshGeometry`, `BoundaryEdgeCache`

## 耦合分析

此文件是整个耦合框架的**类型基础层**：
- `MultiSPMeLayout` 定义了多SPMe状态向量的结构布局，消除了原先的 Dict 键值访问
- `MeshGeometry` 将 Jellyroll 网格的几何拓扑信息（层权重、界面配对、CZM映射）类型化
- `BoundaryEdgeCache` 预计算外边界边列表，消除热模型每步重复扫描网格的开销
- 三个结构体都是不可变的，确保耦合数据在运行时不被意外修改

## 后续变更 (2026-04-07)

- **P4 `_` 前缀函数重命名**：函数名移除 `_` 前缀
- **布局缓存**: `case.multi_spme_layout[...] = ...` Dict 填充 → `case.layout = MultiSPMeLayout(ne, n_chem, nT)` 构造器

## 后续变更 (2026-04-20)

- **MultiSPMeLayout 新增 `areas` 字段**: 预计算每个热单元的面积（网格不变量），消除 CallModel 中每步重复积分的开销
- **新增带 mesh 构造器**: `MultiSPMeLayout(ne, n_chem, nT, mesh_th)` 从网格高斯积分直接计算 `areas`
- **新增 `BoundaryEdgeCache` 结构体**: 存储预计算的边界边对和边长，用于 `apply_convection_bc` 加速
- **新增 `compute_boundary_edge_cache` 函数**: 从网格和外部节点标记中提取去重的外边界边列表
- **MeshGeometry 新增 `boundary_edges` 字段**: 存储边界边缓存，Jellyrollmodel.jl 构造时预计算
- 行数从约 41 行增加到约 95 行
