# CouplingState.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 41
- 路径: `src/CouplingState.jl`
- 来源: 新增（commit ce0cf62），替代 `Dict{String,Any}` 的 `multi_spme_layout`

### 主要结构体列表

| 结构体 | 行号 | 说明 |
|--------|------|------|
| `MultiSPMeLayout` | L9 | 多SPMe状态向量的布局索引（初始化后不可变） |
| `MeshGeometry` | L33 | Jellyroll网格的几何拓扑信息（构建后不可变） |

### 主要函数列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)` | L19 | 便捷构造器：自动计算 range 和 n_total |

## 结构体字段详情

### MultiSPMeLayout
| 字段 | 类型 | 说明 |
|------|------|------|
| `ne` | Int | 热单元数 |
| `n_chem` | Int | 每单元电化学 DOF 数 = Nrn + Nrp + Nel |
| `nT` | Int | 热节点 DOF 数 |
| `n_total` | Int | 全局状态向量总长 = ne*n_chem + nT |
| `chem_range` | UnitRange{Int} | 1:(ne*n_chem) |
| `thermal_range` | UnitRange{Int} | (ne*n_chem+1):(ne*n_chem+nT) |

### MeshGeometry
| 字段 | 类型 | 说明 |
|------|------|------|
| `element_layer` | Vector{Int} | 每个单元的层类型 (1=NE, 2=SP, 3=PE, 4=NCC, 5=PCC) |
| `is_inner_layer` | Vector{Bool} | 是否为内层 |
| `layer_weights` | Matrix{Float64} | ne x 5 层面积权重 [NE, SP, PE, PCC, NCC] |
| `interface_pairs` | Vector{Tuple{Int,Int}} | CZM 界面配对 (top_elem, bot_elem) |
| `czm_element_map` | Dict{Int,Vector{Int}} | 热单元号 → CZM 单元索引向量 |
| `inner_nodes` | Vector{Int} | 内边界节点索引 |
| `outer_nodes` | Vector{Int} | 外边界节点索引 |

## 设计原则

### Fail-Fast 替代 Dict 查询
- 禁止 `haskey`、`get(dict, key, default)`、`try/catch` 吞错误
- 直接访问字段，未初始化就让 Julia 抛异常
- 替代原先 `Dict{String,Any}` 的 `multi_spme_layout` 字典

### 不可变设计
- 使用 `struct`（非 `mutable struct`），字段初始化后不可变
- 避免运行时意外修改布局信息

## 使用位置

| 文件 | 使用方式 |
|------|----------|
| `SetCase.jl` | Case struct 的 `layout` 和 `geometry` 字段类型 |
| `Initialisation.jl` | `ModelInitialisation_MultiSPMe` 中构造 `MultiSPMeLayout` |
| `CallModel.jl` | `CallModel_MultiSPMe` 中使用 `case.layout.ne/n_chem/nT` |
| `Solve.jl` | 外部状态恢复时构造 `MultiSPMeLayout` |
| `Jellyrollmodel.jl` | `setup_thermal2D_mesh` 中构造 `MeshGeometry` |
| `ThermalDistributed.jl` | 通过 `case.geometry.layer_weights` 获取层权重 |
| `CycleData.jl` | 通过 `case.layout` 访问布局信息 |

## 依赖关系

### 该文件依赖哪些其他文件
无外部依赖（纯类型定义和便捷构造器）。

### 哪些文件依赖该文件
- `SetCase.jl` — Case struct 引用 `MultiSPMeLayout` 和 `MeshGeometry` 类型
- `Initialisation.jl` — 构造 `MultiSPMeLayout` 实例
- `CallModel.jl` — 使用 `case.layout` 字段
- `Solve.jl` — 构造 `MultiSPMeLayout` 实例
- `Jellyrollmodel.jl` — 构造 `MeshGeometry` 实例
- `JuBat.jl` — export `MultiSPMeLayout`, `MeshGeometry`

## 耦合分析

此文件是整个耦合框架的**类型基础层**：
- `MultiSPMeLayout` 定义了多SPMe状态向量的结构布局，消除了原先的 Dict 键值访问
- `MeshGeometry` 将 Jellyroll 网格的几何拓扑信息（层权重、界面配对、CZM映射）类型化
- 两个结构体都是不可变的，确保耦合数据在运行时不被意外修改
