# SetMesh.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 494
- 主要结构/函数列表:
  - `mutable struct GaussPoint` - 高斯点数据结构
  - `mutable struct Mesh` - 一维网格结构（L2/L3 单元）
  - `SetMesh()` - 网格构建分发函数（仅支持 L2/L3）
  - `Mesh1D()` - 一维网格构建
  - `PickElement()` - 选取单元子集
  - `CombineMesh()` - 合并多个网格
  - `MultipleMesh()` - 复制/倍增网格
  - `GetGS()` - 高斯积分点计算（仅 L2/L3）
  - `LagrangeBasis()` - 拉格朗日基函数
  - `GSweight()` - 高斯积分权重和点
  - `ShapeFunction1D()` - 一维形状函数及梯度

## Parameters_Design分支
- 行数: 739 (+245, +50%)
- 新增结构/函数列表:
  - `abstract type AbstractCohesiveElement` - CZM 抽象类型
  - `abstract type AbstractDamageState` - 损伤状态抽象类型
  - `mutable struct CohesiveMesh` - 内聚力网格结构（含分层节点、损伤状态等）
  - `Mesh2D()` - 二维 Q4/COH2D4 网格构建
  - `NCweight()` - Newton-Cotes 积分权重
  - `ShapeFunction2D()` - 二维形状函数及梯度（Q4 和 COH2D4）

## 变更详情

### 新增结构体

#### `AbstractCohesiveElement` / `AbstractDamageState`
- 抽象类型，用于多态 CZM 单元和损伤状态
- 定义在 `SetMesh.jl` 而非 `czm.jl`，避免 include 顺序依赖问题

#### `CohesiveMesh` 结构体
- `bulk_mesh::Mesh` - 原始固体网格
- `node::Matrix{Float64}` - 扩展后节点坐标（分层后）
- `nnode::Int64` - 总节点数
- `bulk_element::Matrix{Int64}` - 更新后的固体单元连接
- `cohesive_elements::Vector{AbstractCohesiveElement}` - 内聚力单元
- `n_cohesive::Int64` - 内聚力单元数
- `n_layers::Int64` - 卷绕圈数
- `node_map::Dict{Int64, Vector{Int64}}` - 原节点到分层后节点的映射
- `interface_nodes::Vector{Vector{Tuple{Int64,Int64}}}` - 界面节点对
- `damage_states::Vector{AbstractDamageState}` - 损伤状态向量
- 包含空初始化内部构造函数

### 新增函数

#### `Mesh2D(domain, num, type, gsorder)`
- 构建规则矩形区域的二维 Q4 网格
- 支持 `num = [nx, ny]` 或 `Tuple(nx, ny)`
- 节点编号: 外层 j(y), 内层 i(x)，索引 `j*nnx + i + 1`
- 单元连接: Q4 标准 (1:左下, 2:右下, 3:右上, 4:左上)
- 自动计算高斯积分点

#### `NCweight(order)`
- Newton-Cotes 积分权重和节点
- 支持 order 2-5
- 用于 COH2D4 内聚力单元的积分

#### `ShapeFunction2D(element, type, node, xi, ele_map)`
- Q4 单元: 标准 4 节点双线性形状函数
  - Ni: `(ngs x 4)` 形函数值
  - dNidx: `(ngs x 8)` 梯度 [dN/dx, dN/dy]
  - 通过 Jacobian 变换计算空间梯度
- COH2D4 单元: 退化 4 节点内聚力单元
  - 上表面和下表面中点连线定义局部坐标系
  - 使用法向量 `n_vec` 和切向量 `t_vec` 构建旋转矩阵

### 修改函数

#### `SetMesh()`
- 新增 `type == "Q4"` 和 `type == "COH2D4"` 分支，分发到 `Mesh2D()`

#### `GetGS()`
- 新增 `COH2D4` 分支: 计算内聚力单元的高斯积分
  - 使用 Newton-Cotes 积分 (`NCweight`)
  - 处理退化几何（上下表面中点）
  - 计算 `detJ = L/2`（线段半长）
- 修改通用分支中 `eachindex(w)` 为 `1:size(w,1)` 以兼容二维权重矩阵
- 形状函数调用分发: 根据单元类型选择 `ShapeFunction1D` 或 `ShapeFunction2D`

#### `ShapeFunction1D()`
- 错误信息修正: 从 `mesh.type` 改为 `type`（修正变量名 bug）

## 依赖关系

### 被依赖关系
- `CohesiveMesh` 被 `czm.jl`, `CzmSolve.jl`, `Jellyrollmodel.jl` 使用
- `AbstractCohesiveElement` / `AbstractDamageState` 被 `czm.jl` 中的具体类型继承
- `Mesh2D` 被 `Jellyrollmodel.jl` 中的 `jellyroll_collector_seed_mesh` 调用
- `ShapeFunction2D` 被 `Materialmatrix.jl`, `ThermalDistributed.jl` 等使用
- `NCweight` 被 `czm.jl` 中的 CZM 积分使用
- `Mesh` 结构体被所有模块使用（未修改，仅扩展）

### 依赖关系
- 依赖 `SetParams.jl` 中的 `GaussPoint`（已在此文件定义）
- 无外部依赖

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 是

此文件是整个多物理场耦合框架的基础设施层：
- `CohesiveMesh` + 抽象类型是 CZM 模块的骨架
- `Mesh2D` + `ShapeFunction2D` 是分布式热模型的基础
- `NCweight` 是 CZM 积分的基础
- `COH2D4` 单元支持是内聚力建模的关键

变更性质: 基础设施扩展，从纯 1D 网格扩展到 2D + 内聚力网格，为多物理场耦合提供几何/离散基础。
