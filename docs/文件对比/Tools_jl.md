# Tools.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 18
- 函数列表:
  - `Arrhenius(Eac, T)` - Arrhenius 温度依赖函数
  - `IntV(x, mesh)` - 一维网格上的积分

## Parameters_Design分支
- 行数: 191 (+173, +961%)
- 新增函数列表:
  - `IntQ4(f, x_e, y_e; order)` - Q4 单元积分驱动
  - `identify_boundary_nodes(mesh, param, opt)` - 识别边界节点（内/外）
  - `compute_separation(elem, node, u)` - 计算内聚力单元分离
  - `element_nodal_mean(mesh, nodal_values)` - 节点值到单元均值的映射
  - `thermal2D_volume_average_temperature(mesh, T_nodes)` - 体积加权平均温度
  - `q4_center_gradients(node, elem_nodes)` - Q4 单元中心梯度

## 变更详情

### 原有函数（无变更）
- `Arrhenius(Eac, T)` - 无修改
- `IntV(x, mesh)` - 仅添加末尾 `end`（原文件无换行结尾）

### 新增函数

#### `IntQ4(f, x_e, y_e; order=2)`
- 通用 Q4 单元积分驱动函数
- 使用回调模式: 调用者通过 `do ... end` 块提供积分贡献
- 回调接收: 局部坐标 (xi, eta)、权重 w、空间梯度 (dNdx, dNdy)、detJ
- 用途: 材料矩阵组装、热模型刚度矩阵计算等

#### `identify_boundary_nodes(mesh, param, opt=nothing)`
- 识别网格的内边界和外边界节点
- 支持 `Mesh` 和 `CohesiveMesh` 两种类型
- 使用阿基米德螺旋线参数计算 theta 范围
- 调用 `edge_boundary()`（定义在 `Jellyrollmodel.jl`）
- 返回: `(is_inner, is_outer)` 两个 Bool 向量

#### `compute_separation(elem, node, u)`
- 计算内聚力单元的平均法向/切向分离
- 使用 Newton-Cotes 2 点积分
- 构建旋转矩阵 R 将全局位移转换到局部坐标系
- 返回: `(delta_n, delta_t)` 法向和切向分离

#### `element_nodal_mean(mesh, nodal_values)`
- 将节点场值映射为单元均值
- 对每个单元，取其所有节点值的算术平均
- 用途: 温度场从节点到单元的映射（热-电化学耦合）

#### `thermal2D_volume_average_temperature(mesh, T_nodes)`
- 计算体积加权平均温度
- 使用高斯积分权重和 Jacobian 行列式作为体积权重
- 通过形状函数插值节点温度到积分点
- 回退: 如果 `den <= 0`，使用简单 `mean(T_nodes)`
- 用途: 全局温度监测、集总热模型对比

#### `q4_center_gradients(node, elem_nodes)`
- 计算 Q4 单元中心的形状函数空间梯度
- 在 (xi=0, eta=0) 处求值
- 返回: `(dNdx, dNdy, detJ)` 或 `nothing`（退化单元）
- 用途: 热传导方程的梯度计算

## 依赖关系

### 依赖
- `IntQ4` 依赖 `GSweight`, `LagrangeBasis`（来自 SetMesh.jl）
- `identify_boundary_nodes` 依赖 `edge_boundary`（来自 Jellyrollmodel.jl）
- `compute_separation` 依赖 `NCweight`（来自 SetMesh.jl）
- `thermal2D_volume_average_temperature` 使用 `Statistics.mean`（需要 Statistics 包）
- `q4_center_gradients` 依赖 `LagrangeBasis`（来自 SetMesh.jl）

### 被依赖
- `IntQ4` 被 `Materialmatrix.jl` 使用
- `identify_boundary_nodes` 被热模型边界条件设置使用
- `compute_separation` 被 CZM 求解器使用
- `element_nodal_mean` 被热-电化学耦合使用
- `thermal2D_volume_average_temperature` 被主求解器和循环求解器使用
- `q4_center_gradients` 被热模型使用

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 是

此文件从 2 个简单工具函数扩展为多物理场耦合的核心工具箱：
- `IntQ4` 是 2D 分布式热模型刚度矩阵组装的基础
- `identify_boundary_nodes` 是 Jellyroll 几何热边界条件的基础
- `compute_separation` 是 CZM 损伤演化的基础
- `element_nodal_mean` 是热-电化学耦合（温度传递）的关键桥梁
- `thermal2D_volume_average_temperature` 是全局温度监测的基础
- `q4_center_gradients` 是 2D 热传导方程求解的基础

变更性质: 从"电化学工具集"扩展为"多物理场通用工具集"，新增函数都是为 2D 热模型和 CZM 服务的基础设施。
