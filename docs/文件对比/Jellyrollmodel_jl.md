# Jellyrollmodel.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 562
- 路径: `src/Jellyrollmodel.jl`

### 主要结构体

| 结构体 | 说明 |
|--------|------|
| `JellyrollMesh` | Jellyroll热网格数据结构，包含热网格、合并网格、界面节点对、CZM映射、层信息、极耳节点等 |

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `jellyroll_collector_seed_mesh(param; nθ, gsorder, phase, tol)` | L17 | 基于阿基米德螺旋线生成collector-seeded Q4热网格，识别界面节点对、计算层信息、建立CZM映射 |
| `jellyroll_element_properties(mesh, param)` | L218 | 计算各热单元面积及层面积权重矩阵(ne x 5) [NE, SP, PE, PCC, NCC] |
| `edge_boundary(mesh, nidx, param; which, theta_range, tol)` | L339 | 基于螺旋方程精确识别内/外边界节点 |
| `jellyroll_element_centers(mesh)` | L372 | 计算每个Q4单元的几何中心(ne x 2) |
| `jellyroll_tab_node_indices(mesh, param)` | L389 | 通过弧长二分搜索识别受正/负极耳影响的节点索引 |
| `setup_thermal2D_mesh(case, mesh_data; use_merged)` | L527 | 将JellyrollMesh装入case，自动选择合并/未合并网格，保存界面信息到 `case.geometry` (`MeshGeometry`) |

## 功能描述

本文件实现了Jellyroll（果冻卷）电池的二维螺旋几何建模，是multi-SPMe + distributed2D热模型的核心网格生成模块。主要功能包括：

1. **网格生成**：基于阿基米德螺旋线方程 `r(theta) = a + b*theta`，在内外两条螺旋线之间生成Q4四边形单元。通过等角度采样控制周向分辨率。

2. **界面识别**：检测外螺旋节点与内螺旋节点在下一圈的坐标重合（`interface_pairs`），这些重合点就是相邻卷绕层之间的界面。

3. **层权重计算**：对于每个Q4单元，基于8层结构（PE -> PCC -> PE -> SP -> NE -> NCC -> NE -> SP），计算各材料层的扇形面积权重，用于后续分层热源和等效材料参数计算。

4. **CZM映射**：建立热单元到CZM内聚力单元的映射关系（`czm_element_map`），用于损伤状态对热模型的影响（如界面热阻变化）。

5. **极耳识别**：通过弧长积分 `F(u) = (u*sqrt(u^2+b^2) + b^2*asinh(u/b)) / (2b)` 的二分搜索，精确定位极耳在螺旋线上的角度范围及对应节点。

6. **网格合并**：提供节点合并功能，将重合的界面节点合并为同一节点，生成连续热传导路径的网格。当CZM未启用时使用合并网格。

## 依赖关系

### 该文件依赖
- `src/SetMesh.jl` — `Mesh`结构体、`GetGS`高斯积分函数
- `src/Variables.jl` — `Case`类型定义（`setup_thermal2D_mesh`的参数）
- `src/SetParams.jl` — `Params`结构体（包含PE/NE/SP/PCC/NCC/cell/tab参数）

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("Jellyrollmodel.jl")`（L29）
- `example/testexample.jl` — 调用 `jellyroll_collector_seed_mesh` 和 `setup_thermal2D_mesh`
- `src/Solve.jl` — 通过 `case.mesh["thermal2D"]` 使用生成的网格
- `src/Parallelsolution.jl` — 通过 `case.geometry` 使用界面和层信息

## 耦合分析

本文件是 **multi-SPMe + distributed2D + CZM** 耦合体系的几何基础设施：

- **与multi-SPMe耦合**：`setup_thermal2D_mesh`将界面节点对、层权重、CZM映射等存入 `case.geometry`（`MeshGeometry` struct），供多SPMe求解器使用。每个热单元对应一个独立SPMe模型。

- **与distributed2D热模型耦合**：生成的Q4网格直接作为 `case.mesh["thermal2D"]` 使用。`jellyroll_element_properties` 提供的层权重用于计算各向异性导热系数和分层体积热容。

- **与CZM耦合**：`czm_element_map` 建立热单元到内聚力单元的映射，使损伤状态能影响界面热阻。`interface_pairs` 用于CZM网格创建（`create_czm_mesh`）。当CZM启用时，使用未合并网格保留界面自由度。

- **关键设计决策**：`use_merged` 参数控制网格选择策略——CZM未启用时使用合并网格（消除界面节点重复，保证径向导热连续），CZM启用时使用未合并网格（保留界面自由度，通过界面热阻模型处理层间传热）。

## 后续变更 (2026-04-07)

- **MeshGeometry 构造**: `setup_thermal2D_mesh` 中 6 行 Dict 赋值 → `MeshGeometry` struct 构造 + `layer_weights` 计算
- `case.multi_spme_layout["layer_weights"] = ...` → `case.geometry = MeshGeometry(..., layer_weights, ...)`
- `layer_weights` 存储在 `MeshGeometry` struct 中，提供类型安全访问
- 行数从约 568 行减少到 562 行

## 后续变更 (2026-04-20)

- **预计算边界边缓存**: `setup_thermal2D_mesh` 中在构造 `MeshGeometry` 之前调用 `identify_boundary_nodes` 和 `compute_boundary_edge_cache`
- **MeshGeometry 新增 `boundary_edges` 字段**: 将预计算的 `BoundaryEdgeCache` 存入 `case.geometry.boundary_edges`
- 热模型对流边界条件装配可复用此缓存，消除每步重复扫描网格的开销
- 行数从约 562 行增加到约 569 行
