# CzmMesh.jl

- **源文件**：`src/CzmMesh.jl`
- **行数**：178 行
- **函数/struct 计数**：1 个 struct + 1 个函数
- **职责**：识别 PE-PCC / NE-NCC 真实界面、复制界面节点、重写外层 bulk 连接，并建立 cohesive→热单元映射。
- **相关技术文档**：`md/02_几何与网格.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

### `mutable struct CohesiveElement <: AbstractCohesiveElement` — L3-L12

COH2D4 拓扑：四节点、上下表面节点、长度、`interface_type` 与相邻内/外 bulk 单元。两种 `interface_type` 是参数类型；每个 8 层重复单元仍有四个真实箔–涂层面。

### 跨文件：`CzmSubmesh` / `CohesiveMesh` — `src/SetMesh.jl:L42-L84`

两者为避免 include 顺序问题定义在 `SetMesh.jl`。`CohesiveMesh` 保存拓扑、映射及三项惰性缓存：

- `K_bulk`：线弹性 bulk 刚度
- `cohesive_geom`：纯几何局部标架
- `ws`：预分配装配工作区

缓存随 mesh 对象生灭；mesh 不保存损伤或塑性状态，演化状态归 `MechState`。

## 函数清单

### `create_czm_mesh(czm_submesh, thermal_mesh, param) -> CohesiveMesh` — L33-L178

- 从 `czm_submesh.mesh_bonded` 建立共边→bulk 单元对映射。
- 识别每个周向分段的四个真实 PE-PCC/NCC-NE 界面，归入 `:PE_PCC/:NE_NCC`。
- 按质心半径判定内外层，复制界面节点并重写外层 bulk 连接以允许分离。
- 建立 `cohesive_to_thermal`；温度/SOC 的 bulk 映射由 `CzmSubmesh.thermal_elem_map` 提供。
- 构造时三个缓存保持 `nothing`；首次求解由 `bulk_stiffness/cohesive_geometry/assembly_workspace` 惰性填充。
- 不创建 `DamageState`；`MechState(czm_mesh)` 负责状态初始化。

## 跨文件依赖

- `SetMesh.jl`：`CzmSubmesh`、`CohesiveMesh`、`Mesh`。
- `CouplingState.jl`：`MechState` 与损伤状态。
- `czm.jl`：mesh 缓存访问器和装配。
