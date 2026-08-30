# CouplingState.jl

- **源文件**：`src/CouplingState.jl`
- **行数**：503 行
- **函数/struct 计数**：8 个 struct + 8 个顶层函数
- **职责**：多 SPMe 布局、热网格几何缓存、CZM 几何工作区、力学演化状态、热/SOC→力学载荷映射，以及跨时间步 CZM 更新入口。
- **相关技术文档**：`md/10_参数传递与模块架构.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

### `struct MultiSPMeLayout` — L12-L20

多 SPMe 全局状态向量布局：单元数、每单元化学 DOF、热 DOF、总长度、化学/热区间和预计算单元面积。

### `struct BoundaryEdgeCache` — L50-L53

外边界无向边与无量纲边长；用于热对流边界装配。

### `struct MeshGeometry` — L88-L97

Jellyroll 热网格的层类型、内外层标记、层权重、接口/映射、边界节点和边缓存。

### `struct CohesiveElementGeom` — L108-L118

cohesive 单元的长度、局部标架、旋转矩阵、DOF、上下表面节点和 Gauss 表。该纯几何对象由 `czm.jl` 构建并缓存到 `CohesiveMesh.cohesive_geom`。

### `mutable struct CZMAssemblyWorkspace` — L126-L163

Newton 装配复用的单元/全局缓冲区与稀疏矩阵工作区。实例惰性挂到 `CohesiveMesh.ws`，每轮重写数值。

### `mutable struct DamageState` — L181-L192

逐 cohesive 单元的 `D/D_visc/δ_max_n/δ_max_t/δ_max_eff/fractured/accumulated_damage`。状态由 `MechState` 持有，不属于 mesh。

### `mutable struct MechState` — L201-L208

力学演化状态聚合：

- `u_prev`：上一步收敛位移
- `damage_states`：逐 cohesive 单元损伤
- `plastic_states`：可选 PCC/NCC Gauss 点塑性历史
- `winding_prestress`：可选卷绕预应力场
- `node_ref`：初始螺旋节点快照
- `contact`：SP–涂层接触预留位

求解器只在收敛后原位提交，失败/试探不得修改该结构。

### `mutable struct CZMSnapshot` — L489-L503

CSV 导出的逐步 CZM 快照；位移、损伤、分离和牵引以无量纲形式保存，写出时还原。

## 函数清单

### `MultiSPMeLayout(ne, n_chem, nT)` — L23-L31

建立状态区间与总长度，面积延迟为零向量。

### `MultiSPMeLayout(ne, n_chem, nT, mesh_th)` — L34-L43

额外由 Gauss 权重/Jacobian 预计算热单元面积。

### `compute_boundary_edge_cache(mesh, is_outer)` — L61-L81

抽取、去重外边界边并计算长度。

### `MechState(czm_mesh)` — L211-L216

按 mesh 的 DOF 与 cohesive 数量创建零位移、零损伤初态；其余可选状态为 `nothing`。

### `ensure_node_ref!(case)` — L223-L228

惰性保存初始螺旋节点，之后不重置；用于核心椭圆化相对基准。

### `core_ovalization(czm_mesh, u, ref_node)` — L237-L256

在第一匝窗口计算径向位移，去除 0/1 阶刚体分量后返回 `w_core` 与归一化 `Δ_core`。

### `compute_czm_strain_inputs(case, variables, T_nodes)` — L282-L370

把活动热节点温度和热单元 SOC 经 `CzmSubmesh.thermal_elem_map` 映射为逐力学 bulk 单元的 `dT_czm/Δsoc_p_czm/Δsoc_n_czm`。尺寸、索引和有限性不满足契约时直接失败；不生成细力学节点温度场。

### `update_czm_damage!(case, variables, T_nodes_carry) -> CZMResult` — L390-L476

耦合入口：读取 `case.czm_mesh/case.mech/case.param/case.opt.czm`，计算应变输入，建立可选预应力/塑性状态，并调用：

```julia
solve_czm_step(case.czm_mesh, case.mech, case.param, F_ext, case.opt.czm;
               dT_elem=..., Δsoc_n_elem=..., Δsoc_p_elem=...)
```

返回值或输入非有限、或 `result.converged == false` 时硬失败。求解器已完成损伤/位移收敛提交；J2 塑性在验收后通过 `commit_plastic=true` 提交。

## 缓存与参数边界

- 界面本构/热阻参数直接位于 `param.PCC/param.NCC`；本文件不构造参数缓存。
- 装配缓存随 `CohesiveMesh`，演化状态随 `MechState`。
- `SetCase` 后参数冻结保证 `K_bulk` 新鲜度；BC 每次求解入口现算。

## 已删除接口

`CzmInterfaceParams`、`CzmParamCache`、`compute_czm_params_per_interface`、`CZMAssemblyCache`、`CzmLayout` 及对应 Case 字段均已删除。
