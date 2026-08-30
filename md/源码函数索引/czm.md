# czm.jl

- **源文件**：`src/czm.jl`
- **行数**：866 行
- **函数/struct 计数**：15 个函数，无本文件内 struct
- **职责**：逐层材料映射、cohesive/bulk 内力与切线装配、热化学本征应变载荷，以及随 `CohesiveMesh` 生灭的惰性装配缓存访问。
- **相关技术文档**：`md/06_内聚力模型_CZM.md`、`md/14_粘性正则化.md`、`md/15_颗粒与极片模量区分.md`

## 架构契约

- 装配函数直接读取 `param::Params`；`:PE_PCC` 由 `param.PCC` 提供界面参数，`:NE_NCC` 由 `param.NCC` 提供。
- `SetCase` 归一化后 `param` 冻结。线弹性 bulk 刚度、cohesive 几何标架和工作区分别惰性挂在 `CohesiveMesh.K_bulk/cohesive_geom/ws`；更换 mesh 对象即自然失效。
- BC 不缓存，由求解器每次入口调用 `identify_bc_nodes_czm` 现算。
- `DamageState`、`MechState` 与 `CZMAssemblyWorkspace` 定义在 `CouplingState.jl`；演化状态不属于 mesh。

## 函数清单

### `moduli_of(param, mt) -> (E, ν)` — L31-L47

按 `:PE/:NE/:SP/:PCC/:NCC` 返回层材料模量与泊松比；涂层读取 `E_coat/nu_coat`，连续层读取 `E/nu`，并转换到 CZM 应力空间。

### `eigenstrain_of(param, mt, dT, Δsn, Δsp) -> ε0` — L50-L69

按材料层计算热-化学本征应变。PE 只接收 `Δsp`，NE 只接收 `Δsn`；SP/PCC/NCC 仅接收各自热膨胀项。

### `cohesive_local_frame(czm_mesh, elem)` — L72-L115

根据界面节点计算长度、法向、切向和局部旋转矩阵。

### `assemble_czm_system(czm_mesh, u, param; ...)` — L118-L298

装配 cohesive 切线与内力并返回分离、牵引。循环内通过 `collector_params` 选择 collector，`czm_model` 显式传给双线性本构；尺度因子 `Λ = param.scale.L / param.scale.δ_czm` 在使用点内联。

### `assemble_bulk_stiffness(czm_mesh, param)` — L301-L380

按 `material_type` 直读各层材料参数，装配线弹性 Q4 平面应力 bulk 刚度。

### `gl_element_residual_tangent(...)` — L383-L491

单个 Q4 的 Green-Lagrange/Total-Lagrangian 残差与材料、初应力切线；支持 J2 试探状态与卷绕预应力。

### `assemble_bulk_residual_tangent(czm_mesh, u, param; ...)` — L494-L642

装配几何非线性 bulk 残差与切线。塑性试探态写入局部候选状态，最终提交由求解器控制。

### `assemble_thermal_chemical_load(czm_mesh, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem)` — L645-L713

在线性路径中把逐层本征应变转换为等效节点载荷。

### `assemble_coupled_system(czm_mesh, u, param; ...)` — L716-L764

组合 bulk 与 cohesive 切线和内力；可接收调用方取得的 `K_bulk_cached`、几何缓存和工作区。

### `assemble_coupled_system_full(czm_mesh, u, param, F_ext; ...)` — L767-L798

在耦合装配上加入外载与热化学载荷，形成完整残差。

### `collector_params(param, iface) -> CurrentCollector` — L810-L814

唯一界面分派：`:PE_PCC → param.PCC`、`:NE_NCC → param.NCC`；未知界面直接报错。

### `build_cohesive_geometry(czm_mesh)` — L821-L835

预计算全部 cohesive 单元的局部标架、DOF 和积分表。

### `bulk_stiffness(czm_mesh, param)` — L843-L846

首次调用装配并写入 `czm_mesh.K_bulk`，以后跨步只读。几何非线性/塑性路径不使用该线性缓存。

### `cohesive_geometry(czm_mesh)` — L853-L856

首次调用构建并写入 `czm_mesh.cohesive_geom`。

### `assembly_workspace(czm_mesh)` — L863-L866

首次调用分配并写入 `czm_mesh.ws`；工作区数值每轮覆写。

## 跨文件依赖

- `SetMesh.jl`：`CohesiveMesh` 与三项缓存字段。
- `CouplingState.jl`：`DamageState`、`CohesiveElementGeom`、`CZMAssemblyWorkspace`、塑性状态。
- `Materialmatrix.jl`：双线性牵引-分离本构、切线与 J2 更新。
- `CzmSolve.jl`：消费装配接口并负责收敛提交。

## 已删除接口

`CZMAssemblyCache`、`build_czm_cache`、`ensure_czm_cache`、参数内容哈希和 `case.czm_cache` 均已删除；文档或调用方不得重新引入这些中转层。
