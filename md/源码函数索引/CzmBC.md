# CzmBC.jl

- **源文件**: `src/CzmBC.jl`
- **行数**: 124 行
- **函数/struct 计数**: 0 个 struct + 2 个独立函数
- **职责**: CZM 边界节点识别和 Dirichlet 边界条件施加（相对罚实现）。
- **相关技术文档**: `md/03_边界条件.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

本文件无独立 struct。

## 函数清单

### `apply_bc_czm(K, F; bc_nodes, bc_dofs, bc_vals) -> (K_new, F_new)` — L7-L103

使用相对罚系数对矩阵和载荷施加位移边界条件；支持节点级 `:fixed_x`、`:fixed_y`、`:fixed_xy`，也支持显式 DOF/数值输入。

- 入口契约检查（L8-L62）：方阵、非空、每节点 2 DOF、`F` 长度匹配、`bc_nodes` 与 `bc_dofs/bc_vals` 二选一、索引/类型/重复约束/有限性逐一校验，违反即 `throw`
- 罚系数为**相对罚**（重设计 v2 §6）：`penalty = 1e6 * maximum(abs, diag(K))`（L64-L72），对角尺度非有限正数时抛错（不再有固定 `1e12` 兜底）
- 节点模式按 `:fixed_x/:fixed_y/:fixed_xy` 展开到 DOF（L76-L94）；DOF 模式直接施加（L95-L100）

### `identify_bc_nodes_czm(czm_mesh, param; opt, fix_inner=true) -> (bc_nodes, inner_count, outer_count)` — L105-L124

调用 `identify_boundary_nodes` 获取内外圈掩码；外圈始终固定（`:fixed_xy`），`fix_inner=true` 时同时固定内圈；返回计数用于诊断。

## 跨文件依赖

- `CouplingState.jl`：`CohesiveMesh`
- `Tools.jl`：`identify_boundary_nodes`

## 省略项

无。全部 function 均有独立条目。

### [DEBUG]

无。

### [PLACEHOLDER]

无。（早期版本"全零矩阵时罚系数回退 1e12"的兜底已改为显式抛错。）

### [COMPLEX-CHECK]

无。
