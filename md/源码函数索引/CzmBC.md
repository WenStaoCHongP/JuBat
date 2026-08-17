# CzmBC.jl

- **源文件**: `src/CzmBC.jl`
- **行数**: 63 行
- **函数/struct 计数**: 0 个 struct + 2 个独立函数
- **职责**: CZM 边界节点识别和 Dirichlet 边界条件施加。
- **相关技术文档**: `md/03_边界条件.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

本文件无独立 struct。

## 函数清单

### `apply_bc_czm(K, F; bc_nodes, bc_dofs, bc_vals) -> (K_new, F_new)` — L7-L42

使用相对罚系数对矩阵和载荷施加位移边界条件；支持节点级 `:fixed_x`、`:fixed_y`、`:fixed_xy`，也支持显式 DOF/数值输入。

### `identify_bc_nodes_czm(czm_mesh, param; opt, fix_inner=true)` — L44-L63

调用 `identify_boundary_nodes` 获取内外圈掩码；外圈始终固定，`fix_inner=true` 时同时固定内圈。

## 跨文件依赖

- `CouplingState.jl`：`CohesiveMesh`
- `Tools.jl`：`identify_boundary_nodes`

## 省略项

无。全部 function 均有独立条目。

### [DEBUG]

无。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|---|---|---|
| L13 | 全零矩阵时罚系数回退为 `1e12` | 正常装配路径不触发；异常输入下可能主导条件数 |

### [COMPLEX-CHECK]

无。
