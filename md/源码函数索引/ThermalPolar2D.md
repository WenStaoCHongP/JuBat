# ThermalPolar2D.jl

- **源文件**: `src/ThermalPolar2D.jl`
- **行数**: 120 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 极坐标有限体积法（FVM）圆环热模型矩阵装配，基于径向节点序列与周向离散直接组装稀疏 `M`/`K`/`F`，外边界施加 Robin 对流
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/12_热模型验证方案.md`（圆环精确解）

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ThermalPolar2D_Ring(case::Case, variables::AbstractDict{String,<:Any}, mesh_data)` — L1-L120

极坐标 FVM 圆环热模型装配，返回 `(MT, KT, F)`。

- 输入 `mesh_data.r`（已归一化的径向节点）、`mesh_data.ntheta`（周向离散数）
- 物性：`rho_c_nd = heat_Q / volume`、`k_r_nd = lambda_r`、`k_t_nd = lambda_t`（均无量纲）
- 热源：将 `variables["heat_source_fields"]`（单元值）按节点计数平均得到节点值 `q_node`
- FVM 装配（双重循环 `ir × it`）：
  - 质量矩阵：`M[idx] = (ρc)*·V*`，`V* = 0.5·(r_iph² − r_imh²)·Δθ`
  - 载荷：`F[idx] += q_node·V*`
  - 径向邻居（ir>1 / ir<end）：conductance `a_r = k_r*·r_face·Δθ/Δr`
  - 外边界（ir == end）：施加 Robin `Bi·A_bc·(T − T_amb)`
  - 周向邻居（周期）：`a_t = k_t*·Δr_face/(r·Δθ)`，连接 `it±1`（含周期 wrap）
- 输出：`MT = spdiagm(M)`、`KT = sparse(I,J,V)`、`F`
- 几何校验：`dr_im>0`/`dr_ip>0`/`r_i>0`，违反时 `error(...)`
- 跨文件依赖：仅 `case.param` / `case.mesh` / `variables`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

无。注：双重循环 `ir × it`（L52、L70）天然产生 2 层嵌套，但内部 `if` 均为单层（边界判断），不构成 ≥3 层嵌套。
