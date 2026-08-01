# SPM.jl

- **源文件**: `src/SPM.jl`
- **行数**: 85 行
- **函数/struct 计数**: 3 个独立函数；0 个 struct
- **职责**: Single Particle Model (SPM) 电化学模型主体，仅求解正/负极颗粒扩散，不含电解液浓度/电势方程；与 SPMe 共享颗粒扩散装配流程
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `SPM(case::Case, yt, t; jacobi)` — L1-L31

全局 SPM 求解入口，返回 `(M, K, F, variables)`。

- 先调 `SPM_variables` 计算中间量；若 `mechanicalmodel == "full"`，再调 `Mechanicaloutput` 获取应力耦合系数 `theta_Mn/theta_Mp`
- 颗粒扩散：`ElectrodeDiffusion(param.NE/PE, mesh, nlen, c_gs, theta_M)`，时间尺度归一化 `M .*= ts_*/t0`
- 矩阵装配：`blockdiag(M_np, M_pp)`、`blockdiag(K_np, K_pp)`（无电解液块）
- `jacobi == "constant"` 分支：复用 `param.NE.M_d/K_d`
- 跨文件依赖：`SPM_variables`、`Mechanicaloutput`、`ElectrodeDiffusion`、`SPM_BC`

### `SPM_BC(case, variables)` — L34-L47

构造 SPM 颗粒扩散源项 `F = [flux_np; flux_pp]`。

- 负极表面通量：`flux_np[end] = -j_n·rs^2`
- 正极表面通量：`flux_pp[end] = -j_p·rs^2`
- 跨文件依赖：无（仅数组装配）

### `SPM_variables(case, yt, t)` — L49-L86

计算 SPM 中间变量并填入 `StandardVariables`。

- 电流密度：`j_n = I_app/(as·thickness)`，`j_p = -I_app/(as·thickness)`
- 交换电流：`j0 = k·Arrhenius·sqrt(cn·ce0·|1-cn|)`
- 过电位：`η = 2T·asinh(j/2j0)`；端电压 `V = u_p - u_n + η_p - η_n`
- 熵热项：`u = U(cn) + (T-T0)·dUdT(cn)`
- 温度：从 `case.index["temperature"]` 取，无则用 `cell.T0`
- 跨文件依赖：`StandardVariables`、`Arrhenius`

## 省略项

无（本文件全部 3 个函数均含实质逻辑）。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

无。

---
