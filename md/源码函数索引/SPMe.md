# SPMe.jl

- **源文件**: `src/SPMe.jl`
- **行数**: 291 行
- **函数/struct 计数**: 5 个独立函数；0 个 struct
- **职责**: SPMe (Single Particle Model with electrolyte) 电化学模型主体，计算颗粒扩散 + 电解液扩散的质量/刚度矩阵与源项，提供全局与逐单元两个入口；含力-化学耦合接口（应力影响扩散系数 θ_M）
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `SPMe(case::Case, yt, t; jacobi)` — L1-L35

全局 SPMe 求解入口，返回 `(M, K, F, variables)`。

- 先调 `SPMe_variables` 计算中间量；若 `mechanicalmodel == "full"`，再调 `Mechanicaloutput` 获取应力耦合系数 `theta_Mn/theta_Mp`
- 颗粒扩散：`ElectrodeDiffusion(param.NE/PE, mesh, nlen, c_gs, theta_M)`，时间尺度归一化 `M .*= ts_*/t0`
- 电解液扩散：`ElectrolyteDiffusion(param, mesh_el, nlen, variables)`
- 矩阵装配：`blockdiag(M_np, M_pp, M_el)`
- `jacobi == "constant"` 分支：复用 `param.NE.M_d/K_d`（假设 M、K 不随浓度变化）
- 跨文件依赖：`SPMe_variables`、`Mechanicaloutput`、`ElectrodeDiffusion`、`ElectrolyteDiffusion`、`SPMe_BC`

### `SPMe_element(case, yt_e, t, e; I_e, T_e, jacobi, workspace)` — L37-L93

逐单元（per-element）SPMe 入口，供多 SPMe 并行架构使用；逻辑与 `SPMe` 一致但支持外部电流 `I_e`、温度 `T_e` 注入与 workspace 复用。

- `workspace !== nothing` 时调 `SPMe_variables!`（原位），否则调 `SPMe_variables`（新建）
- 末尾 `variables_e["element index"] = Float64(e)` 便于调试（见 [DEBUG]）
- 跨文件依赖：`SPMe_variables!`、`SPMe_variables`、`Mechanicaloutput`、`ElectrodeDiffusion`、`ElectrolyteDiffusion`、`SPMe_BC`

### `SPMe_variables!(ws, case, yt, t; I_app, T_e)` — L101-L182

`SPMe_variables` 的原位变体：直接写入预分配 workspace Dict，省去 `StandardVariables` 分配。计算逻辑与 `SPMe_variables` 逐行一致。

- 状态提取仅覆写 ws 中已有键（跳过 thermal2D 键）
- 计算：交换电流密度 `j0 = k·Arrhenius·sqrt(cn(1-cn)·ce)`、过电位 `η = 2T·asinh(j/2j0)`、电解液电势降 `dphi_e`（PyBaMM 公式）、端电压 `V_cell = u_p - u_n + η_p - η_n + dphi_e`
- 熵热项：`u = U(cn) + (T-T0)·dUdT(cn)`
- 跨文件依赖：`Arrhenius`、`IntV`、`StandardVariables`（间接）

### `SPMe_BC(case, variables)` — L184-L207

构造 SPMe 的源项向量 `F = [flux_np; flux_pp; flux_el]`。

- 颗粒表面通量：`flux_np[end] = -j_n·rs^2`（负极）、`flux_pp[end] = -j_p·rs^2`（正极）
- 电解液源项：分段系数 `coeff[v_ne] *= (1-t⁺)·as·j_n`，`v_sp *= 0`，`v_pe *= (1-t⁺)·as·j_p`；用 `Assemble1D` 装配
- 跨文件依赖：`Assemble1D`

### `SPMe_variables(case, yt, t; I_app, T_e)` — L209-L291

全局 SPMe 变量计算（非原位）。逻辑同 `SPMe_variables!`，但用 `StandardVariables(case, 1)` 新建 Dict。

- `I_app` 缺省时从 `case.opt.Current` 取
- `T_e` 注入时从 `var_list` 过滤掉 `"temperature"` 键
- 额外：`thermal_distributed = thermal_enabled && thermalmodel == "distributed2D"`；非分布式时写 `variables["thermal2D element current"]`
- 跨文件依赖：`StandardVariables`、`Arrhenius`、`IntV`

## 省略项

无（本文件全部 5 个函数均含实质逻辑，无 trivial getter）。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L89 | `# 8) 可选：添加单元编号到 variables_e（用于调试）`；L90 `variables_e["element index"] = Float64(e)` | 在逐单元模式中将单元编号写入 variables，便于事后定位/排查特定单元；仅诊断用途，对求解无影响 |

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L62 | `if jacobi == "constant" && !isempty(param.NE.M_d) && !isempty(param.NE.K_d)` 三条件 `&&` 链 | 抽 `should_reuse_constant_jacobi(param)` 辅助函数，或缓存判定结果到 `case` |
| L128 | `T = T_e === nothing ? only(yt[case.index["temperature"]]) : T_e` 单行嵌套 `nothing` 检查 + `only` 解包 | 可提取 `_resolve_T(yt, case, T_e)` helper，与 L231-L235 的同义逻辑合并 |

---
