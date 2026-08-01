# Initialisation.jl

- **源文件**: `src/Initialisation.jl`
- **行数**: 153 行
- **函数/struct 计数**: 5 个函数（无独立 struct）
- **职责**: 状态向量初始化——单模型 y0 构造（`ModelInitialisation`）、多 SPMe 并行架构扩展状态向量构造（`ModelInitialisation_MultiSPMe`）、多 SPMe 状态向量切片/更新辅助（`extract_element_state`/`get_thermal_dofs`/`update_state`）
- **相关技术文档**: `md/08_逐单元算法.md`、`md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。状态向量均为 `Vector{Float64}`，布局信息由外部 `MultiSPMeLayout`（定义于 `CouplingState.jl`）描述。

## 函数清单

### `ModelInitialisation(case::Case) -> Vector{Float64}` — L1-L48

为单模型（SPM / SPMe / P2D）构造初始状态向量 `y0`。

- 若 `case.opt.y0` 非空则直接返回用户提供的初值（L44-L45）
- **SPM**（L3-L8）：`y0 = [csn0; csp0]`，浓度初始化为 `param.NE.cs0` / `param.PE.cs0`
- **SPMe**（L9-L16）：`y0 = [csn0; csp0; ce0]`，额外含电解液浓度 `param.EL.ce0`
- **P2D**（L17-L29）：构造 `csn0`、`csp0`、`ce0`，并猜测 `phie0`（L26）、`phis_p`（L27，注释明示"guessed values are not used"）、`phis_n`（L28）；**注意**：`y0` 末尾仅追加 `[csn0; csp0; ce0]`（L29），电位猜测值 `phis_n/phis_p/phie0` 在 L41-L43 才条件性追加
- **热扩展**（L33-L40）：lumped 追加 `T0` 标量；distributed2D 追加 `nT` 个 `T0` 节点温度
- **P2D 电位追加**（L41-L43）：distributed2D 分支后，P2D 模型再追加 `[phis_n; phis_p; phie0]`
- 末尾 `vec(y0)` 展平返回（L47）

### `ModelInitialisation_MultiSPMe(case; initial_soc_distribution=nothing) -> Vector{Float64}` — L49-L113

为多 SPMe 并行架构初始化扩展状态向量。

- 获取 `ne`（热单元数）与 `nT`（热节点数）（L55-L57）
- **临时关闭 thermalmodel 获取纯电化学初值**（L60-L64）：保存 `original_thermalmodel`，设为 `"none"`，调 `ModelInitialisation(case)`，再恢复——避免热 DOF 混入单单元模板
- 计算每单元电化学 DOF 数 `n_chem`（L66）
- **逐单元初始化**（L75-L102）：
  - `initial_soc_distribution === nothing`（均匀）：直接复制 `y0_single_chem`（L78-L80）
  - 非均匀（L82-L101）：按 `soc_e` 调整粒子浓度
    - 负极：`cn_surf_e = NE.cs0 * soc_e`（线性近似，L87）
    - 正极：`cp_surf_e = PE.cs0 * (1 - soc_e)`（L92）
    - 电解液：均匀 `EL.ce0`（L96）
- 追加热场 `T0_nodes`（L105-L108）
- 缓存 `MultiSPMeLayout` 到 `case.layout`（L110）

### `extract_element_state(y, e, layout::MultiSPMeLayout) -> SubArray` — L120-L123

从多 SPMe 全局状态向量中提取单个单元 `e` 的电化学状态。`offset = (e-1) * layout.n_chem`，返回 `y[(offset+1):(offset+n_chem)]`。

### `get_thermal_dofs(y, layout::MultiSPMeLayout) -> SubArray` — L130-L132

从多 SPMe 全局状态向量中提取热场节点温度。返回 `y[layout.thermal_range]`。

### `update_state(y, layout; element_index, element_state, thermal_nodes) -> Vector{Float64}` — L139-L152

更新多 SPMe 全局状态向量（返回新向量，非原位修改）。

- `y_new = copy(y)`（L140）
- `element_index` 非 nothing 时：断言范围与长度（L142-L143），写入对应单元切片（L144-L145）
- `thermal_nodes` 非 nothing 时：断言长度（L148），写入 `thermal_range`（L149）
- 两个参数可同时提供，按顺序应用

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info` / `@warn`。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L26 | P2D 初始电位猜测：`phie0 = -ones * NE.U(cs0)`（L26）；`phis_p = ones * PE.U(PE.cs0) .+ phie0[1]`（L27）；`phis_n = zeros`（L28）。注释明示"guessed values are not used"（跨 L26-L28） | 占位猜测值：实际 P2D 求解器会覆盖这些值；若求解器初始化逻辑变化可能依赖这些猜测。L27 `.+ phie0[1]` 是为了对齐参考电位，语义隐晦 |
| L62 | 临时修改 `case.opt.thermalmodel = "none"` 再恢复以获取纯电化学初值（跨 L62-L64） | 临时 mutation 模式：若 `ModelInitialisation(case)` 抛异常则 `original_thermalmodel` 不会被恢复（无 try-finally 保护）。实践中 `ModelInitialisation` 不太可能抛异常，但模式不健壮 |
| L87 | `cn_surf_e = case.param.NE.cs0 * soc_e`（负极浓度线性近似） | 简化模型：注释明示"假设线性关系"；实际 cs0 是初始表面浓度，乘 soc_e 仅在 soc_e 接近 1 时近似有效。极端 soc_e 下偏差大 |
| L92 | `cp_surf_e = case.param.PE.cs0 * (1.0 - soc_e)`（正极浓度简化） | 简化模型：注释明示"简化模型"；与 L87 对称假设，极端 soc_e 下偏差大 |
| L96 | `ce0_e = ones * case.param.EL.ce0`（电解液假设均匀） | 简化假设：注释明示"假设均匀"；实际循环过程中电解液浓度可能不均匀，初始化时均匀合理 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L3 | `if isempty(case.opt.y0); if case.opt.model == "SPM"; ...; elseif case.opt.model == "SPMe"; ...; elseif case.opt.model == "P2D"; ...; else; error(...); end; if case.opt.thermalmodel == "lumped"; ...; elseif case.opt.thermalmodel == "distributed2D"; ...; end; if case.opt.model == "P2D"; ...; end; else; y0 = case.opt.y0; end`（嵌套 3 层 if + 多个独立 if 块串联，跨 L3-L46） | 抽出 `init_chem_state(case)` / `init_thermal_state(case, y0)` / `init_p2d_potentials(case, y0)` 三个独立函数，主函数组合调用；当前单函数 ~45 行且多模型分支混合 |
| L75 | `for e in 1:ne; offset = ...; if initial_soc_distribution === nothing; y0_chem_all[...] .= y0_single_chem; else; soc_e = ...; cn_surf_e = ...; csn0_e = ...; cp_surf_e = ...; csp0_e = ...; ce0_e = ...; yt_e = [...]; y0_chem_all[...] .= yt_e; end; end`（嵌套 2 层 + 多变量赋值，跨 L75-L102） | 抽出 `init_element_chem(soc_e, param, Nrn, Nrp, Nel)` helper；当前均匀/非均匀分支结构不对称（一个整体拷贝，一个逐字段构造） |
| L139 | `function update_state(y, layout; element_index::Union{Nothing,Int}=nothing, element_state::Union{Nothing,Vector{Float64}}=nothing, thermal_nodes::Union{Nothing,Vector{Float64}}=nothing)`（3 个 Union{Nothing,T} 关键字参数 + 2 个独立 if 块各自断言，跨 L139-L152） | 拆分为 `update_element!(y, layout, idx, state)` 与 `update_thermal!(y, layout, nodes)` 两个函数；当前 Union{Nothing} 模式要求调用者记得传 nothing，易出错 |
