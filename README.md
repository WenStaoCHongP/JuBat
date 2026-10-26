# JuBat
A Julia-based framework for battery modelling. 

Jubat is based on Newman's battery models, including the P2D model, the SPM and the SPMe. All equatins are solved by the finite element method (FEM) with 2nd order elements. Jubat enjoys the benefits from Julia language for speed and 2nd order FEM for accuracy.

Jubat is a compact, readable and fast framework for battery modelling and will be useful for the rapid development of advanced battery management system.

Thermal model and degradation models are under development.

Reference: W. Ai, Y. Liu, Improving the convergence rate of Newman’s battery model using 2nd order finite element method, J. Energy Storage. 67 (2023) 107512. https://doi.org/10.1016/j.est.2023.107512

## 配置取值

配置字符串按大小写精确匹配。

| 配置项 | 合法取值 | 适用范围 |
| --- | --- | --- |
| `opt.model` | `"SPM"`、`"SPMe"`、`"P2D"` | 电化学默认初始化（`opt.y0` 为空） |
| `opt.solveType` | `"Crank-Nicolson"`、`"forward"`、`"backward"` | 电化学差分格式 |

```julia
opt = JuBat.Option()
opt.model = "SPMe"
opt.solveType = "Crank-Nicolson"
```

本分支开发范围为 SPMe + 二维热 + CZM；P2D 虽有初始化分支，但尚未适配分流求解器。
`model="thermal"` 是 `Solve` 的独立纯热入口，需传入 `thermal_variables`，不属于上述电化学初始化。
非法值显式报错，不自动改写或使用默认值替代。

单 SPMe + `distributed2D` 组合已移除；二维 SPMe 耦合需在配置中显式设置 `per_element_spme=true`，并使用 `ModelInitialisation_MultiSPMe` 初始化。
SPMe 的 `none` / `lumped` 热模式及独立 `model="thermal"` 入口保留。

## 多 SPMe 输出

默认初始化将同一套电化学初态复制到各热单元，后续各支路使用自己的状态、温度和电流演化。
多 SPMe 不再对平均状态调用一套代表 SPMe；电化学输出直接收集各支路计算结果，不附加未定义的全局平均电化学量。

| 结果 | 多 SPMe 输出维度与来源 |
| --- | --- |
| `time [s]`、`cell current [A]`、`cell voltage [V]` | 时间向量；总电流与分流求解的公共电压 |
| `temperature [K]` | 体积平均温度时间向量；节点及单元温度场仍分别保留 |
| 正/负颗粒浓度、电解液浓度 `[mol/m^3]` | 内部自由度 × 热单元 × 时间 |
| 表面浓度、交换/界面电流密度、过电位、OCP | 热单元 × 时间；原有键名保留，新增 OCP 和界面电流密度的有量纲输出 |
| 颗粒应力、表面位移 | `mechanicalmodel="full"` 时输出热单元 × 时间；未启用该模型时不发布占位零结果 |
| `thermal2D element branch current [A]` | 各支路实际电流；与原有无量纲 `thermal2D element current` 区分 |

单元维度的顺序与 `case.mesh["thermal2D"].element` 行号一致。例如负极颗粒浓度的 `[:, e, it]`
表示第 `e` 个单元在第 `it` 个记录列的径向浓度分布，表面浓度的 `[e, it]` 表示同一单元和记录列的表面值。
时间层沿用现有显式交错求解约定：首列为初态，后续电化学量、颗粒应力和公共电压来自该步
`CallModel` 的步首已接受状态，电流在 `t_eval` 取样；`time [s]` 标记步末记录时间，节点温度记录步末值。
因此这些电化学输出不能当作重新计算的步末状态，末列浓度也不等同于 `final_state.y` 的浓度。
本批仅保留各单元维度，不改变求解及载荷采样时间层；顶层体积平均温度仍来自同次 `CallModel` 的热场。
内部历史仍为矩阵：每个单元的内部自由度连续排列，时间为列；后处理恢复上述公开维度。
历史记录的完整检查后写入顺序保留。单 SPMe 的键名、单位和输出形状保持原样。

分流接口改为 `solve_branch_currents(case, variables, t, I_total, areas, Te_prev, x_prev; yt_chem, D_elem)`，
移除未使用的代表状态位置参数；必需的逐支路状态仍通过 `yt_chem` 传入。
逐单元浓度历史的存储量随单元数增加；此输出维度变更需调用者显式选择单元，不能沿用旧代表量的线性索引。

## CZM 接口迁移（M7）

- 构网调用改为 `create_czm_mesh(czm_submesh, thermal_mesh)`，移除原先未使用的第三个 `param` 参数。
- `assemble_coupled_system(czm_mesh, u, param; czm_workspace, ...)` 返回 `(K_total, f_int_total, separations, tractions)`；移除原先未生效的 `F_ext`、`F_thermo_chem` 关键字。
- 需要载荷残差时，由调用者计算，或使用保留的 `assemble_coupled_system_full(czm_mesh, u, param, dT_elem, Δsoc_n_elem, Δsoc_p_elem; czm_workspace, F_ext, ...)`。该包装仍构造热化学载荷并计算 `R = F_external + F_thermo_chem - f_int_total`，返回 `(K_total, R, F_thermo_chem, separations, tractions)`；省略 `F_ext` 表示零外载荷。
