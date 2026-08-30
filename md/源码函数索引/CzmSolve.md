# CzmSolve.jl

- **源文件**：`src/CzmSolve.jl`
- **行数**：989 行
- **函数/struct 计数**：1 个 struct + 14 个函数
- **职责**：CZM Newton-Raphson、载荷子步与弧长求解；边界条件现算；试探状态回滚；收敛后向 `MechState` 原位提交。
- **相关技术文档**：`md/06_内聚力模型_CZM.md`、`md/14_粘性正则化.md`、`md/09_分流求解器.md`

## 数据结构

### `mutable struct CZMResult` — L1-L15

单步结果容器：位移、损伤、法/切向牵引与分离、`converged`、迭代数和残差范数。默认构造令 `converged=false`、`residual_norm=Inf`。

## 状态与缓存契约

- 入口接收 `ms::MechState`；各策略先复制 `ms.u_prev` 与 `ms.damage_states` 作为试探态。
- 只有收敛才写回 `ms.u_prev/ms.damage_states`（以及启用塑性时的收敛塑性状态）；失败丢弃试探态，不触碰 `ms`。
- 求解器只返回 `CZMResult`，不克隆或返回新的 `CohesiveMesh`。
- BC 由 `extract_bc_dofs` 每次入口现算。线性 bulk 刚度、cohesive 几何和工作区通过 mesh 访问器惰性获取。

## 函数清单

### `clone_damage_states(damage_states)` — L17-L29

深拷贝损伤历史，供试探/回滚使用。

### `update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta)` — L38-L56

按 cohesive `interface_type` 分组，分别用 `collector_params(param, iface)` 调用 `update_damage`，保持原始单元顺序。

### `extract_bc_dofs(czm_mesh, param; fix_inner)` — L64-L79

调用 `identify_bc_nodes_czm` 现算 Dirichlet DOF/值；BC 不进入缓存。

### `backtrack_line_search!(...)` — L89-L111

basic 路径的回溯线搜索；试验位移上重装残差，返回 `(u_new, R_new_norm, accepted, α_used)`。

### `apply_czm_dirichlet!(u, bc_dofs, bc_vals)` — L113-L118

原位写入 Dirichlet 位移。

### `zero_czm_bc_entries!(v, bc_dofs)` — L120-L125

原位清零受约束条目。

### `fill_czm_result!(result, u, damage_states, separations, tractions)` — L127-L137

把已验收状态写入结果容器。

### `build_arc_length_augmented_matrix(...)` — L139-L149

构造弧长法增广矩阵。

### `spherical_arc_length_correction(...)` — L151-L173

计算球形弧长约束修正；维度、半径、系数或结果非法时显式失败。

### `solve_czm_basic_step(czm_mesh, F_ext, param, ms; ...) -> CZMResult` — L175-L277

单步 Newton + 回溯。冻结损伤求位移，收敛后更新损伤并提交 `ms`；失败返回 `converged=false` 且不修改 `ms`。

### `solve_czm_arc_length_step(czm_mesh, F_ext, param, ms; ...) -> CZMResult` — L280-L483

线性几何路径的弧长增量求解；自适应缩小步长，最终收敛后提交状态。

### `newton_raphson_czm(czm_mesh, F_ext, param, ms; ...) -> CZMResult` — L498-L667

自适应载荷子步 Newton 路径；子步试探可回滚，达到完整载荷且通过最终残差门后才提交。

### `solve_czm_arc_geo_step(czm_mesh, F_ext, param, ms; ...) -> CZMResult` — L679-L935

几何非线性弧长路径；以热/化学本征应变作为载荷参数化对象，先求参考平衡态，再推进到 `λ=1`，验收后提交。

### `solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; ...) -> CZMResult` — L943-L989

统一入口。直接从 `CzmOptions` 展开 `model/iter_method/max_iter/tol/load_steps/arc_length_alpha/viscous_enabled/viscous_tau/fix_inner/geo_nonlinear/j2_plasticity` 并分派：

- `basic` → `solve_czm_basic_step`
- `load_substep` → `newton_raphson_czm`
- `arc_length` 且 `geo_nonlinear=false` → `solve_czm_arc_length_step`
- `arc_length` 且 `geo_nonlinear=true` → `solve_czm_arc_geo_step`

粘性系数由 `opt.czm.viscous_tau / param.scale.t0` 与载荷步长现场计算；未知方法直接报错。

## 已删除接口

`clone_czm_mesh_with_damage`、`param_cache` 参数、`cache` 参数和 `(result, updated_czm_mesh)` 双返回值均已删除。调用方应持有同一 `case.czm_mesh`，并从 `case.mech` 读取已提交状态。
