# CouplingState.jl

- **源文件**: `src/CouplingState.jl`
- **行数**: 750 行
- **函数/struct 计数**: 18 个（10 struct + 8 函数）
- **职责**: 类型安全的状态布局与网格几何定义——CZM 本构参数缓存（`CzmInterfaceParams`/`CzmParamCache`）、多 SPMe 状态向量布局索引（`MultiSPMeLayout`）、CZM 装配缓存与工作区（`CZMAssemblyCache`/`CZMAssemblyWorkspace`）、CZM 跨时间步损伤更新入口（`update_czm_damage!`）、粗热→CZM 节点双线性插值矩阵构造
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/06_内聚力模型_CZM.md`、`docs/planning-with-files/14_力学模块修改/宏观力学模块无量纲化重设计.md`

## 数据结构

### `@with_kw struct CzmInterfaceParams` — L30-L65

单一界面类型（如 `:PE_PCC` 或 `:NE_NCC`）的 CZM 本构参数（归一化后）。20 个字段覆盖 `Materialmatrix.jl` 实际读取的全部字段。

- 体模量与热化学载荷：`E_eff`、`ν`、`α`（L32-L34）
- 无量纲化派生量（重设计 v2）：`Λ = scale.L/scale.δ_czm`、`E_star`（调和平均双材料模量）、`L_ch`（内禀长度）（L37-L39）
- Mode I（法向）：`σ_max`、`K_n`、`δ_0_n`、`δ_c_n`、`G_c`（L42-L46）
- Mode II（切向）：`τ_max`、`K_t`、`δ_0_t`、`δ_c_t`、`G_c_t`（L49-L53）
- BK 混合模式 + 本构选择：`η`、`czm_model`（L56-L57）
- 界面热阻：`h_c0`、`k_air`、`lambda_m`、`beta`、`threshold`（L60-L64）

### `struct CzmParamCache` — L76-L80

按界面类型分组的 CZM 参数缓存（spec §3.5.2）。

- `by_interface::Dict{Symbol, CzmInterfaceParams}`：`:PE_PCC` 与 `:NE_NCC` 两键（L77）
- `param_ref::Params`：保留 param 引用，供 `assemble_bulk_stiffness` 读 `PE/NE.E_coat`（L78）
- `id::UInt64`：内容哈希 `hash((hash(pe_pcc), hash(ne_ncc)))`，Task 4.4 修正——原 `objectid(param)` 漏检原位修改（L79）

### `struct MultiSPMeLayout` — L87-L95

多 SPMe 状态向量的布局索引。初始化后不可变。

- `ne`（热单元数）、`n_chem`（每单元电化学 DOF 数）、`nT`（热节点 DOF 数）、`n_total`（L88-L91）
- `chem_range::UnitRange{Int}`、`thermal_range::UnitRange{Int}`（L92-L93）
- `areas::Vector{Float64}`：预计算的单元面积（网格不变量）（L94）

### `struct BoundaryEdgeCache` — L125-L128

预计算的外边界边列表（网格不变量），用于对流边界条件装配。

- `edges::Vector{Tuple{Int,Int}}`：`(node_a, node_b)` 对，`a < b`（L126）
- `L_edge::Vector{Float64}`：边长（无量纲）（L127）

### `struct MeshGeometry` — L163-L172

Jellyroll 网格的几何拓扑信息。构建后不可变。

- `element_layer::Vector{Int}`：层类型 1=NE/2=SP/3=PE/4=NCC/5=PCC（L164）
- `is_inner_layer::Vector{Bool}`（L165）
- `layer_weights::Matrix{Float64}`：`ne × 5` 层面积权重 `[NE, SP, PE, PCC, NCC]`（L166）
- `interface_pairs::Vector{Tuple{Int,Int}}`：CZM 界面配对（L167）
- `czm_element_map::Dict{Int,Vector{Int}}`：热单元号 → CZM 单元索引向量（L168）
- `inner_nodes`、`outer_nodes`、`boundary_edges`（L169-L171）

### `struct CohesiveElementGeom` — L183-L193

预计算的单个 cohesive 单元几何信息。构建后不变。包含单元长度、法/切向量、旋转矩阵 `R`、全局 DOF 编号 `[8]`、底/顶面节点、Gauss 权重/坐标。

### `mutable struct CZMAssemblyWorkspace` — L201-L238

CZM 每轮 Newton 迭代复用的工作区，避免单元级临时分配。所有中间矩阵/向量运算使用 `mul!` 复用预分配数组。

- 单元级：`u_e[8]`、`K_e[8×8]`、`f_int_e[8]`、`B_global[2×8]`、`B_local[2×8]`、`δ_local[2]`、`BL_dT[8×2]`、`BL_dT_B[8×8]`、`T_vec[2]`、`BLtT[8]`（L203-L212）
- 全局级：`f_int_coh`、`separations`、`tractions`（L214-L216）
- 预分配稀疏矩阵：`K_coh_buf`、`K_coh::SparseMatrixCSC`（L218-L219）
- 构造器 `CZMAssemblyWorkspace(ndof, n_coh)` 内部预算 `nnz_est = max(n_coh * 64, 1)`（L221-L237）

### `mutable struct CZMAssemblyCache` — L248-L266

CZM 求解器的静态/准静态缓存。失效判据基于 `czm_mesh_id` 与 `param_cache_id`——任一变化或 `fix_inner` 切换即重建。挂在 `Case.czm_cache` 上跨时间步复用。

- `K_bulk`、`bulk_dofs`、`cohesive_geom`、`bc_dofs`、`bc_vals`、`ws`、`fix_inner`、`valid`、`czm_mesh_id`、`param_cache_id`（L249-L258）

### `mutable struct CzmLayout` — L273-L277

CZM 求解的布局信息和跨时间步状态。对标电化学的 `MultiSPMeLayout`。

- `n_coh`（cohesive 单元数）、`ndof`（总位移 DOF 数 = 2·nnode）、`u_prev`（上一步位移场，跨时间步持有）（L274-L276）

### `mutable struct CZMSnapshot` — L646-L660

per-step CZM solver state for CSV export。所有物理值以归一化（无量纲）形式存储，denormalization 在 CSV 写出时通过 `case.param.scale` 完成。

- `time_s`（物理时间已还原）、`cycle`、`phase`（L647-L649）
- `displacement`、`damage`、`separation_n/t`、`traction_n/t`（L650-L655）
- `converged`、`iterations`、`residual_norm`、`method`（L656-L659）

## 函数清单

### `MultiSPMeLayout(ne, n_chem, nT)` — L98-L106

便捷构造器：自动计算 `chem_range`、`thermal_range` 和 `n_total`，`areas` 延迟填充为零向量。

### `MultiSPMeLayout(ne, n_chem, nT, mesh_th)` — L109-L118

便捷构造器：接收 mesh 计算单元面积。通过累加 `mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]` 得到每个单元的面积（L113-L116）。`@inbounds` 优化。

### `compute_boundary_edge_cache(mesh, is_outer) -> BoundaryEdgeCache` — L136-L156

从网格和外部节点标记中提取去重的外边界边列表。

- 遍历每个 Q4 单元的 4 条边（L145-L146），仅当两端节点均为外边界节点时收集（L147）
- `(a, b)` 归一为 `a < b` 后用 `Set` 去重（L148-L150）
- 边长 `hypot(dx, dy)`（L152）

### `CzmLayout(czm_mesh::CohesiveMesh)` — L280-L284

便捷构造器：从 `czm_mesh` 初始化 `n_coh` 与 `ndof = 2·nnode`，`u_prev = zeros(ndof)`。

### `compute_czm_params_per_interface(case) -> CzmParamCache` — L302-L401

按界面类型计算 CZM 参数。E_eff 用涂层模量（`PE.E_coat` / `NE.E_coat`），不做全栈均一化。

- 入口断言（L307-L313）：`PE/NE.E_coat > 0`、`scale.σ_czm > 0`、`scale.E_coat > 0`、`σ_max_* > 0`、`G_c_* > 0`、`K_n_* > 0`
- `Λ = scale.L / scale.δ_czm`（L322），断言 `scale.δ_czm > 0`（L321）
- `E_star`：界面双材料等效模量（调和平均）（L330-L331）
- `L_ch`：内禀长度 `E*·G_c/σ_max²`（L332-L333）
- K_0 下界判据（L339-L341）：`δ_0* > 0.1` 时 `@warn`（maxlog=1）
- 分别构造 `:PE_PCC` 与 `:NE_NCC` 的 `CzmInterfaceParams`（L343-L393）
- 内容哈希 `hash((hash(pe_pcc), hash(ne_ncc)))` 用于缓存失效检测（L399，Task 4.4 fix）

### `compute_czm_strain_inputs(case, variables, T_nodes) -> NamedTuple` — L419-L515

计算 CZM 体单元粒度的 `dT`、`Δsoc_p`、`Δsoc_n`，不生成或保留细力学节点温度场。

- `dT_thermal[e] = avg(T_nodes[nodes]) - T0`，通过 `thermal_elem_map` 直接取值到 CZM 单元
- `Δsoc_p/n` 通过 `thermal_elem_map` 直接取值，按 `material_type` 分发（L460-L474）；PCC/NCC/SP 保持 0

### `update_czm_damage!(case, variables, T_nodes_carry) -> result::CZMResult` — L535-L633

更新 CZM 网格的损伤状态（牛顿-拉弗森迭代 + 载荷子步法）。**严格契约版本**：输入或结果含非有限值、或求解未收敛时直接抛错，不静默降级。

- 同步 `czm_model` 选项（L542）
- 计算 per-interface 参数（L545-L546）；特征应变系数不再由此提取——ε₀ 按单元材料经
  `eigenstrain_of(param, mt)` 分层计算（2026-08-29 α/β 分层化，旧 α_eff/β_n/β_p 提取已删除）
- 构建或复用 CZM 缓存（L553，失效判据：objectid(czm_mesh) + param_cache.id + fix_inner）
- 应变输入 `compute_czm_strain_inputs`（L556-L559）
- **输入有限性检查**（L565-L570）：T / soc_n / soc_p / dT / Δsoc_n / Δsoc_p 任一含非有限值即 `throw(ArgumentError)`（早期 @warn 诊断已改为硬失败）
- `u_czm_prev` 从 `czm_layout` 取（L572-L573）
- Viscous regularization（L582-L594）：`β = Δs / (τ_v* + Δs)`，basic=1.0，arc_length/load_substep=1/n_load_steps
- 调用 `solve_czm_step`（L596-L602）
- debug 块（L604-L618）：打印 `max(δ_n)`、`max(δ_eff)`、`converged`、`β`、`D_max` 多个统计
- **结果有限性与收敛检查**（L620-L627）：位移/损伤/分离/牵引/残差非有限或 `!converged` 均 `error()`（替代旧"仅在收敛时提交"的部分提交路径）
- 检查通过后提交损伤状态与位移（L629-L630）

> 旧 6 参数兼容入口 `update_czm_damage!(czm_mesh, czm_params, case, ...)` 已在
> b4c0cde 重构中删除，全仓仅剩本 3 参数版本。


## 省略项

无。所有 struct 与 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L617 | `println("[CZM-Debug] max(δ_n)=$(round(max_delta_n; digits=6)), max(δ_eff)=...")` | 受 `case.opt.debug_coupling` 门控的运行时调试打印（L604-L618）；输出分离量、converged、β、D_max 等关键诊断 |

> 2026-08-18 复核：早期 L532/L540/L593 三处 NaN `@warn` 诊断已改为硬失败
> （输入 `throw(ArgumentError)` L565-L570；结果 `error()` L620-L627），不再计入 DEBUG。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L341 | `@warn "CZM 初始刚度过软：δ_0* > 0.1..." maxlog=1` | 参数验证类告警（非占位）；maxlog=1 抑制重复，不计入 PLACEHOLDER |
| L398 | 注释"原为 objectid(param)，但原位修改 param 字段不改变 objectid，导致漏检" | 说明性注释，无占位代码 |
| L547 | `case.czm_param_cache === nothing && (case.czm_param_cache = czm_param_cache)` | 惰性缓存写入；首次调用后后续不再更新——若 case.param.cohesive 字段运行时变化且未通过 `compute_czm_params_per_interface` 重算则缓存陈旧（依赖 Task 4.4 内容哈希失效机制兜底） |
| L617 | `[CZM-Debug] ... D_max(trrial)=...`（"trrial" 拼写错误） | 拼写笔误，输出字符串不影响逻辑；用户感知 |
| L719 | `abs(detJ) < 1e-20 && break`（Newton 迭代退化保护） | magic number 1e-20；退化单元几何近似奇异时退出，通常合理 |
| L723 | `return ξ, η, abs(ξ) <= 1.0 + 1e-6 && abs(η) <= 1.0 + 1e-6`（容差判定） | magic number 1e-6 容差；数值合理性 |

> 2026-08-18 复核：早期 `push!(V_vals, 1.0)`（Newton 失败回退最近粗热节点）已删除，
> 定位失败改为 `error`（L738-L740）。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L307 | 连续 7 个 `@assert param.X > 0 "..."` 入口断言（参数 positivity 检查，跨 L307-L313） | 抽出 `assert_positive(params, fields...)` helper；当前虽是独立断言，但模式重复且不易维护 |
| L565-L570 | 连续 6 行 `all(isfinite, X) \|\| throw(ArgumentError(...))` 输入有限性检查 | 模式重复；可抽出 `check_all_finite(name, xs)` helper 统一报错消息 |
| L586-L592 | `delta_s = if lowercase(iter_method) == "basic"; 1.0; elseif lowercase(iter_method) in ("arc_length", "arclength", "arc-length"); 1.0 / max(1, n_load_steps); else; 1.0 / max(1, n_load_steps); end`（嵌套 if-elseif-else + 字符串匹配多分支） | 抽出 `compute_delta_s(iter_method, n_load_steps)` helper；arc_length 与 else 分支结果相同可合并 |
| L143 | `for e in 1:ne; nodes = mesh.element[e, :]; for (a, b) in ((...), (...), (...), (...)); (is_outer[a] && is_outer[b]) \|\| continue; key = a < b ? (a, b) : (b, a); key in seen && continue; ...; end; end`（嵌套 2 层 + 多条件 continue，跨 L143-L153） | 抽出 `iter_boundary_edges(mesh, is_outer)` 迭代器 helper，主循环更简洁 |
| L730-L747 | 主循环：结构化父单元定位 + Newton 等参求解 + 越界 `error`（嵌套 2 层 + 多条件） | 可抽出 `interp_weights_in_parent(...)`；相比早期 nearest-10 候选版本已显著简化 |
