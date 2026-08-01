# czm.jl

- **源文件**: `src/czm.jl`
- **行数**: 873 行
- **函数/struct 计数**: 2 个 struct + 11 个独立函数
- **职责**: CZM 本构核心与系统装配——`CohesiveElement` / `DamageState` 数据结构、`create_czm_mesh` 网格构造、`assemble_czm_system` 内聚力装配、`assemble_bulk_stiffness` / `assemble_thermal_chemical_load` 体单元装配、`build_czm_cache` / `ensure_czm_cache` 缓存、`assemble_coupled_system_full` 全耦合残差、`apply_bc_czm` / `identify_bc_nodes_czm` 边界条件
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`md/14_粘性正则化.md`、`md/15_颗粒与极片模量区分.md`

## 数据结构

### `mutable struct CohesiveElement <: AbstractCohesiveElement` — L1-L10

单个 COH2D4 内聚力单元几何与拓扑。

- 字段：`id`、`nodes[4]`、`nodes_bottom[2]`、`nodes_top[2]`、`length`、`interface_type`（`:PE_PCC` / `:NE_NCC`）、`host_outer_elem`、`host_inner_elem`
- `nodes` 顺序 `[n_lo, n_hi, n_hi_copy, n_lo_copy]`（逆时针）；底面对应内层 bulk，顶面对应外层副本

### `mutable struct DamageState <: AbstractDamageState` — L26-L37

每内聚力单元（或高斯点）的损伤历史。

- 字段：`D`（等效损伤 [0,1]）、`D_visc`（粘性正则化有效损伤，`md/14_粘性正则化.md`）、`δ_max_n/t/eff`、`fractured`、`accumulated_damage`
- 内部默认构造 `DamageState()` 全零初始化（L36）

## 函数清单

### `create_czm_mesh(czm_submesh, thermal_mesh, param) -> CohesiveMesh` — L58-L207

基于 CzmSubmesh 构造内聚力网格（spec §4.4）。

- Step 1（L64-L72）：建立 `共边(2 节点) → 单元对` 映射
- Step 2（L74-L101）：遍历共边识别 PE-PCC / NE-NCC 径向界面；用单元质心 `r=hypot(cx,cy)` 判断径向方向（L93-L94）
- Step 3（L103-L175）：节点复制 + 重写外层 bulk 连接（关键：否则分离位移恒为 0）；`sort!(common)` 用 node id 排序避免 `atan(y,x)` ±π 分支切割 bug（L130）
- Step 4（L180-L196）：组装 `CohesiveMesh`，调 `build_thermal_to_czm_interp`
- 末尾正确性自检（L199-L204）：副本坐标一致性、4 节点不重复
- 跨文件依赖：`CohesiveMesh`、`build_thermal_to_czm_interp`、`CzmSubmesh`

### `moduli_of(param, mt::Symbol) -> (E, ν)` — L238-L246

按材料类型查表返回体模量/泊松比，统一到 CZM 应力空间（σ_czm 参考）。

- `:PE` / `:NE` 返回涂层模量 `E_coat` / `nu_coat`（CLAUDE.md §9.4，颗粒 vs 极片区分）
- `:SP` / `:PCC` / `:NCC` 返回连续层 `E` / `nu`
- 双重再缩放因子 `s = scale.E_coat / scale.σ_czm`（>0 时）；σ_czm 未锚定时回退 1.0
- 注：`alphaT` 已从此函数移除（I2-a 修复，见 docstring L234-L236）

### `assemble_czm_system(czm_mesh, u, param_cache; damage_states, geom_cache, ws, visc_beta) -> (K_coh, f_int_coh, separations, tractions)` — L261-L443

内聚力单元全局刚度与内力装配，按 `interface_type` 从 `param_cache.by_interface` 取本构参数。

- 首次调用构建稀疏 sparsity pattern（L285-L310），后续仅 `fill!(nonzeros(K_coh), 0.0)` 清零
- `Λ`：位移空间 → 分离空间换算因子（重设计 v2 §5，L323）；虚功一致性：δ̃ = Λ·B·ũ，切线刚度含一次 Λ
- 单元循环（L315-L440）使用 `mul!` 零分配策略；`bilinear_traction_state` + `bilinear_tangent`（在 `Materialmatrix.jl`）
- 缓存支持：`geom_cache::Vector{CohesiveElementGeom}` / `ws::CZMAssemblyWorkspace`
- 跨文件依赖：`bilinear_traction_state`、`bilinear_tangent`、`NCweight`、`CZMAssemblyWorkspace`

### `assemble_bulk_stiffness(czm_mesh, param_cache) -> K_bulk` — L454-L524

固体 Q4 单元弹性刚度装配，按 `czm_submesh.material_type` 分组取模量。

- 平面应力 D 矩阵 L479-L481
- 调 `IntQ4` 高斯积分；每单元 8×8 块按 DOF 累加
- 缓存此矩阵是 `build_czm_cache` 的最高 ROI（L625）
- 跨文件依赖：`IntQ4`、`moduli_of`

### `assemble_thermal_chemical_load(czm_mesh, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem) -> F_thermo_chem` — L533-L597

热-化学载荷向量装配（初始应变 `ε_0 = α·ΔT + β_n·Δsoc_n + β_p·Δsoc_p`）。

- 按材料类型查表模量（L564）；α_eff / β_n / β_p 为位置参数（跨材料统一）
- 载荷公式：`F = ∫ B^T D ε_0 dΩ`，单元循环 L562-L594

### `build_czm_cache(czm_mesh, param_cache; fix_inner=true) -> CZMAssemblyCache` — L621-L704

构建 CZM 装配缓存（K_bulk、bulk_dofs、cohesive_geom、bc_dofs/vals、ws）。

- 失效判据：`czm_mesh_id = objectid(czm_mesh)` 与 `param_cache.id` 共同决定（L694-L695）
- 步骤 1-6 依次填充 `cache.K_bulk` / `bulk_dofs` / `cohesive_geom` / `bc_dofs` / `ws`
- 跨文件依赖：`assemble_bulk_stiffness`、`identify_bc_nodes_czm`、`CohesiveElementGeom`、`NCweight`

### `ensure_czm_cache(case, czm_mesh, param_cache; fix_inner=true) -> cache` — L719-L729

惰性构建/刷新 `case.czm_cache`，失效条件任一触发即调 `build_czm_cache`。

- 失效条件：`cache === nothing` / `!valid` / mesh objectid 变 / `param_cache.id` 变 / `fix_inner` 切换（L721-L724）

### `assemble_coupled_system(czm_mesh, u, param_cache; F_ext, F_thermo_chem, damage_states, K_bulk_cached, geom_cache, ws, visc_beta) -> (K_total, f_int_total, separations, tractions)` — L737-L769

体刚度 + 内聚力耦合装配（无热化学载荷）。

- `K_total = K_bulk + K_coh`（L763）
- `f_int_total = K_bulk·u + f_int_coh`（L766）
- K_bulk 支持缓存透传避免重复计算

### `assemble_coupled_system_full(czm_mesh, u, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem; F_ext, ...) -> (K_total, R, F_thermo_chem, separations, tractions)` — L777-L810

`assemble_coupled_system` 的扩展版，含热-化学载荷与残差计算。

- 内部调 `assemble_coupled_system` + `assemble_thermal_chemical_load`
- 残差 `R = F_external + F_thermo_chem - f_int_total`（L807）

### `apply_bc_czm(K, F; bc_nodes, bc_dofs, bc_vals) -> (K_new, F_new)` — L817-L852

对 `(K, F)` 应用 Dirichlet 边界条件（相对罚方法）。

- 相对罚 `penalty = 1e6 · max(|diag(K)|)`（重设计 v2 §6，L820-L823），矩阵全零时回退 `1e12`
- 支持两种 BC 形式：`bc_nodes::Dict{Int,Symbol}`（:fixed_x/:fixed_y/:fixed_xy）或 `(bc_dofs, bc_vals)` 对

### `identify_bc_nodes_czm(czm_mesh, param; opt, fix_inner=true) -> (bc_nodes, inner_count, outer_count)` — L854-L873

识别 CZM 网格的边界节点（内圈 + 外圈）。

- 委托 `identify_boundary_nodes` 取 `is_inner` / `is_outer` 掩码
- `fix_inner=true` 时内圈节点固定为 `:fixed_xy`（L862-L865）
- 跨文件依赖：`identify_boundary_nodes`

## 省略项

无。所有函数与 struct 均独立列出。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试用途的 `@info` / `@warn`（参数验证 `@assert` 不计入）。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L355 | `R = [0.0 1.0; 1.0 0.0]`（L < 1e-15 时退化旋转矩阵，对应注释 "退化情况"） | 退化几何下旋转矩阵非物理，但仅在零长度单元（应不存在于合法网格）触发；非数值占位 |
| L658 | `t_vec = [1.0, 0.0]; n_vec = [0.0, 1.0]; R = [0.0 1.0; 1.0 0.0]`（L < 1e-15 时的退化几何回退，在 `build_czm_cache`） | 同 L355：零长度 cohesive 单元的退化方向，物理意义可疑 |
| L823 | `penalty = dmax > 0 ? 1e6 * dmax : 1e12`（矩阵全零时回退固定罚值 1e12） | 注释（L820-L821）明示"回退旧值"，无量纲体系下 1e12 可能主导条件数；正常路径不会触发 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L80 | `if (m1 == :PE && m2 == :PCC) \|\| (m1 == :PCC && m2 == :PE)` 紧跟 `elseif (m1 == :NE && m2 == :NCC) \|\| (m1 == :NCC && m2 == :NE)`（4 个 `&&` 条件 + 2 个 `\|\|`，跨 L80-L86） | 抽出 `classify_interface(m1, m2) -> Union{Nothing,Symbol}` helper，用 `Set((m1,m2)) ∈ [Set([:PE,:PCC]), Set([:NE,:NCC])]` 表驱动 |
| L285 | `if size(K_coh, 1) != ndof \|\| nnz(K_coh) == 0`（首次构建稀疏结构判据，2 个条件） | 简单判据，可保留；或显式 `cache.K_coh_initialized` 标志 |
| L721 | `cache === nothing \|\| !cache.valid \|\| cache.czm_mesh_id != objectid(czm_mesh) \|\| cache.param_cache_id != param_cache.id \|\| cache.fix_inner != fix_inner`（5 个 `\|\|` 条件，L721-L724） | 抽出 `is_cache_valid(cache, czm_mesh, param_cache, fix_inner) -> Bool` 谓词函数，提高可读性 |
