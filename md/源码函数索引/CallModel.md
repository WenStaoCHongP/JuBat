# CallModel.jl

- **源文件**: `src/CallModel.jl`
- **行数**: 269 行
- **函数/struct 计数**: 4 个函数（无独立 struct）
- **职责**: 多物理场装配分发——多 SPMe 模式装配（`CallModel_MultiSPMe`）、单模型/多模型统一入口（`CallModel`）、CZM 损伤到热单元的归约（`map_czm_damage_to_thermal`）、workspace Dict 到轻量 Dict 的提取（`copy_element_results`）
- **相关技术文档**: `md/08_逐单元算法.md`、`md/09_分流求解器.md`、`md/10_参数传递与模块架构.md`

## 数据结构

本文件无独立 struct 定义。`CallModel_MultiSPMe` 返回 `(M, K, F, variables, y_phi)`，其中 `variables` 是 `Dict{String, Union{Array{Float64},Float64}}`，`y_phi` 对 SPMe 为 `Float64[]`。

## 函数清单

### `map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, ne_thermal::Int) -> Vector{Float64}` — L8-L20

按 spec v2 §5.2 把 cohesive 单元损伤 max 归约到粗热单元粒度。

- 初始化 `D_per_thermal = zeros(ne_thermal)`（L9）
- 遍历 cohesive 单元（L10-L18）：读取 `czm_mesh.cohesive_to_thermal[e_coh]` 取对应热单元索引，范围检查（L12），取最大值（L13-L15）
- 无 cohesive 覆盖的热单元保持 0

### `CallModel_MultiSPMe(case, yt, t; jacobi) -> (M, K, F, variables, y_phi)` — L22-L210

多 SPMe 模式的全装配入口。9 个步骤：

1. **验证与状态切片**（L24-L46）：断言 `case.layout` 非空（L24-L25）；提取 `T_nodes`（L40）；用 `@view` 切出每个单元电化学状态 `yt_chem[e]`（L42-L46）
2. **缓存面积**（L48-L56）：复用 `layout.areas`，从节点温度算单元均温 `Te_prev`
3. **分流求解前置**（L58-L79）：归一化总电流 `I_total`（L59）；用单元平均状态生成代表性 `vars_rep`（L71-L76）；保留热场覆盖（L77-L79）
4. **CZM 失效映射**（L81-L99）：D ≥ 0.95 的 cohesive 单元视为断裂，对应热单元加入 `deactivated_elements`（L83-L96）；无 CZM 时为空
5. **渐进式面积损失**（L101-L111）：调 `map_czm_damage_to_thermal`（L104）；超阈值单元的 `debug_coupling` 打印（L106-L110）
6. **分流求解**（L113-L115）：调 `solve_branch_currents(...)` 得到 `I_e`、`Vc`
7. **并行 SPMe 求解**（L117-L138）：预分配线程本地 workspace `ws_pool`（L125）；`Threads.@threads` 调 `SPMe_element`（L129-L137）
8. **装配全局矩阵**（L141-L186）：`M_chem = blockdiag(M_elems...)`（L141）；按 `czm_enabled` 选择 `compute_heat_sources_with_czm` 或 `compute_heat_sources`（L148-L162）；调 `ThermalDistributed2D` + `ThermalDistributed2D_BC` 装配热学矩阵（L181-L185）
9. **合并 variables**（L188-L209）：写入 `cell voltage = Vc`、`temperature = thermal2D_volume_average_temperature(...)`、`thermal2D element voltages`（用于诊断）、四个 timing 字段（L202-L205）；`y_phi = Float64[]`（L207）

### `CallModel(case, yt, t; jacobi) -> (M, K, F, variables, y_phi)` — L211-L238

单模型/多模型统一入口。

- **多 SPMe 委托**（L213-L215）：`case.opt.per_element_spme` 为真时直接 `return CallModel_MultiSPMe(...)`
- **单模型分支**（L218-L230）：按 `case.opt.model` 分发到 `SPM`/`SPMe`/`P2D`/`sP2D`；未实现时 `error`
- **热学 lumped 装配**（L231-L236）：仅当 `thermalmodel == "lumped"` 时追加 1×1 热块；distributed2D 分支不在本函数处理（由 `SPMe`/`SPM` 内部完成）

### `copy_element_results(vars_e) -> Dict` — L247-L269

从 workspace Dict 中提取下游需要的键，返回轻量级独立 Dict。

- 注释明示键分两类（L244-L246）：
  - **计算结果键**（L249-L262）：每次 `=` 赋值创建新对象，引用安全（如 `"cell voltage"`、各种 overpotential / OCP / exchange current density / Gauss 点浓度）
  - **状态提取键**（L263-L268）：通过 `case.index` 原位写入 workspace，注释说"workspace 不再持有引用，直接传递"（如 `"negative particle lithium concentration"`）

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L109 | `println("  [AreaLoss] t=$(t_phys)s \| D_max=... \| 超阈值单元=...")`（在 `case.opt.debug_coupling` 守卫下） | 调试输出：渐进式面积损失触发的单元数与最大损伤值 |
| L164 | `# 保存辅助变量（用于调试）` 注释 | 标注 eta_n_e / eta_p_e / dUdT_e / soc_n / soc_p 等字段为调试用途 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L207 | `y_phi = Float64[]`（多 SPMe 模式下 phi 自由度空数组） | 占位：多 SPMe 架构在 CallModel_MultiSPMe 中不返回额外电位自由度，与 P2D 路径的 y_phi 语义不对称；调用者（Solve.jl L181/L215）需正确处理空数组 |
| L67 | `I_e_prev = hasproperty(case, :I_e_cache) ? case.I_e_cache : nothing`（条件 hasproperty 检查缓存） | 兜底：若 `case` 未提前注入 `I_e_cache` 字段则回退到 nothing，初次调用必然走 nothing 分支；正常但隐式依赖 case 的动态字段 |
| L87 | `if geom !== nothing && hasfield(typeof(geom), :czm_element_map)`（条件检查几何体的 czm_element_map 字段是否存在） | 防御性兜底：若 geometry 类型未定义该字段则跳过整个 deactivated_elements 收集逻辑；意味着 CZM 失效映射仅在 geometry 显式支持时生效 |
| L97 | `deactivated_elements = Int64[]`（跨 L97-L99，无 CZM 或无 czm_mesh 时直接空数组） | 占位：无 CZM 模式下分流求解器收到空 deactivated_elements，行为退化为普通分流 |
| L183 | `t_ratio = 1.0` 然后 `MT = MT .* t_ratio`（注释"统一时间尺度"） | 占位乘子：t_ratio 恒为 1.0，原本为电化学/热时间尺度比，统一后变成无操作；保留是为未来灵活性 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L83 | `if case.opt.czm_enabled && case.czm_mesh !== nothing`（CZM 失效映射入口条件） | 与 L103、L148 重复出现 3 次相同条件 `case.opt.czm_enabled && case.czm_mesh !== nothing`，可提取为本地 `czm_active` 布尔减少重复 |
| L87 | `if geom !== nothing && hasfield(typeof(geom), :czm_element_map); for e in 1:ne; for czm_idx in get(geom.czm_element_map, e, Int64[]); if czm_idx in fractured_czm; push!(...); break; end; end; end; end`（嵌套 4 层：if × 2 + for × 2 + 内 if；跨 L87-L96） | 抽出 `collect_deactivated_elements(geom, ne, fractured_czm)` 函数；当前双层 for + break + 嵌套 if 阅读成本高 |
| L103 | `if case.opt.czm_area_loss_enabled && case.czm_mesh !== nothing && case.czm_mesh.cohesive_to_thermal !== nothing`（3 个 `&&` 链） | 接近 ≥3 阈值但未越界（恰好 3 个）；可考虑抽 `czm_area_loss_ready(case)` helper |
| L161 | `else` 分支调用 `compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas; per_element_spme=true)`（与 L150 的 `compute_heat_sources_with_czm` 二选一） | 双路径不对称：CZM 路径返回完整 variables，非 CZM 路径传 `per_element_spme=true` 关键字，签名不一致；可统一入口 |
