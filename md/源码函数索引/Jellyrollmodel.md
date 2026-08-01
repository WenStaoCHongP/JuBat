# Jellyrollmodel.jl

- **源文件**: `src/Jellyrollmodel.jl`
- **行数**: 691 行
- **函数/struct 计数**: 1 个 struct + 7 个独立函数
- **职责**: Jellyroll 果冻卷电池几何与网格生成——`JellyrollMesh` 数据结构、`jellyroll_collector_seed_mesh` 阿基米德螺旋 collector-seeded 网格、`jellyroll_element_properties` 层面积权重、`jellyroll_tab_node_indices` 极耳节点识别、`setup_thermal2D_mesh` 热网格装配、`build_czm_submesh` 径向 8 层 CZM 子网格
- **相关技术文档**: `md/02_几何与网格.md`

## 数据结构

### `struct JellyrollMesh` — L1-L16

Jellyroll 几何与拓扑的容器，聚合粗热网格（未合并 / 合并）、界面节点对、单元分层信息、极耳节点、CZM 子网格。

- 字段：`thermal2D`（未合并 Mesh）、`thermal2D_merged`（合并重合节点后的 Mesh）、`merge_map`（原节点 → 新节点编号）、`interface_pairs`（界面节点对 `(n_out, n_in)`）、`czm_element_map`（热单元 → CZM 单元 id 列表）、`element_layer`（每单元卷绕层号）、`is_inner_layer`（是否非最外层）、`inner_nodes` / `outer_nodes`（螺旋边界节点）、`pos_tab_nodes` / `neg_tab_nodes`（极耳节点）、`ne` / `nnode`、`czm_submesh`（v3 新增，可为 `nothing`）

## 函数清单

### `jellyroll_collector_seed_mesh(param; nθ, gsorder, phase, tol, nθ_czm) -> JellyrollMesh` — L18-L201

构造 collector-seeded 双螺旋 Q4 网格（内/外螺旋各一条，跨 `[Rin, Rout]`）。

- L21-L31：阿基米德螺旋参数 `a = Rin`、`b = layer/(2π)`；`s_in = 0`、`s_out = layer`（内外螺旋偏移）；`theta0` / `theta1` 由半径边界反解，无有效范围时抛 `error`
- L33-L48：等角度采样——优先用相位对齐的网格点 `phase .+ (k0:k1).*dtheta`，若对齐失败则退化为 `range(theta0, theta1; length=n_theta_eff+1)`
- L50-L67：节点坐标 `r = a + b·θ + s_offset`，x/y 由极坐标转换；前 `n_theta_actual+1` 个为内螺旋，后 `n_theta_actual+1` 个为外螺旋
- L69-L77：单元连接（每段 4 节点 Q4，跨内外螺旋）
- L83-L98：界面节点对识别——双层循环对比坐标 `abs(Δx) < tol && abs(Δy) < tol`，按 `atan(y, x)` 排序
- L101-L116：单元分层 `layer = floor((r_c - Rin)/layer) + 1`；`is_inner_layer` 用外边半径 `< Rout - 0.1·layer` 判定
- L118-L146：构造 `czm_element_map`——建立 `node_to_elem` 反向映射，遍历相邻界面节点对的关联单元
- L148-L149：调 `jellyroll_tab_node_indices` 识别极耳节点
- L151-L193：节点合并——`merged_to` 标记首次出现的节点，重合节点（`dx²+dy² < tol²`）映射到首个；构造 `old_to_new` 重编号、重写单元连接、调 `GetGS` 重算高斯点
- L196-L197：若 `nθ_czm` 非空，调 `build_czm_submesh` 构造细化子网格
- L199-L200：组装 `JellyrollMesh` 并返回
- 跨文件依赖：`GetGS`、`Mesh`、`CzmSubmesh`（`src/SetMesh.jl`）、`jellyroll_tab_node_indices`、`build_czm_submesh`（本文件）

### `jellyroll_element_properties(mesh, param) -> (areas, layer_weights)` — L225-L324

基于螺旋扇形几何计算每单元面积和 5 材料（NE/SP/PE/PCC/NCC）面积权重。

- L229-L234：面积由高斯积分累加 `Σ weight·detJ`
- L240-L249：层序定义（从内到外 8 层，PE→PCC→PE→SP→NE→NCC→NE→SP）+ 各层厚度
- L254-L321：单元循环——计算 4 节点极坐标，取内边半径 `r_in = 0.5·(r1+r4)`，处理角度周期性（`while` 循环把 `dtheta` 归一到 ±π），按扇形面积公式 `A = 0.5·(r_out²-r_in²)·dtheta` 累加各材料，归一化得权重
- L313：总面积为 0 时抛 `error`
- 跨文件依赖：无（纯几何计算）

### `edge_boundary(mesh, nidx, param; which, theta_range, tol) -> Bool` — L346-L368

基于螺旋方程的精确边界节点识别。

- L347-L353：取螺旋参数；`which` 决定偏移 `s_offset`（`:inner` → 0，`:outer` → `layer`，其他抛 `error`）
- L358-L359：反解累计角度 `theta_cum = (r - a - offset)/b`，超出 `[theta_min, theta_max]` 返回 `false`
- L362-L367：理论坐标 `(r_theo·cos, r_theo·sin)` 与实际节点距离 ≤ `tol` 即为边界节点
- 跨文件依赖：无

### `jellyroll_element_centers(mesh) -> Matrix{Float64}` — L379-L382

计算每个 Q4 单元的几何中心 `(ne × 2)`，单行 comprehension 实现。

- 跨文件依赖：无

### `jellyroll_tab_node_indices(mesh, param) -> (pos_idx, neg_idx)` — L396-L499

识别受极耳影响的节点索引（正极耳在内螺旋、负极耳在外螺旋）。

- L397：`@assert mesh.dimension == 2`
- L407-L408：预计算每节点的累计角度 `theta_cum_in` / `theta_cum_out`
- L412-L452：正极耳循环——对每个 `theta_pos`，用 `while` 把 `theta0` 平移到有效范围；二分搜索（L424-L440）求解弧长 `tw` 对应的 `delta_theta`（弧长公式 `F(u) = (u·√(u²+b²) + b²·asinh(u/b))/(2b)`）；筛选满足 `Rin ≤ r ≤ Rout` 且 `theta_start ≤ theta_cum ≤ theta_end` 的节点
- L456-L496：负极耳循环——结构同正极耳，`theta_range = [theta0 - delta_theta, theta0]`
- L498：`unique` 去重后返回
- 跨文件依赖：无

### `setup_thermal2D_mesh(case, mesh_data; use_merged) -> case_new` — L534-L584

将 `JellyrollMesh` 装配到 `case`，返回新的 `case_new`（deepcopy 保持原 case 不变）。

- L535：`deepcopy(case)` 保持原 case 不变
- L537-L549：[v2 修订 2026-07-21] 界面热阻暂禁用，强制 `use_merged = true`（原代码根据 CZM 启用状态自动决定，现注释保留）
- L551-L557：`use_merged=true` 取合并网格 + 空 `interface_pairs`；否则取未合并网格 + 真实 `interface_pairs`
- L561：调 `jellyroll_element_properties` 计算层权重
- L564-L565：调 `identify_boundary_nodes` + `compute_boundary_edge_cache` 预计算边界缓存
- L567-L576：组装 `MeshGeometry`（element_layer / is_inner_layer / layer_weights / interface_pairs / czm_element_map / inner_nodes / outer_nodes / boundary_cache）
- 跨文件依赖：`jellyroll_element_properties`（本文件）、`identify_boundary_nodes`、`compute_boundary_edge_cache`、`MeshGeometry`（外部）

### `build_czm_submesh(param, thermal2D_merged, thermal2D; nθ_czm, gsorder, nθ_thermal, phase) -> CzmSubmesh` — L598-L690

内部辅助函数（不导出）：构造径向 8 层分层 Q4 子网格 + O(1) 解析式 `thermal_elem_map`。

- L599-L601：`thermal2D_merged` 保留在签名中（与 spec §4.1 一致），本函数仅读取 `thermal2D`；预留给未来扩展
- L604-L617：螺旋几何 + 8 层厚度（按层序 PE→PCC→PE→SP→NE→NCC→NE→SP），`@assert sum(thicknesses) ≈ s_total`
- L619-L623：`theta0` / `theta1` 范围裁剪（注释说明省略 `min` 的等价性）
- L625-L634：`n_turns` 用 `ceil` 保证完整覆盖；`n_segments = n_turns · nθ_czm`
- L636-L652：节点生成——`(n_layers+1)` 条螺旋 × `(n_segments+1)` 点
- L654-L685：单元循环——连接 4 节点 Q4，赋材料类型、卷绕圈号 `floor((r_center-a)/s_total)+1`、热单元映射 `clamp(floor((theta_spiral-theta0)/dtheta_thermal)+1, 1, n_thermal)`
- L687-L689：调 `GetGS` 构造高斯点，封装为 `CzmSubmesh`
- 跨文件依赖：`GetGS`、`Mesh`、`CzmSubmesh`（`src/SetMesh.jl`）

## 省略项

无。所有函数与 struct 均独立列出。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L537 | `# ============== [v2 修订 2026-07-21] 界面热阻暂禁用（spec §2.4，先验证 CZM 本构）==============` | 历史决策注释，标记当前禁用界面热阻模型的状态，便于回溯 |
| L542 | `@debug "Auto-selecting thermal mesh" czm_enabled=case_new.opt.czm_enabled use_merged=use_merged`（被注释的原代码） | 结构化调试日志（被注释保留），非临时 `println` |
| L548 | `@debug "Auto-selecting thermal mesh (v2: forced merged)" czm_enabled=case_new.opt.czm_enabled use_merged=use_merged` | 结构化 `@debug`（非 `@info`），默认不输出；标记 v2 强制合并的决策点 |
| L581 | `@debug "Thermal2D mesh setup" ne=ne nnode=nnode n_interface_pairs=n_pairs use_merged=use_merged` | 结构化 `@debug`，非临时调试；网格装配完成后的诊断信息 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L537 | `# 界面热阻暂禁用（spec §2.4，先验证 CZM 本构）` + L547 `use_merged = true   # v2 修订：强制合并网格` | 注释明示"暂禁用"和"v2 修订"，强制使用合并网格绕过界面热阻模型；后续若启用 CZM-热阻耦合需恢复自动选择逻辑 |
| L599 | `thermal2D_merged` 参数保留在签名中但本函数不使用（注释"预留给未来扩展"，跨 L599-L601） | 未使用参数，可能误导调用者认为函数依赖合并网格；spec §4.1 一致性需要 |
| L620 | `# theta1 第二分支 (Rout - a) / b ... 因此此处省略 min 与父函数等价。如未来扩展层序，请恢复 min 形式`（跨 L620-L622） | 显式省略 `min` 依赖于当前层序假设；层序变更时需手动恢复，存在维护风险 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L93 | `if abs(x_out - x_in) < tol && abs(y_out - y_in) < tol`（界面节点对识别的双条件 `&&`） | 仅 2 个条件，未达 ≥3 阈值；可保留。若需简化可用 `hypot(dx, dy) < tol` 单条件替代 |
| L417 | 嵌套 `while` × 2 + `for _ in 1:100` + `for _ in 1:80` 二分搜索（极耳 `delta_theta` 求解，跨 L417-L441 共 25 行） | 算法复杂但非条件链；可抽出 `solve_delta_theta(theta0, tw, a, b) -> delta_theta` 独立函数降低嵌套层级 |
| L448 | `if (Rin - 1e-8 <= r <= Rout + 1e-8) && (theta_start <= theta_cum <= theta_end)`（4 个比较 + 2 个 `&&`，链式比较） | 抽出 `in_annulus(r, Rin, Rout, tol)` 和 `in_angular_window(theta, lo, hi)` 谓词函数提高可读性 |
| L492 | `if (Rin - 1e-8 <= r <= Rout + 1e-8) && (theta_start <= theta_cum <= theta_end)`（同 L448 结构，负极耳分支） | 同 L448 建议；两处重复条件应共用谓词函数 |
