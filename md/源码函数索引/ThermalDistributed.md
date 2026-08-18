# ThermalDistributed.jl

- **源文件**: `src/ThermalDistributed.jl`
- **行数**: 418 行
- **函数/struct 计数**: 8 个独立函数；0 个 struct
- **职责**: 二维分布式热模型主体，组装各向异性导热的质量/刚度/载荷矩阵，提供 2D 平面、2D 圆环（polar）两套网格接口；含对流/侧面/极耳冷却边界条件与分层热源（NE/SP/PE/PCC/NCC）计算
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/03_边界条件.md`、`md/07_界面热阻模型.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ThermalDistributed2D(case::Case, variables::Dict{String,<:Any})` — L1-L47

二维分布式热模型矩阵装配主入口，返回 `(MT, KT, FT)`。

- 通过 `case.geometry.layer_weights` 读取预计算层权重 `fks`（`setup_thermal2D_mesh` 阶段一次性算好，不再现场调 `jellyroll_element_properties`）
- 质量矩阵：`thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)` → `Assemble`
- 刚度矩阵：调 `thermal_anisotropic_conductivity_2d` 得到 `(k_xx, k_xy, k_yy)`，加负号与电化学约定统一；分别装配 `KT_xx/xy/yx/yy` 后求和
- 载荷向量：`q_elem = variables["heat_source_fields"]` → `Assemble1D`
- 跨文件依赖：`Assemble`、`Assemble1D`、`thermal_capacity_weights_2d`、`thermal_anisotropic_conductivity_2d`

### `apply_convection_bc(K, F, mesh, is_outer, case; edge_cache)` — L49-L86

为矩阵 `(K, F)` 叠加 Robin 对流边界条件（外层边界），返回修改后的副本。

- `Bi = scale.h × lambda_r`（无量纲 Biot 数，L50）；`Bi == 0` 时直接返回（L51-L53）
- 使用 1D 2-点 Gauss 积分沿边界边累加 `ke_ij = -Bi·N_i·N_j·J`、`fe_i = Bi·T_amb·N_i·J`（`wt = Bi·w·J`，L74）
- 优先使用预计算 `edge_cache`（`case.geometry.boundary_edges`）；为 `nothing` 时调 `identify_boundary_nodes` + `compute_boundary_edge_cache` 现场计算（L60-L65）
- 跨文件依赖：`identify_boundary_nodes`、`compute_boundary_edge_cache`

### `apply_cool_method(K, F, mesh, case)` — L88-L158

按 `case.opt.cool_method`（"none"/"surface"/"tab"）返回应用冷却边界后的 `(K, F)` 副本。

- "none"：返回原矩阵（L90-L91）
- "surface"（L92-L117）：**经典 Biot 数 `Bi = scale.h`（不乘 lambda_r）**；每高斯点 `wt = conv_factor·wJ`，对节点对 `(i,j)` 累加 `-wt·N_i·N_j`，源项 `+wt·T_amb·N_i`；`conv_factor = 2·Bi/width`
- "tab"（L118-L154）：先 `jellyroll_tab_node_indices` 取极耳节点，按节点弧长加权 `coeff = tab.h·tab.area·w/width`，单点修改 `K[n,n]`、`F[n]`
- 跨文件依赖：`jellyroll_tab_node_indices`

> 早期的原位变体 `apply_convection_bc!` / `apply_cool_method!` 已随 b4c0cde 重构删除，
> 统一为本非原位接口。

### `ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)` — L160-L196

热模型边界条件包装器：单次 `copy`（L162-L164）后依次调用 `apply_convection_bc`（传入 `case.geometry.boundary_edges` 缓存）与 `apply_cool_method`（L192-L195）。

- 见 [PLACEHOLDER] L166-L190：CZM 双向界面热阻耦合代码块已被注释禁用（spec §2.4，v2 修订 2026-07-21），原因：CZM 损伤场与温度场双向耦合会让参数空间与收敛行为同时变化

### `ThermalDistributed2D_Ring(case::Case, variables::Dict{String,Any})` — L198-L248

圆环（polar）网格版本的二维热模型矩阵装配，返回 `(MT, KT, FT)`。

- 体积热容 `rho_c_nd = heat_Q / volume`（无量纲）
- 在每个高斯点用 `atan(gy, gx)` 旋转得到 `(dNdr, dNdtheta)`，分别装配径向 `KT_r`、切向 `KT_t` 刚度
- 不使用层权重（与 `ThermalDistributed2D` 不同）

### `ThermalRing2D_BC(KT, FT, case::Case, outer_nodes, t::Float64)` — L249-L256

圆环网格边界条件包装：基于 `outer_nodes` 构造 `is_outer` 掩码后委托 `apply_convection_bc`（5 参签名，`outer_nodes` 显式传入）。

### `compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme)` — L258-L388

分层热源计算，逐单元计算并写入 `variables`。

- 通过 `fks[e,1..5]` 分配各层权重：NE / SP / PE / PCC / NCC（L355-L366）
- 计算分量（每单元 e，L333-L352）：
  - `Q_rxn_NE/PE = as·|j|·|η|`（反应热）
  - `Q_rev_NE/PE = as·j·T·dUdT`（可逆热/熵热）
  - `Q_ohm_s_NE/PE = I²/(3·σ_eff)`（固相欧姆热）
  - `Q_ohm_e_NE/PE = I²/(3·κ_av)`（电解液相欧姆热）
  - `Q_SP = I²/κ_sp_av`（隔膜欧姆热）
  - `Q_PCC/NCC = I²/(3·σ_pcc/ncc)`（集流体欧姆热；`sigma_pcc/ncc = param.PCC/NCC.sig` 直取，L283-L284，无下限替换）
  - `Q_EL = -I·2·T·(1-t⁺)·(csp_av - csn_av)/ce0`（电解液修正项，L352）
- 总热源 `q_total[e]` 12 项求和（L369）
- 末尾统一缩放：`q .* scale.L³ / cell.volume` 写入 variables 各键（L373-L385）
- 跨文件依赖：`IntV`、`param.EL.kappa`

### `compute_heat_sources_with_czm(case, variables, variables_elems, I_e, T_e, areas, czm_mesh, mesh_data)` — L390-L418

CZM 损伤感知版本的热源计算：先调 `compute_heat_sources`，再将非活跃单元（被 CZM 失效隔离）的热源置零。

- 通过 `get_active_elements(czm_mesh, case.geometry)` 取活跃单元（L394，直接使用 `case.geometry`，无 hasfield 守卫）
- 更新 `variables["heat_source_fields"]`、`variables["active_elements"]`、`variables["total heat source"]`（仅活跃单元求和，L412-L416）
- 签名中的 `mesh_data` 参数当前未在函数体内使用（保留签名兼容）
- 跨文件依赖：`compute_heat_sources`、`get_active_elements`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L166-L190 | 注释块 "界面热阻暂禁用（spec §2.4，v2 修订 2026-07-21）" | CZM-温度场双向耦合被显式禁用，恢复需取消注释并同时恢复 `setup_thermal2D_mesh` 的 `use_merged` 自动逻辑。当前热模型与 CZM 解耦，相关接口未启用 |

> 2026-08-18 复核：早期 `sigma_pcc/ncc = max(param.*.sig, 1e-12)` 的静默下限
> 已移除（现为直取，L283-L284），不再计入 PLACEHOLDER。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L128-L141 | `apply_cool_method` 的 "tab" 分支 `if n_nodes == 1 / elseif i == 1 / elseif i == n_nodes`（嵌套 3+ 层弧长加权） | 抽出独立的 `compute_arc_lengths(tab_nodes, mesh)` 工具函数，主流程只调用一次 |
