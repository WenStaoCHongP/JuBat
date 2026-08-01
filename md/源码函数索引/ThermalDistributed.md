# ThermalDistributed.jl

- **源文件**: `src/ThermalDistributed.jl`
- **行数**: 557 行
- **函数/struct 计数**: 10 个独立函数；0 个 struct
- **职责**: 二维分布式热模型主体，组装各向异性导热的质量/刚度/载荷矩阵，提供 2D 平面、2D 圆环（polar）两套网格接口；含对流/侧面/极耳冷却边界条件与分层热源（NE/SP/PE/PCC/NCC）计算
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/03_边界条件.md`、`md/07_界面热阻模型.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ThermalDistributed2D(case::Case, variables::Dict{String,<:Any})` — L1-L47

二维分布式热模型矩阵装配主入口，返回 `(MT, KT, FT)`。

- 通过 `case.geometry.layer_weights`（或回退到 `jellyroll_element_properties`）获取层权重 `fks`
- 质量矩阵：`thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)` → `Assemble`
- 刚度矩阵：调 `thermal_anisotropic_conductivity_2d` 得到 `(k_xx, k_xy, k_yy)`，加负号与电化学约定统一；分别装配 `KT_xx/xy/yx/yy` 后求和
- 载荷向量：`q_elem = variables["heat_source_fields"]` → `Assemble1D`
- 跨文件依赖：`Assemble`、`Assemble1D`、`thermal_capacity_weights_2d`、`thermal_anisotropic_conductivity_2d`、`jellyroll_element_properties`

### `apply_convection_bc(KT, FT, mesh, is_outer, case; edge_cache)` — L49-L96

为矩阵 `(KT, FT)` 叠加 Robin 对流边界条件（外层边界），返回修改后的副本 `(K, F)`。

- `Bi = scale.h × lambda_r`；`Bi == 0` 时直接返回原矩阵
- 使用 1D 2-点 Gauss 积分沿边界边累加 `ke_ij = -Bi·N_i·N_j·J`、`fe_i = Bi·T_amb·N_i·J`
- `edge_cache === nothing` 时调 `identify_boundary_nodes` + `compute_boundary_edge_cache` 现场计算
- 跨文件依赖：`identify_boundary_nodes`、`compute_boundary_edge_cache`

### `apply_cool_method(KT, FT, mesh, case)` — L99-L174

按 `case.opt.cool_method`（"none"/"surface"/"tab"）返回应用冷却边界后的 `(K, F)` 副本。

- "none"：返回原矩阵副本
- "surface"：每高斯点 `wt = conv_factor·wJ`，对节点对 `(i,j)` 累加 `-wt·N_i·N_j`，源项 `+wt·T_amb·N_i`；`conv_factor = 2·Bi/width`
- "tab"：先 `jellyroll_tab_node_indices` 取极耳节点，按节点弧长加权 `coeff = tab.h·tab.area·w/width`，单点修改 `K[n,n]`、`F[n]`
- 跨文件依赖：`jellyroll_tab_node_indices`

### `apply_convection_bc!(K, F, case; edge_cache)` — L178-L212

`apply_convection_bc` 的原位变体（消除多余 `copy`）；`edge_cache === nothing` 时回退调用非原位版本。

### `apply_cool_method!(K, F, mesh, case)` — L214-L284

`apply_cool_method` 的原位变体，逻辑与原函数逐行一致。

### `ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)` — L286-L323

热模型边界条件包装器：仅做一次 copy 后依次调用原位变体。

- 见 [PLACEHOLDER] L292-L316：CZM 双向界面热阻耦合代码块已被注释禁用（spec §2.4），原因：CZM 损伤场与温度场双向耦合会让参数空间与收敛行为同时变化
- 调 `apply_convection_bc!`（含 edge_cache）→ `apply_cool_method!`

### `ThermalDistributed2D_Ring(case::Case, variables::Dict{String,Any})` — L325-L374

圆环（polar）网格版本的二维热模型矩阵装配，返回 `(MT, KT, FT)`。

- 体积热容 `rho_c_nd = heat_Q / volume`（无量纲）
- 在每个高斯点用 `atan(gy, gx)` 旋转得到 `(dNdr, dNdtheta)`，分别装配径向 `KT_r`、切向 `KT_t` 刚度
- 不使用层权重（与 `ThermalDistributed2D` 不同）

### `ThermalRing2D_BC(KT, FT, case::Case, outer_nodes, t::Float64)` — L376-L383

圆环网格边界条件包装：构造 `is_outer` 掩码后委托 `apply_convection_bc`。

### `compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme)` — L385-L516

分层热源计算（CLAUDE.md §6.3），逐单元计算并写入 `variables`。

- 通过 `fks[e,1..5]` 分配各层权重：NE / SP / PE / PCC / NCC
- 计算分量（每单元 e）：
  - `Q_rxn_NE/PE = as·|j|·|η|`（反应热）
  - `Q_rev_NE/PE = as·j·T·dUdT`（可逆热/熵热）
  - `Q_ohm_s_NE/PE = I²/(3·σ_eff)`（固相欧姆热）
  - `Q_ohm_e_NE/PE = I²/(3·κ_av)`（电解液相欧姆热）
  - `Q_SP = I²/κ_sp_av`（隔膜欧姆热）
  - `Q_PCC/NCC = I²/(3·σ_pcc/ncc)`（集流体欧姆热）
  - `Q_EL = -I·2·T·(1-t⁺)·(csp_av - csn_av)/ce0`（电解液修正项）
- 总热源 `q_total[e]` 见 L497（12 项求和）
- 末尾统一缩放：`q .* scale.L³ / cell.volume` 写入 variables 各键
- 跨文件依赖：`IntV`、`jellyroll_element_properties`、`param.EL.kappa`

### `compute_heat_sources_with_czm(case, variables, variables_elems, I_e, T_e, areas, czm_mesh, mesh_data)` — L518-L556

CZM 损伤感知版本的热源计算：先调 `compute_heat_sources`，再将非活跃单元（被 CZM 失效隔离）的热源置零。

- 通过 `get_active_elements(czm_mesh, geom)` 取活跃单元；`geom` 缺失时全部单元活跃
- 更新 `variables["heat_source_fields"]`、`variables["active_elements"]`、`variables["total heat source"]`
- 跨文件依赖：`compute_heat_sources`、`get_active_elements`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L292 | 注释块 "界面热阻暂禁用（spec §2.4）"（覆盖 L292-L316） | CZM-温度场双向耦合被显式禁用，恢复需取消注释并同时恢复 `setup_thermal2D_mesh` 的 `use_merged` 自动逻辑。当前热模型与 CZM 解耦，相关接口未启用 |
| L411 | `sigma_pcc = max(param.PCC.sig, 1e-12)` | 1e-12 是物理下限，防止除零；若参数集本身 σ 为 0 会被静默替换，可能掩盖参数缺失 |
| L412 | `sigma_ncc = max(param.NCC.sig, 1e-12)` | 同 L411 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L142 | `apply_cool_method` 的 "tab" 分支 `if n_nodes == 1 / elseif i == 1 / elseif i == n_nodes`（嵌套 3+ 层，含 L142-L154） | 抽出独立的 `compute_arc_lengths(tab_nodes, mesh)` 工具函数，主流程只调用一次 |
| L254 | `apply_cool_method!` 中相同的 tab 分支嵌套逻辑（L254-L266） | 与 L142 重复，可共用同一工具函数；同时减少原位/非原位版本间的代码重复 |
| L524 | `if geom !== nothing && hasfield(typeof(geom), :czm_element_map)`（≥2 个 `!== nothing`/`hasfield` 链） | 封装为 `has_czm_element_map(geom)` 谓词函数，提高可读性 |
