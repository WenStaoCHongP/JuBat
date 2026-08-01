# CzmUnitMesh.jl

- **源文件**: `src/CzmUnitMesh.jl`
- **行数**: 108 行
- **函数/struct 计数**: 1 个独立函数
- **职责**: 单元条带（unit strip）CZM 网格生成器，用于验证脚本（8×Q4 + 4×COH2D4，平直几何），底层复用 `create_czm_mesh`
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`test/unit_czm_*.jl`（`unit_czm_bilinear.jl`、`unit_czm_eigenstrain.jl`、`unit_czm_newton.jl`、`unit_czm_strip_mesh.jl`）

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `create_unit_czm_strip(param; width, y0, gsorder) -> (czm_mesh, meta)` — L7-L108

构造平直 8 层单元条带 + 经 `create_czm_mesh` 得到 4×COH2D4 内聚力网格。

- **层序**（L8）：自下而上 `PE → PCC → PE → SP → NE → NCC → NE → SP`
- **几何**（L17-L33）：底边 `y=y0>0`；节点 18 个（9 行 × 2 列），节点 id 沿 x 递增
- **Q4 连接**（L36-L43）：每单元 `[bl, tl, tr, br]`（逆时针，inner=bottom）
- **哑热网格**（L52-L61）：单 Q4 覆盖条带 bbox（含 `pad = 1e-6·max(W,H)` 余量），供 `build_thermal_to_czm_interp` 使用
- **硬断言**（方案 C，L65-L80）：
  - `czm_mesh.n_cohesive == 4`（L66）
  - PE_PCC × 2 + NE_NCC × 2（L68-L69）
  - 每单元 4 节点不重复、副本坐标一致、外层 bulk 含副本节点（L70-L79）
- **meta 字段**（L94-L106）：`y_interfaces`、`bottom_nodes`、`top_nodes_after_czm`、`pcc_nodes`、`ncc_nodes`、`cohesive_ids`、`interface_types`、`layer_materials`、`width`、`heights`
- 跨文件依赖：`GetGS`、`Mesh`、`CzmSubmesh`、`create_czm_mesh`

## 省略项

无。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info`；`@assert` 为单元测试用途的硬断言（方案 C），不计入 DEBUG。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L53 | `pad = 1e-6 * max(W, H)`（哑热网格 bbox 余量） | 数值容差，物理上合理；若验证脚本中 `create_unit_czm_strip` 与 thermal interp 的相对尺度异常，pad 可能掩盖边界效应 |
| L83 | `# 顶层（层 8）在 create_czm_mesh 后：外层若被重写则取 bulk_element[8,:] 上边`（注释解释 `top_nodes_after_czm` 计算逻辑） | 注释说明而非占位代码；逻辑正确但需要读者了解 `create_czm_mesh` 副本机制 |
| L88 | `top_nodes_after = sort([Int(n) for n in e8 if abs(czm_mesh.node[n, 2] - ytop) < 1e-14])`（容差 1e-14） | 紧容差，浮点比较可能漏节点；但节点坐标由 `param.*.thickness` 直接求和，无积分误差，1e-14 应足够 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L70 | `for coh in czm_mesh.cohesive_elements` 内含 4 行 `@assert`（L72-L79），其中 `n_lo_c in outer_nodes && n_hi_c in outer_nodes` 与 `!(n_lo in outer_nodes) && !(n_hi in outer_nodes)` 是 4 个 `&&`/`!` 组合 | 抽出 `verify_cohesive_topology(coh, czm_mesh)` 断言块函数，主流程只调用一次 |
