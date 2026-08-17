# CzmMesh.jl

- **源文件**: `src/CzmMesh.jl`
- **行数**: 182 行
- **函数/struct 计数**: 1 个 struct + 1 个独立函数
- **职责**: CZM 界面单元拓扑与网格构造；识别 PE-PCC / NE-NCC 径向界面、复制界面节点、重写外层 bulk 连接并建立热-CZM 映射。
- **相关技术文档**: `md/02_几何与网格.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

### `mutable struct CohesiveElement <: AbstractCohesiveElement` — L3-L12

单个 COH2D4 内聚力单元的拓扑与几何信息。

- `nodes` 为 `[n_lo, n_hi, n_hi_copy, n_lo_copy]`。
- `interface_type` 为 `:PE_PCC` 或 `:NE_NCC`。
- 两种 `interface_type` 是参数类型；每个 8 层重复单元实际有四个箔–涂层面，完整网格 `n_cohesive = 4 * (length(theta)-1)`。
- `host_outer_elem` / `host_inner_elem` 记录相邻 bulk 单元。

## 函数清单

### `create_czm_mesh(czm_submesh, thermal_mesh, param) -> CohesiveMesh` — L33-L182

从细化的 `CzmSubmesh` 构造 `CohesiveMesh`。

- 建立共边到单元对的映射并识别目标材料界面。
- 每个周向分段识别 PE–PCC/PCC–PE、NE–NCC/NCC–NE 四个真实面，再归入两种 `interface_type`。
- 按径向质心判定内外层。
- 复制界面节点并重写外层 bulk 单元连接，确保可产生分离位移。
- 构造 `cohesive_to_thermal` 与 `thermal_to_czm` 映射。
- 初始化每个 cohesive 单元的 `DamageState` 并执行拓扑自检。

## 跨文件依赖

- `Czm.jl`：`DamageState`
- `CouplingState.jl`：`CohesiveMesh`、`CzmSubmesh`
- `SetMesh.jl`：`Mesh`
- `CouplingState.jl`/映射实现：`build_thermal_to_czm_interp`

## 省略项

无。全部 struct 与 function 均有独立条目。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L55-L61 | 两组材料界面分类使用多重 `&&` / `||` | 后续可抽取表驱动的 `classify_interface`，但需保持当前接口判定语义 |
