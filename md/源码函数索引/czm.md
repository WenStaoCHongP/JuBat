# Czm.jl

- **源文件**: `src/Czm.jl`
- **行数**: 629 行
- **函数/struct 计数**: 1 个 struct + 8 个独立函数
- **职责**: CZM 损伤历史、材料模量映射、内聚力/体单元/热化学载荷装配、装配缓存及完整耦合残差。
- **相关技术文档**: `md/06_内聚力模型_CZM.md`、`md/14_粘性正则化.md`、`md/15_颗粒与极片模量区分.md`

## 数据结构

### `mutable struct DamageState <: AbstractDamageState` — L17-L28

保存等效损伤、粘性有效损伤、历史最大分离、断裂标记和循环累积损伤。默认构造函数将全部状态初始化为零/未断裂。

## 函数清单

### `moduli_of(param, mt) -> (E, ν)` — L57-L65

按材料类型返回模量和泊松比，并用 `scale.E_coat / scale.σ_czm` 将体材料模量统一到 CZM 应力空间。

### `eigenstrain_of(param, mt, dT, Δsn, Δsp) -> ε0` — L67-L89

逐层热-化学本征应变（2026-08-29 α/β 分层化）：`alphaT(mt)·dT + Ω(mt)/3·Δsoc(mt)`，
电极膨胀只作用于本层涂层（NE→Δsn、PE→Δsp），集流体/隔膜仅热应变。取代旧跨层均匀
`α_eff/β_n/β_p` 施加。

### `assemble_czm_system(czm_mesh, u, param_cache; ...)` — L80-L262

装配内聚力刚度、内力、分离和牵引；支持稀疏结构、几何和工作区缓存。

### `assemble_bulk_stiffness(czm_mesh, param_cache)` — L273-L343

按层材料参数装配 Q4 平面应力 bulk 刚度矩阵。

### `assemble_thermal_chemical_load(czm_mesh, param_cache, dT_elem, Δsoc_n_elem, Δsoc_p_elem)` — L352-L416

根据温度与 SOC 变化装配初始应变等效节点载荷；ε₀ 经 `eigenstrain_of` 按单元 `material_type` 分层计算。

### `build_czm_cache(czm_mesh, param_cache; fix_inner=true)` — L440-L523

构建 bulk 刚度、DOF、cohesive 几何、边界 DOF 和装配工作区缓存。

### `ensure_czm_cache(case, czm_mesh, param_cache; fix_inner=true)` — L538-L548

按 mesh identity、参数缓存 id 和 `fix_inner` 状态惰性刷新 `case.czm_cache`。

### `assemble_coupled_system(czm_mesh, u, param_cache; ...)` — L556-L588

组合 bulk 与 cohesive 刚度和内力，不包含热化学载荷。

### `assemble_coupled_system_full(...)` — L596-L629

在耦合装配基础上加入热化学载荷并返回完整残差。

## 跨文件依赖

- `CzmBC.jl`：`identify_bc_nodes_czm`
- `CouplingState.jl`：CZM 缓存、参数缓存和几何工作区类型
- `Materialmatrix.jl`：双线性牵引-分离本构与切线
- `SetMesh.jl`/数值工具：`IntQ4`、`NCweight`

## 省略项

无。全部 struct 与 function 均有独立条目。

### [DEBUG]

无。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|---|---|---|
| L174 | 零长度 cohesive 单元使用固定退化方向 | 非法几何下方向不具物理意义；合法网格不触发 |
| L477 | 缓存构造对零长度单元使用同类固定退化方向 | 与装配路径一致，仅作防御回退 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L104 | 首次构建稀疏结构的双条件判据 | 当前简单可保留，后续可用显式缓存状态 |
| L540-L543 | 缓存有效性由五个条件串联判断 | 后续可抽取 `is_cache_valid` 谓词 |
