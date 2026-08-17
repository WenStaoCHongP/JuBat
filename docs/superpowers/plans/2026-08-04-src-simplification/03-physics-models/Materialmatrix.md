# Materialmatrix.jl Audit Plan

**Status:** ✅ Completed（审计保留） | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。材料矩阵（CZM 参数缓存构造），D-cluster 未命中。

## 现状（427 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `compute_czm_effective_params` | — | 入口 |
| `get_active_elements` | 382 | 活跃单元筛选 |
| 其他 | — | 极片模量加权、CZM 参数 |

**调用者**：`CouplingState.jl:302`（`compute_czm_params_per_interface`）、`czm.jl`

## Audit

- [ ] 通读 427 行
- [ ] grep `try/catch`/`@warn`/`TODO`（spec §7 扫描时未发现）
- [ ] baseline 记录

## Result
无修改。
