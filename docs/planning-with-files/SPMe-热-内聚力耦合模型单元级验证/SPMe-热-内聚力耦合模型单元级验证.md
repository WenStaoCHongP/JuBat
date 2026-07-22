# SPMe-热-内聚力耦合模型单元级验证

## 目标

单元级内聚力网格（8 层 = 8×Q4 + 4×COH2D4）+ 两个验证脚本：

1. 位移 BC，验证双线性牵引-分离律（无电化学/热源）
2. 合成本征应变（ΔT/Δsoc 斜坡），验证界面分离位移（对 1D 解析解）

## 参数

`Jellyroll.jl`（经 `ChooseCell` + `SetCase` 归一化）

## 设计规格（已批准）

完整设计见：

[`docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md`](../../superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md)

### 摘要

| 项 | 选择 |
|----|------|
| 几何 | 平直矩形条带（法向 +y） |
| 网格路径 | 方案 C：复用 `create_czm_mesh` + 硬断言 |
| 生成器 | `src/CzmUnitMesh.jl` → `create_unit_czm_strip` |
| 脚本 1 | `test/unit_czm_bilinear.jl`（Mode I 全曲线 + 卸载 + Mode II） |
| 脚本 2 | `test/unit_czm_eigenstrain.jl`（合成斜坡 + 1D 解析，弹性段） |
| 求解 | 脚本内增量 Newton + `apply_bc_czm(bc_dofs, bc_vals)` |
| 不做 | 真实 SPMe/热步进、改生产本构、环几何 BC |

## 实施计划

[`docs/superpowers/plans/2026-07-23-unit-czm-strip-verification.md`](../../superpowers/plans/2026-07-23-unit-czm-strip-verification.md)

## 状态

- [x] Brainstorming / 设计批准
- [x] 实施计划（writing-plans）
- [x] 实现与验证
