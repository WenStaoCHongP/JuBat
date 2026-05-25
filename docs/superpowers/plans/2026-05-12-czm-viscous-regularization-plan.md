# CZM 粘性正则化实施计划

> **Spec Document:** `docs/superpowers/specs/2026-05-12-czm-viscous-regularization-design.md`
>
> **Context:** 当前 Jellyroll 电-热-CZM 耦合例程在后峰软化区失稳，`D_max` 会增长到约 0.895 后出现残差放大与求解停滞。根因不是 NaN，而是速率无关的双线性内聚律在软化段切线过软，导致全局 Newton / continuation 难以继续跟踪。
>
> **Goal:** 在保持默认行为不变的前提下，引入局部、回滚安全的粘性正则化，使损伤演化和牵引切线在软化区更平滑，从而提高耦合求解的可收敛性。

## Design Principles

1. 正则化要作用在内聚本构/历史变量上，而不是只调外层步长或求解器参数。
2. 默认关闭，确保现有案例结果不变。
3. 粘性历史量只能在非线性求解收敛后提交，失败步必须可回滚。
4. 正则化幅度应能随 `τ_v → 0`（β → 1）退化回当前速率无关模型。
5. 切线刚度中的 `dD/dδ` 项必须显式乘以 β，不能仅替换 D → D_visc。
6. 粘性正则化需配合载荷子步求解器使用；`basic` 模式下 Δs=1.0，正则化效果有限。
7. 先解决单线程、单路径的稳定性，再考虑更广泛的调参或并行。

## Scope

| File                                       | Change | Purpose                                                                    |
| ------------------------------------------ | ------ | -------------------------------------------------------------------------- |
| `src/Option.jl`                          | Modify | 增加 `czm_viscous_enabled` 开关与 `czm_visc_tau` 参数                  |
| `src/SetParams.jl`                       | Modify | τ_v 归一化：`tau_visc = czm_visc_tau / t0`，扩展 `Cohesive` 结构体    |
| `src/czm.jl`                             | Modify | 扩展 `DamageState` 增加 `D_visc` 字段，更新 `clone_damage_states`    |
| `src/Materialmatrix.jl`                  | Modify | `bilinear_traction_state` 使用 D_visc；`bilinear_tangent` 增加 β 因子 |
| `src/CzmSolve.jl`                        | Modify | 计算每子步的 β 和 Δs，传递到本构函数                                     |
| `src/CouplingState.jl`                   | Modify | 传递粘性参数，扩展调试输出（D_eq_max, D_visc_max, β）                     |
| `src/Solve.jl`                           | Verify | 确认 CZM 调用链无需直接修改                                                |
| `src/CycleSolver.jl`                     | Verify | 确认 D_visc 跨周期传递正确                                                 |
| `example/coupled_czm_thermal_example.jl` | Modify | 增加粘性开关与 load_substep 配置，用于回归验证                             |

## Phases

### Phase 1: Formulation Lock-In (已完成)

- [X] 粘性正则化同时作用于牵引力和切线刚度（D_visc 替代 D_eq）
- [X] `τ_v` 归一化方案：物理输入 [s]，归一化 `τ_v* = τ_v / t0`
- [X] `D_eq`（等效损伤）、`D_visc`（粘性有效损伤）、β（松弛因子）的命名和职责已确认
- [X] Δs 定义为载荷分数增量（dimensionless, 0~1），跨三种求解器统一

- **Status:** complete

### Phase 2: Data Model

- [ ] 在 `Option` 中增加 `czm_viscous_enabled::Bool = false` 与 `czm_visc_tau::Float64 = 0.0`
- [ ] 在 `Cohesive` 结构体（`SetParams.jl`）中增加 `tau_visc::Float64 = 0.0`（归一化后的值）
- [ ] 在 `NormaliseParam` 中增加 `param.cohesive.tau_visc = opt.czm_visc_tau / param.scale.t0`
- [ ] 扩展 `DamageState`，增加 `D_visc::Float64` 字段（默认 0.0）
- [ ] 更新 `clone_damage_states` 拷贝 `D_visc`
- [ ] 更新 `DamageState()` 默认构造函数

- **Status:** pending

### Phase 3: Constitutive Update

- [ ] 修改 `bilinear_traction_state`：增加 `beta::Float64=1.0` 参数，牵引力使用 `D_visc` 计算
- [ ] 修改 `bilinear_tangent`：增加 `beta::Float64=1.0` 参数，软化分支 `dD/dδ` 项乘以 `beta`
- [ ] 保持压缩分支和既有双线性规则不变（β 不影响压缩区）
- [ ] 确保 `beta=1.0` 时行为与当前实现完全一致
- [ ] 修改 `assemble_czm_system` 传递 `beta` 到本构函数

- **Status:** pending

### Phase 4: Solver Integration

- [ ] 在 `solve_czm_basic_step` / `newton_raphson_czm` / `solve_czm_arc_length_step` 中计算每子步的 β：
  - `beta = viscous_enabled ? Δs / (tau_visc + Δs) : 1.0`
  - `basic`：Δs = 1.0
  - `load_substep`：Δs = `step_size`
  - `arc_length`：Δs = `delta_lambda`
- [ ] 将 `beta` 传递到 `assemble_coupled_system` → `assemble_czm_system`
- [ ] 在收敛时更新 `D_visc`（松弛公式 + 单调性约束），失败步回滚自动生效
- [ ] 在调试模式下打印 `D_eq_max`、`D_visc_max`、`β` 和收敛结果

- **Status:** pending

### Phase 5: Verification

- [ ] 关闭粘性（`czm_viscous_enabled = false`），确认结果与基线一致
- [ ] 退化测试：`czm_visc_tau = 0`（β=1），确认结果与关闭时一致
- [ ] 使用 `load_substep` 求解器 + 粘性正则化，观察 `D_max ≈ 0.895` 区间是否仍停滞
- [ ] 检查残差峰值、收敛标志、D_visc 单调性
- [ ] 确认 `basic` 模式下粘性正则化效果有限（β ≈ 1）
- [ ] 确认失败步没有推进损伤历史
- [ ] 记录结果到 `docs/planning-with-files/内聚力病态问题解决/progress.md`

- **Status:** pending

## Success Criteria

- 默认关闭时，结果与当前基线完全一致。
- `τ_v = 0` 时，结果与关闭时完全一致（退化验证）。
- 开启粘性 + `load_substep` 后，耦合例程能跨过当前软化停滞区。
- 失败步不会推进损伤历史（D_visc 不变），也不会污染下一次求解。
- `basic` 模式下结果接近无正则化（Δs=1.0 → β≈1）。

## Risks

| Risk                       | Impact                          | Mitigation                                 |
| -------------------------- | ------------------------------- | ------------------------------------------ |
| `τ_v` 过大导致过度增韧  | 物理响应偏硬，损伤推迟          | 从小值（10s）开始测试，观察 D_visc 滞后量  |
| `τ_v` 过小无改善        | 耗时但无效果                    | β 接近 1 时自动退化，不影响正确性         |
| 切线中遗漏 β 因子         | 残差与切线不一致，Newton 不收敛 | Phase 3 明确要求所有 `dD/dδ` 项乘以 β  |
| `basic` 模式下正则化无效 | 用户困惑                        | 文档说明需配合 load_substep 使用           |
| 回滚不完整                 | D_visc 在失败步被半提交         | clone_damage_states 拷贝所有字段，自动安全 |

## Key Formulation Reference

```
Δs 定义:  载荷分数增量 (dimensionless, 0~1)
τ_v*:     τ_v / t0 (无量纲松弛时间)
β:        Δs / (τ_v* + Δs)
D_visc:   D_visc_prev + β * (D_eq - D_visc_prev)
单调性:    D_visc = max(D_visc_prev, D_visc)
牵引力:    T = (1 - D_visc) * K * δ
切线:      dT/dδ = (1 - D_visc) * K - K * δ * β * dD_eq/dδ
退化:      τ_v=0 → β=1 → D_visc=D_eq, 切线=原始
```

## Notes

- 本计划不把"增加步长控制"当作根治方案。
- 粘性正则化应作为本构层稳定化手段，和外层 continuation 配合使用。
- 切线一致性（牵引力用 D_visc，切线中 dD/dδ 乘 β）是实现中最容易出错的地方，必须仔细验证。
