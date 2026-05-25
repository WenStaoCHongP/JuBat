# CZM 粘性正则化设计规格

> 日期: 2026-05-12
> 状态: revised (v2)
> 关联计划: `docs/superpowers/plans/2026-05-12-czm-viscous-regularization-plan.md`

## 1. 问题定义

当前耦合例程中的 CZM 采用速率无关双线性牵引-分离律。调试输出表明，损伤变量 `D_max` 会在后峰软化阶段升至约 0.895，随后 Newton / continuation 失去可跟踪分支，表现为残差放大、`result.converged = false`，以及求解器在最小步长附近停滞。

这类失稳更接近"软化切线过软导致的病态"，而不是数值噪声或 NaN。单纯缩小全局步长或切换求解器不能改变局部本构的病态程度，因此需要在内聚力本构内部加入粘性正则化。

**注意**：当前默认求解器为 `"basic"`（`Option.jl:76`），一次性施加全部载荷，无子步结构。粘性正则化在 `basic` 模式下效果有限（见第 3.4 节）。开启粘性正则化时，建议配合 `"load_substep"` 或 `"arc_length"` 使用。

## 2. 设计目标

1. 让软化段的有效响应更平滑，改善 Newton 雅可比的条件数。
2. 不改变默认关闭状态下的物理和数值结果。
3. 保证损伤历史回滚安全，失败步不能提交状态。
4. 让正则化在 `τ_v → 0` 时退化回当前实现。

## 3. 推荐方案

### 3.1 核心公式

采用一阶滞后型损伤正则化。先由当前分离量和历史变量计算等效损伤 `D_eq`（与当前双线性律完全一致），再通过松弛得到粘性损伤 `D_visc`：

```
D_visc_{n+1} = D_visc_n + β * (D_eq_{n+1} - D_visc_n)
```

其中松弛因子：

```
β = Δs / (τ_v* + Δs)
```

`Δs` 是当前载荷子步的分数增量（dimensionless，取值 0~1），`τ_v*` 是归一化后的粘性松弛时间。为确保损伤单调性：

```
D_visc_{n+1} = max(D_visc_n, D_visc_{n+1})
```

### 3.2 τ_v 的归一化方案

用户在 `Option` 中输入物理松弛时间 `czm_visc_tau`（单位：秒）。

归一化遵循项目的统一时间尺度 `t0 = 3600` s：

```
τ_v* = czm_visc_tau / t0
```

这确保 τ_v* 与求解中使用的无量纲时间一致。在 `SetParams.jl` 中完成此归一化。

**推荐物理值**：`czm_visc_tau` 取 10~100 s 量级（对应 β ≈ 0.01~0.1，取决于 Δs），足以平滑软化段而不显著改变物理响应。

### 3.3 Δs 的跨求解器定义

Δs 统一定义为**载荷分数增量**（dimensionless, 0 to 1）：

| 求解器 | Δs 取值 | 说明 |
|--------|---------|------|
| `basic` | `1.0` | 无子步，全量载荷。β = 1/(τ_v* + 1)，正则化效果弱 |
| `load_substep` | `step_size = 1/n_load_steps` | 默认 step_size = 0.5（n_load_steps=2）。β ≈ Δs/τ_v*，正则化效果显著 |
| `arc_length` | `delta_lambda` | 弧长控制的自适应增量。β ≈ Δs/τ_v*，正则化效果显著 |

**关键推论**：`basic` 模式下 Δs = 1.0，即使 τ_v* 较小（如 0.01），β 也接近 1.0，几乎没有正则化效果。因此粘性正则化需要配合载荷子步求解器使用才能发挥作用。

### 3.4 牵引力修改

牵引力使用 `D_visc` 替代原始 `D`：

```
T_n = (1 - D_visc) * K_n * δ_n     (δ_n ≥ 0)
T_t = K_t * δ_t                      (model1)
T_t = (1 - D_visc) * K_t * δ_t      (mix)
```

效果：D_visc < D_eq（滞后），等效刚度比原始大，软化被推迟。

### 3.5 切线刚度修改（关键）

切线刚度必须同时反映 `D_visc` 的使用和 `dD_visc/dδ` 的贡献。当 β 在 Newton 迭代内视为常数时：

```
dD_visc/dδ = β * dD_eq/dδ
```

因此软化分支的切线变为（以 model1 法向为例）：

```
dT_n/dδ_n = (1 - D_visc) * K_n - K_n * δ_n * β * dD_eq/dδ_n
```

对比原始：
```
dT_n/dδ_n = (1 - D) * K_n - K_n * δ_n * dD_eq/dδ_n
```

两个效应叠加：
1. **D_visc < D_eq** → `(1 - D_visc) > (1 - D_eq)` → 弹性部分刚度更大
2. **β < 1** → 软化项 `K_n * δ_n * β * dD_eq/dδ_n` 被缩小 → 切线下降更慢

净效果：切线刚度更平滑，不再在 D≈0.9 时急剧趋于零。

**实现要求**：`bilinear_tangent` 需要额外接收 `β` 参数（或在内部计算），将所有 `dD_eq/dδ` 项乘以 `β`。

### 3.6 关闭时的退化

当 `czm_viscous_enabled = false` 或 `τ_v = 0` 时：
- `β = 1.0`（因为 τ_v* = 0）
- `D_visc = D_visc + 1.0 * (D_eq - D_visc) = D_eq`
- 切线中的 `β * dD_eq/dδ = dD_eq/dδ`

完全退化为当前实现。

## 4. 状态与回滚语义

### 4.1 状态分类

把 CZM 状态分成三类：

- **等效损伤** `D_eq` + 历史极值量（`δ_max_n`, `δ_max_t`, `δ_max_eff`）：由当前分离量按原始双线性律计算，描述瞬时力学状态。
- **粘性损伤** `D_visc`：用于牵引和切线的有效损伤，通过滞后松弛得到。
- **已提交状态** vs **trial 状态**：仅收敛后提交，失败步回滚。

### 4.2 更新顺序

1. 以当前位移/分离量按原始双线性律计算 `D_eq_trial`。
2. 用 `D_eq_trial` 和上一轮已提交的 `D_visc`，按松弛公式计算 `D_visc_trial`。
3. 在 Newton 子步内使用 `D_visc_trial` 组装残差和切线（切线中 `dD/dδ` 项乘以 β）。
4. 只有当该步收敛后，才把 `D_visc_trial` 和 `D_eq_trial` 提交到正式状态。
5. 若该步失败，回滚到上一个已提交状态（现有的 `clone_damage_states` 机制）。

### 4.3 回滚安全性

粘性状态必须和当前求解步绑定。现有的 `clone_damage_states`（`CzmSolve.jl:17-28`）已经深拷贝 `DamageState` 的所有字段。扩展 `DamageState` 增加 `D_visc` 后，`clone_damage_states` 自然会拷贝它，无需额外回滚逻辑。

三个求解器的回滚路径已确认：
- `solve_czm_basic_step`：line 160 备份，line 221-224 恢复
- `solve_czm_arc_length_step`：line 280 备份，line 411-413 恢复
- `newton_raphson_czm`：line 495 备份，line 571-572 恢复

## 5. 文件级接入点

### `src/Option.jl`

新增两个选项：

```julia
czm_viscous_enabled::Bool = false      # 粘性正则化开关
czm_visc_tau::Float64 = 0.0            # 物理松弛时间 [s]
```

### `src/SetParams.jl`

新增 τ_v 的归一化：

```julia
# 在 NormaliseParam 函数中，cohesive 归一化区域后追加：
param.cohesive.tau_visc = opt.czm_visc_tau / param.scale.t0
```

如果 `Cohesive` 结构体需要新增 `tau_visc` 字段，在此处一并添加。

### `src/czm.jl`

扩展 `DamageState`，增加粘性历史字段：

```julia
mutable struct DamageState <: AbstractDamageState
    D::Float64                     # 等效损伤 D_eq（当前双线性律计算结果）
    D_visc::Float64                # 粘性有效损伤（用于牵引和切线）
    δ_max_n::Float64              # 历史最大法向分离
    δ_max_t::Float64              # 历史最大切向分离
    δ_max_eff::Float64            # 历史最大等效分离
    fractured::Bool               # 是否断裂
    accumulated_damage::Float64   # 累积损伤
end
```

`clone_damage_states` 需要同步拷贝 `D_visc`。

### `src/Materialmatrix.jl`

#### `bilinear_traction_state`

增加 `beta` 参数（可选，默认 1.0）。牵引力计算使用 `D_visc`：

```
T_n = (1 - D_visc) * K_n * δ_n    (当 viscous_enabled 且 β 已知时)
```

当 `beta == 1.0` 或 `viscous_enabled == false` 时，行为与当前完全一致。

#### `bilinear_tangent`

增加 `beta` 参数。在软化分支中，将 `dD_eq/dδ` 替换为 `β * dD_eq/dδ`：

```julia
# 当前代码 (Materialmatrix.jl:219):
dD_dδn = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
dT_dδ[1,1] = (1.0 - D) * K_n - K_n * δ_n * dD_dδn

# 修改后:
dD_dδn = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
dT_dδ[1,1] = (1.0 - D_visc) * K_n - K_n * δ_n * beta * dD_dδn
```

`mix` 模式同理，所有 `dD/dδ` 项均乘以 `beta`。

### `src/CzmSolve.jl`

三个求解器函数需要在调用 `assemble_coupled_system` 前计算当前子步的 β：

```julia
# 计算 β
beta = viscous_enabled ? delta_s / (tau_visc + delta_s) : 1.0
```

并将 `beta` 和 `D_visc` 传递到 `assemble_coupled_system` → `assemble_czm_system` → `bilinear_traction_state` / `bilinear_tangent`。

具体 Δs 来源：
- `basic`：Δs = 1.0
- `load_substep`：Δs = `step_size`（当前子步大小）
- `arc_length`：Δs = `delta_lambda`（当前弧长增量）

### `src/CouplingState.jl`

在 `update_czm_damage!` 中：
1. 从 `case.opt` 读取 `czm_viscous_enabled` 和归一化后的 `tau_visc`
2. 将这些参数传递给 `solve_czm_step`
3. 调试输出中增加 `D_visc_max` 信息

### `src/Solve.jl`

如果 `CallModel_MultiSPMe` 中的 CZM 更新路径需要传递时间步信息，确保 `Δt` 可达。当前 CZM 更新通过 `CouplingState.jl` 间接调用，一般无需直接修改此文件，但需要在实现时确认调用链。

### `src/CycleSolver.jl`

循环求解中 CZM 状态跨周期累积。粘性历史 `D_visc` 需要跨周期保留（与 `D` 和 `δ_max` 一致）。当前的 `final_state` 传递机制应自动覆盖，但实现时需确认。

## 6. 推荐实现约束

1. 只在 CZM 路径内部改动，不改变 SPMe / 热模型的物理方程。
2. 默认关闭，避免影响当前所有已验证案例。
3. τ_v 归一化到无量纲时间尺度 `t0 = 3600` s，与项目统一归一化方案一致。
4. 失败步必须以现有提交状态为准，不能污染后续求解。
5. 切线刚度中的 `dD/dδ` 项必须显式乘以 `β`，不能仅替换 `D → D_visc`。
6. 粘性正则化需要配合载荷子步求解器（`"load_substep"` 或 `"arc_length"`）使用。`"basic"` 模式下 Δs = 1.0，正则化效果有限。
7. 先实现法向主导（model1）的最小版本，再考虑混合模式扩展。

## 7. 验证方式

建议按下面顺序验证：

1. **回归测试**：关闭粘性开关（`czm_viscous_enabled = false`），确认结果与当前基线完全一致。
2. **退化测试**：`czm_visc_tau = 0`（τ_v* = 0 → β = 1），确认结果与关闭时一致。
3. **停滞测试**：开启粘性 + `load_substep`，观察 `D_max ≈ 0.895` 区间是否仍停滞。
4. **诊断输出**：检查 `max(δ_n)`、`max(δ_eff)`、`D_eq_max`、`D_visc_max`、`β`、`result.converged` 的变化趋势。
5. **求解器对比**：对比 `load_substep` 与 `arc_length` 在开启粘性后的收敛稳定性。
6. **回滚验证**：确认失败步没有推进损伤历史（`D_visc` 不变）。
7. **求解器限制确认**：确认 `basic` 模式下粘性正则化效果有限（β ≈ 1），结果接近无正则化。

## 8. 预期效果

- 局部软化段的切线不再过于陡峭（β < 1 缩小软化项）。
- Newton 与 continuation 的步进更平稳（D_visc < D_eq 使有效刚度更大）。
- 耦合场下的求解更容易跨过当前的后峰停滞区。

## 9. 不建议的方案

- 只调小全局步长，不改本构。
- 只换求解器方法，不改损伤演化。
- 在失败后补写历史变量，而不是把 trial / commit 语义做严谨。
- 只替换 `D → D_visc` 而不修改切线中的 `dD/dδ` 项（会导致不一致）。

这些方法都不能从根本上解决软化本构导致的局部病态。
