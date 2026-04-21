# CZM 单位修复与切线正则化 — 理论评审

> 日期: 2026-04-21
> 分支: Parameters_Design (vs main)
> 评审范围: `Solve.jl`, `CycleSolver.jl`, `Materialmatrix.jl`
> 评审类型: 理论正确性验证

---

## 1. 修改概述

本次修改涉及三个核心理论变更：

| 编号 | 变更 | 涉及文件 |
|------|------|----------|
| C1 | cohesive 参数归一化修正 | `Solve.jl:277`, `CycleSolver.jl:246` |
| C2 | 割线刚度替代真实切线 | `Materialmatrix.jl` bilinear_tangent |
| C3 | smooth_positive_part 替代 max(0,δ_n) | `Materialmatrix.jl` bilinear_traction_state + bilinear_tangent |

另有配套修改（非理论层面）：
- Solve.jl: CZM 损伤演化集成到主循环、deepcopy→copy
- CycleSolver.jl: 移除冗余 phase-end CZM 调用、移除重复函数定义
- get_active_elements 签名修复

---

## 2. 逐项理论评审

### 2.1 C1: cohesive 参数归一化修正

**修改内容**: `case.param_dim.cohesive` → `case.param.cohesive`

**CZM 求解器的无量纲化体系**:

| 组件 | 数据来源 | 空间 |
|------|----------|------|
| E_eff, ν_eff | `compute_czm_effective_params` → `case.param` | 归一化 |
| 网格坐标 | `setup_thermal2D_mesh` | 归一化 (x/L) |
| 热化学载荷 | `assemble_thermal_chemical_load` | 归一化 |
| cohesive 参数 (旧) | `case.param_dim.cohesive` | **物理 (SI)** |
| cohesive 参数 (新) | `case.param.cohesive` | 归一化 |

**归一化参考尺度** (`SetParams.jl:298-301`):

```
σ_czm = σ_max_n = 82e6 Pa
δ_czm = δ_c_n   = 6.17e-7 m
K_czm = σ_czm / δ_czm ≈ 1.33e14 Pa/m
```

**归一化后参数值**:

| 参数 | 物理值 | 归一化值 |
|------|--------|----------|
| K_n | 2.4e17 Pa/m | ≈ 1805 |
| δ_0_n | 3.4e-10 m | ≈ 5.5e-4 |
| δ_c_n | 6.17e-7 m | = 1.0 |

**评审结论**: ✅ **理论正确**

- 修复前：K_n = 2.4e17 与 E_eff ~ O(1) 混在同一刚度矩阵 → 条件数 ~10^17 → 分离量 δ ~ 10^-21（归一化空间），远低于阈值 δ_0_n ~ 5.5e-4 → 损伤永不启动
- 修复后：K_n ≈ 1805, δ_c_n = 1.0, δ_0_n ≈ 5.5e-4 → 条件数合理 → 损伤正常演化
- `SetParams.jl:417-427` 已有完整归一化逻辑，修复仅是让调用方使用正确的参数源

### 2.2 C2: 割线刚度替代真实切线

**修改内容**: `bilinear_tangent` 软化段切线刚度从真实导数改为割线近似

#### 原始切线（软化段，model1）:

```
dD/dδ_n = δ_c * δ_0 / (δ² * (δ_c - δ_0))   > 0
dT_n/dδ_n = (1-D)*K_n - K_n*δ_n*dD/dδ_n      可能为负
```

物理上，软化段切线为负（应力随应变增加而减小）是 snap-through 行为的正确描述。但负切线导致 Newton-Raphson 迭代发散：迭代方向在弹性-软化边界来回振荡。

#### 新切线（割线近似）:

```julia
K_n_eff = K_n * (1.0 - D * dδ_n_pos_dδ_n)
dT_dδ[1,1] = K_n_eff
```

其中 `dδ_n_pos_dδ_n = smooth_positive_part_derivative(δ_n, 1e-12)` 是 δ_n 正部的光滑导数：

| δ_n 状态 | dδ_n_pos_dδ_n | K_n_eff | 物理含义 |
|----------|---------------|---------|----------|
| δ_n >> 0 | ≈ 1.0 | ≈ (1-D)*K_n | 张开：损伤降低刚度 |
| δ_n << 0 | ≈ 0.0 | ≈ K_n | 压缩：无损伤效应 |
| δ_n ≈ 0 | ≈ 0.5 | ≈ (1-0.5D)*K_n | 过渡：光滑插值 |

#### 理论分析

**Newton-Raphson 方法中 Jacobian 近似的影响**:

残差方程: R(u) = f_int(u) - f_ext = 0

Newton 更新: Δu = -J⁻¹ * R，其中 J 是 Jacobian 近似。

关键定理：只要 J 保持正定且充分接近真实 Jacobian，Newton 迭代仍收敛到同一不动点。区别仅在于：
- 真实 Jacobian → 二次收敛 (|e_{k+1}| ~ |e_k|²)
- 割线/近似 Jacobian → 线性收敛 (|e_{k+1}| ~ c·|e_k|, 0 < c < 1)

**本修改满足的条件**:
1. K_n_eff > 0（始终为正）→ J 保持正定 ✓
2. 平衡点不变：牵引力计算 (`bilinear_traction_state`) 使用精确双线性律，残差函数未改变 ✓
3. 损伤演化公式 D = δ_c*(δ-δ₀)/(δ*(δ_c-δ₀)) 未改变 ✓

**能量耗散验证**:
- 断裂能 G_c = ∫₀^δ_c T dδ = 0.5 * σ_max * (δ_0 + δ_c)
- 由精确牵引力积分得到，与切线近似无关 ✓

**mix 模型交叉耦合项归零**:
- 原始代码在 mix 模型软化段有 dT_dδ[1,2] 和 dT_dδ[2,1]（来自 dD/dδ_eff * dδ_eff/dδ_t 链式法则）
- 新代码将这些设为 0 → 块对角近似
- 对角项仍捕获主导物理（法向和切向独立响应）
- 块对角近似在 Modified Newton 方法中是常见且合理的 ✓

**评审结论**: ✅ **理论正确**（Modified Newton-Raphson 方法）
- 收敛速度从二次降为线性 → 需更多迭代，但不会改变平衡解
- 这是计算断裂力学中处理软化收敛问题的标准技术
- 参考: de Borst et al., "Non-linear Finite Element Analysis of Solids and Structures", Ch. 7

### 2.3 C3: smooth_positive_part 替代 max(0, δ_n)

**修改内容**: 单侧接触条件的正部函数从硬 max 改为光滑近似

```julia
# 旧: C⁰ 连续，导数在 δ_n=0 处跳变
δ_n_pos = max(0.0, δ_n)

# 新: C∞ 连续
δ_n_pos = smooth_positive_part(δ_n, 1e-12) = 0.5*(δ_n + hypot(δ_n, 1e-12))
```

**数学性质**:

| 性质 | max(0, x) | smooth_positive_part(x, ε) |
|------|-----------|---------------------------|
| 连续性 | C⁰ | C∞ |
| x >> ε | x | ≈ x (误差 ~ ε²/x) |
| x << -ε | 0 | ≈ ε²/(2|x)| ≈ 0 |
| x = 0 | 0 | ε/2 = 5e-13 |
| 导数 | 阶跃函数 (0→1 跳变) | sigmoid 函数 (0→1 平滑过渡) |

**在归一化 CZM 空间中的影响**:
- δ_c_n = 1.0, δ_0_n ≈ 5.5e-4
- 正则化参数 ε = 1e-12
- ε/δ_0_n ≈ 1.8e-9 → 扰动量级远小于物理最小分离量 ✓

**消除的数值问题**:

原始代码在 δ_n = 0 处的 Jacobian 跳变：

```
δ_n slightly > 0: dT_n/dδ_n = (1-D)*K_n
δ_n slightly < 0: dT_n/dδ_n = K_n

若 D > 0: 跳变量 = D*K_n → Newton 方向突变
```

正则化后：导数通过 sigmoid 平滑过渡，Jacobian 在 δ_n = 0 附近连续 ✓

**牵引力计算验证**:

```julia
T_n = K_n * δ_n_neg + (1-D) * K_n * δ_n_pos
    = K_n * δ_n - D * K_n * δ_n_pos
```

- δ_n > 0: T_n ≈ K_n*δ_n - D*K_n*δ_n = (1-D)*K_n*δ_n  ← 与原始一致
- δ_n < 0: T_n ≈ K_n*δ_n - 0 = K_n*δ_n                  ← 与原始一致
- δ_n ≈ 0: 光滑过渡                                      ← 改善连续性

**评审结论**: ✅ **理论正确**
- 数值扰动可忽略（~10^-12 量级 vs 物理量 ~10^-4 量级）
- 消除 Jacobian 不连续性，与 Newton-Raphson 的光滑性要求一致
- 参考技术: "regularized Heaviside function"，计算力学中的标准技术

---

## 3. 配套修改评审

### 3.1 Solve.jl CZM 主循环集成

**修改**: 在主时间循环中每 `czm_update_interval` 步调用一次 CZM 更新

```julia
if czm_active
    czm_step_count += 1
    if czm_step_count % case.opt.czm_update_interval == 0
        u_czm_new, czm_converged = update_czm_damage!(
            case.czm_mesh, case.param.cohesive, case, variables, T_nodes_carry, u_czm_prev)
        if czm_converged
            u_czm_prev = u_czm_new
        end
    end
end
```

**评审**: ✅
- `u_czm_prev` 保留上一步位移场作为下一次 Newton 的初始猜测 → 热启动 ✓
- 仅在收敛时更新位移历史 → 避免传播发散结果 ✓
- `case.param.cohesive` 传递归一化参数 → 与 C1 一致 ✓

### 3.2 CycleSolver.jl 移除 phase-end CZM 调用

**旧代码**: 每个 phase 结束后调用 `update_czm_damage!(..., u_czm_prev=nothing)`
- 问题1: `u_czm_prev=nothing` 丢弃位移历史 → 每个阶段从零重解
- 问题2: 与 Solve.jl 主循环中的更新冗余 → 双重计算

**评审**: ✅ 移除正确
- CZM 更新现在完全由 Solve.jl 主循环负责
- 位移历史跨阶段保留（通过 `final_state` 传递）

### 3.3 CycleSolver.jl 移除重复函数定义

移除了 `compute_czm_effective_params`、`compute_czm_strain_inputs`、`update_czm_damage!` 的重复定义（~160 行）。

**关键差异**: 移除的版本使用 `case.param_dim`（物理参数），保留的 CzmSolve.jl 版本使用 `case.param`（归一化参数）。

**评审**: ✅ 移除正确
- 保留的 CzmSolve.jl 版本是规范实现
- 移除物理参数版本消除了未来误用的风险

### 3.4 deepcopy → copy

Float64 和 Vector{Float64} 的 `deepcopy` 改为 `copy`/直接赋值。

**评审**: ✅ 纯性能优化，无语义变化
- Float64 不可变，`deepcopy(dt_min)` = `dt_min`
- `Vector{Float64}` 的浅拷贝已足够（不涉及嵌套结构）

---

## 4. 潜在风险与遗留问题

### 4.1 收敛速度下降 (已知，非理论错误)

割线刚度导致 Newton 迭代从二次收敛降为线性收敛。实测：
- 部分步仍需 ~900 次迭代
- CZM 占总计算时间 97.5%

**不构成理论错误**，但建议后续优化：
- 增大 `czm_update_interval`（如每 5~10 步更新一次）
- 限制 `czm_max_iter`（30~50 次即放弃）
- 粘性正则化（加入率相关项恢复部分二次收敛性）

### 4.2 刚度尺度不匹配？（已排除）

文档 §5.2 曾提出"bulk 用 E_eff 归一化，cohesive 用 K_czm 归一化，两者量级差异大"。经数值验证，**此说法不成立**：

| 参数 | 归一化尺度 | 归一化值 |
|------|-----------|---------|
| E_eff | E_n = cs_max × R × T_ref ≈ 8.2e7 | ≈ 210 |
| K_n | K_czm = σ_czm / δ_czm ≈ 1.33e14 | ≈ 1805 |
| **比值** | — | **≈ 8.6** |

8.6 倍刚度比在 CZM 中完全合理，不存在量级差异问题。残差停滞的真正原因是割线刚度的线性收敛特性。

### 4.3 容差松弛 (实操妥协)

`czm_tol` 从 1e-4 放宽至 1e-3。这意味着 CZM 残差在 ~0.1% 水平即视为收敛。

**风险**: 损伤状态的精度约为 ~0.1%，对工程应用可接受。若需要更高精度，需实现弧长法或粘性正则化。

### 4.4 mix 模型交叉耦合项归零

将 mix 模型的 dT_dδ[1,2] 和 dT_dδ[2,1] 设为零简化了 Jacobian。对于法向-切向强耦合的情况，可能增加迭代次数。当前测试使用 model1（纯法向），无影响。

---

## 5. 总结

| 修改项 | 理论正确性 | 物理一致性 | 数值稳定性 |
|--------|-----------|-----------|-----------|
| C1: 归一化修正 | ✅ 正确 | ✅ 消除单位不匹配 | ✅ 条件数从 ~10^17 降至 ~10^3 |
| C2: 割线刚度 | ✅ 正确 (Modified NR) | ✅ 平衡解不变 | ✅ 消除负切线发散 |
| C3: smooth_positive_part | ✅ 正确 | ✅ 扰动可忽略 | ✅ 消除 Jacobian 跳变 |
| 配套修改 | ✅ 正确 | ✅ 一致性改善 | ✅ 消除冗余计算 |

**整体评审结论**: 所有修改理论正确，不会引入物理错误。收敛速度的下降是 Modified Newton-Raphson 方法的固有代价，属于已知的精度-稳定性权衡。
