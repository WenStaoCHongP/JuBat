# CZM 单位修复与切线正则化

> 日期: 2026-04-21
> 分支: Parameters_Design
> 涉及文件: `Solve.jl`, `CycleSolver.jl`, `Materialmatrix.jl`, `testexample.jl`

---

## 1. 问题背景

`example/testexample.jl` 启用 `czm_enabled = true` 后，CZM 损伤输出始终为零（D_max = 0）。经排查发现两个层次的问题：

1. **cohesive 参数单位不匹配**：主求解器传递物理参数给归一化框架中的 CZM 求解器，导致刚度矩阵条件数极差、损伤阈值永远达不到。
2. **双线性本构切线不连续**：在软化段（δ > δ_0）切线刚度变负，Newton-Raphson 迭代在弹性-软化边界振荡，无法收敛。

---

## 2. 根因分析

### 2.1 单位不匹配（根本原因）

CZM 求解器内部全部使用归一化量：

| 组件 | 来源 | 单位 |
|------|------|------|
| `E_eff`, `ν_eff` | `compute_czm_effective_params` → `case.param` | 归一化 |
| 网格节点坐标 | `setup_thermal2D_mesh` | 归一化 (x/L) |
| 热化学载荷 | `assemble_thermal_chemical_load` | 归一化 |
| `cohesive_params` | `Solve.jl:277` → `case.param_dim.cohesive` | **物理 (SI)** |

物理 `K_n = 2.4e17 Pa/m` 与归一化 `E_eff ~ O(1)` 混在同一个刚度矩阵中，内聚力刚度比体单元大 ~10¹⁷ 倍。分离量在归一化空间为 ~1e-21，远低于物理阈值 `δ_0_n = 3.4e-10 m`，损伤永远不启动。

**修复**：`SetParams.jl:417-427` 已有完整的 cohesive 归一化逻辑，但调用方未使用。

### 2.2 切线不连续（收敛瓶颈）

`bilinear_tangent` 在 δ_0 处切线从 `+K_n`（弹性）跳到负值（软化），Newton 方向在边界来回振荡。载荷子步法无法解决此本质问题。

---

## 3. 修改内容

### 3.1 cohesive 参数单位修复

**文件**: `src/Solve.jl:277`, `src/CycleSolver.jl:246`

```diff
- case.param_dim.cohesive
+ case.param.cohesive
```

两处调用点均改为传递归一化 cohesive 参数。

### 3.2 切线正则化（割线刚度）

**文件**: `src/Materialmatrix.jl` — `bilinear_tangent` 函数

软化段将可能为负的真实切线替换为正的割线刚度：

```julia
# model1 软化段（原来）:
dT_dδ[1,1] = (1-D)*K_n - K_n*δ_n*dD_dδn   # 可能为负 → Newton 发散

# model1 软化段（现在）:
dT_dδ[1,1] = (1-D)*K_n                       # 始终为正 → Newton 稳定
```

mix 模式同理。弹性段和卸载段不变。牵引力计算（`bilinear_traction_state`）仍用精确双线性律，Newton 方向由割线刚度保证稳定性。

### 3.3 平滑正部函数

**文件**: `src/Materialmatrix.jl` — 新增 `smooth_positive_part` 辅助函数

用 `smooth_positive_part(x, eps) = 0.5*(x + hypot(x, eps))` 替代硬 `max(0, x)`，消除 δ_n 正负切换处的 C⁰ 不连续性，使牵引力-位移关系更光滑。

---

## 4. 验证结果

### 4.1 单位修复效果（basic 方法，czm_tol=1e-4）

| 指标 | 修复前 | 修复后 |
|------|--------|--------|
| D_max | 0.00% | **93.62%** |
| D_mean | 0.00% | **0.57%** |
| δ_max_n | 0 | **1.49e-6 m** |
| 温度范围 | 298.15 ~ 321.59 K | 298.15 ~ 313.18 K |

修复后损伤值合理：少数界面损伤严重（D_max ~94%），大部分界面损伤轻微（D_mean ~0.6%），符合物理预期。

### 4.2 切线正则化 + load_substep + czm_tol=1e-3

| 指标 | 值 |
|------|-----|
| D_max | **99.59%** |
| D_mean | **3.22%** |
| δ_max_n | **2.06e-5 m** |
| 温度范围 | 298.15 ~ 311.77 K |
| 容量 | 4.85 Ah |
| CZM 耗时占比 | 97.5%（383s） |

### 4.3 残差收敛对比

| 方法 | 残差量级 | 收敛状态 |
|------|---------|---------|
| basic + 物理参数 | 9.7 → 3166（递增发散） | 全部不收敛 |
| load_substep + 物理参数 | 0.01 ~ 0.4 | 全部不收敛 |
| load_substep + 归一化参数（原始切线） | 0.01 ~ 0.4 | 全部不收敛 |
| load_substep + 归一化参数 + 割线刚度 | **0.001 ~ 0.002** | 大幅改善，接近 tol |
| load_substep + 割线刚度 + tol=1e-3 | **≤ 0.01** | 大部分收敛 |

---

## 5. 遗留问题与后续方向

### 5.1 CZM 计算耗时（97.5%）

load_substep 方法在部分步仍需 ~900 次 Newton 迭代。割线刚度牺牲了二次收敛性，改为线性收敛。可能的改善方向：

- **减少调用频率**：增大 `czm_update_interval`（如每 5~10 步更新一次）
- **限制最大迭代**：将 `czm_max_iter` 从 100 降至 30~50，避免在不可收敛的步上浪费
- **粘性正则化**：在软化段加入率相关项 `η·dδ/dt`，可恢复部分二次收敛性

### 5.2 完全消除 stalled 警告

当前残差稳定在 ~0.01 水平，部分步仍无法达到 1e-3。根因是割线刚度近似（Modified Newton）只有线性收敛，而非刚度尺度不匹配。

**刚度尺度验证**（2026-04-21 复核）：

| 参数 | 归一化尺度 | 归一化值 |
|------|-----------|---------|
| E_eff | E_n = cs_max × R × T_ref ≈ 8.2e7 | ≈ 210 |
| K_n | K_czm = σ_czm / δ_czm ≈ 1.33e14 | ≈ 1805 |
| **比值 K_n / E_eff** | — | **≈ 8.6** |

8.6 倍的刚度比在 CZM 建模中完全合理（cohesive 应比 bulk 刚以模拟弹性阶段的刚性连接），不存在量级差异问题。

**消除 stalled 的正确方向**：

- **粘性正则化**：在软化段加入率相关项 `η·dδ/dt`，恢复光滑切线和二次收敛
- **指数软化律**：替代双线性律（消除 δ₀ 处 C¹ 不连续的折点）
- **弧长法**（`arc_length`）：直接跟踪软化路径，绕过 snap-through 困难

---

## 6. 修改文件清单

| 文件 | 变更 | 说明 |
|------|------|------|
| `src/Solve.jl:277` | 1 行 | `param_dim.cohesive` → `param.cohesive` |
| `src/CycleSolver.jl:246` | 1 行 | `param_dim.cohesive` → `param.cohesive` |
| `src/Materialmatrix.jl` | ~30 行 | 新增 `smooth_positive_part`；`bilinear_traction_state` 和 `bilinear_tangent` 软化段改用割线刚度 |
