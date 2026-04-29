# 网格敏感性分析：统计指标体系设计

> **Date:** 2026-04-29
>
> **Status:** Draft
>
> **Supersedes:** 部分 `2026-04-22-grid-sensitivity-analysis-design.md` 的 §3 Metric Definitions
>
> **Source discussion:** 用户要求将网格敏感性分析的对比指标从"偶然点值"改为统计值

---

## 1. 背景与动机

现有网格敏感性分析脚本（`example/网格敏感性/` Script 2-5）中的所有对比指标都基于**单点值**（snapshot）：

- 电化学 Track：`V_end` 端电压点值差、`T_peak` 峰值点值差、`max|dT/dt|` 极值点值差
- 热学 Track：`T_max` 终态点值差
- CZM Track：`D_max` 终态点值差

**问题**：点值对局部波动高度敏感，无法反映整条曲线或空间场的整体收敛质量。

**目标**：将所有对比指标统一替换为基于 **RMSPE**（相对均方根百分比误差）的统计量，覆盖时间历程、空间场和曲线面积三个维度。

---

## 2. 核心误差公式

### 2.1 RMSPE（相对均方根百分比误差）

$$\epsilon_{\text{RMSPE}}(y, y_{\text{ref}}) = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(\frac{y_i - y_{\text{ref},i}}{|y_{\text{ref},i}|}\right)^2} \times 100\%$$

### 2.2 零点保护

当 $|y_{\text{ref},i}| < \delta$ 时跳过该采样点，其中：

$$\delta = 10^{-3} \times \max(|y_{\text{ref}}|)$$

报告中需标注跳过率。若跳过率超过 50%，该指标的 RMSPE 不可信，需改用绝对误差。

### 2.3 空间场 RMSPE（时间平均）

$$\epsilon_{\text{spatial}} = \frac{1}{N_t}\sum_{k=1}^{N_t}\sqrt{\frac{1}{N_n}\sum_{n=1}^{N_n}\left(\frac{T_n^{(k)} - T_{n,\text{ref}}^{(k)}}{|T_{n,\text{ref}}^{(k)}|}\right)^2} \times 100\%$$

### 2.4 曲线面积偏差

$$\epsilon_A = \frac{|A - A_{\text{ref}}|}{|A_{\text{ref}}|} \times 100\%$$

其中 $A$ 为曲线的梯形积分面积。

### 2.5 验收标准

所有 RMSPE 指标 < 5%。

---

## 3. 各 Track 指标定义

### 3.1 电化学 Track

| 指标名称 | 公式 | 数据来源 |
|----------|------|---------|
| 电压曲线 RMSPE | $\epsilon_{\text{RMSPE}}(V(t), V_{\text{ref}}(t))$ | `cell voltage [V]` 时间序列 |
| 温度曲线 RMSPE | $\epsilon_{\text{RMSPE}}(T(t), T_{\text{ref}}(t))$ | `temperature [K]` 时间序列 |
| 温度梯度曲线 RMSPE | $\epsilon_{\text{RMSPE}}(dT/dt(t), dT/dt_{\text{ref}}(t))$ | 温度时间序列差分 |

**注意**：dT/dt 在差分前需要确保两组时间步对齐。若时间步不一致，需先将参考解插值到候选解的时间网格上。

### 3.2 热学 Track

| 指标名称 | 公式 | 数据来源 |
|----------|------|---------|
| 峰值温度曲线 RMSPE | $\epsilon_{\text{RMSPE}}(T_{\max}(t), T_{\max,\text{ref}}(t))$ | 每时间步取 max(T_nodes) |
| 温度范围曲线 RMSPE | $\epsilon_{\text{RMSPE}}(T_{\text{range}}(t), T_{\text{range,ref}}(t))$ | 每时间步取 max-min |
| 空间场 RMSPE | $\epsilon_{\text{spatial}}$（§2.3） | T_hist 节点温度时间序列 |

**说明**：T_range 在温度场接近均匀时可能接近零，此时零点保护会触发跳过。需在输出中标注。

### 3.3 CZM Track

| 指标名称 | 公式 | 数据来源 |
|----------|------|---------|
| 损伤演化 RMSPE | $\epsilon_{\text{RMSPE}}(D_{\max}(t), D_{\max,\text{ref}}(t))$ | `czm D_max` 时间序列 |
| 断裂数演化 RMSPE | $\epsilon_{\text{RMSPE}}(n_f(t), n_{f,\text{ref}}(t))$ | `czm n_fractured` 时间序列 |
| 分离曲线 RMSPE | $\epsilon_{\text{RMSPE}}(\delta_n(t), \delta_{n,\text{ref}}(t))$ | `czm δ_max_n [m]` 时间序列 |
| 牵引-分离面积偏差 | $\epsilon_A$（§2.4） | traction-separation 数据 |

### 3.4 能量守恒检查

能量残余 $\epsilon_R(t)$ 保留瞬时值形式（能量残余本身是瞬时物理量），但新增：

$$\epsilon_{R,\text{RMSPE}} = \text{RMSPE}(\epsilon_R(t), \mathbf{0})$$

即能量残余的 RMS 值占电功的比例，作为整体守恒质量的统计度量。

---

## 4. 工具函数

所有脚本共用以下 Julia 工具函数：

```julia
"""
    rmspe(y, y_ref; rel_tol=1e-3)

计算相对均方根百分比误差。跳过 |y_ref| < rel_tol * max(|y_ref|) 的点。
"""
function rmspe(y, y_ref; rel_tol=1e-3)
    threshold = rel_tol * maximum(abs.(y_ref))
    mask = abs.(y_ref) .> threshold
    count(mask) == 0 && return NaN
    return sqrt(mean(((y[mask] .- y_ref[mask]) ./ y_ref[mask]).^2)) * 100
end

"""
    spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)

计算空间场 RMSPE 的时间平均值。
T_hist: (nnode × nt) 矩阵
"""
function spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)
    nt = size(T_hist, 2)
    errs = [rmspe(T_hist[:,k], T_ref_hist[:,k]; rel_tol) for k in 1:nt]
    return mean(filter(!isnan, errs))
end

"""
    area_error(x, y, x_ref, y_ref)

计算归一化曲线面积偏差（梯形积分）。
"""
function area_error(x, y, x_ref, y_ref)
    A  = abs(trapz(x, y))
    Ar = abs(trapz(x_ref, y_ref))
    Ar == 0 && return NaN
    return abs(A - Ar) / Ar * 100
end

"""
    trapz(x, y)

简单梯形积分。
"""
function trapz(x, y)
    return sum(0.5 .* diff(x) .* (y[2:end] .+ y[1:end-1]))
end

"""
    align_to_ref(t_cand, y_cand, t_ref)

将候选解插值到参考解的时间网格上，返回对齐后的 y_cand_aligned。
"""
function align_to_ref(t_cand, y_cand, t_ref)
    return linear_interpolation(t_cand, y_cand).(t_ref)
end
```

---

## 5. 脚本改动清单

### 5.1 Script 2（电化学）

**改动范围**：仅后处理逻辑，仿真配置不变。

**汇总表**：
- 原 `V_end [V]` + `V err%` → 新 `V(t) RMSPE%`
- 原 `T_peak [K]` + `T err%` → 新 `T(t) RMSPE%`
- 原 `dTdt_max [K/s]` + `dTdt err%` → 新 `dT/dt RMSPE%`

**图**：
- 图 1-2（电压/温度曲线对比）：不变
- 图 3（收敛误差）：横轴不变，纵轴改为 RMSPE%

### 5.2 Script 3（热学）

**改动范围**：仅后处理逻辑。

**汇总表**：
- 原 `T_max err%` → 新 `T_max(t) RMSPE%`
- 新增 `T_range(t) RMSPE%`
- 新增 `Spatial RMSPE%`

**图**：
- 图 3（收敛误差）：改为 RMSPE%
- 新增图：空间场 RMSPE 随 nθ 变化

### 5.3 Script 4（CZM）

**改动范围**：仅后处理逻辑。

**汇总表**：
- 原 `D_max err%` → 新 `D_max(t) RMSPE%`
- 新增 `n_frac(t) RMSPE%`
- 新增 `δ_max_n(t) RMSPE%`
- 新增 `Traction-Sep 面积偏差%`

**图**：
- 图 4（收敛误差）：改为 log-scale RMSPE

### 5.4 Script 5（能量守恒）

**改动范围**：后处理逻辑。

- $\epsilon_R(t)$ 瞬时曲线保留
- 新增能量残余的 RMSPE 统计值

---

## 6. 不变的部分

以下内容不受本次改动影响：

- 仿真配置（网格候选值、模型选择、时间步策略）
- 物理判据（Biot 数分析、cohesive 特征长度、l_c 计算）
- 验收阈值（统一 5%）
- 脚本 1（特征长度计算脚本）

---

## 7. 参考文件

- [2026-04-22-grid-sensitivity-analysis-design.md](2026-04-22-grid-sensitivity-analysis-design.md) — 原始设计规格（§3 指标定义被本文档取代）
- [2026-04-22-grid-sensitivity-analysis-plan.md](../plans/2026-04-22-grid-sensitivity-analysis-plan.md) — 实施计划
- [findings.md](../../planning-with-files/网格敏感性分析/findings.md) — 发现记录
