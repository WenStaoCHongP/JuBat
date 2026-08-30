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

**注意**：
- dT/dt 在差分前需要确保两组时间步对齐。若时间步不一致，需先用 `align_to_ref` 将参考解插值到候选解的时间网格上。
- 所有脚本使用 `opt.dtType = "auto"`（自适应时间步），不同网格分辨率会产生不同的时间步序列。插值可能引入额外误差，需在报告中标注。
- 若候选解提前触发截止电压（时间范围短于参考解），RMSPE 只计算两组解共同覆盖的时间范围。

### 3.1.1 与旧规格的关系

旧规格 `2026-04-22-grid-sensitivity-analysis-design.md` §3 定义了以下指标，本文档中未保留：

- **角变化收敛**（旧 §3.2）：由于 Biot 数分析已表明热场近乎轴对称（$Bi_t \approx 0.004$），角变化本身极小，RMSPE 的零点保护会大量跳过，不具备统计意义。改由空间场 RMSPE 间接覆盖。
- **应力峰值**（旧 §3.3）：应力场由扩散应力 + 热应力驱动，非外部加载。其空间分布不均匀，单点峰值意义有限。建议在实施时作为可选补充指标加入。
- **损伤起始时间**（旧 §3.3）：起始时刻是一个事件时间，RMSPE 不适用。改为在 D_max(t) 曲线的 RMSPE 中间接体现。
- **载荷-位移曲线偏差**（旧 §3.3）：保留牵引-分离面积偏差作为替代。载荷-位移曲线需要纯机械模型施加位移边界条件，与电池仿真的实际驱动方式（扩散应力+热应力）不同，优先级较低。

### 3.2 热学 Track

| 指标名称 | 公式 | 数据来源 |
|----------|------|---------|
| 峰值温度曲线 RMSPE | $\epsilon_{\text{RMSPE}}(T_{\max}(t), T_{\max,\text{ref}}(t))$ | 每时间步取 max(T_nodes)，T_nodes 来自 `result.T_hist` |
| 温度范围曲线 RMSPE | $\epsilon_{\text{RMSPE}}(T_{\text{range}}(t), T_{\text{range,ref}}(t))$ | 每时间步取 max(T_nodes) - min(T_nodes) |
| 空间场 RMSPE | $\epsilon_{\text{spatial}}$（§2.3） | `result["thermal2D temperature at nodes [K]"]`（nnode × nt 矩阵） |

**说明**：T_range 在温度场接近均匀时可能接近零，此时零点保护会触发跳过。需在输出中标注。

### 3.3 CZM Track

| 指标名称 | 公式 | 数据来源 |
|----------|------|---------|
| 损伤演化 RMSPE | $\epsilon_{\text{RMSPE}}(D_{\max}(t), D_{\max,\text{ref}}(t))$ | `czm D_max` 时间序列 |
| 断裂数演化 RMSPE | $\epsilon_{\text{RMSPE}}(n_f(t), n_{f,\text{ref}}(t))$ | `czm n_fractured` 时间序列 |
| 分离曲线 RMSPE | $\epsilon_{\text{RMSPE}}(\delta_n(t), \delta_{n,\text{ref}}(t))$ | `czm δ_max_n [m]` 时间序列 |
| 牵引-分离面积偏差 | $\epsilon_A$（§2.4） | 选取 D_max 峰值单元的 traction-separation 数据（见下方说明） |

**牵引-分离面积偏差的单元选择策略**：CZM 输出 `czm traction normal [Pa]` 和 `czm separation normal [m]` 为 (n_coh × nt) 矩阵。选取 D_max 达到最大值的那个单元（即最终时刻 D 值最大的单元），提取其 traction-separation 时间序列构建曲线。若多个单元共享最大 D 值，取其中索引最小的一个。不同 nθ 的 CZM 网格单元数不同，无法逐单元对应，因此只比较"峰值损伤单元"的局部响应。

### 3.4 能量守恒检查

能量残余 $\epsilon_R(t)$ 保留瞬时值形式（能量残余本身是瞬时物理量），但新增归一化 RMS 残余：

$$\epsilon_{R,\text{rms}} = \frac{\sqrt{\frac{1}{N}\sum_{i=1}^{N}R(t_i)^2}}{|W_{\text{elec}}(t_{\text{end}})|} \times 100\%$$

其中 $R(t)$ 是能量平衡残余，$W_{\text{elec}}(t_{\text{end}})$ 是仿真结束时的累积电功。这**不是** RMSPE（因为参考值为零时 RMSPE 无定义），而是绝对残余的 RMS 归一化。

---

## 4. 工具函数

所有脚本共用以下 Julia 工具函数（需 `using Statistics`）：

```julia
using Statistics

"""
    rmspe(y, y_ref; rel_tol=1e-3) -> (rmspe_val, skip_rate)

计算相对均方根百分比误差。跳过 |y_ref| < rel_tol * max(|y_ref|) 的点。
返回 (RMSPE值, 跳过率)。若跳过率 > 50%，调用者应改用绝对误差。
"""
function rmspe(y, y_ref; rel_tol=1e-3)
    threshold = rel_tol * maximum(abs.(y_ref))
    mask = abs.(y_ref) .> threshold
    skip_rate = 1.0 - count(mask) / length(y_ref)
    count(mask) == 0 && return (NaN, 1.0)
    val = sqrt(mean(((y[mask] .- y_ref[mask]) ./ y_ref[mask]).^2)) * 100
    return (val, skip_rate)
end

"""
    spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)

计算空间场 RMSPE 的时间平均值。
T_hist: (nnode × nt) 矩阵，来自 result["thermal2D temperature at nodes [K]"]
"""
function spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)
    nt = size(T_hist, 2)
    errs = Float64[]
    for k in 1:nt
        val, skip = rmspe(T_hist[:,k], T_ref_hist[:,k]; rel_tol)
        isnan(val) || push!(errs, val)
    end
    isempty(errs) && return 0.0  # 所有时步均无空间变化 → 误差为零
    return mean(errs)
end

"""
    area_error(x, y, x_ref, y_ref)

计算归一化曲线面积偏差（梯形积分）。
x 和 y 至少需要 2 个点。
"""
function area_error(x, y, x_ref, y_ref)
    length(x) < 2 && return NaN
    length(x_ref) < 2 && return NaN
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
不依赖外部包，手写线性插值。超出 t_cand 范围的值用端点值填充。
"""
function align_to_ref(t_cand, y_cand, t_ref)
    return [let
        idx = searchsortedfirst(t_cand, t)
        if idx == 1
            y_cand[1]
        elseif idx > length(t_cand)
            y_cand[end]
        else
            frac = (t - t_cand[idx-1]) / (t_cand[idx] - t_cand[idx-1])
            y_cand[idx-1] + frac * (y_cand[idx] - y_cand[idx-1])
        end
    end for t in t_ref]
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
- 新增能量残余的归一化 RMS 统计值 $\epsilon_{R,\text{rms}}$（见 §3.4 公式）

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
- [findings.md](../../planning-with-files/07_网格敏感性分析/findings.md) — 发现记录
