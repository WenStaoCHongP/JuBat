# 网格收敛性分析：基于 GCI 的标准框架设计

> **Date:** 2026-06-01
>
> **Status:** Draft
>
> **Supersedes:** `2026-04-29-grid-sensitivity-statistical-metrics-design.md`
>
> **Source discussion:** 用户要求基于学术文献（Roache 1998, ASME V&V 10-2019, Turon et al. 2007）重新设计网格收敛性分析方法论

---

## 1. 背景与动机

### 1.1 旧方法的问题

旧规格 `2026-04-29-grid-sensitivity-statistical-metrics-design.md` 使用 **RMSPE**（相对均方根百分比误差）作为核心收敛指标。经文献调研，存在以下问题：

1. **RMSPE 不是 FEM 收敛性分析的标准指标**。FEM 文献和 ASME V&V 10-2019 标准使用 L2 范数、H1 半范数和 GCI (Grid Convergence Index)。
2. **零点除法问题**。RMSPE 在参考值接近零时需要跳过采样点，跳过率超过 50% 时指标不可信。L2 范数用全局范数归一化，不存在此问题。
3. **缺少 Richardson 外推和 GCI**。Roache (1998) 提出的 GCI 是量化网格不确定性的行业标准，可给出置信区间。
4. **CZM 指标不够稳健**。D_max 的逐点 RMSPE 对损伤前沿位置偏移敏感。Turon et al. (2007) 推荐以积分量（断裂能耗散）为主要指标。
5. **缺少收敛阶分析**。观测收敛阶 $p$ 与理论阶数比较是验证数值方法实现正确性的关键手段。

### 1.2 目标

采用严格学术标准重新设计网格收敛性分析框架，使其：

- 符合 ASME V&V 10-2019 标准要求
- 提供 GCI 量化网格不确定性
- 计算观测收敛阶并与理论值比较
- 适用于论文发表（J. Power Sources、J. Electrochem. Soc. 等）

### 1.3 新旧文件关系

| 旧文件 | 新文件 | 说明 |
|--------|--------|------|
| `example/网格敏感性/` | `example/网格敏感性_v2/` | 新目录，完全独立 |
| `0_rmspe_utils.jl` | `0_convergence_utils.jl` | 工具函数重写 |
| `2_electrochemical_mesh_sensitivity.jl` | `2_electrochemical_convergence.jl` | Script 2 重写 |
| `3_thermal_mesh_sensitivity.jl` | `3_thermal_convergence.jl` | Script 3 重写 |
| `4_czm_mesh_sensitivity.jl` | `4_czm_convergence.jl` | Script 4 重写 |
| `5_energy_conservation_check.jl` | `5_energy_conservation.jl` | Script 5 重写 |
| `1_cohesive_characteristic_length.jl` | `1_cohesive_characteristic_length.jl` | 从原版复制，不变 |

---

## 2. 文献基础

### 2.1 Richardson 外推与 GCI

**核心参考文献**：

- **Roache, P.J. (1998)**. *Verification and Validation in Computational Science and Engineering*. Hermosa Publishers.
- **Celik, I.B. et al. (2008)**. "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications." *ASME J. Fluids Eng.*, 130(7), 078001.
- **Schwer, L. (2008)**. "Is Your Mesh Refined Enough? Estimating Discretization Error via GCI." *DynaMORE*.

**GCI 定义**：

$$\text{GCI}_{\text{fine}} = \frac{F_s \cdot |\varepsilon_{21}|}{r_{21}^p - 1} \times 100\%$$

其中：
- $\varepsilon_{21} = (f_2 - f_1) / f_1$ 是相邻两级网格的标量结果相对变化
- $r_{21} = h_2 / h_1$ 是网格细化比（$h$ 为特征单元尺寸）
- $p$ 是观测收敛阶
- $F_s$ 是安全系数

**安全系数**：

| 网格级数 | $F_s$ |
|----------|-------|
| 2 级 | 3.0 |
| 3+ 级 | 1.25 |

$F_s = 1.25$ 代表 95% 置信区间。

### 2.2 ASME V&V 10-2019

**核心要求**：

1. 代码验证必须在解验证之前完成
2. 至少 3 个网格级别（粗、中、细）
3. 使用 Richardson 外推确定观测收敛阶
4. GCI 是量化网格不确定性的标准指标
5. 积分量（应变能、反力、总热通量）优先于点值用于收敛评估
6. 应力奇异点处不期望点态收敛

### 2.3 CZM 网格收敛性

**核心参考文献**：

- **Turon, A. et al. (2007)**. "An engineering solution for mesh size effects in the simulation of delamination using cohesive zone models." *Eng. Fract. Mech.*, 74(10), 1665-1682.

**关键结论**：

- 双线性牵引-分离律需要至少 **3 个单元**跨越 cohesive 过程区长度 $l_{cz}$
- $l_{cz} = E' G_{Ic} / \sigma_{\max}^2$
- 积分量（断裂能耗散）是最可靠的收敛指标
- 网格过粗时，可降低 $\sigma_{\max}$ 同时保持 $G_{Ic}$ 不变来扩大 $l_{cz}$

---

## 3. 核心误差度量

### 3.1 离散范数

**L2 范数**（全局误差测度）：

$$\|e\|_{L_2} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(f_i - f_{\text{ref},i})^2}$$

**归一化 L2 范数**（跨物理场可比较）：

$$\|e\|_{L_2,\text{rel}} = \frac{\|f - f_{\text{ref}}\|_{L_2}}{\|f_{\text{ref}}\|_{L_2}}$$

其中 $\|f_{\text{ref}}\|_{L_2} = \sqrt{\frac{1}{N}\sum f_{\text{ref},i}^2}$。

**与 RMSPE 的关键区别**：归一化 L2 范数的分母是整个参考场的 L2 范数（永远非零），而非逐点参考值（可能为零）。

**H1 半范数**（梯度/通量误差测度）：

$$|e|_{H_1} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}(\nabla f_i - \nabla f_{\text{ref},i})^2}$$

**L∞ 范数**（最坏情况局部误差）：

$$\|e\|_{\infty} = \max_i |f_i - f_{\text{ref},i}|$$

### 3.2 观测收敛阶

对于恒定细化比 $r$：

$$p = \frac{\ln|(f_3 - f_2)/(f_2 - f_1)|}{\ln(r)}$$

对于非恒定细化比（本项目情况），遵循 Celik et al. (2008) 的迭代求解法。

定义：

$$\alpha = \frac{f_3 - f_2}{f_2 - f_1}$$

迭代方程（固定点迭代）：

$$p_{k+1} = \frac{1}{\ln(r_{21})} \left| \ln|\alpha| + \ln\left(\frac{r_{21}^{p_k} - 1}{r_{32}^{p_k} - 1}\right) \right|$$

迭代步骤：

1. 取初始值 $p_0 = \ln|\alpha| / \ln(r_{21})$（恒定比近似）
2. 固定点迭代：$p_{k+1} = \frac{1}{\ln(r_{21})} |\ln|\alpha| + \ln\frac{r_{21}^{p_k} - 1}{r_{32}^{p_k} - 1}|$
3. 收敛判据：$|p_{k+1} - p_k| < 10^{-6}$
4. 最大迭代次数：100（防止发散）
5. 若发散（$p_k > 20$ 或 NaN），回退到恒定比近似值 $p_0$

### 3.3 等效单元尺寸

对于 Jellyroll 螺旋网格，每级网格的等效 $h$ 定义为：

$$h = \sqrt{\frac{A_{\text{total}}}{N_{\text{elem}}}}$$

其中 $A_{\text{total}}$ 是热网格总面积，$N_{\text{elem}}$ 是单元数。这是 ASME V&V 10 推荐的等效 $h$ 定义。

### 3.4 GCI 适用范围说明

**GCI 仅适用于标量值**。对于每个 Track，从每个网格级别提取一个标量结果（如 $V_{\text{end}}$、$T_{\max}$、$E_{\text{frac}}(t_{\text{end}})$），然后用这些标量值计算 GCI。

**L2 范数适用于曲线和场量**。对于时间序列（如 $V(t)$）和空间场（如 $T(r, \theta)$），计算 L2 误差范数量化整体差异。

两者不可混用：不要对 L2 范数值序列计算 GCI，也不要对曲线做 GCI 计算。

### 3.5 渐近收敛检查

$$\text{Ratio} = \frac{\text{GCI}_{23}}{r_{21}^p \cdot \text{GCI}_{12}}$$

若 Ratio $\in [0.8, 1.2]$，表明解在渐近收敛区间内，GCI 值可靠。

---

## 4. 各 Track 指标定义

### 4.1 电化学 Track (Script 2)

**网格变量**：`{Nn, Ns, Np}` 四级

| 配置 | 标签 |
|------|------|
| (40, 20, 40) | 参考解（最细） |
| (20, 10, 20) | Level 2 |
| (20, 5, 20) | Level 3 |
| (10, 5, 10) | Level 4（最粗） |

**收敛指标**：

| 指标 | 范数类型 | 物理量 | GCI 标量 |
|------|----------|--------|----------|
| 电压曲线误差 | $\|e\|_{L_2,\text{rel}}$ | `cell voltage [V]` 时间序列 | $V_{\text{end}}$ |
| 温度曲线误差 | $\|e\|_{L_2,\text{rel}}$ | `temperature [K]` 时间序列 | $T_{\text{peak}}$ |
| dT/dt 曲线误差 | $\|e\|_{L_2}$（绝对值） | dT/dt 差分序列 | $\max|dT/dt|$ |

**时间序列对齐**：不同网格产生不同时间步序列。使用线性插值将所有候选解对齐到参考解的时间网格上（保留原版 `align_to_ref` 逻辑）。

**理论收敛阶**：对于 Crank-Nicolson 时间离散 + 1D FEM 空间离散（二阶），L2 收敛阶 $O(h^2)$。

### 4.2 热学 Track (Script 3)

**网格变量**：`nθ ∈ {20, 40, 80, 160}`

**物理基础**：均匀体积热源 + 表面冷却的纯热模型。

**收敛指标**：

| 指标 | 范数类型 | 物理量 | GCI 标量 |
|------|----------|--------|----------|
| 空间温度场 L2 | $\|e\|_{L_2,\text{rel}}$（时间平均） | T_hist (nnode × nt) | $T_{\max}(t_{\text{end}})$ |
| T_max(t) 曲线误差 | $\|e\|_{L_2,\text{rel}}$ | $T_{\max}(t)$ 曲线 | — |
| T_range(t) 曲线误差 | $\|e\|_{L_2}$（绝对值） | $\Delta T(t)$ 曲线 | — |

**空间场插值策略**：不同 nθ 的网格节点数不同。需将粗网格场插值到细网格参考节点上：

1. 将所有节点坐标从极坐标 $(r, \theta)$ 转换为笛卡尔坐标 $(x, y)$：$x = r\cos\theta$, $y = r\sin\theta$
2. 取最细网格 (nθ=160) 的笛卡尔节点坐标作为参考
3. 对每级粗网格，使用反距离加权 (IDW) 插值到参考节点
4. 在插值后的场上计算 L2 范数

**使用笛卡尔坐标的原因**：
- 极坐标下 $\theta$ 具有周期性（$\theta=0$ 和 $\theta=2\pi$ 物理上相邻），直接在 $(r, \theta)$ 空间计算距离会导致边界附近邻居搜索错误
- $r$（米）和 $\theta$（弧度）量纲不匹配，会导致距离计算偏向 $\theta$ 方向
- 笛卡尔坐标天然消除了周期性和量纲问题

IDW 公式：

$$\hat{T}(x, y) = \frac{\sum_j w_j T_j}{\sum_j w_j}, \quad w_j = \frac{1}{d_j^2}$$

其中 $d_j = \sqrt{(x - x_j)^2 + (y - y_j)^2}$ 是目标点到第 $j$ 个源节点的笛卡尔距离，搜索半径取最近 4 个节点。

**理论收敛阶**：本项目使用 Q4（4 节点双线性四边形单元），多项式阶数 $k=1$。L2 收敛阶 $O(h^2)$，H1 收敛阶 $O(h^1)$。

**注**：参考文献 Ai & Liu (2023) 的 "2nd order FEM" 指的是方法达到二阶收敛率 $O(h^2)$，而非使用二次单元（Q8/Q9）。

### 4.3 CZM Track (Script 4)

**网格变量**：由 cohesive 特征长度 $l_c$ 确定的 nθ 四级

**物理基础**：全耦合模型（SPMe + distributed2D thermal + CZM），热-化学载荷驱动损伤。

**收敛指标**：

| 指标 | 范数类型 | 物理量 | GCI 标量 | 优先级 |
|------|----------|--------|----------|--------|
| **断裂能耗散** | — | $E_{\text{frac}}(t) = \sum_e G_c \cdot l_e \cdot D_e(t)$ | $E_{\text{frac}}(t_{\text{end}})$ | 主 |
| D_max(t) 曲线误差 | $\|e\|_{L_2,\text{rel}}$ | `czm D_max` 时间序列 | $D_{\max}(t_{\text{end}})$ | 辅 |
| 断裂数曲线误差 | $\|e\|_{L_2}$（绝对值） | `czm n_fractured` 时间序列 | $n_f(t_{\text{end}})$ | 辅 |
| 牵引-分离面积偏差 | $\epsilon_A$（保留） | 峰值损伤单元的 T-δ 曲线 | — | 辅 |

**断裂能耗散计算**：

$$E_{\text{frac}}(t) = \sum_{e=1}^{N_{\text{coh}}} G_c \cdot l_e \cdot D_e(t)$$

其中 $l_e$ 是第 $e$ 个 cohesive 单元的长度，$D_e(t)$ 是其损伤值。这是积分量，单调收敛，不受损伤前沿位置偏移影响。

**单位说明**：$G_c$ [J/m²] × $l_e$ [m] × $D_e$ [-] = [J/m]，即单位深度（2D 模型假设平面应变）的断裂能。不同网格级别的 $E_{\text{frac}}$ 比较时使用相对值（GCI），单位不影响收敛性分析。

**牵引-分离面积偏差的单元选择策略**（保留旧设计）：选取 $D_{\max}$ 达到最大值的单元的 traction-separation 时间序列。

### 4.4 能量守恒检查 (Script 5)

保持独立脚本，单网格 (nθ=80) 运行。

**简化方案**：

$$R(t) = Q_{\text{gen}}(t) - \Delta E_{\text{th}}(t) - Q_{\text{loss}}(t) - E_{\text{frac}}(t)$$

**指标**：

| 指标 | 公式 | 说明 |
|------|------|------|
| 瞬时相对误差 | $\varepsilon_R(t) = |R(t)| / |Q_{\text{gen}}(t)| \times 100\%$ | 保留 |
| 归一化 RMS 残余 | $\varepsilon_{R,\text{rms}} = \sqrt{\frac{1}{N}\sum R(t_i)^2} / |W_{\text{elec}}(t_{\text{end}})| \times 100\%$ | 保留 |

**说明**：能量守恒检查是物理验证，不是网格收敛性指标。标注清楚两者关系。

---

## 5. 工具函数接口

所有脚本共用 `0_convergence_utils.jl`，需 `using Statistics, LinearAlgebra`。

### 5.1 离散范数函数

```julia
"""
    l2_norm(f, f_ref) -> Float64

绝对 L2 范数：sqrt(mean((f - f_ref)^2))
"""
function l2_norm(f, f_ref)
    return sqrt(mean((f .- f_ref).^2))
end

"""
    l2_rel_norm(f, f_ref) -> Float64

归一化 L2 范数：||f - f_ref||_L2 / ||f_ref||_L2
分母 ||f_ref||_L2 = sqrt(mean(f_ref^2))，对非零参考场永远有定义。
"""
function l2_rel_norm(f, f_ref)
    norm_ref = sqrt(mean(f_ref.^2))
    norm_ref == 0 && return NaN
    return l2_norm(f, f_ref) / norm_ref
end

"""
    max_norm(f, f_ref) -> Float64

L∞ 最大范数：max|f - f_ref|
"""
function max_norm(f, f_ref)
    return maximum(abs.(f .- f_ref))
end
```

### 5.2 GCI 框架函数

```julia
"""
    observed_order(f1, f2, f3, r21, r32) -> Float64

计算观测收敛阶 p，使用 Celik et al. (2008) 的迭代法处理非恒定细化比。
f1: 最细网格标量值, f2: 中间, f3: 最粗
r21 = h2/h1, r32 = h3/h2

迭代方程：
  p_{k+1} = |ln|alpha| + ln((r21^p_k - 1)/(r32^p_k - 1))| / ln(r21)
其中 alpha = (f3 - f2) / (f2 - f1)
"""
function observed_order(f1, f2, f3, r21, r32)
    denom = f2 - f1
    numer = f3 - f2
    if abs(denom) < 1e-15 || abs(numer) < 1e-15
        return NaN  # 无法确定收敛阶
    end
    s = sign(denom / numer)  # 检查单调收敛
    if s < 0
        return NaN  # 非单调收敛，GCI 不适用
    end

    alpha = numer / denom

    # 初始值：恒定比近似
    p = abs(log(abs(alpha)) / log(r21))

    # 固定点迭代（Celik et al. 2008）
    for _ in 1:100
        r21p = r21^p
        r32p = r32^p
        if r21p - 1 < 1e-15 || r32p - 1 < 1e-15
            break
        end
        correction = log((r21p - 1) / (r32p - 1))
        p_new = abs(log(abs(alpha)) + correction) / log(r21)
        if abs(p_new - p) < 1e-6
            return p_new
        end
        p = p_new
        if p > 20 || isnan(p)
            # 发散，回退到恒定比近似
            return abs(log(abs(alpha)) / log(r21))
        end
    end
    return p
end

"""
    compute_gci(f_fine, f_coarse, r; p, Fs=1.25) -> Float64

计算 Grid Convergence Index [%]。
f_fine: 细网格标量值, f_coarse: 粗网格标量值
r: 细化比 = h_coarse / h_fine
p: 收敛阶
Fs: 安全系数（3+ 级网格用 1.25）
"""
function compute_gci(f_fine, f_coarse, r; p, Fs=1.25)
    abs(f_fine) < 1e-15 && return NaN
    epsilon = (f_coarse - f_fine) / f_fine
    denom = r^p - 1
    denom <= 0 && return NaN
    return Fs * abs(epsilon) / denom * 100
end

"""
    asymptotic_check(gci_12, gci_23, r21, p) -> Float64

渐近收敛检查。返回 GCI_23 / (r21^p * GCI_12)。
比值接近 1.0 表明解在渐近收敛区间内。
"""
function asymptotic_check(gci_12, gci_23, r21, p)
    gci_12 <= 0 && return NaN
    return gci_23 / (r21^p * gci_12)
end

"""
    effective_h(mesh) -> Float64

计算等效单元尺寸 h = sqrt(A_total / N_elem)。
接收 thermal2D mesh 对象（case.mesh["thermal2D"]），需有 .area 和 .element 字段。
对于电化学 Track（无热网格），h 可从 (1/Nn + 1/Np) 的倒数近似。
"""
function effective_h(mesh)
    A_total = sum(mesh.area)
    N_elem = size(mesh.element, 1)
    return sqrt(A_total / N_elem)
end
```

### 5.3 场对齐与插值函数

```julia
"""
    align_to_ref(t_cand, y_cand, t_ref)

将候选解插值到参考解的时间网格上，返回对齐后的 y_cand_aligned。
线性插值，超出范围用端点值填充。
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

"""
    interpolate_to_ref_field(coarse_vals, coarse_r, coarse_theta,
                             ref_r, ref_theta; k=4)

使用反距离加权 (IDW) 将粗网格场插值到参考节点上。
所有距离计算在笛卡尔坐标下进行，以正确处理 θ 的周期性和 r-θ 量纲不匹配。

coarse_vals: 粗网格节点值 (n_coarse,)
coarse_r, coarse_theta: 粗网格节点极坐标
ref_r, ref_theta: 参考网格节点极坐标
k: 最近邻数量
"""
function interpolate_to_ref_field(coarse_vals, coarse_r, coarse_theta,
                                  ref_r, ref_theta; k=4)
    # 极坐标 → 笛卡尔坐标
    coarse_x = coarse_r .* cos.(coarse_theta)
    coarse_y = coarse_r .* sin.(coarse_theta)
    ref_x = ref_r .* cos.(ref_theta)
    ref_y = ref_r .* sin.(ref_theta)

    n_ref = length(ref_r)
    n_coarse = length(coarse_r)
    result = similar(coarse_vals, n_ref)
    k_eff = min(k, n_coarse)

    for i in 1:n_ref
        # 笛卡尔距离
        dists = sqrt.(
            (coarse_x .- ref_x[i]).^2 .+
            (coarse_y .- ref_y[i]).^2
        )
        # 取最近 k 个
        idx = partialsortperm(dists, 1:k_eff)
        w = 1.0 ./ (dists[idx].^2 .+ 1e-30)  # 加小量防止除零
        result[i] = sum(w .* coarse_vals[idx]) / sum(w)
    end
    return result
end
```

### 5.4 辅助函数（保留）

```julia
"""
    trapz(x, y)

梯形积分。
"""
function trapz(x, y)
    return sum(0.5 .* diff(x) .* (y[2:end] .+ y[1:end-1]))
end

"""
    area_error(x, y, x_ref, y_ref)

归一化曲线面积偏差 [%]。
"""
function area_error(x, y, x_ref, y_ref)
    length(x) < 2 && return NaN
    length(x_ref) < 2 && return NaN
    A  = abs(trapz(x, y))
    Ar = abs(trapz(x_ref, y_ref))
    Ar == 0 && return NaN
    return abs(A - Ar) / Ar * 100
end
```

---

## 6. 脚本结构

### 6.1 目录布局

```
example/网格敏感性_v2/
├── 0_convergence_utils.jl          # 工具函数
├── 1_cohesive_characteristic_length.jl  # 从原版复制
├── 2_electrochemical_convergence.jl     # Script 2
├── 3_thermal_convergence.jl             # Script 3
├── 4_czm_convergence.jl                 # Script 4
└── 5_energy_conservation.jl             # Script 5
```

### 6.2 通用收敛分析流程（每个 Script）

```
1. 运行 4 级网格 → 收集结果
2. 计算每级网格的等效 h (effective_h)
3. 对每对标量值 (f_fine, f_coarse)：
   a. 计算细化比 r = h_coarse / h_fine
   b. 用三组值计算观测收敛阶 p (observed_order)
   c. 计算 GCI_fine = Fs * |epsilon| / (r^p - 1) * 100%
4. 对曲线/场量：
   a. 时间序列对齐 (align_to_ref)
   b. 空间场插值 (interpolate_to_ref_field)，仅热学 Track
   c. 计算 L2_rel / L2 / L∞ 范数
5. 渐近收敛检查
6. 输出：
   a. 收敛误差表 (h, L2_rel, L∞, p_obs, GCI, asymptotic ratio)
   b. GCI 汇总表 (level pair, r, epsilon, p, GCI_fine)
   c. 收敛图 (log-log，含理论斜率线)
   d. 物理量曲线对比图
```

### 6.3 收敛图标准格式

所有收敛图遵循统一格式：

- **x 轴**：$\log(h)$，等效单元尺寸
- **y 轴**：$\log(\|e\|_{L_2})$ 或 $\log(\|e\|_{L_2,\text{rel}})$
- **叠加**：理论收敛阶斜率线 $O(h^p)$
- **标注**：观测收敛阶 $p_{\text{obs}}$，与理论值 $p_{\text{theory}}$ 的比较

**理论收敛阶参考**（本项目单元类型与离散方法）：

| 物理场 | 单元类型 | L2 理论阶 | H1 理论阶 | 说明 |
|--------|----------|-----------|-----------|------|
| 2D 热场 | Q4 ($k=1$) | $O(h^2)$ | $O(h^1)$ | 标准 FEM：$L_2 = O(h^{k+1})$ |
| 1D 电化学 | 线性 FEM + CN | $O(h^2)$ | $O(h^1)$ | Crank-Nicolson 时间 + FEM 空间 |
| CZM | 依赖问题 | 观测值为准 | — | 无通用理论，以观测收敛阶为准 |

### 6.4 输出目录

```
output/mesh_convergence/
├── echem_voltage_convergence.png
├── echem_temperature_convergence.png
├── echem_error_convergence.png      # log-log 收敛图
├── thermal_temperature_convergence.png
├── thermal_spatial_error.png        # log-log 收敛图
├── thermal_peak_convergence.png
├── czm_damage_evolution.png
├── czm_fracture_energy.png          # 新增：断裂能耗散图
├── czm_error_convergence.png        # log-log 收敛图
├── energy_components.png
├── energy_residual.png
└── energy_relative_error.png
```

---

## 7. 验收标准

### 7.1 网格收敛性

| 指标类型 | 阈值 | 说明 |
|----------|------|------|
| GCI_fine | < 5% | 网格不确定性量化 |
| 观测收敛阶 $p$ | 接近理论值 | 确认数值方法实现正确 |
| 渐近检查比值 | $\in [0.8, 1.2]$ | 解在渐近收敛区间内 |
| L2_rel | 单调递减（随 $h$ 减小） | 收敛行为合理 |

### 7.2 物理验证

| 指标 | 阈值 | 说明 |
|------|------|------|
| 能量残余 $\varepsilon_R$ | < 1% | 物理守恒 |
| 归一化 RMS 残余 $\varepsilon_{R,\text{rms}}$ | < 5% | 整体能平衡 |

### 7.3 报告要求

每个 Track 需输出：

1. **收敛误差表**：每级网格的 $h$, $\|e\|_{L_2,\text{rel}}$, $\|e\|_{\infty}$, $p_{\text{obs}}$, GCI, 渐近比值
2. **GCI 汇总表**：每对相邻级别的 $r$, $\varepsilon$, $p$, GCI_fine
3. **收敛图**：log-log 误差图（含理论斜率线）
4. **物理量曲线对比图**
5. **结论**：推荐网格参数及对应的 GCI 值

---

## 8. 与旧规格的映射关系

| 旧指标 (2026-04-29) | 新指标 (本文档) | 变化说明 |
|---------------------|-----------------|----------|
| `rmspe(y, y_ref)` | `l2_rel_norm(f, f_ref)` | 范数替换 RMSPE，消除零点除法 |
| `spatial_rmspe_over_time()` | 空间场 L2_rel + IDW 插值 | 增加跨网格插值，不再跳过 |
| `area_error()` | `area_error()`（保留） | 不变 |
| `align_to_ref()` | `align_to_ref()`（保留） | 不变 |
| V(t) RMSPE% | V(t) L2_rel% + GCI on $V_{\text{end}}$ | 增加 GCI |
| T(t) RMSPE% | T(t) L2_rel% + GCI on $T_{\text{peak}}$ | 增加 GCI |
| dT/dt RMSPE% | dT/dt L2（绝对值） | 去掉百分比，保留绝对值 |
| Spatial RMSPE% | 空间场 L2_rel% + GCI on $T_{\max}$ | 增加 IDW 插值和 GCI |
| T_range(t) RMSPE% | T_range(t) L2（绝对值） | 去掉百分比 |
| D_max(t) RMSPE% | D_max(t) L2_rel% + **断裂能耗散 GCI** | 断裂能升为主指标 |
| n_frac(t) RMSPE% | n_frac(t) L2（绝对值） + GCI | 增加 GCI |
| δ_max_n(t) RMSPE% | δ_max_n(t) L2_rel% | 范数替换 |
| 面积偏差% | 面积偏差%（保留） | 不变 |
| — | 观测收敛阶 $p$ | **新增** |
| — | GCI [%] | **新增** |
| — | 渐近收敛检查 | **新增** |
| — | log-log 收敛图 + 理论斜率线 | **新增** |

---

## 9. 不变的部分

以下内容不受本次重新设计影响：

- 仿真配置（模型选择、电流函数、边界条件）
- 物理判据（Biot 数分析、cohesive 特征长度、$l_c$ 计算）
- 网格候选值（电化学 (40,20,40) 等、热学 nθ={20,40,80,160}、CZM nθ 由 $l_c$ 决定）
- Script 1（cohesive 特征长度计算脚本）
- 时间步策略（自适应 dtType="auto"）
- 高斯积分阶（gsorder=2）
- 能量守恒检查的物理公式

---

## 10. 参考文献

1. Roache, P.J. (1998). *Verification and Validation in Computational Science and Engineering*. Hermosa Publishers.
2. Celik, I.B. et al. (2008). "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications." *ASME J. Fluids Eng.*, 130(7), 078001.
3. Schwer, L. (2008). "Is Your Mesh Refined Enough? Estimating Discretization Error via GCI." *DynaMORE Forum*.
4. ASME V&V 10-2019. *Standard for Verification and Validation in Computational Solid Mechanics*.
5. Turon, A. et al. (2007). "An engineering solution for mesh size effects in the simulation of delamination using cohesive zone models." *Eng. Fract. Mech.*, 74(10), 1665-1682.
6. Riedelheimer et al. (2014). "Application of Grid Convergence Index in FE Computation."
7. NASA GRC. "Examining Spatial (Grid) Convergence." https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
8. Ai, W., Liu, Y. (2023). "Improving the convergence rate of Newman's battery model using 2nd order finite element method." *J. Energy Storage*, 67, 107512.
