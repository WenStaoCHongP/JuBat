# 对流边界条件正负号最终验证

## 代码约定确认

### 第132行的约定

```julia
/**
 * 重要约定：返回的 KT 已包含负号，即 M dT/dt = KT T + F
 */
```

**标准形式**：
$$
M \frac{dT}{dt} + K_{\text{standard}} T = F
$$

**代码约定**：
$$
M \frac{dT}{dt} = KT \cdot T + F
$$

**转换关系**：
$$
KT = -K_{\text{standard}}
$$

## 理论推导验证

### 1. 标准热传导方程

物理方程：
$$
\rho c \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + q
$$

等价于：
$$
\rho c \frac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = q
$$

### 2. 对流边界条件

在边界 $\Gamma$ 上：
$$
-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})
$$

### 3. 弱形式推导

使用测试函数 $w$，在 $\Omega$ 上积分：
$$
\int_\Omega \rho c \frac{\partial T}{\partial t} w \, d\Omega - \int_\Omega \nabla \cdot (k \nabla T) w \, d\Omega = \int_\Omega q w \, d\Omega
$$

分部积分（Green公式）：
$$
\int_\Omega \rho c \frac{\partial T}{\partial t} w \, d\Omega + \int_\Omega k \nabla T \cdot \nabla w \, d\Omega - \int_{\Gamma} k \frac{\partial T}{\partial n} w \, d\Gamma = \int_\Omega q w \, d\Omega
$$

代入边界条件 $-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})$：
$$
\int_\Omega \rho c \frac{\partial T}{\partial t} w \, d\Omega + \int_\Omega k \nabla T \cdot \nabla w \, d\Omega + \int_{\Gamma} h(T - T_{\text{amb}}) w \, d\Gamma = \int_\Omega q w \, d\Omega
$$

展开对流项：
$$
\int_\Omega \rho c \frac{\partial T}{\partial t} w \, d\Omega + \int_\Omega k \nabla T \cdot \nabla w \, d\Omega + \int_{\Gamma} h T w \, d\Gamma - \int_{\Gamma} h T_{\text{amb}} w \, d\Gamma = \int_\Omega q w \, d\Omega
$$

移项整理：
$$
\int_\Omega \rho c \frac{\partial T}{\partial t} w \, d\Omega + \left[\int_\Omega k \nabla T \cdot \nabla w \, d\Omega + \int_{\Gamma} h T w \, d\Gamma\right] = \int_\Omega q w \, d\Omega + \int_{\Gamma} h T_{\text{amb}} w \, d\Gamma
$$

### 4. 离散化

$$
M \dot{T} + K_{\text{standard}} T = F_{\text{total}}
$$

其中：
- **质量矩阵**：$M = \int \rho c N_i N_j \, d\Omega$
- **导热刚度**：$K_{\text{cond}} = \int k \nabla N_i \cdot \nabla N_j \, d\Omega$（**正号**）
- **对流刚度**：$K_{\text{conv}} = \int_\Gamma h N_i N_j \, d\Gamma$（**正号**）
- **标准刚度**：$K_{\text{standard}} = K_{\text{cond}} + K_{\text{conv}}$（**正号**）

### 5. 代码约定转换

代码约定：$M \dot{T} = KT \cdot T + F$

对比标准形式：$M \dot{T} + K_{\text{standard}} T = F$

得到：
$$
KT = -K_{\text{standard}} = -(K_{\text{cond}} + K_{\text{conv}})
$$

因此：
- **导热项**（代码）：$KT_{\text{cond}} = -K_{\text{cond}} = -\int k \nabla N \cdot \nabla N$（**负号**）
- **对流项**（代码）：$KT_{\text{conv}} = -K_{\text{conv}} = -\int h N N$（**负号**）

**结论**：按代码约定，所有刚度矩阵项都应该是**负号** ✅

## 代码实现检查

### 1. 导热刚度矩阵（基准，第291行）

```julia
weights = -λ_iso_nd .* (wJ ./ L_th^2)  // 权重为负
KT_x = Assemble(Vi, Vj, dNdx, dNdx, weights, nnode)
KT_y = Assemble(Vi, Vj, dNdy, dNdy, weights, nnode)
```

**系数分析**：
- $\lambda_{\text{iso}} > 0$（热导率，正数）
- $wJ > 0$（权重×雅可比，正数）
- $weights = -\lambda \cdot wJ < 0$（**负数**）

**刚度矩阵**：
$$
KT_{\text{cond}}[i,j] = \sum weights \cdot \nabla N_i \cdot \nabla N_j < 0
$$

**符号**：**负号** ✅

### 2. 外圈对流边界（第396-401行）

```julia
for (s, w) in zip(s_vals, w_vals)
    N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
    wt = Bi * w * (J / L_th)  // 正数
    
    // 加负号与体内扩散项统一
    ke11 += -wt * N1 * N1  // 计算：-（正数）= 负数
    ke12 += -wt * N1 * N2  // 计算：-（正数）= 负数
    ke22 += -wt * N2 * N2  // 计算：-（正数）= 负数
    fe1 += wt * T_amb * N1  // 正数
    fe2 += wt * T_amb * N2  // 正数
end

KT[a, a] += ke11  // 装配：KT += 负数
```

**系数分析**：
- $Bi = h \cdot L / k > 0$（Biot数，正数）
- $w > 0$（高斯权重，正数）
- $J = L/2 > 0$（雅可比，正数）
- $wt = Bi \cdot w \cdot J / L_{th} > 0$（**正数**）
- $ke11 = -wt \cdot N_1 \cdot N_1 < 0$（**负数**）

**刚度矩阵**：
$$
KT[a,a] += ke11 < 0
$$

**符号**：**负号** ✅

### 3. 表面冷却（修正后，第547行）

```julia
for g in 1:ngs
    wt = conv_factor * wJ[g]  // 正数
    for i, j in element_nodes
        KT[ni, nj] -= wt * Ni_g * Nj_g  // 装配：KT -= 正数
    end
end
```

**系数分析**：
- $conv\_factor = Bi_z / L_{th}^2 > 0$（正数）
- $wJ[g] > 0$（高斯权重×雅可比，正数）
- $wt = conv\_factor \cdot wJ > 0$（**正数**）
- $Ni_g, Nj_g > 0$（形函数值，正数）

**刚度矩阵**：
$$
KT[ni, nj] -= wt \cdot N_i \cdot N_j
$$

即：
$$
KT[ni, nj] = KT[ni, nj] - wt \cdot N_i \cdot N_j < KT_{\text{原}}
$$

**效果**：减去一个正数，KT变得更负

**符号**：**负号** ✅

### 4. 极耳冷却（修正后，第630行）

```julia
for (i, n) in enumerate(tab_nodes)
    weight = arc_lengths[i] / total_arc_length  // 正数
    coeff = h_tab * tab_area * weight / (H * k_th * L_th)  // 正数
    
    KT[n, n] -= coeff  // 装配：KT -= 正数
end
```

**系数分析**：
- $h_{\text{tab}} > 0$（对流系数，正数）
- $A_{\text{tab}} > 0$（极耳面积，正数）
- $weight > 0$（弧长权重，正数，且 $\sum weight = 1$）
- $H, k_{th}, L_{th} > 0$（正数）
- $coeff = \frac{h_{\text{tab}} \cdot A_{\text{tab}} \cdot weight}{H \cdot k_{th} \cdot L_{th}} > 0$（**正数**）

**刚度矩阵**：
$$
KT[n, n] -= coeff
$$

即：
$$
KT[n, n] = KT[n, n] - coeff < KT_{\text{原}}
$$

**效果**：减去一个正数，KT变得更负

**符号**：**负号** ✅

## 符号一致性总结

| 项 | 代码操作 | 系数符号 | 刚度矩阵变化 | 最终符号 | 位置 |
|---|---------|---------|------------|---------|------|
| **导热** | `KT += weights * ...` | `weights < 0` | 减小 | **负** ✅ | L291 |
| **外圈对流** | `KT += ke` | `ke < 0` | 减小 | **负** ✅ | L405 |
| **表面冷却** | `KT -= wt * ...` | `wt > 0` | 减小 | **负** ✅ | L547 |
| **极耳冷却** | `KT -= coeff` | `coeff > 0` | 减小 | **负** ✅ | L630 |

**结论**：所有项的最终符号**完全一致**，都是**负号** ✅

## 物理意义验证

### 能量方程

代码形式：
$$
M \frac{dT}{dt} = KT \cdot T + F
$$

其中 $KT < 0$（对角元素为负）。

移项：
$$
M \frac{dT}{dt} - KT \cdot T = F
$$

即：
$$
M \frac{dT}{dt} + (-KT) \cdot T = F
$$

因为 $KT < 0$，所以 $-KT > 0$。

### 稳态分析

稳态时 $\dot{T} = 0$：
$$
KT \cdot T = -F
$$

展开（无内热源，$F = F_{\text{conv}}$）：
$$
KT \cdot T = -F_{\text{conv}}
$$

即：
$$
(KT_{\text{cond}} + KT_{\text{conv}}) T = -F_{\text{conv}}
$$

因为 $KT < 0$（负号），$F_{\text{conv}} > 0$（正号），所以：
$$
(-|KT|) T = -F_{\text{conv}}
$$

$$
|KT| T = F_{\text{conv}}
$$

**物理意义**：
- 左边：散热（导热+对流）
- 右边：环境温度驱动

如果 $T > T_{\text{amb}}$：
- $|KT| T > F_{\text{conv}}$（散热大于驱动）
- 温度下降，趋向 $T_{\text{amb}}$ ✅

### 温度响应

假设：
- 初始温度：$T_0 = 350$ K（高于环境）
- 环境温度：$T_{\text{amb}} = 298$ K
- 无内热源

**瞬态方程**：
$$
M \frac{dT}{dt} = KT \cdot (T - T_{\text{amb}})
$$

因为 $T > T_{\text{amb}}$ 且 $KT < 0$：
$$
\frac{dT}{dt} = \frac{KT}{M} (T - T_{\text{amb}}) < 0
$$

**结论**：温度下降 ✅

**最终状态**：$T \to T_{\text{amb}}$ ✅

## 数值稳定性

### 对角占优

刚度矩阵对角元素：
$$
KT_{ii} = KT_{\text{cond},ii} + KT_{\text{conv},ii} + KT_{\text{z-cool},ii}
$$

所有项都是负数：
$$
KT_{ii} < 0
$$

**效果**：
- 矩阵 $-KT$ 的对角元素为正
- 对角占优增强
- 数值稳定性提高 ✅

### 系统矩阵

时间离散（如隐式Euler）：
$$
\left(\frac{M}{\Delta t} - KT\right) T^{n+1} = \frac{M}{\Delta t} T^n + F
$$

系统矩阵：
$$
A = \frac{M}{\Delta t} - KT
$$

因为 $M > 0$（对角）且 $KT < 0$（对角元素负）：
$$
A_{ii} = \frac{M_{ii}}{\Delta t} - KT_{ii} > 0
$$

**效果**：系统矩阵对角占优，数值稳定 ✅

## 与外圈对流的完全一致性

### 实现方式对比

**外圈对流**（边界线积分）：
```julia
wt = Bi * w * J / L_th  // 正数
ke11 = -wt * N1 * N1     // 负数
KT[a,a] += ke11          // 装配负数 → KT减小
```

**表面冷却**（xy平面面积分）：
```julia
wt = Bi_z * wJ / L_th^2  // 正数
KT[ni,nj] -= wt * Ni * Nj  // 减去正数 → KT减小
```

**极耳冷却**（节点集中）：
```julia
coeff = h * A * w / (H * k * L)  // 正数
KT[n,n] -= coeff                  // 减去正数 → KT减小
```

### 等价性验证

三种方式都是：
1. 计算一个**正数**系数（wt, coeff）
2. 使其对KT的贡献为**负**（通过 `-wt` 或 `-=`）
3. 最终效果：KT变得更负

**完全一致** ✅

## 最终结论

### 符号正确性

✅ **导热项**：负号（代码约定）  
✅ **外圈对流**：负号（已有实现）  
✅ **表面冷却**：负号（已修正）  
✅ **极耳冷却**：负号（已修正）  

### 物理正确性

✅ **能量平衡**：散热项正确  
✅ **温度响应**：高温冷却到环境温度  
✅ **稳态解**：收敛到 $T_{\text{amb}}$  

### 数值正确性

✅ **矩阵性质**：对角占优  
✅ **系统稳定**：数值收敛  
✅ **一致性**：所有对流项符号统一  

---

**最终确认**：符号问题已**完全正确** ✅
