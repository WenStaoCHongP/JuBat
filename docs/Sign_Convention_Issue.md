# 对流边界条件正负号不一致问题

## 问题发现

检查三种冷却方式的符号，发现**不一致**：

### 1. 外圈对流边界（第396-401行）

```julia
for (s, w) in zip(s_vals, w_vals)
    N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
    wt = Bi * w * (J / L_th)
    
    // 加负号与体内扩散项统一
    ke11 += -wt * N1 * N1  // ❌ 负号
    ke12 += -wt * N1 * N2  // ❌ 负号
    ke22 += -wt * N2 * N2  // ❌ 负号
    fe1 += wt * T_amb * N1  // ✅ 正号
    fe2 += wt * T_amb * N2  // ✅ 正号
end

KT[a, a] += ke11  // 负值
KT[a, b] += ke12  // 负值
```

**刚度矩阵**：**负号** ❌  
**载荷向量**：正号 ✅

### 2. 表面冷却（第547行）

```julia
KT[ni, nj] += wt * Ni_g * Nj_g  // ✅ 正号
FT[ni] += wt * T_amb_nd * Ni_g  // ✅ 正号
```

**刚度矩阵**：**正号** ✅  
**载荷向量**：正号 ✅

### 3. 极耳冷却（第630行）

```julia
KT[n, n] += coeff              // ✅ 正号
FT[n] += coeff * T_amb_nd      // ✅ 正号
```

**刚度矩阵**：**正号** ✅  
**载荷向量**：正号 ✅

## 符号不一致！

| 冷却方式 | KT符号 | FT符号 | 代码位置 |
|---------|--------|--------|---------|
| 外圈对流 | **负号** ❌ | 正号 | L396-401 |
| 表面冷却 | **正号** ✅ | 正号 | L547 |
| 极耳冷却 | **正号** ✅ | 正号 | L630 |

**问题**：外圈对流的刚度矩阵用了负号，但表面/极耳冷却用了正号！

## 代码约定分析

### 约定说明（第132行）

```julia
**重要约定**：返回的 KT 已包含负号，即 M dT/dt = KT T + F
```

### 约定解释

标准FEM形式：
$$
M \frac{dT}{dt} + K_{\text{standard}} T = F
$$

代码约定形式：
$$
M \frac{dT}{dt} = KT \cdot T + F
$$

对比得：
$$
KT = -K_{\text{standard}}
$$

**含义**：代码中的 `KT` 是**负的**标准刚度矩阵。

## 理论推导

### 标准热传导方程

$$
\rho c \frac{\partial T}{\partial t} + \nabla \cdot (-k \nabla T) = q
$$

其中 $\nabla \cdot (-k \nabla T) = -\nabla \cdot (k \nabla T)$

标准形式：
$$
\rho c \frac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = q
$$

### 对流边界条件

边界上：
$$
-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})
$$

### 弱形式（分部积分后）

$$
M \dot{T} + \underbrace{\int_\Omega k \nabla T \cdot \nabla w \, d\Omega}_{K_{\text{cond}} T} + \underbrace{\int_{\Gamma} h T w \, d\Gamma}_{K_{\text{conv}} T} = \int_\Omega q w \, d\Omega + \int_{\Gamma} h T_{\text{amb}} w \, d\Gamma
$$

标准刚度矩阵：
$$
K_{\text{standard}} = K_{\text{cond}} + K_{\text{conv}} = \int k \nabla N \cdot \nabla N + \int_\Gamma h N N
$$

### 代码约定转换

$$
KT = -K_{\text{standard}} = -K_{\text{cond}} - K_{\text{conv}}
$$

因此：
- 导热项（代码）：$KT_{\text{cond}} = -\int k \nabla N \cdot \nabla N$（**负号**）
- 对流项（代码）：$KT_{\text{conv}} = -\int h N N$（**负号**）

**结论**：按代码约定，所有刚度矩阵贡献都应该是**负号**！

## z方向冷却的符号

### 物理方程

电池上下表面冷却（z方向）：
$$
\rho c \frac{\partial T}{\partial t} = \nabla_{xy} \cdot (k \nabla_{xy} T) + q - \frac{2h}{H}(T - T_{\text{amb}})
$$

整理：
$$
\rho c \frac{\partial T}{\partial t} + \left[-\nabla_{xy} \cdot (k \nabla_{xy} T) + \frac{2h}{H} T\right] = q + \frac{2h}{H} T_{\text{amb}}
$$

标准形式：
$$
M \dot{T} + K_{\text{cond}} T + K_{\text{z-cool}} T = F + F_{\text{z-cool}}
$$

其中：
- $K_{\text{cond}} = \int k \nabla N \cdot \nabla N$（正号）
- $K_{\text{z-cool}} = \int \frac{2h}{H} N N$（正号）

### 代码约定

$$
KT = -(K_{\text{cond}} + K_{\text{z-cool}})
$$

因此：
- $KT_{\text{cond}} = -\int k \nabla N \cdot \nabla N$（**负号**）
- $KT_{\text{z-cool}} = -\int \frac{2h}{H} N N$（**负号**）

**结论**：z方向冷却也应该是**负号**！

## 当前代码的问题

### 正确的实现

✅ **外圈对流边界**（L396-401）：
```julia
ke11 += -wt * N1 * N1  // ✅ 负号，符合约定
```

❌ **表面冷却**（L547）：
```julia
KT[ni, nj] += wt * Ni_g * Nj_g  // ❌ 应该是 -= 或 wt 本身为负
```

❌ **极耳冷却**（L630）：
```julia
KT[n, n] += coeff  // ❌ 应该是 -= 或 coeff 本身为负
```

## 错误的后果

### 物理意义

标准形式：$M \dot{T} + K T = F$
- $K > 0$：系统稳定（散热）
- $K < 0$：系统不稳定（增热）

代码形式：$M \dot{T} = KT \cdot T + F$
- $KT < 0$：系统稳定（对应 $K > 0$）
- $KT > 0$：系统不稳定（对应 $K < 0$）

**如果表面/极耳冷却用正号**：
- 相当于 $K_{\text{z-cool}} < 0$（错误）
- 物理上变成了**增热**而非散热！
- 温度会**发散**而非收敛到 $T_{\text{amb}}$

### 数值影响

稳态时 $\dot{T} = 0$：
$$
KT \cdot T = -F
$$

如果 $KT$ 的符号错误：
- 导热项：$KT_{\text{cond}} < 0$（正确，散热）
- z冷却项：$KT_{\text{z-cool}} > 0$（**错误，增热**）

**总刚度矩阵**：
$$
KT_{\text{total}} = KT_{\text{cond}} + KT_{\text{z-cool}}
$$

如果 $KT_{\text{z-cool}}$ 很大且为正，可能导致：
- $KT_{\text{total}} > 0$（系统不稳定）
- 温度求解发散
- 数值不收敛

## 修正方案

### 方案1：修改表面冷却和极耳冷却

```julia
// 表面冷却（修正）
KT[ni, nj] -= wt * Ni_g * Nj_g  // 改为减号
FT[ni] += wt * T_amb_nd * Ni_g  // 载荷向量保持正号

// 极耳冷却（修正）
KT[n, n] -= coeff  // 改为减号
FT[n] += coeff * T_amb_nd  // 载荷向量保持正号
```

### 方案2：系数中包含负号

```julia
// 表面冷却（修正）
wt = -conv_factor * wJ[g]  // 系数为负
KT[ni, nj] += wt * Ni_g * Nj_g  // 装配用加号

// 极耳冷却（修正）
coeff = -h_tab * tab_area * weight / (H * k_th * L_th)  // 系数为负
KT[n, n] += coeff  // 装配用加号
```

### 推荐方案

**方案1更清晰**：
- 装配时明确减号，与外圈对流一致
- 易于理解和维护

## 验证方法

### 1. 符号检查

所有对流/冷却项的刚度矩阵贡献应该是**负号**：
```julia
@assert KT[i,i] < 0 "对流项应为负号（代码约定）"
```

### 2. 物理测试

设置简单场景：
- 初始温度：$T_0 = 350$ K
- 环境温度：$T_{\text{amb}} = 298$ K
- 无内热源：$q = 0$

**预期行为**：温度应逐渐降低到 $T_{\text{amb}}$

**如果符号错误**：温度会升高（发散）

### 3. 稳态测试

稳态时 $KT \cdot T = -F$：
```julia
KT_total = KT_cond + KT_conv + KT_surface + KT_tab
@assert all(diag(KT_total) .< 0) "总刚度对角元应为负"
```

## 总结

### 问题

1. ❌ **表面冷却**：刚度矩阵用了**正号**（应该是负号）
2. ❌ **极耳冷却**：刚度矩阵用了**正号**（应该是负号）
3. ✅ **外圈对流**：刚度矩阵用了**负号**（正确）

### 后果

- 物理意义错误：冷却变成了增热
- 数值不稳定：可能发散
- 与外圈对流不一致

### 修正

**方案1（推荐）**：装配时改为减号
```julia
KT[i, j] -= wt * Ni * Nj  // 减号
```

**方案2**：系数中包含负号
```julia
wt = -conv_factor * wJ[g]  // 负系数
KT[i, j] += wt * Ni * Nj   // 加号
```

---

**紧急修正建议**：表面冷却和极耳冷却的刚度矩阵贡献应改为**负号**！
