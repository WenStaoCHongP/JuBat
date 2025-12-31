# 表面冷却归一化问题分析

## 问题发现

对比外边界对流和表面冷却的h处理，发现**归一化不一致**。

## SetParams.jl 中的定义

### Scale结构（第180行）

```julia
h_th::Float64 = 0  # Biot number reference (h*L_th/k_th)
```

**说明**：`h_th` 是**无量纲Biot数**，不是有量纲的对流系数！

### 归一化处理（第264行）

```julia
h_ref = param_dim.cell.h * L_th / k_ref  # Biot number (dimensionless h)
param_dim.scale.h_th = h_ref
```

**公式**：
$$
h_{\text{th}} = \frac{h \cdot L_{\text{th}}}{k_{\text{th}}} = Bi
$$

**含义**：
- 输入：`param_dim.cell.h` [W/(m²·K)]（有量纲）
- 输出：`param_dim.scale.h_th` [无量纲]（Biot数）

## 外边界对流的使用

### 代码（第356行）

```julia
function _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    scale = case.param_dim.scale
    Bi = scale.h_th  // ← 使用无量纲Biot数
    
    Bi == 0 && return
    
    L_th = scale.L_th
    ...
    wt = Bi * w * (J / L_th)  // ← Bi已经是无量纲
    ...
end
```

### 无量纲分析

对流项系数：
$$
\text{coeff} = \frac{Bi}{L_{\text{th}}} = \frac{h \cdot L_{\text{th}} / k_{\text{th}}}{L_{\text{th}}} = \frac{h}{k_{\text{th}}}
$$

**正确** ✅

## 表面冷却的使用（问题所在）

### 代码（第499-512行）

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    # 获取参数
    h_surface = hasproperty(case.opt, :h_surface) ? case.opt.h_surface : 10.0  # ← 有量纲！W/(m²·K)
    H = hasproperty(case.param_dim.cell, :height) ? case.param_dim.cell.height : case.param_dim.cell.width
    
    scale = case.param_dim.scale
    k_th = scale.k_th
    L_th = scale.L_th
    T_ref = scale.T_ref
    T_amb_nd = case.param_dim.cell.T_amb / T_ref
    
    # 体积散热系数：2h/H
    vol_coeff = 2.0 * h_surface / H  // ← h_surface是有量纲的
    
    # 无量纲Biot数：Bi_z = 2h*L_th^2 / (H*k_th)
    Bi_z = vol_coeff * L_th^2 / k_th  // ← 自己计算Biot数
    
    conv_factor = Bi_z / L_th^2  // ← 再除以L_th^2
    
    for g in 1:ngs
        wt = conv_factor * wJ[g]  // ← 使用conv_factor
        ...
    end
end
```

### 问题分析

**问题1：h的来源不一致**
- 外边界对流：使用 `scale.h_th`（**无量纲**Biot数）
- 表面冷却：使用 `opt.h_surface`（**有量纲** W/(m²·K)）

**问题2：重复归一化？**

让我们追踪 `conv_factor` 的量纲：

1. `h_surface` [W/(m²·K)]（有量纲）
2. `vol_coeff = 2*h_surface/H` [W/(m³·K)]
3. `Bi_z = vol_coeff * L_th^2 / k_th` 
   $$
   = \frac{2h_{\text{surface}} \cdot L_{\text{th}}^2}{H \cdot k_{\text{th}}}
   $$
   [无量纲]

4. `conv_factor = Bi_z / L_th^2`
   $$
   = \frac{2h_{\text{surface}}}{H \cdot k_{\text{th}}}
   $$
   [1/m²]

5. `wt = conv_factor * wJ[g]`
   - `wJ[g]` 的量纲是什么？

### wJ的量纲

`wJ[g] = weight[g] * detJ[g]`
- `weight[g]`：高斯权重（无量纲）
- `detJ[g]`：雅可比行列式 [m²]（从参数空间到物理空间）

所以 `wJ[g]` [m²]

因此：
$$
wt = conv\_factor \cdot wJ = \frac{2h_{\text{surface}}}{H \cdot k_{\text{th}}} \cdot [m^2] = \frac{2h_{\text{surface}} \cdot [m^2]}{H \cdot k_{\text{th}}}
$$

量纲：$[W/(m^2 \cdot K)] \cdot [m^2] / ([m] \cdot [W/(m \cdot K)]) = [无量纲]$ ✓

看起来量纲是对的，但是...

### 与外边界对流的对比

**外边界对流**：
```julia
Bi = scale.h_th  // 无量纲
wt = Bi * w * (J / L_th)
```

$$
wt = \frac{h \cdot L_{\text{th}}}{k_{\text{th}}} \cdot w \cdot \frac{J}{L_{\text{th}}} = \frac{h \cdot J \cdot w}{k_{\text{th}}}
$$

其中 $J$ 是边长的一半 [m]，$w$ 是高斯权重（无量纲）。

**表面冷却**：
```julia
h_surface = 10.0  // 有量纲 [W/(m²·K)]
conv_factor = 2*h_surface / (H*k_th)
wt = conv_factor * wJ
```

$$
wt = \frac{2h_{\text{surface}}}{H \cdot k_{\text{th}}} \cdot wJ
$$

### 问题所在

**关键问题**：`h_surface` 是**有量纲**的，但 `scale.h_th` 是**无量纲**的！

如果用户在 Jellyroll.jl 中设置了：
```julia
cell.h = 10.0  # W/(m²·K)
```

那么经过 NormaliseParam 后：
```julia
scale.h_th = cell.h * L_th / k_th  # 无量纲Biot数
```

**但表面冷却直接从 `opt.h_surface` 读取有量纲的值**，而不是使用已经归一化的 `scale.h_th`！

## 正确的做法

### 方案1：使用无量纲Biot数（推荐）

表面冷却应该像外边界对流一样，使用已经归一化的Biot数：

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    scale = case.param_dim.scale
    
    # 方案1a：使用现有的scale.h_th（外圈对流的Biot数）
    Bi_outer = scale.h_th
    
    # 或者方案1b：定义专门的z方向Biot数
    h_surface_dim = hasproperty(case.param_dim.cell, :h_surface) ? 
                    case.param_dim.cell.h_surface : 
                    case.param_dim.cell.h
    Bi_surface = h_surface_dim * scale.L_th / scale.k_th
    
    H = case.param_dim.cell.width
    L_th = scale.L_th
    
    # 体积散热系数（无量纲）
    # q_vol* = (2h/H) * (T - T_amb) / (k_th/L_th)
    #        = (2h*L_th)/(H*k_th) * (T - T_amb)
    #        = (2*Bi*L_th/H) * (T - T_amb)
    vol_coeff_nd = 2.0 * Bi_surface / H  # 注意：这里Bi已经包含了L_th
    
    # 不需要再除以L_th^2，因为Bi已经包含了L_th
    conv_factor = vol_coeff_nd / L_th  # [1/m]
    
    for g in 1:ngs
        wt = conv_factor * wJ[g]  # wJ[g]的量纲是[m²]
        ...
    end
end
```

等等，我需要重新理解无量纲化...

### 重新分析无量纲化

**物理方程**：
$$
\rho c \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + q - \frac{2h}{H}(T - T_{\text{amb}})
$$

**无量纲化**：
- $T^* = T / T_{\text{ref}}$
- $x^* = x / L_{\text{th}}$
- $t^* = t \cdot (k_{\text{th}} / (\rho c_{\text{ref}} L_{\text{th}}^2))$
- $q^* = q \cdot L_{\text{th}}^2 / (k_{\text{th}} T_{\text{ref}})$

**无量纲方程**：
$$
\frac{\partial T^*}{\partial t^*} = \nabla^{*2} T^* + q^* - \frac{2h L_{\text{th}}^2}{H k_{\text{th}}} (T^* - T_{\text{amb}}^*)
$$

定义：
$$
Bi_z = \frac{2h L_{\text{th}}^2}{H k_{\text{th}}}
$$

**弱形式（无量纲）**：
$$
\int M \frac{\partial T^*}{\partial t^*} N_i \, d\Omega^* + \int \nabla^* T^* \cdot \nabla^* N_i \, d\Omega^* + \int Bi_z T^* N_i \, d\Omega^* = \int q^* N_i \, d\Omega^* + \int Bi_z T_{\text{amb}}^* N_i \, d\Omega^*
$$

其中 $d\Omega^* = d\Omega / L_{\text{th}}^2$。

**离散化**：
$$
K_{ij} = -\int \nabla^* N_i \cdot \nabla^* N_j \, d\Omega^* - \int Bi_z N_i N_j \, d\Omega^*
$$

$$
= -\int \frac{1}{L_{\text{th}}^2} \nabla N_i \cdot \nabla N_j \, d\Omega - \int Bi_z \frac{1}{L_{\text{th}}^2} N_i N_j \, d\Omega
$$

$$
= -\frac{1}{L_{\text{th}}^2} \int (\nabla N_i \cdot \nabla N_j + Bi_z N_i N_j) \, d\Omega
$$

### 正确的实现

```julia
# 获取有量纲的h
h_surface = case.param_dim.cell.h_surface  # [W/(m²·K)]

# 计算无量纲Biot数
Bi_z = 2.0 * h_surface * L_th^2 / (H * k_th)

# 系数（包含在积分中）
conv_factor = Bi_z / L_th^2  # [1/m²]

for g in 1:ngs
    wt = conv_factor * wJ[g]  # wJ[g] = weight * detJ，detJ的量纲是[m²]
    # wt是无量纲的
    KT[i,j] -= wt * Ni * Nj
end
```

这样看起来是对的！

但问题在于：**h_surface从哪里来？**

## 真正的问题

**外边界对流**：
- 从 `param_dim.cell.h` 获取有量纲值
- 在 `NormaliseParam` 中归一化为 `scale.h_th`
- 使用时直接读取 `scale.h_th`

**表面冷却**：
- 从 `opt.h_surface` 获取有量纲值（**没有归一化**！）
- 使用时自己计算Biot数

### 不一致之处

1. **参数来源不同**：
   - 外边界：`param_dim.cell.h` → `scale.h_th`
   - 表面冷却：`opt.h_surface`（直接从运行时选项读取）

2. **归一化时机不同**：
   - 外边界：在 `NormaliseParam` 中预先归一化
   - 表面冷却：在使用时临时归一化

3. **量纲混乱**：
   - 如果用户误以为 `opt.h_surface` 是无量纲的，结果就错了
   - 如果用户想用与外边界相同的h，需要设置两次

## 解决方案

### 方案1：统一使用无量纲Biot数（推荐）

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    scale = case.param_dim.scale
    
    # 使用统一的Biot数定义
    # 如果opt中有h_surface_Bi，使用它；否则使用scale.h_th
    if hasproperty(case.opt, :h_surface_Bi)
        Bi_surface = case.opt.h_surface_Bi  # 无量纲Biot数
    else
        # 回退到scale.h_th（与外边界相同）
        Bi_surface = scale.h_th
    end
    
    H = case.param_dim.cell.width
    L_th = scale.L_th
    
    # z方向Biot数
    Bi_z = 2.0 * Bi_surface * L_th / H
    
    conv_factor = Bi_z / L_th^2
    
    for g in 1:ngs
        wt = conv_factor * wJ[g]
        KT[i,j] -= wt * Ni * Nj
    end
end
```

### 方案2：从param_dim读取有量纲值，统一归一化

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    scale = case.param_dim.scale
    
    # 从param_dim读取有量纲的h（而非opt）
    h_surface = hasproperty(case.param_dim.cell, :h_surface) ?
                case.param_dim.cell.h_surface :
                case.param_dim.cell.h  # 回退到外边界的h
    
    H = case.param_dim.cell.width
    k_th = scale.k_th
    L_th = scale.L_th
    
    # 计算z方向Biot数
    Bi_z = 2.0 * h_surface * L_th^2 / (H * k_th)
    
    conv_factor = Bi_z / L_th^2
    
    for g in 1:ngs
        wt = conv_factor * wJ[g]
        KT[i,j] -= wt * Ni * Nj
    end
end
```

## 问题总结

### 当前实现的问题

1. ❌ **参数来源不一致**：
   - 外边界从 `param_dim.cell.h` 读取
   - 表面冷却从 `opt.h_surface` 读取

2. ❌ **归一化不统一**：
   - 外边界使用预归一化的 `scale.h_th`
   - 表面冷却自己临时归一化

3. ❌ **默认值问题**：
   - `opt.h_surface` 默认10.0，但没说清楚量纲
   - 容易混淆是有量纲还是无量纲

4. ❌ **用户体验差**：
   - 如果想让表面冷却和外边界用相同的h，需要设置两次

### 推荐修正

**方案2更好**：从 `param_dim.cell` 读取有量纲值，保持与外边界一致。

```julia
// 修正后
h_surface = hasproperty(case.param_dim.cell, :h_surface) ?
            case.param_dim.cell.h_surface :
            case.param_dim.cell.h  // 回退到外边界的h

Bi_z = 2.0 * h_surface * L_th^2 / (H * k_th)
```

这样：
- ✅ 参数来源统一（都从 `param_dim.cell`）
- ✅ 归一化方式一致（都使用 h * L / k 的形式）
- ✅ 默认行为合理（回退到外边界的h）
- ✅ 用户友好（只需设置一次h）
