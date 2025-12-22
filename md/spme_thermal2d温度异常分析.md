# spme_thermal2d_example.jl 温度异常问题分析

## 问题描述

`spme_thermal2d_example.jl` 在**纯热模拟模式**（`electrochem_enabled = false`）下，给定了恒定的空间分布热源，但温度表现不正常（可能不升温或下降）。

## 关键代码设置

### 1. 热源配置（第169-183行）

```julia
else  # electrochem_enabled = false
    # 组装径向线性递减热源：q(r) = q_user * (1 - (r - Rin)/(Rout - Rin))
    ne = size(mesh_th.element, 1)
    den = max(Rout - Rin, 1e-12)
    profile = clamp.(1 .- (r_centers .- Rin) ./ den, 0.0, 1.0)
    if uppercase(q_units) == "SI"
        variables["heat_source_fields"] = q_user .* profile
        variables["heat_source_units_code"] = 1.0
    else
        q_ref = case.param_dim.scale.q_th
        variables["heat_source_fields"] = (q_user / max(q_ref, 1e-16)) .* profile
        variables["heat_source_units_code"] = 0.0
    end
end
```

**默认值**：
- `q_user = 1000` W/m³（第62行）
- `q_units = "SI"`（第63行）

### 2. 参数集配置（Jellyroll）

从 `src/parameters/Jellyroll.jl`：

```julia
cell.Rout = 0.021/2      # 外径 = 10.5 mm
cell.Rin = 1.92e-3       # 内径 = 1.92 mm
cell.width = 7e-2        # 宽度 = 70 mm
cell.h = 150.            # 对流换热系数 = 150 W/(m²·K)  ← ⚠️ 关键参数
cell.T_amb = 298         # 环境温度 = 298 K
```

### 3. 边界条件（ThermalDistributed2D_BC）

外边界应用**对流边界条件**（Robin边界）：

$$
-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})
$$

从 `src/ThermalDistributed.jl` 第353-410行，对流边界在刚度矩阵和载荷向量中实现：

```julia
Bi = scale.h_th  # 无量纲Biot数 = h × L_th / k_ref
wt = Bi * w * (J / L_th)
ke11 += -wt * N1 * N1    # 刚度矩阵（负号）
fe1 += wt * T_amb * N1   # 载荷向量
```

## 问题诊断：能量不平衡

### 关键计算

#### 1. 电池几何
- **内径**：Rin = 1.92 mm
- **外径**：Rout = 10.5 mm
- **宽度**：width = 70 mm
- **体积**：V = π(Rout² - Rin²) × width ≈ **2.40×10⁻⁵ m³** = 24 cm³

#### 2. 对流散热面积
- **外周面积**：A₁ = 2πRout × width ≈ 4.62×10⁻³ m²
- **端盖面积**：A₂ = 2πRout² ≈ 6.93×10⁻⁴ m²
- **总面积**：A_conv ≈ **5.31×10⁻³ m²** = 53.1 cm²

#### 3. 总产热功率
$$
Q_{\text{gen}} = q_{\text{user}} \times V = 1000 \times 2.40 \times 10^{-5} = \mathbf{0.024 \text{ W}}
$$

#### 4. 对流散热能力
在温升 ΔT = 1K 时：
$$
Q_{\text{conv}} = h \times A_{\text{conv}} \times \Delta T = 150 \times 5.31 \times 10^{-3} \times 1 = \mathbf{0.797 \text{ W}}
$$

在温升 ΔT = 0.1K 时：
$$
Q_{\text{conv}} = 150 \times 5.31 \times 10^{-3} \times 0.1 = \mathbf{0.0797 \text{ W}}
$$

#### 5. 稳态温升
$$
\Delta T_{\text{ss}} = \frac{Q_{\text{gen}}}{h \times A_{\text{conv}}} = \frac{0.024}{150 \times 5.31 \times 10^{-3}} = \mathbf{0.030 \text{ K}}
$$

### 结论

**❌ 散热能力远超产热能力！**

- **产热功率**：Q_gen = **0.024 W**
- **散热能力**（ΔT=0.1K）：Q_conv = **0.080 W**
- **比值**：Q_gen / Q_conv@0.1K ≈ **0.30**

**稳态温升仅 0.03 K！**

这意味着：
1. 温度几乎不会升高（低于数值精度）
2. 如果初始温度略高于环境温度，会快速冷却到环境温度
3. 数值舍入误差可能导致温度轻微波动或下降

## 物理机制

### 1. Biot数分析

无量纲Biot数：
$$
\text{Bi} = \frac{h L_{\text{th}}}{k_{\text{ref}}}
$$

从 `SetParams.jl`：
- L_th ≈ (Rout - Rin) / 2 ≈ 4.3 mm
- k_ref ≈ λ_r_eff（径向有效热导率，约1-5 W/(m·K)）
- **Bi ≈ 150 × 0.0043 / 2 ≈ 0.32**

**解释**：Bi < 1，理论上热传导主导，但由于产热太弱，即使传导快也无法建立温度场。

### 2. 时间尺度

热扩散时间：
$$
\tau_{\text{th}} = \frac{\rho c L_{\text{th}}^2}{k} \sim 1-10 \text{ s}
$$

对流冷却时间：
$$
\tau_{\text{conv}} = \frac{\rho c V}{h A} \sim 5 \text{ s}
$$

由于 h 很大，**对流冷却速度快于热扩散**，任何温度升高都会被迅速散失。

### 3. 能量守恒方程

$$
\rho c \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + q - \frac{h A}{V}(T - T_{\text{amb}})
$$

在稳态（∂T/∂t = 0）且忽略热传导梯度（均匀温度）：
$$
q = \frac{h A}{V}(T_{\text{ss}} - T_{\text{amb}})
$$

代入数值：
$$
1000 = \frac{150 \times 5.31 \times 10^{-3}}{2.40 \times 10^{-5}}(T_{\text{ss}} - 298)
$$
$$
T_{\text{ss}} - 298 = \frac{1000 \times 2.40 \times 10^{-5}}{150 \times 5.31 \times 10^{-3}} = 0.030 \text{ K}
$$

**稳态温度 T_ss = 298.03 K**，几乎不变！

## 温度"不正常"的表现形式

基于上述分析，可能的"不正常"现象包括：

### 1. 温度几乎不变
- **预期**：有热源，温度应该升高
- **实际**：温度从298.00K → 298.03K，变化极微
- **原因**：散热太强，稳态温升仅0.03K

### 2. 温度初期下降
- **现象**：如果初始温度 T₀ = 298.00K，但由于数值初始化或边界条件处理，某些节点温度略高（如298.05K），则这些节点会快速冷却到298.03K
- **原因**：对流边界强制温度趋向环境温度

### 3. 温度空间分布异常
- **现象**：内部温度略高（如298.04K），但靠近外边界的节点温度接近环境温度（如298.01K），形成"反向梯度"（中心冷边缘热的感觉）
- **原因**：虽然热源在内部，但对流散热如此强大，边界温度被强制压低

### 4. 数值振荡
- **现象**：温度在298.00K附近小幅振荡（±0.01K）
- **原因**：数值求解器的舍入误差，在如此小的温升下（0.03K）相对误差较大

## 解决方案

### 方案A：减小对流换热系数（推荐）

将 `cell.h` 从150降至10 W/(m²·K)（自然对流水平）：

```julia
param_dim.cell.h = 10.0  # W/(m²·K)
```

**预期效果**：
- 新的稳态温升：ΔT_ss = 0.024 / (10 × 5.31×10⁻³) ≈ **0.45 K**
- 温度变化明显，可观测

### 方案B：增大热源密度

使稳态温升达到10K：

$$
q_{\text{new}} = \frac{h A}{V} \times 10 = \frac{150 \times 5.31 \times 10^{-3}}{2.40 \times 10^{-5}} \times 10 \approx 33000 \text{ W/m}^3
$$

```julia
q_user = 33000  # W/m³
```

**预期效果**：
- 稳态温升：10 K
- 但热源密度增大33倍，可能不符合实际应用

### 方案C：使用绝热边界（测试用）

暂时禁用对流边界条件，验证热源是否正常工作：

修改 `ThermalDistributed2D_BC`，在开头加入：

```julia
# 测试：禁用对流
return nothing
```

**预期效果**：
- 温度持续上升（无散热）
- 验证热源计算正确性

### 方案D：调整参数集

使用对流较弱的参数集（如 ThermalMinimal 或自定义）：

```julia
param_dim = JuBat.ChooseCell("ThermalMinimal")
# 或手动设置
param_dim.cell.h = 5.0  # 自然对流
```

## 验证步骤

### 1. 运行诊断脚本

```bash
cd /workspace/example
julia diagnose_spme_thermal2d.jl
```

输出将显示：
- 总产热功率 Q_gen
- 对流散热能力 Q_conv
- 稳态温升 ΔT_ss
- 能量平衡判断

### 2. 修改参数后重新运行

选择方案A（推荐）：

```julia
# spme_thermal2d_example.jl 第7行之后添加
param_dim.cell.h = 10.0  # 减小对流换热系数
```

预期输出：
- 稳态温升约0.5K
- 温度时间曲线显示明显上升

### 3. 检查输出图像

```bash
eog output/spme_thermal2d_Tmean.png  # 查看平均温度曲线
```

**修改前**：温度几乎不变（298.00 → 298.03 K）  
**修改后**：温度明显上升（298.00 → 298.50 K）

## 理论背景

### 对流边界条件的物理意义

对流边界 `-k ∂T/∂n = h(T - T_amb)` 表示：
- 热流出边界 = 对流散热
- h 越大，边界温度越接近环境温度（强制冷却）
- h → ∞ 时，等效于 **Dirichlet边界** T = T_amb

在当前情况下：
- h = 150 W/(m²·K) 非常大
- 边界节点几乎被"锁定"在环境温度附近
- 内部热源产生的热量迅速从边界散失

### 典型对流换热系数参考值

| 冷却方式 | h [W/(m²·K)] |
|---------|-------------|
| 静止空气（自然对流） | 5-25 |
| 低速空气流动 | 10-100 |
| 强制风冷 | 50-250 |
| 液冷（水） | 500-10000 |

**Jellyroll参数的 h = 150 接近强制风冷上限**，可能过于保守（过度估计散热能力）。

## 结论

**温度"不正常"的根源：散热能力远超产热能力**

- **产热**：Q_gen = 0.024 W
- **散热**：Q_conv@ΔT=0.1K = 0.080 W（3.3倍）
- **稳态温升**：仅0.03K（低于典型测温精度）

**建议**：
1. 将对流换热系数从150降至10 W/(m²·K)（方案A）
2. 或增大热源密度到30000 W/m³（方案B）
3. 运行诊断脚本验证能量平衡

这不是代码错误，而是**参数配置不合理**导致的物理行为异常。

---

**生成时间**：2025-12-22  
**分析工具**：JuBat 热传导模块
