# 归一化问题发现记录

**创建日期**: 2026-03-23

---

## 1. 问题发现

### 1.1 现象描述

在二维分布式热模型中，边界温度始终接近环境温度，导致：
1. 温度梯度过大（T_vol - T_edge ≈ 20 K）
2. 散热效率过低（~15-20% vs 集总模型的 80-90%）
3. 无法达到稳态（净热源持续为正）

### 1.2 根因定位

**代码位置**: `src/ThermalDistributed.jl` 第53行

```julia:src/ThermalDistributed.jl
function apply_convection_bc(KT, FT, mesh, is_outer, case)
    K = copy(KT)
    F = copy(FT)

    Bi = case.param.cell.h  # ❌ 错误！这不是 Biot 数
```

**问题**: `case.param.cell.h` 是**传热系数的无量纲形式**（用于集总模型），而不是**Biot 数**（用于边界条件）。

---

## 2. 参数归一化分析

### 2.1 两个 "h" 的定义

| 参数 | 定义位置 | 公式 | 用途 |
|------|----------|------|------|
| `scale.h` | SetParams.jl:296 | `h*L/lambda_r` | 分布式边界条件的 Biot 数 |
| `param.cell.h` | SetParams.jl:394 | `h*A*T_ref/P_ref` | 集总模型的传热系数 |

### 2.2 代码对比

**`scale.h` 的定义**（SetParams.jl:296）：
```julia:src/SetParams.jl
param_dim.scale.h = param_dim.cell.h * param_dim.scale.L / param_dim.cell.lambda_r  # Biot 数
```

**`param.cell.h` 的定义**（SetParams.jl:394）：
```julia:src/SetParams.jl
param.cell.h = param_dim.cell.h * param_dim.cell.area * param.scale.T_ref / param.scale.phi / param.scale.I_typ
```

### 2.3 量纲分析

**`scale.h`（Biot 数）**：
$$Bi = \frac{h \cdot L}{\lambda_r}$$

单位验证：
$$[Bi] = \frac{[W/(m^2 \cdot K)] \cdot [m]}{[W/(m \cdot K)]} = [-]$$

**`param.cell.h`（传热系数无量纲）**：
$$h^* = \frac{h \cdot A \cdot T_{ref}}{P_{ref}}$$

单位验证：
$$[h^*] = \frac{[W/(m^2 \cdot K)] \cdot [m^2] \cdot [K]}{[W]} = [-]$$

---

## 3. 边界条件分析

### 3.1 正确的无量纲边界条件

**有量纲形式**：
$$-k \frac{\partial T}{\partial n} = h(T - T_{amb})$$

**无量纲形式**：
$$-k^* \frac{\partial T^*}{\partial n^*} = Bi \cdot (T^* - T_{amb}^*)$$

其中：
- $k^* = k \cdot L \cdot T_{ref} / P_{ref}$（无量纲导热率）
- $Bi = h \cdot L / k$（Biot 数）
- $T^* = T / T_{ref}$（无量纲温度）

### 3.2 边界积分公式

**有限元弱形式**：
$$\int_{\Gamma} h \cdot N_i \cdot N_j \, ds = \int_{\Gamma} Bi \cdot N_i \cdot N_j \, ds^*$$

其中 $ds^* = ds/L$ 是无量纲弧长。

**代码实现**：
```julia
wt = Bi * w * J  # w 是高斯权重，J 是无量纲 Jacobian
```

### 3.3 错误原因

如果使用 `param.cell.h`（量级 ~10^4）而不是 `scale.h`（量级 ~0.001），会导致：

1. **边界积分过大**：散热项被放大
2. **边界温度过低**：强制边界温度接近 T_amb
3. **物理意义错误**：混淆了集总参数和分布式参数

---

## 4. 数值验证

### 4.1 Jellyroll 参数（来自代码）

```julia
# 从 parameters/Jellyroll.jl
h = 10.0                  # 对流换热系数 [W/(m²·K)]
lambda_r = 1.30           # 径向等效导热率 [W/(m·K)]
L = 1.56e-4               # 电极堆叠厚度 [m]
A = 0.0048                # 冷却面积 [m²]
T_ref = 298.0             # 参考温度 [K]
phi = 0.0257              # 电势尺度 [V]
I_typ = 0.0147            # 电流尺度 [A]
P_ref = phi * I_typ       # 功率参考 [W]
```

### 4.2 计算对比

| 参数 | 公式 | 数值 |
|------|------|------|
| `scale.h` | h*L/lambda_r | **0.0012** |
| `param.cell.h` | h*A*T_ref/P_ref | **3.76×10^4** |
| 比值 | param.cell.h / scale.h | **3.13×10^7** |

### 4.3 物理意义

- `scale.h = 0.0012`：对流热阻 >> 导热热阻，边界温度应接近内部温度
- `param.cell.h = 3.76×10^4`：被错误地用于边界条件，导致散热过大

---

## 5. 集总模型 vs 分布式模型

### 5.1 集总模型

**散热公式**：
$$P_{out} = h \cdot A \cdot (T - T_{amb})$$

**无量纲形式**：
$$P_{out}^* = h^* \cdot (T^* - T_{amb}^*)$$

其中 $h^* = h \cdot A \cdot T_{ref} / P_{ref}$，即 `param.cell.h`。

**使用位置**：`src/Thermal.jl`

### 5.2 分布式模型

**边界条件**：
$$-k \frac{\partial T}{\partial n} = h(T - T_{amb})$$

**无量纲形式**：
$$-k^* \frac{\partial T^*}{\partial n^*} = Bi \cdot (T^* - T_{amb}^*)$$

其中 $Bi = h \cdot L / k$，即 `scale.h`。

**使用位置**：`src/ThermalDistributed.jl`

### 5.3 关键区别

| 特征 | 集总模型 | 分布式模型 |
|------|----------|------------|
| 温度 | 单一值 T | 分布场 T(x,y) |
| 散热 | 基于总体积平均温度 | 基于边界节点温度 |
| 参数 | h*A*T_ref/P_ref | h*L/k |
| 物理意义 | 传热功率/电功率 | 对流热阻/导热热阻 |

---

## 6. 从第一性原理推导边界条件传热系数（2026-03-23 更新）

### 6.1 有量纲边界条件

侧面对流边界条件（Robin 边界）：

$$-k \frac{\partial T}{\partial n} = h(T - T_{amb})$$

### 6.2 有限元弱形式中的边界积分

边界积分项：

$$\int_{\Gamma} h(T - T_{amb}) N_i \, ds = \int_{\Gamma} h \, T \, N_i \, ds - \int_{\Gamma} h \, T_{amb} \, N_i \, ds$$

对应：
- **刚度矩阵贡献**：$K_{ij} += \int_{\Gamma} h \, N_i \, N_j \, ds$
- **载荷向量贡献**：$F_i += \int_{\Gamma} h \, T_{amb} \, N_i \, ds$

### 6.3 无量纲化推导

**关键变量无量纲化**：
- $ds = L \cdot ds^*$（边界弧长）
- $T = T_{ref} \cdot T^*$（温度）

**边界积分无量纲化**：

$$\int_{\Gamma} h(T - T_{amb}) N_i \, ds = \int_{\Gamma^*} h \cdot T_{ref}(T^* - T_{amb}^*) \cdot N_i \cdot L \, ds^*$$

$$= L \cdot T_{ref} \cdot h \int_{\Gamma^*} (T^* - T_{amb}^*) N_i \, ds^*$$

**热流密度无量纲化**（基于能量守恒，统一使用 $P_{ref}$）：

$$q_s^* = q_s \cdot \frac{L^2}{P_{ref}}$$

因此边界热流：

$$q_s^* = h(T - T_{amb}) \cdot \frac{L^2}{P_{ref}} = \frac{h \cdot T_{ref} \cdot L^2}{P_{ref}}(T^* - T_{amb}^*)$$

### 6.4 边界条件传热系数定义

**定义边界条件无量纲传热系数**：

$$\boxed{h_{conv}^* = \frac{h \cdot T_{ref} \cdot L^2}{P_{ref}}}$$

**与现有参数的关系**：

$$h_{conv}^* = \frac{h \cdot T_{ref} \cdot L^2}{P_{ref}} = \frac{h \cdot L}{\lambda_r} \cdot \frac{\lambda_r \cdot L \cdot T_{ref}}{P_{ref}} = Bi \times k_r^*$$

即：

$$h_{conv}^* = \text{scale.h} \times \text{param.cell.lambda\_r}$$

### 6.5 数值验证

**Jellyroll 参数**：
- h = 10 W/(m²·K)
- L = 1.56×10⁻⁴ m
- T_ref = 298 K
- P_ref = 3.78×10⁻⁴ W
- λ_r = 1.30 W/(m·K)

**直接计算**：

$$h_{conv}^* = \frac{10 \times 298 \times (1.56 \times 10^{-4})^2}{3.78 \times 10^{-4}} = \frac{10 \times 298 \times 2.43 \times 10^{-8}}{3.78 \times 10^{-4}} \approx \mathbf{0.19}$$

**通过 Biot 数计算**：

$$Bi = \frac{10 \times 1.56 \times 10^{-4}}{1.30} = 0.0012$$

$$k_r^* = \frac{1.30 \times 1.56 \times 10^{-4} \times 298}{3.78 \times 10^{-4}} = 160$$

$$h_{conv}^* = 0.0012 \times 160 = 0.19 \quad \checkmark$$

### 6.6 与错误用法的对比

| 参数 | 公式 | 数值 | 用途 |
|------|------|------|------|
| `param.cell.h` (错误) | h·A·T_ref/P_ref | 3.78×10⁴ | 集总模型 |
| `scale.h` | h·L/λ_r | 0.0012 | Biot 数 |
| **$h_{conv}^*$ (正确)** | h·T_ref·L²/P_ref | **0.19** | **边界条件** |

**差异来源**：$A / L^2 = 0.0048 / (1.56 \times 10^{-4})^2 = 2 \times 10^5$

---

## 7. 修正方案

### 7.1 代码修改

**位置**：`src/ThermalDistributed.jl` 第53行

```julia
# 修正前（错误）
Bi = case.param.cell.h  # ❌ 这是集总模型参数，不是边界条件参数

# 修正后（正确）
h_conv_nd = case.param_dim.scale.h * case.param.cell.lambda_r
# 或者等价地：
# h_conv_nd = case.param_dim.cell.h * param_dim.scale.T_ref * param_dim.scale.L^2 / (param_dim.scale.phi * param_dim.scale.I_typ)
```

### 7.2 文档中需要删除的错误引用

以下文档中将 `param.cell.h` 称为 "Biot 数" 是错误的：

1. `md/05_热模型_二维分布式.md`：第53行注释 `Bi = case.param.cell.h  # Biot 数（统一能量尺度）`
2. 相关分析文档中的类似错误引用

**注意**：`scale.h` 才是真正的 Biot 数，`param.cell.h` 是集总模型传热系数。

### 7.3 后续工作

1. ✅ 从第一性原理推导正确的边界条件参数
2. [ ] 修改 `ThermalDistributed.jl` 代码
3. [ ] 删除文档中的错误引用
4. [ ] 运行验证脚本
5. [ ] 更新技术文档

---

## 8. 参考资料

- `md/01_参数定义与归一化.md`：归一化方案文档
- `md/05_热模型_二维分布式.md`：热模型文档
- `docs/thermal_verify/2D模型侧面散热低估问题分析.md`：问题分析文档
