# 统一时间尺度热归一化 - 代码实现计划

> **注意**: 关于热参数归一化的描述已被 `docs/superpowers/specs/2026-03-13-heat-Q-normalization-fix-design.md` 取代。
> 如有矛盾，以设计文档为准。

> **创建日期**: 2026-03-13
> **基于设计文档**: `2026-03-12-unified-timescale-normalization-design.md`
> **目标**: 将代码从传统热尺度切换到统一能量尺度

---

## 1. 设计原则

### 1.1 核心原则

1. **不新增结构体字段**：直接复用现有 `Scale` 结构体字段
2. **不保留传统热尺度**：完全切换到统一能量尺度
3. **字段重定义**：现有字段的物理意义改变，但名称保持兼容

4. **保持边界条件兼容**：`h` 作为 Biot 数保留，因为边界条件代码依赖它

### 1.2 字段重定义映射

| 字段 | 原定义（传统热尺度） | 新定义（统一能量尺度） |
|------|---------------------|------------------------|
| `scale.lambda` | $k_{th} = \lambda_r$（径向等效导热率） | $k_{th} = P_{ref} / (L \cdot T_{ref})$（导热率尺度参数） |
| `scale.h` | $Bi = h \cdot L / k_{th}$（Biot 数） | **保持不变**（边界条件需要） |
| `scale.t_th` | $t_{th} = (\rho c)_{th} L^2 / k_{th}$（热扩散时间） | **删除**（统一使用 `t0`） |
| `scale.q` | $q_{th} = k_{th} T_{ref} / L^2$（参考热源） | $q_{th} = P_{ref} / L^3$（热源尺度参数） |
| `scale.L_th` | 特征热长度 | **删除**（统一使用 `L`） |

### 1.3 统一能量尺度参数

```julia
# 功率参考（不存储，仅用于计算）
P_ref = phi * I_typ    # 功率参考 [W]

# 尺度参数重定义
lambda = P_ref / (L * T_ref)         # 导热率尺度 [W/(m·K)]
h     ← Biot 数 = h_cell * L / lambda_r  # 传热系数尺度（无量纲 Biot 数）
q     = P_ref / L^3                    # 热源尺度 [W/m³]
t0    = 3600                           # 统一时间尺度 [s]
```

### 1.4 关键决策

**`h` 保持为 Biot 数**:
- 边界条件代码依赖 Biot 数
- Biot 数计算方式不变： $Bi = h \cdot L / \lambda_r$
- 这与传统方案一致，不修改边界条件代码

**体积热容归一化**:
- 使用统一能量尺度
- 公式：$(\rho c)^* = \rho c \cdot L^3 \cdot T_{ref} / (t_0 \cdot P_{ref})$
- 数学等价于：
  ```
  (\rho c)^* = \frac{(\rho c)_{th} \cdot L^3}{t_0 \cdot P_{ref}} = \rho c \cdot L^2 / (k_{th} \cdot t_{th})
    = \rho c \cdot L^2 / (k_{th}) \cdot T_{ref}
    = \rho c \cdot L^2 \cdot T_{ref} / (t_0 \cdot P_{ref})
  ```
- 注意：这里的 `rho` 实际存储的是体积热容 $(\rho c)^*$，而非无量纲密度

- `t0` 直接作为统一时间尺度，代码简洁

- **保持向后兼容**：如果将来需要回退，可以快速恢复

- **单位一致性**：新方案避免了 $(\rho c)^* \neq \rho^* \cdot c^*$ 的混淆

**代码实现**:
```julia
# ChooseCell 中（统一能量尺度参数计算）
P_ref = scale.phi * scale.I_typ
lambda = P_ref / (scale.L * scale.T_ref)
h = cell.h * scale.L / cell.lambda_r      # Biot 数
q = P_ref / scale.L^3

# 设置尺度参数
param_dim.scale.lambda = lambda
param_dim.scale.h = h                     # Biot 数
param_dim.scale.q = q
# t_th 不再需要，统一使用 scale.t0
# L_th 不再需要，统一使用 scale.L
```

---

## 2. 修改任务清单

### 任务 1: 更新 Scale 结构体 - 删除 L_th 和 t_th 字段

**文件**: `src/SetParams.jl`

**修改位置**: 约第 202-211 行

**原代码**:
```julia
    # --- Thermal scaling  ---
    L_th::Float64 = 0          # characteristic thermal length
    k_th::Float64 = 0          # reference thermal conductivity
    t_th::Float64 = 0          # thermal diffusion time scale rho_c_th L_th^2 / k_th
    q_th::Float64 = 0          # reference volumetric heat source k_th*T_ref/L_th^2
    h_th::Float64 = 0          # Biot number reference (h*L_th/k_th)
```

**新代码**:
```julia
    # --- Thermal scaling (统一能量尺度) ---
    lambda::Float64 = 0          # 导热率尺度参数 P_ref/(L*T_ref) [W/(m·K)]
    q::Float64 = 0               # 热源尺度参数 P_ref/L^3 [W/m³]
    h::Float64 = 0               # Biot 数 = h_cell*L/lambda_r (边界条件)
```

**验证**: 结构体能正常编译，默认值正确。

---

### 任务 2: 更新 ChooseCell 函数 - 热尺度计算

**文件**: `src/SetParams.jl`

**修改位置**: 约第 291-299 行

**原代码**:
```julia
    #Thermal scaling (统一使用 scale.L 作为长度尺度)
    param_dim.scale.k_th = param_dim.cell.lambda_r
    param_dim.scale.q_th = param_dim.scale.k_th * param_dim.scale.T_ref / param_dim.scale.L^2
    param_dim.scale.t_th = param_dim.cell.rho * param_dim.cell.heat_Q * param_dim.scale.L^2 / param_dim.scale.k_th
    param_dim.scale.h_th = param_dim.cell.h * param_dim.scale.L / param_dim.scale.k_th  # Biot number
```

**新代码**:
```julia
    # === 统一能量尺度热参数 ===
    # 功率参考（不存储，仅用于计算）
    P_ref = param_dim.scale.phi * param_dim.scale.I_typ

    # 热尺度参数（统一能量尺度）
    param_dim.scale.lambda = P_ref / (param_dim.scale.L * param_dim.scale.T_ref)
    param_dim.scale.h = param_dim.cell.h * param_dim.scale.L / param_dim.cell.lambda_r  # Biot 数
    param_dim.scale.q = P_ref / param_dim.scale.L^3
    # t_th 不再需要，统一使用 scale.t0
    # L_th 不再需要，统一使用 scale.L
```

**验证**:
- `P_ref = phi * I_typ` 计算正确
- `lambda` 新公式计算正确
- `h` 仍为 Biot 数（边界条件兼容）
- `q` 新公式计算正确
- `t_th` 和 `L_th` 已删除

---

### 任务 3: 更新 NormaliseParam - 电极热参数

**文件**: `src/SetParams.jl`

**修改位置**: 约第 333-334 行（PE）和 357-358 行（NE）

**原代码**:
```julia
    param.PE.lambda = param_dim.PE.lambda / param.scale.k_th
    param.PE.rho = (param_dim.PE.rho * param_dim.PE.heat_Q) / param.scale.rho_c_th
```

**新代码**:
```julia
    # 正极热参数（统一能量尺度）
    # k* = k / lambda, (ρc)* = ρc * L^2 / (lambda * t0)
    param.PE.lambda = param_dim.PE.lambda / param.scale.lambda
    param.PE.rho = (param_dim.PE.rho * param_dim.PE.heat_Q) * param.scale.L^2 / (param.scale.lambda * param.scale.t0)
```

**同样修改 NE**:
```julia
    # 负极热参数（统一能量尺度）
    param.NE.lambda = param_dim.NE.lambda / param.scale.lambda
    param.NE.rho = (param_dim.NE.rho * param_dim.NE.heat_Q) * param.scale.L^2 / (param.scale.lambda * param.scale.t0)
```

**注意**:
- 直接计算：$L^2 / (\lambda \cdot t_0)$
- 这等价于 $\rho c / (\rho c)_{th}$
- 注：`param.PE.rho` 实际存储的是无量纲体积热容 $(\rho c)^*$，而非无量纲密度

**单元测试验证**:
```julia
# 验证尺度参数
P_ref = scale.phi * scale.I_typ
@assert scale.lambda ≈ P_ref / (scale.L * scale.T_ref)
@assert scale.h ≈ cell.h * scale.L / cell.lambda_r
@assert scale.q ≈ P_ref / scale.L^3
```

**验证**:
- `param.PE.lambda` 和 `param.PE.rho` 的无量纲值正确
- `param.PE.rho` 经计算后应为体积热容 $(\rho c)^*$
- 单元测试通过
- 温度场结果与修改前一致（误差 < 1%)
- 无量纲参数为 O(1) 量级
- 电压曲线一致
- CZM 损伤演化一致
