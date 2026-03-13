# heat_Q 参数归一化修正设计

> **创建日期**: 2026-03-13
> **状态**: 待审核
> **问题**: `heat_Q` 定义为比热容 (J/kg/K)，但归一化公式缺少密度因子，被当作体积热容处理
> **取代**: 本设计取代 `docs/plans/2026-03-13-code-implementation-plan.md` 中关于热参数归一化的描述

---

## 0. 范围与先决条件

### 0.1 范围

**包含**:
- `src/SetParams.jl` 中的 Scale struct 和归一化公式
- `src/Materialmatrix.jl` 中的体积热容计算
- `md/01_参数定义与归一化.md` 技术文档

**不包含**:
- `src/ThermalPolar2D.jl` - 使用不同的热尺度方法 (`L_th`, `k_th` 等)，不在本次修改范围内
- 电化学参数归一化 - 本次仅修改热参数

### 0.2 实现先决条件

```
[步骤 1] 添加 Scale.rho 字段
    ↓
[步骤 2] 修改归一化公式 (PE/NE/SP/PCC/NCC)
    ↓
[步骤 3] 修改 Materialmatrix.jl 体积热容计算
    ↓
[步骤 4] 更新文档
    ↓
[步骤 5] 验证测试
```

**关键**: 步骤 1 必须先完成，否则步骤 2 的代码会因 `scale.rho` 不存在而报错。

---

## 1. 问题描述

### 1.1 当前状态

**参数定义** (`SetParams.jl` 注释):
- `heat_Q -- specific heat capacities [J/(kg K)]` - 比热容
- `rho -- density [kg/m³]` - 密度

**当前归一化代码** (SetParams.jl:336):
```julia
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

**技术文档** (01_参数定义与归一化.md) 声称的公式:
$$c^* = c \cdot \rho_{ref} \cdot \frac{L^3 \cdot T_{ref}}{t_0 \cdot P_{ref}}$$

**问题**: 代码实现中**缺少 `scale.rho` 因子**，导致公式不一致。

### 1.2 代码层面的问题

1. **Scale struct 缺少 `rho` 字段**
   - `SetParams.jl:292` 使用 `param_dim.scale.rho = param_dim.cell.rho`
   - 但 Scale struct 定义中**没有** `rho` 字段
   - 这是一个运行时错误

2. **Materialmatrix.jl 错误使用 `param.layer.rho`**
   - 第 20 行: `rho_c_e[e] = fks[e, 1] * param.NE.rho + ...`
   - 将 `param.layer.rho` 当作体积热容使用
   - 但 `param.layer.rho` 实际存储的是无量纲密度 $\rho^*$

### 1.3 影响范围

- 所有使用热容矩阵的计算 (`Materialmatrix.jl`)
- Scale struct 定义 (`SetParams.jl`)
- 无量纲参数的量级
- 文档与代码的不一致导致维护困难

---

## 2. 设计方案

### 2.1 核心决策

| 决策项 | 选择 | 理由 |
|--------|------|------|
| 归一化策略 | 密度与比热容分离 | 物理意义清晰，便于参数研究 |
| 体积热容计算 | 运行时计算 $(\rho c)^* = \rho^* \cdot c^*$ | 不新增字段，改动最小 |
| 文档修改 | 全面清理，统一术语 | 删除矛盾描述，避免混淆 |

### 2.2 统一术语表

| 物理量 | 有量纲单位 | 无量纲公式 | 代码变量 |
|--------|-----------|------------|----------|
| 密度 | $\rho$ [kg/m³] | $\rho^* = \rho / \rho_{ref}$ | `param.layer.rho` |
| 比热容 | $c$ [J/(kg·K)] | $c^* = c \cdot \rho_{ref} \cdot L^3 \cdot T_{ref} / (t_0 \cdot P_{ref})$ | `param.layer.heat_Q` |
| 体积热容 | $\rho c$ [J/(m³·K)] | $(\rho c)^* = \rho^* \cdot c^*$ | **运行时计算** |

---

## 3. 代码修改

### 3.1 归一化公式修正

**文件**: `src/SetParams.jl`

**修改位置**: `NormaliseParam` 函数

**修正后公式**:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
# 其中 P_ref = phi * I_typ
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                  (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

**适用层**: PE, NE, SP, PCC, NCC

### 3.2 使用体积热容的代码

在需要体积热容 $(\rho c)^*$ 的地方，使用:
```julia
rho_c_star = param.layer.rho * param.layer.heat_Q
```

---

## 4. 文档修改

### 4.1 `01_参数定义与归一化.md` 修改要点

1. **删除矛盾的"体积热容"字段描述**
   - 体积热容不是独立字段，而是运行时计算值

2. **修正比热容归一化公式描述**
   - 确保文档公式与代码一致

3. **明确字段语义**
   - `param.layer.rho` = 无量纲密度 $\rho^*$
   - `param.layer.heat_Q` = 无量纲比热容 $c^*$
   - $(\rho c)^* = \rho^* \cdot c^*$ 在使用时计算

### 4.2 `2026-03-13-code-implementation-plan.md` 修正

删除或修正关于 `param.PE.rho` 存储"体积热容"的错误描述。

---

## 5. 验证方案

### 5.1 单元验证

```julia
# 归一化后应满足: (ρc)* = ρ* · c*
ρ_star = param.PE.rho
c_star = param.PE.heat_Q
ρc_star = ρ_star * c_star

# 应等价于直接归一化体积热容:
P_ref = scale.phi * scale.I_typ
ρc_direct = (param_dim.PE.rho * param_dim.PE.heat_Q) *
            scale.L^3 * scale.T_ref / (scale.t0 * P_ref)

@assert ρc_star ≈ ρc_direct
```

### 5.2 回归测试

运行现有热模型验证算例，确保:
- 温度场结果与修改前一致（误差 < 0.1%）
- 电压曲线一致
- 无量纲参数为 O(1) 量级

### 5.3 验收标准

| 检查项 | 标准 |
|--------|------|
| 归一化公式正确 | $(\rho c)^* = \rho^* \cdot c^*$ 单位验证通过 |
| 代码修改 | PE/NE/SP/PCC/NCC 五层全部修正 |
| 文档一致性 | 术语统一，无矛盾描述 |
| 回归测试 | 温度误差 < 0.1% |

---

## 6. 代码修改详情

### 6.1 Scale Struct 修改

**文件**: `src/SetParams.jl`

**位置**: Scale struct 定义（约第 175 行）

**添加字段**:
```julia
# 在 Scale struct 中添加
rho::Float64 = 0          # 密度尺度 = 电池平均密度 [kg/m³]
```

### 6.2 归一化公式修正

**文件**: `src/SetParams.jl`

**修改位置**: `NormaliseParam` 函数

**当前代码 (错误)**:
```julia
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

**修正后**:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
# 其中 P_ref = phi * I_typ
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                  (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

**适用层**: PE, NE, SP, PCC, NCC (共 5 处)

### 6.3 体积热容使用位置修改

**文件**: `src/Materialmatrix.jl`

**位置**: `thermal_capacity_weights_2d` 函数（第 20 行）

**当前代码 (错误)**:
```julia
rho_c_e[e] = fks[e, 1] * param.NE.rho + fks[e, 2] * param.SP.rho +
             fks[e, 3] * param.PE.rho + fks[e, 4] * param.PCC.rho +
             fks[e, 5] * param.NCC.rho
```

**修正后**:
```julia
# 体积热容 (ρc)* = ρ* · c*
rho_c_e[e] = fks[e, 1] * (param.NE.rho * param.NE.heat_Q) +
             fks[e, 2] * (param.SP.rho * param.SP.heat_Q) +
             fks[e, 3] * (param.PE.rho * param.PE.heat_Q) +
             fks[e, 4] * (param.PCC.rho * param.PCC.heat_Q) +
             fks[e, 5] * (param.NCC.rho * param.NCC.heat_Q)
```

---

## 7. 实现清单

### 7.1 代码修改（按顺序执行）

| 顺序 | 任务 | 文件 | 先决条件 |
|------|------|------|----------|
| 1 | Scale struct 添加 `rho` 字段 | `src/SetParams.jl` | 无 |
| 2 | PE 热参数归一化公式 | `src/SetParams.jl` | 任务 1 |
| 3 | NE 热参数归一化公式 | `src/SetParams.jl` | 任务 1 |
| 4 | SP 热参数归一化公式 | `src/SetParams.jl` | 任务 1 |
| 5 | PCC 热参数归一化公式 | `src/SetParams.jl` | 任务 1 |
| 6 | NCC 热参数归一化公式 | `src/SetParams.jl` | 任务 1 |
| 7 | 体积热容计算 | `src/Materialmatrix.jl` | 任务 2-6 |

### 7.2 文档修改

- [ ] 修改 `md/01_参数定义与归一化.md` - 删除矛盾描述，统一术语
- [ ] 修改 `docs/plans/2026-03-13-code-implementation-plan.md` - 修正关于 `rho` 的错误描述

### 7.3 验证

- [ ] 运行回归测试验证（温度误差 < 0.1%）

---

## 8. 与相关文档的关系

### 8.1 取代关系

本设计文档**取代** `docs/plans/2026-03-13-code-implementation-plan.md` 中关于热参数归一化的所有描述。如有矛盾，以本文档为准。

### 8.2 不在范围内

`src/ThermalPolar2D.jl` 使用独立的热尺度参数（`scale.L_th`, `scale.k_th`, `scale.rho_c_th` 等），与本次修改的统一能量尺度方案不同。如需统一，应作为后续独立任务处理。
