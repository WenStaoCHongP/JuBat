# CZM 渐进式有效面积损失机制设计

**日期**: 2026-05-26
**状态**: 已批准
**前置文档**: `2026-05-26-czm-damage-convergence-analysis-design.md`

---

## 1. 背景与动机

在单次放电仿真中，D_max 收敛于 0.96，无法达到 1.0（完全断裂）。根本原因：

1. δ_c/δ₀ = 1806，损伤公式 D ≈ 1 - δ₀/δ 在 δ=8.5nm 时就达到 D=0.96，但断裂需要 δ=617nm
2. 0 < D < 0.99 的单元在电化学上完全等价于完好单元——存在"死区"
3. 不存在足够强的正反馈推动 δ 从 8.5nm 增长到 617nm

**本设计目标**：在现有 CZM-电化学-热耦合框架内，引入渐进式有效面积损失机制，消除死区并增强正反馈。

---

## 2. 有效面积公式

当 D > D_threshold 时，热单元的有效电化学面积按线性缩减：

```
D ≤ D_threshold:  A_eff[e] = A_e
D > D_threshold:  A_eff[e] = A_e × (1 - D[e]) / (1 - D_threshold)
D = 1.0:          A_eff[e] = 0
```

**参数**：
- `D_threshold`：面积开始缩减的损伤阈值，默认 0.83（可调）
- 公式在 D = D_threshold 处连续（A_eff = A_e）

**D_threshold 物理含义**：
- D = 0.83 对应 δ ≈ δ₀/(1-D) = 0.34nm/0.17 ≈ 2.0 nm
- 此时界面脱粘开始显著影响电接触

**有效面积对应值**：

| D | (1-D)/(1-0.83) | A_eff / A_e |
|---|----------------|-------------|
| ≤ 0.83 | 1.0 | 100% |
| 0.90 | 0.59 | 59% |
| 0.96 | 0.24 | 24% |
| 0.99 | 0.06 | 6% |
| 1.00 | 0.0 | 0% |

---

## 3. CZM → 热单元 D 映射

每个内聚力单元直接对应一个热单元（1-to-1 映射）。

**映射方式**：
- 利用 `geometry.czm_element_map`
- D_elem[e] = czm_mesh.damage_states[对应 czm_idx].D
- 无需取 max 或 mean——直接 1-to-1

**实现**：新增辅助函数 `map_czm_damage_to_thermal`

---

## 4. 四个修改位置

### 4.1 位置 1：D 映射（`CallModel.jl`，新增函数）

新增函数将 CZM damage_states 映射为热单元级别的 D 数组：

```julia
function map_czm_damage_to_thermal(czm_mesh, geometry, ne)
    D_elem = zeros(ne)
    for e in 1:ne
        czm_indices = get(geometry.czm_element_map, e, Int64[])
        if !isempty(czm_indices)
            # 1-to-1: 取第一个（也是唯一一个）
            D_elem[e] = czm_mesh.damage_states[czm_indices[1]].D
        end
    end
    return D_elem
end
```

### 4.2 位置 2：分流求解器权重（`Parallelsolution.jl`）

**当前代码**（约第 360 行）：
```julia
w = areas ./ sum(areas)
```

**修改为**：
```julia
A_eff = areas .* effective_area_factor.(D_elem, D_threshold)
w = A_eff ./ sum(A_eff)
```

其中 `effective_area_factor(D, D_th)` 实现为：
```julia
function effective_area_factor(D::Float64, D_threshold::Float64)
    D ≤ D_threshold && return 1.0
    return (1.0 - D) / (1.0 - D_threshold)
end
```

**影响**：
- D > D_threshold 的单元在分流求解器中权重降低
- 分配到这些单元的电流减少
- 剩余电流自动重分配到未损伤/低损伤单元

### 4.3 位置 3：SPMe 内部电流密度（`CallModel.jl` / `SPMe.jl`）

SPMe 模型计算电流密度时使用 A_eff 而非 A_e：

```
j_e = I_e / A_eff[e]
```

**修改方式**：
- 将 A_eff 数组传递给 `CallModel_MultiSPMe`
- 在调用 `SPMe_element` 时传入 `A_eff[e]` 作为有效面积参数
- SPMe 内部用此面积计算电流密度

**影响**：
- 电流密度在并联电路中仍对所有单元相等：j = I_total / sum(A_eff)
- 但 j > j_original（因为 sum(A_eff) < sum(A)）
- 系统级电流密度上升

### 4.4 位置 4：热源计算（`ThermalDistributed.jl`）

**不需要额外缩放因子**。热源由 SPMe 内部计算：
- Q_rxn ∝ j（反应热）
- Q_ohm ∝ j²（焦耳热，二次方放大）
- Q_rev ∝ j（可逆热）

当 SPMe 正确使用 A_eff 计算电流密度后，热源自动反映增大的 j。

### 4.5 保留 D ≥ 0.99 硬截止

`get_fractured_elements`（`Materialmatrix.jl:365`）的完全失活机制保留作为安全阀：
- D ≥ 0.99 时 A_eff ≈ 6% A_e，接近完全失活
- 完全失活（I_e = 0, q = 0）仍由原有逻辑处理

---

## 5. 反馈链分析

### 5.1 系统级焦耳热反馈

在并联电路模型中，所有活跃单元的电流密度相同：

```
j = I_total / sum(A_eff)
```

当某些单元 D > 0.83 → sum(A_eff) ↓ → j 全局上升。

**正反馈链**：
```
D_e↑ → A_eff_e↓ → sum(A_eff)↓
    → j = I_total/sum(A_eff) ↑ (所有单元)
    → Q_ohm ∝ j² ↑ (焦耳热二次方放大)
    → ΔT↑ (全局温升)
    → 热应变↑ → F_thermo_chem↑ → δ↑ → D↑ (所有单元)
```

**反馈强度估算**（20 个等面积单元，1 个单元 D=0.96）：

| 指标 | 改前 | 改后 (D=0.96) | 变化 |
|------|------|--------------|------|
| sum(A_eff) | 20A | 19.24A | -3.8% |
| j (所有单元) | I/20A | I/19.24A | +3.8% |
| Q_ohm (∝ j²) | 基准 | ×1.078 | +7.8% |
| Q_rxn (∝ j) | 基准 | ×1.038 | +3.8% |

当更多单元同时损伤时反馈更强。

### 5.2 局部热障反馈（已有机制增强）

已损伤单元的界面热阻增加（`compute_gap_conductance`）：

```
D_e↑ → h_eff_e↓ → 局部 ΔT_interface↑
    → 局部热应变额外↑ → δ_e↑ → D_e↑ (单元级正反馈)
```

在系统级 j 上升的背景下，局部热障效应被放大（更多热流通过受损界面）。

### 5.3 渐进式容量损失

```
D > 0.83 → A_eff↓ → 可参与反应的活性面积↓
    → 该单元可提取锂量↓
    → 并联模型总容量↓
```

不再有 D=0~0.99 的死区。任何 D > 0.83 的损伤都立即反映为容量能力的下降。

---

## 6. 预期效果

### 6.1 单次放电

- D 收敛值应高于 0.96（因系统级焦耳热反馈增强）
- 收敛值取决于 D_threshold、单元数、C-rate 等参数

### 6.2 多循环

- 损伤累积加速（每循环的增量反馈更强）
- 更快达到 D ≥ 0.99 触发完全失活
- 容量衰减曲线从第一个循环起就有渐进式下降

### 6.3 与现有机制的兼容性

- D < D_threshold (0.83)：无面积缩减，行为与改前完全一致
- D_threshold < D < 0.99：渐进式面积缩减（新机制）
- D ≥ 0.99：完全失活（原有机制保留作为安全阀）

---

## 7. 新增参数

在 `Option` 结构体（`src/Option.jl`）中新增：

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `czm_area_loss_enabled` | Bool | false | 启用渐进式面积损失 |
| `czm_area_loss_threshold` | Float64 | 0.83 | 面积开始缩减的 D 阈值 |

---

## 8. 代码修改汇总

| 文件 | 修改内容 | 复杂度 |
|------|----------|--------|
| `src/Option.jl` | 新增 `czm_area_loss_enabled`, `czm_area_loss_threshold` 字段 | 低 |
| `src/CallModel.jl` | 新增 `map_czm_damage_to_thermal` 函数；调用时计算 D_elem 并传递 | 中 |
| `src/Parallelsolution.jl` | `solve_branch_currents` 接收 D_elem，计算 A_eff 用于权重 | 中 |
| `src/CallModel.jl` / `SPMe.jl` | SPMe 使用 A_eff 计算电流密度 | 中 |
| `src/ThermalDistributed.jl` | 无需修改（热源由 SPMe 自动计算） | — |
| `src/Materialmatrix.jl` | 新增 `effective_area_factor` 辅助函数 | 低 |

---

## 9. 验证方案

1. **单元测试**：`effective_area_factor` 在 D=0, D=0.83, D=0.96, D=1.0 处的返回值
2. **单次放电**：对比启用/禁用面积损失的 D_max 收敛值
3. **多循环**：验证渐进式容量衰减（SOH 曲线不再是阶跃式）
4. **守恒检查**：电流守恒 sum(w .* I_e) = I_total 仍然成立
5. **回归**：D_threshold=1.0 时行为与改前完全一致（面积不缩减）
