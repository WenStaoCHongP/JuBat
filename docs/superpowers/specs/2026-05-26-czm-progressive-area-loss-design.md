# CZM 渐进式有效面积损失机制设计

**日期**: 2026-05-26
**状态**: 审阅修正后
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
- 精确公式反解：δ = δ₀·δ_c / (δ_c - D·(δ_c - δ₀))
- D = 0.83 时：δ = 0.34e-9 × 617e-9 / (617e-9 - 0.83 × 616.66e-9) = 2.095e-9 / 105.83e-9 ≈ 1.98 nm
- 约等于近似值 δ₀/(1-D) = 2.0 nm（因 δ_c >> δ₀，近似误差 <5%）
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

**核心假设：一个内聚力单元只影响对应的内侧热单元，一对一。**

因此，热单元的损伤仅由其对应的那个 CZM 单元决定，无需取 maximum 或 mean。

依据：
- `create_czm_mesh`（`czm.jl:87-110`）中，每个 CZM 单元的底面节点 `(n_in_1, n_in_2)` 恰好对应唯一一个内侧热单元
- `get_active_elements`（`Materialmatrix.jl:381`）只影响 `is_inner_layer[e] = true` 的热单元

**映射实现**：遍历内侧热单元，直接取其对应 CZM 单元的 D 值。

```julia
function map_czm_damage_to_thermal(czm_mesh, geometry, ne)
    D_elem = zeros(ne)
    for e in 1:ne
        czm_indices = get(geometry.czm_element_map, e, Int64[])
        if !isempty(czm_indices) && geometry.is_inner_layer[e]
            # 1-to-1: 内侧热单元只对应一个 CZM 单元，直接取其 D
            D_elem[e] = czm_mesh.damage_states[czm_indices[1]].D
        end
    end
    return D_elem
end
```

---

## 4. 修改位置

### 4.1 位置 1：D 映射（`CallModel.jl`，新增函数）

如第 3 节所述，新增 `map_czm_damage_to_thermal` 函数。仅遍历内侧热单元，直接取其对应 CZM 单元的 D（1-to-1）。

### 4.2 位置 2：分流求解器权重（`Parallelsolution.jl`）

**这是唯一的干预点。** 仅修改分流求解器的面积权重，后续效应全部自动传导。

**当前代码**（约第 360 行）：
```julia
w = areas ./ sum(areas)
```

**修改为**：
```julia
if czm_area_loss_enabled
    D_elem = map_czm_damage_to_thermal(czm_mesh, geometry, ne)
    A_eff = areas .* effective_area_factor.(D_elem, D_threshold)
    w = A_eff ./ sum(A_eff)
else
    w = areas ./ sum(areas)
end
```

其中 `effective_area_factor(D, D_th)` 实现为：
```julia
function effective_area_factor(D::Float64, D_threshold::Float64)
    D ≤ D_threshold && return 1.0
    return (1.0 - D) / (1.0 - D_threshold)
end
```

**影响传导链**（全部自动，无需额外修改）：

```
w[e]↓ (受损单元权重降低)
  → Newton 迭代调整 I_e 分配：
    受损单元 I_e↓，健康单元 I_e↑
  → SPMe: j_n = I_e / (as × thickness) 自动变化
    受损单元 j_n↓，健康单元 j_n↑
  → 热源: 两条路径均自动响应
    欧姆热: Q_ohm = I_local² / (3σ) → 健康单元↑
    反应热: Q_rxn = as × |j_n| × |η| → 健康单元↑
  → 温度场 → 热应变 → CZM 驱动力 → δ → D
```

### 4.3 SPMe 和热源：无需修改

**SPMe 电流密度**（`SPMe.jl:111`）：
```julia
j_n = I_app / param.NE.as / param.NE.thickness  # 体积电流密度 (A/m³)
```
这是体积电流密度，输入 `I_app` 来自分流求解器输出的 `I_e[e]`。权重修改后 `I_e[e]` 自动变化，`j_n` 随之变化。**不需要额外修改 SPMe。**

**热源计算**（`ThermalDistributed.jl:455-471`）使用两条混合路径：

| 热源分量 | 计算方式 | 对 I_e 变化的响应 |
|----------|----------|-------------------|
| Q_ohm_s (固相欧姆) | `I_local² / (3σ)` | 直接响应 I_e² |
| Q_ohm_e (液相欧姆) | `I_local² / (3κ)` | 直接响应 I_e² |
| Q_SP (隔膜) | `I_local² / κ_sp` | 直接响应 I_e² |
| Q_PCC/NCC (集流体) | `I_local² / (3σ)` | 直接响应 I_e² |
| Q_rxn (反应热) | `as × \|j_n\| × \|η\|` | 通过 j_n 间接响应 |
| Q_rev (可逆热) | `as × j_n × T × dU/dT` | 通过 j_n 间接响应 |

其中 `I_local = I_e[e]`。权重修改后 `I_e[e]` 变化 → 两条路径均自动响应。**不需要额外修改热源计算。**

### 4.4 保留 D ≥ 0.99 硬截止

`get_fractured_elements`（`Materialmatrix.jl:365`）的完全失活机制保留作为安全阀：
- D ≥ 0.99 时 A_eff ≈ 6% A_e，接近完全失活
- 完全失活（I_e = 0, q = 0）仍由原有逻辑处理

---

## 5. 反馈链分析

### 5.1 通过分流求解器的反馈机制

分流求解器的 Newton 迭代求解：
- `V_e(I_e[e]) = V`（所有活跃单元等电压）
- `sum(w[e] × I_e[e]) = I_total`（电流守恒）

当受损单元 w[e] 减小，Newton 迭代调整 V 和所有 I_e[e] 的分配：

**正反馈链**：
```
D_e↑ → A_eff_e↓ → w[e]↓
    → Newton 迭代：受损单元 I_e↓，健康单元 I_e↑
    → 健康单元 I_local²↑ → Q_ohm↑ (欧姆热二次方放大)
    → 健康单元温升↑ → 热应变↑ → CZM 驱动力↑
    → δ↑ → D↑ (所有单元，包括已损伤的)
```

**反馈强度**取决于：
- 受损单元数量和严重程度
- 健康单元对额外电流的热响应灵敏度
- 具体的 Newton 迭代结果（非线性，无法简化为线性公式）

**定量分析需通过实际仿真验证**，不建议使用简化公式估算。

### 5.2 局部热障反馈（已有机制增强）

已损伤单元的界面热阻增加（`compute_gap_conductance`）：

```
D_e↑ → h_eff_e↓ → 局部 ΔT_interface↑
    → 局部热应变额外↑ → δ_e↑ → D_e↑ (单元级正反馈)
```

在健康单元 I_e 增加的背景下（更多热生成），局部热障效应被放大。

### 5.3 渐进式容量损失

```
D > D_threshold → w[e]↓ → I_e[e]↓
    → 受损单元承载电流减少 → 可提取电荷减少
    → 健康单元承载更多电流 → 更快达到电压截止
    → 总放电时间缩短 → capacity = I_total × duration↓
```

不再有 D=0~0.99 的死区。任何 D > D_threshold 的损伤都通过电压响应间接反映为容量损失。

---

## 6. 预期效果

### 6.1 单次放电

- D 收敛值应高于 0.96（因正反馈增强）
- 收敛值取决于 D_threshold、单元数、C-rate 等参数
- 具体数值需通过仿真验证

### 6.2 多循环

- 损伤累积加速（每循环的增量反馈更强）
- 更快达到 D ≥ 0.99 触发完全失活
- 容量衰减曲线从第一个循环起就有渐进式下降

### 6.3 与现有机制的兼容性

- `czm_area_loss_enabled = false`：行为与改前完全一致（回归保证）
- D < D_threshold (0.83)：无面积缩减
- D_threshold < D < 0.99：渐进式面积缩减（新机制）
- D ≥ 0.99：完全失活（原有机制保留作为安全阀）

---

## 7. 新增参数

在 `Option` 结构体（`src/Option.jl`）中新增：

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `czm_area_loss_enabled` | Bool | false | 启用渐进式面积损失 |
| `czm_area_loss_threshold` | Float64 | 0.83 | 面积开始缩减的 D 阈值（可调） |

**D_threshold 调参建议**：
- 建议范围：0.70–0.90
- 较低值（0.70）：更早启动面积缩减，反馈更强，但可能过度惩罚轻微损伤
- 较高值（0.90）：仅在严重脱粘时启动，反馈更温和
- 推荐通过参数扫描对比（0.75, 0.80, 0.83, 0.85, 0.90）确定最优值

---

## 8. 代码修改汇总

| 文件 | 修改内容 | 复杂度 |
|------|----------|--------|
| `src/Option.jl` | 新增 `czm_area_loss_enabled`, `czm_area_loss_threshold` 字段 | 低 |
| `src/CallModel.jl` | 新增 `map_czm_damage_to_thermal` 函数；计算 D_elem 传递给分流求解器 | 中 |
| `src/Parallelsolution.jl` | `solve_branch_currents` 新增 `D_elem` 参数，计算 A_eff 用于权重 | 中 |
| `src/Materialmatrix.jl` | 新增 `effective_area_factor` 辅助函数 | 低 |
| `src/SPMe.jl` | **无需修改**（I_e 自动反映权重变化） | — |
| `src/ThermalDistributed.jl` | **无需修改**（I_local 自动反映 I_e 变化） | — |

---

## 9. 验证方案

1. **单元测试**：`effective_area_factor` 在 D=0, D=0.83, D=0.96, D=1.0 处的返回值
2. **映射测试**：`map_czm_damage_to_thermal` 对 1-to-many 情况取 maximum 的正确性
3. **单次放电**：对比启用/禁用面积损失的 D_max 收敛值
4. **多循环**：验证渐进式容量衰减（SOH 曲线不再是阶跃式）
5. **能量守恒**：总产热功率与电功率输入的平衡检查
6. **回归**：`czm_area_loss_enabled = false` 时行为与改前完全一致
7. **衔接测试**：D 从 0.98 → 0.99 时渐进面积损失与原有硬截止的平滑过渡
8. **参数扫描**：D_threshold 在 0.75–0.90 范围内的灵敏度分析
