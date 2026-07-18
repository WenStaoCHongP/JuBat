# CZM 损伤收敛分析：D_max = 0.96 的理论评估

**日期**: 2026-05-26
**状态**: 已完成
**范围**: 单次放电工况下，双线性内聚力模型损伤收敛机制与反馈链分析

---

## 1. 问题描述

在单次放电仿真中，观察到最大损伤 D_max 收敛于 0.96 而非达到 1.0（完全断裂）。
本文档从理论上分析：

1. D_max 收敛于 0.96 的数学和物理原因
2. 当前电化学-热-内聚力耦合模型是否存在反馈机制推动损伤继续增长
3. D = 0.96 是否导致容量损失

---

## 2. CZM 关键参数

来自 `src/parameters/Jellyroll.jl:148-153`：

| 参数 | 符号 | 值 |
|------|------|------|
| 最大牵引应力 | σ_max | 82 MPa |
| 界面刚度 | K_n | 2.4×10¹⁷ Pa/m |
| 损伤起始分离 | δ₀ = σ_max/K_n | **0.34 nm** |
| 断裂能 | G_c | 25.3 J/m² |
| 临界分离位移 | δ_c = 2G_c/σ_max | **617 nm** |
| 刚度比 | δ_c/δ₀ | **~1806** |

---

## 3. 定量临界值计算

### 3.1 D = 0.96 时的分离位移反算

双线性损伤公式：

$$D = \frac{\delta_c(\delta - \delta_0)}{\delta(\delta_c - \delta_0)}$$

由于 δ_c >> δ₀ (1806×)，当 δ >> δ₀ 时可近似为：

$$D \approx 1 - \frac{\delta_0}{\delta}$$

反算 D = 0.96：

$$\delta_{0.96} \approx \frac{\delta_0}{1 - D} = \frac{0.34\text{ nm}}{0.04} \approx 8.5 \text{ nm}$$

精确解：δ* = 8.43 nm。

### 3.2 损伤-位移对应表

| D 值 | 对应 δ | 与 δ_c 的差距 | 倍率 |
|------|--------|---------------|------|
| 0.96 | 8.5 nm | 608.5 nm | δ_c/δ = 72× |
| 0.99 | 34 nm | 583 nm | δ_c/δ = 18× |
| 1.00 | 617 nm | 0 | — |

**核心发现**：从 D = 0.96 到 D = 1.0，分离位移需要增加约 **72 倍**。这是 δ_c/δ₀ 比值极大（~1806）的直接后果。

---

## 4. 反馈链增益分析

### 4.1 路径 A：热阻反馈（正反馈，增益极弱）

**反馈链**：
```
D↑ → h_eff↓ → R_th↑ → ΔT_interface↑ → Δε_thermal↑ → F_thermo_chem↑ → δ↑ → D↑
```

**实现位置**：`src/Materialmatrix.jl` (`compute_gap_conductance`)、`src/ThermalDistributed.jl` (`ThermalDistributed2D_BC`)

**D = 0.96 时的界面导热系数**：

完好界面 (D = 0)：
```
h_eff ≈ h_c0 + k_air/(2βλ_m) = 1e7 + 1.86e5 ≈ 1.02×10⁷ W/(m²·K)
```

损伤界面 (D = 0.96, δ = 8.5 nm)：
```
h_eff = h_c0(1-D) + k_air/(δ + 2βλ_m)
      = 1e7 × 0.04 + 0.026/(8.5e-9 + 140e-9)
      = 4.0e5 + 1.75e5
      ≈ 5.75×10⁵ W/(m²·K)
```

**降低约 17.5 倍**。

界面温差估算（1C 放电）：
- 体热源 q ≈ 10⁴–10⁵ W/m³
- 界面热流 q_interface ≈ q × t_layer ≈ 10⁴ × 100e-6 = 1 W/m²
- ΔT_interface ≈ 1 / 5.75e5 ≈ **2×10⁻⁶ K**

**结论**：界面热阻增加导致的温差增量在微开尔文量级，远不足以产生有意义的热应变增量。

### 4.2 路径 B：本体刚度软化（自洽效应）

在软化段，CZM 牵引力：
```
T = (1-D) × K_n × δ
```

D = 0.96 时，有效刚度降为 4% K_n。这是 CZM 求解器已在平衡方程中考虑的自洽效应，不提供额外驱动力。

### 4.3 路径 C：电流重分配（D < 0.99 时未激活）

**代码依据** (`src/Materialmatrix.jl:365`, `src/CallModel.jl:63`):
```julia
# get_fractured_elements: state.fractured || state.D >= 0.99
```

D = 0.96 的单元：
- 不被 deactivate
- 仍承载电流
- 仍产生热源
- 无电流集中效应

**此路径在 D < 0.99 时完全不激活。**

### 4.4 路径 D：热源屏蔽（仅在断裂时激活）

**代码依据** (`src/ThermalDistributed.jl:536-541`):
```julia
for e in 1:ne
    if !is_active[e]
        q_total[e] = 0.0
    end
end
```

通过 `get_active_elements` → `get_fractured_elements` (D ≥ 0.99) 触发。D = 0.96 不受影响。

### 4.5 路径 E：化学应变（恒定速率，无增量驱动）

**代码依据** (`src/czm.jl:397-451`, `src/CouplingState.jl:262`):
```
ε₀ = α_eff × ΔT + β_n × Δsoc_n + β_p × Δsoc_p
```

单次放电中，SOC 单调下降，Δsoc 变化率恒定。化学应变不提供增量驱动。

### 4.6 反馈链汇总

| 反馈路径 | 增益量级 | 是否激活 | 结论 |
|----------|----------|----------|------|
| A: 热阻反馈 | ~10⁻⁶ K 级温差 | 激活但极弱 | 不足以推动 |
| B: 刚度软化 | 自洽效应 | 已包含在求解器中 | 非驱动因素 |
| C: 电流重分配 | 有效 | **未激活** (需 D ≥ 0.99) | 关键路径缺失 |
| D: 热源屏蔽 | 有效 | **未激活** (需 D ≥ 0.99) | 关键路径缺失 |
| E: 化学应变 | 恒定速率 | 激活 | 无增量驱动 |

---

## 5. D = 0.96 与容量损失

### 5.1 当前模型的容量损失机制

**代码依据** (`src/PostProcessing.jl:251-257`):
```julia
current_soh = initial_capacity > 0 ?
    cycle_result.capacity_discharge / initial_capacity : 1.0
```

容量损失仅通过单元失活实现：
1. `get_fractured_elements` 判定 D ≥ 0.99 的单元
2. 失活单元电流归零，热源归零
3. 总电流重分配至剩余活跃单元
4. 活跃单元数量减少 → 可放电容量降低

### 5.2 D = 0.96 对容量的影响

**直接影响：无**。D = 0.96 的单元：
- 仍承载电流（未被 deactivate）
- 仍产生热源（未被屏蔽）
- 电化学状态完整参与

**间接影响：可忽略**。界面热阻增加 17.5 倍导致的微开尔文级温差对反应速率的影响可忽略。

---

## 6. 收敛机制总结

### 6.1 为什么 D_max 收敛于 0.96

在单次放电过程中：
1. 温度趋于稳态：热生成与散热达到平衡，ΔT 增长率 → 0
2. SOC 变化率恒定：恒流放电下 Δsoc 变化率大致不变
3. CZM 驱动力稳定：F_thermo_chem 不再显著增长
4. 分离位移收敛：δ 稳定于 δ* ≈ 8.5 nm，对应 D ≈ 0.96

### 6.2 根本原因

当前模型在 D = 0.96 到 D = 0.99 之间存在一个 **"死区"**：
- 没有耦合机制能提供足够的驱动力跨越 8.5 nm → 617 nm 的间隙
- 唯一的断裂触发条件（D ≥ 0.99 → 单元失活 → 电流重分配 → 热集中）被阻隔在这个死区之后

### 6.3 真实电池中可能存在的反馈（当前模型未实现）

1. **SEI 生长加速**：界面损伤 → 局部 SEI 生长加速 → 体积膨胀 → 界面应力增加
2. **电解液浸润性下降**：界面脱粘 → 电解液难以浸润 → 反应不均匀 → 电流密度集中
3. **锂沉积**：界面损伤 → 局部电流密度不均匀 → 锂枝晶生长 → 进一步损伤
4. **渐进式容量衰减**：即使未完全断裂，部分脱粘也会降低有效反应面积

---

## 7. D = 1 后受流面积减小的传导分析

### 7.1 从断裂到容量损失的完整传导链

```
D ≥ 0.99 (断裂判定)
    │  get_fractured_elements (Materialmatrix.jl:365)
    ▼
映射到热单元 → deactivated_elements (CallModel.jl:63-74)
    │  遍历 czm_element_map，将有断裂 CZM 的热单元标记为 deactivated
    ▼
分流求解器电流重分配 (Parallelsolution.jl:394-408)
    │  active_mask 排除 deactivated 单元
    │  I_e[deactivated] = 0
    │  归一化: sf = I_total / sum(w[active] .* I_e[active])
    │  → 每个活跃单元电流增大
    ▼
SPMe 电化学响应 (各单元独立求解)
    │  活跃单元 I_e ↑ → 过电位 ↑ → 端电压 ↓
    ▼
端电压截止 (Solve.jl:357-364)
    │  V_cell < v_l → termination_reason = "voltage_cutoff_low"
    │  → 放电提前终止，duration 缩短
    ▼
容量计算 (PostProcessing.jl:187)
    capacity = I_total × duration / 3600
    → duration 更短 → 容量更小
```

### 7.2 面积权重的关键细节

分流求解器中，面积权重 `w = areas ./ sum(areas)` 基于 **全部** 单元面积计算一次（含未来将被 deactivate 的单元）。

当 n_deact 个单元断裂后，剩余活跃单元承载的总电流约束：

```
sum(w[e] × I_e[e] for e in active) = I_total
```

对于面积相等的 n 个单元，断裂 n_deact 个后，每个活跃单元电流：

$$I_{e,\text{active}} = I_{\text{total}} \times \frac{n}{n - n_{\text{deact}}}$$

| 断裂比例 | 电流密度增幅 | 预期容量影响 |
|----------|-------------|-------------|
| 10% | +11% | 小幅下降 |
| 30% | +43% | 明显下降 |
| 50% | +100% | 严重下降 |

### 7.3 容量损失是"间接"的

关键代码（`PostProcessing.jl:187`）：
```julia
capacity = abs(I_current) * duration / 3600.0
```

这里 `I_current` 是恒定的总电流（如 5A），`duration` 是实际放电持续时间。模型 **不直接** 计算 `活跃面积 × SOC = 可提取锂量`，而是通过以下间接路径：

1. 活跃单元电流密度 ↑ → Butler-Volmer 过电位 ↑ → 单元端电压 ↓
2. 公共端电压 V（所有活跃单元的并联电压）下降更快
3. `V_cell < V_lower` 更早触发 → `duration` 更短
4. `capacity = I_total × duration` 更小

**结论：受流面积减小的影响通过电压响应的间接路径得到了体现，物理上是正确的。**

### 7.4 当前模型的两个关键限制

#### 限制 1：二元开关，无渐变过渡

```
D < 0.99  →  单元完全活跃（100% 面积参与）
D ≥ 0.99  →  单元完全失活（0% 面积参与）
```

代码依据（`Materialmatrix.jl:365`）：
```julia
if state.fractured || state.D >= 0.99
    push!(fractured, i)
```

在真实电池中，部分脱粘（0 < D < 0.99）也应该导致：
- 接触电阻增加（有效导通面积减小）
- 反应面积减小（部分活性物质失联）
- 局部电流密度不均匀

**当前模型在这些中间状态下不体现任何容量效应。** D = 0.96 与 D = 0 的单元在电化学上完全等价。

#### 限制 2：面积权重不随断裂更新

```julia
w = areas ./ sum(areas)  # 始终用全部面积做分母
```

即使 50% 的单元断裂，`sum(areas)` 仍然是原始总面积。归一化通过 `sf = I_total / sum(w[active] .* I_e[active])` 保证总电流约束，但权重分母不变。

这在物理上意味着：**面积参考是固定的（初始总电极面积），断裂单元的面积被"短路"掉而非从总面积中扣除。**

对于并联模型来说，这种处理等价于：
```
I_e[active] = I_total × w[e] / sum(w[active])
```

即活跃单元按原始面积比例分担总电流，在并联电路假设下是正确的。

### 7.5 "陷阱锂"效应

当单元断裂（D ≥ 0.99）被 deactivate 时：
- 其 SPMe 状态被冻结（I_e = 0 → 无电化学反应）
- 单元内剩余的嵌入锂无法被提取
- 这部分锂代表直接的容量损失

模型 **通过电压响应间接体现了这一效应**：可提取的总锂量减少 → 各活跃单元 SOC 消耗更快 → 电压下降更快 → 放电时间更短 → 容量更小。

但如果需要直接量化"因断裂损失的活性锂量"，当前模型没有提供这种计算。可以考虑增加：
```
Q_trapped = sum(A_e × (c_s,e - c_s,empty) × F / 3600 for e in deactivated)
```

### 7.6 D = 1 后受流面积分析总结

| 问题 | 回答 |
|------|------|
| D = 1 后受流面积减小是否体现？ | **是**，通过 deactivate → 电流重分配 → 电压下降加速 → 提前截止 → 容量减小 |
| 体现方式是直接的还是间接的？ | **间接的**，通过电压响应而非面积计算 |
| D = 0.96 时受流面积是否减小？ | **否**，0 < D < 0.99 的单元仍然 100% 活跃 |
| "陷阱锂"是否被量化？ | **否**，仅通过放电时间缩短间接反映 |
| 面积权重是否随断裂更新？ | **否**，w 分母始终为总面积，但归一化保证了正确性 |

---

## 8. 代码引用索引

| 机制 | 文件 | 函数/行号 |
|------|------|-----------|
| 双线性本构 | `src/Materialmatrix.jl` | `bilinear_traction_state` (68-161) |
| 损伤更新 | `src/Materialmatrix.jl` | `update_damage` (293-307) |
| 断裂判定 | `src/Materialmatrix.jl` | `get_fractured_elements` (362, D ≥ 0.99) |
| 活跃单元 | `src/Materialmatrix.jl` | `get_active_elements` (375) |
| 界面导热 | `src/Materialmatrix.jl` | `compute_gap_conductance` (310) |
| CZM 调度 | `src/CouplingState.jl` | `update_czm_damage!` (329-430) |
| 有效参数 | `src/CouplingState.jl` | `compute_czm_effective_params` (228) |
| 应变输入 | `src/CouplingState.jl` | `compute_czm_strain_inputs` (262) |
| 热-化学载荷 | `src/czm.jl` | `assemble_thermal_chemical_load` (397-451) |
| CZM 求解入口 | `src/CzmSolve.jl` | `solve_czm_step` (662) |
| 电流重分配 | `src/CallModel.jl` | CZM 失效处理 (60-78) |
| 分流求解 | `src/Parallelsolution.jl` | `solve_branch_currents` (358) |
| 热源屏蔽 | `src/ThermalDistributed.jl` | `compute_heat_sources_with_czm` (512-550) |
| SOH 计算 | `src/PostProcessing.jl` | `_update_soh_and_capacity!` (251-260) |
| CZM 参数 | `src/parameters/Jellyroll.jl` | cohesive (146-171) |
