# CZM 耦合修复与单步损伤演化

> 日期: 2026-04-20
> 更新时间: 2026-04-20
> 分支: Parameters_Design
> 涉及文件: 12 个 src/ 文件, +459 / -228 行

---

## 1. 问题背景

`example/testexample.jl` 和 `example/coupled_czm_thermal_example.jl` 在 `opt.czm_enabled = true` 时无法运行，暴露出 CZM（内聚力模型）与主求解器之间多处接口断裂。此外，CZM 损伤仅在循环求解（`solve_cycling`）的阶段末更新，单次放电（`Solve`）中损伤不演化。

---

## 2. 修改总览

### 2.1 API 修正（example 脚本 → src 接口对齐）

| 问题 | 修复 |
|------|------|
| `param_dim.gap_conductance.h_0` 不存在 | → `param_dim.cohesive.h_c0`，间隙导热参数在 `Cohesive` 结构体中 |
| `mesh_data.Jellyroll_czm` 不存在 | → `mesh_data.thermal2D`，`JellyrollMesh` 无 `Jellyroll_czm` 字段 |
| `czm_enabled = true` 但未创建 CZM 网格 | 在 `testexample.jl` 中添加 `create_czm_mesh` 调用 |
| `compute_all_gap_conductances(czm_mesh, param_dim.gap_conductance)` | → `param_dim.cohesive` |

### 2.2 数据流修复（模块间参数传递）

| 位置 | 问题 | 修复 |
|------|------|------|
| `Solve.jl:93` | `initial_state` 为 `Dict` 时 `vec()` 崩溃 | 添加 Dict 分支，提取 `"y"` 键 |
| `CycleSolver.jl:142` | CZM 网格未挂载到 `case.czm_mesh` | `solve_phase` 中 `case.czm_mesh = czm_mesh` |
| `ThermalDistributed.jl:498` | `get_active_elements` 收到 `Mesh` 而非 `MeshGeometry` | 从 `case.geometry` 获取 `is_inner_layer` / `czm_element_map` |
| `Materialmatrix.jl:346` | `mesh_data.ne` 字段不存在 | 改用 `length(mesh_data.is_inner_layer)` |
| `ThermalDistributed.jl:523` | `active_elements` 为 `Vector{Int64}` 无法写入 `Dict{String, Union{Float64, Array{Float64}}}` | 转换为 `Float64.(active_elements)` |

### 2.3 CZM 损伤演化集成到 Solve 主循环

**核心变更**：将 CZM 损伤更新从"仅循环求解阶段末"扩展到"单次 Solve 每步"。

```
调用链（修改前）:
  Solve.jl 主循环 → CallModel → [无 CZM 损伤更新]
  CycleSolver.jl → solve_phase → Solve → 阶段末 update_czm_damage!

调用链（修改后）:
  Solve.jl 主循环 → CallModel → 每步 update_czm_damage!
  CycleSolver.jl → solve_phase → Solve → 每步自动更新（复用同一逻辑）
```

**具体修改**：

| 文件 | 内容 |
|------|------|
| `CzmSolve.jl` | 从 `CycleSolver.jl` 迁入 `compute_czm_effective_params`、`compute_czm_strain_inputs`、`update_czm_damage!` 三个函数，确保在 `Solve.jl` 之前加载 |
| `CycleSolver.jl` | 删除上述三个函数（161 行），改为调用 `CzmSolve.jl` 中的定义 |
| `Solve.jl` | 主循环中每步接受后按 `czm_update_interval` 间隔调用 `update_czm_damage!`，并追踪耗时 |

### 2.4 CZM 摘要统计记录与输出

使 CZM 损伤数据在每步可追踪、可绘图。

| 文件 | 新增内容 |
|------|----------|
| `Variables.jl` | `StandardVariables` 新增: `czm D_max`, `czm D_mean`, `czm δ_max_n`, `czm δ_mean_n`, `czm n_fractured` |
| `CallModel.jl` | `CallModel_MultiSPMe` 中从 `czm_mesh.damage_states` 计算摘要统计写入 `variables` |
| `PostProcessing.jl` | 输出 `czm D_max`, `czm D_mean`, `czm δ_max_n [m]`, `czm δ_mean_n [m]`, `czm n_fractured` |

### 2.5 计时器修复

CZM 损伤更新（`update_czm_damage!`）在 Solve 主循环中调用但未计入 `timing_totals`，导致计时统计与墙钟时间严重不符。

**修复**：在 `Solve.jl` 中用 `time_ns()` 包裹 `update_czm_damage!` 调用，累加到 `timing_totals["czm"]`。

**效果**（60s 仿真，18 步）：

| 模块 | 修复前 | 修复后 |
|------|--------|--------|
| SPMe 求解 | 63.34% | 0.66% |
| CZM 模型 | 1.89% | **98.90%** |
| 墙钟 / 合计 | 19s / 10s（不符） | 225s / 221s（吻合） |

### 2.6 其他小修复

| 文件 | 修改 |
|------|------|
| `Assemble.jl` | `KJ = deepcopy(KI)` → `KJ = zeros(Int64, length(KI))`，消除不必要的深拷贝 |
| `CouplingState.jl` | `MultiSPMeLayout` 新增 `areas` 字段及构造器；新增 `BoundaryEdgeCache` 结构体和 `compute_boundary_edge_cache` 函数 |
| `MeshGeometry` | 新增 `boundary_edges` 字段 |
| `Jellyrollmodel.jl` | `setup_thermal2D_mesh` 中填充 `boundary_cache` |
| `Initialisation.jl` | 线程局部工作区使用精简型 `create_element_workspace` |
| `ThermalDistributed.jl` | 大量重构：预缓存网格数据、BC 缓存、精简线程工作区等 |

---

## 3. 文件变更明细

```
 src/Assemble.jl           |   2 +-
 src/CallModel.jl          |  10 ++-
 src/CouplingState.jl      |  58 +++++++++++++-
 src/CycleSolver.jl        | 166 ++--------------------------------
 src/CzmSolve.jl           | 165 +++++++++++++++++++++++++++++++++
 src/Initialisation.jl     |   2 +-
 src/Jellyrollmodel.jl     |   8 +-
 src/Materialmatrix.jl     |   2 +-
 src/PostProcessing.jl     |   9 ++-
 src/Solve.jl              |  59 +++++++++++---
 src/ThermalDistributed.jl | 199 ++++++++++++++++++++++++++++++++------
 src/Variables.jl          |   7 +-
 12 files changed, 459 insertions(+), 228 deletions(-)
```

---

## 4. 已知问题与后续优化

1. **CZM 性能瓶颈**：`solve_czm_step` 每步 ~12s（Newton-Raphson + 20 载荷子步 × 801 CZM 单元），占 98.9% 耗时。可通过 `czm_update_interval` 增大间隔、减少 `czm_load_steps`、或向量化 CZM 求解器缓解。
2. **损伤值在短时仿真中为零**：60s 放电损伤尚未萌生，需更长仿真时间或循环工况才能观察到损伤演化。

---

## 5. 理论正确性评审

### 5.1 循环求解中 CZM 损伤被双重更新（本次引入 · 需修复）

**问题**：本次修改使 `Solve` 主循环每步调用 `update_czm_damage!`，但 `CycleSolver.solve_phase` 在 `Solve` 返回后**再次**调用 `update_czm_damage!`（L158-168）。当 `solve_cycling` 运行时，每个阶段末尾 CZM 损伤被更新两次。

```julia
# CycleSolver.jl L157-168（阶段末冗余调用）
if czm_mesh !== nothing && czm_params !== nothing
    y_end = final_state["y"]
    ...
    update_czm_damage!(czm_mesh, czm_params, case, vars_end, T_nodes_nd, nothing)
    #                              传入 u_czm_prev = nothing ← 从零开始！
end
```

**风险**：虽然 CZM 损伤只能单调递增（不会因冗余调用产生"负损伤"），但冗余调用传入 `u_czm_prev = nothing` 导致 CZM 从零位移重解，而非延续 Solve 最后一步的位移场。理论上可能产生不一致的应力-位移解，且浪费 ~12s/次的计算。

**建议**：`Solve` 已内置 CZM 损伤演化后，移除 `CycleSolver.solve_phase` 中的阶段末冗余调用。

### 5.2 扩散应变系数 β 与 cs_max 的关系（旧版问题已修复）

**位置**：`CzmSolve.jl` → `compute_czm_effective_params` + `czm.jl:363`

**分析**：

旧版分析成立的前提是：`CzmSolve.jl` 直接从 `param_dim` 读取 `Omega`，而 `Δsoc_n` / `Δsoc_p` 已经是无量纲 SOC 差。这样会少乘一个 `cs_max`，导致量纲不闭合。

当前代码已改为与 `SetParams.NormaliseParam` 一致：`Omega` 在归一化阶段已经乘上 `cs_max`，`CzmSolve.jl` 也改为使用 `case.param` 的归一化参数。因此当前路径下，下面这个应变写法是自洽的。

CZM 初始应变公式（`czm.jl:363`）：

```julia
ε_0 = α_eff × dT + β_n × Δsoc_n + β_p × Δsoc_p
```

其中：
- `β_n = Ω_n / 3 = (3.1e-6 × 33133) / 3 ≈ 0.0342` [-]
- `Δsoc_n = mean(cs/cs_max) − cs0/cs_max` [无量纲]
- `β_n × Δsoc_n` 量纲 = 无量纲应变

**对应公式**为：

$$\varepsilon_\text{diff} = \frac{\Omega}{3} \cdot \Delta c = \frac{\Omega}{3} \cdot c_\text{max} \cdot \Delta\!\left(\frac{c}{c_\text{max}}\right)$$

在当前归一化链路中，可等价写成 `β = (Ω/3) × c_s,max`，因为 `Ω` 已经在 `SetParams.NormaliseParam` 中乘过 `cs_max`。

**结论**：

- 旧版 `param_dim` 路径下，这个问题成立。
- 当前 `SetParams` + `CzmSolve` 路径下，这个问题已修复，不应再把 `β` 视为缺少 `cs_max`。

**数值对照**（以 Jellyroll NE 参数为例）：

| 项 | 当前代码 | 正确值 |
|----|----------|--------|
| β_n | 0.0342 | 0.0342 |
| Δsoc ~ 0.5 时的扩散应变 | 0.017 | 0.017 |

**说明**：原先“低估约 4 个数量级”的结论只适用于旧版 `param_dim` 路径；在当前归一化链路下不再成立。

### 5.3 热应变 α × ΔT 量纲（正确）

- `α_eff` = 加权平均热膨胀系数 [-]（`SetParams.NormaliseParam` 中按 `T_ref` 归一化）
- `dT_elem` = 无量纲温度差 `T* - T0*` [-]（当前 `CzmSolve.jl` 与归一化链路一致）
- `α_eff × dT_elem` = 无量纲应变 ✓

### 5.4 CZM 准静态假设合理性（正确）

每步调用 CZM 求解器采用 Newton-Raphson 求解力学平衡，忽略惯性效应。对电池电化学-热-力耦合问题，力学响应速率远快于热/电化学过程，准静态假设物理上合理。

### 5.5 `czm_update_interval` 的近似影响（可接受）

当 `interval > 1` 时，中间步使用陈旧的损伤状态计算间隙导热边界条件。这是可控近似，误差取决于损伤演化速率与步长之比。对于损伤缓慢累积的场景（典型电池工况），影响可忽略。

### 5.6 评估结论

| 编号 | 问题 | 严重性 | 来源 | 状态 |
|------|------|--------|------|------|
| 5.1 | 循环求解双重 CZM 更新 | 中 | 本次引入 | 待修复 |
| 5.2 | 扩散应变缺少 cs_max | **高** | 旧版路径 | 已修复 |
| 5.3 | 热应变 α×ΔT | — | 正确 | ✓ |
| 5.4 | 准静态假设 | — | 合理 | ✓ |
| 5.5 | update_interval 近似 | 低 | 可接受 | ✓ |
