# 电化学-热-内聚力耦合实现现状

本文档按当前代码实现更新，用于说明已经落地的能力、仍处于半接线状态的部分，以及与原始设计稿的差异。

---

## 一、当前总体状态

当前代码已经具备以下基础能力：

- 果冻卷二维热网格生成，以及未合并/合并节点两套热网格并行返回。
- 基于热网格坐标重合关系创建 CZM 界面单元。
- 基于双线性牵引-分离律的 CZM 求解与损伤演化。
- 在二维热边界条件装配阶段，根据 CZM 损伤状态附加层间界面导热耦合。
- 循环仿真中记录 SOH，并在 SOH 低于阈值或断裂比例过高时提前终止。

当前代码尚未完全实现或尚未接入主流程的部分：

- 原设计稿中的独立 `GapConductance` 结构并不存在，间隙导热参数已并入 `Cohesive`。
- `create_czm_mesh_from_thermal`、`ThermalDistributed2D_with_CZM`、`solve_branch_currents_with_damage` 等设计接口未落地为同名正式接口。
- `czm_update_interval` 选项已经加入 `Option`，但当前主流程未使用该参数做“每若干时间步更新一次损伤”的控制。
- 失效单元退出电化学反应的支撑代码已经存在，但 `variables["deactivated_elements"]` 目前没有在主流程中自动写入，因此“断裂即退出”的链路尚未完全闭环。

---

## 二、当前代码中的关键设计决策

| 项目 | 当前实现 |
|------|----------|
| CZM 启用开关 | `Option.czm_enabled`，默认值为 `true` |
| 损伤判据 | 以 `DamageState.D` 和 `DamageState.fractured` 共同表征；辅助函数中将 `D >= 0.99` 视为断裂 |
| CZM 模式 | `Option.czm_model` 支持 `"model1"` 和 `"mix"`；并非仅保留设计稿中的 Mode I-only 开关 |
| 间隙导热参数位置 | 存在于 `SetParams.jl` 的 `Cohesive` 结构内 |
| 热耦合位置 | 在 `ThermalDistributed2D_BC` 中直接叠加界面导热刚度 |
| 热网格选择 | `setup_thermal2D_mesh` 默认在 `czm_enabled=false` 时选 merged 网格，在 `czm_enabled=true` 时选未合并网格 |
| 循环终止 | SOH 阈值终止；或断裂单元超过 50% 时提前终止 |
| 损伤更新时间 | 当前在 `solve_phase` 的阶段末调用 `_update_czm_damage!`，不是每个时间步 |

---

## 三、数据结构与参数

### 3.1 `Option` 中已经存在的 CZM 相关字段

当前 `src/Option.jl` 中的主要字段如下：

```julia
czm_enabled::Bool = true
czm_update_interval::Int64 = 1
czm_soh_threshold::Float64 = 0.8
czm_inner_exit_only::Bool = true
czm_iter_method::String = "load_substep"
czm_max_iter::Int64 = 100
czm_tol::Float64 = 1e-4
czm_load_steps::Int64 = 20
czm_arc_length_alpha::Float64 = 1.0
```

说明：

- `czm_update_interval` 目前只定义了选项，尚未在主求解循环中消费。
- `czm_inner_exit_only` 表达了“仅内圈退出”的设计意图，但当前退出逻辑仍属于部分接线状态。

### 3.2 `Cohesive` 中已包含界面导热参数

当前 `src/SetParams.jl` 没有 `GapConductance` 结构，而是在 `Cohesive` 中直接加入了界面导热相关参数：

```julia
h_c0::Float64 = 1e7
k_air::Float64 = 0.026
lambda_m::Float64 = 70e-9
beta::Float64 = 1.0
threshold::Float64 = 70e-9
```

对应的 Jellyroll 默认值在 `src/parameters/Jellyroll.jl` 中已经给定。

### 3.3 `CohesiveMesh` 当前结构

当前 `src/czm.jl` 中的 `CohesiveMesh` 没有设计稿里提出的 `inner_elements`、`outer_elements`、`element_active`、`fractured_interface` 等字段。当前实现主要依赖：

- `cohesive_elements`
- `damage_states`
- `node_map`
- `interface_nodes`

体积单元与 CZM 单元的映射关系不存储在 `CohesiveMesh` 内，而是由 `jellyroll_collector_seed_mesh` 返回的 `czm_element_map` 和 `is_inner_layer` 保存在 `case.multi_spme_layout` 中。

---

## 四、网格与界面映射

### 4.1 `jellyroll_collector_seed_mesh` 当前返回内容

当前 `src/Jellyrollmodel.jl` 中，`jellyroll_collector_seed_mesh(param_dim; nθ, gsorder, phase, tol)` 会同时返回：

```julia
(
    thermal2D = mesh_unmerged,
    thermal2D_merged = thermal2D_merged,
    Jellyroll_czm = mesh_unmerged,
    merge_map = merge_map,
    interface_pairs = interface_pairs,
    czm_element_map = czm_element_map,
    element_layer = element_layer,
    is_inner_layer = is_inner_layer,
    inner_nodes = inner_nodes,
    outer_nodes = outer_nodes,
    pos_tab_nodes = pos_tab_nodes,
    neg_tab_nodes = neg_tab_nodes,
    ne = ne,
    nnode = nnode
)
```

与原设计稿不同：

- 没有 `create_czm` 关键字参数。
- 没有把 `czm_data` 直接挂在 mesh 对象上。
- CZM 相关映射信息通过返回值元组显式传出。

### 4.2 `setup_thermal2D_mesh` 的当前行为

当前正式接口名为：

```julia
setup_thermal2D_mesh(case, mesh_data; use_merged=nothing)
```

注意：

- 该函数返回新的 `case`，不是原地修改版本，因此没有 `setup_thermal2D_mesh!`。
- 当 `use_merged === nothing` 时，默认策略是：
  - `czm_enabled=false` 时使用 `thermal2D_merged`
  - `czm_enabled=true` 时使用 `thermal2D`
- 同时会把 `interface_pairs`、`czm_element_map`、`is_inner_layer` 等信息写入 `case.multi_spme_layout`。

### 4.3 `create_czm_mesh` 的当前方式

当前仍使用：

```julia
create_czm_mesh(thermal_mesh, param_dim; tol=1e-8)
```

实现方法是：

1. 在未合并热网格中寻找坐标重合的内外节点对。
2. 按极角排序后，将相邻节点对连接成界面单元。
3. 生成 `CohesiveElement` 数组与对应的 `DamageState`。

代码中不存在 `create_czm_mesh_from_thermal` 这个替代接口。

---

## 五、当前间隙导热模型

### 5.1 代码中的实际公式

当前 `compute_gap_conductance(D, δ_n, cohesive)` 定义在 `src/Materialmatrix.jl`，实现为：

```julia
delta = max(δ_n, 0.0)

D_sep = if delta <= delta0
    0.0
elseif delta < delta_c
    (delta - delta0) / (delta_c - delta0)
else
    1.0
end

D_clamped = clamp(max(D, D_sep), 0.0, 0.9999)
two_beta_lambda = 2.0 * beta * lambda_m

denom = if delta < delta0
    h_c0 + k_air / (delta + two_beta_lambda)
elseif delta < threshold
    h_c0 * (1.0 - D_clamped) + k_air / (delta + two_beta_lambda)
else
    h_c0 * (1.0 - D_clamped) + k_air / (delta + delta0)
end

h_eff = 1.0 / max(denom, 1e-30)
```

这意味着当前实现不是设计稿中最初给出的简化串联热阻公式，而是：

- 同时考虑损伤变量 `D` 与法向分离 `δ_n`。
- 使用 `D_sep` 对分离位移驱动的退化做补充。
- 引入 `lambda_m`、`beta`、`threshold` 来区分不同间隙尺度区间。

### 5.2 已实现的辅助接口

当前已经导出以下接口：

```julia
compute_gap_conductance(D, δ_n, cohesive)
compute_element_gap_conductance(czm_mesh, elem_idx, cohesive)
compute_all_gap_conductances(czm_mesh, cohesive)
```

---

## 六、热方程中的耦合方式

### 6.1 当前耦合入口

代码没有 `ThermalDistributed2D_with_CZM` 这个独立入口。当前做法是直接在：

```julia
ThermalDistributed2D_BC(KT, FT, case, t)
```

中处理 CZM 界面热耦合。

### 6.2 当前装配流程

当 `case.opt.czm_enabled == true` 时：

1. 从 `case.multi_spme_layout["czm_mesh"]` 读取 CZM 网格。
2. 若不存在，则尝试基于当前 `thermal2D` 网格即时创建。
3. 遍历每个 `cohesive_element`。
4. 读取 `DamageState.D` 与 `DamageState.δ_max_n`。
5. 用 `compute_gap_conductance` 计算 `h_eff`。
6. 按 `coeff = h_eff * czm_elem.length / (k_th * L_th)` 组装节点对耦合项。

当前写入矩阵的形式为：

```julia
K[nb, nb] -= coeff
K[nb, nt] += coeff
K[nt, nb] += coeff
K[nt, nt] -= coeff
```

说明：

- 这部分逻辑已经接入主热求解流程。
- 代码中没有独立的 `_add_interface_thermal_coupling!` 辅助函数。
- 代码中也没有早先文档里的 `_apply_interlayer_thermal_conductance!` 实现。

---

## 七、电流重分配与失效单元退出

### 7.1 已实现部分

当前 `solve_branch_currents_newton(...; deactivated_elements=nothing)` 已经支持传入失效单元列表：

- 失效单元在 `active_mask` 中被标记为非活跃。
- 非活跃单元电流被置零。
- 总电流约束重新只对活跃单元归一化。
- 结果中会记录：
  - `thermal2D n_active_elements`
  - `thermal2D n_cutoff_elements`
  - `thermal2D n_deactivated_elements`

### 7.2 尚未闭环部分

主流程中虽然会尝试从 `variables["deactivated_elements"]` 读取退出单元列表，但当前代码库中没有看到自动写入该字段的路径。因此当前状态更准确地说是：

- 分流求解器已经具备“失效单元退出”的能力。
- 但 CZM 断裂结果尚未稳定传递到该输入通道。

### 7.3 热源屏蔽的当前状态

`src/ThermalDistributed.jl` 中存在 `heatQ_Source_with_czm(...)`，会基于 `get_active_elements(czm_mesh, mesh_data)` 将非活跃单元热源置零；但当前主流程默认仍调用 `heatQ_Source(...)`，因此这部分属于“已实现辅助函数、尚未接入主线”。

---

## 八、循环求解与 SOH

### 8.1 当前 `CyclingResult` 字段

当前循环结果结构包含：

```julia
cycle_idx
capacity_charge
capacity_discharge
coulombic_efficiency
D_max
D_mean
n_fractured
T_max
soh
cycle_results
final_czm_mesh
initial_capacity
soh_terminated
```

与旧设计稿不同：

- 当前没有 `termination_reason` 字段。
- 当前没有 `n_active_elements` 历史字段。

### 8.2 当前 SOH 逻辑

SOH 由 `_update_soh_and_capacity!` 计算：

```julia
current_soh = cycle_result.capacity_discharge / initial_capacity
```

其中 `initial_capacity` 在第一个有效放电循环后确定。

### 8.3 当前终止逻辑

`_check_cycle_termination(...)` 当前实现为：

- 当 `current_soh <= case.opt.czm_soh_threshold` 且循环数大于 1 时终止。
- 当 `cycle_result.n_fractured > 0.5 * czm_mesh.n_cohesive` 时提前终止。

因此，当前代码并不是“所有内聚力单元断裂才终止”，而是“断裂比例超过 50% 就提前停”。

### 8.4 当前损伤更新频率

当前 `solve_phase(...)` 在阶段求解完成后，会尝试调用 `_update_czm_damage!` 更新损伤。因此当前实现应描述为：

- 充/放/静置每个阶段结束时进行一次 CZM 更新。
- 不是设计稿里写的“每个时间步更新一次”。

---

## 九、当前可直接使用的接口

### 9.1 网格与几何

```julia
jellyroll_collector_seed_mesh(param_dim; nθ=360, gsorder=2, phase=0.0, tol=1e-8)
setup_thermal2D_mesh(case, mesh_data; use_merged=nothing)
jellyroll_element_properties(mesh, param_dim)
jellyroll_tab_node_indices(mesh, param_dim)
edge_boundary(mesh, nidx, param_dim; which=:inner/:outer, theta_range=nothing, tol=1e-4)
```

### 9.2 CZM 与导热

```julia
create_czm_mesh(thermal_mesh, param_dim; tol=1e-8)
compute_gap_conductance(D, δ_n, cohesive)
compute_element_gap_conductance(czm_mesh, elem_idx, cohesive)
compute_all_gap_conductances(czm_mesh, cohesive)
get_fractured_elements(czm_mesh)
get_active_elements(czm_mesh, mesh_data)
```

### 9.3 循环求解

```julia
solve_phase(case, phase_type, t_max, I_current, V_limit, initial_state; ...)
solve_cycling(case, cycle_opt, czm_mesh; verbose=true, save_detailed=false)
```

---

## 十、建议的文档口径

为避免后续文档与实现再次脱节，建议统一采用以下表述：

1. 将本文视为“实现现状说明”，而不是“未执行的开发计划”。
2. 对于已经存在但未接线的能力，明确标注为“辅助函数已实现，主流程未接入”。
3. 对于未来计划新增的接口，不再在文档中写成已经存在的函数名。

目前最接近真实代码状态的总结是：

- 界面热阻耦合已经实现并接入二维热边界条件装配。
- CZM 损伤演化与循环 SOH 监控已经实现。
- 断裂导致的电化学退出和热源屏蔽处于部分实现状态，仍需要把 CZM 结果稳定写入主流程变量。