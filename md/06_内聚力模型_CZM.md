# 内聚力模型 (CZM)

本文档描述果冻卷电池层间界面的内聚力模型 (Cohesive Zone Model) 实现，用于模拟界面损伤的演化。

---

## 1. 模型概述

### 1.1 物理背景

果冻卷电池在循环过程中，层间界面会受到应力的作用而产生损伤。内聚力模型用于描述这种界面损伤的起始和演化过程。

**主要功能**：
- 模拟层间界面的力学行为
- 描述损伤的起始和扩展
- 预测电池的循环寿命

### 1.2 主要假设

- 界面破坏发生在预设的界面单元上
- 损伤是渐进的，从完好到完全断裂
- 损伤状态影响界面的力学和热学性能

---

## 2. 本构模型

### 2.1 双线性牵引-分离律

**弹性阶段** (δ < δ₀)：

$$
T = K \cdot \delta
$$

其中 K 为界面刚度，δ 为分离位移。

**软化阶段** (δ₀ ≤ δ < δ_c)：

$$
T = (1 - D) \cdot K \cdot \delta
$$

**损伤变量定义**：

$$
D = \frac{\delta_c(\delta - \delta_0)}{\delta(\delta_c - \delta_0)}
$$

**完全断裂** (δ ≥ δ_c)：

$$
T = 0, \quad D = 1
$$

### 2.2 参数定义

| 参数 | 符号 | 说明 |
|------|------|------|
| 界面刚度 | K_n, K_t | 法向/切向刚度 |
| 荷载强度 | T_max | 最大牵引应力 |
| 初始分离 | δ₀ = T_max/K | 开始损伤的分离位移 |
| 临界分离 | δ_c = 2G_c/T_max | 完全断裂的分离位移 |
| 断裂能 | G_c | 界面断裂所需能量 |

### 2.3 BK 混合模式准则

**混合模式临界分离**：

$$
\delta_c^{mix} = \left[(G_n/G_c^I)^\alpha + (G_t/G_c^{II})^\alpha\right]^{1/\alpha} \cdot \delta_c^I
$$

### 2.4 bilinear_traction_state 函数

**功能**：计算双线性本构的牵引力和切线模量

```julia
function bilinear_traction_state(δ_n, δ_t, δ0_n, δc_n, δ0_t, δc_t, K_n, K_t, D_prev, mixed_mode)
```

**返回**：
- T_n, T_t：法向/切向牵引力
- D_tan_n, D_tan_t：法向/切向切线模量
- D：损伤变量

---

## 3. 有限元实现

### 3.1 单元数据结构

**CohesiveElement**：
- `nodes`：节点编号 (4节点界面单元)
- `length`：界面长度
- `normal`：法向量

**DamageState**：
- `D`：损伤变量 (0-1)
- `δ_max_n`：历史最大法向分离
- `δ_max_t`：历史最大切向分离
- `fractured`：是否已完全断裂

### 3.2 CZMResult 结构体

**定义**：

```julia
mutable struct CZMResult
    displacement::Vector{Float64}    # 位移场 (ndof)
    damage::Vector{Float64}          # 损伤变量 (n_coh)
    traction_n::Vector{Float64}      # 法向牵引力 (n_coh)
    traction_t::Vector{Float64}      # 切向牵引力 (n_coh)
    separation_n::Vector{Float64}    # 法向分离 (n_coh)
    separation_t::Vector{Float64}    # 切向分离 (n_coh)
    converged::Bool                  # 是否收敛
    iterations::Int64                # 迭代次数
    residual_norm::Float64           # 最终残差范数
end
```

**构造函数**：

```julia
CZMResult(ndof::Int, n_coh::Int) = new(
    zeros(ndof), zeros(n_coh), zeros(n_coh), zeros(n_coh),
    zeros(n_coh), zeros(n_coh), false, 0, Inf)
```

### 3.3 单元刚度矩阵

**分离位移计算**：

$$
\delta = B \cdot u_{element}
$$

其中 B 为位移-分离关系矩阵。

**切线刚度矩阵**：

$$
K_{coh} = \int B^T D_{tan} B \, dA
$$

其中 \(D_{tan}\) 为本构切线模量。

**内力向量**：

$$
f_{int} = \int B^T T \, dA
$$

### 3.4 系统组装

**全局刚度矩阵**：

$$
K_{global} = \sum_e K_{coh,e}
$$

**全局内力向量**：

$$
f_{global} = \sum_e f_{int,e}
$$

---

## 4. 求解方法

### 4.1 求解框架

所有迭代方法通过统一入口 `solve_czm_step` 分发：

```julia
result, updated_czm_mesh = solve_czm_step(
    czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev;
    iter_method="basic",   # "basic" | "load_substep" | "arc_length"
    max_iter=50, tol=1e-8, n_load_steps=10, visc_beta=1.0, ...
)
```

**残差格式**：

$$
R = F_{applied} - f_{int}(u)
$$

其中 $F_{applied}$ 为外力向量（含热-化学载荷），$f_{int}$ 为内力向量。

### 4.2 损伤更新策略（隐式冻结）

**所有三种方法采用相同的损伤更新策略**：在一个时间步内，NR 迭代过程中损伤状态保持冻结（使用时间步开始时的 $D$），仅当位移场收敛后才统一更新损伤。

```
进入时间步 → damage_states = D_begin
  NR 迭代 (全过程用 D_begin 计算刚度 K 和内力 f_int)
  → 位移 u 收敛
  → damage_states = update_damage(D_begin, separations(u_converged))
退出时间步 → D_end 用于下一时间步
```

**设计原因**：CZM 软化段的临界分离位移很小，即使微量载荷即可导致损伤从 0 跳至近 1.0。若在 NR 迭代或载荷子步中途更新损伤，刚度矩阵会骤降，导致后续迭代/子步在近乎零刚度的系统上发散。冻结损伤保证刚度矩阵在整个求解过程中一致，与隐式时间积分的物理假设（一个时间步内损伤为常数）相符。

### 4.3 迭代方法

#### 4.3.1 basic — 基本 Newton-Raphson

```julia
function solve_czm_basic_step(czm_mesh, F_ext, ..., u_prev; ...)
```

- 从 `u_prev`（上一步收敛位移）出发
- 一次性施加全量载荷 $F_{ext} + F_{thermo-chem}$
- NR 迭代 + 回溯线搜索（最多 8 次减半）
- 收敛后更新损伤

**适用场景**：时间步内载荷增量较小（标准时间步进），是默认且最稳健的方法。

#### 4.3.2 load_substep — 载荷子步法

```julia
function newton_raphson_czm(czm_mesh, F_ext, ..., param; n_load_steps=10, u0=u_prev, ...)
```

将载荷从当前平衡态逐步推进到目标态：

$$
F_{applied} = f_{int}^{ref} + \alpha \cdot F_{\Delta}
$$

其中：
- $f_{int}^{ref} = f_{int}(u_{prev})$：上一步平衡态的内力
- $F_{\Delta} = F_{target} - f_{int}^{ref}$：载荷增量（通常很小）
- $\alpha \in [0, 1]$：载荷进度，从 0 到 1 逐步推进

**自适应子步控制**：
- 初始步长：`step_size = 1 / n_load_steps`
- 子步收敛 → 步长 ×1.25（加速）
- 子步失败 → 步长 ×0.5（回退），直至 `step_size_min`
- 每个子步内：NR 迭代 + 线搜索（8 次减半）

**收敛判据**：
- 子步收敛：`R_norm < tol × 10`
- 整体收敛：`load_progress ≥ 1 - ε` 且 `R_norm < tol × 100`

**适用场景**：需要更细粒度载荷控制的工况；对 basic 已能收敛的标准仿真不会带来额外精度提升，但代价更高（约 4-5 倍 CZM 耗时）。

#### 4.3.3 arc_length — Crisfield 弧长法

```julia
function solve_czm_arc_length_step(czm_mesh, F_ext, ..., u_prev; n_load_steps=10, arc_length_alpha=1.0, ...)
```

在载荷子步的基础上引入 Crisfield 圆柱弧长约束：

$$
\| \Delta u \|^2 = l_{arc}^2
$$

**求解过程**（每个子步）：
1. 预测步：$\Delta u_{pred} = K_{bc}^{-1} F_{\Delta} \cdot \Delta\lambda_{pred}$，计算弧长目标 $l_{arc}$
2. 校正步（NR 迭代）：求解两个线性系统
   - $\Delta u_R = K_{bc}^{-1} R_{bc}$（残差校正）
   - $\Delta u_F = K_{bc}^{-1} F_{\Delta,bc}$（载荷方向）
3. 通过弧长约束的二次方程求解载荷参数增量 $\Delta\lambda$
4. 选择与预测方向最接近的根

**适用场景**：snap-through / snap-back 等强非线性问题。

### 4.4 三种方法性能对比

在 testexample.jl（1C 放电，nθ=16，336 单元，368 时间步）上的典型表现：

| 指标 | basic | load_substep | arc_length |
|------|-------|-------------|------------|
| CZM 耗时 | ~5.5 s (35%) | ~22 s (68%) | ~22 s (67%) |
| 总耗时 | ~25 s | ~43 s | ~41 s |
| 物理结果 | 基准 | 与 basic 一致 | 与 basic 一致 |

对于标准时间步进仿真，`basic` 是推荐的首选方法。

### 4.5 热-化学载荷

**功能**：计算由温度变化和 SOC 变化引起的等效节点力

```julia
function assemble_thermal_chemical_load(
    czm_mesh::CohesiveMesh,
    E_eff::Float64,
    ν_eff::Float64,
    α_eff::Float64,          # 热膨胀系数
    β_n::Float64,            # 负极化学膨胀系数
    β_p::Float64,            # 正极化学膨胀系数
    dT_elem::Vector{Float64},
    Δsoc_n_elem::Vector{Float64},
    Δsoc_p_elem::Vector{Float64}
)
```

**初始应变公式**：

$$
\varepsilon_0 = \alpha \cdot \Delta T + \beta_n \cdot \Delta soc_n + \beta_p \cdot \Delta soc_p
$$

**等效节点力**：

$$
F_{thermo-chem} = \int B^T D \varepsilon_0 \, d\Omega
$$

### 4.6 完整耦合系统组装

```julia
function assemble_coupled_system(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    E_eff::Float64,
    ν_eff::Float64,
    cohesive_params::Cohesive;
    damage_states=nothing,
    K_bulk_cached=nothing,
    geom_cache=nothing,
    ws=nothing,
    visc_beta=1.0
)
```

**返回**：
- `K_total`：总刚度矩阵
- `f_int_total`：内力向量
- `separations`：分离位移
- `tractions`：牵引力

---

## 5. 与热的耦合

### 5.1 界面导热系数模型

**间隙导热系数**：

$$
h_{eff} = \frac{1}{h_{c0}(1-D) + \frac{k_{air}}{\delta + 2\beta\lambda_m}}
$$

其中：
- h_c0：完好界面导热系数
- k_air：空气导热系数
- λ_m：微观粗糙度参数
- β：间隙尺度系数

### 5.2 损伤对导热的影响

**影响机制**：
1. 损伤增加 → 接触面积减少
2. 间隙增大 → 热阻增加
3. 导热系数降低

**有限元装配**：

```julia
K[nb,nb] -= h_eff * l / (k_ref * L_th)
K[nb,nt] += h_eff * l / (k_ref * L_th)
K[nt,nb] += h_eff * l / (k_ref * L_th)
K[nt,nt] -= h_eff * l / (k_ref * L_th)
```

其中 l 为界面单元长度，nb, nt 为界面两侧节点。

---

## 6. 循环求解中的 CZM

### 6.1 SOH 监控

**SOH 定义**：

$$
SOH = C_{discharge} / C_{initial}
$$

**终止条件**：
- SOH ≤ SOH_threshold（默认 0.8）
- 断裂单元比例 > 50%

### 6.2 损伤更新

**更新频率**：由 opt.czm_update_interval 控制
- 默认每个充/放电阶段结束更新一次
- 非每个时间步更新

**更新流程**：
1. 收集当前应力状态
2. 计算各界面单元的分离位移
3. 更新损伤变量
4. 记录最大损伤和断裂数量

### 6.3 失效单元处理

**电化学退出**：
- 断裂单元不再参与电化学反应
- 分流求解时排除断裂单元

**热源屏蔽**：
- 断裂单元热源置零
- 通过 `heatQ_Source_with_czm` 函数实现

---

## 7. 后处理

### 7.1 czm_output_to_variables 函数

**功能**：将 CZM 求解结果写入 variables 字典

```julia
function czm_output_to_variables(
    czm_mesh::CohesiveMesh, 
    result::CZMResult, 
    variables::Dict{String, Union{Array{Float64}, Float64}}
)
```

**输出变量**：

| 键 | 类型 | 说明 |
|----|------|------|
| `czm displacement x` | Vector | x 方向位移 |
| `czm displacement y` | Vector | y 方向位移 |
| `czm damage` | Vector | 各单元损伤变量 |
| `czm traction normal` | Vector | 法向牵引力 |
| `czm traction tangent` | Vector | 切向牵引力 |
| `czm separation normal` | Vector | 法向分离 |
| `czm separation tangent` | Vector | 切向分离 |
| `czm max damage` | Float64 | 最大损伤 |
| `czm mean damage` | Float64 | 平均损伤 |
| `czm fractured elements` | Float64 | 断裂单元数 |

---

## 8. 代码位置

| 功能 | 文件 | 函数 |
|------|------|------|
| CZMResult 结构体 | src/CzmSolve.jl | `CZMResult` (1-15) |
| basic NR 求解 | src/CzmSolve.jl | `solve_czm_basic_step` (154) |
| 载荷子步法求解 | src/CzmSolve.jl | `newton_raphson_czm` (486) |
| 弧长法求解 | src/CzmSolve.jl | `solve_czm_arc_length_step` (261) |
| 统一求解入口 | src/CzmSolve.jl | `solve_czm_step` (662) |
| 损伤克隆/还原 | src/CzmSolve.jl | `clone_damage_states`, `clone_czm_mesh_with_damage` |
| BC 处理 | src/CzmSolve.jl | `extract_bc_dofs`, `apply_czm_dirichlet!` |
| 回溯线搜索 | src/CzmSolve.jl | `backtrack_line_search!` (83) |
| 弧长增广矩阵 | src/CzmSolve.jl | `build_arc_length_augmented_matrix` (142) |
| CZM 调度入口 | src/CouplingState.jl | `update_czm_damage!` |
| 本构模型 | src/Materialmatrix.jl | `bilinear_traction_state`, `update_damage` |
| 单元矩阵 | src/czm.jl | `cohesive_element_matrices` |
| 系统组装 | src/czm.jl | `assemble_coupled_system` |
| 热-化学载荷 | src/czm.jl | `assemble_thermal_chemical_load` |
| CZM 网格创建 | src/czm.jl | `create_czm_mesh` |
| 间隙导热 | src/Materialmatrix.jl | `compute_gap_conductance`, `compute_element_gap_conductance` |
| 损伤统计 | src/czm.jl | `get_damage_statistics`, `check_fracture_criterion` |
| 损伤更新（循环） | src/CycleSolver.jl | `_update_czm_damage!` |
| SOH 计算 | src/CycleSolver.jl | `_update_soh_and_capacity!` |
| 循环求解 | src/CycleSolver.jl | `solve_cycling` |
