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

### 1.3 界面类型、真实面与离散单元

三种计数层级严格分开：

| 层级 | 数量/定义 | 含义 |
|---|---|---|
| 本构/材料类型 | `N_type_coh = 2` | `:PE_PCC`、`:NE_NCC` |
| 真实 cohesive 面 | `N_face_coh_per_repeat = 4` | PE(A)–PCC、PCC–PE(B)、NE(A)–NCC、NCC–NE(B) |
| 离散单元 | `N_elem_coh = 4 * N_seg` | `N_seg = length(theta) - 1`，按整条螺旋周向分段计数 |

`Φ` / `phi_pairs` 只表示跨匝 outer/inner 节点配对，不表示 cohesive 面或单元数。SP–涂层接触面属于独立接触模型，也不计入上述 cohesive 数量。

### 1.4 模块文件分工（b4c0cde 重构后）

| 文件 | 职责 |
|------|------|
| `src/czm.jl` | 本构装配：`moduli_of`、`eigenstrain_of`、`assemble_czm_system`、`assemble_bulk_stiffness`、`assemble_thermal_chemical_load`、`assemble_coupled_system(_full)`、网格惰性缓存访问器（`bulk_stiffness`/`cohesive_geometry`/`assembly_workspace`）、`collector_params` |
| `src/CzmMesh.jl` | 网格拓扑：`CohesiveElement`、`create_czm_mesh`（共边识别界面、节点复制、外层 bulk 连接重写、`cohesive_to_thermal` 映射） |
| `src/CzmBC.jl` | 边界条件：`apply_bc_czm`（罚函数 Dirichlet）、`identify_bc_nodes_czm`（内外圈 `:fixed_xy`） |
| `src/CzmSolve.jl` | 非线性求解：`CZMResult`、`solve_czm_basic_step`、`newton_raphson_czm`、`solve_czm_arc_length_step`、统一入口 `solve_czm_step`、`update_damage_per_interface`、`backtrack_line_search!` |
| `src/CzmPostProcess.jl` | 统计与损伤管理：`get_damage_statistics`、`check_fracture_criterion`、`reset_damage_states!`、`accumulate_cycle_damage!`、`czm_output_to_variables`（2026-08-30 起接收 `damage_states` 向量/MechState） |
| `src/CzmUnitMesh.jl` | 验证用 8 层平直条带网格 `create_unit_czm_strip`（4 个 COH2D4 + 硬断言） |
| `src/CouplingState.jl` | 状态与耦合：`MechState`（演化状态聚合）、`DamageState`、`update_czm_damage!`、`compute_czm_strain_inputs` |
| `src/Materialmatrix.jl` | 本构数学：`bilinear_traction_state`、`bilinear_traction`、`bilinear_tangent`、`update_damage`、间隙导热系列 |

结构体宿主：`CzmSubmesh` / `CohesiveMesh` / 抽象类型定义在 `src/SetMesh.jl:25-76`（避免 include 顺序问题）。

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

参数按**界面类型**独立给出（2026-08-30 重构：挂 `CurrentCollector`——`param.PCC` 承载
PE-PCC 界面、`param.NCC` 承载 NE-NCC 界面；`Cohesive` 结构已删除）：每套
`σ_max / K_n / δ_0 / G_c / δ_c`（法向）与 `τ_max / K_t / δ_0_t / G_c_t / δ_c_t`（切向）
加 `eta`（BK 指数）与界面热阻五参数（`h_c0/k_air/lambda_m/beta/threshold`）。

| 参数 | 符号 | 说明 |
|------|------|------|
| 界面刚度 | K_n, K_t | 法向/切向刚度（每界面独立） |
| 荷载强度 | σ_max / τ_max | 最大法向/切向牵引应力（每界面独立） |
| 初始分离 | δ₀ = T_max/K | 开始损伤的分离位移 |
| 临界分离 | δ_c = 2G_c/T_max | 完全断裂的分离位移 |
| 断裂能 | G_c | 界面断裂所需能量（每界面独立） |
| 极片弹性模量 | E_coat | PE/NE 极片（涂层）宏观弹性模量，逐层体刚度（`moduli_of`）与宏观应力 |

### 2.3 BK 混合模式准则

**混合模式临界分离**：

$$
\delta_c^{mix} = \left[(G_n/G_c^I)^\alpha + (G_t/G_c^{II})^\alpha\right]^{1/\alpha} \cdot \delta_c^I
$$

### 2.4 bilinear_traction_state 函数

**功能**：计算双线性本构的牵引力、损伤变量和更新后的损伤状态

```julia
function bilinear_traction_state(δ_n::Float64, δ_t::Float64,
    damage_state::DamageState, ip::CurrentCollector, czm_model::String; visc_beta::Float64=1.0)
```

（定义于 `src/Materialmatrix.jl`；第 4 参 `ip` 是**界面参数宿主**——`param.PCC` 或
`param.NCC`（2026-08-30 重构，`collector_params(param, iface)` 分派）；第 5 参
`czm_model`（"model1"/"mix"）来源 `opt.czm.model`，显式传入。）

**返回**：
- T_n, T_t：法向/切向牵引力
- D_eq：等效损伤变量
- new_state::DamageState：更新后的损伤状态（含 D_visc、δ_max_eff、accumulated_damage、fractured）

**粘性正则化**：牵引力使用 $D_{visc}$ 而非 $D_{eq}$：$T_n = (1 - D_{visc}) K_n \delta_n$

---

## 3. 有限元实现

### 3.1 单元数据结构

**CohesiveElement**（`src/CzmMesh.jl:3-12`）：
- `id`：单元编号
- `nodes`：节点编号 (4节点界面单元)
- `nodes_bottom` / `nodes_top`：底面/顶面节点（顺序一致，用于装配节点对）
- `length`：界面长度
- `interface_type`：`:PE_PCC` 或 `:NE_NCC`，表示参数类型而非真实面编号
- `host_outer_elem` / `host_inner_elem`：真实面两侧的 Q4 单元编号

**DamageState**（`src/czm.jl:17`）：
- `D`：损伤变量 (0-1)
- `δ_max_n`：历史最大法向分离
- `δ_max_t`：历史最大切向分离
- `D_visc`：粘性正则化损伤变量
- `δ_max_eff`：等效最大分离位移
- `accumulated_damage`：累积损伤
- `fractured`：是否已完全断裂

### 3.2 CZMResult 结构体

**定义**（`src/CzmSolve.jl:1-15`）：

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
K_{coh} = \int_L B^T D_{tan} B \, dL
$$

其中 \(D_{tan}\) 为本构切线模量。

**内力向量**：

$$
f_{int} = \int_L B^T T \, dL
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

**装配缓存（2026-08-30 重构）**：缓存随网格生灭——`CohesiveMesh` 惰性字段
`K_bulk`/`cohesive_geom`/`ws`，经访问器 `bulk_stiffness(czm_mesh, param)` /
`cohesive_geometry(czm_mesh)` / `assembly_workspace(czm_mesh)` 首次调用构建、
跨步只读复用；**对象身份即失效判据**（重建网格即丢缓存）。新鲜度由参数冻结契约保证
（SetCase 归一化后不得改 `param`）。BC 不缓存（每次求解入口现算）。
`CZMAssemblyWorkspace` 提供 `mul!` 复用的工作区；`CohesiveElementGeom` 缓存单元几何量。

---

## 4. 求解方法

### 4.1 求解框架

所有迭代方法通过统一入口 `solve_czm_step` 分发（2026-08-30 终态签名）：

```julia
result = solve_czm_step(
    czm_mesh, ms::MechState, param::Params, F_ext, czm_opt::CzmOptions;
    dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,
    eigenstrain=nothing, mech_state=nothing, prestress=nothing)
# -> CZMResult（收敛后在 ms 上原位提交损伤/位移；失败不触碰 ms）
```

求解配置（iter_method/max_iter/tol/load_steps/arc_length_alpha/viscous_*/model/
fix_inner/geo_nonlinear/j2_plasticity）全部从 `czm_opt` 展开，不再逐 kwarg 传递。

热-化学应变系数不再逐参量传入：ε₀ 按单元 `material_type` 经
`eigenstrain_of(param, mt)` 分层计算（α=该层 `alphaT`，β=该层 `Ω/3`，
2026-08-29 α/β 同批分层化）。

**参数直读（2026-08-30 重构）**：求解/装配函数统一接收 `param::Params`，热路径按
单元 `interface_type` 分派 `collector_params(param, iface)`（`:PE_PCC`→`param.PCC`、
`:NE_NCC`→`param.NCC`）；`CzmParamCache`/内容哈希失效机制已删除，失效语义唯一化为
对象身份。

**默认方法**：函数签名默认 `iter_method="load_substep"`；`Option.czm_iter_method`
默认 `"basic"`（`src/Option.jl:73`），主循环调用时显式传入。

**残差格式**：

$$
R = F_{applied} - f_{int}(u)
$$

其中 $F_{applied}$ 为外力向量（含热-化学载荷），$f_{int}$ 为内力向量。

### 4.2 损伤更新策略（隐式冻结）

**所有三种方法采用相同的损伤更新策略**：在一个时间步内，NR 迭代过程中损伤状态保持冻结（使用时间步开始时的 $D$），仅当位移场收敛后才统一更新损伤。

```
进入时间步 → damage_states = clone(ms.damage_states)（试探态，ms 不受影响）
  NR 迭代 (全过程用 D_begin 计算刚度 K 和内力 f_int)
  → 位移 u 收敛
  → ms.damage_states = update_damage_per_interface(D_begin, separations(u_converged), param, czm_model)
退出时间步 → 收敛则 ms.damage_states/ms.u_prev 原位提交；失败则丢弃试探态
```

**设计原因**：CZM 软化段的临界分离位移很小，即使微量载荷即可导致损伤从 0 跳至近 1.0。若在 NR 迭代或载荷子步中途更新损伤，刚度矩阵会骤降，导致后续迭代/子步在近乎零刚度的系统上发散。冻结损伤保证刚度矩阵在整个求解过程中一致，与隐式时间积分的物理假设（一个时间步内损伤为常数）相符。

### 4.3 迭代方法

#### 4.3.1 basic — 基本 Newton-Raphson

```julia
function solve_czm_basic_step(czm_mesh, F_ext, param, ms::MechState; ...)
```
（`src/CzmSolve.jl`）

- 从 `ms.u_prev`（上一步收敛位移）出发
- 一次性施加全量载荷 $F_{ext} + F_{thermo-chem}$
- NR 迭代 + 回溯线搜索（最多 8 次减半）
- 收敛后更新损伤

**适用场景**：时间步内载荷增量较小（标准时间步进），是 `Option` 层面的默认方法且最稳健。

#### 4.3.2 load_substep — 载荷子步法

```julia
function newton_raphson_czm(czm_mesh, F_ext, param, ms::MechState; ...)
```
（`src/CzmSolve.jl`）

将载荷从当前平衡态逐步推进到目标态：

$$
F_{applied} = f_{int}^{ref} + \alpha \cdot F_{\Delta}
$$

其中：
- $f_{int}^{ref} = f_{int}(ms.u_{prev})$：上一步平衡态的内力
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
function solve_czm_arc_length_step(czm_mesh, F_ext, param, ms::MechState; ...)
```
（`src/CzmSolve.jl`）

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
    param::Params,
    dT_elem::Vector{Float64},
    Δsoc_n_elem::Vector{Float64},
    Δsoc_p_elem::Vector{Float64}
)
```
（`src/czm.jl`；α/β 分层分辨率——由 `eigenstrain_of(param, material_type[e])` 按层取
`alphaT` 与 `Ω/3`，不再跨层统一传入）

**初始应变公式**（逐层，mt 为该单元材料类型）：

$$
\varepsilon_0^{(mt)} = \alpha_T^{(mt)} \cdot \Delta T + \tfrac{\Omega^{(mt)}}{3} \cdot \Delta soc^{(mt)}
$$

电极膨胀只作用于本层涂层（NE→Δsoc_n、PE→Δsoc_p），集流体/隔膜只有热应变
（Jellyroll 参数集 SP/PCC/NCC.alphaT 显式置零）。

**等效节点力**：

$$
F_{thermo-chem} = \int B^T D \varepsilon_0 \, d\Omega
$$

### 4.6 完整耦合系统组装

```julia
function assemble_coupled_system(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param::Params;
    F_ext=nothing,
    F_thermo_chem=nothing,
    damage_states=nothing,
    K_bulk_cached=nothing,
    geom_cache=nothing,
    ws=nothing,
    visc_beta=1.0
)
```
（`src/czm.jl:546`；一步完成"组装 + 残差"的变体是 `assemble_coupled_system_full`，`:586`）

**返回**：
- `K_total`：总刚度矩阵
- `f_int_total`：内力向量
- `separations`：分离位移
- `tractions`：牵引力

### 4.7 模量层级（极片 vs 颗粒）

CZM 应变驱动所需的体模量必须使用**极片（涂层）宏观弹性模量**（经 `moduli_of` 逐层取
`E_coat`/连续层 `E`），而**非**颗粒层面的 `PE.E`/`NE.E`。两者物理尺度相差约两个数量级
（典型值：颗粒 ~1e10 Pa，极片 ~5e8 Pa），错误使用会造成刚度矩阵系统性偏差。

2026-08-30 重构后 E_eff 概念与按界面参数缓存（含派生量 Λ、E_star、L_ch）已整体删除：
装配直读 `param`，模量经 `moduli_of` 按材料类型逐层给出；Λ 使用点内联
`scale.L/δ_czm`。全叠合厚度加权仅保留在参考尺度 `scale.E_coat` 的定义中。

| 按界面损伤更新 | src/CzmSolve.jl | `update_damage_per_interface` (52) |
| 损伤克隆 | src/CzmSolve.jl | `clone_damage_states`（克隆链 `clone_czm_mesh_with_damage` 已删除） |
| BC 处理 | src/CzmSolve.jl | `extract_bc_dofs` (79), `apply_czm_dirichlet!` (130) |
| 回溯线搜索 | src/CzmSolve.jl | `backtrack_line_search!` (107) |
| 弧长增广矩阵 | src/CzmSolve.jl | `build_arc_length_augmented_matrix` (156) |
| CZM 调度入口 | src/CouplingState.jl | `update_czm_damage!` |
| 热-化学应变输入 | src/CouplingState.jl | `compute_czm_strain_inputs` |
| 界面参数分派 | src/czm.jl | `collector_params(param, iface)`（`:PE_PCC`→PCC、`:NE_NCC`→NCC） |
| 本构模型 | src/Materialmatrix.jl | `bilinear_traction_state` (68), `bilinear_tangent` (182), `update_damage` (293) |
| 内聚力单元装配 | src/czm.jl | `assemble_czm_system` (79) |
| 体刚度装配 | src/czm.jl | `assemble_bulk_stiffness` (268) |
| 系统组装 | src/czm.jl | `assemble_coupled_system` (546), `assemble_coupled_system_full` (586) |
| 热-化学载荷 | src/czm.jl | `assemble_thermal_chemical_load` (347) |
| 网格惰性缓存 | src/czm.jl | `bulk_stiffness` / `cohesive_geometry` / `assembly_workspace` |
| CZM 网格创建 | src/CzmMesh.jl | `create_czm_mesh` (33) |
| 验证条带网格 | src/CzmUnitMesh.jl | `create_unit_czm_strip` (7) |
| CZM 边界条件 | src/CzmBC.jl | `apply_bc_czm` (7), `identify_bc_nodes_czm` (105) |
| 间隙导热 | src/Materialmatrix.jl | `compute_gap_conductance` (329), `compute_element_gap_conductance` (359) |
| 有效面积因子 | src/Materialmatrix.jl | `effective_area_factor` (424) |
| 损伤统计 | src/CzmPostProcess.jl | `get_damage_statistics` (12), `check_fracture_criterion` (37) |
| CZM 后处理 | src/CzmPostProcess.jl | `czm_output_to_variables` (99) |
| SOH 计算 | src/CyclePostProcess.jl | `update_soh_and_capacity!` (154) |
| 循环求解 | src/CycleSolver.jl | `solve_cycling` (209) |
