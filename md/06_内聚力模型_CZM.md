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

### 4.1 Newton-Raphson 迭代

**迭代格式**：

$$
F(u) = f_{int}(u) - F_{ext} - F_{thermo-chem}
$$

**更新公式**：

$$
u^{(k+1)} = u^{(k)} - J^{-1} \cdot F
$$

**函数签名**：

```julia
function newton_raphson_czm(
    czm_mesh::CohesiveMesh, 
    F_ext::Vector{Float64}, 
    E_eff::Float64, 
    ν_eff::Float64, 
    cohesive_params::Cohesive, 
    param_dim;
    α_eff::Float64=0.0,      # 热膨胀系数
    β_n::Float64=0.0,        # 负极化学膨胀系数
    β_p::Float64=0.0,        # 正极化学膨胀系数
    dT_elem=nothing,         # 单元温度变化
    Δsoc_n_elem=nothing,     # 负极 SOC 变化
    Δsoc_p_elem=nothing,     # 正极 SOC 变化
    max_iter::Int=50, 
    tol::Float64=1e-8, 
    u0=nothing, 
    n_load_steps::Int=10
)
```

**返回**：
- `result::CZMResult`：求解结果
- `new_czm_mesh::CohesiveMesh`：更新损伤状态后的网格

### 4.2 热-化学载荷

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

### 4.3 完整耦合系统组装

```julia
function assemble_coupled_system_full(
    czm_mesh::CohesiveMesh, 
    u::Vector{Float64},
    E_eff::Float64, 
    ν_eff::Float64,
    α_eff::Float64, 
    β_n::Float64, 
    β_p::Float64,
    cohesive_params::Cohesive,
    dT_elem::Vector{Float64},
    Δsoc_n_elem::Vector{Float64},
    Δsoc_p_elem::Vector{Float64};
    F_ext=nothing,
    damage_states=nothing
)
```

**返回**：
- `K_total`：总刚度矩阵
- `R`：残差向量
- `separations`：分离位移
- `tractions`：牵引力

### 4.4 迭代方法选项

| 方法 | 说明 | 适用场景 |
|------|------|----------|
| `"basic"` | 基本 Newton-Raphson | 简单问题，单步载荷 |
| `"load_substep"` | 载荷子步法（默认） | 一般问题，提高收敛性 |
| `"arc_length"` | 弧长法 | snap-through/snap-back 问题 |

**载荷子步法**：
- 将总载荷分成多个子步
- 每个子步内使用 Newton 迭代
- 提高收敛性但增加计算量

### 4.5 收敛性改进

**线搜索**：当残差下降缓慢时启用

```julia
Δu_norm = norm(Δu)
max_Δu = 1e-6
if Δu_norm > max_Δu
    Δu = Δu * (max_Δu / Δu_norm)
end
```

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

| 功能 | 文件 | 函数/行号 |
|------|------|-----------|
| CZMResult 结构体 | src/CzmSolve.jl | `CZMResult` (1-15) |
| Newton-Raphson 求解 | src/CzmSolve.jl | `newton_raphson_czm` (30-120) |
| CZM 单步求解 | src/CzmSolve.jl | `solve_czm_step` |
| 后处理 | src/CzmSolve.jl | `czm_output_to_variables` (496-514) |
| 本构模型 | src/Materialmatrix.jl | `bilinear_traction_state` |
| 单元矩阵 | src/czm.jl | `cohesive_element_matrices` |
| 系统组装 | src/czm.jl | `assemble_czm_system`, `assemble_coupled_system` |
| 完整耦合组装 | src/czm.jl | `assemble_coupled_system_full` (427-443) |
| 热-化学载荷 | src/czm.jl | `assemble_thermal_chemical_load` (349-403) |
| CZM 网格创建 | src/czm.jl | `create_czm_mesh` |
| 间隙导热 | src/Materialmatrix.jl | `compute_gap_conductance`, `compute_element_gap_conductance` |
| 损伤统计 | src/czm.jl | `get_damage_statistics`, `check_fracture_criterion` |
| 损伤更新 | src/CycleSolver.jl | `_update_czm_damage!` |
| SOH 计算 | src/CycleSolver.jl | `_update_soh_and_capacity!` |
| 循环求解 | src/CycleSolver.jl | `solve_cycling` |
