# SPMe + 二维分布式热 + CZM 耦合仿真：函数依赖分析

> 聚焦范围：Jellyroll 型号、SPMe 模型、distributed2D 热模型、CZM 内聚力模型
> 目的：为后续代码简化提供函数调用关系和依赖指导

---

## Part A: 执行流程总览

### A1. 初始化阶段

```
用户代码
│
├── ChooseCell("Jellyroll")                          [SetParams.jl:231]
│   └── include("parameters/Jellyroll.jl")           ← 加载 Jellyroll 物理参数
│   └── 计算统一参考尺度 (scale 字段)
│   → 输出: param_dim::Params (有量纲参数)
│
├── Option()                                         [Option.jl:35]
│   └── 设置: model="SPMe", thermal_enabled=true, thermalmodel="distributed2D"
│   └── 设置: per_element_spme=true, czm_enabled=true, mechanicalmodel="full"
│   → 输出: opt::Option
│
├── SetCase(param_dim, opt)                          [SetCase.jl:1]
│   └── NormaliseParam(param_dim)                    [SetParams.jl:308]
│   └── 分支 opt.model=="SPMe": 创建颗粒网格 + 电解液网格
│   └── 分支 opt.thermalmodel=="distributed2D": 设置温度DOF索引
│   → 输出: case::Case (含 param, opt, mesh, index)
│
├── jellyroll_collector_seed_mesh(param_dim; nθ, gsorder)  [Jellyrollmodel.jl:17]
│   └── 阿基米德螺旋线几何生成
│   └── 识别界面节点对、分层信息、极耳节点
│   └── GetGS()                                      ← 高斯积分
│   → 输出: JellyrollMesh 结构体
│
├── setup_thermal2D_mesh(case, mesh_data)             [Jellyrollmodel.jl]
│   └── 将 JellyrollMesh 写入 case.mesh["thermal2D"]
│   └── 创建 MultiSPMeLayout(case.layout)
│   └── 创建 MeshGeometry(case.geometry)
│   → 修改: case (添加 mesh, layout, geometry)
│
└── [可选] create_czm_mesh(case.mesh["thermal2D"], param_dim)  [czm.jl:56]
    └── 坐标重合检测界面节点对
    └── 创建 CohesiveElement 列表和 DamageState 列表
    → 输出: czm_mesh::CohesiveMesh
    → 挂载: case.czm_mesh = czm_mesh
```

**关键开关影响：**
- `opt.model == "SPMe"` → 创建电解液网格和索引 (SetCase.jl:26-43)
- `opt.thermalmodel == "distributed2D"` → 设置温度DOF索引 (SetCase.jl:82-85)
- `opt.czm_enabled == true` → 需要在初始化后创建 czm_mesh

---

### A2. 求解阶段 — 主循环

```
Solve(case; initial_state, return_final_state, ...)   [Solve.jl:1]
│
├── [初始化]
│   ├── 分支 opt.model=="thermal": 纯热求解路径 (跳过)
│   ├── multi_spme_enabled = opt.per_element_spme
│   ├── 分支 multi_spme_enabled:
│   │   └── ModelInitialisation_MultiSPMe(case)
│   ├── 分支 !multi_spme_enabled:
│   │   └── ModelInitialisation(case)
│   └── theta = 0.5 (Crank-Nicolson)
│
├── [初始步]
│   └── CallModel(case, y0, t, jacobi="update")       ← 首次系统矩阵计算
│       └── 分支 opt.per_element_spme:
│           └── CallModel_MultiSPMe(...)               ← 多SPMe路径 ⭐
│
├── [时间步循环] while t <= t_end
│   │
│   ├── CallModel(case, y_old, t, jacobi="update")    [CallModel.jl:172]
│   │   └── 分支 opt.per_element_spme → CallModel_MultiSPMe
│   │
│   ├── 时间积分: y_new = Mt \ (Kt * y_old + Ft)
│   │
│   ├── 分支 multi_spme_enabled: 提取温度自由度 T_nodes
│   │
│   ├── ErrorEstimation(case, y_old, y_new, dt_min/dt)
│   │   └── 分支 opt.model=="SPMe": norm(y_new-y_old)/norm(y_old)*coeff
│   │
│   ├── 自适应时间步调整
│   │
│   ├── [CZM损伤更新] 分支 czm_active && step % czm_update_interval == 0:
│   │   └── update_czm_damage!(case, variables, T_nodes_carry)  [CouplingState.jl:331]
│   │
│   ├── [截止电压检测] 分支 multi_spme_enabled:
│   │   └── 检查 n_cutoff_elements, V_cell
│   │
│   └── [整体电压截止]
│       └── V_cell < v_l 或 V_cell > v_h → break
│
├── PostProcessing(case, variables_hist, v)            [PostProcessing.jl]
│
└── [附加输出]
    ├── 耗时统计 (SPMe/branch/thermal/CZM)
    ├── 截止信息
    └── 分支 return_final_state: 保存 final_state
```

---

### A3. CallModel_MultiSPMe 内部流程（核心耦合函数）

```
CallModel_MultiSPMe(case, yt, t; jacobi)              [CallModel.jl:1]
│
├── [1] 解析状态向量
│   ├── get_thermal_dofs(yt, case.layout)              ← 提取温度场 T_nodes
│   └── 提取每个单元的电化学状态 yt_chem[e]
│
├── [2] 提取单元面积和温度
│   └── Te_prev[e] = mean(T_nodes[element_nodes])
│
├── [3] 分流求解
│   ├── I_total = opt.Current(t * t0) / I_typ
│   ├── SPMe_variables(case, yt_representative, t)     ← 代表性电化学计算
│   ├── 分支 opt.czm_enabled: 获取 deactivated_elements
│   └── solve_branch_currents(case, variables, ..., deactivated_elements)
│       ├── compute_prefactors(variables, param, ...)
│       ├── compute_all_coefficients(ne, Te_prev, ...)
│       ├── detect_cutoff_elements(coeffs, ...)
│       └── newton_iteration(I_e, V, ne, w, I_total, ...)
│       → 输出: variables, I_e, Vc
│
├── [4] 并行求解每单元 SPMe  (Threads.@threads)
│   └── SPMe_element(case, yt_chem[e], t, e; I_e, T_e, jacobi, workspace)
│       ├── SPMe_variables!(ws, case, yt_e, t; I_app, T_e)  ← 原位计算
│       ├── 分支 opt.mechanicalmodel=="full":
│       │   └── Mechanicaloutput(case, variables_e)
│       │       └── Calstressdisp(param.NE/PE, mesh, cs, T)
│       ├── ElectrodeDiffusion(param.NE/PE, mesh, nlen, cs_gs, theta_M)
│       ├── ElectrolyteDiffusion(param, mesh_el, nlen, variables_e)
│       └── SPMe_BC(case, variables_e)
│       → 输出: M_e, K_e, F_e, variables_e
│
├── [5] 装配电化学全局矩阵
│   ├── M_chem = blockdiag(M_elems...)
│   ├── K_chem = blockdiag(K_elems...)
│   └── F_chem = vcat(F_elems...)
│
├── [6] 计算逐单元热源
│   ├── 分支 opt.czm_enabled && czm_mesh !== nothing:
│   │   └── compute_heat_sources_with_czm(case, variables, ...)
│   │       └── compute_heat_sources(case, ...)         ← 基本热源
│   │       └── get_active_elements(czm_mesh, geom)     ← CZM活跃单元
│   │       └── 非活跃单元热源置零
│   └── 分支 !czm_enabled:
│       └── compute_heat_sources(case, variables, ..., per_element_spme=true)
│
├── [7] 装配热学矩阵
│   ├── ThermalDistributed2D(case, variables)           [ThermalDistributed.jl:1]
│   │   ├── thermal_capacity_weights_2d(param, fks, ...)
│   │   ├── thermal_anisotropic_conductivity_2d(param, fks, ...)
│   │   └── Assemble() / Assemble1D()
│   │   → 输出: MT, KT, FT
│   └── ThermalDistributed2D_BC(KT, FT, case, t)
│       ├── 分支 opt.czm_enabled: 间隙热阻修正 (compute_gap_conductance)
│       ├── apply_convection_bc!(K, F, case)            ← 外边界对流
│       │   └── 分支 cool_method: 无/不执行
│       └── apply_cool_method!(K, F, mesh, case)
│           ├── 分支 "surface": 全表面对流
│           ├── 分支 "tab": 极耳冷却 (jellyroll_tab_node_indices)
│           └── 分支 "none": 不处理
│       → 输出: KT_bc, FT_bc
│
├── [8] 全局拼装
│   ├── M = blockdiag(M_chem, MT)
│   ├── K = blockdiag(K_chem, KT)
│   └── F = [F_chem; FT]
│
└── [9] 输出
    ├── variables["cell voltage"] = Vc
    ├── variables["temperature"] = thermal2D_volume_average_temperature(...)
    └── 返回: M, K, F, variables, y_phi=[]
```

---

### A4. CZM 损伤更新流程

```
update_czm_damage!(case, variables, T_nodes_carry)   [CouplingState.jl:331]
│
├── compute_czm_effective_params(case)
│   └── 计算有效弹性模量 E_eff、泊松比 ν_eff、热膨胀系数 α_eff、扩散系数 β_n/β_p
│
├── ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)
│   ├── 分支 缓存无效/参数变化:
│   │   └── build_czm_cache(czm_mesh, E_eff, ν_eff, param)
│   │       ├── assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)
│   │       ├── 预计算 CohesiveElementGeom
│   │       └── identify_bc_nodes_czm(czm_mesh, param)
│   └── 返回: cache::CZMAssemblyCache
│
├── compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)
│   └── 提取 dT_elem, Δsoc_n_elem, Δsoc_p_elem
│
├── [粘性正则化] 分支 opt.czm_viscous_enabled && tau_visc > 0:
│   └── visc_beta = delta_s / (tau_visc + delta_s)
│
├── solve_czm_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev; ...)
│   ├── 分支 iter_method=="basic":
│   │   └── solve_czm_basic_step(...)
│   │       └── Newton迭代 + backtrack_line_search!
│   ├── 分支 iter_method=="load_substep":
│   │   └── newton_raphson_czm(...)
│   │       └── 自适应载荷步进 + 回溯线搜索
│   └── 分支 iter_method=="arc_length":
│       └── solve_czm_arc_length_step(...)
│           └── Crisfield弧长法
│   → 内部调用:
│       assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, ...)
│       assemble_czm_system(czm_mesh, u, cohesive_params; ...)
│       assemble_bulk_stiffness(...) [或使用缓存]
│       bilinear_traction_state(δ_n, δ_t, damage_state, params)
│       bilinear_tangent(δ_n, δ_t, damage_state, params)
│       apply_bc_czm(K, R; bc_dofs, bc_vals)
│       update_damage(damage_states, separations, params)
│   → 输出: result::CZMResult, updated_czm_mesh
│
└── 分支 result.converged:
    └── czm_mesh.damage_states = updated_czm_mesh.damage_states
    └── case.czm_layout.u_prev = result.displacement
```

---

### A5. 数据流概要

```
┌──────────────┐     I_total      ┌──────────────────┐
│   Opt.Current │ ──────────────→  │  分流求解器       │
└──────────────┘                   │ solve_branch_currents│
                                   └────────┬─────────┘
                                            │ I_e[e] (逐单元电流)
                                            ▼
┌──────────────┐  yt_chem[e], I_e, T_e  ┌──────────────────┐
│   状态向量    │ ─────────────────────→  │  SPMe_element    │
│   yt (全局)  │                         │  (每单元独立)     │
└──────┬───────┘                         └────────┬─────────┘
       │ T_nodes                                   │ M_e, K_e, F_e
       │ (温度场)                                  │ variables_e
       │                                           │ (含 eta, j, soc, OCV)
       │                                           ▼
       │                              ┌──────────────────────┐
       │                              │  热源计算             │
       │                              │  compute_heat_sources │
       │                              │  compute_heat_sources │
       │                              │  _with_czm            │
       │                              └──────────┬───────────┘
       │                                         │ q_total[e]
       │                                         ▼
       │                              ┌──────────────────────┐
       │                              │ ThermalDistributed2D  │
       │                              │ + BC (含CZM间隙热阻)  │
       │                              └──────────┬───────────┘
       │                                         │ MT, KT, FT
       │                                         ▼
       │                              ┌──────────────────────┐
       │                              │  全局拼装 blockdiag   │
       │                              │  M=blk(M_chem,MT)    │
       │                              │  K=blk(K_chem,KT)    │
       │                              └──────────┬───────────┘
       │                                         │
       ▼                                         ▼
┌──────────────────────────────────────────────────────┐
│              Solve 主循环: 时间积分                     │
│   y_new = Mt \ (Kt * y_old + Ft)                     │
│   → 更新 T_nodes, yt_chem                            │
└──────────────────────┬───────────────────────────────┘
                       │ T_nodes_carry
                       ▼
              ┌──────────────────────┐
              │ CZM 损伤更新          │
              │ update_czm_damage!   │
              │ → D, δ, displacement │
              └──────────┬───────────┘
                         │ 损伤状态 D
                         ▼
              ┌──────────────────────┐
              │ 反馈到热模型          │
              │ - 间隙热阻 h_eff(D)  │
              │ - 活跃单元屏蔽       │
              │ - 分流截止处理       │
              └──────────────────────┘
```

---

## Part B: 模块函数清单

### B1. SetParams.jl — 参数定义与归一化

#### 结构体定义

| 结构体 | 行号 | 用途 |
|--------|------|------|
| `Params` | 217 | 顶层参数容器 (PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale, cohesive) |
| `Electrode` | 37 | 电极参数 (含力学: E, nu, alphaT, Omega) |
| `Separator` | 71 | 隔膜参数 |
| `CurrentCollector` | 81 | 集流体参数 |
| `Electrolyte` | 89 | 电解液参数 (函数型: De, kappa, dlnf_dlnc) |
| `Cell` | 101 | 电池参数 (含 Jellyroll: Rin, Rout, lambda_r/t, Nr_th, Nθ_th) |
| `Tab` | 133 | 极耳参数 |
| `Cohesive` | 148 | CZM参数 (法向/切向/混合模式/界面热阻/粘性) |
| `Scale` | 177 | 统一归一化尺度 (含热/力学尺度) |

#### 函数

##### `ChooseCell(CellType::String)` → `Params` [行231]
- **分支** `CellType=="Jellyroll"`: include("parameters/Jellyroll.jl")
- **输出**: param_dim::Params (有量纲)
- **副作用**: 计算 scale 中所有参考量 (L, t0, phi, sig, P_ref, σ_czm, δ_czm 等)

##### `NormaliseParam(param_dim::Params)` → `Params` [行308]
- **输入**: param_dim (有量纲参数)
- **输出**: param::Params (无量纲参数)
- **关键归一化**:
  - 热参数: `heat_Q* = heat_Q × ρ_ref × L³ × T_ref / (t0 × P_ref)`
  - CZM参数: `σ_max_n* = σ_max_n / σ_czm`, `δ_0_n* = δ_0_n / L`
  - 粘性: 外部通过 `SetCase` 中 `opt.czm_visc_tau / param.scale.t0` 设置 `tau_visc`

---

### B2. Option.jl — 配置选项

#### 结构体 `Option` [行35]

**耦合路径关键开关:**

| 字段 | 类型 | 默认值 | 耦合路径值 | 影响范围 |
|------|------|--------|-----------|---------|
| `model` | String | "SPM" | "SPMe" | SetCase, CallModel, SPMe |
| `thermal_enabled` | Bool | false | true | 变量初始化 |
| `thermalmodel` | String | "none" | "distributed2D" | SetCase, Solve, CallModel |
| `per_element_spme` | Bool | false | true | Solve初始化, CallModel |
| `czm_enabled` | Bool | false | true | CallModel, Solve, BC |
| `mechanicalmodel` | String | "none" | "full" | SPMe, SPMe_element |
| `czm_iter_method` | String | "basic" | "basic" | solve_czm_step |
| `czm_update_interval` | Int | 1 | 1 | Solve主循环 |
| `cool_method` | String | "tab" | "tab"/"surface" | ThermalDistributed2D_BC |
| `czm_model` | String | "model1" | "model1"/"mix" | CZM牵引力模型 |
| `czm_viscous_enabled` | Bool | false | true/false | 粘性正则化 |
| `czm_visc_tau` | Float64 | 0.0 | 10~100s | 松弛时间 |

#### 结构体 `CycleOption` [行5]

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `n_cycles` | Int | 50 | 循环次数 |
| `I_charge/discharge` | Float64 | 5.0 | 充放电电流 |
| `V_upper/lower` | Float64 | 4.2/2.5 | 截止电压 |
| `reset_T_each_cycle` | Bool | true | 每循环重置温度 |

---

### B3. SetCase.jl — 案例构建

#### `SetCase(param_dim::Params, opt::Option, y0=[])` → `Case` [行1]
- **调用**: `NormaliseParam(param_dim)` → param
- **调用**: `SetMesh(...)` → 颗粒/电解液网格
- **调用**: `PickElement(mesh_el, ...)` → 分层电解液网格
- **分支** `opt.model=="SPMe"`: 创建电解液网格和索引
- **分支** `opt.thermalmodel=="distributed2D"`: 设置温度DOF索引
- **粘性归一化**: `param.cohesive.tau_visc = opt.czm_visc_tau / param.scale.t0`
- **输出**: `Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, nothing, nothing)`

#### 结构体 `Case` [行100]
- `param_dim::Params` — 有量纲参数
- `param::Params` — 无量纲参数
- `opt::Option` — 选项
- `mesh::Dict{String, Mesh}` — 网格集合
- `index::Dict{String, Array{Int64}}` — DOF索引
- `layout::Union{Nothing, MultiSPMeLayout}` — 多SPMe布局
- `geometry::Union{Nothing, MeshGeometry}` — 网格几何
- `czm_mesh::Union{Nothing, CohesiveMesh}` — CZM网格
- `czm_cache::Union{Nothing, CZMAssemblyCache}` — CZM装配缓存
- `czm_layout::Union{Nothing, CzmLayout}` — CZM布局+u_prev

---

### B4. CouplingState.jl — 状态布局与耦合辅助

#### 结构体

| 结构体 | 行号 | 用途 |
|--------|------|------|
| `MultiSPMeLayout` | 9 | 多SPMe状态向量布局 (ne, n_chem, nT, chem_range, thermal_range, areas) |
| `BoundaryEdgeCache` | 48 | 预计算的外边界边列表 |
| `MeshGeometry` | 86 | 网格几何拓扑 (element_layer, layer_weights, interface_pairs, czm_element_map, boundary_edges) |
| `CohesiveElementGeom` | 106 | 预计算cohesive单元几何 (length, n_vec, t_vec, R, dofs) |
| `CZMAssemblyWorkspace` | 124 | CZM每轮Newton迭代复用工作区 |
| `CZMAssemblyCache` | 170 | CZM跨时间步缓存 (K_bulk, cohesive_geom, bc_dofs, ws) |
| `CzmLayout` | 193 | CZM布局+跨时间步u_prev |

#### 函数

##### `MultiSPMeLayout(ne, n_chem, nT, mesh_th)` → `MultiSPMeLayout` [行31]
- **输入**: 热单元数, 化学DOF数, 热节点数, 热网格
- **输出**: 布局结构体 (含预计算areas)

##### `compute_boundary_edge_cache(mesh, is_outer)` → `BoundaryEdgeCache` [行59]
- **输入**: 网格, 外边界标记
- **输出**: 去重的外边界边列表和边长

##### `compute_czm_effective_params(case)` → `(E_eff, ν_eff, α_eff, β_n, β_p)` [行227]
- **输入**: case
- **输出**: 有效弹性模量, 泊松比, 热膨胀系数, 扩散应变系数

##### `compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)` → `(dT_elem, Δsoc_n_elem, Δsoc_p_elem)` [行264]
- **输入**: case, 变量字典, CZM网格, 温度场
- **从 variables 提取**: `"thermal2D element soc_n/p"`
- **输出**: 每单元温度变化, SOC变化

##### `update_czm_damage!(case, variables, T_nodes_carry)` → `(u_czm, converged)` [行331]
- **调用**: `compute_czm_effective_params`, `ensure_czm_cache`, `compute_czm_strain_inputs`
- **调用**: `solve_czm_step` (核心求解)
- **分支** `opt.czm_viscous_enabled`: 计算 visc_beta
- **分支** `result.converged`: 提交损伤状态和位移
- **输出**: 位移场, 收敛标志
- **副作用**: 更新 case.czm_mesh.damage_states, case.czm_layout.u_prev

---

### B5. CallModel.jl — 模型调度

##### `CallModel(case, yt, t; jacobi)` → `(M, K, F, variables, y_phi)` [行172]
- **分支** `opt.per_element_spme`: 委托给 CallModel_MultiSPMe ⭐
- **分支** `opt.model=="SPMe"`: 调用 SPMe(case, yt, t)
- **分支** `opt.thermalmodel=="lumped"`: 追加热矩阵
- **输出**: 全局M/K/F矩阵, variables字典, y_phi=[]

##### `CallModel_MultiSPMe(case, yt, t; jacobi)` → `(M, K, F, variables, y_phi)` [行1]
- **调用链**: 见 A3 节详细流程图
- **关键调用**: `solve_branch_currents`, `SPMe_element` (并行), `compute_heat_sources[_with_czm]`, `ThermalDistributed2D`, `ThermalDistributed2D_BC`
- **输出**: blockdiag(M_chem, MT), blockdiag(K_chem, KT), [F_chem; FT], variables

##### `copy_element_results(vars_e)` → `Dict{String,...}` [行208]
- **输入**: 单元变量字典
- **输出**: 轻量级独立拷贝 (含14个关键键)

---

### B6. SPMe.jl — SPMe 模型

##### `SPMe_element(case, yt_e, t, e; I_e, T_e, jacobi, workspace)` → `(M_e, K_e, F_e, variables_e)` [行37]
- **调用**: `SPMe_variables!(ws, ...)` 或 `SPMe_variables(case, ...)`
- **分支** `opt.mechanicalmodel=="full"`: `Mechanicaloutput(case, variables_e)`
- **调用**: `ElectrodeDiffusion(param.NE/PE, mesh, nlen, cs_gs, theta_M)`
- **调用**: `ElectrolyteDiffusion(param, mesh_el, nlen, variables_e)`
- **调用**: `SPMe_BC(case, variables_e)` → F_e
- **输出**: blockdiag(M_np, M_pp, M_el), blockdiag(K_np, K_pp, K_el), F_e, variables_e

##### `SPMe_variables!(ws, case, yt, t; I_app, T_e)` → `ws` [行102]
- **输入**: workspace字典, case, 状态向量, 时间, 电流(可选), 温度(可选)
- **计算**: j_n/j_p, j0_n/j0_p (交换电流), eta_n/eta_p (过电位), dphi_e (电解液电位差), V_cell
- **写入 ws**: cell voltage, exchange current, overpotential, OCP, temperature, etc.

##### `SPMe_variables(case, yt, t; I_app, T_e)` → `variables` [行210]
- 同上逻辑但创建新 Dict（非原位版本）

##### `SPMe_BC(case, variables)` → `F::Vector{Float64}` [行185]
- **输入**: case, variables
- **计算**: 粒子扩散通量 + 电解液源项
- **输出**: 源向量 F = [flux_np; flux_pp; flux_el]

---

### B7. Mechanical.jl — 力学耦合

##### `Mechanicaloutput(case, variables)` → `variables` [行1]
- **分支** `opt.model=="SPMe"`:
  - **调用**: `Calstressdisp(param.NE, mesh_n, c_n, T)` → 应力/位移
  - **调用**: `Calstressdisp(param.PE, mesh_p, c_p, T)` → 应力/位移
  - **修正**: eta_n/p, V_cell (应力耦合修正)
  - **输出**: 更新 variables 中的应力、位移、扩散耦合系数
- **开关**: `opt.mechanicalmodel == "full"`

##### `Calstressdisp(electrode, mesh, cs, T)` → `(σ_r_center, σ_θ_surf, disp_surf, θ_M, cs_gs)` [行112]
- **输入**: 电极参数, 网格, 浓度, 温度
- **输出**: 中心径向应力, 表面切向应力, 表面位移, 应力-扩散耦合系数, 高斯点浓度

##### `thermal_diffusion_stress_2D(case, variables)` → `variables` [行165]
- **输入**: case, variables (含 T_nodes, soc_n/p)
- **调用**: `identify_boundary_nodes(mesh, param, opt)` → 边界节点
- **输出**: 更新 variables 中的应力/位移场 (2D平面应力)

---

### B8. Parallelsolution.jl — 分流求解器

##### `solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements)` → `(variables, I_e, V)` [行358]
- **调用**: `compute_prefactors(variables, param, mesh_ne, mesh_pe)`
- **调用**: `compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref)`
- **调用**: `detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)`
- **调用**: `initialize_currents(ne, w, I_total, x_prev)`
- **调用**: `newton_iteration(I_e, V, ne, w, I_total, coeffs; active_mask)`
- **分支** `deactivated_elements !== nothing`: CZM失效单元电流置零
- **输出**: 更新 variables (含 element current, cutoff info, OCV), I_e, V

##### `compute_prefactors(variables, param, mesh_ne, mesh_pe)` → `NamedTuple` [行23]
- **从 variables 提取**: 颗粒表面浓度, 电解液浓度
- **输出**: prefactor_n/p, csn_av/csp_av, OCP参考值, 固相电导

##### `compute_element_coefficients(e, T_e, param, prefactors, T_ref)` → `NamedTuple` [行59]
- **输出**: C1 (OCV项), C2 (温度项), alpha_p/n (Butler-Volmer系数), C5 (欧姆电阻)

##### `detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)` → `CutoffInfo` [行161]
- **输出**: 活跃掩码, 截止单元列表, OCV分布

##### `newton_iteration(I_e, V, ne, w, I_total, coeffs; active_mask)` → `(V, converged, last_iter)` [行200]
- **分支** `active_mask`: 仅对活跃单元求解
- **调用**: `branch_voltage(coeff, I)`, `branch_dVdI(coeff, I)`, `line_search(...)`

---

### B9. ThermalDistributed.jl — 二维分布式热模型

##### `ThermalDistributed2D(case, variables)` → `(MT, KT, FT)` [行1]
- **输入**: case (含 mesh["thermal2D"]), variables (含 heat_source_fields)
- **调用**: `thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)`
- **调用**: `thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy)`
- **调用**: `Assemble(...)` — 质量矩阵、刚度矩阵
- **输出**: MT (质量), KT (刚度, 含各向异性), FT (载荷, 含热源)

##### `ThermalDistributed2D_BC(KT, FT, case, t)` → `(K, F)` [行286]
- **分支** `opt.czm_enabled`: CZM间隙热阻修正
  - **调用**: `compute_gap_conductance(D, δ_n, param.cohesive)` → h_eff
- **调用**: `apply_convection_bc!(K, F, case; edge_cache)` — 外边界对流
- **调用**: `apply_cool_method!(K, F, mesh, case)`:
  - **分支** `"none"`: 不处理
  - **分支** `"surface"`: 全表面Biot数对流
  - **分支** `"tab"`: 极耳冷却 (jellyroll_tab_node_indices)
- **输出**: 修正后的 K, F

##### `compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme)` → `variables` [行379]
- **分支** `per_element_spme && variables_elems !== nothing`: 从 variables_elems[e] 提取
- **计算**: Q_rxn, Q_rev, Q_ohm (分NE/PE/SP/PCC/NCC各层)
- **输出**: q_total, 各层热源, 写入 variables

##### `compute_heat_sources_with_czm(case, variables, variables_elems, I_e, T_e, areas, czm_mesh, mesh_data)` → `variables` [行512]
- **调用**: `compute_heat_sources(case, ...)` — 基本热源
- **调用**: `get_active_elements(czm_mesh, geom)` — 活跃单元
- **效果**: 非活跃单元热源置零
- **输出**: 更新 variables

---

### B10. czm.jl — CZM 本构模型与网格

#### 结构体

| 结构体 | 行号 | 用途 |
|--------|------|------|
| `CohesiveElement` | 1 | 内聚力单元 (id, nodes[4], nodes_bottom[2], nodes_top[2], length, layer_idx) |
| `DamageState` | 24 | 损伤状态 (D, D_visc, δ_max_n/t/eff, fractured, accumulated_damage) |
| `CohesiveMesh` | (Mutable) | CZM网格 (bulk_mesh, node, bulk_element, cohesive_elements, damage_states, interface_nodes) |

##### `create_czm_mesh(thermal_mesh, param_dim; tol)` → `CohesiveMesh` [行56]
- **输入**: 热网格 (Q4), 参数, 坐标容差
- **算法**: 检测外螺旋与内螺旋的坐标重合节点对
- **输出**: CohesiveMesh (含 cohesive_elements, damage_states)

##### `assemble_czm_system(czm_mesh, u, cohesive_params; damage_states, geom_cache, ws, visc_beta)` → `(K_coh, f_int_coh, separations, tractions)` [行156]
- **核心**: 遍历cohesive单元, 计算B矩阵、分离位移、牵引力、切线刚度
- **调用**: `bilinear_traction_state(δ_n, δ_t, damage_state, params; visc_beta)`
- **调用**: `bilinear_tangent(δ_n, δ_t, damage_state, params; visc_beta)`
- **使用**: CZMAssemblyWorkspace (避免分配)

##### `assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)` → `K_bulk::SparseMatrixCSC` [行331]
- **输入**: CZM网格, 有效弹性模量, 泊松比
- **输出**: 固体Q4单元刚度矩阵

##### `assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n/p_elem)` → `F_thermo_chem::Vector{Float64}` [行397]
- **输入**: CZM网格, 材料参数, 单元温度变化, SOC变化
- **计算**: ε_0 = α*ΔT + β_n*Δsoc_n + β_p*Δsoc_p
- **输出**: 热-化学载荷向量

##### `build_czm_cache(czm_mesh, E_eff, ν_eff, param)` → `CZMAssemblyCache` [行463]
- **调用**: `assemble_bulk_stiffness`, `identify_bc_nodes_czm`
- **预计算**: bulk DOF映射, cohesive单元几何, 边界条件, 工作区
- **输出**: CZMAssemblyCache

##### `ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)` → `cache` [行551]
- **分支** 缓存无效/参数变化: 调用 `build_czm_cache`
- **输出**: 有效的 CZMAssemblyCache

##### `assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; ...)` → `(K_total, f_int_total, separations, tractions)` [行562]
- **调用**: `assemble_bulk_stiffness` (或使用缓存), `assemble_czm_system`
- **输出**: K_total = K_bulk + K_coh, f_int_total = f_int_bulk + f_int_coh

##### `assemble_coupled_system_full(...)` → `(K_total, R, F_thermo_chem, separations, tractions)` [行585]
- **调用**: `assemble_coupled_system`, `assemble_thermal_chemical_load`
- **输出**: K_total, 残差 R = F_ext + F_thermo_chem - f_int_total

##### `apply_bc_czm(K, F; bc_nodes/bc_dofs)` → `(K_new, F_new)` [行608]
- **方法**: 罚函数法 (penalty=1e12)

##### `identify_bc_nodes_czm(czm_mesh, param; opt)` → `(bc_nodes, inner_count, outer_count)` [行642]
- **调用**: `identify_boundary_nodes(czm_mesh, param, opt)`
- **输出**: 内/外边界节点 → :fixed_xy

---

### B11. CzmSolve.jl — CZM 求解器

#### 结构体
- `CZMResult` [行1]: displacement, damage, traction_n/t, separation_n/t, converged, iterations, residual_norm

##### `solve_czm_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev; ...)` → `(result, new_czm_mesh)` [行639]
- **分支** `iter_method=="basic"`: `solve_czm_basic_step`
- **分支** `iter_method=="load_substep"`: `newton_raphson_czm`
- **分支** `iter_method=="arc_length"`: `solve_czm_arc_length_step`
- **公共参数**: α_eff, β_n/β_p, dT_elem, Δsoc_n/p_elem, max_iter, tol, cache, visc_beta

##### `solve_czm_basic_step(...)` → `(result, new_czm_mesh)` [行154]
- **算法**: 标准 Newton-Raphson + backtrack_line_search!
- **调用**: `extract_bc_dofs`, `assemble_thermal_chemical_load`, `assemble_coupled_system`
- **调用**: `apply_bc_czm`, `backtrack_line_search!`, `update_damage`

##### `newton_raphson_czm(...)` → `(result, new_czm_mesh)` [行473]
- **算法**: 自适应载荷步进 + 回溯线搜索
- **特点**: 将总载荷分为 n_load_steps 个子步, 失败时自动减半步长

##### `solve_czm_arc_length_step(...)` → `(result, new_czm_mesh)` [行261]
- **算法**: Crisfield圆柱弧长法
- **特点**: 弧长约束 ||Δu||² = target², 二次方程求根选最近预测方向

##### `update_damage(damage_states, separations, cohesive_params; visc_beta)` → `damage_states` [CzmPostProcess.jl]
- **分支** `visc_beta < 1.0`: 粘性正则化 D_visc = β*D + (1-β)*D_visc_prev
- **更新**: δ_max_n/t/eff, D, fractured 标志

##### `clone_damage_states(damage_states)` → `Vector{DamageState}` [行17]
- **输出**: 深拷贝

##### `clone_czm_mesh_with_damage(czm_mesh, damage_states)` → `CohesiveMesh` [行31]
- **输出**: 共享几何/拓扑, 替换 damage_states

---

### B12. Solve.jl — 主求解器

##### `Solve(case; initial_state, return_final_state, ...)` → `result::Dict` [行1]
- **分支** `opt.model=="thermal"`: 纯热求解 (跳过耦合路径)
- **分支** `per_element_spme`: `ModelInitialisation_MultiSPMe(case)` / `ModelInitialisation(case)`
- **分支** `initial_state !== nothing`: 从外部状态恢复
- **主循环**: 时间步进 + CallModel + CZM更新 + 截止检测
- **输出**: PostProcessing结果 + 耗时统计 + 截止信息 + final_state (可选)

##### `ErrorEstimation(case, y_old, y_new, coeff)` → `Float64` [行423]
- **分支** `opt.model=="SPMe"`: norm(y_new-y_old)/norm(y_old)*coeff

##### `RecordMatrix!(case, M, K)` → `case` [行413]
- **分支** `opt.jacobi=="constant"`: 缓存 M/K 矩阵

---

### B13. Jellyrollmodel.jl — Jellyroll 几何

##### `jellyroll_collector_seed_mesh(param; nθ, gsorder, phase, tol)` → `JellyrollMesh` [行17]
- **算法**: 阿基米德螺旋线 r(θ) = a + bθ, collector-seeded网格
- **输出**: JellyrollMesh (含 thermal2D mesh, interface_pairs, czm_element_map, element_layer, layer_weights, tab_nodes)

##### `setup_thermal2D_mesh(case, mesh_data)` → `case` (修改)
- **功能**: 将 JellyrollMesh 数据写入 case, 创建 layout/geometry

##### `jellyroll_element_properties(mesh, param)` → `(element_layer, layer_weights)`
- **输出**: 每个单元的层类型和5层面积权重 [NE, SP, PE, PCC, NCC]

##### `jellyroll_tab_node_indices(mesh, param)` → `(pos_idx, neg_idx)`
- **输出**: 正/负极耳对应的节点索引

##### `thermal2D_volume_average_temperature(mesh, T_nodes)` → `Float64`
- **输出**: 体积加权平均温度

---

### B14. CycleSolver.jl — 循环求解器

##### `solve_cycling(case, cycle_opt)` → `CyclingResult`
- **循环**: 充电 → 静置 → 放电 → 静置
- **每阶段**: `solve_phase(case, ...)`
- **跨周期**: final_state 传递, damage 累积
- **分支** `reset_T_each_cycle`: 温度重置
- **分支** `soh <= czm_soh_threshold`: 终止

##### `solve_phase(case; phase_type, I_func, t_end, dt, initial_state)` → `PhaseResult`
- **调用**: `Solve(case; initial_state, return_final_state=true)`
- **输出**: PhaseResult (含 final_state)

---

## 附录: 完整开关→函数映射表

| 开关条件 | 影响的函数 | 跳过的分支 |
|----------|-----------|-----------|
| `opt.model == "SPMe"` | SetCase, CallModel, SPMe, ErrorEstimation | SPM, P2D, sP2D 分支 |
| `opt.per_element_spme` | Solve, CallModel → CallModel_MultiSPMe | 单模型路径 |
| `opt.thermalmodel == "distributed2D"` | SetCase, Solve, ThermalDistributed2D | lumped, ring2D, ring2D_polar |
| `opt.czm_enabled` | CallModel, Solve, ThermalDistributed2D_BC | 无CZM路径 |
| `opt.mechanicalmodel == "full"` | SPMe_element → Mechanicaloutput | 无力学耦合 |
| `opt.cool_method == "tab"` | apply_cool_method! → jellyroll_tab_node_indices | surface, none |
| `opt.czm_iter_method` | solve_czm_step → 3种求解器 | 其他方法 |
| `opt.czm_model == "model1"` | bilinear_traction_state (仅Mode I) | mix (Mode I+II) |
| `opt.czm_viscous_enabled` | update_czm_damage! → visc_beta | 无粘性正则化 |
