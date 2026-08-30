# 力学结构体四层重构实施计划（Mechanics Struct Refactor）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 按 spec 删除 `Cohesive`/`CzmInterfaceParams`/`CzmParamCache`/`CzmLayout`/`CZMAssemblyCache`，界面参数挂 `CurrentCollector`，缓存挂 `CohesiveMesh`，演化状态聚合 `MechState`，`opt.czm` 嵌套收敛——物理零变更，基线 v7 逐位一致。

**Architecture:** 大爆炸单批：Task 2–8 构成一个连续的 break-fix 弧（中间态不可加载、**不得**以恢复旧符号来"修复"加载错误），Task 9 全量验证门，Task 10 文档，最后单次提交。纯重构无新行为——既有 34 测试 + 基线 v7 即验收标准，不新增测试文件。

**Tech Stack:** Julia 1.11.2（`D:/Julia-1.11.2/bin/julia.exe`），`@with_kw`（Parameters.jl），SparseArrays。

**Spec:** `docs/superpowers/specs/2026-08-30-mechanics-struct-refactor-design.md`（字段映射表 §6、契约 §5、验证 §7——执行者必须同时持有 spec）。

## Global Constraints

- 环境固定：`GKSwstype=100`、`JULIA_NUM_THREADS=1`、`--startup-file=no --project=.`
- 验收 = testexample 与基线 `testexample-20260830T005629+0800` **逐位一致**（电压 4.0367→3.9438 V、温度 298.15~299.00 K、环向应力 −1.3952~+3.8570 MPa、分离 9.4407e-13 m、19 步、零损伤）+ runtests 34/34 + couple_example 四 PNG 哈希 `540fe42f/ecdc9f58/2c29e35b/5124dec3`
- **参数冻结契约**：SetCase 归一化后不改 `param`（写入 AGENTS，Task 10）
- 新代码风格：新函数无下划线前缀/感叹号后缀（除既有 `!` 惯例外——`update_czm_damage!`/`commit!` 沿用）；不新增防御性断言（AGENTS 9.7）
- commit 守卫拦截 `git commit`：最终提交步骤需用户授权执行
- **前置条件**：工作树已有两个未提交已过门批次（旋转清理、α/β 分层化）；本计划执行前应先提交它们（用户授权），否则 diff 混流
- 界面↔集流体映射：`:PE_PCC → param.PCC`、`:NE_NCC → param.NCC`（层序中 PCC 两面皆 PE、NCC 两面皆 NE）

---

### Task 1: Option 层收敛（opt.czm 嵌套，全库读写点机械替换）

**Files:**
- Modify: `src/Option.jl`（新增 `CzmOptions`，删 20 个 `czm_*` 平铺字段）
- Modify: 所有 `opt.czm_*` / `case.opt.czm_*` 读写点（`grep -rln "czm_" src/ test/ example/` 枚举）

**Interfaces:**
- Produces: `CzmOptions` struct（字段与默认值见下）；`Option.czm::CzmOptions`。后续任务以 `case.opt.czm.<field>` 读取。
- Consumes: 无（首个任务）。

- [ ] **Step 1: 定义 CzmOptions 并替换 Option 字段**

`src/Option.jl` 中删除 20 个 `czm_*` 字段（`czm_model/czm_enabled/czm_update_interval/czm_soh_threshold/czm_inner_exit_only/czm_fix_inner/czm_iter_method/czm_max_iter/czm_tol/czm_load_steps/czm_arc_length_alpha/czm_viscous_enabled/czm_visc_tau/czm_area_loss_enabled/czm_area_loss_threshold/czm_geo_nonlinear/czm_winding_prestress/czm_j2_plasticity/czm_continuous_feedback/czm_friction_mu`），新增：

```julia
@with_kw mutable struct CzmOptions
    enabled::Bool = false                # 是否启用CZM损伤模型
    model::String = "model1"             # "model1" or "mix"
    update_interval::Int64 = 1           # 损伤更新间隔（时间步数）
    soh_threshold::Float64 = 0.8         # SOH终止阈值
    inner_exit_only::Bool = true         # 断裂时仅内圈单元退出电化学反应
    fix_inner::Bool = true               # 边界：true=内外圈均固定
    iter_method::String = "basic"        # "basic" | "load_substep" | "arc_length"
    max_iter::Int64 = 100
    tol::Float64 = 1e-4
    load_steps::Int64 = 2                # 载荷子步数（load_substep）
    arc_length_alpha::Float64 = 1.0      # 弧长法系数
    viscous_enabled::Bool = false
    viscous_tau::Float64 = 0.0           # 物理松弛时间 [s]
    area_loss_enabled::Bool = false
    area_loss_threshold::Float64 = 0.83
    geo_nonlinear::Bool = false          # GL/TL 残差 + K_G（Batch 2）
    winding_prestress::Bool = false      # 卷绕预应力 σ₀（Batch 2'）
    j2_plasticity::Bool = false          # PCC/NCC J2（Batch 3）
    continuous_feedback::Bool = false    # 损伤–电–热反馈（Batch 6）
    friction_mu::Float64 = 0.10          # SP Coulomb 摩擦（Batch 8 预留）
end
```

`Option` struct 增加一行：`czm::CzmOptions = CzmOptions()`。

- [ ] **Step 2: 全库机械替换读写点**

按映射 `opt.czm_X → opt.czm.X`、`case.opt.czm_X → case.opt.czm.X`（注意 `czm_model→czm.model`、`czm_visc_tau→czm.viscous_tau`）。逐一处理：

```bash
grep -rn "opt\.czm_\|\.czm_model\|\.czm_visc_tau" src/ test/ example/ --include="*.jl"
```

已知热点：`src/Solve.jl`（czm_update_interval/geo_nonlinear/debug 块）、`src/CouplingState.jl`（update_czm_damage! 读取全部求解选项、`czm_params.czm_model = case.opt.czm_model` 反向突变行**删除**——model 从此只活在 opt）、`src/CycleSolver.jl`（soh_threshold/update_interval）、`src/SetCase.jl`（`param.cohesive.tau_visc = opt.czm_visc_tau / param.scale.t0` 行**删除**，visc_beta 改在 CouplingState 从 `case.opt.czm.viscous_tau / param.scale.t0` 现算）、`src/Parallelsolution.jl`、`src/Variables.jl`、全部 `test/test_czm_*.jl` 与 `example/*.jl` 中 `opt.czm_X =` 赋值。

- [ ] **Step 3: 验证**

```bash
grep -rn "opt\.czm_\|\.czm_model\b\|czm_visc_tau" src/ test/ example/ --include="*.jl"
# 期望：无输出（czm_model 只允许出现在 opt.czm.model 与 bilinear 显式参数——后者 Task 3 建）
julia --startup-file=no -t 1 -e 'include("src/JuBat.jl"); using .JuBat; o = JuBat.Option(); @assert o.czm.iter_method == "basic"; println("OK")'
```

（此时 `param.cohesive` 仍在——Task 2 处理；若 Step 3 加载失败且报错指向 `param.cohesive`，那是 Task 2 范围，先确认报错非本任务引入再继续。）

---

### Task 2: 参数层——CurrentCollector 承载界面，Cohesive 删除

**Files:**
- Modify: `src/SetParams.jl`（CurrentCollector +16 字段；删 `Cohesive` struct 与 `Params.cohesive`；NormaliseParam 归一化块；ChooseCell 锚点）
- Modify: `src/parameters/Jellyroll.jl`（界面参数从 `cohesive.X_pe_pcc = …` 改写到 `PCC.X = …`）
- Modify: `src/JuBat.jl`（导出表删 `CzmInterfaceParams, CzmParamCache`——本任务连同结构体一起删）

**Interfaces:**
- Produces: `CurrentCollector` 含界面字段（名：`σ_max, K_n, δ_0, G_c, δ_c, τ_max, K_t, δ_0_t, G_c_t, δ_c_t, eta, h_c0, k_air, lambda_m, beta, threshold`；默认值 `0.0`，唯 `eta=1.0, h_c0=1e7, k_air=0.026, lambda_m=70e-9, beta=1.0, threshold=70e-9`）；归一化后量纲契约同 spec §4.1。
- Consumes: Task 1 的 `opt.czm.model/viscous_tau`。

- [ ] **Step 1: CurrentCollector 增加界面字段**

在 `CurrentCollector`（`src/SetParams.jl:86` 附近，`sigma_y/H` 塑性字段之后）追加：

```julia
    # 界面本构（该集流体与其相邻涂层之间的 COH 界面；法向无后缀、切向 _t）
    σ_max::Float64 = 0.0    # 最大法向牵引 [Pa]
    K_n::Float64 = 0.0      # 法向初始刚度 [Pa/m]
    δ_0::Float64 = 0.0      # 法向损伤起始位移 [m]
    G_c::Float64 = 0.0      # 法向断裂能 [J/m²]
    δ_c::Float64 = 0.0      # 法向临界位移 [m]
    τ_max::Float64 = 0.0    # Mode II 最大切向牵引 [Pa]
    K_t::Float64 = 0.0
    δ_0_t::Float64 = 0.0
    G_c_t::Float64 = 0.0
    δ_c_t::Float64 = 0.0
    eta::Float64 = 1.0      # BK 混合模式指数 [-]
    # 界面热阻（md/07）
    h_c0::Float64 = 1e7
    k_air::Float64 = 0.026
    lambda_m::Float64 = 70e-9
    beta::Float64 = 1.0
    threshold::Float64 = 70e-9
```

- [ ] **Step 2: 删除 Cohesive 与归一化/锚点迁移**

删除 `mutable struct Cohesive`（SetParams.jl:161-198）、`Params` 的 `cohesive::Cohesive` 字段、NormaliseParam 中 `param.cohesive.*` 归一化块（SetParams.jl:503 起整块）、`ChooseCell` 的 cohesive 锚点块（:342-355）。替换为：

```julia
# ChooseCell 锚点（PE-PCC 界面 = PCC 为锚）
if param_dim.PCC.σ_max == 0 || param_dim.PCC.G_c == 0
    @warn "[ChooseCell] PCC.σ_max 或 PCC.G_c = 0; scale.σ_czm/δ_czm/G_czm/K_czm 无法锚定。"
end
param_dim.scale.σ_czm = param_dim.PCC.σ_max
if param_dim.PCC.σ_max > 0 && param_dim.PCC.G_c > 0
    param_dim.scale.δ_czm = 2 * param_dim.PCC.G_c / param_dim.PCC.σ_max
else
    param_dim.scale.δ_czm = param_dim.scale.L
end
param_dim.scale.G_czm = param_dim.scale.σ_czm * param_dim.scale.δ_czm
param_dim.scale.K_czm = param_dim.scale.σ_czm / param_dim.scale.δ_czm
```

```julia
# NormaliseParam 界面字段归一化（逐字段显式赋值，沿用既有 per-field 惯用式；锚点见 ChooseCell）
param.PCC.σ_max = param_dim.PCC.σ_max / param_dim.scale.σ_czm
param.PCC.τ_max = param_dim.PCC.τ_max / param_dim.scale.σ_czm
param.PCC.K_n = param_dim.PCC.K_n / param_dim.scale.K_czm
param.PCC.K_t = param_dim.PCC.K_t / param_dim.scale.K_czm
param.PCC.δ_0 = param_dim.PCC.δ_0 / param_dim.scale.δ_czm
param.PCC.δ_c = param_dim.PCC.δ_c / param_dim.scale.δ_czm
param.PCC.δ_0_t = param_dim.PCC.δ_0_t / param_dim.scale.δ_czm
param.PCC.δ_c_t = param_dim.PCC.δ_c_t / param_dim.scale.δ_czm
param.PCC.G_c = param_dim.PCC.G_c / param_dim.scale.G_czm
param.PCC.G_c_t = param_dim.PCC.G_c_t / param_dim.scale.G_czm
param.PCC.eta = param_dim.PCC.eta
param.PCC.h_c0 = param_dim.PCC.h_c0
param.PCC.k_air = param_dim.PCC.k_air
param.PCC.lambda_m = param_dim.PCC.lambda_m
param.PCC.beta = param_dim.PCC.beta
param.PCC.threshold = param_dim.PCC.threshold
param.NCC.σ_max = param_dim.NCC.σ_max / param_dim.scale.σ_czm
param.NCC.τ_max = param_dim.NCC.τ_max / param_dim.scale.σ_czm
param.NCC.K_n = param_dim.NCC.K_n / param_dim.scale.K_czm
param.NCC.K_t = param_dim.NCC.K_t / param_dim.scale.K_czm
param.NCC.δ_0 = param_dim.NCC.δ_0 / param_dim.scale.δ_czm
param.NCC.δ_c = param_dim.NCC.δ_c / param_dim.scale.δ_czm
param.NCC.δ_0_t = param_dim.NCC.δ_0_t / param_dim.scale.δ_czm
param.NCC.δ_c_t = param_dim.NCC.δ_c_t / param_dim.scale.δ_czm
param.NCC.G_c = param_dim.NCC.G_c / param_dim.scale.G_czm
param.NCC.G_c_t = param_dim.NCC.G_c_t / param_dim.scale.G_czm
param.NCC.eta = param_dim.NCC.eta
param.NCC.h_c0 = param_dim.NCC.h_c0
param.NCC.k_air = param_dim.NCC.k_air
param.NCC.lambda_m = param_dim.NCC.lambda_m
param.NCC.beta = param_dim.NCC.beta
param.NCC.threshold = param_dim.NCC.threshold
```

（`eta/h_c0/k_air/lambda_m/beta/threshold` 无因次或已在 L 空间，直接拷贝——与现行 Cohesive 归一化块对这些字段的处理一致。）

- [ ] **Step 3: 参数集迁移（Jellyroll.jl）**

`cohesive.X_pe_pcc = v` → `PCC.X = v`、`cohesive.X_ne_ncc = v` → `NCC.X = v`；共用字段写两份：

```julia
# PCC：PE-PCC 界面（数值沿用现 cohesive._pe_pcc 值）
PCC.σ_max = …; PCC.K_n = …; PCC.δ_0 = …; PCC.G_c = …; PCC.δ_c = …
PCC.τ_max = …; PCC.K_t = …; PCC.δ_0_t = …; PCC.G_c_t = …; PCC.δ_c_t = …
PCC.eta = 1.45; PCC.h_c0 = …; PCC.k_air = …; PCC.lambda_m = …; PCC.beta = …; PCC.threshold = …
# NCC：NE-NCC 界面（同构；eta/热阻五参数与 PCC 相同值）
NCC.σ_max = …; …; NCC.eta = 1.45; NCC.h_c0 = …; …
# cohesive = Cohesive() 整块删除
```

其余参数集（LGM50/Ring/Enertech/Northrop）未赋 cohesive，不动。同时删除 `src/JuBat.jl` 中 `CzmInterfaceParams/CzmParamCache` 导出及 struct 定义本身（CouplingState.jl:30-79——`compute_czm_params_per_interface` 仍被引用，函数体与调用点的删除在 Task 4，此处**先保留函数但让其编译**：它读 `param.cohesive.*` 会断——故本任务内直接把 `CzmInterfaceParams/CzmParamCache/compute_czm_params_per_interface` 定义与 `case.czm_param_cache` 字段一并删除，断点留给 Task 4 接）。

- [ ] **Step 4: 验证（符号残留检查）**

```bash
grep -rn "Cohesive(\|\.cohesive\b\|cohesive\.\|CzmInterfaceParams\|CzmParamCache\|czm_param_cache" src/ --include="*.jl"
# 期望：仅剩 Task 3/4/5 将改写的消费点清单（czm.jl/CzmSolve.jl/CouplingState.jl/Mechanical.jl/CycleSolver.jl/SetCase.jl），
# 不再有任何 struct 定义或参数集赋值
```

（本任务后加载必然失败——break-fix 弧开始，继续 Task 3。）

---

### Task 3: 本构与间隙导热签名迁移（Materialmatrix.jl）

**Files:**
- Modify: `src/Materialmatrix.jl`（bilinear_*、update_damage、compute_gap_conductance）
- Modify: `src/CouplingState.jl`（`compute_czm_params_per_interface` 已在 Task 2 删除；K_0 下界 @warn 块的宿主处理）

**Interfaces:**
- Produces（Task 4/7/8 消费）:
  - `bilinear_traction_state(δ_n::Float64, δ_t::Float64, damage_state::DamageState, ip::CurrentCollector, czm_model::String; visc_beta::Float64=1.0) -> (T_n, T_t, D_eq, new_state)`
  - `bilinear_tangent(δ_n, δ_t, damage_state, ip::CurrentCollector, czm_model::String; visc_beta=1.0) -> dT_dδ::Matrix{Float64}`
  - `update_damage(damage_states, separations, ip::CurrentCollector, czm_model::String; visc_beta=1.0) -> Vector{DamageState}`
  - `compute_gap_conductance(D::Float64, δ_n::Float64, ip::CurrentCollector, param::Params) -> h_eff`
- Consumes: Task 2 的 CurrentCollector 字段（`ip.δ_0/ip.δ_c/ip.δ_0_t/ip.δ_c_t/ip.eta/ip.K_n/ip.K_t/ip.σ_max/τ_max/G_c/G_c_t`、热阻五参数）。

- [ ] **Step 1: bilinear_traction_state 改造**

签名与字段读取映射：`params.δ_0_n→ip.δ_0`、`params.δ_c_n→ip.δ_c`、`params.δ_0_t/δ_c_t` 同名、`params.η→ip.eta`、`params.K_n/K_t/σ_max` 同名；`czm_model` 从第 5 位置参数读（原 `params.czm_model`）。函数体逻辑零变更，例：

```julia
function bilinear_traction_state(δ_n::Float64, δ_t::Float64, damage_state::DamageState,
                                 ip::CurrentCollector, czm_model::String; visc_beta::Float64=1.0)
    K_n = ip.K_n; K_t = ip.K_t
    δ_0_n = ip.δ_0; δ_c_n = ip.δ_c
    δ_0_t = ip.δ_0_t; δ_c_t = ip.δ_c_t
    η = ip.eta
    # …函数体逐行保留，仅上述变量来源变化…
```

- [ ] **Step 2: bilinear_tangent / update_damage 同构改造**（同 Step 1 映射）。

- [ ] **Step 3: compute_gap_conductance 改造**

```julia
function compute_gap_conductance(D::Float64, δ_n::Float64, ip::CurrentCollector, param::Params)
    # Λ 内联（spec §4.1）：分离空间（δ_czm 归一）→ 热模型长度空间（L 归一）
    inv_Λ = param.scale.δ_czm / param.scale.L
    delta0 = ip.δ_0 * inv_Λ
    delta_c = ip.δ_c * inv_Λ
    delta = max(δ_n, 0.0) * inv_Λ
    D_clamped = clamp(D, 0.0, 0.9999)
    two_beta_lambda = 2.0 * ip.beta * ip.lambda_m
    h_eff = if delta < delta0
        ip.h_c0 + ip.k_air / (delta + two_beta_lambda)
    elseif delta < ip.threshold
        ip.h_c0 * (1.0 - D_clamped) + ip.k_air / (delta + two_beta_lambda)
    else
        ip.h_c0 * (1.0 - D_clamped) + ip.k_air / (delta + delta0)
    end
    h_eff > 0 || error("compute_gap_conductance: zero or negative conductance")
    return h_eff
end
```

（其余逐行保留原逻辑；`compute_element_gap_conductance/compute_all_gap_conductances` 增加 `param` 透传参数。）`CouplingState.jl` 原 K_0 下界 @warn（δ_0*>0.1）迁移至 `NormaliseParam` 末尾（读 `param.PCC.δ_0 / param.NCC.δ_0`）。

- [ ] **Step 4: 验证**

```bash
grep -n "params\.\|CzmInterfaceParams" src/Materialmatrix.jl
# 期望：无输出（本文件不再出现 params::CzmInterfaceParams 任何痕迹）
```

---

### Task 4: 装配与求解链直读 param，Λ 内联，参数缓存删除

**Files:**
- Modify: `src/czm.jl`（`assemble_czm_system/assemble_coupled_system/assemble_coupled_system_full/assemble_bulk_stiffness/assemble_thermal_chemical_load/assemble_bulk_residual_tangent/gl_element…调用`，~43 处 `param_cache`）
- Modify: `src/CzmSolve.jl`（~39 处 `param_cache`，签名 `param_cache::CzmParamCache` → `param::Params`）
- Modify: `src/CouplingState.jl`（`update_czm_damage!` 主干）
- Modify: `src/Mechanical.jl`（若残留 param_cache 引用）
- Modify: `src/SetCase.jl`（Case 字段：删 `czm_param_cache`、`czm_cache`、`czm_layout`——`czm_layout` 字段删除放 Task 6，本任务先删 `czm_param_cache`）

**Interfaces:**
- Produces: 装配族签名 `assemble_czm_system(czm_mesh, u, param::Params; damage_states=nothing, geom_cache=nothing, ws=nothing, visc_beta=1.0, czm_model="model1", …)`（其余 kwarg 同现状）；`assemble_thermal_chemical_load(czm_mesh, param, dT, Δsn, Δsp)`（批次①已定型）；`solve_czm_step` 族签名本任务只做 `param_cache→param` 替换，**完整收敛签名在 Task 6**。
- Consumes: Task 3 的 bilinear 族新签名。

- [ ] **Step 1: assemble_czm_system 热路径改造**

关键片段（循环外取 Λ、循环内取 ip）：

```julia
function assemble_czm_system(czm_mesh::CohesiveMesh, u::Vector{Float64}, param::Params;
        damage_states=nothing, geom_cache=nothing, ws=nothing,
        visc_beta::Float64=1.0, czm_model::String="model1")
    # …
    Λ = param.scale.L / param.scale.δ_czm   # 内联，不再存字段
    # …循环内：
        iface = czm_mesh.cohesive_elements[i].interface_type
        ip = iface === :PE_PCC ? param.PCC : param.NCC
        # …
        T_n, T_t, _, _ = bilinear_traction_state(δ_n, δ_t, damage_state, ip, czm_model; visc_beta=visc_beta)
        dT_dδ = bilinear_tangent(δ_n, δ_t, damage_state, ip, czm_model; visc_beta=visc_beta)
```

其余逐行保留（`params.Λ` 的两处使用点都换成内联 Λ；`param_cache.param_ref` 全部改 `param`）。`czm_model` 由调用链从 `opt.czm.model` 传入（Task 6 完成接线前，装配函数带默认值 "model1" 保证编译）。

- [ ] **Step 2: 其余装配/求解函数同构替换**

`assemble_bulk_stiffness(czm_mesh, param)`、`assemble_thermal_chemical_load(czm_mesh, param, …)`（已是）、`assemble_coupled_system(_full)`、`assemble_bulk_residual_tangent(czm_mesh, u, param; …)`、`update_damage_per_interface(czm_mesh, damage_states, separations, param, czm_model; visc_beta)`、CzmSolve 五个求解函数 `param_cache::CzmParamCache, param` 两参 → 单 `param::Params`。机械规则：

```bash
grep -n "param_cache" src/czm.jl src/CzmSolve.jl src/CouplingState.jl src/Mechanical.jl src/SetCase.jl
# 逐处：param_cache → param（签名）；param_cache.param_ref → param（取参）
```

`build_czm_cache/ensure_czm_cache/CZMAssemblyCache` 定义与调用点本任务**暂留**（Task 5 删），但其 `param_cache` 参数改传 `param`。

- [ ] **Step 3: update_czm_damage! 主干适配**

删除 `czm_param_cache = compute_czm_params_per_interface(case)` 块与 `case.czm_param_cache` 读写；`ensure_czm_cache` 调用暂改为传参形态（Task 5 删）；`α_eff/β` 提取已不存在（批次①）；`visc_beta` 改为：

```julia
czm_opt = case.opt.czm
visc_beta = 1.0
if czm_opt.viscous_enabled && czm_opt.viscous_tau > 0.0
    delta_s = lowercase(czm_opt.iter_method) == "basic" ? 1.0 : 1.0 / max(1, czm_opt.load_steps)
    visc_beta = (czm_opt.viscous_tau / param.scale.t0) / ((czm_opt.viscous_tau / param.scale.t0) + delta_s)
end
```

- [ ] **Step 4: 验证**

```bash
grep -rn "param_cache\|by_interface\|compute_czm_params_per_interface\|\.Λ\b" src/ --include="*.jl"
# 期望：仅剩 Task 5/6 将处理的 ensure_czm_cache/build_czm_cache/CZMAssemblyCache 内部（互相引用成孤岛，无外部消费）
```

---

### Task 5: 网格缓存挂载 + BC 现算，缓存机制删除

**Files:**
- Modify: `src/CzmMesh.jl`（CohesiveMesh 新增 3 字段）
- Modify: `src/czm.jl`（删 `CZMAssemblyCache/build_czm_cache/ensure_czm_cache`；K_bulk/geom/ws 惰性构建）
- Modify: `src/CzmBC.jl` / BC 提取路径（每次 solve 入口现算）
- Modify: `src/CzmSolve.jl`（`extract_bc_dofs(czm_mesh, param; fix_inner)` 现算）

**Interfaces:**
- Produces: `CohesiveMesh.K_bulk::Union{Nothing,SparseMatrixCSC}`、`.cohesive_geom::Union{Nothing,Vector{CohesiveElementGeom}}`、`.ws::Union{Nothing,CZMAssemblyWorkspace}`（惰性：`=== nothing` 时构建一次）；`bulk_stiffness(czm_mesh, param)`、`cohesive_geometry(czm_mesh)`、`assembly_workspace(czm_mesh)` 三个惰性访问器。失效 = 对象身份。
- Consumes: Task 4 的 `assemble_bulk_stiffness(czm_mesh, param)`。

- [ ] **Step 1: CohesiveMesh 增加缓存字段**（`src/CzmMesh.jl`，mutable struct 内）：

```julia
    K_bulk::Union{Nothing, SparseMatrixCSC{Float64, Int64}} = nothing        # 惰性装配缓存（弹性路径）
    cohesive_geom::Union{Nothing, Vector{CohesiveElementGeom}} = nothing     # 纯几何标架（gs 同款）
    ws::Union{Nothing, CZMAssemblyWorkspace} = nothing                       # 预分配工作区
```

（若 CohesiveMesh 用裸 `mutable struct` 无默认值构造，则在 `create_czm_mesh` 组装处显式置 `nothing`。）

- [ ] **Step 2: 三个惰性访问器 + 删除旧缓存机制**

```julia
function bulk_stiffness(czm_mesh::CohesiveMesh, param::Params)
    czm_mesh.K_bulk === nothing && (czm_mesh.K_bulk = assemble_bulk_stiffness(czm_mesh, param))
    return czm_mesh.K_bulk
end
function cohesive_geometry(czm_mesh::CohesiveMesh)
    czm_mesh.cohesive_geom === nothing && (czm_mesh.cohesive_geom = build_cohesive_geometry(czm_mesh))
    return czm_mesh.cohesive_geom
end
function assembly_workspace(czm_mesh::CohesiveMesh)
    czm_mesh.ws === nothing && (czm_mesh.ws = CZMAssemblyWorkspace(2 * czm_mesh.nnode, czm_mesh.n_cohesive))
    return czm_mesh.ws
end
```

`build_cohesive_geometry` = 现 `build_czm_cache` 第 3 段（czm.jl:757-775）抽出。删除 `CZMAssemblyCache` struct、`build_czm_cache`、`ensure_czm_cache`、`case.czm_cache` 字段及全部引用；`JuBat.jl` 导出表删 `build_czm_cache, ensure_czm_cache`。

- [ ] **Step 3: 调用点切换**

`cache.K_bulk → bulk_stiffness(czm_mesh, param)`、`cache.cohesive_geom → cohesive_geometry(czm_mesh)`、`cache.ws → assembly_workspace(czm_mesh)`；`extract_bc_dofs(czm_mesh, param; cache=cache)` 改为内部每次调 `identify_bc_nodes_czm(czm_mesh, param; fix_inner=fix_inner)`（fix_inner 从 `czm_opt.fix_inner` 传入求解函数）。geo_nl 路径继续 `K_bulk_cached=nothing` 语义（不触碰 `czm_mesh.K_bulk`）。

- [ ] **Step 4: 验证**

```bash
grep -rn "CZMAssemblyCache\|ensure_czm_cache\|build_czm_cache\|czm_cache" src/ test/ example/ --include="*.jl"
# 期望：无输出
julia --startup-file=no -t 1 -e 'include("src/JuBat.jl")'
# 期望：仍可能因 czm_layout（Task 6）失败；报错不得指向本任务符号
```

---

### Task 6: MechState 聚合 + 提交语义 + 克隆链删除 + Case 收敛

**Files:**
- Modify: `src/CouplingState.jl`（新 `MechState`；删 `CzmLayout`；`update_czm_damage!` 完整收敛）
- Modify: `src/czm.jl` / `src/CzmSolve.jl`（damage_states 从 mesh 迁 ms；删 `clone_czm_mesh_with_damage`；求解返回 `result::CZMResult` 单值）
- Modify: `src/SetCase.jl`（Case：删 `czm_layout`，增 `mech`；SetCase 末尾在 czm_mesh 建立后 `case.mech = MechState(czm_mesh)`——注意 czm_mesh 在 `setup_thermal2D_mesh` 后由调用方挂载，mech 创建放 `update_czm_damage!` 首次惰性初始化：`case.mech === nothing && (case.mech = MechState(czm_mesh))`）
- Modify: `src/Solve.jl`、`src/CycleSolver.jl`、`src/PostProcessing.jl`、`src/CsvExport.jl`、`src/Variables.jl`（`case.czm_layout.X → case.mech.X`、`czm_mesh.damage_states → case.mech.damage_states`、`ensure_node_ref!(case) → case.mech.node_ref` 惰性块）

**Interfaces:**
- Produces:
  - `mutable struct MechState`（spec §4.3 字段表，含 `contact::Nothing`）
  - `solve_czm_step(czm_mesh, ms::MechState, param, F_ext, czm_opt::CzmOptions; dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing) -> CZMResult`（**终态签名**；basic/load_substep/arc 三个内部分发同构）
  - 求解器收敛后原位提交：`ms.u_prev = u`、`ms.damage_states = new_states`；失败不触碰 ms
- Consumes: Task 1 `CzmOptions`、Task 4/5 装配与惰性缓存。

- [ ] **Step 1: MechState 定义与 CzmLayout 删除**

```julia
mutable struct MechState
    u_prev::Vector{Float64}
    node_ref::Union{Nothing, Matrix{Float64}}
    damage_states::Vector{DamageState}
    plastic_states::Union{Nothing, Matrix{PlasticState}}
    winding_prestress::Union{Nothing, Vector{Tuple{Float64,Float64,Float64}}}
    contact::Nothing
end
function MechState(czm_mesh::CohesiveMesh)
    MechState(zeros(2 * czm_mesh.nnode), nothing,
              [DamageState() for _ in 1:czm_mesh.n_cohesive], nothing, nothing, nothing)
end
```

`CzmLayout/init_czm_layout` 删除；`create_czm_mesh` 不再初始化 `czm_mesh.damage_states`（字段从 CohesiveMesh 删除）。

- [ ] **Step 2: 求解器提交语义**

各求解函数：入口 `damage_states = copy.(ms.damage_states)`（或现有 clone_damage_states）作局部试探态；**收敛后** `ms.damage_states = damage_states; ms.u_prev = copy(u)`；塑性 `commit_plastic` 路径原位写 `ms.plastic_states`（时机不变）。返回值 `result, updated_czm_mesh → result`；删除 `clone_czm_mesh_with_damage`。调用方（`update_czm_damage!`、CycleSolver、测试）改单返回值。

- [ ] **Step 3: update_czm_damage! 终态收敛**（spec §4.5 全貌）：

```julia
res = solve_czm_step(case.czm_mesh, case.mech, case.param, F_ext, case.opt.czm;
                     dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem)
```

调试块与有限性检查逐行保留；`core_ovalization` 的 `ref_node` 改由 `case.mech.node_ref`（惰性：首次为 nothing 时 `case.mech.node_ref = copy(case.czm_mesh.node)`）。

- [ ] **Step 4: 验证（弧闭合——首次全量加载）**

```bash
grep -rn "CzmLayout\|czm_layout\|clone_czm_mesh_with_damage\|czm_mesh.damage_states\|czm_param_cache\|czm_cache\b" src/ test/ example/ --include="*.jl"
# 期望：无输出
julia --startup-file=no -t 1 -e 'include("src/JuBat.jl"); using .JuBat; println("load OK")'
# 期望：load OK
```

---

### Task 7: 测试套适配（8 个 CZM 测试 + runtests 编排）

**Files:**
- Modify: `test/test_czm_mech_core.jl`、`test_czm_geometric_stiffness.jl`、`test_czm_geo_c1.jl`、`test_czm_arc_geo.jl`、`test_czm_j2_integration.jl`、`test_czm_multicycle_c4lite.jl`、`test_czm_winding_prestress.jl`、`test/unit_czm_eigenstrain.jl`、`test/unit_czm_newton.jl`（若引用缓存/布局）

**Interfaces:**
- Consumes: Task 1–6 全部终态签名。
- Produces: 无。

- [ ] **Step 1: 机械模式适配**（四类样板，逐文件套用）：

```julia
# (a) Option 赋值：opt.czm_geo_nonlinear = true  →  opt.czm.geo_nonlinear = true
# (b) 缓存获取：cache = JuBat.ensure_czm_cache(case, czm_mesh, pc)  →  删除；K_bulk/geom/ws 经惰性访问器或直接不传
# (c) 求解调用：solve_czm_step(czm_mesh, F_ext, pc, param, u0; kw...) →
#     ms = JuBat.MechState(czm_mesh); res = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; dT_elem=…, …)
#     多步推进时复用同一 ms（u_prev/损伤在 ms 内累积）
# (d) param_cache/by_interface：
#     pc = JuBat.compute_czm_params_per_interface(case) → 删除
#     cache.by_interface[:PE_PCC].σ_max → case.param.PCC.σ_max
#     cache.by_interface[:NE_NCC].δ_0_n → case.param.NCC.δ_0
```

`unit_czm_eigenstrain.jl` 的 `retune_czm_for_damage!` 重写（调谐归一化 param 实例）：

```julia
function retune_czm_for_damage!(param; σ_scale::Float64, δc_over_δ0::Float64)
    # PCC（PE-PCC 界面）
    σ_pcc = param.PCC.σ_max * σ_scale
    τ_pcc = param.PCC.τ_max * σ_scale
    δ0_pcc = param.PCC.δ_c / δc_over_δ0
    δ0_pcc_t = param.PCC.δ_c_t / δc_over_δ0
    param.PCC.σ_max = σ_pcc
    param.PCC.K_n = σ_pcc / δ0_pcc
    param.PCC.δ_0 = δ0_pcc
    param.PCC.G_c = 0.5 * σ_pcc * param.PCC.δ_c
    param.PCC.τ_max = τ_pcc
    param.PCC.K_t = τ_pcc / δ0_pcc_t
    param.PCC.δ_0_t = δ0_pcc_t
    param.PCC.G_c_t = 0.5 * τ_pcc * param.PCC.δ_c_t
    # NCC（NE-NCC 界面）
    σ_ncc = param.NCC.σ_max * σ_scale
    τ_ncc = param.NCC.τ_max * σ_scale
    δ0_ncc = param.NCC.δ_c / δc_over_δ0
    δ0_ncc_t = param.NCC.δ_c_t / δc_over_δ0
    param.NCC.σ_max = σ_ncc
    param.NCC.K_n = σ_ncc / δ0_ncc
    param.NCC.δ_0 = δ0_ncc
    param.NCC.G_c = 0.5 * σ_ncc * param.NCC.δ_c
    param.NCC.τ_max = τ_ncc
    param.NCC.K_t = τ_ncc / δ0_ncc_t
    param.NCC.δ_0_t = δ0_ncc_t
    param.NCC.G_c_t = 0.5 * τ_ncc * param.NCC.δ_c_t
end
retune_czm_for_damage!(param; σ_scale=σ_scale, δc_over_δ0=δc_over_δ0)
```

（注意：param 冻结契约针对生产路径；测试内调谐参数是既有测试做法的延续，加一行注释说明。`assemble_thermal_chemical_load(czm_mesh, cache, …)` → `(czm_mesh, param, …)`。）

- [ ] **Step 2: 验证**

```bash
for f in test_czm_mech_core test_czm_geometric_stiffness test_czm_geo_c1 test_czm_arc_geo test_czm_j2_integration test_czm_multicycle_c4lite test_czm_winding_prestress unit_czm_eigenstrain; do
  GKSwstype=100 JULIA_NUM_THREADS=1 "D:/Julia-1.11.2/bin/julia.exe" --startup-file=no --project=. test/$f.jl 2>&1 | tail -2
done
# 期望：全部 PASS（断言集不变——纯 API 适配，数值应与适配前一致）
```

---

### Task 8: 示例脚本适配

**Files:**
- Modify: `example/testexample.jl`（opt.czm_* 赋值——若 Task 1 未覆盖）、`example/couple_example.jl`、`example/czm_cycle_example.jl`、`example/coupled_czm_thermal_example.jl`（L244 `czm_param_cache.by_interface[…]` → `param.PCC/param.NCC` + `compute_element_gap_conductance` 新签名）、`example/内聚力验证/verify_czm_per_interface.jl`（同 Task 7(d) 模式）、`example/jellyroll_stress_displacement.jl`（若有引用）

**Interfaces:** Consumes Task 1–6 终态。Produces: 无。

- [ ] **Step 1: 逐脚本适配**

```bash
grep -rln "czm_param_cache\|by_interface\|opt\.czm_\|czm_layout\|ensure_czm_cache\|\.damage_states" example/ --include="*.jl"
# 逐文件按 Task 7 的 (a)-(d) 模式替换
```

`coupled_czm_thermal_example.jl` 间隙导热段示例：

```julia
h_eff_all = [JuBat.compute_element_gap_conductance(czm_mesh, i, JuBat.collector_params(case.param, czm_mesh.cohesive_elements[i].interface_type), case.param)
             for i in 1:czm_mesh.n_cohesive]
```

（`collector_params(param, iface) = iface === :PE_PCC ? param.PCC : param.NCC`——Task 4 中在 czm.jl 定义并导出的公共小助手。）

- [ ] **Step 2: 验证**

```bash
grep -rn "czm_param_cache\|by_interface\|opt\.czm_\b\|ensure_czm_cache\|CzmLayout" example/ --include="*.jl"
# 期望：无输出
```

---

### Task 9: 全量验证门（v7 逐位一致 + 34/34 + PNG）

**Files:** 无代码改动（失败则回到对应 Task 修复后重跑本门）。

- [ ] **Step 1: 全套测试**

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 "D:/Julia-1.11.2/bin/julia.exe" --startup-file=no --project=. test/runtests.jl 2>&1 | tail -3
# 期望：34/34 PASS, exit 0
```

- [ ] **Step 2: testexample 快门 vs 基线 v7 逐位比对**

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 "D:/Julia-1.11.2/bin/julia.exe" --startup-file=no --project=. example/testexample.jl
# 期望（对照 Simplify/baseline/testexample/metrics.toml）：
#   4.0367 → 3.9438 V；0.0833 Ah；298.15~299.00 K；19 步；D 0.0000%；分离 9.4407e-13 m；
#   环向 −1.3952~+3.8570 MPa；切向 −0.3739~+0.3952 MPa；exit 0
# 任一指标漂移 = 回对应 Task 定位，禁止"数值接近"放行
```

- [ ] **Step 3: couple_example PNG 门**

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 "D:/Julia-1.11.2/bin/julia.exe" --startup-file=no --project=. example/couple_example.jl
sha256sum output/couple_example/final_*.png
# 期望四哈希：540fe42f… / ecdc9f58… / 2c29e35b… / 5124dec3…（v7 记录）
```

- [ ] **Step 4: 三重门全绿后进入 Task 10；否则定位修复并从 Step 1 重跑。**

---

### Task 10: 文档同步与提交

**Files:**
- Modify: `AGENTS.md` §9.4（重写 CZM 入口描述：CurrentCollector 承载界面、参数冻结契约、MechState/缓存挂载）、§5.3 选项表（opt.czm 嵌套）
- Modify: `md/01_参数定义与归一化.md`、`md/06_内聚力模型_CZM.md`、`md/07_界面热阻模型.md`、`md/15_颗粒与极片模量区分.md`、`md/对照/06_CZM对照.md`、`md/源码函数索引/`（czm/CzmSolve/CouplingState/CzmMesh/Mechanical/SetParams/Option 七篇）
- Modify: `src/parameters/Jellyroll.jl` 顶部注释（若提及 cohesive 块）

- [ ] **Step 1: 文档更新**（核心语义一句话：界面属性挂集流体、参数 SetCase 后冻结、缓存随网格、状态在 mech、Λ 内联）。

- [ ] **Step 2: 提交（需用户逐次授权 commit 守卫）**

```bash
git add -A
git commit -m "refactor(mech): four-layer struct refactor per 2026-08-30 spec

- interface params on CurrentCollector (PCC<->:PE_PCC, NCC<->:NE_NCC); Cohesive/CzmInterfaceParams/CzmParamCache deleted
- czm_model/tau_visc single-source in opt.czm (CzmOptions nested, 20 fields)
- assembly caches lazily on CohesiveMesh (K_bulk/geom/ws); BC computed per solve; content-hash invalidation deleted
- MechState aggregates evolution state (damage/plastic/prestress/u_prev/node_ref + contact slot); clone chain deleted; converged in-place commit
- Lambda inlined scale.L/δ_czm; verification: runtests 34/34, testexample bit-identical to v7, 4 PNG hashes unchanged"
```

---

## Self-Review 记录

1. **Spec 覆盖**：spec §4.1→Task 2/3/4；§4.2→Task 5；§4.3→Task 6；§4.4→Task 1/6；§4.5→Task 6 Step 3；§5 契约→Task 10（AGENTS）+Task 6（提交语义）；§6 映射表→Task 1/2/7(d)；§7 验证→Task 9；§8 顺序→任务序；§9 风险→前置条件。无遗漏。
2. **占位符扫描**：参数集迁移（Task 2 Step 3）中 `PCC.σ_max = …` 为"沿用现值"的映射指令而非占位——执行者从现行 `cohesive.σ_max_pe_pcc` 值逐字搬入，已显式说明。
3. **类型一致性**：`ip::CurrentCollector`、`czm_model::String` 第 5 参、`solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; 3 载荷)`、`collector_params(param, iface)` 在 Task 3/4/6/8 间一致；`opt.czm.model/viscous_tau` 与 Task 1 字段名一致。
