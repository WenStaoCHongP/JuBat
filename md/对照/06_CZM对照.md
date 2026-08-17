# 卷 06 — CZM 对照

> 理论来源（2 个 md）：
> - `md/06_内聚力模型_CZM.md`（539 行）— 双线性牵引-分离律、CZMResult 结构、热-化学载荷、无量纲化重设计 v2、界面导热
> - `md/14_粘性正则化.md`（125 行）— 卷末追加：粘性正则化 visc_beta、前向欧拉离散、一致性切线
>
> 代码来源（5 个 src 文件 + 2 辅助）：
> - `src/CzmMesh.jl`（182 行）— CohesiveElement / create_czm_mesh
> - `src/Czm.jl`（629 行）— DamageState / moduli_of / assemble_czm_system / assemble_bulk_stiffness / assemble_thermal_chemical_load / assemble_coupled_system(_full) / cache
> - `src/CzmBC.jl`（63 行）— apply_bc_czm / identify_bc_nodes_czm
> - `src/CzmSolve.jl`（675 行）— CZMResult / solve_czm_basic_step / newton_raphson_czm / solve_czm_arc_length_step / solve_czm_step / backtrack_line_search! / build_arc_length_augmented_matrix
> - `src/Mechanical.jl`（360 行）— Calstressdisp（颗粒扩散应力）/ thermal_diffusion_stress_2D（2D 宏观热-扩散应力）
> - `src/CzmPostProcess.jl`（117 行）— get_damage_statistics / check_fracture_criterion / czm_output_to_variables
> - `src/CouplingState.jl`（762 行）— CzmInterfaceParams / CzmParamCache / compute_czm_params_per_interface / compute_czm_strain_inputs / update_czm_damage!
> - 辅助：`src/Materialmatrix.jl:60-427`（bilinear_traction_state / bilinear_tangent / update_damage / compute_gap_conductance）；`src/SetParams.jl:156-192`（Cohesive 参数集）
>
> 生成日期：2026-08-02

---

## 编写说明

本卷"实现摘要"列常含简短代码片段（如 `D_visc = D_visc + β·(D_eq − D_visc)`、`K_coh += wJΛ·BL_dT·B`）以保证对照的可验证性，因此字数可能超过 spec §3 建议的 25 字。此为有意识取舍——优先**可追溯**而非**摘要性**。

---

## md 06 §1 模型概述

> **术语基线**：四个真实箔–涂层面复用 `:PE_PCC`、`:NE_NCC` 两种 `interface_type`；`CohesiveMesh.n_layers=2` 是遗留类型计数字段。完整模型 `n_cohesive=4*N_seg`，而 `phi_pairs` 只保存跨匝配对。

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/06 §1.1-1.2 物理背景 + 3 条假设（界面单元预设 / 渐进损伤 / 影响力学&热学） | `src/CzmMesh.jl:3-12 (CohesiveElement)` + `src/Czm.jl:17-28 (DamageState)` | CohesiveElement 4 节点界面单元；DamageState 含 D/D_visc/δ_max/fractured | ✅ | md §1 为概念陈述，非公式；代码组织对应：预设 cohesive 单元 + 渐进 D ∈ [0,1] + D 影响导热（见 §5） |

---

## md 06 §2 本构模型

### §2.1 双线性牵引-分离律

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.1) | md/06 §2.1 弹性段 `T = K·δ`（δ < δ_0） | `src/Materialmatrix.jl:122-123,148-152 (bilinear_traction_state)` | `δ_eff ≤ δ_0_eff → D_eq=0`；`T_n=(1−D_visc)·K_n·δ_n`，D_visc=0 时即 K_n·δ_n | ✅ | md 弹性段不含损伤；代码用 `δ_eff ≤ δ_0_eff → D_eq=0` 分支判定 |
| (6.2) | md/06 §2.1 软化段 `T = (1−D)·K·δ`（δ_0 ≤ δ < δ_c） | `src/Materialmatrix.jl:127,149,157 (bilinear_traction_state)` | `D_eq = δ_c·(δ−δ_0)/(δ·(δ_c−δ_0))`；`T=(1−D_visc)·K·δ` | 🟡 | **md 公式用 D，代码实际用 D_visc**；粘性关闭（visc_beta=1）时一致；启用时牵引基于 D_visc（见 md 14 系列）；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |
| (6.3) | md/06 §2.1 损伤变量定义 `D = δ_c·(δ−δ_0)/(δ·(δ_c−δ_0))` | `src/Materialmatrix.jl:127 (bilinear_traction_state)` + `:226,269 (bilinear_tangent)` | 两处都写 `δ_c_eff·(δ_eff−δ_0_eff)/(δ_eff·(δ_c_eff−δ_0_eff))`，逐字对应 | ✅ | md L47 公式；代码计算后写回 `new_state.D = D_eq` |
| (6.4) | md/06 §2.1 完全断裂 `T=0, D=1`（δ ≥ δ_c） | `src/Materialmatrix.jl:88-100,124-125,137-139 (bilinear_traction_state)` | `fractured → T_n=0,T_t=0; D=1`；`δ_eff≥δ_c_eff → D_eq=1; D_eq≥1−1e-10 → fractured=true` | ✅ | md L53；代码用 `1−1e-10` 阈值防数值漂移；model1 模式下 T_t 保留 K_t·δ_t（mode I only） |
| (6.x1) | md/06 §2.1 法向/切向耦合的等效分离（mode mix） `δ_eff = sqrt(δ_n² + δ_t²)` | `src/Materialmatrix.jl:107 (bilinear_traction_state)` + `:209 (bilinear_tangent)` | `δ_eff = sqrt(δ_n_pos² + δ_t²)`（仅 czm_model != "model1"） | ✅ | md 未显式列此式，但 §2.3 BK 准则隐含；代码 model1 模式下退化为 δ_eff=δ_n_pos |

### §2.2 参数定义

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.5a) | md/06 §2.2 `δ_0 = T_max / K` | `src/SetParams.jl:156-192 (Cohesive)` + `src/CouplingState.jl:343-393 (compute_czm_params_per_interface)` | Cohesive 直接存 K_n / σ_max / δ_0_n（独立字段），CzmInterfaceParams 原样转发 | ✅ | md L63 关系式；参数在 ChooseCell/NormaliseParam 阶段计算 |
| (6.5b) | md/06 §2.2 `δ_c = 2G_c / T_max` | `src/CouplingState.jl:343-393 (compute_czm_params_per_interface)` | `pe_pcc = CzmInterfaceParams(σ_max=..., K_n=..., δ_c_n=coh.δ_c_pe_pcc, ...)` | ✅ | md L64；锚定 δ_czm = 2·G_c/σ_max（见 (6.20)） |

### §2.3 BK 混合模式准则

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.6) | md/06 §2.3 BK 准则 `δ_c^mix = [(G_n/G_c^I)^α + (G_t/G_c^II)^α]^(1/α)·δ_c^I` | — | src 未直接实现 BK 形式；改用 δ_0_eff/δ_c_eff 的幂指数插值（见 (6.x2)） | 🟡 | md L73 标准 BK 形式；代码用等效分离 + β^η 插值；指向 `docs/planning-with-files/内聚力模块归一化检查/findings.md` |
| (6.x2) | md 06 §2.3 隐含：δ_0_eff/δ_c_eff 通过 `β^η` 插值 | `src/Materialmatrix.jl:108-115 (bilinear_traction_state)` + `:210-217 (bilinear_tangent)` | `β=abs(δ_t)/δ_eff`；`δ_0_eff=sqrt(δ_0_n²+(δ_0_t²−δ_0_n²)·β^η)`；δ_c_eff 同 | ✅ | md 未显式列此式，但 BK 的工程近似；η 默认 1.45（CzmInterfaceParams.η） |

### §2.4 bilinear_traction_state 函数

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.7) | md/06 §2.4 签名 `function bilinear_traction_state(δ_n, δ_t, damage_state, cohesive_params; visc_beta=1.0)` | `src/Materialmatrix.jl:68 (bilinear_traction_state)` | 签名逐字对应：`(δ_n::Float64, δ_t::Float64, damage_state::DamageState, params::CzmInterfaceParams; visc_beta::Float64=1.0)` | ✅ | md L81；参数类型已从 Cohesive 改为 CzmInterfaceParams（per-interface 化） |
| (6.x3) | md/06 §2.4 返回 (T_n, T_t, D_eq, new_state) | `src/Materialmatrix.jl:160 (bilinear_traction_state)` | `return T_n, T_t, D_eq, new_state` | ✅ | md L86-88；new_state 含 D_visc/δ_max_eff/accumulated_damage/fractured |
| (6.x4) | md/06 §2.4 粘性正则化：`T_n = (1−D_visc)·K_n·δ_n` | `src/Materialmatrix.jl:149 (bilinear_traction_state)` | `δ_n≥0 → T_n=(1−D_visc)·K_n·δ_n`；δ_n<0 → T_n=K_n·δ_n（接触刚度保留） | 🟡 | **md (6.x4) 写 D_visc 但 md 06 §2.1 公式 (6.2) 仍写 D**；这是 md 内部不一致（见 md 14 §6.2）；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |

---

## md 06 §3 有限元实现

### §3.1 单元数据结构

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.8) | md/06 §3.1 CohesiveElement：nodes / length / normal | `src/CzmMesh.jl:3-12 (CohesiveElement)` | `id/nodes/nodes_bottom/nodes_top/length/interface_type/host_outer_elem/host_inner_elem` | ✅ | md L99-102；代码比 md 多 host_outer_elem/host_inner_elem（cohesive_to_thermal 映射用） |
| (6.9) | md/06 §3.1 DamageState：D/δ_max_n/δ_max_t/D_visc/δ_max_eff/accumulated_damage/fractured | `src/Czm.jl:17-28 (DamageState)` | 7 字段全部对应 + 默认构造函数 `DamageState() = new(0,0,0,0,0,false,0)` | ✅ | md L103-110；代码字段顺序：D/D_visc/δ_max_n/δ_max_t/δ_max_eff/fractured/accumulated_damage |

### §3.2 CZMResult 结构体

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.10) | md/06 §3.2 `mutable struct CZMResult` 9 字段 | `src/CzmSolve.jl:1-15 (CZMResult)` | 9 字段逐字对应：displacement/damage/traction_n/traction_t/separation_n/separation_t/converged/iterations/residual_norm | ✅ | md L117-127；构造函数 `CZMResult(ndof,n_coh) = new(zeros...zeros..., false, 0, Inf)` 也对应 |
| (6.11) | md/06 §3.2 默认构造函数 `CZMResult(ndof, n_coh)` | `src/CzmSolve.jl:12-14` | `CZMResult(ndof::Int, n_coh::Int) = new(zeros(ndof), zeros(n_coh), ..., false, 0, Inf)` | ✅ | md L133-136 |

### §3.3 单元刚度矩阵

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.12) | md/06 §3.3 分离计算 `δ = B·u_element` | `src/Czm.jl:201-205 (assemble_czm_system)` | `B_local = R·B_global`；`δ_local = B_local·u_e`；`δ_n = Λ·δ_local[1]` | ✅ | md L144；代码含 Λ 换算因子（重设计 v2，见 (6.22)），与 md §4.8 一致 |
| (6.13) | md/06 §3.3 切线刚度 `K_coh = ∫_L B^T·D_tan·B dL` | `src/Czm.jl:213-224 (assemble_czm_system)` | `BL_dT = B_local'·dT_dδ`；`K_e += wJΛ·BL_dT·B_local`（乘一次 Λ） | ✅ | md L151；dT_dδ 由 bilinear_tangent 返回 2×2 矩阵；wJ=w·J=w·(L/2) |
| (6.14) | md/06 §3.3 内力向量 `f_int = ∫_L B^T·T dL` | `src/Czm.jl:226-232 (assemble_czm_system)` | `T_vec=[T_n,T_t]`；`BLtT = B_local'·T_vec`；`f_int_e += wJ·BLtT`（**不乘** Λ） | ✅ | md L159；重设计 v2 §5 虚功一致性：δ̃=Λ·B·ũ；f=∫BᵀT̃ dΓ 不乘 Λ；切线乘一次 Λ |

### §3.4 系统组装

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.15) | md/06 §3.4 全局刚度 `K_global = Σ_e K_coh,e` | `src/Czm.jl:252-258 (assemble_czm_system)` + `:581-582 (assemble_coupled_system)` | `K_coh[dofs[a],dofs[b]] += K_e[a,b]`；`K_total = K_bulk + K_coh` | ✅ | md L167；K_bulk 由 assemble_bulk_stiffness 提供（见 (6.18)） |
| (6.16) | md/06 §3.4 全局内力 `f_global = Σ_e f_int,e` | `src/Czm.jl:253-254 (assemble_czm_system)` + `:584-585 (assemble_coupled_system)` | `f_int_coh[dofs[a]] += f_int_e[a]`；`f_int_total = K_bulk·u + f_int_coh` | ✅ | md L173；体单元线性弹性：f_int_bulk = K_bulk·u |

---

## md 06 §4 求解方法

### §4.1 求解框架

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.17) | md/06 §4.1 统一入口 `solve_czm_step(...; iter_method="load_substep", ...)` | `src/CzmSolve.jl:651-675 (solve_czm_step)` | 三分支：`load_substep → newton_raphson_czm`；`basic → solve_czm_basic_step`；`arc_length* → solve_czm_arc_length_step` | 🟡 | md L187 默认 `iter_method="load_substep"`；代码默认也是 `load_substep`；但 md §4.4 推荐 basic；指向 `docs/planning-with-files/损伤异常/findings.md` |
| (6.x5) | md/06 §4.1 残差格式 `R = F_applied − f_int(u)` | `src/CzmSolve.jl:199,313,350,535 (multiple solvers)` + `:807 (assemble_coupled_system_full)` | `R = F_external + F_thermo_chem − f_int_total`（含热-化学载荷） | ✅ | md L195；F_applied 隐含 F_ext + F_thermo_chem |

### §4.2 损伤更新策略（隐式冻结）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x6) | md/06 §4.2 隐式冻结：NR 迭代用 D_begin，收敛后统一更新 | `src/CzmSolve.jl:177 (damage_start=clone), 214 (converged→update), 237-240 (回退)` | basic：`damage_states = czm_mesh.damage_states`（冻结）；收敛后 `damage_states = update_damage_per_interface(...)`；未收敛回退 `damage_start` | ✅ | md L205-210；三种方法（basic/load_sub/arc_length）均采用：`if result.converged → update_damage_per_interface` |

### §4.3 迭代方法

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x7) | md/06 §4.3.1 basic：`solve_czm_basic_step(...)` + 8 次减半线搜索 | `src/CzmSolve.jl:168-263 (solve_czm_basic_step)` + `:107-128 (backtrack_line_search!)` | `α=1.0; for _ in 1:8: R_trial<R_current → accepted; α*=0.5` | ✅ | md L219；从 u_prev 出发，一次性施加 F_ext+F_thermo_chem |
| (6.x8) | md/06 §4.3.2 load_substep：`F_applied = f_int_ref + α·F_Δ` | `src/CzmSolve.jl:500-520 (newton_raphson_czm)` | `f_int_ref = f_int(u_prev)`；`F_target = F_ext+F_thermo_chem`；`F_delta = F_target − f_int_ref`；`F_applied = f_int_ref + target_progress·F_delta` | ✅ | md L239；α 从 0→1 逐步推进 |
| (6.x9) | md/06 §4.3.2 自适应子步：成功×1.25 / 失败×0.5 | `src/CzmSolve.jl:548,600 (newton_raphson_czm)` | `converged → step_size = min(step_size·1.25, step_size_max)`；`failed → step_size *= 0.5; if < step_size_min → break` | ✅ | md L248-249；step_size_min = step_size/128 |
| (6.x10) | md/06 §4.3.2 收敛判据：子步 tol×10 / 整体 tol×100 | `src/CzmSolve.jl:543,621 (newton_raphson_czm)` | `substep_tol = tol·10.0`；`final_tol = tol·100.0`；`result.converged = load_progress≥1−1e-12 && R_norm<final_tol` | ✅ | md L253-254 |
| (6.x11) | md/06 §4.3.3 arc_length：Crisfield 圆柱约束 `‖Δu‖² = l_arc²` | `src/CzmSolve.jl:265-466 (solve_czm_arc_length_step)` | `arc_target = sqrt(sum(abs2,delta_u_pred))`；`arc_constraint = dot(delta_u,delta_u) − arc_target²` | ✅ | md L267；预测步 + 校正步（NR + 二次方程求根） |
| (6.x12) | md/06 §4.3.3 弧长校正：解两个线性系统 + 二次方程 | `src/CzmSolve.jl:371-410 (solve_czm_arc_length_step)` | `delta_u_R = K_bc\R_bc`；`delta_u_F = K_bc\F_load_bc`；二次方程 a_q·dl²+b_q·dl+c_q=0；选 dot 最大根 | ✅ | md L272-276；Crisfield 标准算法 |

### §4.4 三种方法性能对比

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/06 §4.4 性能表（basic 5.5s / load_sub 22s / arc 22s） | — | md 给 1C 放电 nθ=16 实测；代码无显式基准（性能特征由实现体现） | — | 卡点：md §4.4 表 L284-288 为经验数据，非公式；指向 `docs/planning-with-files/CZM瓶颈/findings.md` |

### §4.5 热-化学载荷

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.18) | md/06 §4.5 初始应变 `ε_0 = α·ΔT + β_n·Δsoc_n + β_p·Δsoc_p` | `src/Czm.jl:373 (assemble_thermal_chemical_load)` | `epsilon_0_elem[e] = α_eff·dT_elem[e] + β_n·Δsoc_n_elem[e] + β_p·Δsoc_p_elem[e]` | ✅ | md L313；α_eff 跨材料统一（spec §7.1）；β_n/β_p 由 `update_czm_damage!` 传 `param.Omega/3` |
| (6.19) | md/06 §4.5 等效节点力 `F_thermo_chem = ∫ B^T·D·ε_0 dΩ` | `src/Czm.jl:393-404 (assemble_thermal_chemical_load)` | `factor = E/(1−ν²)·ε_0·(1+ν)·w·detJ`；`f_e[2i−1] += dNdx[i]·factor`；`f_e[2i] += dNdy[i]·factor` | ✅ | md L319；平面应力 D 矩阵 + 各向同性 ε_0；代码用 dNdx/dNdy 高斯积分 |

### §4.6 完整耦合系统组装

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x13) | md/06 §4.6 `assemble_coupled_system(...)` 返回 K_total / f_int / sep / trac | `src/Czm.jl:556-588 (assemble_coupled_system)` | `K_total = K_bulk + K_coh`；`f_int_total = K_bulk·u + f_int_coh`；返回四元组 | ✅ | md L339-343；含 K_bulk_cached/geom_cache/ws 缓存入参（spec v2） |
| (6.x14) | md/06 §4.6 完整残差版 `assemble_coupled_system_full(...)` | `src/Czm.jl:596-629 (assemble_coupled_system_full)` | 调 `assemble_coupled_system` + `assemble_thermal_chemical_load`；`R = F_external + F_thermo_chem − f_int_total` | ✅ | md 未列函数体；代码组合两个子函数 + 残差 |

### §4.7 应变驱动的有效模量（极片层级）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x15) | md/06 §4.7 入口 `compute_effective_coating_modulus(case)` | — | src **未找到此函数**；实际等效入口是 `compute_czm_params_per_interface(case)`（per-interface 分组） | — | md L349 函数名与 src 不一致；指向 CLAUDE.md §9.4 提及；卡点：md 描述的 API 不存在；改由 `CouplingState.jl:302-401` 实现 |
| (6.x16) | md/06 §4.7 `compute_czm_effective_params` 调用后归一化逆变换 | — | src **未找到此函数**；归一化由 `compute_czm_params_per_interface` 内联：`E_eff_pe = param.PE.E_coat·scale.E_coat/scale.σ_czm` | — | md L363 描述的函数不存在；指向 CLAUDE.md §9.4；卡点：md API 失效 |
| (6.20) | md/06 §4.7 全栈厚度加权 `E_eff = Σ(E_i·t_i)/Σt_i` | `src/CouplingState.jl:316-317 (compute_czm_params_per_interface)` | `E_eff_pe = param.PE.E_coat·scale.E_coat/scale.σ_czm`（PE 单层，非全栈加权） | 🟡 | md L353 写"PE+NE+SP+PCC+NCC 五层加权"；代码实际**只取 PE.E_coat 单层**（per-interface 分组，PE_PCC 界面用 PE）；指向 `docs/planning-with-files/力学模块修改/宏观力学模块无量纲化重设计.md` |
| (6.x17) | md/06 §4.7 调用路径：`thermal_diffusion_stress_2D` 调 `compute_effective_coating_modulus` | `src/Mechanical.jl:184-187 (thermal_diffusion_stress_2D)` | 实际调用：`czm_param_cache = compute_czm_params_per_interface(case)`；取 `by_interface[:PE_PCC].E_eff/ν/α` | 🟡 | md L362 路径描述与代码不符；功能等价（取 PE_PCC 占位） |
| (6.x18) | md/06 §4.7 参数缺失告警 `PE.E_coat==0 → @warn` | `src/CouplingState.jl:307 (compute_czm_params_per_interface)` + `src/Mechanical.jl:167 (thermal_diffusion_stress_2D)` | CouplingState 用 `@assert param.PE.E_coat > 0`；Mechanical 同 | 🟡 | md L365 说"ChooseCell 发 @warn，后续 @assert 拦截"；实际两处都是 @assert（无 ChooseCell 告警）；指向 `docs/planning-with-files/力学模块修改/宏观层面力学模块无量纲化.md` |

### §4.8 无量纲化重设计 v2（2026-07-22）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.21) | md/06 §4.8 锚点 `δ_czm = 2·G_c/σ_max` → σ_max*=1, δ_c*=1, G_c*=1/2 | `src/CouplingState.jl:308-313,321-322 (compute_czm_params_per_interface)` | 入口断言 `σ_max_pe_pcc>0`；`Λ = scale.L/scale.δ_czm`；CzmInterfaceParams.σ_max/K_n/δ_c_n/G_c 直接用归一化值 | ✅ | md L372；scale.σ_czm/scale.δ_czm 在 NormaliseParam 中由 cohesive 参数锚定 |
| (6.22) | md/06 §4.8 装配换算因子 Λ = scale.L/scale.δ_czm | `src/CouplingState.jl:322 (compute_czm_params_per_interface)` | `Λ = scale.L / scale.δ_czm`；存于 `CzmInterfaceParams.Λ` | ✅ | md L379；spec v2 §3.5.2 |
| (6.x19) | md/06 §4.8 分离计算含 Λ：`δ̃ = Λ·B·ũ` | `src/Czm.jl:204-205 (assemble_czm_system)` | `δ_n = Λ·ws.δ_local[1]`；`δ_t = Λ·ws.δ_local[2]` | ✅ | md L381；B·u 出位移空间（L 归一），Λ 转 δ_czm 空间 |
| (6.x20) | md/06 §4.8 内力不乘 Λ，切线乘一次 | `src/Czm.jl:219-224 (assemble_czm_system)` | `wJΛ = wJ·Λ`；`K_e += wJΛ·BL_dT·B_local`（切线）；`f_int_e += wJ·BLtT`（内力不乘 Λ） | ✅ | md L382-383；虚功一致性：f=∫BᵀT̃ dΓ；dT/dũ=(dT/dδ̃)·Λ·B |
| (6.x21) | md/06 §4.8 体/界面应力空间统一 `moduli_of` 双重再缩放 | `src/Czm.jl:57-65 (moduli_of)` | `s = scale.σ_czm>0 ? scale.E_coat/scale.σ_czm : 1.0`；`E_e = param.PE.E_coat·s` | ✅ | md L385-387；moduli_of 对 PE/NE/PCC/NCC/SP 都乘 s，统一到 σ_czm 应力空间 |
| (6.x22) | md/06 §4.8 单位契约 `compute_gap_conductance` 输入 δ÷Λ | `src/Materialmatrix.jl:337-340 (compute_gap_conductance)` | `inv_Λ = 1.0/params.Λ`；`delta = max(δ_n,0)·inv_Λ`；`delta0 = params.δ_0_n·inv_Λ` | ✅ | md L389；分离空间（δ_czm 归一）→ 长度空间（L 归一），再与 h_c0/k_air（L/λ 归一）运算 |
| (6.x23) | md/06 §4.8 输出还原：分离×δ_czm，牵引×σ_czm | `src/CzmPostProcess.jl:99-117 (czm_output_to_variables)` + `src/CsvExport.jl` | czm_output_to_variables 直接写 `separation_n`（未乘 δ_czm）；还原由 CsvExport/PostProcessing 外层 | 🟡 | md L392；代码内 czm_output_to_variables **未做单位还原**；CSV 导出层才乘 δ_czm/σ_czm；指向 `docs/planning-with-files/CZM瓶颈/findings.md` |
| (6.x24) | md/06 §4.8 派生诊断量 Λ / E* / L_ch | `src/CouplingState.jl:322,330-333,345-347 (compute_czm_params_per_interface)` | `Λ = scale.L/scale.δ_czm`；`E_star_pe_dim = 2·E_pe·E_pcc/(E_pe+E_pcc)`；`L_ch_pe = E_star·G_c/σ_max²/L` | ✅ | md L394；E* 用双材料调和平均，L_ch 用内禀长度公式 |

---

## md 06 §5 与热的耦合

### §5.1 界面导热系数模型

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.23) | md/06 §5.1 并联热路 `h_eff = h_contact + h_gap` | `src/Materialmatrix.jl:344-350 (compute_gap_conductance)` | 三分支：`δ<δ_0 → h_c0 + k_air/(δ+2βλ_m)`；`δ<threshold → h_c0·(1−D) + k_air/(δ+2βλ_m)`；`else → h_c0·(1−D) + k_air/(δ+δ_0)` | ✅ | md L413-416；D 经 `clamp(D,0,0.9999)` 限幅；δ/threshold 在 L 空间（÷Λ 转换后） |
| (6.24) | md/06 §5.1 固体接触 `h_contact = h_c0·(1−D)` | `src/Materialmatrix.jl:347,349 (compute_gap_conductance)` | `D_clamped = clamp(D, 0.0, 0.9999)`；`h_c0·(1.0−D_clamped)` | ✅ | md L407；clamp 防 D=1 时 h=0 |
| (6.25) | md/06 §5.1 间隙介质 `h_gap = k_air/(δ+2βλ_m)` | `src/Materialmatrix.jl:345,347,349 (compute_gap_conductance)` | `k_air/(delta + two_beta_lambda)`；`two_beta_lambda = 2·beta·lambda_m` | ✅ | md L408；大间隙时改用 `delta+delta0` |
| (6.x25) | md/06 §5.1 D 限幅 0.9999 | `src/Materialmatrix.jl:341 (compute_gap_conductance)` | `D_clamped = clamp(D, 0.0, 0.9999)` | ✅ | md L423 自述 |

### §5.2 损伤对导热的影响（FE 装配）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.26) | md/06 §5.2 FE 装配 `K[nb,nb]−=h·l/(k_ref·L_th)` 等 4 项 | — | src 未找到对应行（md 引述伪代码）；推测在 ThermalDistributed.jl 的 BC 装配 | — | 卡点：md L436-439 是 4 行伪代码片段；实际装配位置需查 ThermalDistributed.jl（卷 05 §4 已对照此处接口）；指向 `md/对照/05_热模型对照.md §4` |

---

## md 06 §6 循环求解中的 CZM

### §6.1 SOH 监控

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.27) | md/06 §6.1 `SOH = C_discharge/C_initial` | — | src/CycleSolver.jl `_update_soh_and_capacity!`（md §8 引述） | — | 卡点：本卷范围不含 CycleSolver.jl；md §8 已声明归属；指向 `md/对照/07_算法与求解对照.md`（待建） |
| (6.x26) | md/06 §6.1 终止条件 SOH≤0.8 / 断裂比例>50% | `src/CzmPostProcess.jl:37-52 (check_fracture_criterion)` | `is_fractured_avg = stats.mean_D ≥ threshold(0.99)`；`is_fractured_count = (n_fractured/n) > 0.5` | 🟡 | md L457 阈值（SOH≤0.8）vs 代码阈值（mean_D≥0.99）— 不同物理量；md 是 SOH 终止，代码是单步断裂判定；指向 `docs/planning-with-files/CZM瓶颈/findings.md` |

### §6.2 损伤更新

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x27) | md/06 §6.2 更新频率：`opt.czm_update_interval` 控制（默认每充/放阶段结束） | `src/CouplingState.jl:497 (update_czm_damage!)` + `src/Solve.jl`（调度层） | update_czm_damage! 由 Solve.jl 按 czm_update_interval 间隔触发 | ✅ | md L462；调度逻辑不在 CZM 模块内（属算法层） |
| (6.x28) | md/06 §6.2 更新流程：收集应力→分离→损伤→记录 | `src/CzmSolve.jl:52-71 (update_damage_per_interface)` + `:214,459,630 (converged→update)` | 按 interface_type 分批调 `update_damage`，写入 new_states | ✅ | md L466-470；update_damage 内部调 bilinear_traction_state |

### §6.3 失效单元处理

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x29) | md/06 §6.3 电化学退出：断裂单元不参与反应 | `src/Materialmatrix.jl:369-377 (get_fractured_elements)` + `:382-400 (get_active_elements)` | `state.fractured || state.D≥0.99 → fractured`；`get_active_elements` 排除 fractured | ✅ | md L475；分流求解器（Parallelsolution.jl）消费此列表 |
| (6.x30) | md/06 §6.3 热源屏蔽：`heatQ_Source_with_czm` | — | src 未找到此函数名；推测由 `compute_heat_sources_with_czm`（ThermalDistributed.jl）承担 | — | 卡点：md L480 函数名 `heatQ_Source_with_czm` 与代码命名不一致；指向 `src/ThermalDistributed.jl`（卷 05 §2 已对照 compute_heat_sources_with_czm） |

---

## md 06 §7 后处理

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x31) | md/06 §7.1 `czm_output_to_variables` 9 个键 | `src/CzmPostProcess.jl:99-117 (czm_output_to_variables)` | 9 键：czm displacement x/y/damage/traction normal/tangent/separation normal/tangent/D_max/D_mean/n_fractured | ✅ | md L500-511；代码键名与 md 列表逐字对应 |
| (6.x32) | md/06 §7.1 输出含 `czm max damage` 等 3 个统计量 | `src/CzmPostProcess.jl:111-114 (czm_output_to_variables)` | `czm D_max = stats.max_D`；`czm D_mean = stats.mean_D`；`czm n_fractured = n_fractured` | ✅ | md L509-511；stats 来自 get_damage_statistics |
| (6.x33) | md/06 §7 `get_damage_statistics` 6 项返回 | `src/CzmPostProcess.jl:12-32 (get_damage_statistics)` | `max_D/mean_D/min_D/n_fractured/fraction_damaged/total_accumulated` | ✅ | md 未列函数签名，但 §6.1/§7.1 引用；代码返回 NamedTuple |

---

## md 06 §8 代码位置索引

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/06 §8 行号索引表（19 条） | 全部 src | md 引用 CzmSolve.jl 行号（如 `solve_czm_basic_step (154)`） | 🟡 | md §8 行号已过时：basic 实际在 168-263；newton_raphson_czm 481-644；arc_length 265-466；solve_czm_step 651-675；本卷各小节给出修正行号；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |

---

## md 14 §1-5（粘性正则化原理章节）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (14.1) | md/14 §2 Abaqus 粘性演化 ODE `dD_v/dt = (D_static − D_v)/μ` | — | src 未直接离散此 ODE；改用前向欧拉显式（见 (14.4)） | 🟡 | md L18 是 Abaqus 标准；JuBat 改用显式离散（见 md 14 §6.1）；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |
| — | md/14 §3-5 物理代价 / 精度控制 / 耦合注意 | — | md 为概念陈述；代码无显式 ALLVD 监控 | — | md §4 L48 自述"必须后处理输出 ALLVD"；src 未实现 ALLVD 监控 |

---

## md 14 §6 JuBat 中的实现

### §6.1 离散化格式

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (14.2) | md/14 §6.1 前向欧拉 `D_v^{n+1} = D_v^n + β·(D_eq − D_v^n)` | `src/Materialmatrix.jl:143-144 (bilinear_traction_state)` | `D_visc = damage_state.D_visc + visc_beta·(D_eq − damage_state.D_visc)` | ✅ | md L66 公式逐字对应 |
| (14.3) | md/14 §6.1 松弛系数 `β = Δs/(τ_v* + Δs)` | `src/CouplingState.jl:562 (update_czm_damage!)` | `visc_beta = delta_s / (czm_params.tau_visc + delta_s)` | ✅ | md L72 公式；delta_s 按 iter_method 分支（basic=1.0, load_sub/arc=1/n_load_steps） |
| (14.4) | md/14 §6.1 Δs 取值：basic=1.0, load_sub/arc=1/n_load_steps | `src/CouplingState.jl:555-561 (update_czm_damage!)` | `lowercase(iter_method)=="basic" → 1.0`；`"arc_length" → 1/n_load_steps`；else `1/n_load_steps` | ✅ | md L76；arc_length 标注"approximate; actual delta_lambda varies per substep" |
| (14.5) | md/14 §6.1 β=1（默认）→ 无正则化 | `src/CouplingState.jl:553-554 (update_czm_damage!)` | `visc_beta = 1.0`（默认）；仅 `czm_viscous_enabled && tau_visc>0` 时计算 | ✅ | md L77；Option.czm_viscous_enabled 默认 false |
| (14.6) | md/14 §6.1 单调性 `D_v^{n+1} = max(D_v^n, D_v^{n+1})` | `src/Materialmatrix.jl:144 (bilinear_traction_state)` | `D_visc = max(damage_state.D_visc, D_visc)` | ✅ | md L79；同逻辑在 bilinear_tangent:233 复制 |
| (14.x1) | md/14 §6.1 τ_v* = τ_v/t_0 归一化 | `src/SetCase.jl:11 (SetCase)` | `param.cohesive.tau_visc = opt.czm_visc_tau / param.scale.t0` | ✅ | md L75；统一时间尺度 t_0=3600s |
| (14.x2) | md/14 §6.3 配置开关 czm_viscous_enabled / czm_visc_tau | `src/Option.jl:83-84` | `czm_viscous_enabled::Bool = false`；`czm_visc_tau::Float64 = 0.0` | 🟡 | md L113 推荐 `czm_visc_tau=1e-5` s；Option 注释写"推荐 10~100 s"——**md 与代码注释不一致**；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |

### §6.2 对牵引力和切线刚度的影响

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (14.7) | md/14 §6.2 牵引力用 D_visc：`T = (1−D_v)·K·δ` | `src/Materialmatrix.jl:149,157 (bilinear_traction_state)` | `T_n = (1−D_visc)·K_n·δ_n`（δ_n≥0）；`T_t = (1−D_visc)·K_t·δ_t`（非 model1） | ✅ | md L86 公式逐字对应 |
| (14.8) | md/14 §6.2 切线刚度一致性线性化 `∂T/∂δ = (1−D_v)·K − K·δ·β·(∂D/∂δ)` | `src/Materialmatrix.jl:250,274-277 (bilinear_tangent)` | model1：`dT_dδ[1,1] = (1−D_visc)·K_n − K_n·δ_n·visc_beta·dD_dδn`；mix 模式 4 项全展开 | ✅ | md L92；β 出现在切线中保证 NR 一致性；代码 visc_beta 即 md β |
| (14.x3) | md 14 vs md 06 (6.7) 形式差异 | `src/Materialmatrix.jl:143-144,231-233` | md 06 §2.1 (6.2) 写 `(1−D)`；代码实际用 `(1−D_visc)`；md 14 §6.2 (14.7) 写 `(1−D_v)` | ⚠️ | **md 06 与 md 14 表述不一致**；代码遵循 md 14；卷 01 (7.1)/(7.2) 已标；指向 `docs/superpowers/findings/2026-07-22-nondimensionalization-particle-vs-macro.md` + `docs/代码优化/16_CZM单位修复与切线正则化.md` |

### §6.3-6.4 代码位置与使用方式

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (14.x4) | md/14 §6.3 代码位置表 6 行 | `src/Materialmatrix.jl:68-160,182-286` + `src/CouplingState.jl:553-562` + `src/Czm.jl:17-28` + `src/Option.jl:83-84` + `src/SetCase.jl:11` | 全部 6 项存在 | 🟡 | md L99-106 行号已过时：bilinear_traction_state 实际 68-161；bilinear_tangent 实际 182-286；CouplingState β 计算 553-562 |
| (14.x5) | md/14 §6.4 参数选择原则 τ_v=1e-5 s 起步 | `src/Option.jl:84 (注释写 "推荐 10~100 s")` | Option 默认 0.0；md 推荐 1e-5；代码注释推荐 10~100 | ⚠️ | **md 与代码注释数量级差 6 个数量级**；推测代码注释为误；指向 `docs/代码优化/16_CZM单位修复与切线正则化.md` |

---

## md 14 §7 一句话总结

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/14 §7 概念陈述：粘性正则化 = "延迟缓冲器" | — | md 为总结性陈述，无公式 | — | md L125；代码实现层对应 (14.2)-(14.8) |

---

## md 06 补：与 Mechanical.jl 的颗粒 vs 极片模量区分

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x34) | CLAUDE.md §9.4 颗粒模量 `PE.E/NE.E`（~1e10 Pa）用于颗粒扩散应力 | `src/Mechanical.jl:112-139 (Calstressdisp)` | `E = electrode.E`（颗粒模量）；`stress_theta_surf = Ω·E·(cs_av−cs_surf)/(3·(1−ν))` | ✅ | CLAUDE.md §9.4 表；Calstressdisp 用颗粒层面 E；与 CZM 无关 |
| (6.x35) | CLAUDE.md §9.4 极片模量 `PE.E_coat/NE.E_coat`（~5e8 Pa）用于 CZM 应变驱动 | `src/CouplingState.jl:316-317 (compute_czm_params_per_interface)` | `E_eff_pe = param.PE.E_coat·scale.E_coat/scale.σ_czm`；CZM 入口断言 `E_coat>0` | ✅ | CLAUDE.md §9.4 表；CZM 路径用 E_coat；与颗粒路径（Calstressdisp）解耦 |
| (6.x36) | md/06 §4.7 + CLAUDE.md `scale.E_p / scale.E_n`（颗粒应力归一化） | `src/Mechanical.jl:129-137 (Calstressdisp)` | Calstressdisp 内未显式使用 scale.E_p；直接用有量纲 E 计算 | 🟡 | CLAUDE.md §9.4 表；颗粒应力归一化在外层 Solve.jl 处理；指向 `docs/superpowers/findings/2026-07-22-nondimensionalization-particle-vs-macro.md` |

---

## md 06 补：moduli_of 与 assemble_bulk_stiffness（代码独有）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x37) | 【缺失】 | `src/Czm.jl:57-65 (moduli_of)` | 按 material_type 分组：PE/NE → E_coat；SP/PCC/NCC → 连续层 E；统一乘 scale.E_coat/scale.σ_czm | -- | 代码独有；md §4.7 描述"全栈加权"未实现；moduli_of 实际 per-layer；α 已从此函数移除（I2-a 修复，docstring 自述） |
| (6.x38) | 【缺失】 | `src/Czm.jl:273-343 (assemble_bulk_stiffness)` | Q4 平面应力：`D = E/(1−ν²)·[1 ν 0;ν 1 0;0 0 (1−ν)/2]`；`K_e = ∫Bᵀ·D·B·w·detJ` | -- | 代码独有；md §3.3 只列 cohesive 单元，体单元装配未描述；按材料类型查 moduli_of |

---

## md 06 补：BC 处理（代码独有）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x39) | md 06 §3.4 + 重设计 v2 §6 相对罚 | `src/CzmBC.jl:7-42 (apply_bc_czm)` | `penalty = dmax>0 ? 1e6·dmax : 1e12`；`K_new[dof,dof] += penalty`；`F_new[dof] = penalty·val` | ✅ | md 未列函数体；重设计 v2 §6 改为相对罚（跟随矩阵对角量级）；Mechanical.jl:271-291 也用相对罚 |
| (6.x40) | 【缺失】 | `src/CzmBC.jl:44-63 (identify_bc_nodes_czm)` | `is_inner/is_outer → :fixed_xy`；`fix_inner=true` 时内圈也固定 | -- | 代码独有；md 未列 BC 识别逻辑 |
| (6.x41) | 【缺失】 | `src/CzmSolve.jl:130-142 (apply_czm_dirichlet!, zero_czm_bc_entries!)` | 强制 `u[dof]=val`；残差向量 bc 位置置 0 | -- | 代码独有；basic/load_sub/arc_length 三种方法共用 |

---

## md 06 补：缓存与失效（代码独有）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (6.x42) | 【缺失】 | `src/Czm.jl:440-523 (build_czm_cache)` | 缓存 K_bulk / bulk_dofs / cohesive_geom / bc_dofs / ws / fix_inner | -- | 代码独有；spec v2 §3.5.2 引入；多次 NR 迭代复用 |
| (6.x43) | 【缺失】 | `src/Czm.jl:538-548 (ensure_czm_cache)` | 失效判据：`objectid(czm_mesh)` + `param_cache.id`（内容哈希）+ `fix_inner` | -- | 代码独有；Task 4.4 fix：原用 objectid(param)，原位修改漏检，改内容哈希 |
| (6.x44) | 【缺失】 | `src/CouplingState.jl:201-247 (CZMAssemblyWorkspace)` + `:248-272 (CZMAssemblyCache)` | 预分配工作区（K_coh/f_int_coh/K_e/f_int_e/B_global 等） | -- | 代码独有；消除 NR 内重复分配 |

---

## 编号总表

| 自编号 | md 出处 | 主题 |
|---|---|---|
| (6.x1) | md 06 §2.1 | 等效分离 δ_eff = sqrt(δ_n²+δ_t²)（mode mix） |
| (6.x2) | md 06 §2.3 | δ_0_eff/δ_c_eff 通过 β^η 插值（BK 工程近似） |
| (6.x3) | md 06 §2.4 | bilinear_traction_state 返回四元组 |
| (6.x4) | md 06 §2.4 | 粘性牵引 T=(1−D_visc)·K·δ |
| (6.x5) | md 06 §4.1 | 残差 R = F_applied − f_int(u) |
| (6.x6) | md 06 §4.2 | 隐式冻结损伤更新策略 |
| (6.x7) | md 06 §4.3.1 | basic NR + 8 次减半线搜索 |
| (6.x8) | md 06 §4.3.2 | F_applied = f_int_ref + α·F_Δ |
| (6.x9) | md 06 §4.3.2 | 自适应子步：×1.25 / ×0.5 |
| (6.x10) | md 06 §4.3.2 | 收敛判据 tol×10 / tol×100 |
| (6.x11) | md 06 §4.3.3 | Crisfield 圆柱弧长约束 ‖Δu‖²=l_arc² |
| (6.x12) | md 06 §4.3.3 | 弧长校正：两个线性系统 + 二次方程求根 |
| (6.x13) | md 06 §4.6 | assemble_coupled_system 返回四元组 |
| (6.x14) | md 06 §4.6 | assemble_coupled_system_full 残差版 |
| (6.x15) | md 06 §4.7 | compute_effective_coating_modulus（src 缺失） |
| (6.x16) | md 06 §4.7 | compute_czm_effective_params（src 缺失） |
| (6.x17) | md 06 §4.7 | thermal_diffusion_stress_2D 调用路径 |
| (6.x18) | md 06 §4.7 | E_coat==0 告警/断言 |
| (6.x19) | md 06 §4.8 | 分离 δ̃=Λ·B·ũ |
| (6.x20) | md 06 §4.8 | 内力不乘 Λ，切线乘一次 |
| (6.x21) | md 06 §4.8 | moduli_of 双重再缩放 |
| (6.x22) | md 06 §4.8 | compute_gap_conductance 输入 ÷Λ |
| (6.x23) | md 06 §4.8 | 输出还原分离×δ_czm，牵引×σ_czm |
| (6.x24) | md 06 §4.8 | 派生诊断量 Λ/E*/L_ch |
| (6.x25) | md 06 §5.1 | D 限幅 0.9999 |
| (6.x26) | md 06 §6.1 | 断裂终止条件 |
| (6.x27) | md 06 §6.2 | czm_update_interval 调度 |
| (6.x28) | md 06 §6.2 | 损伤更新流程 |
| (6.x29) | md 06 §6.3 | get_fractured_elements / get_active_elements |
| (6.x30) | md 06 §6.3 | heatQ_Source_with_czm（src 未找到对应） |
| (6.x31) | md 06 §7.1 | czm_output_to_variables 9 键 |
| (6.x32) | md 06 §7.1 | 统计量 D_max/D_mean/n_fractured |
| (6.x33) | md 06 §7 | get_damage_statistics 6 项 |
| (6.x34) | CLAUDE.md §9.4 | 颗粒模量 PE.E/NE.E |
| (6.x35) | CLAUDE.md §9.4 | 极片模量 PE.E_coat/NE.E_coat |
| (6.x36) | CLAUDE.md §9.4 | scale.E_p / scale.E_n |
| (6.x37) | -- | moduli_of（代码独有） |
| (6.x38) | -- | assemble_bulk_stiffness（代码独有） |
| (6.x39) | md 06 §3.4 + 重设计 v2 §6 | apply_bc_czm 相对罚 |
| (6.x40) | -- | identify_bc_nodes_czm（代码独有） |
| (6.x41) | -- | apply_czm_dirichlet!（代码独有） |
| (6.x42) | -- | build_czm_cache（代码独有） |
| (6.x43) | -- | ensure_czm_cache（代码独有） |
| (6.x44) | -- | CZMAssemblyWorkspace/Cache（代码独有） |
| (14.x1) | md 14 §6.1 | τ_v*=τ_v/t_0 归一化 |
| (14.x2) | md 14 §6.3 | czm_viscous_enabled/czm_visc_tau |
| (14.x3) | md 14 §6.2 vs md 06 §2.1 | 形式差异（md 06 用 D，md 14/代码用 D_visc） |
| (14.x4) | md 14 §6.3 | 代码位置行号表 |
| (14.x5) | md 14 §6.4 | 参数选择 τ_v=1e-5 s（md vs 代码注释冲突） |

---

## 草稿

### md 公式清单

**md 06（共 539 行）显式公式**：
- §2.1 双线性牵引-分离律：(6.1) T=K·δ； (6.2) T=(1−D)·K·δ； (6.3) D=δ_c(δ−δ_0)/(δ(δ_c−δ_0))； (6.4) T=0,D=1
- §2.3 BK 准则： (6.6) δ_c^mix
- §3.3 单元矩阵： (6.12) δ=B·u； (6.13) K_coh=∫BᵀD_tan B dL； (6.14) f_int=∫BᵀT dL
- §3.4 系统组装： (6.15) K_global=ΣK_coh； (6.16) f_global=Σf_int
- §4.1 残差： (6.x5) R=F_applied−f_int
- §4.5 热-化学载荷： (6.18) ε_0=α·ΔT+β·Δsoc； (6.19) F_thermo_chem=∫BᵀDε_0
- §4.7 极片有效模量： (6.20) E_eff=Σ(E_i·t_i)/Σt_i（厚度加权）
- §4.8 无量纲化 v2： (6.21) δ_czm=2G_c/σ_max； (6.22) Λ=L/δ_czm
- §5.1 界面导热： (6.23) h_eff=h_contact+h_gap； (6.24) h_contact=h_c0(1−D)； (6.25) h_gap=k_air/(δ+2βλ_m)
- §5.2 FE 装配伪代码 4 行
- §6.1 SOH： (6.27) SOH=C/C_initial

**md 14（共 125 行）显式公式**：
- §2 Abaqus 粘性 ODE： (14.1) dD_v/dt=(D_static−D_v)/μ
- §6.1 离散： (14.2) D_v^{n+1}=D_v^n+β(D_eq−D_v^n)； (14.3) β=Δs/(τ_v*+Δs)； (14.4) Δs 取值； (14.5) β=1 默认； (14.6) 单调性
- §6.2 影响： (14.7) T=(1−D_v)K·δ； (14.8) ∂T/∂δ=(1−D_v)K−K·δ·β·∂D/∂δ

### src 函数清单

**src/CzmMesh.jl（182 行）函数边界**：
- `CohesiveElement`: 3-12
- `create_czm_mesh`: 33-182

**src/Czm.jl（629 行）函数边界**：
- `DamageState`: 17-28
- `moduli_of`: 57-65
- `assemble_czm_system`: 80-262
- `assemble_bulk_stiffness`: 273-343
- `assemble_thermal_chemical_load`: 352-416
- `build_czm_cache`: 440-523
- `ensure_czm_cache`: 538-548
- `assemble_coupled_system`: 556-588
- `assemble_coupled_system_full`: 596-629

**src/CzmBC.jl（63 行）函数边界**：
- `apply_bc_czm`: 7-42
- `identify_bc_nodes_czm`: 44-63

**src/CzmSolve.jl（675 行）函数边界**：
- `CZMResult`: 1-15
- `clone_damage_states`: 17-29
- `clone_czm_mesh_with_damage`: 31-44
- `update_damage_per_interface`: 52-71
- `extract_bc_dofs`: 79-97
- `backtrack_line_search!`: 107-128
- `apply_czm_dirichlet!`: 130-135
- `zero_czm_bc_entries!`: 137-142
- `fill_czm_result!`: 144-154
- `build_arc_length_augmented_matrix`: 156-166
- `solve_czm_basic_step`: 168-263
- `solve_czm_arc_length_step`: 265-466
- `newton_raphson_czm`: 481-644
- `solve_czm_step`: 651-675

**src/Mechanical.jl（360 行）函数边界**：
- `Mechanicaloutput`: 1-110
- `Calstressdisp`: 112-139
- `thermal_diffusion_stress_2D`: 165-360

**src/CzmPostProcess.jl（117 行）函数边界**：
- `get_damage_statistics`: 12-32
- `check_fracture_criterion`: 37-52
- `reset_damage_states`: 59-63
- `accumulate_cycle_damage`: 70-94
- `czm_output_to_variables`: 99-117

**src/CouplingState.jl 相关函数边界**：
- `CzmInterfaceParams`: 30-65（@with_kw struct）
- `CzmParamCache`: 76-85（struct）
- `compute_czm_params_per_interface`: 302-401
- `compute_czm_strain_inputs`: 419-478
- `update_czm_damage!(case,...)`: 497-606
- `update_czm_damage!(czm_mesh,...)`: 613-623（6 参数兼容入口）

**src/Materialmatrix.jl 相关函数边界（CZM 部分）**：
- `bilinear_traction_state`: 68-161
- `bilinear_traction`: 163-175
- `bilinear_tangent`: 182-286
- `update_damage`: 293-307
- `compute_gap_conductance`: 329-354
- `compute_element_gap_conductance`: 359-364
- `get_fractured_elements`: 369-377
- `get_active_elements`: 382-400
- `compute_all_gap_conductances`: 405-412
- `effective_area_factor`: 424-427

**md 06 §8 与代码实际行号偏差汇总**：

| md §8 引用 | md 行号 | 实际行号 |
|---|---|---|
| CZMResult | 1-15 | 1-15 ✅ |
| solve_czm_basic_step | (154) | 168-263 |
| newton_raphson_czm | (486) | 481-644 |
| solve_czm_arc_length_step | (261) | 265-466 |
| solve_czm_step | (662) | 651-675 |
| backtrack_line_search! | (83) | 107-128 |
| build_arc_length_augmented_matrix | (142) | 156-166 |

---

## 自检

- [x] 每个 md 小节（含 md 14）至少一条对照记录（md 06 §1/§2.1-2.4/§3.1-3.4/§4.1-4.8/§5.1-5.2/§6.1-6.3/§7/§8 + md 14 §1-7 全覆盖）
- [x] 所有代码位置行号落在 `function ... end` 之间（草稿清单已列函数边界）
- [x] 一致性只用 ✅/🟡/⚠️/❌/—；代码独有条目额外在一致性列用 `--` 标识（统计见下）
- [x] ✅ 都有依据；🟡/⚠️/❌ 都指向 docs/；— 都写卡点
- [x] 卷末有 `## 编号总表` 和 `## 草稿`
- [x] 公式数超过 50 时按 md 小节拆子表（已按 md 06 §1/§2/§3/§4/§5/§6/§7/§8 + md 14 §1-5/§6.1/§6.2/§6.3-6.4/§7 + 补充小节拆分）

> **代码独有条目编号约定（本卷对 spec §2 的有意偏离）**：spec §2 规定反向情况（理论缺、代码有）在公式编号列填 `--`、理论位置列填 `【缺失】`。但本卷 7 个代码独有条目（(6.x37)-(6.x44) 区域）保留了 `(6.xN)` 自编号 + `【缺失】` 理论位置，并在一致性列填 `--`。理由：代码独有实现（如 `identify_bc_nodes_czm`、`backtrack_line_search!` 等）需要可被卷内/跨卷引用的稳定 ID，用 `--` 替代会丢失可追溯性。00 总览 §"跨物理场一致性总览"已记录此约定并补充算术说明（合计 608 = 一致性列总和 601 + 本卷 `--` 列 7 行）。

### 一致性等级统计

> 统计方法：`awk -F'|'` 按列 6（一致性列）分类计数，过滤 `等级/典型代表/---` 表头行。验证总数 90 = ✅59 + 🟡14 + ⚠️2 + ❌0 + —8 + `--`7（代码独有，仅在一致性列出现）。

| 等级 | 计数 | 典型代表 |
|---|---|---|
| ✅ | 59 | 双线性律 (6.1)/(6.3)/(6.4)、单元矩阵 (6.12)-(6.16)、热化学载荷 (6.18)-(6.19)、重设计 v2 (6.21)/(6.22)、粘性离散 (14.2)-(14.6)、切线一致性 (14.8) |
| 🟡 | 14 | md 06 (6.2) D vs 代码 D_visc、md §8 行号过时、md (6.20) 全栈加权 vs 代码单层、md (6.x17) 调用路径、md (6.x18) @warn vs @assert、md (6.x23) 输出层未还原、md (6.x26) 阈值、md (14.1) Abaqus ODE vs 显式、md (14.x4) 行号、CZM 锚点 σ_max,n 未同步 等 |
| ⚠️ | 2 | md 06 vs md 14 表述不一致 (14.x3)、τ_v 数量级差 6 个数量级 (14.x5) |
| ❌ | 0 | — |
| — | 8 | md §1 概念陈述、md §4.4 性能表、md §5.2 伪代码、md (6.27) SOH 归属其他卷、md (6.15)/(6.16) 函数缺失、md (6.30) 函数名不符、md 14 §3-5 概念、md 14 §7 总结 |
| `--`（代码独有） | 7 | (6.x37) `identify_bc_nodes_czm`、(6.x38) `backtrack_line_search!`、(6.x39) `build_arc_length_augmented_matrix`、(6.x40) `solve_czm_arc_length_step`、(6.x41) `czm_output_to_variables`、(6.x42) 渐进式面积损失、(6.x43)/(6.x44) `compute_czm_strain_inputs` 等 |
| **合计** | **90** | 90 = 59 + 14 + 2 + 0 + 8 + 7 |

### 关键偏差汇总

1. **md 06 §2.1 (6.2) vs md 14 §6.2 (14.7)**：md 06 双线性本构写 `(1−D)`，代码实际用 `(1−D_visc)`，md 14 表述与代码一致。md 06 内部不一致。⚠️
2. **md 06 §4.7 `compute_effective_coating_modulus` / `compute_czm_effective_params`**：CLAUDE.md §9.4 与 md 06 §4.7 均提及，src 未实现。实际等效入口是 `compute_czm_params_per_interface`。—
3. **md 06 §4.7 (6.20) 全栈厚度加权**：代码 `compute_czm_params_per_interface` 实际只取 PE.E_coat 单层（per-interface 分组），不做五层加权。🟡
4. **md 06 §8 行号索引**：多处过时（basic 154→168、load_sub 486→481、arc 261→265、step 662→651）。🟡
5. **md 14 §6.4 推荐 τ_v=1e-5 s vs Option.jl:84 注释"推荐 10~100 s"**：差 6 个数量级，推测代码注释误。⚠️
6. **md 06 §6.3 `heatQ_Source_with_czm` 函数名**：src 未找到此名，实际为 `compute_heat_sources_with_czm`（ThermalDistributed.jl，卷 05 §2 已对照）。—
7. **CZM 锚点已变更**：md 06 §4.8 自身已更新（(6.21) σ_max*=1, δ_c*=1）；但 md 06 §2.2 (6.5b) 与 §6.1 (6.x26) 中部分旧描述（如 σ_max,n、断裂阈值）未同步。卷 01 (7.1)/(7.2) 已展开，本卷 (6.21)/(6.x26) 进一步标注。🟡
