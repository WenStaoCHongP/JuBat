# 堆芯塌陷力学建模 Batch 3（PCC/NCC J2 塑性 / C2）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现 `czm_j2_plasticity=true` 的 PCC/NCC 平面应力一致 J2 返回映射（Simo–Hughes），接入 GL 单元与 `mech_state` 槽位（Batch 1 预留位的首个消费者），塑性状态跨时间步持久化且收敛才提交；默认关零漂移；通过 spec §7 Batch 3 测试与 C2/L1 门。

**Architecture:** 三任务。Task 1 新建 `src/CzmPlasticity.jl`（`PlasticState` 类型、`foil_params_of` 物理箔参数、`return_mapping_plane_stress` 4×4 Newton 返回映射 + 一致切线，切线由收敛系统雅可比线性化精确导出）+ 单元测试；Task 2 把塑性接入 `gl_element_residual_tangent`/`assemble_bulk_residual_tangent(mech_state, plasticity=true)`/求解器链路/`update_czm_damage!`（`CzmLayout.plastic_states` 持久化 + 收敛后提交）+ 集成测试；Task 3 三道门禁与记录。

**Tech Stack:** Julia 1.11.2；`Test`/`LinearAlgebra`/`SparseArrays`。

**Spec:** `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md`（v1.3）§3.3、§4.1/§4.2、§7 Batch 3 行；参数依据：`docs/planning-with-files/30_堆芯塌陷力学建模/宏观电池模型_正负极极片力学参数_文献整理.md` §10.7（联网核实记录见 findings"Batch 3 参数核实"节）。

## Global Constraints

（同 Batch 2 计划：运行环境、强制行为基线 PNG SHA `4ba6207c…e932`、兼容性契约、错误处理 AGENTS 9.7、每 Task 一提交。）

### 三项设计决策（执行契约）

**D-B3-0（物理箔本构切换，opt-in）**：现状 PCC/NCC 用 500 MPa 软等效层（`src/parameters/Jellyroll.jl:99,109`），若直接叠加 σ_y=60/200 MPa 屈服应变达 ~12%，塑性永不激活。故 `czm_j2_plasticity=true` 时 PCC/NCC 切换为物理箔参数集（文献整理 §10.7，联网核实区间内）：Al（PCC）`E=70 GPa, ν=0.33, σ_y=60 MPa`；Cu（NCC）`E=110 GPa, ν=0.34, σ_y=200 MPa`；`H=0`（理想塑性，敏感性扫描 0–2 GPa）。新字段 `PCC/NCC.E_foil/nu_foil/sigma_y/H` 带默认值加入 `CurrentCollector`（SetParams.jl:86）；默认关时 PCC/NCC 仍走 500 MPa 等效层路径，零漂移。
**D-B3-1（塑性需 geo_nl）**：`plasticity=true && geo_nl=false` → `error`（线弹性常刚度路径无逐 GP 应力评估，linear+plastic 组合无消费者；工况 C 塑性与几何非线性同开）。
**D-B3-2（收敛才提交）**：装配对已提交塑性状态做纯函数试算（不原地变异）；求解收敛后由 `update_czm_damage!` 以收敛位移调用一次 `commit_plastic!` 写回 `CzmLayout.plastic_states`——满足 spec"trial 非收敛不提交"的原子性（无需求解内 clone/回滚）。

### 归一化契约

新应力类字段一律 `/σ_czm` 归一（`NormaliseParam` 追加四行；σ_czm=8.2e7 Pa）；塑性路径内部不再乘 `moduli_of` 的 `E_coat/σ_czm` 链（`foil_params_of` 直接返回 σ_czm 单位的 `(E, ν, σ_y, H)`）。`mech_state` = `Matrix{PlasticState}`（`[ne, 4]`，2×2 高斯点；仅 PCC/NCC 行被消费）。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/CzmPlasticity.jl` | 新建 | `PlasticState`、`foil_params_of`、`return_mapping_plane_stress`、`commit_plastic!`、`clone_plastic_states` |
| `src/JuBat.jl` | 修改 | include 新文件 |
| `src/SetParams.jl` | 修改 | `CurrentCollector` 四新字段；`NormaliseParam` 归一化 |
| `src/czm.jl` | 修改 | `gl_element_residual_tangent` 塑性分支；`assemble_bulk_residual_tangent` 消费 `mech_state`/`plasticity`；`assemble_coupled_system` 透传 |
| `src/CzmSolve.jl` | 修改 | `plasticity`/`mech_state` 关键字贯穿（与 geo_nl 同模式） |
| `src/CouplingState.jl` | 修改 | `CzmLayout.plastic_states`；`update_czm_damage!` 接线 + 提交 + σ_y 断言 |
| `test/unit_czm_j2.jl` | 新建 | 本构级测试（spec §7 Batch 3：单轴屈服/卸载、KKT、耗散非负、一致切线 FD） |
| `test/test_czm_j2_integration.jl` | 新建 | 集成测试（C2-lite：强本征应变下 PCC/NCC 屈服、塑性关退回 C1 逐位、耗散非负） |

---

## Task 1: 本构层——类型、箔参数、平面应力 J2 返回映射（TDD）

`return_mapping_plane_stress(e_mech, C, σ_y, H, eps_p, κ) -> (σ, C_ep, Δeps_p, Δκ)`：
未知量 `(σ11,σ22,σ12,Δγ)` 解 `σ = C(e_mech − eps_p − Δγ n(σ))` 与 `f = σ̄ − σ_y − H(κ+Δγ) = 0`（σ̄² = σ11²+σ22²−σ11σ22+3σ12²），4×4 Newton（弹性试算起步，tol 1e-12，≤50 步）；弹性步（f≤0）直接返回 `σ=C·e_mech, C_ep=C`。一致切线 `C_ep = dσ/de_mech` 由收敛系统的线性化 `J[dσ;dΔγ] = −∂R/∂e·de` 精确给出（算法一致性切线，FD 可精确复现）。单轴特例解析可验：加载 `σ11 = σ_y+Hκ`、卸载斜率 E、永久应变累计。

测试（`test/unit_czm_j2.jl`，5 testset）：① 弹性步不变（σ=Ce、C_ep=C）；② 单轴加载屈服后停留屈服面（f=0 精确）+ 卸载永久应变 + 再加载沿弹性斜率至屈服面折返；③ KKT/耗散（σ:Δεᵖ ≥ 0、f≤0、Δγ≥0）；④ 一致切线 FD（塑性流动区扰动 e_mech 三分量，C_ep vs 中心差分 rtol 1e-6）；⑤ 剪切+双轴组合步应力点回面。

- [ ] Step 1: 写失败测试（上述 5 testset， foil 参数取归一化前的物理值直接构造 C=平面应力矩阵）
- [ ] Step 2: 运行确认 `UndefVarError: PlasticState`
- [ ] Step 3: 实现 `src/CzmPlasticity.jl`（完整代码：struct + 返回映射 + 线性化切线 + commit/clone；`CurrentCollector` 四字段 `E_foil=70e9/nu_foil=0.33/sigma_y=0.0/H=0.0` 为通用默认，Jellyroll.jl 显式设 PCC 70/0.33/60e6/0、NCC 110/0.34/200e6/0；NormaliseParam 四行 `/σ_czm`）
- [ ] Step 4: 测试通过；`-e 'include("src/JuBat.jl")'` 干净加载
- [ ] Step 5: 提交 `feat(czm): 平面应力一致 J2 返回映射与物理箔参数（Batch 3 本构层）`

## Task 2: 单元/装配/求解器接入与状态持久化（TDD）

`gl_element_residual_tangent` 增加关键字 `plastic=nothing`（`(σ_y, H, eps_p_vec::Vector{NTuple{3}}, κ_vec::Vector)`，每 GP 一项）：PCC/NCC 且 plastic≠nothing 时，`e_mech = E_vec − ε₀v − eps_p[gp]`，`(σ, C_ep, …) = return_mapping(...)`，σ 替代弹性 S、C_ep 替代 D_mat 进材料切线与 K_G 的 Ŝ；可选 `commit_to` 写回试算状态。`assemble_bulk_residual_tangent`：`plasticity=true` 要求 `geo_nl && mech_state isa Matrix{PlasticState}`（否则 error），PCC/NCC 行取 `foil_params_of`，把每单元 4-GP 状态切片传入；`assemble_coupled_system`/`solve_czm_step`/basic/newton 增 `plasticity::Bool=false, mech_state=nothing` 透传（与 geo_nl 同模式，含线搜索）。`CzmLayout` 增 `plastic_states::Union{Nothing,Matrix{PlasticState}}=nothing`；`update_czm_damage!`：`case.opt.czm_j2_plasticity` 时断言 σ_y>0（指明材料层）、惰性初始化状态矩阵、求解后以 `result.displacement` 调 `commit_plastic!`。

测试（`test/test_czm_j2_integration.jl`，4 testset）：① 塑性关 + geo 关 → 与 Batch 1 逐位（回归锚）；② geo 开 + 塑性关 → 与 Batch 2 C1 参照一致；③ 全开 + 强本征应变（压缩 ε₀≈−3e-3，超 Cu 箔屈服应变 200/110000≈1.8e-3）：PCC/NCC 单元 κ>0、其余层 κ=0、求解收敛；④ 跨两次求解状态累计（同载荷重解 κ 不变——max-history 型不可逆性由提交时序保证）+ 耗散非负抽查。

- [ ] Step 1: 写失败测试（复用 `build_c1_fixture` 模式，nθ=8）
- [ ] Step 2: 确认失败（plasticity 关键字不存在）
- [ ] Step 3: 按"架构"段实现（czm.jl 塑性分支 ~40 行、CzmSolve 透传、CouplingState ~25 行）
- [ ] Step 4: 定向测试全绿 + `test_czm_geo_c1.jl`/`test_czm_geometric_stiffness.jl`/`test_czm_mech_core.jl` 回归绿
- [ ] Step 5: 提交 `feat(czm): J2 塑性接入 GL 单元与求解链路（mech_state 消费、收敛提交、CzmLayout 持久化）`

## Task 3: 三道门禁与批次记录

- [ ] 全套 28/28；探针快照逐位；testexample 冻结指标 + PNG SHA 一致
- [ ] `Simplify/baseline.md` 追加行；progress 追加 Batch 3 小节（D-B3-0/1/2 执行情况、参数核实引用、C2 证据）；index 更新
- [ ] 提交 `docs(baseline): Batch 3 门禁记录（C2 达成）`

---

## 自评审记录

1. **spec 覆盖**：§3.3 全条目（平面应力一致返回映射/C_ep/参数默认/状态/耗散断言/回滚→D-B3-2 原子性等价）→ Task 1/2；§4.1 Czm.jl PlasticState、CouplingState 持久化 ✓；§4.2 PlasticState 定义（eps_p NTuple{3}, κ）✓、MechHistory 仍延后（u_committed 已由 czm_layout.u_prev 承担，delta_core Batch 5——偏差登记第 5 项）；§7 Batch 3 测试五项 → Task 1 四项 + Task 2 集成；C2/L1 判据 → Task 2 testset ③。
2. **决策登记**：D-B3-0（箔本构切换，文献 §10.7 依据）、D-B3-1（塑性需 geo_nl）、D-B3-2（收敛提交原子性）、 MechHistory 继续延后（第 5 项偏差）。
3. **参数核实**（联网，findings 记录）：Al 15 μm 硬态 UTS 150–290 MPa、σ_y≈0.9UTS（软态 40–150），60 MPa 为 Shi 2026 保守起始 ✓；ED Cu σ_y 108–441 MPa，200 MPa 居中 ✓；薄硬态箔硬化近耗尽 → H=0 合理，扫描 0–2 GPa。
4. **类型一致性**：`PlasticState(eps_p::NTuple{3,Float64}, kappa::Float64)`；`mech_state::Matrix{PlasticState}`；返回映射签名在测试/单元/装配三处一致。
5. **无占位符**：返回映射的 4×4 Newton 实现要点（残差、雅可比、线性化切线公式）在 Task 1 Step 3 给出完整规格；测试断言值由解析特例（单轴）给出。
