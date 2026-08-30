# 力学结构体根本性重构设计（Mechanics Struct Refactor）

- **日期**: 2026-08-30
- **状态**: 已与用户逐节确认（brainstorming 流程），待 spec 审阅
- **性质**: 纯重构，物理零变更；验收标准 = 与基线 `testexample-20260830T005629+0800`（v7）逐位一致
- **迁移策略**: 大爆炸单批落地（用户选定），批内按层推进、单次验证门

---

## 1. 背景与动机

当前力学侧存在四层结构体，经 2026-08-29 会话逐字段审查确认以下事实：

1. `Cohesive`（SetParams.jl，28 字段平铺后缀）→ `CzmInterfaceParams`（23 字段）存在 17 字段改名复制，构造体 ~50 行纯字段搬运；`E_eff/ν/E_star/L_ch` 零消费（分层 `moduli_of` 已旁路）；
2. `Λ = scale.L/δ_czm` 全局量重复存储于两个界面实例；
3. `czm_model` 存在 Option↔参数双写（`CouplingState.jl` 反向突变）；`tau_visc` 存在 Option→参数拷贝（`SetCase.jl`）；
4. 快照失效靠内容哈希补丁（Task 4.4 修复史），而装配缓存 `CZMAssemblyCache` 的失效键本来就绑定 `objectid(czm_mesh)`；
5. 演化状态散布于 `CzmLayout`（u_prev/node_ref/prestress/plastic_states）与 `czm_mesh.damage_states`，损伤更新走 `clone_czm_mesh_with_damage` 克隆链。

用户否决了"Cohesive 嵌套化"过渡方案，指定按**电化学侧结构体惯用法**（属性挂物理部件、原地归一化、派生量挂原处如 `M_d`、求解器直读字段）根本性重构全部四层。

## 2. 目标与非目标

**目标**（用户确认，按优先序）：
1. 单一事实来源：每物理量只存一处；
2. 电化学范式一致性：`CurrentCollector` 承载界面属性（与 `PE.E_coat`、`PCC.sigma_y` 同构）；
3. 为 SP–涂层接触（AGENTS §9.8）预留结构扩展位；
4. API 收敛：Case 挂载 4→2、求解签名 17 参数→6、Option 收敛为嵌套子结构。

**非目标**：
- 不改任何物理模型与数值结果（验收即逐位一致）；
- 不做输入/视图类型级分离——采用**单类型原地归一化**（用户明确选定，接受输入身份被归一化覆盖，与电化学一致）；
- 电化学侧结构体不动（未来若统一另立批次）；
- 不实现接触本身（仅预留位）。

## 3. 决策记录

| 决策点 | 结论 |
|---|---|
| 范围 | 四层全含：参数、缓存/状态、网格/拓扑、Option/Case |
| 分离语义 | 单类型原地归一化（电化学式），放弃类型级分离 |
| 参数组织 | 界面参数挂 `CurrentCollector`（用户提案）；`:PE_PCC`↔`PCC`、`:NE_NCC`↔`NCC` 1:1（层序 PCC 两面皆 PE、NCC 两面皆 NE，无歧义） |
| 状态/缓存归属 | 方案 3：装配缓存挂 `CohesiveMesh`（gs 先例、对象身份即失效判据）；演化状态聚合 `MechState` |
| 迁移 | 大爆炸单批 |

**损伤-刚度解耦事实**（支撑缓存设计）：损伤只改变 K_coh 切线——每次 Newton 迭代现算、从不缓存；被缓存的 K_bulk 只依赖（冻结的）参数与网格，弹性路径恒定、geo_nl/塑性路径禁用之。

## 4. 设计

### 4.1 参数层：CurrentCollector 承载界面，Cohesive 消亡

`CurrentCollector` 新增 16 字段（量纲输入）：

```julia
# 界面本构（该集流体与其相邻涂层之间的 COH 界面；法向无后缀、切向带 _t）
σ_max, K_n, δ_0, G_c, δ_c            # Mode I
τ_max, K_t, δ_0_t, G_c_t, δ_c_t      # Mode II
eta                                    # BK 混合模式指数
# 界面热阻（md/07）
h_c0, k_air, lambda_m, beta, threshold
```

- `Cohesive` struct、`Params.cohesive` 字段、`CzmInterfaceParams`、`CzmParamCache`、`compute_czm_params_per_interface` 全部删除；`E_eff/ν/E_star/L_ch` 死字段随之消亡；
- `czm_model`、`tau_visc` 删除（Option 已有等价字段，见 4.4）；`update_czm_damage!` 的反向突变与 `SetCase` 的 `tau_visc` 拷贝链消失；
- Λ 不再存储：使用点内联 `param.scale.L / param.scale.δ_czm`（czm.jl 装配、compute_gap_conductance 两处）；
- `NormaliseParam` 归一化块（现 SetParams.jl:503+ 的类比）：`σ_max/τ_max ÷ σ_czm`、`K_n/K_t ÷ K_czm`、`δ_0/δ_c/δ_0_t/δ_c_t ÷ δ_czm`、`G_c/G_c_t ÷ G_czm`；`eta/h_c0/k_air/lambda_m/beta/threshold` 沿用现状归一化（无因次或 L 空间，与今天一致）；
- 锚点：`scale.σ_czm = PCC.σ_max`（量纲阶段）、`scale.δ_czm = 2·PCC.G_c/PCC.σ_max`；缺参回退（`scale.δ_czm = scale.L`）与 `ChooseCell` @warn 改查 `PCC.σ_max/PCC.G_c`；
- 热路径分派：`ip = iface === :PE_PCC ? param.PCC : param.NCC`；`bilinear_traction_state/bilinear_tangent` 改收 collector 实例 + `czm_model` 作为显式参数（本构分支不再从参数结构体读模型选择）；`compute_gap_conductance` 改收 `(D, δ_n, ip, param)`（Λ 由 param.scale 内联）；
- 共用值（eta、五个热阻参数）在 PCC/NCC 各写一份，由参数集显式赋相同值——电化学惯用法接受实例间重复（`PE.rho/NE.rho` 先例），比隐藏的结构体复制诚实。

### 4.2 网格/缓存层：CohesiveMesh 承载装配缓存

`CohesiveMesh` 新增字段：

```julia
K_bulk::Union{Nothing, SparseMatrixCSC}   # 惰性：首次线性装配时构建，此后跨步只读
cohesive_geom::Union{Nothing, Vector{CohesiveElementGeom}}  # 纯几何派生标架（gs 同款）
ws::Union{Nothing, CZMAssemblyWorkspace}  # 预分配工作区
```

- **BC 不缓存**：`identify_bc_nodes_czm` 每次 solve 入口现算（O(nnode)，成本可忽略），消灭 `fix_inner` 作为缓存键的问题；
- `CZMAssemblyCache`、`ensure_czm_cache`、`param_cache.id` 内容哈希全部删除；
- 失效语义唯一化：**对象身份**——重建 mesh 即丢弃缓存；复用 mesh 即缓存有效。K_bulk 新鲜度由 §5 参数冻结契约保证；
- `damage_states` 迁出网格（→ MechState），`clone_czm_mesh_with_damage` 克隆链整体删除；
- `CohesiveElement/CzmSubmesh/phi_pairs` 拓扑类型保留（本会话审查无结构性问题）。

### 4.3 状态层：MechState 聚合全部演化状态

```julia
mutable struct MechState
    u_prev::Vector{Float64}               # 上步收敛位移
    node_ref::Union{Nothing, Matrix{Float64}}  # 初始螺旋快照（Δ_core 冻结基准）
    damage_states::Vector{DamageState}    # ← czm_mesh 迁入
    plastic_states::Union{Nothing, Matrix{PlasticState}}  # ← CzmLayout 迁入
    winding_prestress::Union{Nothing, Vector{Tuple{Float64,Float64,Float64}}}  # ← CzmLayout 迁入（winding_prestress_field 返回型）
    contact::Nothing                      # AGENTS §9.8 SP–涂层接触预留位
end
```

- `CzmLayout` 消亡；`case.mech` 随 SetCase（czm_enabled）创建；
- **提交语义**：求解器在收敛后**原位**提交损伤/塑性/位移到 ms（与塑性 D-B3-2"收敛提交"同时机），失败回滚 = 不触碰 ms（试探态为局部变量）；求解器返回值从 `(result, updated_czm_mesh)` 变为 `result::CZMResult`。

### 4.4 Option/Case 层：收敛挂载

20 个 `czm_*` 平铺字段收敛为 `opt.czm::CzmOptions` 嵌套子结构（默认值沿用现状）：

```
enabled, model, update_interval, soh_threshold, inner_exit_only, fix_inner,
iter_method, max_iter, tol, load_steps, arc_length_alpha,
viscous_enabled, viscous_tau,
area_loss_enabled, area_loss_threshold,
geo_nonlinear, winding_prestress, j2_plasticity,
continuous_feedback, friction_mu（Batch 8 预留）
```

- `opt.mechanicalmodel` 保留顶层（模型选择器家族：none/full）；
- Case 挂载 4→2：`case.czm_mesh`（拓扑+缓存）+ `case.mech`；`czm_param_cache/czm_cache/czm_layout` 消亡；
- `visc_beta` 直接从 `opt.czm.viscous_tau` 计算（原 `param.cohesive.tau_visc` 链删除）。

### 4.5 调用链（改后全貌）

```julia
# Solve.jl CZM 更新块（现 ~40 行准备代码 → ）：
res = solve_czm_step(case.czm_mesh, case.mech, param, F_ext, opt.czm;
                     dT=dT_elem, Δsn=Δsoc_n, Δsp=Δsoc_p)
# 求解器内部：收敛即提交 ms；返回 CZMResult
```

## 5. 契约

1. **参数冻结契约**（写入 AGENTS）：`SetCase` 归一化后不得修改 `param` 任何字段。K_bulk/界面缓存/标架的新鲜度由此保证（电化学 `M_d` 同款敞口，从内容哈希补丁改为声明式契约）。
2. **失效语义**：缓存随宿主对象生灭（换 mesh = 缓存作废）；不存在任何手动失效调用。
3. **提交语义**：演化状态只在收敛后写入 ms；任何试探/回滚不触碰 ms。

## 6. 字段迁移映射表

| 旧（`cohesive.X`） | 新位置 |
|---|---|
| `σ_max_pe_pcc / σ_max_ne_ncc` | `PCC.σ_max / NCC.σ_max` |
| `K_n_pe_pcc / K_n_ne_ncc` | `PCC.K_n / NCC.K_n` |
| `δ_0_pe_pcc / δ_0_ne_ncc` | `PCC.δ_0 / NCC.δ_0` |
| `G_c_pe_pcc / G_c_ne_ncc` | `PCC.G_c / NCC.G_c` |
| `δ_c_pe_pcc / δ_c_ne_ncc` | `PCC.δ_c / NCC.δ_c` |
| `τ_max/K_t/δ_0_t/G_c_t/δ_c_t` 同构 ×2 | `PCC.* / NCC.*`（切向后缀 `_t` 保留） |
| `eta` | `PCC.eta` + `NCC.eta`（各一份） |
| `h_c0/k_air/lambda_m/beta/threshold` | `PCC.* + NCC.*`（各一份） |
| `czm_model` | 删除 → `opt.czm.model`（读点：`bilinear_*` 显式参数） |
| `tau_visc` | 删除 → `opt.czm.viscous_tau` |
| `CzmInterfaceParams.Λ` | 删除 → 使用点内联 `scale.L/δ_czm` |
| `CzmInterfaceParams.E_eff/ν/E_star/L_ch` | 删除（零消费） |
| `CzmInterfaceParams.α` | 已于批次①删除（α/β 分层化，基线 v7） |
| `Option.czm_*`（20 字段） | `opt.czm.*`（同名去前缀，`czm_visc_tau→viscous_tau` 等） |
| `case.czm_param_cache/czm_cache/czm_layout` | 删除 → `case.czm_mesh`（缓存内嵌）+ `case.mech` |
| `czm_mesh.damage_states` | `case.mech.damage_states` |

## 7. 验证计划（大爆炸单批，单次验证门）

1. **逐位一致门**：`testexample.jl` 快门与基线 v7 全指标逐位一致（电压/温度/应力范围/分离 9.4407e-13 m/19 步/零损伤）；
2. **全套测试**：`runtests.jl` 34/34（8 个 CZM 测试文件适配新 API；`verify_czm_per_interface` 的 `retune` 辅助改为调谐 PCC/NCC 实例）；
3. **绘图门**：`couple_example.jl` 适配后运行，四张 PNG SHA-256 与 v7 记录一致（`540fe42f/ecdc9f58/2c29e35b/5124dec3`）；
4. **契约冒烟**：geo_nl/plasticity/prestress 路径各跑既有专项测试确认 K_bulk 未被错误消费；
5. 任一数值漂移 = 定位或回滚，不得带差异进入下一批。

## 8. 波及面与实施顺序

代码：`SetParams.jl`、`parameters/Jellyroll.jl`（及其他参数集）、`czm.jl`、`CzmSolve.jl`、`CouplingState.jl`、`CzmMesh.jl`、`Mechanical.jl`、`Materialmatrix.jl`、`Option.jl`、`SetCase.jl`、`Solve.jl`、`CycleSolver.jl`、`Variables.jl`、`JuBat.jl`（导出表）；测试 8 文件；示例 3+ 脚本；文档 AGENTS §9.4、md/01/06/15/07、对照/06、函数索引 6 篇。粗估 90–120 处调用点。

批内推进顺序（单一验证门在最后）：① Option/参数层（字段迁移+归一化+锚点）→ ② 删除参数缓存与 Λ 内联 → ③ 网格缓存挂载+BC 现算 → ④ MechState+提交语义+克隆链删除 → ⑤ 调用点适配（src/测试/示例）→ ⑥ 文档 → ⑦ 验证门。

## 9. 风险与开放问题

1. **与堆芯塌陷 Batch 1+ 并发协调**：同一批文件（CouplingState/czm/CzmSolve）正在按 2026-08-20 spec 演进；本重构落地前需确认 Batch 1 规格评审状态，避免两条线冲突（建议：本重构先落，或冻结其一）；
2. **参数冻结契约的执行**：语言层面不可强制（原地字段写无法拦截），靠文档+code review；如未来需要硬保证，再评估 debug 模式采样断言（本期不做）；
3. `bilinear_*` 增加 `czm_model` 参数后 Materialmatrix.jl 的调用面变化需同步检查间隙导热/损伤批更新路径。
