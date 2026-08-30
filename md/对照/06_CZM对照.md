# 卷 06 — CZM 理论与当前实现对照

> 更新日期：2026-08-30
> 理论来源：`md/06_内聚力模型_CZM.md`、`md/07_界面热阻模型.md`、`md/14_粘性正则化.md`、`md/15_颗粒与极片模量区分.md`
> 实现来源：`src/Materialmatrix.jl`、`src/czm.jl`、`src/CzmSolve.jl`、`src/CouplingState.jl`、`src/CzmMesh.jl`、`src/SetMesh.jl`、`src/SetParams.jl`、`src/Option.jl`

本文以 2026-08-30 力学四层结构重构后的接口为准。旧的 `Cohesive`、
`CzmInterfaceParams`、`CzmParamCache`、`CZMAssemblyCache`、`CzmLayout` 和
`case.czm_param_cache/czm_cache/czm_layout` 已删除，不再作为现行实现依据。

---

## 1. 结构与归属

| 层 | 当前宿主 | 内容 | 契约 |
|---|---|---|---|
| 参数 | `Params` 中的 `PCC/NCC::CurrentCollector` | 连续层、塑性、CZM 界面本构和热阻字段 | `SetCase` 归一化后冻结 |
| 选项 | `Option.czm::CzmOptions` | 20 个模型/求解/耦合选项 | 模型与粘性只在此处单一持有 |
| 拓扑/缓存 | `CohesiveMesh` | bulk/cohesive 拓扑、映射、`K_bulk/cohesive_geom/ws` | 缓存随 mesh 对象生灭；BC 不缓存 |
| 演化状态 | `MechState` | `u_prev/damage/plastic/prestress/node_ref/contact` | 只在收敛后提交，失败不触碰 |
| 单步结果 | `CZMResult` | 位移、损伤、牵引、分离和收敛诊断 | 返回结果，不返回克隆 mesh |

参数分派唯一规则：

```julia
collector_params(param, :PE_PCC) === param.PCC
collector_params(param, :NE_NCC) === param.NCC
```

层序中 PCC 两面均邻接 PE、NCC 两面均邻接 NE，因此四个真实箔–涂层面复用两种
`interface_type`。完整离散计数仍为 `4 * (length(theta) - 1)`；
`CohesiveMesh.n_layers == 2` 是遗留的本构类型计数，不是物理层数。

---

## 2. 字段迁移对照

| 旧位置 | 当前位置 |
|---|---|
| `cohesive.σ_max_pe_pcc / σ_max_ne_ncc` | `PCC.σ_max / NCC.σ_max` |
| `cohesive.K_n_pe_pcc / K_n_ne_ncc` | `PCC.K_n / NCC.K_n` |
| `cohesive.δ_0_pe_pcc / δ_0_ne_ncc` | `PCC.δ_0 / NCC.δ_0` |
| `cohesive.G_c_pe_pcc / G_c_ne_ncc` | `PCC.G_c / NCC.G_c` |
| `cohesive.δ_c_pe_pcc / δ_c_ne_ncc` | `PCC.δ_c / NCC.δ_c` |
| 两套 `τ_max/K_t/δ_0_t/G_c_t/δ_c_t` | `PCC.* / NCC.*`（切向后缀 `_t`） |
| `cohesive.eta` | `PCC.eta` 与 `NCC.eta` 各一份 |
| `cohesive.h_c0/k_air/lambda_m/beta/threshold` | `PCC.*` 与 `NCC.*` 各一份 |
| `cohesive.czm_model` | `opt.czm.model` |
| `cohesive.tau_visc` | `opt.czm.viscous_tau` |
| `CzmInterfaceParams.Λ` | 使用点内联 `scale.L / scale.δ_czm` |
| `czm_mesh.damage_states` | `case.mech.damage_states` |
| `case.czm_param_cache/czm_cache/czm_layout` | `case.czm_mesh` + `case.mech` |

`E_eff/ν/E_star/L_ch` 参数缓存字段已经消亡。bulk 材料模量由 `moduli_of(param, mt)`
逐层读取；全叠合厚度加权只用于参考尺度 `scale.E_coat`，不作为单元本构模量。

---

## 3. 牵引–分离本构

| 理论项 | 当前实现 | 一致性与说明 |
|---|---|---|
| 弹性段 `T = Kδ` | `bilinear_traction_state` | `δ_eff ≤ δ_0_eff` 时静态损伤为 0 |
| 软化段 `T = (1-D_v)Kδ` | `bilinear_traction_state` | 粘性关闭时 `D_v = D`；开启时按 md/14 使用粘性损伤 |
| 完全断裂 `D=1, T=0` | `bilinear_traction_state` | `fractured` 历史不可逆 |
| 一致切线 | `bilinear_tangent` | 显式接收 `ip::CurrentCollector`、`czm_model` 与 `visc_beta` |
| 历史更新 | `update_damage` / `update_damage_per_interface` | 按界面分组并保持单元顺序 |

现行本构接口：

```julia
bilinear_traction_state(δ_n, δ_t, damage_state, ip::CurrentCollector,
                        czm_model::String; visc_beta=1.0)
bilinear_tangent(δ_n, δ_t, damage_state, ip::CurrentCollector,
                 czm_model::String; visc_beta=1.0)
update_damage(damage_states, separations, ip::CurrentCollector,
              czm_model::String; visc_beta=1.0)
```

`czm_model == "model1"` 为 Mode-I 主导分支；`"mix"` 使用法/切向等效分离与
`eta` 插值。模型选择不从材料结构反向读取。

---

## 4. 无量纲化与有限元装配

CZM 锚点来自量纲阶段的 PE-PCC 界面：

```text
scale.σ_czm = PCC.σ_max
scale.δ_czm = 2 * PCC.G_c / PCC.σ_max
scale.G_czm = scale.σ_czm * scale.δ_czm
scale.K_czm = scale.σ_czm / scale.δ_czm
```

`NormaliseParam` 对 PCC/NCC 分别执行：牵引除以 `σ_czm`、刚度除以 `K_czm`、
分离除以 `δ_czm`、断裂能除以 `G_czm`。`eta/beta` 不变；热阻参数转入
`scale.L` 空间。

| 装配内容 | 函数 | 当前数据源 |
|---|---|---|
| 层材料模量 | `moduli_of` | PE/NE 涂层 `E_coat/nu_coat`；SP/PCC/NCC 连续层 `E/nu` |
| 层本征应变 | `eigenstrain_of` | 各层 `alphaT`；PE 接 `Δsp`、NE 接 `Δsn` |
| cohesive 内力/切线 | `assemble_czm_system` | `collector_params(param, interface_type)` |
| 线弹性 bulk 刚度 | `assemble_bulk_stiffness` | `param` 逐层直读 |
| GL/TL bulk 残差/切线 | `assemble_bulk_residual_tangent` | 当前位移、逐层本征应变、可选塑性/预应力 |
| 热化学等效载荷 | `assemble_thermal_chemical_load` | 父热单元温度/SOC 映射 |
| 完整耦合 | `assemble_coupled_system(_full)` | bulk + cohesive + 外载/本征应变载荷 |

分离从长度空间转到 CZM 分离空间时，`assemble_czm_system` 在使用点内联
`Λ = param.scale.L / param.scale.δ_czm`。间隙导热反向使用
`param.scale.δ_czm / param.scale.L`。Λ 不存入任何结构或缓存。

---

## 5. 缓存、BC 与状态提交

| 对象/动作 | 当前实现 |
|---|---|
| 线弹性 bulk 刚度 | `bulk_stiffness` 首次写 `czm_mesh.K_bulk` |
| cohesive 几何 | `cohesive_geometry` 首次写 `czm_mesh.cohesive_geom` |
| 装配工作区 | `assembly_workspace` 首次写 `czm_mesh.ws` |
| 缓存失效 | 重建 mesh；无手动失效和内容哈希 |
| 边界条件 | `extract_bc_dofs` 每次求解入口调用 `identify_bc_nodes_czm` |
| 损伤/位移试探 | 从 `MechState` 深拷贝为局部状态 |
| 收敛提交 | 求解策略原位写回 `ms.damage_states/ms.u_prev` |
| 失败回滚 | 丢弃局部试探态，`ms` 保持不变 |

`K_bulk` 的有效性依赖参数冻结：`SetCase` 归一化后不得修改 `param`。若参数或
拓扑改变，应重建 Case/mesh。

---

## 6. 求解链

统一入口：

```julia
solve_czm_step(czm_mesh, case.mech, param, F_ext, opt.czm;
               dT_elem=..., Δsoc_n_elem=..., Δsoc_p_elem=...) -> CZMResult
```

| `opt.czm.iter_method` | 分派 | 备注 |
|---|---|---|
| `basic` | `solve_czm_basic_step` | 单步 Newton + 回溯 |
| `load_substep` | `newton_raphson_czm` | 自适应载荷子步 |
| `arc_length`，线性几何 | `solve_czm_arc_length_step` | 位移/载荷弧长 |
| `arc_length`，几何非线性 | `solve_czm_arc_geo_step` | 本征应变载荷参数化到 `λ=1` |

`visc_beta` 在统一入口由 `opt.czm.viscous_tau / param.scale.t0` 与当前载荷步尺度计算。
`update_czm_damage!` 对输入/输出有限性和 `result.converged` 进行硬门控；未收敛结果不会
被当作有效工程状态继续传播。

---

## 7. 界面热阻与耦合边界

现行独立接口为：

```julia
compute_gap_conductance(D, δ_n, ip::CurrentCollector, param::Params)
compute_element_gap_conductance(damage_states, elem_idx, ip, param)
compute_all_gap_conductances(damage_states, ip, param)
```

界面热阻字段位于 `PCC/NCC`，但 `ThermalDistributed2D_BC` 中的界面热阻装配块目前
仍被禁用；默认热网格使用合并节点。损伤对热模型的现行影响仅限已有的热源屏蔽/
有效面积路径，不应把独立 `compute_gap_conductance` 接口解释为已启用双向界面热阻耦合。

`CzmSubmesh.phi_pairs` 仅提供匝间节点配对；当前未装配 SP–涂层单边接触、穿透约束或
Coulomb stick-slip。`opt.czm.friction_mu` 仍是预留字段，不能据此宣称已实现摩擦接触。

---

## 8. 输出与后处理

- `CZMResult` 保存归一化位移、损伤、牵引、分离和求解诊断。
- `case.mech.damage_states` 是跨步/跨循环的已提交损伤真源。
- `czm_output_to_variables` 与 CSV/PostProcessing 层负责结果键和单位还原。
- 层分辨宏观应力由 `compute_macro_stress/export_macro_stress` 使用
  `case.mech.u_prev` 在线恢复；CZM-off 工具 `thermal_diffusion_stress_2D` 在
  `mesh_bonded` 域独立求解。

公共输出键、单位与 `CyclingResult.soh` 字段契约不因结构重构改变。

---

## 9. 当前实现自检清单

- [x] 界面参数只在 `PCC/NCC::CurrentCollector`。
- [x] `opt.czm::CzmOptions` 是 CZM 配置唯一宿主。
- [x] Λ 在装配/导热使用点内联。
- [x] 装配缓存随 `CohesiveMesh`，BC 每次入口现算。
- [x] 演化状态聚合在 `MechState`，收敛提交、失败不触碰。
- [x] 求解统一入口返回单个 `CZMResult`。
- [x] 文档明确界面热阻装配与 SP–涂层接触尚未启用。
