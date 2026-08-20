# 堆芯塌陷力学建模（工况 C）设计规格

- **日期**：2026-08-20
- **状态**：待用户评审（评审通过前不进入代码实现）
- **任务记录**：`docs/planning-with-files/堆芯塌陷力学建模/`
- **理论基线**：`Theory/` v3.0（各向异性表述已按 Batch 0' 修订：逐层各向同性 + 叠层等效）
- **前置任务**：`理论与模型差异评审`（差距识别）、`界面术语统一`（计数契约）、Batch 0'（Theory/03 修订，已完成并推送）

---

## 1. 背景与目标

Theory v3.0 的工况 C（堆芯塌陷）要求"显式逐层非线性 bulk + 几何刚度 + 箔塑性 + CZM + 多圈状态 + 闭合反馈"的最小能力包；当前代码是"逐层各向同性小应变线弹性 bulk + 箔–涂层 CZM + 顺序耦合代理反馈"。本设计在不破坏任何现有数值行为的前提下，分批补齐该非线性核心，最终通过 C1/C2 全量与 C4-lite 验收，获得堆芯塌陷模拟能力（完整 C3/C4 宣称等 Batch 8 接触摩擦完成后解锁）。

**非目标（当前批次明确不做）**：

- 不做单层本征正交各向异性/横观各向异性本构（已废除，见 Theory/03 注 2.1b）。
- 不做 SP–涂层单边接触与 Coulomb 摩擦（Batch 8 后置；此前 SP–电极界面按完美粘结假设）。
- 不做 5×5 monolithic Jacobian（仅在 operator-splitting 验证不足时再立项）。
- 不改变默认（所有新开关关闭时）的任何数值结果、公开 API 与结果键。

## 2. 已确认决策（用户拍板，不可在实现中变更）

| # | 决策 |
|---|---|
| D1 | 各向异性只存在于叠层结构层面：每个材料各向同性（`E_i, ν_i`），切向并联/法向串联产生等效各向异性；Theory 已同步修订 |
| D2 | SP–电极界面先按完美接触（完美粘结）：层内 SP 面共享节点，跨匝 `phi_pairs` 绑定；接触/摩擦后置 Batch 8 |
| D3 | `czm_enabled` 为唯一主开关；新能力全部为默认关的子选项；`czm_enabled=true` 单独使用与现状行为完全一致 |
| D4 | 求解先 operator-splitting（Theory 07 §6.6.3 策略 B + 同时间步外层固定点）；‖ΔD‖>0.05 且外层失效时才考虑 monolithic |
| D5 | 新反馈（连续 `R_contact(D)/A_eff(D)`、界面热阻恢复）全部 opt-in，默认路径保留现有代理 |
| D6 | 双文档：本规格在 `docs/superpowers/specs/`，任务三件套在 `docs/planning-with-files/堆芯塌陷力学建模/` |
| D7 | 每批遵守 AGENTS 9.6/9.7 基线门禁（Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no`；按记录精度判等） |

## 3. 物理与数学模型

### 3.1 本构（逐层各向同性 + 叠层等效）

- 工况 C 每层用各向同性平面应力刚度 `C_i = E_i/(1-ν_i²)·[1 ν 0; ν 1 0; 0 0 (1-ν)/2]`（Theory/03 式 2.36），现有 `moduli_of` 查表路径不变。
- 工况 R 的叠层等效矩阵（式 2.37–2.38）仅为理论参考；本设计所有新路径均为显式分层（工况 C 路径），不实现等效矩阵凝聚。
- 参数完整性：`E_i/ν_i` 缺失或非法（`E≤0`、`ν∉(-1,1)`）→ `error` 并指明材料层；不借用、不置零。

### 3.2 几何非线性（Batch 2，C1）

- von Kármán 运动学（Theory/02 式 1.30a）：`ε_ss = ∂u_s/∂s + u_n/ρ_s + ½(∂u_n/∂s)²`；剪应变含二次项；曲率 `κ_ss = -∂²u_n/∂s²` 进入弯曲能。
- 总 Lagrangian 实现：Q4 上 Green-Lagrange 应变 = 线性部分 + 二次膜项；内力 `f_int = ∫Bᵀσ`，切线 = 材料切线 + 几何刚度 `K_G(σ)`（当前应力状态的初应力矩阵，含曲率修正项）。
- `czm_geo_nonlinear=false` 时严格走现有线性 B 矩阵代码路径（不是同一公式取零极限），保证基线零漂移。
- 不叠加 CLT `D` 项（Theory/03 Step 4）。

### 3.3 PCC/NCC J2 塑性（Batch 3，C2/L1）

- 仅 `i∈{PCC,NCC}`；屈服 `f = √(3/2)‖s‖ - (σ_y0 + Hκ)`；平面应力 `σ_zz=0` 径向返回映射 + 一致切线 `C_ep`（标准 plane-stress J2 算法）。
- 参数：`param.PCC.sigma_y`（Al 默认 60 MPa）、`param.NCC.sigma_y`（Cu 默认 200 MPa）、`param.PCC.H`/`param.NCC.H`（默认 0，理想塑性；>0 时各向同性硬化）。`ChooseCell` 按 AGENTS 9.4 模式：工况 C 塑性开关开而缺 `sigma_y` → `@warn` + 求解入口 `assert` 拦截。
- 高斯点状态：`eps_p`（3 分量）、`κ`；trial 状态非收敛不提交（复用 `clone_damage_states` 快照/回滚模式）。
- 塑性耗散 `σ:ε̇_p ≥ 0` 为测试断言。

### 3.4 Φ 跨匝完美粘结（Batch 4'，小批次）

- `czm_phi_bond=true` 时在机械子网格构建阶段合并 `phi_pairs` 中的重合节点对（外匝节点并入内匝节点索引，连接表重写）：零穿透、零分离、零滑移。
- 默认关时 `phi_pairs` 行为不变（仅拓扑记录）。这是 `phi_pairs` 的第一个求解器消费者。
- 声明义务：完美粘结抑制 SP 滑移先皱机理（Theory 04:174、07:209）；结果输出必须携带 `collapse_approx = "phi_perfect_bond"` 标记，直至 Batch 8。

### 3.5 Δ_core 与塌陷状态（Batch 5，C4-lite）

**Δ_core 离散定义（本规格冻结）**：

1. 取机械子网格**最内条螺旋边界**（第 1 条层界螺旋）的全部节点在**当前（变形后）构型**的坐标 `(x_i, y_i)`；
2. 最小二乘拟合圆（圆心 `(x_c, y_c)`、半径 `r_fit`）；
3. `Δ_core = max_i |r_i - r_fit| / r_fit`（归一化），同时输出绝对幅值 `w_core = max_i |r_i - r_fit|`（m）。

跨圈演化输出 `Delta_core [-]`、`core wrinkle amplitude [m]`；多圈单调恶化（允许物理合理的非单调）为 C4-lite 判据。

**多圈状态持久化**（挂在 Case 层，仿 `czm_layout.u_prev`）：

- `u`/当前几何（更新参考构型，**禁止每圈重置完美圆**）；
- 逐高斯点 `eps_p/κ`；cohesive `D`（现有 `damage_states`）；
- 每半循环结束时提交；失败步全部回滚。

**路径跟踪**：把 `solve_czm_arc_length_step`（Crisfield 圆柱弧长）从"cohesive 载荷参数"扩展为"全机械残差"（含 K_G、塑性、Φ 绑定贡献）；负切线/分岔附近步长缩减（下限 `step/128` 后报错终止，不伪造收敛）。

### 3.6 损伤–电–热反馈（Batch 6，opt-in）

- `czm_continuous_feedback=true`：分支电压/过电位残差加入 `R_contact(D)/A_eff(D)`（替换仅失效置零与阈值面积权重两条代理路径；代理路径保持默认可用）；总电流守恒约束不变。
- 界面热阻：按 `cohesive_to_thermal` N:1 物理归约恢复 `k_n(D)` 装配（当前被注释禁用的代码段重写启用），并加 `Q_contact = j²R_contact/A_eff`；D=0 时必须与现模型完全一致。
- 同时间步外层固定点：`I/T/u/D` 联合迭代，收敛判据为残差与状态增量双阈值；不收敛 → `error`，不提交部分状态。

## 4. 架构与组件

### 4.1 文件级变更表

| 文件 | 批次 | 变更 |
|---|---|---|
| `src/Option.jl` | 1 | 新增默认关子选项（§5 表）；不改现有字段语义 |
| `src/Czm.jl` | 1-3 | 新增 `assemble_bulk_residual_tangent`（bulk 残差/切线，三槽位之一）；`PlasticState` 类型；`assemble_coupled_system_full` 改为调用新入口（开关全关时逐句等价） |
| `src/CzmSolve.jl` | 2,5 | Newton/弧长子步内每次迭代重组切线（不再用常量 `K_bulk` 缓存；弹性+无塑性+无 K_G 时仍走缓存快路径）；弧长扩展为全机械残差 |
| `src/CouplingState.jl` | 1,3,5,6 | `CZMAssemblyCache` 增加参考构型与机械状态字段；塑性状态随 Case 持久化；外层固定点迭代器 |
| `src/Jellyrollmodel.jl` | 4' | `czm_phi_bond=true` 时 phi_pairs 节点合并 |
| `src/CycleSolver.jl` | 5 | phase 交接扩展：传递机械状态与当前几何，不重置完美圆 |
| `src/parameters/Jellyroll.jl` | 3 | `PCC/NCC.sigma_y`、`PCC/NCC.H`（默认 0/未设） |
| `src/SPMe.jl`、`src/ThermalDistributed.jl` | 6 | 连续反馈入口（opt-in）；热阻装配恢复 |
| `test/` | 各批 | 新测试文件见 §7 |

### 4.2 关键新接口（内部，不进公共文档）

```julia
# 残差/切线分列入口（contact 槽位 Batch 8 前为空实现）
assemble_bulk_residual_tangent(czm_mesh, u, param_cache, mech_state; geo_nl, plasticity)
    -> (f_int_bulk::Vector{Float64}, K_tangent::SparseMatrixCSC)

mutable struct PlasticState
    eps_p::NTuple{3,Float64}   # (ss, nn, sn) 塑性应变
    kappa::Float64             # 等效塑性应变
end

struct MechHistory        # Case 层持久化
    plastic::Matrix{PlasticState}     # [elem, gp] 仅 PCC/NCC 有效
    u_committed::Vector{Float64}      # 当前几何（参考构型更新源）
    delta_core::Float64
end
```

### 4.3 数据流（工况 C 开启时）

```
每 czm_update_interval 时间步:
  CallModel (电-热) → T, soc 场
      ↓ thermal_to_czm 插值
  CZM 力学步 (Newton/弧长):
      assemble_bulk_residual_tangent   ← u, PlasticState, 参考构型
      + assemble_czm_system            ← DamageState
      [+ phi_bond 绑定贡献]
      → 收敛: 提交 u, eps_p/κ, D, Δ_core;失败: 回滚全部 trial 状态
      ↓ (Batch 6 起, opt-in) R_contact(D)/A_eff(D), k_n(D)
  外层固定点: 重复 CallModel↔CZM 直至 ‖ΔD‖,‖Δu‖,‖ΔT‖ < tol
```

## 5. 配置接口（`src/Option.jl` 新增字段）

| 字段 | 类型/默认 | 批次 | 语义 |
|---|---|---|---|
| `czm_geo_nonlinear` | `Bool=false` | 2 | von Kármán 残差 + K_G |
| `czm_j2_plasticity` | `Bool=false` | 3 | PCC/NCC J2（要求 `sigma_y` 参数存在） |
| `czm_phi_bond` | `Bool=false` | 4' | Φ 跨匝完美粘结（节点合并） |
| `czm_continuous_feedback` | `Bool=false` | 6 | 连续损伤–电–热反馈 + 界面热阻 |
| `czm_friction_mu` | `Float64=0.10` | 8（预留） | SP Coulomb 摩擦系数（当前无消费者） |

兼容性契约：以上全 false 时，行为、结果键、`tools/czm_baseline_probe.jl` 三方法快照与现状逐指标一致。

## 6. 错误处理约定（继承 AGENTS 9.7）

- 缺参（`E/ν/sigma_y`）→ `error` 指明材料层；不默认、不置零、不截断。
- 非收敛：子步回滚并减步长；至下限仍失败 → `error` 终止并输出诊断（残差史、活跃集、步长），不提交部分状态。
- 状态提交原子性：`u/εp/κ/D/Δ_core` 同步提交或同步回滚。
- 不引入任何新的静默回退分支。

## 7. 测试与验收门

| 批次 | 新测试 | 验收门 |
|---|---|---|
| 1 | `test_czm_mech_core.jl`：开关全关时新入口 ≡ `K_bulk*u`（按记录精度）；缓存不变量 | `testexample.jl` 基线一致；`czm_baseline_probe.jl` 快照不变 |
| 2 | `test_czm_geometric_stiffness.jl`：切线有限差分、刚体运动零应变、单元 patch、自由膨胀零应力、K_G→0 退化回 Batch 1 冻结解 | **C1** |
| 3 | `unit_czm_j2.jl`：单轴屈服/卸载、KKT、耗散非负、一致切线有限差分、trial 回滚 | **C2/L1**：强约束+完美圆仅预屈曲场；塑性关退回 C1 |
| 4' | `test_czm_phi_bond.jl`：绑定对位移一致、默认关时网格逐位不变 | 基线不变 |
| 5 | `test_czm_multicycle_state.jl`：跨圈状态持久化、失败回滚、Δ_core 计算（构造解析圆+扰动算例） | **C4-lite**：多圈 Δ_core 与 D 联合增长（带 phi_perfect_bond 声明） |
| 6 | `test_czm_feedback.jl`：D=0 退化一致、D↑ 方向正确、外层不收敛报错 | 反馈方向正确；总电流/热功率守恒 |
| 7 | 汇总：耗散非负、守恒、网格/步长敏感性、Φ 完整性、运行时间冒烟 | 全部门禁通过后才允许 C4-lite 宣称 |
| 8 | （后置）接触测试组 | **C3+C5**；解锁完整宣称 |

每批均跑 `example/testexample.jl` 冻结门禁；`unit_czm_eigenstrain.jl` 既有失败（58/60）在 Batch 1 开工时登记为前置状态，不作为新回归。

## 8. 风险与停止条件

| 风险 | 缓解/停止条件 |
|---|---|
| 基线漂移 | 任一强制指标不一致 → 停止回退当批 |
| 切线实现错误 | 有限差分校验失败 → 不进下一批 |
| 弧长在分岔点反复失败 | 步长减至下限 → 报错终止（可诊断），不伪造收敛 |
| 性能（每迭代重组切线） | 弹性+无塑性+无 K_G 走缓存快路径；Batch 7 加运行时间冒烟预算；超预算先优化再继续 |
| 完美粘结误导解读 | 结果键携带 `collapse_approx` 声明；宣称门禁在 Batch 8 前禁止"SP 先皱/完整塌陷"表述 |

## 9. 批次依赖与执行顺序

```
0'(done) → 1 → 2(C1) → 3(C2) → 4' → 5(C4-lite) → 6 → 7 → [8(C3/C5, 后置另立]
```

每批一次提交（现有提交风格），同步更新 planning 三件套与总索引；规格评审通过是 Batch 1 代码工作的前置条件。

## 10. 自评审记录（2026-08-20）

- 无 TBD/占位符；Δ_core 定义已冻结（§3.5）。
- 与 Theory 修订版（逐层各向同性）完全一致；无本征各向异性残留要求。
- 与用户决策 D1–D7 逐条对应；所有新能力 opt-in，默认零漂移。
- 范围单一（力学补充建模），可由本规格 + 既有代码直接派生实现计划。
