# 堆芯塌陷力学建模 Batch 5（多圈状态/Δ_core/弧长，C4-lite）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐 spec §3.5 核心内容——Δ_core 位移基计算与跨圈输出、多圈状态持久化审计、弧长法 geo_nl 解锁（λ 缩放本征应变增量）——达成 **C4-lite**：多圈 Δ_core 与 D 联合增长（携带 SP 粘结/Φ 自由近似声明）；不可达则执行 D10 有界敏感性探针后停止评审。

**用户修订（2026-08-22）**：**Batch 4' 不做**——保持 SP 与活性材料的完美粘结（层内共享节点，现状），Φ 跨匝界面维持自由（节点不合并）；spec §9 的 4' 行留空待后续决策；C4-lite 声明相应为 `collapse_approx = "sp_perfect_bond_phi_free"`（如实反映：SP–电极层内粘结、Φ 缝自由）。

**Architecture:** 三任务。Task 1 Δ_core（§3.5 冻结定义）+ 滤波单测 + 输出键 + 多圈快照；Task 2 弧长 geo 路径（λ-本征应变斜坡，见 D-B5-2）；Task 3 多圈集成测试 + C4-lite 门 + 三道门禁（v2）。

**Spec:** spec v1.4 §3.5（Δ_core 冻结定义/多圈持久化/路径跟踪）、§7 Batch 5 行、§8 D10 行；Theory §6.10。

## Global Constraints

（同前：环境；基线 v2（探针 v2 表、PNG `b31ffb49…`）；默认零漂移（`czm_phi_bond=false` 网格逐位、其余开关默认关）；AGENTS 9.7；每 Task 一提交；§9.9 输出。）

### 设计决策（执行契约）

**D-B5-1（Δ_core 冻结实现）**：`core_ovalization(czm_mesh, u, ref_node) -> (w_core, Δ_core)`——Γ_in,free = 内螺旋（螺旋 1）**第一匝窗口**节点（`1:nθ`，同探针方法）；`u_n,i = u·r̂`（r̂ 取初始螺旋节点向径，**始终相对初始螺旋**）；窗口内离散 Fourier 去除 n=0（均值）与 n=1（刚体平移）后 `w_core = max|ũ_n|`；`Δ_core = w_core/r_ref`，`r_ref = a`（`cell.Rin`，11.111 归一）。初始节点快照存 `CzmLayout.node_ref`（首次更新时惰性拷贝，永不重置）。输出键 `Delta core [-]`、`core wrinkle amplitude [m]`（写入 variables，CallModel 同站点）+ `collapse_approx = "sp_perfect_bond_phi_free"`（geo_nl 开启时输出；**用户修订：4' 不做**，SP–电极层内完美粘结、Φ 缝自由，至 Batch 8 重估）。
**D-B5-2（已撤回，用户批准 2026-08-22：完整 Crisfield 柱面弧长）**：解锁 `arc_length + geo_nl`，实现 Theory §6.10 完整列式：λ 缩放本征应变增量（`ε*(λ) = ε_ref + λ·Δε*`，`ε_ref` 快照于 `CzmLayout.eigenstrain_prev`）；增广 Newton `Δu = Δu_r + Δλ·Δu_t`（`Δu_r = −K⁻¹R`、`Δu_t = K⁻¹f̂`，`f̂ = ∫BᵀC·Δε* dA` 等效载荷）；柱面约束 `ΔuᵀΔu = Δl²` → Δλ 二次方程 (6.90–6.91)，取与上步切向内积为正的根，无实根 Δl/2 重试；越过极值点（λ 可降）推进至 λ≥1 提交；步长自适应，`step/128` 下限报错终止（不伪造收敛）；**损伤 max-history 与弧长回溯的交互按 §6.10 语义**（D 不回退、重试不幂等，诊断记录重试次数）。撤回理由：D13 探针实测预应力参考态近临界（μ≈5e-3），C4-lite 目标现象即极值点——λ-斜坡会在临界点停滞并阻断 C4-lite。
**D-B5-3（多圈持久化审计而非新建）**：`solve_phase` 已跨相位携带 `czm_mesh`（u_prev/damage/plastic 均在其上与 `czm_layout`），不重置完美圆（旧 bug 已修，代码注释在案）；本批补：`node_ref`/`eigenstrain_prev` 快照字段 + 审计测试锁死"跨相位/跨圈不重置"契约。
**D-B5-4（C4-lite 判据与 D10 分支）**：判据 = 多圈（≥3 相位）Δ_core 与 D_max 联合增长（允许物理合理的非单调，取趋势断言）；输出必须携带 `collapse_approx` 声明。不可达 → D10：预应力量级（0.5/1×）× 本征应变幅值（0.5/1/2×，≤6 组）有界探针 → 停止评审记录结论，不降级验收门。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/czm.jl` | 修改 | `core_ovalization`；geo 分支 eigenstrain 可混参（ε_ref+λΔε*） |
| `src/CzmSolve.jl` | 修改 | 解锁 arc_length+geo；λ-本征应变斜坡（去 F_tc 外载） |
| `src/CouplingState.jl` | 修改 | `CzmLayout.node_ref/eigenstrain_prev`；update_czm_damage!：快照、Δ_core 计算与输出、collapse_approx 标记 |
| `src/CallModel.jl` | 修改 | variables 写 `Delta core [-]`/`core wrinkle amplitude [m]` |
| `test/test_czm_delta_core.jl` | 新建 | 滤波单测（构造含 0/1 阶污染的已知场）+ r_ref 归一 + 持久化审计 |
| `test/test_czm_arc_geo.jl` | 新建 | 弧长 geo 收敛 ≈ basic（弹性区）、λ→1、越界步长报错 |
| `test/test_czm_multicycle_c4lite.jl` | 新建 | 多相位 Δ_core/D 联合增长（C4-lite）或 D10 探针分支 |

---

## Task 1: Δ_core 与多圈快照（TDD）

测试（`test_czm_delta_core.jl`，4 testset）：① 滤波精确性：构造 `u_n(θ) = 0.3+0.5cos θ+0.2cos 2θ+0.1cos 5θ`（放大 1e-4）位移场 → max|ũ_n| ≥ 0.25×max|u_n| 且纯 0/1 阶场 → w≈0（机器零级）；② Δ_core 归一（÷a）；③ `node_ref` 恒定：两次 update 后 Δ_core 基准不变（相对初始螺旋）；④ 持久化审计：跨两次相位调用 u_prev/damage_states 连续（第二次初值 == 第一次终值）。
实现：`core_ovalization`（第一匝窗口 DFT，同探针法）；`CzmLayout.node_ref/eigenstrain_prev` 两新字段；`update_czm_damage!` 尾部计算并写 variables；`CallModel` 键透传 + `collapse_approx = "sp_perfect_bond_phi_free"`（geo_nl 开启时）。
提交 `feat(czm): Φ 缝默认完美粘结（节点合并入构造，v1.5）+ Δ_core 位移基计算与多圈快照`。

测试（`test_czm_delta_core.jl`，4 testset）：① 滤波精确性：构造 `u_n(θ) = 0.3+0.5cos θ+0.2cos 2θ+0.1cos 5θ`（放大 1e-4）位移场 → `w_core ≈ 0.2+0.1`（幅值）量级断言（max|ũ| ≥ 0.25·max|u_n|，且纯 0/1 阶场 → w≈0 机器零）；② Δ_core 归一（÷a）；③ `node_ref` 持久：两次 update_czm_damage! 后 Δ_core 基准不变（相对初始螺旋）；④ 持久化审计：跨两次 solve_phase 调用 u_prev/damage_states 连续（不重置断言：第二次初值 == 第一次终值）。
实现：`core_ovalization`（单匝窗口 DFT 同探针法）；`CzmLayout` 两新字段；`update_czm_damage!` 尾部计算并写 variables；`CallModel` 键透传。
提交 `feat(czm): Φ 缝默认完美粘结（节点合并入构造，v1.5）+ Δ_core 位移基计算与多圈快照`。

## Task 2: 弧长 geo 解锁——完整 Crisfield 柱面弧长（TDD）

测试（`test_czm_arc_geo.jl`，3 testset）：① `arc_length + geo_nl + eigenstrain`（小载荷）收敛且位移 ≈ basic 方法解（rtol 1e-6，弹性区）；② λ 达 1（`result.iterations` 有限、残差 < tol）；③ 越极值点能力：构造软化解（人工降低 σ_max 使 cohesive 早期软化）验证 λ 中途下降仍推进至 λ≥1 并提交失稳后构型。④ `geo_nl=false + arc_length` 行为不变（回归锚）。
实现：`solve_czm_arc_length_step` geo 分支：F_tc 不装配；`eigenstrain_mix(t) = (α,βn,βp, ε_ref + t·(ε_tot−ε_ref))` 逐子步混参装配（新增 `mix_eigenstrain` 帮助函数）；目标 R = F_ext − f_int；收敛/步长/终止沿用现结构（§6.10 终止约定已在）；dispatch 删除 geo_nl 拒绝分支（plasticity/prestress/phi_bond 透传同 load_substep 模式）。`update_czm_damage!`：`eigenstrain_prev` 快照（几何无关，首次建立后每相位更新）。
提交 `feat(czm): 弧长 geo 路径解锁——完整 Crisfield 柱面弧长（§6.10，可越极值点）`。

## Task 3: C4-lite 多圈集成 + 三道门禁（v2）

- [ ] `test_czm_multicycle_c4lite.jl`：nθ=8 夹具，`geo+plastic+prestress(0.2×)` 全开（Φ 默认粘结，v1.5），模拟 3 相位（Δsoc −0.3 → 静置 → 再放电幅值递增 20%）：断言 Δ_core 相位间总体增长（允许一次回落 ≤10%）、D_max 不减、输出含 `collapse_approx = "phi_perfect_bond"`。**不可达即转 D10**：探针（预应力 0.5/1× × 本征应变 0.5/1/2×）打印矩阵 → 测试改 `@test_broken` + findings 记录，**不降级**。
- [ ] 全套；探针与 testexample **v3 重冻结**（Φ 粘结拓扑变化 + 决策声明；v1/v2/v3 档案同存）。
- [ ] `Simplify/baseline.md` 行、progress（C4-lite 证据或 D10 结论）、index；spec 若 C4-lite 达成 → §7 Batch 5 行勾稽注记（v1.4.1 小修，不动决策表）。
- [ ] 提交 `feat/test(czm): C4-lite 多圈集成与门禁记录`。

---

## 自评审记录

1. **spec 覆盖**：§3.5 Δ_core/输出键/多圈持久化 → Task 1；§3.5 路径跟踪 + §4.1 CzmSolve 行 → Task 2（含 D-B5-2 偏离）；§7 Batch 5 三测试（滤波单测/持久化/弧长回归）→ Task 1①/1④/2③；C4-lite/D10 → Task 3。**§3.4（4'）按用户修订不实现**——spec §9 序列留空，C4-lite 声明如实标注 Φ 自由。
2. **决策登记**：**用户二次修订（spec v1.5）：Φ 缝默认完美粘结，4' 取消，基线 v3 重冻结**；D-B5-1（Δ_core 窗口=第一匝、基准=初始螺旋恒定）；D-B5-2（**原偏离已撤回**——用户批准完整 Crisfield，理由：预应力近临界 + C4-lite 目标即极值点）；D-B5-3（持久化审计不重建）；D-B5-4（C4-lite 判据与 D10 分支）。
3. **类型一致性**：`core_ovalization(czm_mesh,u,ref_node)`、`mix_eigenstrain(ref,tot,t)` 命名一致；新键名与 spec §3.5 逐字（collapse_approx 除外——用户修订后的近似声明）。
4. **风险**：C4-lite 可达性未知（D10 兜底已定义）；弧长在小载荷弹性区应与 basic 一致（rtol 1e-6 偏紧则放宽至 1e-4 并记录）；32 个测试文件的运行时长增长。
