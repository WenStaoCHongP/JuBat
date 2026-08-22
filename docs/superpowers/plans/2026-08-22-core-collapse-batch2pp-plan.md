# 堆芯塌陷力学建模 Batch 2''（D13 网格探针）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐薄层径向细分网格能力（opt-in），新建 `tools/czm_mesh_probe.jl` 执行 D13 四组扫描（`nθ∈{80,360}` × 细分 `{1,3}`），以临界载荷因子（Cholesky 二分）与主模态的网格敏感性落定宣称边界，并把结论写回 spec（v1.4）。

**Architecture:** 三任务。Task 1 径向细分（`build_czm_submesh`/`jellyroll_collector_seed_mesh` 增 `thin_subdiv::Int=1`——径向层表展开，SP/PCC/NCC 各分 N 等厚子层；同材内部边界被 `create_czm_mesh` 的材料突变识别自动豁免）+ 切线分离装配（`gl_element`/`assemble_bulk_residual_tangent` 增 `split_KG` 可选返回，默认路径逐位不动）+ 网格拓扑测试。Task 2 探针脚本：预应力参考态下 `K_mat`/`K_G` 分离装配 → Cholesky 二分求失稳载荷因子 μ_crit + 逆幂迭代主模态 + Δ_core-lite（spec §3.5 D8 的探针本地实现）→ `output/czm_mesh_probe/` 输出。Task 3 四组扫描 + D13 判定 + spec v1.4 宣称边界 + 三道门禁（v2 基线）。

**Tech Stack:** Julia 1.11.2；sparse `\` + `cholesky`（无 Arpack 依赖，二分+逆幂足够）；CSV/Printf。

**Spec:** spec v1.3 §3.8（D13）、§3.5（Δ_core 定义）、§7 Batch 2'' 行；代码事实：`build_czm_submesh`（`src/Jellyrollmodel.jl:601-687`）径向层表结构、`create_czm_mesh` 界面按材料突变识别（`src/CzmMesh.jl:50-83`）。

## Global Constraints

（同前：运行环境；基线 **v2**（探针 v2 表、testexample PNG `b31ffb49…`）；默认零漂移——`thin_subdiv=1`/`split_KG=false` 为默认；AGENTS 9.7；每 Task 一提交；§9.9 输出到 `output/czm_mesh_probe/`。）

### 设计决策（执行契约）

**D-B2''-1（细分的层表展开法）**：细分只改 `build_czm_submesh` 的 `layer_thicknesses/material_sequence` 输入表（薄层 `:SP/:PCC/:NCC` 各拆 `thin_subdiv` 个等厚子层、继承材料类型），下游（节点螺旋、单元循环、φ 配对、winding_turn、thermal_elem_map）对 `n_layers` 泛化——`n_layers=8` 硬编码改 `length(layer_thicknesses)`。界面识别按材料突变 ⟹ 细分不产生新 cohesive/Φ 拓扑。cohesive 计数恒为 `4·n_segments`。
**D-B2''-2（切线分离的位级安全）**：`gl_element_residual_tangent(...; split_KG=false)`——`true` 时返回 `(f_e, K_mat_e, K_G_e)` 三元组（分别累积，与现和值逐位一致由测试断言）；`false` 路径原样（bitwise）。`assemble_bulk_residual_tangent` 透传，`split_KG=true` 时返回 `(f, K_mat, K_G)`（仅 geo_nl 路径支持，线弹性路径传 true 即 error）。
**D-B2''-3（临界载荷因子的 Cholesky 二分）**：预应力参考态（Batch 2' σ₀ 场，u=0）下装配 `(K_mat, K_G)`，BC 缩减后二分 μ：`cholesky(K_mat − μ·K_G)` 成功=PD。μ_crit = PD 丧失点（精度：区间收敛至 1e-3 相对）。主模态 = μ_crit−ε 处逆幂迭代（≤50 步）。四组共用同一 σ₀ 场函数与参考定义，仅网格不同。
**D-B2''-4（Δ_core-lite 探针本地实现）**：按 spec §3.5/D8 在探针内实现（Γ_in,free 第一匝内边界节点的 u_n 相对初始螺旋、Fourier 去 0/1 阶取 max 残差、除 r_ref）——Batch 5 生产实现仍按其批次交付，探针版仅供敏感性对比。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/Jellyrollmodel.jl` | 修改 | 两函数增 `thin_subdiv=1`；层表展开；n_layers 泛化 |
| `src/czm.jl` | 修改 | `gl_element`/`assemble_bulk_residual_tangent` 增 `split_KG` |
| `test/test_czm_thin_subdiv.jl` | 新建 | 细分拓扑（子层单元数/材料继承/cohesive 恒 4·n_seg/Φ 重合断言/默认关逐位）+ split_KG 位级一致 |
| `tools/czm_mesh_probe.jl` | 新建 | D13 探针（四组扫描、μ_crit、主模态、Δ_core-lite；输出 `output/czm_mesh_probe/`） |
| spec v1.4 / planning 三件套 | 修改 | D13 判定与宣称边界落定 |

---

## Task 1: 径向细分与切线分离（TDD）

- [ ] Step 1: 失败测试 `test/test_czm_thin_subdiv.jl`（4 testset）：① `thin_subdiv=3` 时 bulk 单元数 = `(5+3·3)·n_seg`（5 厚层+9 薄子层）、材料类型逐层继承、cohesive 恒 `4·n_seg`、Φ 配对坐标重合；② `thin_subdiv=1`（默认）与现状网格逐位（node/element/material_type `==`）；③ `split_KG=true` 返回 `(K_mat+K_G) == K_total`（同一 u 与 σ₀，geo 路径，isbits 逐位）；④ 线弹性路径传 `split_KG=true` → error。
- [ ] Step 2: 确认失败（无 `thin_subdiv`/`split_KG` 关键字）。
- [ ] Step 3: 实现（层表展开 ~10 行 + n_layers 泛化 + split_KG 双累积分支 ~15 行）。
- [ ] Step 4: 4/4 通过 + 回归（`test_czm_mech_core`/`geo_c1`/`geometric_stiffness`/`winding_prestress`/`j2_integration` 零失败）+ 提交 `feat(mesh/czm): 薄层径向细分 thin_subdiv 与切线分离装配 split_KG（默认零漂移）`。

## Task 2: D13 探针脚本

- [ ] Step 1: `tools/czm_mesh_probe.jl`：四组 `(nθ, thin_subdiv) ∈ {80,360}×{1,3}`；每组：归一化网格 → σ₀ 场（Batch 2'）→ `split_KG` 装配（u=0, geo_nl=true, prestress）→ BC 缩减（`extract_bc_dofs`）→ Cholesky 二分 μ_crit → 主模态（逆幂）→ 模态周向阶数 n（对 Γ_in 节点模态位移做离散 Fourier 主峰识别）→ Δ_core-lite（在 μ=0.5·μ_crit 载荷下线性求解 K_mat·u = λ·(−K_G·u_ref)? **简化**：直接用主模态幅值归一化对比，Δ_core-lite 仅在 μ=0.5 μ_crit 线性响应下计算一次）→ 行输出。
- [ ] Step 2: 运行四组（预计每组 1–6 min），落盘 `output/czm_mesh_probe/probe_results.csv`（列：nθ, subdiv, ndof, μ_crit, mode_n, mode_amplitude_ratio, delta_core_lite, seconds）+ 控制台汇总表。
- [ ] Step 3: 提交 `feat(tools): D13 网格探针（μ_crit Cholesky 二分 + 主模态 + Δ_core-lite，四组扫描）`。

## Task 3: D13 判定、spec v1.4 与门禁

- [ ] Step 1: 判定（按 spec §3.8）：μ_crit 与 mode_n 四组间相对变化 <10% 且 mode_n 一致（预计 n=2 椭圆化）→ "接受叠层级边界"；否则 → 单元技术改造立项建议。结论与数据表写入 findings/progress。
- [ ] Step 2: spec v1.4：§3.8 增"探针实测结论"小节（数据+判定+宣称边界语句），头部版本行。
- [ ] Step 3: 三道门禁（v2）：全套 30/30；探针 v2 快照逐位；testexample v2（电/热指标+PNG `b31ffb49…`）。
- [ ] Step 4: `Simplify/baseline.md` 行 + progress/index 更新 + 提交 `docs(spec/probe): D13 网格探针结论与宣称边界落定（spec v1.4）`。

---

## 自评审记录

1. **spec 覆盖**：§3.8 全部（细分能力/4 组/μ_crit/Δ_core/模态/判定）→ Task 1/2/3；§7 Batch 2'' 验收门（网格敏感性结论+宣称边界）→ Task 3 Step 1-2。
2. **决策登记**：D-B2''-1（层表展开，下游泛化）；D-B2''-2（split 位级安全）；D-B2''-3（Cholesky 二分——无外部特征值依赖，PD 丧失点即线性屈曲因子，预应力参考态固定）；D-B2''-4（Δ_core-lite 探针本地，生产版 Batch 5）。
3. **类型一致性**：`thin_subdiv::Int=1`/`split_KG::Bool=false` 命名与既有 kwarg 风格一致；探针仅消费公开接口。
4. **风险**：nθ=360 网格 ndof ~20 万，稀疏分解二分 ~40 次因子分解（每次秒级）——单组上限 ~6 min 可接受；若 μ_crit 不存在（K_G 全稳定，预应力以拉为主）→ 探针改用"特征应变压缩参考态"（ε₀=−0.01 全场）作为 K_G 源，脚本内两参考态都跑并记录（判定以压缩参考态为准——屈曲物理以压缩为准）。
5. **无占位符**：层表展开与 split 双累积的实现规格完整；探针算法（二分区间 [0, μ_max·] 初始倍增、逆幂 50 步、Fourier 主峰）均有明确算法描述与输出契约。
