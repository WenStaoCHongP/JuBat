# 堆芯塌陷力学建模 Batch 0''+1 实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 消解 Theory 中 8 项阻塞实现的内部矛盾，并在 `src/czm.jl` 建立 bulk 残差/切线统一入口 `assemble_bulk_residual_tangent`，使后续几何非线性、塑性、预应力批次有唯一理论参考和唯一装配入口，且开关全关时行为逐位不变。

**Architecture:** 两段式。前段（Batch 0''）是纯文档批次：按代码实参重算 Theory 几何参数表并标注字段来源，统一层编号与 `A_eff` 量纲约定，消除 `07` 的同号小节、`K_uu` 双定义，补齐弧长法与 Φ 约束施加方式说明，加 `κ_ss` 与 C⁰ Q4 不兼容的实现注记。后段（Batch 1）在 `src/czm.jl` 插入一个只有线弹性槽位的 `assemble_bulk_residual_tangent`，把 `assemble_coupled_system` 的 `K_bulk*u` 计算改走该入口；`geo_nl`/`plasticity`/`mech_state` 三个未实现槽位传非默认值即 `error`，不静默降级。

**Tech Stack:** Julia 1.11.2；`Test` 标准库；`SparseArrays`/`LinearAlgebra`；Markdown + LaTeX（Theory 文档）；ripgrep（文档门禁校验）。

## Global Constraints

以下约束逐字来自 spec `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md` 与 `AGENTS.md`，每个 Task 的验收隐含包含本节。

- **运行环境**：Julia 1.11.2，单线程（`JULIA_NUM_THREADS=1`），`GKSwstype=100`，`--startup-file=no`。
- **强制行为基线**（AGENTS 9.6）：入口 `example/testexample.jl`，命令 `& 'D:\Julia-1.11.2\bin\julia.exe' --startup-file=no example\testexample.jl`；必须 exit code 0；网格/步数与 `Simplify/baseline/testexample/README.md` 冻结表全部科学结果在脚本打印精度下完全一致；`output/testexample/testexample_results.png`（a2caecc 起 testexample 按 AGENTS §9.9 输出到 `output/<脚本名>/` 子目录，路径见 README:42）的 SHA-256 应与 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932` 一致。任一科学指标不一致 → 停止并回退当批。
- **兼容性契约**（spec §5 v1.3）：`czm_geo_nonlinear`/`czm_winding_prestress`/`czm_j2_plasticity`/`czm_phi_bond`/`czm_continuous_feedback` 全 false 时，行为、结果键、`tools/verify_czm_standalone.jl` 三方法快照与现状逐指标一致。
- **错误处理**（spec §6，继承 AGENTS 9.7）：缺参 → `error` 指明材料层；不默认、不置零、不截断；不引入任何新的静默回退分支；未实现的槽位传入非默认值一律 `error`，不得静默忽略。
- **网格约定**（AGENTS 9.3）：`nθ` 同时控制热网格与力学/CZM 网格周向分辨率，力学网格直接继承热网格实际角节点，**不得**新增 `nθ_czm` 或独立角度数组。
- **几何数值唯一来源**（D14）：Batch 0'' 之后，Theory 中一切几何数值以 `src/parameters/Jellyroll.jl` 与 `src/Jellyrollmodel.jl` 为唯一来源，表中每行必须标注对应代码字段名。
- **planning 文件归属**（AGENTS 9.5）：进度与发现写入 `docs/planning-with-files/30_堆芯塌陷力学建模/`，并同步更新 `docs/planning-with-files/index.md`。
- **提交粒度**（spec §9）：每批一次提交，沿用现有中文提交风格（`docs(...)`/`refactor(...)`/`test(...)`）。
- **测试套件现状**：`test/runtests.jl` 全套现为 **22/22 通过**（`e117fd2` 已修复 `unit_czm_eigenstrain.jl` 两条既有失败断言，60/60）。spec §7 脚注所载"58/60 既有失败登记"已过期作废；本计划所有批次的全套测试门禁为**全绿**，任何失败都视为新回归。

## 范围与后续计划

本计划只覆盖 spec §9 执行顺序中的 `0''` 与 `1`：

```
0'(done) → [本计划: 0''(理论修订, 8 项阻塞) → 1(力学核边界)] → 2(C1) → 2'(预应力) → 2''(网格探针) → 3(C2) → 4' → 5(C4-lite) → 6 → 7 → [8 后置]
```

后续批次单独立计划，原因是其内容尚未唯一确定，现在写入即为占位符：

- Batch 2/2''：单元技术是否需要改造取决于 D13 网格探针结论（spec §3.8）。
- Batch 6：残差中是否写入 `(1−D)⁻²` 取决于 D12 三选一的专项推导结论（spec §3.6）。
- Batch 7：并入 5 项文档级修订，须在 Batch 0''–6 的实际修订结果之后才能定稿。

### 与 spec §4.1/§4.2 的批次归属与签名偏差（已记录，供评审）

1. **`src/Option.jl` 六个子选项字段在本计划 Batch 1 全部加入**（与 §4.1 一致）。理由：§5 的兼容性契约"全 false 时逐指标一致"需要字段存在才可测试，Task 5 为此写了显式的默认值契约测试。
2. **`PlasticState`、`MechHistory` 类型与 `CZMAssemblyCache` 的参考构型/机械状态字段不在本计划实现**（§4.1 将其标为 Batch 1-3 / 1,3,5,6）。理由：Batch 1 无消费者，提前定义即死类型，违反 YAGNI 且无法测试。`assemble_bulk_residual_tangent` 的 `mech_state` 参数按 §4.2 签名保留为尾置可选位置参数，Batch 3 引入 `PlasticState` 时无需改签名。
3. **`assemble_bulk_residual_tangent` 的关键字参数比 spec §4.2 冻结签名多一个 `K_bulk_cached`**（v1.1 登记）。理由：Batch 1 接线必须在不重复装配的前提下透传既有 `assemble_coupled_system` 的 `K_bulk_cached` 快路径（零漂移）；§4.2 标注"内部，不进公共文档"，故以实现签名为准。Batch 2/3 实现者以本计划签名为准（4 位置参数 + 3 关键字参数）。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `tools/verify_czm_standalone.jl` | 不修改（直接运行） | Batch 1 快照门禁来源：三方法 × 8 载荷水平收敛对比（方案 B；原 `czm_baseline_probe.jl` 已被 `2bf2ac7` 删除，v1.1 修订，见 Task 1） |
| `docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md` | 已建立 | 冻结 `verify_czm_standalone.jl` 快照（2026-08-21，HEAD `e117fd2`），Batch 1 及后续批次的回归基准 |
| `tools/theory_geometry_recompute.jl` | 新建 | 只读脚本：从 `ChooseCell("Jellyroll")` 重算螺旋几何表并打印代码字段来源（D14 防漂移） |
| `Theory/07_弱形式与求解.md` | 修改 | 消除同号 §6.4.5 与 KKT 残留、统一 `K_uu`、新增 §6.10 柱面弧长法、补 Φ 约束施加方式 |
| `Theory/02_几何与运动学.md` | 修改 | 几何数值按代码重算（`46.6\|0.132` 15 处 + `284` 7 处 + `\approx 6` 5 处）、γ 量级联动订正、`κ_ss` 与 C⁰ Q4 不兼容实现注记 |
| `Theory/01_符号与守恒律公理.md` | 修改 | 几何数值（`46.6\|0.132` 7 处 + `284` 3 处）、层编号双约定声明、`A_eff` 量纲 |
| `Theory/03_本构理论.md` | 修改 | 式 (2.68) `A_eff` 改为无量纲面积分数 |
| `Theory/04_CZM.md` | 修改 | 几何数值（仅 `284` 2 处：`:34`、`:49`；无 `46.6/0.132` 命中） |
| `Theory/09_附录A_符号表.md` | 修改 | `A_eff` 量纲与关系式同步 |
| `src/Option.jl` | 修改 | 新增 6 个默认关子选项字段 |
| `src/czm.jl` | 修改 | 新增 `assemble_bulk_residual_tangent`；`assemble_coupled_system` 改走该入口 |
| `test/test_czm_option_defaults.jl` | 新建 | 子选项默认值契约（spec §5） |
| `test/test_czm_mech_core.jl` | 新建 | 新入口 ≡ `K_bulk*u`、缓存不变量、未实现槽位报错（spec §7 Batch 1） |

`test/runtests.jl` 自动发现 `test/` 下所有 `.jl` 文件并以独立进程运行，新增测试文件无需注册。

---

## Task 1: 复核 CZM 基线快照（verify_czm_standalone，方案 B）

**v1.1 修订背景**：本 Task 原为"修复 `czm_baseline_probe.jl` 的 4 处过期 API 并冻结快照"。计划评审（2026-08-21）对照 `49ec452` 历史版本核实该 4 处诊断属实，但该探针已被 `2bf2ac7` 作为过时脚本**整体删除**，修复对象不存在。用户决策采用**方案 B**：Batch 1 快照门禁改用 `tools/verify_czm_standalone.jl`（同覆盖 basic/load_substep/arc_length 三方法 × 8 个 Δsoc_n 载荷水平；`2bf2ac7` 已修复其 `build_czm_cache` 签名并实测通过），不恢复旧探针。spec 已随 v1.3 同步替换 §5/§7 引用。

基线快照 `docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md` 已于 2026-08-21 冻结（HEAD `e117fd2`，计划 v1.1 修订时执行）。本 Task 在开工时**复核其可复现性**，不改任何代码。

**Files:**
- 无修改。只读运行 `tools/verify_czm_standalone.jl`，对照既有快照 `docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md`。

**Interfaces:**
- Consumes: `JuBat.solve_czm_step(czm_mesh, F_ext, param_cache, param, u_prev; α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem, max_iter, tol, n_load_steps, iter_method, cache)`（`src/CzmSolve.jl:651`）；`JuBat.get_damage_statistics(czm_mesh)`（`src/CzmPostProcess.jl:12`，返回 `max_D`/`mean_D`/`n_fractured`）。
- Produces: 无新文件；复核结论记入 progress.md。Task 7 Step 6 与后续所有批次按该快照比对。

- [ ] **Step 1: 重跑工具，逐位比对冻结快照**

Run（与冻结时同环境）: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. tools/verify_czm_standalone.jl`

Expected: 退出码 0。输出中（i）网格统计行（Nodes/Bulk/Cohesive）、（ii）有效参数行（E_eff/ν_eff/α_eff/β_n/β_p）、（iii）8 行收敛对比表的每个 `OK/FAIL it D r` 字段、（iv）Summary 各方法行——全部与 `baseline_czm_standalone.md` 冻结表**在打印精度下逐位一致**。

两点纪律（与冻结快照时的记录一致）：

1. **`FAIL` 条目是冻结行为**。`2bf2ac7` 实测为三方法 7/8 载荷水平收敛；不收敛的条目原样冻结，**不得**调参（改 tol/max_iter/n_load_steps）使其收敛后当作通过。若冻结表存在 FAIL 行而重跑变 OK（或反向），同样视为行为漂移，须停下定位。
2. 工具在 `tools/verify_czm_standalone.jl:134` 重建网格时传 `param_dim`（与 `:66` 模板用的 `case.param` 不一致），且 `:66/:134` 用未合并的 `mesh_data.thermal2D`——这是工具自身的既有状态，**本 Task 不修**（改它会移动基线、混淆批次归因）；已登记为 findings 待办，见 Step 2。CZM 求解只用 bulk/cohesive 拓扑，不受 thermal 网格合并影响，门禁有效性不因此受损。

- [ ] **Step 2: 确认工具既有瑕疵已在 findings 登记（不改代码）**

该登记已于 2026-08-21（计划 v1.1 修订）写入 `docs/planning-with-files/30_堆芯塌陷力学建模/findings.md`"verify_czm_standalone.jl 既有瑕疵登记"节（`:66`/`:134` 传参不一致、用未合并 `thermal2D`、处置约定）。本步只核对条目仍在；若缺失按 findings 同名节补齐。

- [ ] **Step 3: 复核结论记入 progress.md**

在 `docs/planning-with-files/30_堆芯塌陷力学建模/progress.md` 追加一行复核记录（日期、HEAD、比对结果一致/不一致）。本 Task 不产生代码变更，无独立提交；若 Step 1 比对不一致，**停止**，先定位漂移来源（src 改动 or 工具改动 or 环境），不得进入 Task 2。


---

## Task 2: Theory/07 求解章四项阻塞矛盾

对应 spec §9 Batch 0'' 的 ①②③④。这四项都在 `Theory/07_弱形式与求解.md`，同一评审视角，合为一个 Task。

**Files:**
- Modify: `Theory/07_弱形式与求解.md:156-157`、`:245`、`:551-557`、`:573`、`:677`，并在 §6.9 之后追加 §6.10

**Interfaces:**
- Consumes: 无（纯文档）
- Produces: 唯一的 `K_uu` 定义（式 6.8）、唯一的损伤更新框架表述（位移阈值代数更新）、`§6.10` 柱面弧长法（式 6.88–6.91，供 Batch 5 弧长实现引用）、Φ 约束施加方式说明（供 Batch 4' 引用）

- [ ] **Step 1: 确认两处矛盾现状**

Run: `rg -n "§6\.4\.5" "Theory/07_弱形式与求解.md"`
Expected: 两行同号小节 —— `551:#### §6.4.5 §6.4 小结` 与 `555:#### §6.4.5 本节小结`。

Run: `rg -n "K_\{uu\}=" "Theory/07_弱形式与求解.md"`
Expected: 两行不同定义 —— `156`（含 `K_G(\sigma)`，无 cohesive 项）与 `677`（含 cohesive 项，无 `K_G`）。

- [ ] **Step 2: 删除残留的 KKT 版小结（第 555–557 行）**

删除这三行（含其后空行，保留第 559 行的 `---`）：

```markdown
#### §6.4.5 本节小结

损伤残差 $R_D$ 由 KKT 条件 (6.38) 或其互补问题形式 (6.42)/(6.43) 给出。损伤演化为准静态不可逆过程，由 CZM TSL 的能量释放率 $Y$ 与韧度 $G_c$ 的互补条件决定。Jacobian 分块 $K_{Du}$（式 6.44）、$K_{DD}$（正则化）、$K_{DT}=K_{D c_s}=K_{D\phi}=0$（式 6.45）。
```

保留第 551–553 行的 `#### §6.4.5 §6.4 小结`，并把其标题改为 `#### §6.4.5 本节小结`，使小节号唯一、标题风格与全文一致。

被删段落与保留段落互斥：前者称 `K_DD` 需正则化、损伤由 KKT 互补给出；后者（也是 §6.4 正文 (6.38')/(6.42')/(6.45') 的实际结论）称 `K_DD=1` 恒等、无需正则化。实现必须以后者为准。

- [ ] **Step 3: 清除 §6.5.1 的 KKT 引用（第 573 行）**

把

```markdown
其中 $x\in\{T,c_s,c_e\}$（含时间导数的场量）。位移 $u$ 与电势 $\phi$ 为准静态（无 $\dot u,\dot\phi$），损伤 $D$ 由 KKT 代数约束（无 $\dot D$ 显式，通过互补条件）。
```

替换为

```markdown
其中 $x\in\{T,c_s,c_e\}$（含时间导数的场量）。位移 $u$ 与电势 $\phi$ 为准静态（无 $\dot u,\dot\phi$）；损伤 $D$ 由 §6.4 的位移阈值代数更新给出（max-history 式 (6.38') + 代数残差 (6.42')），无 $\dot D$ 显式，**不使用** KKT 互补条件。
```

**v1.1 扩充（评审发现：另有三处现行框架的 KKT 表述，原 Step 3 未覆盖而 Step 9 门禁会失败，一并清除）**：

`Theory/07:19`（§1.1 未知量表），把

```markdown
| $D,\dot D$，KKT 条件 | §3.3 式 (3.28)–(3.30) | 损伤未知量 $R_D$ |
```

替换为

```markdown
| $D,\dot D$，位移阈值代数更新 | §3.3 式 (3.28)–(3.30) | 损伤未知量 $R_D$ |
```

`Theory/07:80`（§1.1 方程表），把

```markdown
| $R_D=0$ | $D$ | §3.3 式 (3.30) 损伤演化 + KKT | 内聚力一致性 |
```

替换为

```markdown
| $R_D=0$ | $D$ | §3.3 式 (3.30) 损伤演化（max-history 代数更新） | 内聚力一致性 |
```

`Theory/07:1162`（本章小结第 4 条；位于 §6.10 插入点之后，与 Step 7 同区域），把

```markdown
4. **§6.4 损伤演化变分**：KKT 条件（式 6.38）、互补问题形式（式 6.42–6.43）、$K_{Du},K_{DD}$（式 6.44, 6.66）；
```

替换为

```markdown
4. **§6.4 损伤演化变分**：位移阈值代数更新（式 6.38' max-history、式 6.42' 代数残差，v2.3 重写）、$K_{Du},K_{DD}$（式 6.44, 6.66）；
```

- [ ] **Step 4: 统一 `K_uu` 定义（式 6.8 补全三项）**

第 155–158 行，把

```markdown
$$
K_{uu}=\int_\Omega B^T C^{\star} B\,dV+K_G(\sigma),\qquad
C^{\star}\in\{C_{\mathrm{eff}}^{\mathrm{ps}},\,C^{ep}\} \tag{6.8}
$$
```

替换为

```markdown
$$
K_{uu}=\underbrace{\int_\Omega B^T C^{\star} B\,dV}_{\text{材料切线}}
+\underbrace{K_G(S)}_{\text{几何刚度（工况 C）}}
+\underbrace{\int_{\Gamma_{\mathrm{coh}}}N_\delta^T\frac{\partial\mathbf{T}}{\partial\boldsymbol{\delta}}N_\delta\,dS}_{\text{cohesive 切线}},\qquad
C^{\star}\in\{C_{\mathrm{eff}}^{\mathrm{ps}},\,C^{ep}\} \tag{6.8}
$$
```

`K_G` 的自变量由 `\sigma` 改为第二类 PK 应力 `S`，与 spec §3.2 冻结的完全 Green-Lagrange 全 Lagrangian 列式功共轭一致（D9）。

- [ ] **Step 5: 让式 (6.54) 复述而非另立定义**

第 674–679 行，把

```markdown
**$K_{uu}$（力学对角块，式 6.8）**：

$$
K_{uu}=\int_\Omega B^T C_{\mathrm{eff}}^{\mathrm{ps}} B\,dV+\int_{\Gamma_{\mathrm{coh}}}N_\delta^T\frac{\partial\mathbf{T}}{\partial\boldsymbol{\delta}}N_\delta\,dS \tag{6.54}
$$
```

替换为

```markdown
**$K_{uu}$（力学对角块）**：与式 (6.8) 同一定义，此处仅复述，不另立形式：

$$
K_{uu}=\int_\Omega B^T C^{\star} B\,dV+K_G(S)+\int_{\Gamma_{\mathrm{coh}}}N_\delta^T\frac{\partial\mathbf{T}}{\partial\boldsymbol{\delta}}N_\delta\,dS\quad(\text{同式 }6.8) \tag{6.54}
$$

工况 R 取 $K_G=0$、$C^{\star}=C_{\mathrm{eff}}^{\mathrm{ps}}$。
```

保留 `\tag{6.54}` 不动：第 1071 行仍在引用该编号，改号会产生新的悬空引用。第 1071 行把 (6.54) 说成 `C_eff^ps` 的分量来源，这是引用错位，属文档级问题，**本批不修**，按 Step 8 记入 Batch 7 清单。

- [ ] **Step 6: 补 Φ 跨匝约束的施加方式**

第 245 行段落以 `表 7.1 力学内部反馈在工况 C 增加 $\varepsilon^p,\kappa$ 路径（§7）。` 结尾。在该段之后插入一个新段：

```markdown
**Φ 跨匝配对的施加方式（实现约定）**：§1.1.5 的 Φ 伪周期映射 $(s,t_{\mathrm{repeat}})\leftrightarrow(s+L_{\mathrm{turn},j},0)$ **不以附加弱形式项（Lagrange 乘子或罚项）进入式 (6.4)**，而是在离散层通过**节点合并 / 主从自由度消元**施加：配对节点共享同一位移自由度，约束在装配前即被精确满足，故残差与 Jacobian 中不出现对应的乘子块。因此式 (6.4)–(6.8) 与 §6.6 的 Jacobian 分块结构不含 Φ 项，这不是遗漏。合并的前置条件是配对节点坐标重合（`CzmSubmesh.phi_pairs` 两端节点坐标之差在几何容差内），实现侧须在合并前断言该条件，不重合即 `error`，不做就近吸附。该路径为 opt-in（`czm_phi_bond`，默认关）；关闭时两侧为独立自由度，界面自由。
```

- [ ] **Step 7: 新增 §6.10 柱面弧长法**

理论现状：全库只有 `Theory/04_CZM.md:1277` 一句提及弧长，无任何列式。Batch 5 要实现路径跟踪，必须先有唯一参考。

在 §6.9 接口传递小节的末尾、下一个顶级分节之前，插入下面整节。新编号从 (6.88) 起，不与现有最大编号 (6.87) 冲突：

```markdown
### §6.10 路径跟踪：柱面弧长法（Crisfield）

极限点（临界载荷处切线奇异）附近，载荷控制与位移控制的 Newton 迭代均会失败：前者在极限点后无解，后者在多个自由度同时失稳时不唯一。堆芯塌陷的目标现象（内圈不圆度突增）正是极限点/分岔行为，故须引入路径跟踪。

**载荷参数化**：本模型的驱动量不是外力，而是热-化学本征应变。因此弧长参数 $\lambda$ 缩放**本征应变增量**，而非外力向量：

$$
\boldsymbol{\varepsilon}^*(\lambda)=\boldsymbol{\varepsilon}^{*,n}+\lambda\,\Delta\boldsymbol{\varepsilon}^*,\qquad
\Delta\boldsymbol{\varepsilon}^*=\alpha_{\mathrm{eff}}\Delta T+\beta_{\mathrm{eff}}\Delta\bar\theta \tag{6.88}
$$

若同时开启卷绕预应力（自平衡初始应力场 $\sigma_0$），$\sigma_0$ 不参与 $\lambda$ 缩放：它是构型的一部分，而非加载路径的一部分。

**增广方程组**：未知量扩展为 $(\Delta u,\Delta\lambda)$，在 Newton 修正的基础上附加柱面弧长约束：

$$
\begin{cases}
K_{uu}\,\Delta u=-R_u(u,\lambda)+\Delta\lambda\,\hat{f}\\[2pt]
\Delta u^T\Delta u=\Delta l^2
\end{cases} \tag{6.89}
$$

其中 $\hat f=\partial(-R_u)/\partial\lambda=\int_\Omega B^T C^{\star}\Delta\boldsymbol{\varepsilon}^*\,dV$ 是本征应变增量的等效载荷向量，$\Delta l$ 为弧长半径。式 (6.89) 第二式取**柱面**形式（Crisfield，不含 $\Delta\lambda^2$ 项），避免为 $\lambda$ 选择量纲权重。

**预测-修正**：迭代 $k$ 处把修正分解为 $\Delta u^{(k)}=\Delta u_r^{(k)}+\Delta\lambda^{(k)}\Delta u_t^{(k)}$，其中 $\Delta u_r=-K_{uu}^{-1}R_u$（残差修正）、$\Delta u_t=K_{uu}^{-1}\hat f$（切线预测）。代入约束得 $\Delta\lambda$ 的二次方程

$$
a_2(\Delta\lambda)^2+a_1\Delta\lambda+a_0=0 \tag{6.90}
$$

$$
a_2=\Delta u_t^T\Delta u_t,\quad
a_1=2(\bar u+\Delta u_r)^T\Delta u_t,\quad
a_0=(\bar u+\Delta u_r)^T(\bar u+\Delta u_r)-\Delta l^2 \tag{6.91}
$$

$\bar u$ 为本增量步内已累积的位移修正。两根中取与上一步切线方向内积为正者（避免回溯已走过的路径）；实根不存在时按 $\Delta l\leftarrow\Delta l/2$ 重试。

**终止约定**：$\Delta l$ 减至下限仍无实根或仍不收敛时**报错终止并输出诊断**（残差史、$\lambda$ 史、步长史），不外推、不取单根近似、不把未收敛状态作为收敛提交。这与 §6.8 的收敛判据一致：路径跟踪失败是可诊断的物理/离散信息（可能正是分岔点密集），不是需要绕过的数值噪声。

**与损伤更新的次序**：$\lambda$ 只缩放本征应变；损伤 $D$ 仍按 §6.4 的位移阈值代数更新在每次迭代内更新，其 max-history 不可逆性使弧长回溯时 $D$ 不回退，故弧长步长缩减后的重试路径上 $D$ 单调不减——这是物理正确的（损伤不可逆），但也意味着弧长重试不是幂等操作，诊断输出须记录重试次数。
```

- [ ] **Step 8: 记录 Batch 7 待办**

在 `docs/planning-with-files/30_堆芯塌陷力学建模/findings.md` 末尾追加：

```markdown
### Batch 0'' 执行中新增的 Batch 7 文档级待办

- `Theory/07:1071` 把式 (6.54) 说成 $C_{\mathrm{eff}}^{\mathrm{ps}}$ 分量（$C_{nn},C_{ss},C_{sn},C_s$）的来源，而 (6.54) 实际是 $K_{uu}$ 的复述。Batch 0'' 保留 (6.54) 编号以免产生新悬空引用，引用错位本身并入 Batch 7 修订。
```

- [ ] **Step 9: 校验四项已消解**

Run: `rg -n "§6\.4\.5" "Theory/07_弱形式与求解.md"`
Expected: 恰好 1 行。

Run: `rg -n "KKT" "Theory/07_弱形式与求解.md"`
Expected: 仅出现在历史沿革表述中（§6.4 的 v2.3 修订说明块、"与 v2.2 KKT 形式对比"类措辞）。v1.1 清理后 `:19`、`:80`、`:573`、`:1162` 四处现行框架表述均已改为位移阈值代数更新措辞。逐条人工确认每处命中属于历史沿革；出现任何"损伤由 KKT 约束/互补条件给出"的现行表述即未完成。

Run: `rg -n "K_\{uu\}=" "Theory/07_弱形式与求解.md"`
Expected: 2 行（156 与 677 区域），两者三项完全一致，其中一处带 `(\text{同式 }6.8)`。

Run: `rg -n "Crisfield|柱面弧长" "Theory/07_弱形式与求解.md"`
Expected: 至少 2 行命中，位于新增 §6.10。

- [ ] **Step 10: 提交**

```bash
git add "Theory/07_弱形式与求解.md" "docs/planning-with-files/30_堆芯塌陷力学建模/findings.md"
git commit -m "docs(theory): 07 消除同号小节与 K_uu 双定义，补柱面弧长法与 Φ 约束施加方式"
```

---

## Task 3: 几何参数表按代码重算并标注字段来源（D14）

对应 spec §9 Batch 0'' 的 ⑧。Theory 的几何数值疑因把 21700 的直径当半径，导致匝数、弧长、螺距角全线偏离。命中分布（2026-08-21 评审实测，v1.1 订正原"共 38 处"）：`46.6|0.132` 共 22 处（`01`:7、`02`:15、`04`:0）；`284` 共 18 处（`01`:3、`02`:7、`03`:3、`04`:2、`06`:2、`CLAUDE.md`:1，其中 `06:69` 为 ρ≈2846 豁免）；`02` 另有 `\approx 6` 5 处。计数仅供对账，最终以 Step 7 的零命中门禁为准。

**Files:**
- Create: `tools/theory_geometry_recompute.jl`
- Modify: `Theory/02_几何与运动学.md`（`46.6/0.132/6 m` 类命中 15+5 处 + `284` 7 处）、`Theory/01_符号与守恒律公理.md`（`46.6|0.132` 7 处 + `284` 3 处）、`Theory/04_CZM.md`（仅 `284` 2 处：`:34`、`:49`，无 `46.6/0.132` 命中）、`Theory/03_本构理论.md`（`284` 与层厚典型值 3 处：`:44`、`:77`、`:331`）、`Theory/06_热源.md:281`（`284` 1 处）、`Theory/CLAUDE.md:17`（`284` 1 处）

`Theory/03_本构理论.md` 也被 Task 4 修改，但两者改的是不同行（Task 3 改层厚与 `t_repeat`，Task 4 改式 2.68 的 `A_eff`），按顺序执行不冲突。

**Interfaces:**
- Consumes: `JuBat.ChooseCell("Jellyroll")` → `param.cell.layer`、`param.cell.Rin`、`param.cell.Rout`、`param.{PE,NE,SP,PCC,NCC}.thickness`
- Produces: `tools/theory_geometry_recompute.jl` 的打印表，Theory 几何表每行标注的代码字段名与之一一对应；Batch 2''/5 的网格与 Δ_core 讨论直接引用这些值

- [ ] **Step 1: 建立可复现的重算脚本**

创建 `tools/theory_geometry_recompute.jl`：

```julia
# tools/theory_geometry_recompute.jl
#
# 从代码实参重算 Jellyroll 螺旋几何，供 Theory/01,02,04 几何表引用（D14）。
# 只读诊断脚本：不改任何求解路径，不写文件。
#
# 运行方式: julia --startup-file=no --project=. tools/theory_geometry_recompute.jl
#
# 每一行输出都标注对应代码字段名，Theory 表格必须逐行引用该字段名，
# 使几何数值只有一个来源，避免再次漂移。

using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

p = JuBat.ChooseCell("Jellyroll")

layer = p.cell.layer
Rin = p.cell.Rin
Rout = p.cell.Rout
b = layer / (2pi)

# theta 上界与 Jellyrollmodel.jl:40 的 theta1 完全一致（s_in=0, s_out=cell.layer）
theta_end = min((Rout - Rin - layer) / b, (Rout - Rin) / b)
N_turns = theta_end / (2pi)
# 阿基米德螺旋弧长，b << r 下取 ∫ r dθ
L_spiral = Rin * theta_end + b * theta_end^2 / 2
r_avg = L_spiral / theta_end
L_turn = 2pi * r_avg
r_out_mid = Rin + b * theta_end

println("=== 层厚（src/parameters/Jellyroll.jl） ===")
@printf("PE.thickness   = %6.2f um\n", p.PE.thickness * 1e6)
@printf("NE.thickness   = %6.2f um\n", p.NE.thickness * 1e6)
@printf("SP.thickness   = %6.2f um\n", p.SP.thickness * 1e6)
@printf("PCC.thickness  = %6.2f um\n", p.PCC.thickness * 1e6)
@printf("NCC.thickness  = %6.2f um\n", p.NCC.thickness * 1e6)

println("\n=== 螺旋几何 ===")
@printf("cell.layer (t_repeat)      = %.6e m = %.2f um\n", layer, layer * 1e6)
@printf("cell.Rin   (a)             = %.6e m = %.3f mm\n", Rin, Rin * 1e3)
@printf("cell.Rout                  = %.6e m = %.3f mm\n", Rout, Rout * 1e3)
@printf("b = cell.layer/(2pi)       = %.6e m/rad\n", b)
@printf("theta_end (Jellyrollmodel.jl:40) = %.4f rad\n", theta_end)
@printf("N_turns = theta_end/(2pi)  = %.3f\n", N_turns)
@printf("L_spiral = int r dtheta    = %.5f m\n", L_spiral)
@printf("r_avg = L_spiral/theta_end = %.6e m = %.3f mm\n", r_avg, r_avg * 1e3)
@printf("L_turn = 2pi*r_avg         = %.6e m = %.2f mm\n", L_turn, L_turn * 1e3)
@printf("L_turn/L_spiral            = %.3f %% (= 1/N_turns = %.3f %%)\n",
    100 * L_turn / L_spiral, 100 / N_turns)

println("\n=== 螺距角 gamma = t_repeat/(2*pi*r) ===")
for (nm, r) in (("r = Rin  ", Rin), ("r = r_avg", r_avg), ("r = r_out", r_out_mid))
    g = layer / (2pi * r)
    @printf("%s : gamma = %.5e rad = %.4f deg, gamma^2 = %.3e\n", nm, g, g * 180 / pi, g^2)
end

println("\n=== 周向单元长与长宽比（每层厚度方向仅 1 个 Q4）===")
for nth in (40, 80, 360)
    Lc = L_turn / nth
    @printf("n_theta=%3d: elem_len=%7.2f um | SP/NCC %5.1f:1 | PCC %5.1f:1 | PE %4.1f:1 | NE %4.1f:1\n",
        nth, Lc * 1e6, Lc / p.SP.thickness, Lc / p.PCC.thickness,
        Lc / p.PE.thickness, Lc / p.NE.thickness)
end
```

- [ ] **Step 2: 运行脚本，核对权威值**

Run: `julia --startup-file=no --project=. tools/theory_geometry_recompute.jl`

Expected（本计划撰写时实测，用于交叉校验；若不一致说明参数已变更，以脚本输出为准并同步更新本计划）：

```
cell.layer (t_repeat)      = 3.736000e-04 m = 373.60 um
cell.Rin   (a)             = 1.920000e-03 m = 1.920 mm
cell.Rout                  = 1.015000e-02 m = 10.150 mm
b = cell.layer/(2pi)       = 5.946029e-05 m/rad
theta_end                  = 132.1285 rad
N_turns                    = 21.029
L_spiral                   = 0.77271 m
r_avg                      = 5.848200e-03 m = 5.848 mm
L_turn                     = 3.674532e-02 m = 36.75 mm
L_turn/L_spiral            = 4.755 %
gamma: Rin 1.7744 deg / r_avg 0.5825 deg / r_out 0.3485 deg
gamma^2: 9.591e-04 / 1.034e-04 / 3.699e-05
n_theta= 80: elem_len=459.32 um | SP/NCC 38.3:1 | PCC 28.7:1 | PE 6.1:1 | NE 5.4:1
n_theta=360: elem_len=102.07 um | SP/NCC  8.5:1 | PCC  6.4:1 | PE 1.4:1 | NE 1.2:1
```

`n_theta=80` 一行与 spec §3.8 的长宽比表（薄层 29–38:1、涂层约 6:1）一致，可作为脚本正确性的旁证。

- [ ] **Step 3: 建立替换对照表**

按下表在三个 Theory 文件中逐处替换。旧值全部错误，不是精度差异。

| 量 | 旧值（错） | 新值（代码） | 标注的代码字段 |
|---|---|---|---|
| $t_{\mathrm{repeat}}$ | $\approx 284\;\mu$m | $=373.6\;\mu$m | `cell.layer` |
| $t_{\mathrm{PE}}$ | $\approx 60\;\mu$m | $=75.6\;\mu$m | `PE.thickness` |
| $t_{\mathrm{NE}}$ | $\approx 60\;\mu$m | $=85.2\;\mu$m | `NE.thickness` |
| $t_{\mathrm{SP}}$ | $\approx 12\;\mu$m | $=12.0\;\mu$m | `SP.thickness` |
| $t_{\mathrm{PCC}}$ | $\approx 12\;\mu$m | $=16.0\;\mu$m | `PCC.thickness` |
| $t_{\mathrm{NCC}}$ | $\approx 8\;\mu$m | $=12.0\;\mu$m | `NCC.thickness` |
| $a$（内半径） | $\sim 3$ mm | $=1.92$ mm | `cell.Rin` |
| $R_{\mathrm{out}}$ | （未列） | $=10.15$ mm | `cell.Rout` |
| $r_{\mathrm{avg}}$ | $\sim 10$ mm | $=5.848$ mm | 派生 `L_spiral/theta_end` |
| $\theta_{\mathrm{end}}$ | $46.6\cdot 2\pi$ | $=132.13$ rad $=21.03\cdot 2\pi$ | `Jellyrollmodel.jl:40` `theta1` |
| $N_{\mathrm{turns}}$ | $46.6$ | $=21.03$ | 派生 `theta_end/(2π)` |
| $L_{\mathrm{spiral}}$ | $\approx 6$ m | $=0.7727$ m | 派生 $\int r\,d\theta$ |
| $L_{\mathrm{turn}}$ | $\approx 0.132$ m | $=0.03675$ m | 派生 $2\pi r_{\mathrm{avg}}$ |
| 一匝占总弧长 | $\sim 2.2\%$（$1/46.6$） | $=4.755\%$（$1/21.03$） | 派生 |
| $\gamma$ 范围 | $[0.26°,0.86°]$ | $[0.3485°,1.7744°]$ | 派生 $t_{\mathrm{repeat}}/(2\pi r)$ |
| $\mathcal{O}(\gamma^2)$ | $\sim 10^{-5}$ | $\sim 10^{-4}$（最内匝达 $9.6\times10^{-4}$） | 派生 |

**派生比值**（整档改变，必须一并重算，不能只换分子分母数字）：

| 位置 | 旧值（错） | 新值 |
|---|---|---|
| `Theory/01:238` $L_{\mathrm{spiral}}/t_{\mathrm{repeat}}$ | $6\,\mathrm{m}/284\,\mu\mathrm{m}\sim 2\times 10^4\sim O(10^5)$ | $0.7727\,\mathrm{m}/373.6\,\mu\mathrm{m}\approx 2.07\times 10^3\sim O(10^3)$ |
| `Theory/02:338` 弧长 vs SPMe 厚度向 | $6\,\mathrm{m}/100\,\mu\mathrm{m}\sim 6\times 10^4$ | $0.7727\,\mathrm{m}/100\,\mu\mathrm{m}\approx 7.7\times 10^3\sim O(10^4)$ |

`Theory/01:238` 旧文本自身也不自洽（先写 $2\times10^4$ 又写 $O(10^5)$），改写时把量级标注统一为 $O(10^3)$。两处比值虽降一个量级，但"尺度分离成立"的结论不变（$10^3$–$10^4$ 仍是充分分离），不要改结论、只改数值与量级符号。

- [ ] **Step 4: 改 `Theory/02_几何与运动学.md`（`46.6/0.132/6 m` 25 处 + `284` 7 处）**

Run: `rg -n "46\.6|284|0\.132|6\\\\,\\\\mathrm\{m\}|\\\\approx 6" "Theory/02_几何与运动学.md"`

按 Step 3 对照表逐处替换。除数值外，还须处理下列联动：

1. **第 86–91 行几何参数表**：在表头追加一列 `代码字段`，每行填 Step 3 表最后一列；量级列（`\mathcal{O}(\cdot)`）随新值订正（如 $L_{\mathrm{turn}}$ 由 $\mathcal{O}(10^{-1})$ 改为 $\mathcal{O}(10^{-2})$ m）。
2. **第 112–113、118 行 $\gamma$ 推导**：重代入 $t_{\mathrm{repeat}}=373.6\,\mu$m、$a=1.92$ mm、$r_{\mathrm{avg}}=5.848$ mm。第 118 行是 v2.3 的一次订正说明，**不要删除该说明块**，改为二次订正：原 v2.2 称 $\gamma\approx 0.01°$、v2.3 订为 $[0.26°,0.86°]$（当时误用 $a=3$ mm、$r_{\mathrm{avg}}=10$ mm、$t_{\mathrm{repeat}}=284\,\mu$m），本次按代码实参订为 $[0.3485°,1.7744°]$，$\mathcal{O}(\gamma^2)$ 由 $10^{-5}$ 订为 $10^{-4}$。等价性论证结论不变，仅"等价精度"再降一个量级。
3. **第 300 行**：一匝占比 $2.2\%\to 4.755\%$，$1/46.6\to 1/21.03$。
4. **第 321 行**：五个层厚典型值全部按代码值改写，并把 `故 $t_{\mathrm{repeat}}\approx 284\;\mu\mathrm{m}$` 改为 `故 $t_{\mathrm{repeat}}=2(t_{\mathrm{PE}}+t_{\mathrm{NE}}+t_{\mathrm{SP}})+t_{\mathrm{PCC}}+t_{\mathrm{NCC}}=373.6\;\mu\mathrm{m}$（代码字段 \`cell.layer\`）`。
5. **第 447、453、705、707、709、725、768 行**：`46.6\cdot 2\pi` → `21.03\cdot 2\pi`。这些是坐标范围声明，不涉及推导结论。
6. **第 63、72、88、329 行**：`L_{\mathrm{spiral}}\approx 6` m → `=0.7727` m；第 329 行"尺度 $\sim 6$ m"→"尺度 $\sim 0.8$ m"。注意该处论证是"$x$（$\sim 100\,\mu$m）与 $s$ 是两个独立嵌套尺度"，$0.8$ m 与 $100\,\mu$m 仍差 4 个量级，论证成立，不需改结论。

- [ ] **Step 5: 改 `Theory/01`、`Theory/04`**

Run: `rg -n "46\.6|284|0\.132" "Theory/01_符号与守恒律公理.md" "Theory/04_CZM.md"`

按 Step 3 两张表替换。另有两处需要额外处理：

1. `01:18` 的 B 方案坐标系声明：`θ_unfold ∈ [0, 46.6·2π]` → `[0, 21.03·2π]`；`n ∈ [0, t_repeat]` 处补注 `（t_repeat = 373.6 μm，代码字段 cell.layer）`。
2. `01:238` 按"派生比值"表重算，量级符号统一为 $O(10^3)$。

`04:34` 与 `04:49` 的 `t_repeat ≈ 284 μm` 是"1 单元/匝的径向分辨率"论证的一部分：径向分辨率随 `t_repeat` 增大到 373.6 μm，而 `04:49` 的热扩散特征长度 $\ell_T\sim1$–$10$ mm 不变，$\ell_T \gg t_{\mathrm{repeat}}$ 仍成立，**结论不变**，只改数值。

- [ ] **Step 6: 改 `Theory/03`、`Theory/06`、`Theory/CLAUDE.md`**

Run: `rg -n "284" "Theory/03_本构理论.md" "Theory/06_热源.md" "Theory/CLAUDE.md"`
Expected: `03` 3 行（`:44`、`:77`、`:331`）、`06` 1 行（`:281`）、`CLAUDE.md` 1 行（`:17`）。

- `03:44` 与 `03:331` 同时含五个层厚典型值和 `t_repeat`，全部按 Step 3 对照表替换。这两处原文称"与 spec §1.3 一致"/"基于 LGM50 典型值"——改数值后该出处已不成立，把出处改为 `代码字段 PE.thickness/NE.thickness/SP.thickness/PCC.thickness/NCC.thickness（src/parameters/Jellyroll.jl）`。`Jellyroll.jl` 的层厚注释虽标 Chen2020，但与 `03` 原文引用的 LGM50 整理值不同，以代码为准（D14）。
- `03:77` 除把 `≈ 284 μm` 改为 `= 373.6 μm` 外，句中"亦等于 spec §1.3 中的 $t_{\mathrm{repeat}}\approx 284\,\mu\mathrm{m}$"的出处引用随 `:44`/`:331` 同标准改写为代码字段 `cell.layer`（`src/parameters/Jellyroll.jl`）——旧出处数值已失效，不得保留（v1.1 补充）。
- `06:281` 是热正反馈量级估计中的 $L=t_{\mathrm{repeat}}$，改数值；该估计只给量级，$284\to373.6$ μm 不改变量级结论，不要改结论。
- `CLAUDE.md:17` 改为 `n ∈ [0, t_repeat] = 373.6 μm（代码字段 cell.layer）`。

- [ ] **Step 7: 全库门禁校验**

Run: `rg -n "46\.6|0\.132" Theory/`
Expected: 0 命中。

Run: `rg -n "284" Theory/`
Expected: 恰好 1 命中 —— `Theory/06_热源.md:69` 的 `ρ≈2846 kg/m³`（等效密度实测值，与 `t_repeat` 无关）。这是已知豁免项，须在 findings 中记明；出现任何其他命中即视为漏改。

Run: `rg -n "2\\\\times 10\^4|6\\\\times 10\^4" Theory/`
Expected: 0 命中（两处派生比值已重算）。

Run: `rg -n "代码字段" "Theory/02_几何与运动学.md"`
Expected: ≥1 命中，且第 86–91 行几何表每一行都带字段名。

- [ ] **Step 8: 提交**

```bash
git add tools/theory_geometry_recompute.jl "Theory/01_符号与守恒律公理.md" "Theory/02_几何与运动学.md" "Theory/03_本构理论.md" "Theory/04_CZM.md" "Theory/06_热源.md" Theory/CLAUDE.md
git commit -m "docs(theory): 几何参数表按代码实参重算并标注字段来源（21.03 匝/0.7727 m/373.6 um）"
```

---

## Task 4: 层编号、`A_eff` 量纲、`κ_ss` 实现注记

对应 spec §9 Batch 0'' 的 ⑤⑥⑦。三项都是"同一符号两种含义"造成的误读风险，合为一个 Task。

**Files:**
- Modify: `Theory/01_符号与守恒律公理.md:18` 后（层编号声明）、`:649`（`A_eff`）
- Modify: `Theory/03_本构理论.md:1196`（式 2.68）
- Modify: `Theory/09_附录A_符号表.md:63`（`A_eff`）
- Modify: `Theory/02_几何与运动学.md:536` 后（`κ_ss` 注记）

**Interfaces:**
- Consumes: `CzmSubmesh.material_type` 的 5 个 Symbol（`:PE`、`:NE`、`:SP`、`:PCC`、`:NCC`，见 `src/czm.jl:40-45` 的 `moduli_of` 文档）；8 层物理层序（`src/Jellyrollmodel.jl` 的 `n_layers = 8`）
- Produces: 唯一的层编号约定（后续所有 Theory/spec 表述引用）；`A_eff` 唯一为无量纲面积分数（Batch 6 反馈实现直接依据）

- [ ] **Step 1: 确认三处现状**

Run: `rg -n "第 \\$?i\\$? 层|第 i 层" Theory/`
Expected: 多处命中，且无法从上下文判定指的是 5 种材料类型还是 8 层物理层序。

Run: `rg -n "A_\{\\mathrm\{eff\}\}" "Theory/01_符号与守恒律公理.md" "Theory/09_附录A_符号表.md" "Theory/03_本构理论.md"`
Expected: `01:649` 与 `09:63` 单位列为 `m²`；`03:1196` 为 $A_{\mathrm{eff}}(D)=A_0(1-D)$。

- [ ] **Step 2: 插入层编号双约定声明**

`Theory/01_符号与守恒律公理.md`，在第 18 行"B 方案坐标系声明"段之后插入：

```markdown
**层编号双约定声明（Batch 0'' 统一）**：本理论中"层"有两种互不相同的编号，全文出现"第 $i$ 层"时必须按本声明区分，不得混用：

1. **材料类型编号**（5 种）：$i\in\{\mathrm{PE},\mathrm{NE},\mathrm{SP},\mathrm{PCC},\mathrm{NCC}\}$，对应代码 `CzmSubmesh.material_type` 的 5 个 Symbol。本构参数（$E_i,\nu_i,\alpha_i$）按材料类型查表，与其在卷绕重复单元中出现几次无关。
2. **物理层序编号**（8 层）：$i=1,\ldots,8$ 沿厚度 $n$ 方向遍历一个卷绕重复单元，层序为 PE → PCC → PE → SP → NE → NCC → NE → SP，对应代码 `n_layers = 8`。厚度加权、叠层等效（式 2.37–2.38）、逐层积分一律用此编号。

两者不可互推：PE 与 NE 各在 8 层序中出现两次，SP 出现两次，故 $\sum_{i=1}^{8}t_i=t_{\mathrm{repeat}}$ 而 $\sum_{i\in\text{5 类}}t_i\ne t_{\mathrm{repeat}}$。**误读后果**：按材料类型编号把"第 4 层"当作 PCC（5 类顺序的第 4 个），而按物理层序"第 4 层"是 SP，两者本构相差一个数量级。凡涉及厚度求和、逐层积分、叠层等效的表述，一律显式写"物理层序第 $i$ 层（$i=1..8$）"；凡涉及本构查表的表述，一律显式写"材料类型 $i$"。
```

- [ ] **Step 3: 统一 `A_eff` 为无量纲面积分数**

`Theory/01_符号与守恒律公理.md:649`，把

```markdown
| $A_{\mathrm{eff}}(D)$ | 损伤相关有效反应面积 | m² | §4.5, R-EC-1 | $A_{\mathrm{eff}}(D) = A_0(1-D)$ |
```

替换为

```markdown
| $A_{\mathrm{eff}}(D)$ | 损伤相关有效反应面积**分数**（相对无损面积） | —（无量纲） | §4.5, R-EC-1 | $A_{\mathrm{eff}}(D) = 1-D$；绝对面积按需写作 $A_0\,A_{\mathrm{eff}}$ |
```

`Theory/09_附录A_符号表.md:63`，把

```markdown
| $A_{\mathrm{eff}}(D)$ | 损伤相关有效反应面积 | m² | §3.5, 式(3.97) | $A_{\mathrm{eff}}=A_0(1-D)$ |
```

替换为

```markdown
| $A_{\mathrm{eff}}(D)$ | 损伤相关有效反应面积**分数**（相对无损面积） | —（无量纲） | §3.5, 式(3.97) | $A_{\mathrm{eff}}=1-D$；绝对面积按需写作 $A_0\,A_{\mathrm{eff}}$ |
```

`Theory/03_本构理论.md:1196-1198`，把

```markdown
A_{\mathrm{eff}}(D)=A_0(1-D) \tag{2.68}
```

及其后的说明句改为

```markdown
A_{\mathrm{eff}}(D)=1-D \tag{2.68}
```

说明句改为：

```markdown
其中 $A_{\mathrm{eff}}$ 是**无量纲**有效面积分数（相对无损面积 $A_0$ 归一，$A_{\mathrm{eff}}(0)=1$）；需要绝对面积时显式写 $A_0\,A_{\mathrm{eff}}$。$D\uparrow$ → $A_{\mathrm{eff}}\downarrow$ → 电化学反应总面积 $Q\downarrow$（标量，进入 SPMe 源项）与界面热导 $k_n^{\mathrm{eff}}\downarrow$。
```

量纲自查：$\eta_{\mathrm{eff}}=\eta-j\,R_{\mathrm{contact}}/A_{\mathrm{eff}}$ 中 $j$ 为界面反应电流密度（A/m²）、$R_{\mathrm{contact}}$ 为面电阻率（Ω·m²）、$A_{\mathrm{eff}}$ 无量纲 → 结果为 V ✓。这与 spec §3.6 的实现约定一致。

**注意**：本步只统一量纲约定，**不触碰** $R_{\mathrm{contact}}=R_0/(1-D)$ 与 $A_{\mathrm{eff}}=1-D$ 联合产生的 $(1-D)^{-2}$ 是否双计——那是 D12，延后至 Batch 6 专项推导，本批不得预先选定任一候选。

- [ ] **Step 4: 加 `κ_ss` 与 C⁰ Q4 不兼容的实现注记**

`Theory/02_几何与运动学.md`，在第 536 行段落（以"——对应内圈法向皱褶模态。"结尾）之后插入：

```markdown
> **实现注记（D9，Batch 0''）**：式 (1.30a) 的曲率 $\kappa_{ss}=-\partial^2 u_n/\partial s^2$ 含位移二阶导，要求 $C^1$ 连续的插值（Kirchhoff 板/梁类单元）。JuBat 的力学离散是**双线性 $C^0$ Q4 平面单元**，节点间只保证位移连续、不保证转角连续，$\partial^2 u_n/\partial s^2$ 在单元内恒为常数、跨单元不连续，因此 **$\kappa_{ss}$ 不可作为实现公式**。
>
> 实现侧改用**物理 $(x,y)$ 坐标上的完全 Green-Lagrange 全 Lagrangian 列式**：$E=\tfrac12(F^TF-I)$ 取全部二次项，本构为逐层各向同性 St. Venant–Kirchhoff（与第二类 PK 应力 $S$ 功共轭），几何刚度取标准初应力形式 $K_G(S)=\int G^T\hat S\,G\,dV$。弯曲效应不通过显式曲率自由度引入，而由**物理坐标网格几何本身**携带（每层厚度方向的 Q4 在螺旋上天然弯曲）；逐层各向同性使全程无需材料局部标架。
>
> 因此式 (1.30a) 及其 $\kappa_{ss}$ 在本文中**降格为物理动机引用**：它说明中等转动下膜–弯耦合为何产生几何刚度与分岔，但不参与离散。完全 GL 在中等转动极限下涵盖 (1.30a)，且对低阶周向模态（$n=2$ 椭圆化，即 $\Delta_{\mathrm{core}}$ 主模态）无 Donnell 型截断误差。
>
> **能力边界**：每层厚度方向仅 1 个全积分双线性 Q4，薄层（SP/NCC 12 μm、PCC 16 μm）在 $n_\theta=80$ 时长宽比达 28.7–38.3:1，剪切闭锁会显著高估弯曲刚度。故本列式可解叠层级低阶不圆度，**不足以解析单层起皱**；宣称边界由 Batch 2'' 网格探针（spec §3.8）确定。
```

- [ ] **Step 5: 校验三项已消解**

Run: `rg -n "层编号双约定声明" "Theory/01_符号与守恒律公理.md"`
Expected: 1 命中。

Run: `rg -n "A_\{\\mathrm\{eff\}\}.{0,80}m²" Theory/`
Expected: 0 命中。

Run: `rg -n "C\^0 Q4|C⁰ Q4|不可作为实现公式" "Theory/02_几何与运动学.md"`
Expected: ≥1 命中。

- [ ] **Step 6: 跑全套测试确认文档批次无回归**

Batch 0'' 是纯文档批次，理论上不影响代码。但 `test/` 中有若干前置条件测试会读取文档/参数，须确认。

Run: `julia --startup-file=no --project=. test/runtests.jl`
Expected: 全套 22/22 通过（`e117fd2` 起 `unit_czm_eigenstrain.jl` 已修复，Global Constraints"测试套件现状"）。把完整的通过清单记入下一步的 progress.md。

- [ ] **Step 7: 更新 planning 三件套与总索引，并提交**

在 `docs/planning-with-files/30_堆芯塌陷力学建模/progress.md` 追加 Batch 0'' 小节，逐条记录 8 项的处置结果（含 Task 2 Step 8 转入 Batch 7 的 1 项）、Step 6 的测试清单、以及 `tools/theory_geometry_recompute.jl` 的实测输出。

在 `docs/planning-with-files/index.md` 更新"堆芯塌陷力学建模"任务行的状态与时间。

```bash
git add "Theory/01_符号与守恒律公理.md" "Theory/02_几何与运动学.md" "Theory/03_本构理论.md" "Theory/09_附录A_符号表.md" "docs/planning-with-files/30_堆芯塌陷力学建模/progress.md" "docs/planning-with-files/index.md"
git commit -m "docs(theory): 统一层编号双约定与 A_eff 量纲，补 kappa_ss 与 C0 Q4 不兼容实现注记"
```

**Batch 0'' 完成门**：spec §7 要求"8 项阻塞矛盾全部消解；全库无 `46.6 匝/6 m/0.132 m/284 μm` 残留；几何表每行带代码字段名"。三项由 Task 2 Step 9、Task 3 Step 7、Task 4 Step 5 的 grep 共同覆盖（`284` 只允许 `Theory/06:69` 的 `ρ≈2846` 一处豁免）。未全绿不得进入 Task 5。

---

## Task 5: Option 子选项字段与默认关契约

**Files:**
- Modify: `src/Option.jl:87`（在 `czm_area_loss_threshold` 之后）
- Test: `test/test_czm_option_defaults.jl`

**Interfaces:**
- Consumes: `JuBat.Option()`（`@with_kw` 关键字构造器）
- Produces: `Option` 新字段 `czm_geo_nonlinear::Bool`、`czm_winding_prestress::Bool`、`czm_j2_plasticity::Bool`、`czm_phi_bond::Bool`、`czm_continuous_feedback::Bool`、`czm_friction_mu::Float64`，全部默认关。Batch 2/2'/3/4'/6 分别消费前 5 个。

- [ ] **Step 1: 写失败测试**

创建 `test/test_czm_option_defaults.jl`：

```julia
using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §5 兼容性契约：堆芯塌陷力学建模的所有新能力 opt-in，默认关。
# 本测试锁定"默认关"这一契约本身——任何把默认值改为开的提交都会在此失败。
@testset "堆芯塌陷子选项默认全关（spec §5）" begin
    opt = JuBat.Option()

    @test opt.czm_geo_nonlinear == false
    @test opt.czm_winding_prestress == false
    @test opt.czm_j2_plasticity == false
    @test opt.czm_phi_bond == false
    @test opt.czm_continuous_feedback == false

    # Batch 8 预留：当前无消费者，只锁定默认值
    @test opt.czm_friction_mu == 0.10
end

@testset "czm_enabled 是唯一主开关（spec §5）" begin
    # czm_enabled=true 且子选项全关时，子选项仍为 false——
    # 主开关不得隐式打开任何新能力
    opt = JuBat.Option()
    opt.czm_enabled = true

    @test opt.czm_enabled == true
    @test opt.czm_geo_nonlinear == false
    @test opt.czm_winding_prestress == false
    @test opt.czm_j2_plasticity == false
    @test opt.czm_phi_bond == false
    @test opt.czm_continuous_feedback == false
end

@testset "子选项可显式开启（关键字构造与字段赋值两条路径）" begin
    opt_kw = JuBat.Option(czm_geo_nonlinear=true)
    @test opt_kw.czm_geo_nonlinear == true
    @test opt_kw.czm_j2_plasticity == false

    opt_set = JuBat.Option()
    opt_set.czm_j2_plasticity = true
    @test opt_set.czm_j2_plasticity == true
    @test opt_set.czm_geo_nonlinear == false
end
```

- [ ] **Step 2: 运行测试，确认失败**

Run: `julia --startup-file=no --project=. test/test_czm_option_defaults.jl`
Expected: FAIL，报 `type Option has no field czm_geo_nonlinear`（`UndefVarError`/`ErrorException`），因为字段还不存在。

- [ ] **Step 3: 加字段**

`src/Option.jl`，在第 87 行 `czm_area_loss_threshold::Float64 = 0.83` 之后、`end` 之前插入：

```julia
    # 堆芯塌陷力学建模（spec 2026-08-20 §5）：全部 opt-in，默认关时行为零漂移
    czm_geo_nonlinear::Bool = false        # 完全 Green-Lagrange TL 残差 + 标准初应力 K_G（Batch 2）
    czm_winding_prestress::Bool = false    # 卷绕预应力初始应力场 σ₀(r)，缺参即 error（Batch 2'）
    czm_j2_plasticity::Bool = false        # PCC/NCC 平面应力一致 J2 返回映射，缺 sigma_y 即 error（Batch 3）
    czm_phi_bond::Bool = false             # Φ 跨匝完美粘结（phi_pairs 节点合并）（Batch 4'）
    czm_continuous_feedback::Bool = false  # 连续损伤–电–热反馈 + 界面热阻（Batch 6）
    czm_friction_mu::Float64 = 0.10        # SP Coulomb 摩擦系数（Batch 8 预留，当前无消费者）
```

- [ ] **Step 4: 运行测试，确认通过**

Run: `julia --startup-file=no --project=. test/test_czm_option_defaults.jl`
Expected: 三个 testset 全部 PASS，无 Fail/Error。

- [ ] **Step 5: 提交**

```bash
git add src/Option.jl test/test_czm_option_defaults.jl
git commit -m "feat(option): 新增堆芯塌陷力学六个默认关子选项及默认值契约测试"
```

---

## Task 6: `assemble_bulk_residual_tangent` 统一入口

Batch 1 的核心交付。只实现线弹性槽位；`geo_nl`、`plasticity`、`mech_state` 三个未实现槽位传非默认值即 `error`（spec §6，AGENTS 9.7：不静默降级）。

**Files:**
- Modify: `src/czm.jl`（在 `assemble_bulk_stiffness` 之后、`assemble_thermal_chemical_load` 之前插入新函数）
- Test: `test/test_czm_mech_core.jl`

**Interfaces:**
- Consumes: `JuBat.assemble_bulk_stiffness(czm_mesh::CohesiveMesh, param_cache::CzmParamCache) -> SparseMatrixCSC{Float64,Int64}`（`src/czm.jl:268`）；`JuBat.build_czm_cache(czm_mesh, param_cache; fix_inner) -> CZMAssemblyCache`（字段 `K_bulk`、`cohesive_geom`、`ws`）
- Produces:

```julia
assemble_bulk_residual_tangent(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache,
    mech_state=nothing;
    geo_nl::Bool=false,
    plasticity::Bool=false,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,
) -> (f_int_bulk::Vector{Float64}, K_tangent::SparseMatrixCSC{Float64, Int64})
```

  `mech_state` 为尾置可选位置参数，与 spec §4.2 的四位置参数签名一致；Batch 3 引入 `PlasticState`/`MechHistory` 时无需改签名。Task 7 的 `assemble_coupled_system` 消费本入口。

- [ ] **Step 1: 写失败测试**

创建 `test/test_czm_mech_core.jl`：

```julia
using Test
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec 2026-08-20-core-collapse-mechanics-design.md §4.2 / §7 Batch 1
# 开关全关时，新入口必须与既有 K_bulk*u 路径逐位等价；未实现的槽位必须报错。

function build_mech_core_fixture(; nθ::Int=40)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    param_cache = JuBat.compute_czm_params_per_interface(case)
    return case, param_cache
end

# 确定性的非零位移场：正弦扰动，避免全零掩盖 f_int 差异
function make_test_u(nnode::Int)
    u = zeros(Float64, 2 * nnode)
    for n in 1:nnode
        u[2*n-1] = 1e-6 * sinpi(2 * n / nnode)
        u[2*n]   = 1e-6 * cospi(2 * n / nnode)
    end
    return u
end

@testset "新入口与 K_bulk*u 逐位等价（线弹性槽位）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    K_ref = JuBat.assemble_bulk_stiffness(czm_mesh, param_cache)
    f_ref = K_ref * u

    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache)

    # 同一装配同一乘法，必须是逐位相等而非近似相等
    @test f_int == f_ref
    @test K_tan == K_ref
    @test length(f_int) == 2 * czm_mesh.nnode
    @test size(K_tan) == (2 * czm_mesh.nnode, 2 * czm_mesh.nnode)
    @test !any(isnan, f_int)
    @test !any(isnan, K_tan)
end

@testset "零位移给零内力（线弹性无预应力）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u0 = zeros(Float64, 2 * czm_mesh.nnode)

    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u0, param_cache)

    @test all(iszero, f_int)
    @test nnz(K_tan) > 0
end

@testset "缓存不变量：传入 K_bulk_cached 时直接复用同一对象" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    cache = JuBat.ensure_czm_cache(case, czm_mesh, param_cache)
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; K_bulk_cached=cache.K_bulk)

    # === 不是"数值相等"，而是同一对象：证明没有重复装配
    @test K_tan === cache.K_bulk
    @test f_int == cache.K_bulk * u
end

@testset "切线对称性（线弹性 + 平面应力各向同性）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    _, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache)

    @test norm(K_tan - transpose(K_tan), Inf) ≤ 1e-12 * norm(K_tan, Inf)
end

@testset "未实现槽位必须报错，不得静默降级（AGENTS 9.7 / spec §6）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; geo_nl=true)
    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; plasticity=true)
    # mech_state 的消费者在 Batch 3 引入；此前传入非 nothing 即报错
    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache, :dummy_state)
end

@testset "位移向量长度不符必须报错，不得截断或补零" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh

    u_short = zeros(Float64, 2 * czm_mesh.nnode - 1)
    @test_throws DimensionMismatch JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u_short, param_cache)

    u_long = zeros(Float64, 2 * czm_mesh.nnode + 3)
    @test_throws DimensionMismatch JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u_long, param_cache)
end
```

- [ ] **Step 2: 运行测试，确认失败**

Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl`
Expected: FAIL，报 `UndefVarError: assemble_bulk_residual_tangent not defined`（`JuBat` 模块中无此符号）。

- [ ] **Step 3: 实现新入口**

`src/czm.jl`，在 `assemble_bulk_stiffness` 的 `end`（第 338 行）之后、`assemble_thermal_chemical_load` 的文档字符串之前插入：

```julia
"""
    assemble_bulk_residual_tangent(czm_mesh, u, param_cache, mech_state=nothing;
                                  geo_nl=false, plasticity=false, K_bulk_cached=nothing)
        -> (f_int_bulk, K_tangent)

bulk 残差/切线的统一入口（spec 2026-08-20-core-collapse-mechanics-design.md §4.2）。

三个槽位，当前只实现第一个：

1. **线弹性**（`geo_nl=false, plasticity=false`）：`f_int = K_bulk*u`，`K_tangent = K_bulk`，
   与既有 `assemble_bulk_stiffness` 路径逐位等价。
2. **几何非线性**（`geo_nl=true`）：完全 Green-Lagrange TL + 标准初应力 `K_G`，Batch 2。
3. **J2 塑性**（`plasticity=true`）：PCC/NCC 平面应力一致返回映射，Batch 3；届时经
   `mech_state` 传入 `PlasticState`。

未实现的槽位传入非默认值一律 `error`——静默走线弹性会让上层误以为几何非线性/塑性已生效
（AGENTS 9.7）。

`mech_state` 按 spec §4.2 保留为尾置可选位置参数，使 Batch 3 引入塑性状态时无需改签名。
"""
function assemble_bulk_residual_tangent(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache,
    mech_state=nothing;
    geo_nl::Bool=false,
    plasticity::Bool=false,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing
)
    ndof = 2 * czm_mesh.nnode
    length(u) == ndof || throw(DimensionMismatch(
        "assemble_bulk_residual_tangent: u 长度为 $(length(u))，应为 $ndof " *
        "(2 × nnode=$(czm_mesh.nnode))"))

    geo_nl && error(
        "assemble_bulk_residual_tangent: geo_nl=true 尚未实现（Batch 2：完全 " *
        "Green-Lagrange TL + K_G）。不静默退回线弹性。")
    plasticity && error(
        "assemble_bulk_residual_tangent: plasticity=true 尚未实现（Batch 3：PCC/NCC " *
        "平面应力一致 J2 返回映射）。不静默退回线弹性。")
    mech_state === nothing || error(
        "assemble_bulk_residual_tangent: mech_state 的消费者在 Batch 3 引入，" *
        "当前必须传 nothing，收到 $(typeof(mech_state))。")

    K_tangent = K_bulk_cached !== nothing ? K_bulk_cached :
                assemble_bulk_stiffness(czm_mesh, param_cache)
    f_int_bulk = K_tangent * u

    return f_int_bulk, K_tangent
end
```

- [ ] **Step 4: 运行测试，确认通过**

Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl`
Expected: 6 个 testset 全部 PASS，无 Fail/Error。

若 `@test K_tan === cache.K_bulk` 失败，说明 `K_bulk_cached` 被复制而非透传；检查是否误写成 `copy(K_bulk_cached)` 或 `sparse(K_bulk_cached)`。

- [ ] **Step 5: 确认模块可干净加载**

Run（PowerShell 下用单引号包裹 Julia 表达式，避免内层双引号被吞）：`julia --startup-file=no --project=. -e 'include("src/JuBat.jl")'`
Expected: 无警告、无错误地完成加载（新符号不与既有定义冲突）。

- [ ] **Step 6: 提交**

```bash
git add src/czm.jl test/test_czm_mech_core.jl
git commit -m "feat(czm): 新增 assemble_bulk_residual_tangent 统一入口（线弹性槽位，未实现槽位报错）"
```

---

## Task 7: `assemble_coupled_system` 接线与 Batch 1 门禁

**Files:**
- Modify: `src/czm.jl:558-577`（`assemble_coupled_system` 函数体）
- Test: `test/test_czm_mech_core.jl`（追加 testset）

**Interfaces:**
- Consumes: Task 6 的 `assemble_bulk_residual_tangent(czm_mesh, u, param_cache; K_bulk_cached=...)`
- Produces: `assemble_coupled_system` 与 `assemble_coupled_system_full` 的返回值与现状逐位一致；bulk 残差/切线自此只有一个装配入口

- [ ] **Step 1: 追加等价性测试**

在 `test/test_czm_mech_core.jl` 末尾追加：

```julia
@testset "assemble_coupled_system 接线后逐位等价（Batch 1 零漂移）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    cache = JuBat.ensure_czm_cache(case, czm_mesh, param_cache)

    # 参照解：显式按"bulk 刚度 + cohesive"两路各自装配后相加
    K_bulk_ref = JuBat.assemble_bulk_stiffness(czm_mesh, param_cache)
    K_coh_ref, f_coh_ref, _, _ = JuBat.assemble_czm_system(
        czm_mesh, u, param_cache;
        damage_states=czm_mesh.damage_states,
        geom_cache=cache.cohesive_geom, ws=cache.ws)

    K_tot, f_tot, seps, tracts = JuBat.assemble_coupled_system(
        czm_mesh, u, param_cache;
        damage_states=czm_mesh.damage_states,
        K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom, ws=cache.ws)

    @test K_tot == K_bulk_ref + K_coh_ref
    @test f_tot == K_bulk_ref * u + f_coh_ref
    @test length(seps) == czm_mesh.n_cohesive
    @test length(tracts) == czm_mesh.n_cohesive
end

@testset "assemble_coupled_system_full 残差构成不变" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)
    ne = size(czm_mesh.bulk_element, 1)

    cache = JuBat.ensure_czm_cache(case, czm_mesh, param_cache)
    pe = param_cache.by_interface[:PE_PCC]
    α_eff = pe.α
    β_n = case.param.NE.Omega / 3.0
    β_p = case.param.PE.Omega / 3.0
    dT_elem = fill(1e-4, ne)
    Δsoc_n = fill(1e-4, ne)
    Δsoc_p = fill(1e-4, ne)

    K_tot, R, F_tc, _, _ = JuBat.assemble_coupled_system_full(
        czm_mesh, u, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n, Δsoc_p;
        damage_states=czm_mesh.damage_states,
        K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom, ws=cache.ws)

    _, f_int = JuBat.assemble_coupled_system(
        czm_mesh, u, param_cache;
        damage_states=czm_mesh.damage_states,
        K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom, ws=cache.ws)

    # R = F_ext + F_thermo_chem - f_int，F_ext 缺省为零；加 0.0 不改变浮点值
    @test R == F_tc - f_int
    @test !any(isnan, K_tot)
    @test !any(isnan, R)
end
```

若 `R == F_tc - f_int` 只在末位不符（差值范数约 `eps`），不要改成 `≈` 放过：这意味着 `assemble_czm_system` 的结果依赖复用的 `ws` 工作区状态（同参数两次调用不等值），本身是需要记入 findings 的发现，先定位再决定。

- [ ] **Step 2: 运行测试，确认当前已通过**

Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl`
Expected: 全部 PASS。这两个 testset 描述的是**接线前后都必须成立**的等价关系，因此现在就应通过。它们的作用是在 Step 3 改动 `assemble_coupled_system` 内部实现时守住外部行为——先确认基准为绿，改完再确认仍为绿。

若此步不通过，说明参照解写错了（例如漏传 `ws` 导致 workspace 状态不同），先修测试再继续，不要改 `src/`。

- [ ] **Step 3: 把 `assemble_coupled_system` 改走新入口**

`src/czm.jl`，把第 560–575 行

```julia
    # 固体刚度（使用缓存或重新计算）
    K_bulk = K_bulk_cached !== nothing ? K_bulk_cached : assemble_bulk_stiffness(czm_mesh, param_cache)

    # 内聚力刚度和内力（使用几何缓存和工作区，透传 param_cache）
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, param_cache; damage_states=damage_states,
        geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

    # 固体内力（线性弹性：f_int = K * u）
    f_int_bulk = K_bulk * u

    # 总刚度矩阵
    K_total = K_bulk + K_coh

    # 总内力 = 固体内力 + 内聚力内力
    f_int_total = f_int_bulk + f_int_coh
```

替换为

```julia
    # 固体残差与切线（统一入口，spec §4.2）。Batch 1 只有线弹性槽位，
    # 与原 K_bulk*u 逐位等价；Batch 2/3 的 K_G 与塑性从此处接入。
    f_int_bulk, K_bulk = assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; K_bulk_cached=K_bulk_cached)

    # 内聚力刚度和内力（使用几何缓存和工作区，透传 param_cache）
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, param_cache; damage_states=damage_states,
        geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

    # 总刚度矩阵
    K_total = K_bulk + K_coh

    # 总内力 = 固体内力 + 内聚力内力
    f_int_total = f_int_bulk + f_int_coh
```

`assemble_coupled_system_full` 无需改动：它通过 `assemble_coupled_system` 间接走新入口。

`assemble_coupled_system` 开头的 `ndof` 局部变量在改动前后都未被使用（`F_ext`/`F_thermo_chem` 两个关键字参数同样是既有的未消费参数）。这属于先前遗留，**本 Task 不清理**——零漂移批次只做接线，无关清理留给后续简化批次，以免混淆基线比对的归因。

- [ ] **Step 4: 运行定向测试**

Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl`
Expected: 8 个 testset 全部 PASS。

Run: `julia --startup-file=no --project=. test/test_assemble_coupled_system.jl`
Expected: PASS（既有测试，验证签名与形状未变）。

- [ ] **Step 5: 跑全套测试**

Run: `julia --startup-file=no --project=. test/runtests.jl`
Expected: 全套 22/22 通过——不得出现任何失败项（v1.1 起 eigenstrain 豁免作废，见 Global Constraints）。

- [ ] **Step 6: 基线快照门禁（verify_czm_standalone，方案 B）**

Run（与冻结时同环境）: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. tools/verify_czm_standalone.jl`
Expected: 网格统计行（Nodes/Bulk/Cohesive）、有效参数行（E_eff/ν_eff/α_eff/β_n/β_p）、8 行收敛对比表的每个 `OK/FAIL it D r` 字段、Summary 三行的全部数值，与 `docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md` 冻结表在打印精度下逐位一致。

任一数值不同（含 OK↔FAIL 翻转）即停止：新入口在开关全关时应当是同一算式，出现差异说明接线改变了语义（例如 `K_bulk_cached` 未透传导致重新装配，或求和次序变化）。定位后修复，不得放行；不得调参使冻结的 FAIL 条目"变好"后放行。

- [ ] **Step 7: 强制行为基线门禁**

Run: `& 'D:\Julia-1.11.2\bin\julia.exe' --startup-file=no example\testexample.jl`（`GKSwstype=100`，单线程）

Expected: exit code 0；`thermal elements = 1682`、`thermal nodes = 1763`、`result time steps = 19`、`initial voltage = 4.0367 V`、`final voltage = 3.9438 V`、`voltage drop = 0.0929 V`、`final capacity = 0.0833 Ah`、`minimum temperature = 298.15 K`、`maximum temperature = 299.00 K`、`final CZM D_max = 0.0000%`、`final CZM D_mean = 0.0000%`、`maximum normal separation = 1.2557e-14 m`、`fractured elements = 0`、`CZM converged updates = 19 / 19`，全部与 `Simplify/baseline/testexample/README.md` 冻结表一致。

Run: `Get-FileHash output\testexample\testexample_results.png -Algorithm SHA256`
Expected: `4BA6207C3CCF92DA5E37349EE335CF21A10A50B46A14CDA13DE95EEFA6CAE932`（大小写不敏感）。

任一科学指标不一致 → 停止并回退本 Task（Global Constraints）。

- [ ] **Step 8: 更新 planning 与基线记录**

在 `docs/planning-with-files/30_堆芯塌陷力学建模/progress.md` 追加 Batch 1 小节：新入口签名、`assemble_coupled_system` 接线方式、Step 5/6/7 三道门禁的实测结果、以及与 spec 的偏差记录（`PlasticState`/`MechHistory`/cache 字段延后至消费批次；`K_bulk_cached` 签名扩充，见"与 spec §4.1/§4.2 的批次归属与签名偏差"）。

在 `Simplify/baseline.md` 的批次记录表追加一行：日期、"堆芯塌陷 Batch 1：bulk 残差/切线统一入口"、行数变化、定向测试结果、"所有科学指标及 PNG SHA-256 完全一致"、PASS。

在 `docs/planning-with-files/index.md` 更新任务行状态与时间。

- [ ] **Step 9: 提交**

```bash
git add src/czm.jl test/test_czm_mech_core.jl Simplify/baseline.md "docs/planning-with-files/30_堆芯塌陷力学建模/progress.md" "docs/planning-with-files/index.md"
git commit -m "refactor(czm): assemble_coupled_system 改走 bulk 残差/切线统一入口（零漂移）"
```

---

## 自评审记录

**1. spec 覆盖检查**

| spec 条目 | 覆盖 Task |
|---|---|
| §9 Batch 0'' ①（`07` 同号 §6.4.5 + KKT 引用） | Task 2 Step 2–3 |
| §9 Batch 0'' ②（`K_uu` 双定义） | Task 2 Step 4–5 |
| §9 Batch 0'' ③（弧长法小节） | Task 2 Step 7 |
| §9 Batch 0'' ④（Φ 约束施加方式） | Task 2 Step 6 |
| §9 Batch 0'' ⑤（`κ_ss` 与 C⁰ Q4） | Task 4 Step 4 |
| §9 Batch 0'' ⑥（层编号双约定） | Task 4 Step 2 |
| §9 Batch 0'' ⑦（`A_eff` 量纲） | Task 4 Step 3 |
| §9 Batch 0'' ⑧（几何表重算 + 字段溯源，D14） | Task 3 全部 |
| §7 Batch 0'' 验收门（8 项消解 + 无残留 + 字段名） | Task 2 Step 9、Task 3 Step 7、Task 4 Step 5 |
| §7 Batch 0'' 验收门（全套测试无回归） | Task 4 Step 6 |
| §4.1 `src/Option.jl` Batch 1 | Task 5 |
| §4.2 `assemble_bulk_residual_tangent` | Task 6 |
| §4.1 `src/czm.jl` Batch 1（`assemble_coupled_system` 改走新入口） | Task 7 Step 3 |
| §7 Batch 1 新测试（新入口 ≡ `K_bulk*u`、缓存不变量） | Task 6 Step 1、Task 7 Step 1 |
| §7 Batch 1 验收门（`testexample` 基线一致） | Task 7 Step 7 |
| §7 Batch 1 验收门（基线快照不变；v1.3 起 `verify_czm_standalone.jl`） | Task 1（复核）、Task 7 Step 6（比对）；快照 `baseline_czm_standalone.md` 已于 2026-08-21 冻结 |
| §5 兼容性契约（子选项全 false 时一致） | Task 5、Task 7 Step 6–7 |
| §6 错误处理（不静默降级） | Task 6 Step 1/3 的槽位报错与维度检查 |

**已知缺口（有意为之，理由见"与 spec §4.1/§4.2 的批次归属与签名偏差"）**：`PlasticState`、`MechHistory`、`CZMAssemblyCache` 的参考构型/机械状态字段不在本计划实现，随消费批次引入。

**计划外新增的交付（v1.1 更新）**：`tools/theory_geometry_recompute.jl`（D14"防止再次漂移"的可复现实现，不改求解路径）；`docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md`（方案 B 基线冻结，2026-08-21 随 v1.1 修订执行）。原计划的 `czm_baseline_probe.jl` 修复交付因对象被 `2bf2ac7` 删除而取消，由 `verify_czm_standalone.jl` 快照替代。

**2. 占位符扫描**

无 TBD / TODO / "类似 Task N" / "适当处理" / 无代码的测试步骤。每个改动步骤都给出完整的替换前后文本或完整函数体。Task 1 Step 7 与 Task 3 Step 2 的 `<...>` 是**待实测填入的数值**，不是待决定的设计——填写规则已在同一步给出。

**3. 类型与签名一致性**

- `assemble_bulk_residual_tangent` 的签名在 Task 6 的 Interfaces、文档字符串、实现、以及 Task 7 Step 3 的调用处完全一致（4 个位置参数 + 3 个关键字参数）。
- `build_czm_cache(czm_mesh, param_cache; fix_inner)`：Task 1 Step 3 的修复与 `src/czm.jl:435` 实际签名一致。
- `bilinear_traction_state(δ_n, δ_t, ::DamageState, ::CzmInterfaceParams)`：Task 1 Step 5 传 `pe`（`CzmInterfaceParams`），与 `src/Materialmatrix.jl:68` 一致。
- `create_czm_mesh(::CzmSubmesh, ::Mesh, param)`：Task 1 Step 4 与 Task 6 fixture 均传归一化 `case.param`。
- Option 字段名在 Task 5 的测试、实现、以及 spec §5 表格中三处一致（`czm_geo_nonlinear`、`czm_winding_prestress`、`czm_j2_plasticity`、`czm_phi_bond`、`czm_continuous_feedback`、`czm_friction_mu`）。
- `get_damage_statistics` 返回字段 `max_D`/`mean_D`/`n_fractured` 与探针打印一致（`src/CzmPostProcess.jl:12`）。

---

## 修订记录

### v1.1（2026-08-21，计划评审后修订）

评审结论与用户决策（方案 B）回填，共 8 项；评审确认无误的部分（src/Theory 行号与签名、几何重算数值、基线冻结值、TDD 步骤）未改动：

1. **Task 1 重写（方案 B）**：`czm_baseline_probe.jl` 已被 `2bf2ac7` 删除，原"修复 4 处 API 漂移"交付取消；Batch 1 门禁改用 `verify_czm_standalone.jl`（spec v1.3 同步替换 §5/§7 引用）。基线快照 `baseline_czm_standalone.md` 已于本修订时冻结（HEAD `e117fd2`，实测三方法各 7/8 收敛、10.0 水平 FAIL 原样冻结）；Task 1 改为复核可复现性。工具既有瑕疵（`:66`/`:134` 传参不一致、用未合并 `thermal2D`）登记 findings，不修改不改基线。
2. **PNG 基线路径更正**：`a2caecc`（AGENTS §9.9）后 testexample 输出于 `output/testexample/testexample_results.png`；Global Constraints 与 Task 7 Step 7 两处同步更正（原路径已过期，且存在哈希迁移前旧文件造成假通过的风险）。
3. **测试套件预期更正**：`e117fd2` 已修复 `unit_czm_eigenstrain.jl`（60/60），全套 22/22；删除"既有失败登记"豁免，Task 4 Step 6、Task 7 Step 5 门禁改为全绿。
4. **Task 2 KKT 清理扩充**：评审发现 `07:19`、`07:80`、`07:1162`（本章小结）三处现行框架 KKT 表述不在原 Step 3 范围，而 Step 9 自设门禁（"仅历史沿革"）会失败；三处补入 Step 3，Step 9 预期同步改写。
5. **Task 3 计数订正**：实测 `46.6|0.132` 22 处（01:7、02:15、04:0）+ `284` 18 处（06:69 豁免）+ `\approx 6` 5 处；原"25/11/2、共 38 处"与实际不符（04 的"2 处"与"284 2 处"为同一组，疑笔误）。零命中门禁不变，计数仅供参考。
6. **`03:77` 出处联动**：该行"亦等于 spec §1.3 中的 t_repeat≈284 μm"出处随 `:44`/`:331` 同标准改为代码字段溯源。
7. **偏差登记补第三项**：`assemble_bulk_residual_tangent` 比 spec §4.2 冻结签名多 `K_bulk_cached` 关键字参数（零漂移快路径所需，内部接口），章节标题相应改为"与 spec §4.1/§4.2 的批次归属与签名偏差"。
8. **评审方法说明**：上述各项均对当前 HEAD（`e117fd2`）逐项核实——含 `src/czm.jl` 引用代码逐字比对、Theory 行号逐行核对、几何数值独立复算、基线冻结表比对；非推测性修订。
