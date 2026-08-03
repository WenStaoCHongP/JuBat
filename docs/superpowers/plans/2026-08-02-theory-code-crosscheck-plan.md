# 理论代码对照检查文档 实施计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 `md/对照/` 下生成 8 卷覆盖全物理场的理论 ↔ 代码对照文档，作为 `czm-refactor` 分支重构前的现状基线快照。

**Architecture:** 公式级粒度、表格驱动（每条对照记录 6 列）；分卷按物理场切分；每卷严格按 md 小节顺序，先列"公式清单 + 函数清单"草稿再合并对照表，杜绝边读边写。

**Tech Stack:** 纯 Markdown 文档；Julia 源码只读；用 Grep `^function\s+\w+` 定位函数边界。

**Spec:** `docs/superpowers/specs/2026-08-02-theory-code-crosscheck-design.md`

---

## 文件结构

新建目录 `md/对照/`，包含 8 个文件：

| 文件 | 责任 |
|---|---|
| `md/对照/00_总览与索引.md` | 入口：图例、命名规则、跨卷导航、统计汇总、验证现状索引、跨卷公式编号总索引 |
| `md/对照/01_参数与归一化对照.md` | md 01 + md 15 ↔ `SetParams.jl`, `parameters/*.jl`, `Materialmatrix.jl` |
| `md/对照/02_几何与网格对照.md` | md 02 ↔ `SetMesh.jl`, `Jellyrollmodel.jl`, `ring.jl`, `CzmUnitMesh.jl` |
| `md/对照/03_边界条件对照.md` | md 03 ↔ `ThermalDistributed.jl` 中 BC 段、`Solve.jl` 电化学 BC |
| `md/对照/04_电化学SPMe对照.md` | md 04 ↔ `SPMe.jl`, `Electrode*.jl`, `Electrolyte*.jl`, `SPM.jl`, `P2D.jl` |
| `md/对照/05_热模型对照.md` | md 05 + md 07 ↔ `Thermal.jl`, `ThermalDistributed.jl`, `ThermalPolar2D.jl` |
| `md/对照/06_CZM对照.md` | md 06 + md 14 ↔ `czm.jl`, `CzmSolve.jl`, `Mechanical.jl`, `CzmPostProcess.jl`, `CouplingState.jl` |
| `md/对照/07_算法与求解对照.md` | md 08 + 09 + 10 ↔ `Solve.jl`, `CallModel.jl`, `Parallelsolution.jl`, `CycleSolver.jl`, `SetCase.jl`, `Initialisation.jl`, `Variables.jl`, `Option.jl`, `Assemble.jl` |

**不修改任何 src/md 源文件**。所有产出仅限 `md/对照/` 下。

---

## 通用约定（每个卷任务都遵守）

**4 步工作流**（每个卷任务都按此推进，不要跳号）：

1. **读 md** — 打开本卷对应的 md 文档，全文通读；在新建的 `md/对照/NN_xxx.md` 文件末尾 `## 草稿` 小节列出每小节的公式清单（公式编号 + 一句话内容）
2. **读 src** — 用 `Grep "^function\\s+\\w+" path/to/file.jl` 拿到本卷所有 src 文件的函数边界（function 行 + end 行），同样写入草稿区
3. **逐小节对照** — 按 md 小节顺序，每个小节生成一张 markdown 表（6 列模板见 spec §3）；完成后卷末追加 `## 编号总表`（汇总本卷所有自编号 `(N.xM)` → md 出处映射，无自编号则写"无"）
4. **自检 + 提交** — 行号精度、一致性等级、模板合规三查后提交；草稿区**保留不删**作为附录

**禁止行为**：
- ❌ 边读边写（必须先列草稿再合并）
- ❌ 重写或修订 md 公式（属理论审计，本任务外）
- ❌ 修改任何 src 代码
- ❌ 修复发现的真实偏差（仅记录，另开计划）

**行号查证命令**（每个有疑问的函数都用此核对）：
```bash
# 查看 src 中某文件所有 function 边界
```
实际工具调用：`Grep pattern="^function\\s+\\w+" path="<file>" output_mode="content" -n=true`

**自检清单**（提交前每条都要确认）：
- [ ] 每个 md 小节至少有一条对照记录（或显式"本节无公式"）
- [ ] 所有代码位置行号落在 `function ... end` 之间，不越界
- [ ] 一致性列只用了 ✅/🟡/⚠️/❌/— 五种标记
- [ ] ✅ 条目都有依据；🟡/⚠️/❌ 都指向 docs/ 笔记或 issue；— 都写了卡点
- [ ] 代码独有实现用 `--` + `【缺失】` 混在主表
- [ ] 实现摘要 ≤ 25 字
- [ ] 卷末有"草稿"小节保留中间产物
- [ ] 卷末有"编号总表"（若有自编号）

---

## Chunk 1: 任务 1 — 00 总览与索引骨架

### Task 1: 建立 `md/对照/00_总览与索引.md` 骨架

**Files:**
- Create: `md/对照/00_总览与索引.md`

- [ ] **Step 1: 创建目录与文件骨架**

写入以下内容（统计字段先占位 `TBD`，待各卷完成后回填）：

```markdown
# 理论代码对照检查文档 — 总览与索引

**目的**：JuBat 项目 `czm-refactor` 分支重构前的现状基线快照。覆盖全物理场（电化学 SPMe / 二维热 / CZM / 耦合 / 几何边界 / 参数归一化），公式级粒度。

**生成日期**：2026-08-02
**对应 spec**：`docs/superpowers/specs/2026-08-02-theory-code-crosscheck-design.md`

---

## 命名规则

本子目录使用 `NN_主题_对照.md` 命名，独立于 `md/00_文档索引与命名规范.md` 的 01–15 理论编号体系，**不更新**该索引文档。

---

## 一致性图例

| 标记 | 含义 |
|---|---|
| ✅ 一致 | 公式与代码数学等价 |
| 🟡 单位/归一化差异 | 数学形式一致，归一化常数/尺度因子/单位换算偏差 |
| ⚠️ 形式差异 | 同物理但实现形式不同（离散/迭代/近似） |
| ❌ 不一致 | 已确认偏差、漏项、错符号、错边界 |
| — 待核 | 暂时无法判定 |

---

## 分卷导航

| 卷 | 文件 | 覆盖 md | 覆盖 src |
|---|---|---|---|
| 01 | [参数与归一化](./01_参数与归一化对照.md) | md 01 + 15 | SetParams.jl, parameters/*.jl, Materialmatrix.jl |
| 02 | [几何与网格](./02_几何与网格对照.md) | md 02 | SetMesh.jl, Jellyrollmodel.jl, ring.jl, CzmUnitMesh.jl |
| 03 | [边界条件](./03_边界条件对照.md) | md 03 | ThermalDistributed.jl BC, Solve.jl 电化学 BC |
| 04 | [电化学 SPMe](./04_电化学SPMe对照.md) | md 04 | SPMe.jl, Electrode*.jl, Electrolyte*.jl, SPM.jl, P2D.jl |
| 05 | [热模型](./05_热模型对照.md) | md 05 + 07 | Thermal.jl, ThermalDistributed.jl, ThermalPolar2D.jl |
| 06 | [CZM](./06_CZM对照.md) | md 06 + 14 | czm.jl, CzmSolve.jl, Mechanical.jl, CzmPostProcess.jl, CouplingState.jl |
| 07 | [算法与求解](./07_算法与求解对照.md) | md 08 + 09 + 10 | Solve.jl, CallModel.jl, Parallelsolution.jl, CycleSolver.jl, SetCase.jl, Initialisation.jl, Variables.jl, Option.jl, Assemble.jl |

**不进对照的文件**：install.jl, JuBat.jl, Tools.jl, CycleData.jl, CsvExport.jl, PostProcessing.jl（纯基础设施/输出，无公式对应）

---

## 跨物理场一致性总览

> 待各卷完成后回填。

| 卷 | 公式总数 | ✅ | 🟡 | ⚠️ | ❌ | — | 代码独有 (`--`) |
|---|---|---|---|---|---|---|---|
| 01 参数与归一化 | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 02 几何与网格 | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 03 边界条件 | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 04 电化学 SPMe | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 05 热模型 | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 06 CZM | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 07 算法与求解 | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| **合计** | **TBD** | **TBD** | **TBD** | **TBD** | **TBD** | **TBD** | **TBD** |

---

## 验证现状索引

> md 11–13 验证方案不进对照卷，在此列出脚本路径与执行状态。

### md 11 电化学验证方案
- TBD（待回填：脚本路径、执行状态、对应卷条目反向链接）

### md 12 热模型验证方案
- `example/热模块验证/thermal_verify.jl` — TBD
- `example/热模块验证/thermal_error_source_analysis.jl` — TBD
- `example/热模块验证/thermal_equivalent_lumped_compare.jl` — TBD

### md 13 耦合验证方案
- TBD

---

## 跨卷公式编号总索引

> 待各卷完成后回填。汇总所有自编号（`(N.xM)`）的 md 出处。

TBD

---

## 相关 docs/ 笔记反向索引

> 列出本对照文档反向引用过的 `docs/` 下笔记，便于追溯。

TBD
```

- [ ] **Step 2: 验证目录创建成功**

Run: `ls md/对照/`
Expected: 列出 `00_总览与索引.md`

**注**：此时 01–07 卷文件均未创建，骨架中"分卷导航"的相对链接（`./01_xxx.md` 等）为**预期死链**，待 Task 9 回填后验证。

- [ ] **Step 3: 提交**

```bash
git add "md/对照/00_总览与索引.md"
git commit -m "docs(check): 新增理论代码对照卷 00 总览与索引骨架"
```

---

## Chunk 2: 任务 2–4 — 卷 01/02/03

### Task 2: 卷 01 — 参数与归一化对照

**Files:**
- Create: `md/对照/01_参数与归一化对照.md`
- Read (md): `md/01_参数定义与归一化.md`, `md/15_颗粒与极片模量区分.md`
- Read (src): `src/SetParams.jl`, `src/parameters/*.jl`, `src/Materialmatrix.jl`

**覆盖范围说明**：
- md 01 §1 物理参数定义 → `SetParams.jl` 的 Electrode/Cell/Cohesive/Scale/Params struct
- md 01 §2 电化学归一化 → `SetParams.jl` 的 `NormaliseParam`
- md 01 §3 热学归一化 → 同上
- md 01 §4 电流归一化 → 同上
- md 15 颗粒 vs 极片模量 → `SetParams.jl`（E vs E_coat）、`Materialmatrix.jl`（实际函数清单：`thermal_capacity_weights_2d`, `thermal_anisotropic_conductivity_2d`, `bilinear_traction`, `bilinear_tangent`, `update_damage`, `compute_gap_conductance`, `effective_area_factor` 等）。**注意**：CLAUDE.md 与 spec 提到的 `compute_effective_coating_modulus(case)` 在 src 中**不存在**（已被内联或未实现）——发现此现象时，对照条目一致性填 `— 待核`，备注"CLAUDE.md §9.4 提及但 src 未找到实现，疑似待补"。

- [ ] **Step 1: 读 md 全文，列公式清单到草稿**

打开 `md/01_参数定义与归一化.md` 与 `md/15_颗粒与极片模量区分.md`，在文件末尾 `## 草稿` 小节先写：
```
### md 公式清单（待合并）
- §X.X eq.(1.1) ...
- §X.X eq.(1.2) ...
- ...
```

- [ ] **Step 2: Grep src 函数边界，列函数清单到草稿**

对 `src/SetParams.jl`、`src/Materialmatrix.jl` 执行：
`Grep pattern="^function\\s+\\w+|^struct\\s+\\w+" path="<file>" output_mode="content" -n=true`

补充 `src/parameters/*.jl` 中相关 struct。在草稿区写：
```
### src 函数清单（待合并）
- SetParams.jl:LINE-LINE (NormaliseParam)
- SetParams.jl:LINE-LINE (struct Electrode)
- Materialmatrix.jl:LINE-LINE (compute_effective_coating_modulus)
- ...
```

- [ ] **Step 3: 逐小节合并对照表**

按 md 01 小节顺序，每小节生成一张 6 列表。表格模板见 spec §3 示例。md 15 同理在卷末追加。完成后按"通用约定 4 步工作流 Step 3"在卷末追加 `## 编号总表`。

- [ ] **Step 4: 自检**

按"通用约定 → 自检清单"逐条核对。

- [ ] **Step 5: 提交**

```bash
git add "md/对照/01_参数与归一化对照.md"
git commit -m "docs(check): 新增理论代码对照卷 01 参数与归一化"
```

---

### Task 3: 卷 02 — 几何与网格对照

**Files:**
- Create: `md/对照/02_几何与网格对照.md`
- Read (md): `md/02_几何与网格.md`
- Read (src): `src/SetMesh.jl`, `src/Jellyrollmodel.jl`, `src/ring.jl`, `src/CzmUnitMesh.jl`

**覆盖范围**：
- 阿基米德螺旋线 `r(θ) = a + bθ` → `Jellyrollmodel.jl`（`jellyroll_collector_seed_mesh` 等）
- collector-seeded 网格 → `Jellyrollmodel.jl`
- COH2D4 单元 → `CzmUnitMesh.jl`
- 圆环网格 → `ring.jl`
- CohesiveMesh 结构 → `SetMesh.jl` 或 `CzmUnitMesh.jl`

- [ ] **Step 1: 读 md 全文，列公式清单到草稿**

- [ ] **Step 2: Grep src 函数边界，列函数清单到草稿**

对所有 4 个 src 文件执行 `Grep pattern="^function\\s+\\w+|^struct\\s+\\w+" -n=true`。

- [ ] **Step 3: 逐小节合并对照表**

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/02_几何与网格对照.md"
git commit -m "docs(check): 新增理论代码对照卷 02 几何与网格"
```

---

### Task 4: 卷 03 — 边界条件对照

**Files:**
- Create: `md/对照/03_边界条件对照.md`
- Read (md): `md/03_边界条件.md`
- Read (src): `src/ThermalDistributed.jl`（仅 BC 相关函数）、`src/Solve.jl`（电化学 BC 段）

**注意**：边界条件在 src 中**跨文件散落**。先用 Grep 定位所有 BC 相关函数，再回 md 比对。建议搜索关键词：`_BC`, `apply_convection`, `apply_cool_method`, `apply_bc`。

- [ ] **Step 1: 读 md 全文，列公式清单到草稿**

- [ ] **Step 2: Grep BC 函数边界**

```text
Grep pattern="^function\\s+\\w+(BC|bc|convection|cool)" path="src/" output_mode="content" -n=true
Grep pattern="apply_convection_bc|apply_cool_method|_BC" path="src/ThermalDistributed.jl" output_mode="content" -n=true
```
补充 `src/Solve.jl` 中的电化学 BC 应用点。

- [ ] **Step 3: 逐小节合并对照表**

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/03_边界条件对照.md"
git commit -m "docs(check): 新增理论代码对照卷 03 边界条件"
```

---

## Chunk 3: 任务 5–7 — 卷 04/05/06

### Task 5: 卷 04 — 电化学 SPMe 对照

**Files:**
- Create: `md/对照/04_电化学SPMe对照.md`
- Read (md): `md/04_电化学模型_SPMe.md`
- Read (src): `src/SPMe.jl`, `src/ElectrodeDiffusion.jl`, `src/ElectrodePotential.jl`, `src/ElectrolyteDiffusion.jl`, `src/ElectrolytePotential.jl`, `src/SPM.jl`, `src/P2D.jl`

**覆盖范围**：
- 颗粒扩散方程 → `ElectrodeDiffusion.jl`、`SPMe.jl:37-93 (SPMe_element)` 中扩散段
- 电解液守恒 → `ElectrolyteDiffusion.jl`, `ElectrolytePotential.jl`
- Butler-Volmer → `SPMe.jl` 内 BV 段
- 端电压 → `SPMe.jl:SPMe_variables!` 等
- 机械耦合（仅电化学侧输入） → 散落

- [ ] **Step 1: 读 md 全文，列公式清单到草稿**

- [ ] **Step 2: Grep 函数边界（7 个 src 文件）**

```text
Grep pattern="^function\\s+\\w+" path="src/SPMe.jl" output_mode="content" -n=true
Grep pattern="^function\\s+\\w+" path="src/ElectrodeDiffusion.jl" output_mode="content" -n=true
...（其余 5 文件同理）
```

- [ ] **Step 3: 逐小节合并对照表**

注意：本卷公式数预估 30–40，**可能超过 50**。若超过，按 spec §3 "长表拆分"规则：每 md 小节一张子表，前加 `### 小节名` 三级标题。

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/04_电化学SPMe对照.md"
git commit -m "docs(check): 新增理论代码对照卷 04 电化学SPMe"
```

---

### Task 6: 卷 05 — 热模型对照

**Files:**
- Create: `md/对照/05_热模型对照.md`
- Read (md): `md/05_热模型_二维分布式.md`, `md/07_界面热阻模型.md`
- Read (src): `src/Thermal.jl`, `src/ThermalDistributed.jl`, `src/ThermalPolar2D.jl`

**覆盖范围**：
- 能量方程 → `ThermalDistributed.jl`、`ThermalPolar2D.jl`
- 分层热源（Q_rxn / Q_rev / Q_ohm） → `ThermalDistributed.jl:385-516 (compute_heat_sources)` 等
- 各向异性导热 → `Thermal.jl`
- 极坐标 FVM → `ThermalPolar2D.jl`
- 界面热阻/损伤耦合 → `ThermalDistributed.jl` 中相关段（参考 `docs/thermal_verify/findings.md`）

- [ ] **Step 1: 读 md 全文（2 个文件），列公式清单**

- [ ] **Step 2: Grep 函数边界（3 个 src 文件）**

- [ ] **Step 3: 逐小节合并对照表**

md 07 在卷末追加。

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/05_热模型对照.md"
git commit -m "docs(check): 新增理论代码对照卷 05 热模型"
```

---

### Task 7: 卷 06 — CZM 对照

**Files:**
- Create: `md/对照/06_CZM对照.md`
- Read (md): `md/06_内聚力模型_CZM.md`, `md/14_粘性正则化.md`
- Read (src): `src/czm.jl`, `src/CzmSolve.jl`, `src/Mechanical.jl`, `src/CzmPostProcess.jl`, `src/CouplingState.jl`

**覆盖范围**：
- 双线性牵引-分离律 → `czm.jl`
- CZMResult 结构 → `czm.jl`
- 热-化学载荷 → `czm.jl:assemble_thermal_chemical_load`、`CzmSolve.jl` 中 `assemble_coupled_system`（配套装配路径）
- Newton-Raphson 求解 → `CzmSolve.jl:168-263 (solve_czm_basic_step)`
- 粘性正则化 → `CzmSolve.jl`（visc_beta 参数）、`czm.jl`
- 应力计算 → `Mechanical.jl`
- 后处理 → `CzmPostProcess.jl`
- 耦合状态 → `CouplingState.jl`

- [ ] **Step 1: 读 md 全文（2 个文件），列公式清单**

- [ ] **Step 2: Grep 函数边界（5 个 src 文件）**

注意：`czm.jl` 与 `CzmSolve.jl` 都很长（663 / 686 行），函数多。务必先用 Grep 拿到全部边界再对照。

- [ ] **Step 3: 逐小节合并对照表**

md 14 粘性正则化在卷末追加。**关注 spec §3 示例已标注的 ⚠️ 形式差异**——md 14 在 md 06 公式 (6.7) 中的体现情况。

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/06_CZM对照.md"
git commit -m "docs(check): 新增理论代码对照卷 06 CZM"
```

---

## Chunk 4: 任务 8–9 — 卷 07 + 回填

### Task 8: 卷 07 — 算法与求解对照

**Files:**
- Create: `md/对照/07_算法与求解对照.md`
- Read (md): `md/08_逐单元算法.md`, `md/09_分流求解器.md`, `md/10_参数传递与模块架构.md`
- Read (src): `src/Solve.jl`, `src/CallModel.jl`, `src/Parallelsolution.jl`, `src/CycleSolver.jl`, `src/SetCase.jl`, `src/Initialisation.jl`, `src/Variables.jl`, `src/Option.jl`, `src/Assemble.jl`

**覆盖范围**：
- 多 SPMe 并行架构、状态向量设计 → `Solve.jl`, `Variables.jl`
- 分层热源算法 → `Solve.jl`、关联到卷 05
- Newton-Raphson 分流 → `Parallelsolution.jl:358 (solve_branch_currents)`（**注意**：CLAUDE.md §7.5 提到的 `solve_branch_currents_newton` 在 src 中不存在，实际函数名为 `solve_branch_currents`；配套函数 `newton_iteration:200`、`line_search:297`、`detect_cutoff_elements:161`）
- 截止电压检测 → `Solve.jl`
- CZM 失效处理 → `Solve.jl`、`CycleSolver.jl`
- Case/variables 结构 → `SetCase.jl`, `Variables.jl`
- CycleSolver → `CycleSolver.jl`
- Option/CycleOption → `Option.jl`
- 矩阵组装 → `Assemble.jl`

- [ ] **Step 1: 读 md 全文（3 个文件），列公式清单**

- [ ] **Step 2: Grep 函数边界（9 个 src 文件）**

- [ ] **Step 3: 逐小节合并对照表**

md 08/09/10 按 spec §3 小节级拆分；若公式数 < 50 可各为一张大表。md 10 偏架构图，可能公式少，主要是数据流对应——可放宽为"功能点 ↔ 代码模块"映射，但仍用 6 列模板。

- [ ] **Step 4: 自检**

- [ ] **Step 5: 提交**

```bash
git add "md/对照/07_算法与求解对照.md"
git commit -m "docs(check): 新增理论代码对照卷 07 算法与求解"
```

---

### Task 9: 回填 00 总览与索引

**Files:**
- Modify: `md/对照/00_总览与索引.md`

- [ ] **Step 1: 统计各卷一致性计数**

打开卷 01–07，逐卷统计：公式总数、✅ / 🟡 / ⚠️ / ❌ / — / `--` 各几条。把数字记到草稿。

- [ ] **Step 2: 回填"跨物理场一致性总览"表**

把上一步的 TBD 全部替换为实际数字，"合计"行求和。

- [ ] **Step 3: 回填"验证现状索引"**

打开 `example/热模块验证/`（用 `Glob pattern="example/**/*.jl"`）确认现有验证脚本路径与状态。md 11/12/13 各列脚本路径 + 一句话状态（如"已修正 2026-03-15"、"未实现"、"通过"）。把对应 md 的章节反向链接到对照卷相关条目。

- [ ] **Step 4: 回填"跨卷公式编号总索引"**

扫描卷 01–07 末尾的"编号总表"小节，把所有自编号 `(N.xM)` 汇总到此表，格式：
```
| 自编号 | md 出处 | 主题 |
|---|---|---|
| (4.x1) | md/04 §2.1（位于 eq.(4.12) 与 eq.(4.13) 之间） | 颗粒扩散系数温度修正 |
```

- [ ] **Step 5: 回填"相关 docs/ 笔记反向索引"**

扫描所有卷的"备注"列，列出所有被引用过的 `docs/` 下文件，去重排序。

- [ ] **Step 6: 自检**

- [ ] **Step 7: 提交**

```bash
git add "md/对照/00_总览与索引.md"
git commit -m "docs(check): 回填理论代码对照卷 00 总览与索引统计"
```

- [ ] **Step 8: 最终验收**

确认 `md/对照/` 下有 8 个文件；`ls md/对照/` 应列出：
```
00_总览与索引.md
01_参数与归一化对照.md
02_几何与网格对照.md
03_边界条件对照.md
04_电化学SPMe对照.md
05_热模型对照.md
06_CZM对照.md
07_算法与求解对照.md
```

`git log --oneline md/对照/` 应有 9 条提交（00 骨架 + 7 卷 + 00 回填）。

---

## 完成判据

- [ ] `md/对照/` 下 8 个文件全部存在
- [ ] 每个对照卷的每个 md 小节至少有一条对照记录（或显式"本节无公式"）
- [ ] 所有代码位置行号经 Grep 核对落在 `function ... end` 之间
- [ ] 一致性列只用 ✅/🟡/⚠️/❌/— 五种标记
- [ ] 00 总览的统计、验证索引、跨卷编号索引、docs 反向索引全部回填
- [ ] 9 条提交按顺序落地

## 风险与边界

- **发现 md 公式错误**：不在本任务修复，仅在备注列指向"建议在 md 审计任务中处理"
- **发现代码 bug**：不在本任务修复，仅记录为 ❌ 不一致并指向 docs/ 已有笔记或新建 issue
- **md 公式无编号**：用自编号 `(N.xM)`，并在卷首"编号说明"声明
- **某 md 小节确实无公式**：显式标注"本节无公式（仅文字说明/架构图）"
- **src 中找不到对应实现**：填 `--` + `【缺失】`（代码独有）或反之记为"md 独有"，一致性填 `—` 待核
