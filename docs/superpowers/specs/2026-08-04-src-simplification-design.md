# JuBat src/ 简化设计（Reduce 模式）

**日期**：2026-08-04
**作者**：Reduce 模式 brainstorming 产物
**状态**：待 executing-plans 消化
**关联文档**：
- `Simplify/reducing-ai-code-bloat-plan.md`（原始 bloat 治理计划，本 spec 对其做本地化与具体化）
- `md/00_代码风格规范.md`、`md/00_文档索引与命名规范.md`
- `CLAUDE.md`、`AGENTS.md`

**口径说明（统一）**：除非特别注明，"src/" 行数与文件数均**仅指 `src/*.jl`**（36 个文件 / 11,184 行）；`src/parameters/*.jl`（5 个文件 / 711 行）单独统计。这是为避免基线口径混乱。

---

## 0. 摘要与决策

### 0.1 目标

依据 `Simplify/reducing-ai-code-bloat-plan.md` 的方法论，对 `src/*.jl` 全部 36 个文件（共 11,184 行）+ `src/parameters/*.jl` 5 个文件（共 711 行）执行一次 **Reduce 模式**：只做简化、删除、合并，不引入新功能。

### 0.2 五项已对齐决策

| # | 决策 | 选择 |
|---|---|---|
| 1 | Reduce 范围 | **只产出计划，不动代码**。执行交给后续 writing-plans → executing-plans |
| 2 | 计划粒度 | **逐文件展开**（每文件一节，列删除/合并/保留/风险/前置测试） |
| 3 | 优先级准则 | **重复 + dead code 驱动** |
| 4 | 文件边界 | **选择性允许新建**（如 `VariableKeys.jl`，但必须在 spec 中给出理由） |
| 5 | 总体骨架 | **Strategy B：依赖层级自底向上** —— 参数/变量 → 几何/网格 → 物理模型 → 求解器 → 后处理 |

### 0.3 不变量

- **活跃公共 API 不删减**：`JuBat.jl` 的 `export` 表中，**经 grep 确认仍被外部调用者（example/、test/、用户脚本）使用的**不删减；**dead code 的 export 允许同步移除**（如 SPM/P2D/ThermalLumped 经 grep 确认无调用者时，对应 export 一并移除）。所有 export 变更必须在 PR 描述中列出，并在 CHANGELOG 标 breaking change
- 数值核心（High-risk-leave-alone 桶，见 §4）只记录不动手
- 所有 `example/*.jl`（minimal / SPMe_Thermal / czm_cycle / testexample / jellyroll_stress_displacement）必须跑通且结果在 1% 内一致

---

## 1. 对原始 bloat 计划的评审

### 1.1 计划的优点（保留）

- **Build mode / Reduce mode 分离**：JuBat 当前 `czm-refactor` 分支上混改文件、未提交测试与文档并存，正是混模式的证据
- **Phase 0 基线快照**：可量化（行数 / 重复率 / 复杂度 / dead code）
- **Phase 1 Map-before-cut**：本 spec §6-§9 即此产物
- **Phase 3 小批量可验证 + characterization test 先行**
- **§4 提示词模式**：实用，本 spec §5 吸收

### 1.2 计划的不足（本项目错配）

| 原计划条款 | 错配 | 本 spec 的修正 |
|---|---|---|
| §3 Phase 0 推荐工具链（eslint/ruff/vulture/knip/jscpd） | 语言完全错配（JuBat 是 Julia） | §2 给出 Julia 侧替代与脚本策略 |
| §3 Phase 2 三桶分类（Delete/Consolidate/Leave） | 缺风险维度 | §4 新增第四桶 High-risk-leave-alone |
| §3 Phase 3 假设测试已存在 | JuBat 大文件测试覆盖不足 | 每文件标前置测试；§5 characterization test 要求 |
| §5 Pre-merge checklist 通用 | 未结合 JuBat 归一化/字段层级陷阱 | §11 给出 JuBat 专属版 |
| §2 未定义退出标准与切换机制 | Reduce 何时收尾不明 | §3 明确退出阈值与 baseline 跟踪 |
| §6 未列 Julia 生态技能 | — | §2 + §5 补 |
| 全文无"耦合接口边界"概念 | 应沿耦合数据流切，而非按文件大小切 | Strategy B 的层级即耦合接口边界 |

---

## 2. 工具链本地化（替代原 §3 Phase 0）

### 2.1 Julia 生态替代品

| 原计划推荐 | Julia 侧替代 | 备注 |
|---|---|---|
| eslint / ruff（lint） | `StaticLint.jl`、`Jet.jl` | 仅类型与基本 lint，规则远少于 eslint |
| ts-prune / knip（dead code） | **无等价物** | 手写 grep 脚本 |
| jscpd（duplication） | **无等价物** | 手写 grep + AST 比对 |
| cloc（size trends） | `wc -l` + `git log --stat` | 够用 |
| sonarqube | **无** | — |
| vulture / radon（Python） | **无** | — |

### 2.2 实际策略

由于 Julia 生态缺工具，本次 Reduce 采用：

1. **grep 脚本为主**，放在 `Simplify/scripts/`：
   - `dead_code_grep.sh`：批量 grep 每个顶层函数名的调用者
   - `dup_finder.sh`：基于正则模式（`^function\s+\w+`、`variables\[\"...\"\]`、`try.*catch`）找重复簇
   - `baseline.sh`：跑 `wc -l src/*.jl` 并写入 `Simplify/baseline.md`
2. **行数 / 重复簇数 / dead code 数** 三项作为可量化指标，写入 `Simplify/baseline.md`，每完成一个文件更新一行
3. **版本控制策略**：`Simplify/scripts/` 与 `Simplify/baseline.md` 纳入 git；`baseline.md` 每文件更新一行即一次 commit（保留历史轨迹，便于回看）

### 2.3 基线快照

执行 executing-plans 之前必须先固化基线。`Simplify/baseline.md`（待 executing 阶段生成）应包含：

```
## Baseline 2026-08-04
- src/*.jl：36 文件 / 11,184 行
- src/parameters/*.jl：5 文件 / 711 行（单独统计，主线 Reduce 仅涉及验证是否 dead）
- variables["..."] 使用数：477 跨 17 文件
- 硬编码键表数：3 处（Variables.jl / CallModel.jl / CouplingState.jl）
- 重复簇：D1-D10（见 spec §6）
- 兜底 try/catch：~13 处（见 spec §7 类别 A）
- 兼容入口：~9 处（见 spec §7 类别 C）
```

**统计口径**：
- 行数以 `wc -l` 为准，**仅 `src/*.jl`**（不含 `src/parameters/`）
- "硬编码键表数"指 `StandardVariables`、`create_element_workspace`、`copy_element_results` 三处**显式键名字面量列表**，不含散落在代码中的 `variables["..."]` 访问点
- "兜底 try/catch" 仅指 §7.1 类别 A 列出的；不含类别 C 中 `try/finally` 资源清理

---

## 3. 退出标准

### 3.1 量化指标

- **行数**：`src/*.jl` 总行数 −10% 以上（基线 11,184 → 目标 ≤ 10,066）
- **重复簇清零**：D1、D3、D4 必须消除；D2 收敛到 `VariableKeys.jl` 单一来源
- **Dead code**：§9 嫌疑表经 grep 验证后，确认无调用者的全部删除
- **兜底清零**：§7 类别 A（错误吞没 try/catch）全部删除或改 rethrow + 注释；类别 B（静默降级）全部改 fail-fast
- **兼容入口验证**：§7 类别 C 每条 grep 验证；无外部调用者的删除
- **止损条件**：若 executing-plans 阶段累计产出超过 12 个 PR 仍未达标，暂停并回到 brainstorming 重新评估策略

### 3.2 不变量验证

- 所有 `example/*.jl` 跑通，结果与基线快照在 1% 内一致
- `test/unit_czm_*` 全套通过
- 新增 characterization test 覆盖 czm.jl / CzmSolve.jl / CouplingState.jl / ThermalDistributed.jl 各至少一个 happy path
- 活跃公共 API（§0.3 定义）的 export 不删减

### 3.3 节奏（来自原 §7）

| 频率 | 活动 |
|---|---|
| 每完成一个文件 | 更新 `Simplify/baseline.md` 对应行 |
| 每个 PR | §11 Pre-merge checklist |
| 退出阈值达成 | 本次 Reduce 收尾，归档 spec |

---

## 4. 四桶 triage（替代原 §3 Phase 2）

每文件必须显式标桶。

| 桶 | 处理 | JuBat 适用 |
|---|---|---|
| **Delete** | 删除 | dead code（SPM / P2D / ThermalLumped / ring 验证后；类别 A try/catch；类别 C 兼容入口无调用者） |
| **Consolidate** | 合并 | D1-D9 重复簇；类别 B 改 fail-fast；不一致约定归一 |
| **Leave alone** | 不动 | 重复密度低、工作正常、改动收益小 |
| **High-risk-leave-alone** | 不动，仅记录 | 数值核心求解体 |

### 4.1 High-risk-leave-alone 清单

以下文件/函数即便发现重复也**只记录不动手**：

- `czm.jl`：`assemble_czm_system` (261-454)
- `CzmSolve.jl`：`solve_czm_basic_step` (168-480)
- `SPMe.jl`：`SPMe_element` (37) 内部数值
- `ThermalDistributed.jl`：`ThermalDistributed2D` (1) 求解体
- `Parallelsolution.jl`：`solve_branch_currents` (358) 主迭代
- `Materialmatrix.jl`：CZM 本构函数（`bilinear_traction` 等）
- `Mechanical.jl`：`Calstressdisp` (112) 颗粒扩散应力

理由：数值代码对边界条件、量纲、归一化极其敏感，"consolidate" 的潜在收益远小于回归风险。如必须动，需要先写 characterization test 锁定行为。

---

## 5. Reduce 模式执行守则（吸收原 §4 + §6）

### 5.1 提示词模式（原 §4）

执行 executing-plans 时，每个文件的 Reduce 任务应在 prompt 中包含：

- **Root-cause fix**："不要新增分支/标志来处理此情况；找出为何现有逻辑不覆盖，直接修。"
- **No net-new abstractions**："新建函数/类/文件前，先检查是否能复用现有的；若新建，一行说明为何旧代码不能复用。"
- **Deletion-first**："函数经过多次打补丁已膨胀。仅按当前需求重写，不保留任何无法指出活跃调用者的旧分支/注释/标志。"
- **Bound the diff**："改动仅限本任务直接相关文件。无关问题另行列出，不内联修复。"
- **Force a plan before code**："先给短计划：改什么、删什么、为什么。等确认再写 diff。"
- **Explicit reduce-mode**："你在清理模式，不是功能模式。不新增功能。唯一目标：减重复、删死码、简化控制流，所有现有测试保持通过。"
- **Self-review**："以严格的简洁性代码评审者身份复审你的 diff，指出任何不严格必要的新增。"

### 5.2 必备技能（原 §6）

| 技能 | 在 JuBat 的具体体现 |
|---|---|
| Reading diffs critically | 检查每次 diff 的 scope，而非仅正确性 |
| Static analysis tooling | §2.2 的 grep 脚本 |
| Characterization testing | 改动 czm.jl / CzmSolve.jl 等前先写行为快照测试 |
| Architectural judgment | 区分"值得合并的重复" vs "巧合重复" |
| Session/prompt discipline | Build 与 Reduce 严格分离，本 spec 即 Reduce 会话 |
| Reading dependency maps | Strategy B 的层级依赖此判断 |

---

## 6. 跨文件重复主表（核心索引）

> 所有结论带 file:line 证据。Strategy B 各层引用此表。

| # | 重复簇 | 涉及文件:行 | 性质 | 估计冗余 |
|---|---|---|---|---|
| **D1** | 循环求解双实现：`solve_phase_with_export` / `solve_cycle_with_export` 几乎复制 `solve_phase` / `solve_cycling` | CycleData.jl:32, 268 vs CycleSolver.jl:118, 217 | 大量拷贝 | ~400 行 |
| **D2** | 三处硬编码变量键表：`StandardVariables` + `create_element_workspace` + `copy_element_results` 各维护一份字符串键 | Variables.jl:1（`StandardVariables`）、:157（`create_element_workspace`）；CallModel.jl:247（`copy_element_results`）；CouplingState.jl 散布 | 同步风险 | 三处**显式键名列表**共 ~30 个键名重复定义（非所有字符串字面量） |
| **D3** | ThermalDistributed 双胞胎对：`apply_convection_bc`/`!` 和 `apply_cool_method`/`!` | ThermalDistributed.jl:49 vs 178；99 vs 214 | 双胞胎函数 | ~90 行 |
| **D4** | 热源双实现：`compute_heat_sources` 与 `compute_heat_sources_with_czm` 路径分裂 | ThermalDistributed.jl:385 vs 518 | 部分重叠 | ~70 行 |
| **D5** | 后处理/导出五重奏：PostProcessing / CycleSolver `_postprocess_*` / CsvExport / CycleData `export_*_csv` / CzmPostProcess | PostProcessing.jl:1, 159, 219, 236, 254；CsvExport.jl 全文；CycleData.jl:425；CzmPostProcess.jl:99 | 责任分散 | 多处 |
| **D6** | `variables["..."]` 字符串键 477 处 | 17 文件（机械 73 / 后处理 75 / P2D 59 / 变量 68 / 热 41 / SPMe 30 …） | 单点变更成本高 | 不冗余但脆弱 |
| **D7** | `E_coat` 模量入口分散：两个独立 `@assert param.PE.E_coat > 0` | CouplingState.jl:307, 309；Mechanical.jl:165 入口 | 检查重复 | ~10 行 |
| **D8** | CZM 损伤更新双调用路径（已部分修复） | CycleSolver.jl（docs 显示已删） | 历史问题 | 已清 |
| **D9** | 多重 method overload：`MultiSPMeLayout` 双构造、`update_czm_damage!` 双方法（CouplingState.jl）；`assemble_coupled_system` + `_full`（czm.jl） | CouplingState.jl:98/109；497/613；**czm.jl:737/777** | 需 triage 是 overload 还是 fork | triage 标准已定（见 §10.4 czm.jl 与 §10.5 CouplingState.jl） |
| **D10** | 单函数小文件 | ring.jl（1 fn）、ThermalPolar2D.jl（1 fn）、CzmUnitMesh.jl（1 fn）、install.jl（5 行） | 拆分过度 | 整合潜力 |

---

## 7. 兜底处理与打补丁清单（用户指正后补扫）

### 7.1 类别 A：错误吞没 try/catch（必须删除或改 rethrow）

| 文件:行 | 模式 | 处理 |
|---|---|---|
| `CsvExport.jl:98-194` | 7 个 `_write_*` 各被 `try/catch e @warn ... end` 包裹 | 默认删 try/catch，让 I/O 错误抛出；如个别确需保留，注释写明恢复保证 |
| `Mechanical.jl:294-297` | `U_M = try ... catch e @warn "using zero displacement" end` | 删除；力学失败应抛错而非静默 0 位移 |
| `CzmSolve.jl:220, 319, 372, 381` | 4 处线性求解 `try \" catch fallback` | **逐处判定**（见 §7.1.1） |
| `Solve.jl:273-310, 418-424` | 步进失败恢复 + 末尾清理 | diff 评估：步进失败路径是否掩盖根因 |

**统一原则**：默认让错误抛出；仅在能证明"该错误可恢复且恢复策略明确"时保留，且必须注释写明恢复保证。

#### 7.1.1 CzmSolve.jl 四处 try/catch 逐处判定

| 行 | 当前 fallback | spec 决策 | 理由 |
|---|---|---|---|
| 220 | `Δu = try K\F catch zeros end` | **保留 + 加注释** | 奇异/病态矩阵时返回零位移增量，避免一步发散毁掉整个 transient；恢复策略明确（下一步重新尝试），但必须在注释中写明"返回零增量、依赖下一步恢复" |
| 319 | `tangent = try ... catch zeros end` | **保留 + 加注释** | 切线刚度计算失败时用零矩阵，配合 `backtrack_line_search!` 收敛机制；恢复路径明确 |
| 372 | `delta_u_R = try ... catch zeros end` | **保留 + 加注释** | 弧长法 R 分支求解失败时返回零，与 381 配对 |
| 381 | `delta_u_F = try ... catch zeros end` | **保留 + 加注释** | 弧长法 F 分支求解失败时返回零，与 372 配对 |

四处均属"数值核心中已设计好的恢复路径"，**保留但每处加注释**："本 try/catch 是 [奇异矩阵/弧长法分支] 的恢复路径，返回零让下一步重试；移除会导致 transient 一步发散即崩溃"。

### 7.2 类别 B：静默降级 @warn + 继续（改 fail-fast）

| 文件:行 | 当前行为 | 处理 |
|---|---|---|
| `CouplingState.jl:540` | NaN 时重置 `czm_mesh` 为 zeros | 改 `@assert` 或 `error()` 显式失败 |
| `Solve.jl:123` | "外部状态长度不匹配，回退到模型初始化" | 改 `@assert length(y0)==expected` |
| `PostProcessing.jl:276, 284` | SOH/断裂超阈值 @warn 终止 | **保留**（合理运行时终止） |
| `Variables.jl:244` | 时间步超预分配动态扩容 | **保留**（合理运行时行为） |
| `CouplingState.jl:340, 532, 593` | @warn 刚度过软 / NaN / Solve 问题 | **保留**（合理警告） |

### 7.3 类别 C：向后兼容入口（grep 验证后处理）

| 文件:行 | 兼容代码 | 处理 |
|---|---|---|
| `CouplingState.jl:611` | `update_czm_damage!` 6 参数兼容入口 | grep 验证；无外部调用者则删 |
| `SetCase.jl:114` | `Case` 5 参数兼容构造器 | 同上 |
| `czm.jl:192` | `interface_nodes = [[]]` 旧字段 | grep 是否仍被读 |
| `czm.jl:232` | δ_czm 缺失回退为 1 | 评估必要性 |
| `czm.jl:821` | 矩阵空/全零回退旧值 | 数值兜底，单独评估 |
| `Materialmatrix.jl:327` | "旧方案 δ_czm = L 时行为不变" | 与 czm.jl:232 同根 |
| `SetParams.jl:339` | "cohesive 数据缺失时回退 scale.L" | 同根 |
| `CouplingState.jl:37, 452, 553` | "向后兼容""默认 1.0" | 逐一评估 |
| `Option.jl:55, 64` | `solveType` 的 forward/backward 分支 | grep 实际取值；未走的分支可删 |

### 7.4 类别 D：占位 TODO（不在本次 Reduce 范围）

- `parameters/Jellyroll.jl:160-189`：7 处 `TODO 用户提供实测值`
- 这是数据缺口，应在 `md/01_参数定义与归一化.md` 跟踪，**不属于代码简化**

### 7.5 类别 E：补丁累积痕迹

每处问"现在还需要吗"：

- `CycleSolver.jl:168` 注释 `# 旧版在此处调用 update_czm_damage!...` —— 旧路径已删，注释可清
- `Solve.jl:17, 132` `elseif case.opt.solveType == "backward"` 出现两次 —— 验证 forward/backward 是否真被走
- `czm.jl:821` 注释 `# 主导条件数。矩阵为空/全零时回退旧值`

---

## 8. 不一致约定清单（盲点 1 补扫）

### 8.1 命名前缀 5 种方言

| 前缀风格 | 代表 | 文件 |
|---|---|---|
| `Set*` PascalCase | `SetMesh`, `SetCase`, `SetParams` | SetMesh/SetCase/SetParams.jl |
| `create_*` snake_case | `create_czm_mesh`, `create_element_workspace`, `create_unit_czm_strip` | czm/Variables/CzmUnitMesh.jl |
| `setup_*` snake_case | `setup_thermal2D_mesh` | Jellyrollmodel.jl |
| `build_*` snake_case | `build_czm_cache`, `build_czm_submesh`, `build_thermal_to_czm_interp`, `build_arc_length_augmented_matrix` | czm/Jellyrollmodel/CzmSolve/CouplingState.jl |
| `Model*` / `Initialize*` | `ModelInitialisation`, `ModelInitialisation_MultiSPMe`, `initialize_currents` | Initialisation/Parallelsolution.jl |

`md/00_代码风格规范.md §1.2` 规定"功能前缀+功能描述"，但未规定何时用 `Set` vs `create` vs `setup` vs `build`。**Reduce 阶段不统一改名**（公共 API 变更风险高），仅在 spec 记录；如未来要统一，应有独立 PR。

### 8.2 其他不一致

- **英式 vs 美式拼写**：`ModelInitialisation`（英式）vs `initialize_currents`（美式）
- **错误处理 95 处跨 19 文件**：`@assert/error/@warn/throw` 混用，密度差异大（CouplingState 14 处 vs CycleSolver 1 处）。**Reduce 阶段不改风格**，但 §7 类别 A/B 的具体改动会顺带减少部分 `@warn`
- **`case.param.scale` vs `case.param_dim.scale`**：125 处跨 12 文件。**Reduce 阶段不强行统一**，但在改每个文件时抽查：数值计算应用 `case.param.scale`（已归一化），结果还原应用 `case.param_dim.scale`（物理单位）

---

## 9. Dead code 嫌疑表

| 嫌疑 | 文件:行 | 验证方法 | 处理 |
|---|---|---|---|
| P2D 全套（5 函数 287 行） | P2D.jl 全文 | `grep -r "P2D(" --include="*.jl"` 排除 `P2D_` 前缀误匹配 | 无内部调用者则整文件删 + 移除 JuBat.jl export |
| SPM 全套（85 行） | SPM.jl 全文 | `grep -r "\bSPM(" --include="*.jl"` 排除 SPMe | 同上 |
| `ThermalLumped` | Thermal.jl:1 | `grep -r "ThermalLumped"` | 无调用者则删 |
| `ring_mesh` + `ThermalPolar2D_Ring` | ring.jl:7；ThermalPolar2D.jl:1 | `grep -r "ring_mesh\|ThermalPolar2D_Ring"` | 无调用者则删两个文件 |
| `RecordMatrix!`、`ErrorEstimation` | Solve.jl:442, 452 | `grep -r "RecordMatrix!\|ErrorEstimation"` | 无调用者则删 |
| `clone_damage_states`、`clone_czm_mesh_with_damage` | CzmSolve.jl:17, 31 | grep | 无调用者则删 |
| `Assemble1D` | Assemble.jl:28 | grep | 同上 |
| `install.jl` | 全文 5 行 | `grep "include(\"install"` in JuBat.jl | 未 include 则删 |
| `czm_mesh.interface_nodes` 旧字段 | czm.jl:192 | `grep "interface_nodes"` | 无读取者则删 |

---

## 10. 逐文件计划（Strategy B 五层）

### 10.1 第 1 层：参数与变量

#### `SetParams.jl`（522 行）
- **桶**：Leave alone
- **删除**：无（struct 定义文件，删字段破坏公共 API）
- **合并**：`NormaliseParam` (350) 已是统一入口，保留
- **保留**：`ChooseCell` (249)、所有 `mutable struct`
- **风险**：低；归一化逻辑分散在单函数内可读性差，但动它风险高
- **前置测试**：现有 `test/` 覆盖归一化（`md/11_电化学验证方案.md`）
- **重复密度**：低

#### `Variables.jl`（277 行）—— **D2 核心**
- **桶**：Consolidate
- **删除**：`StandardVariables` (1) 与 `create_element_workspace` (157) 两份独立字符串键字面量表
- **合并**：统一到 `VariableKeys.jl` 的 `const ELEMENT_KEYS / GLOBAL_KEYS` 常量；`StandardVariables` 内部调用 `create_element_workspace` 拼装全局字典
- **保留**：`Variable_update!` (236)
- **风险**：中；变更影响 17 文件 477 处 `variables["..."]` 访问
- **前置测试**：随机跑一次 `minimal_example.jl`，把所有 variables 键 dump 出来作快照
- **新建文件**：`src/VariableKeys.jl`，集中所有键常量
  - **理由**：键名是跨 17 文件的契约，单一来源是根因修复（与"选择性允许新建"决策一致）

#### `Option.jl`（88 行）
- **桶**：Leave alone
- **删除**：无（验证 `solveType` forward/backward 是否被走，未走可清分支但非必须）
- **合并**：`CycleOption` + `PhaseType` 理论上仅 CycleSolver 用，但已 export，保持原位
- **保留**：原样
- **风险**：极低

#### `parameters/{Jellyroll, Enertech, LGM50, Northrop, Ring}.jl`
- **桶**：Delete 候选（除 Jellyroll）
- **删除**：grep 验证 `ChooseCell("Enertech"|"LGM50"|"Northrop"|"Ring")` 调用者；无则删整个文件
- **保留**：`Jellyroll.jl`（主线必保）
- **风险**：低
- **前置测试**：保留的每个参数集都要 `ChooseCell(...)` 跑通

#### `CallModel.jl`（268 行）—— **D2 成员**
- **桶**：Consolidate
- **删除**：`copy_element_results` (247) 中 12 个键字面量
- **合并**：改用 `VariableKeys.jl` 的 `ELEMENT_RESULT_KEYS` 元组循环构造字典
- **保留**：`CallModel` (211)、`CallModel_MultiSPMe` (22)、`map_czm_damage_to_thermal` (8)
- **风险**：中；`CallModel_MultiSPMe` 是核心调度
- **前置测试**：`test/unit_czm_*` 部分覆盖

#### `SetCase.jl`（116 行）
- **桶**：Leave alone（验证 5 参数兼容构造器）
- **删除**：`Case` 5 参数兼容构造器 (114) grep 验证后处理
- **保留**：`SetCase` (1)、`Case` struct (100)
- **风险**：极低

---

### 10.2 第 2 层：几何与网格

#### `SetMesh.jl`（767 行）
- **桶**：High-risk-leave-alone
- **理由**：基础设施，动一处影响所有模型
- **行动**：不动，仅记录已审查

#### `Jellyrollmodel.jl`（690 行）
- **桶**：High-risk-leave-alone
- **理由**：几何决定下游一切
- **行动**：不动

#### `CzmUnitMesh.jl`（108 行）—— **D10**
- **桶**：Leave alone
- **理由**：仅一个函数 `create_unit_czm_strip` (7)，但用于 unit test 隔离
- **行动**：保留独立，spec 标注"D10 评估后保留——unit test 隔离有合理理由"

#### `ring.jl`（70 行）—— **D10 + Dead code 嫌疑**
- **桶**：Delete 候选
- **行动**：先 grep `ring_mesh` 调用者；仅 `example/` 用则保留；无调用者则删

---

### 10.3 第 3 层：物理模型

#### `SPMe.jl`（290 行）
- **桶**：High-risk-leave-alone（数值核心）
- **删除/合并**：仅对比 `SPMe_variables!` (101) vs `SPMe_variables` (209) 是否双胞胎；如是，合并为单一函数 + `mutate` 标志
- **保留**：`SPMe` (1)、`SPMe_element` (37)、`SPMe_BC` (184)
- **风险**：高
- **前置测试**：`md/11_电化学验证方案.md`

#### `SPM.jl`（85 行）—— **Dead code 嫌疑**
- **桶**：Delete 候选
- **行动**：grep `SPM(` 排除 SPMe 误匹配；无内部调用者整文件删 + 移除 export

#### `P2D.jl`（287 行）—— **Dead code 嫌疑**
- **桶**：Delete 候选
- **行动**：grep `P2D(` 排除 `P2D_` 前缀；无内部调用者整文件删 + 移除 export

#### `Electrode{Diffusion, Potential}.jl`（18+19 行）和 `Electrolyte{Diffusion, Potential}.jl`（30+23 行）
- **桶**：Consolidate 候选（文件粒度）
- **倾向决策**：**保留独立**。理由：(1) 合并收益仅 ~20 行；(2) 每个文件是独立 FEM 单元算子，与 `Assemble.jl`（仅 40 行的拼装工具）语义不同；(3) 跨文件 include 已在 `JuBat.jl` 固化，迁移成本 > 收益
- **风险**：低

#### `czm.jl`（873 行）
- **桶**：High-risk-leave-alone（数值核心）
- **D9 triage**（与 CouplingState.jl 同节同标准）：
  - `assemble_coupled_system` (737) + `assemble_coupled_system_full` (777)：grep 两函数名；**仅其一被调用 → 删另一**；**两版都被调用 → 对比内部逻辑差异**：≤5 行差异合并为 `full::Bool=false` 关键字；>5 行差异保留双实现并注释差异点
  - `build_czm_cache` (621) vs `ensure_czm_cache` (719)：grep + diff 验证是否 wrapper 关系；如 `ensure_czm_cache` 仅是 `build_czm_cache` + memoize，合并为 `ensure_czm_cache`（保留 memoize 行为）
- **保留**：`CohesiveElement` (1)、`DamageState` (26)、`create_czm_mesh` (58)、`assemble_czm_system` (261-454)
- **风险**：高（数值核心，仅在不破坏 `unit_czm_newton.jl` 快照的前提下做 triage）
- **前置测试**：`test/unit_czm_bilinear.jl` / `unit_czm_newton.jl`

#### `Mechanical.jl`（360 行）—— **D7 + 类别 A**
- **删除（类别 A）**：`Mechanical.jl:294-297` 的 try/catch + 零位移降级
- **合并（D7）**：`thermal_diffusion_stress_2D` (165) 入口 `@assert param.PE.E_coat > 0` 与 CouplingState.jl:307 重复 —— 抽 `assert_E_coat(case)` helper，放 SetParams.jl 末尾
- **保留**：`Mechanicaloutput` (1)（对比 PostProcessing 是否能合并由 executing 决定）、`Calstressdisp` (112)
- **风险**：中

#### `Materialmatrix.jl`（427 行）
- **桶**：High-risk-leave-alone（CZM 本构敏感）
- **行动**：不动

#### `ThermalDistributed.jl`（557 行）—— **D3 + D4 核心**
- **桶**：Consolidate
- **删除**：验证后删 `ThermalDistributed2D_Ring` (325) 与 `ThermalRing2D_BC` (376)（若 ring.jl 也删）
- **合并 D3**（**统一为 `!` 后缀版本**——JuBat 风格规范 §7.2 规定修改参数的函数用 `!` 后缀；非 `!` 版本仅在不修改入参时使用）：
  - **方案**：删 `apply_convection_bc` (49)，保留 `apply_convection_bc!` (178) 作为唯一实现；非 `!` 版本的调用点改为显式 `K2 = copy(K); apply_convection_bc!(K2, F2, ...)`（grep 验证非 `!` 版本调用点数量；如 ≤3 处，迁移成本可接受）
  - 同理 `apply_cool_method` (99) → `apply_cool_method!` (214)
  - **签名**（合并后）：`apply_convection_bc!(K, F, mesh, is_outer, case; edge_cache=nothing)`、`apply_cool_method!(K, F, mesh, case)`
- **合并 D4**：
  - **方案**：删 `compute_heat_sources_with_czm` (518)，保留 `compute_heat_sources` (385) 作为唯一实现，新增关键字 `czm_data::Union{Nothing, NamedTuple}=nothing`
  - **签名**（合并后）：`compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme=false, czm_data=nothing)`
  - `czm_data` 为 `NamedTuple`（字段 `mesh`、`damage` 等）或 `nothing`；为 `nothing` 时等价于原 `compute_heat_sources` 行为
  - **风险点**：需 diff 验证两函数除 CZM 分支外的主体逻辑是否完全一致；如不一致则保留双实现并标注差异
- **保留**：`ThermalDistributed2D` (1)、`ThermalDistributed2D_BC` (286)
- **风险**：中-高（数值核心 + 结构改动）
- **前置测试**：`md/12_热模型验证方案.md` 圆环精确解 + `example/热模块验证/thermal_verify.jl`

#### `Thermal.jl`（80 行）—— **Dead code 嫌疑**
- **桶**：Delete 候选
- **行动**：grep `ThermalLumped`；无调用者删

#### `ThermalPolar2D.jl`（120 行）—— **D10**
- **桶**：依赖 ring.jl 决策
- **行动**：若 ring.jl 删则本文件也删；否则保留

---

### 10.4 第 4 层：求解器

#### `Solve.jl`（471 行）
- **桶**：Consolidate（边缘）+ 类别 A/B
- **删除**：
  - 类别 A：`Solve.jl:273-310, 418-424` 的 try/catch diff 评估
  - 类别 B：`Solve.jl:123` "回退到模型初始化" 改 `@assert`
  - `RecordMatrix!` (442)、`ErrorEstimation` (452) grep 验证后处理
- **保留**：`Solve` (1) 函数签名不变
- **风险**：高（主求解器）
- **前置测试**：所有 example

#### `CzmSolve.jl`（675 行）
- **桶**：High-risk-leave-alone（数值核心）+ 类别 A
- **删除**：`clone_damage_states` (17)、`clone_czm_mesh_with_damage` (31) grep 验证
- **类别 A 处理**：4 处 `try \" catch fallback` (220, 319, 372, 381) **逐处评估**，能改 rethrow 则改；必要的保留并注释
- **合并**：`apply_czm_dirichlet!` (130) + `zero_czm_bc_entries!` (137) diff 后决定
- **保留**：`CZMResult` (1)、`solve_czm_basic_step` (168-480)、`backtrack_line_search!` (107)、`fill_czm_result!` (144)
- **风险**：高
- **前置测试**：`test/unit_czm_newton.jl`

#### `Parallelsolution.jl`（453 行）
- **桶**：High-risk-leave-alone（数值核心）
- **保留**：`solve_branch_currents` (358) 公共 API 不动
- **风险**：高（已知 cutoff 误触发问题，见 CLAUDE.md §8.3）
- **前置测试**：需新增 characterization test（端电压曲线快照）

#### `CycleSolver.jl`（546 行）—— **D1 保留方**
- **桶**：Leave alone（D1 删除在 CycleData.jl 一侧）
- **保留**：`solve_phase` (118)、`solve_cycling` (217)、`compute_cs0_from_soc` (495)、`apply_initial_soc!` (533)
- **理由**：CLAUDE.md 明确 `solve_cycling` 已审过决定保持内联
- **风险**：中
- **前置测试**：`example/czm_cycle_example.jl`、`example/coupled_czm_thermal_example.jl`

#### `Initialisation.jl`（152 行）
- **桶**：Consolidate（验证双胞胎）
- **删除**：`extract_element_state` (120)、`get_thermal_dofs` (130)、`update_state` (139) grep 验证
- **合并**：`ModelInitialisation` (1) vs `ModelInitialisation_MultiSPMe` (54) —— diff 确认是否双胞胎
- **风险**：中

#### `Assemble.jl`（40 行）
- **桶**：Leave alone
- **删除**：`Assemble1D` (28) grep 验证
- **保留**：原样

---

### 10.5 第 5 层：后处理与导出

#### `CycleData.jl`（623 行）—— **D1 删除方**
- **桶**：Consolidate（最大删除动作）
- **删除**：
  - `solve_phase_with_export` (32-268, ~236 行) —— 整体删除，复制了 CycleSolver.solve_phase
  - `solve_cycle_with_export` (268-425, ~157 行) —— 整体删除
- **合并**：将"导出钩子"作为可选 callback 注入 `CycleSolver.solve_phase` / `solve_cycling`。**callback 设计**：
  - **签名**：`export_callback(step_data::TimeStepData, cycle::Int, phase::PhaseType) -> Nothing`
  - **关键字参数**：`solve_phase(...; export_callback::Union{Nothing,Function}=nothing)`；`solve_cycling(...; export_callback=nothing)`
  - **调用点**：在 `solve_phase` 内每个时间步成功后（拿到 `step_data` 之后）调用；callback 为 `nothing` 时跳过
  - **异常处理**：callback 抛异常**不捕获**——直接终止求解（与 Reduce 模式 fail-fast 原则一致）；用户应在 callback 内部 `try/catch` 自己处理 I/O 错误
  - **行为等价性**：原 `solve_phase_with_export(case, ...; export_interval=k)` 等价于 `solve_phase(case, ...; export_callback=(sd,c,p) -> iszero(c*step_count % k) && push!(records, sd))`
- **保留**：`TimeStepData` (7)、`CycleExportData` (23)、`export_cycle_data_to_csv` (425)、`load_cycle_data_from_csv` (514)
- **风险**：中（~400 行删除）
- **前置测试**：必须先跑 `solve_cycle_with_export` 的 happy path 作快照，确保 callback 化后输出一致

#### `CsvExport.jl`（636 行）—— **D5 + 类别 A**
- **桶**：Consolidate
- **删除（类别 A）**：7 个 `_write_*` 的 try/catch 包裹
- **合并**：8 个 `_write_*` 函数模板化 —— 抽 `write_table(path, rows, cols)` helper，每个 `_write_*` 仅声明列定义
- **保留**：`export_cycling_csv` (88) 公共 API、`CsvExportOptions` (20)
- **风险**：中
- **前置测试**：`example/czm_cycle_example.jl` 跑通 + CSV 输出文件 diff

#### `PostProcessing.jl`（350 行）—— **D5**
- **桶**：Consolidate（边缘）
- **合并**：`PostProcessing` (1) 与 `Mechanical.Mechanicaloutput` 职责重叠 —— **倾向合并**：`Mechanicaloutput` 内联进 `PostProcessing` 作为末尾的力学分支，删除 `Mechanicaloutput` 名字；保持 `PostProcessing` 作为唯一入口（grep 验证 `Mechanicaloutput` 调用点，迁移至 `PostProcessing`）
- **保留**：7 个 `_postprocess_*` / `_print_*` helper 短小，不合并
- **风险**：中

#### `CouplingState.jl`（762 行）—— **D2 + D9 + 类别 B**
- **桶**：Consolidate
- **删除（类别 B）**：`CouplingState.jl:540` NaN 重置改 fail-fast
- **合并 D9 triage**（**判定标准**：grep 两版的调用点；**两版都有外部调用者 = overload（保留）**；**只有一版有调用者 = fork（删未用版）**）：
  - `MultiSPMeLayout` 双构造 (98, 109)：grep `MultiSPMeLayout(` 所有调用点；如两版签名都被调用，保留双构造并加注释说明差异；如仅一版被调用，删另一版
  - `update_czm_damage!` 双方法 (497, 613)：同上 grep；spec 已知 613 行注释"自动构建 CzmLayout 并委托给 3 参数版本"——如 6 参数版无外部调用者，删 6 参数版
  - **注意**：`assemble_coupled_system` (737) + `_full` (777) 在 **czm.jl** 而非本文件，其 triage 见 §10.4 czm.jl 段
- **合并 D2 部分**：CouplingState.jl 内散布的字符串键字面量（与 D2 主表对齐，约 50 处 `variables["..."]` 访问与键字面量）改用 `VariableKeys.jl` 的常量；**与 Variables.jl / CallModel.jl 的 ~30 个键名列表保持单一来源**
- **保留**：所有 struct
- **风险**：高（CZM 装配核心）
- **前置测试**：`test/unit_czm_*`

#### `CzmPostProcess.jl`（117 行）—— **D5**
- **桶**：Consolidate（验证）
- **合并**：5 个函数与 CouplingState 的 CZM 函数对比，重叠则合并入 CouplingState
- **风险**：低

#### `Tools.jl`（191 行）
- **桶**：Leave alone
- **删除**：`q4_center_gradients` (178) grep 验证
- **保留**：8 个独立 helper 不合并
- **风险**：低

#### `install.jl`（5 行）
- **桶**：Delete 候选
- **行动**：grep `include("install` in JuBat.jl；未 include 则删

#### `JuBat.jl`（89 行）
- **桶**：Leave alone（同步表）
- **删除**：dead code 删除后同步移除 export 表对应名字
- **合并**：include 顺序如调整，同步更新
- **风险**：低

---

## 11. JuBat 专属 Pre-merge Checklist

### 11.1 通用项（原 §5 + 兜底/补丁）

- [ ] 删除/合并是否引入新的字符串键字面量？（如是，加进 `VariableKeys.jl`）
- [ ] 数值核心（§4.1）是否被改动？如否，跳过下两项
- [ ] 数值改动是否跑了对应的 `md/11-13` 验证方案脚本？
- [ ] 新增 try/catch 是否必要？如必要，注释中是否写明恢复保证？
- [ ] 新增 @warn 是否伴随降级？如伴随，是否已改为 fail-fast？
- [ ] 新增"兼容入口""旧字段"是否 grep 过无调用者？
- [ ] `variables["..."]` 新访问点是否走 `VariableKeys.jl` 常量？

### 11.2 JuBat 归一化项

- [ ] `scale.L^2` / `scale.q` / `scale.E_coat` 归一化是否正确使用？
- [ ] `PE.E`（颗粒）vs `PE.E_coat`（极片）是否混用？
- [ ] 数值计算用 `case.param.scale`，结果还原用 `case.param_dim.scale` 是否一致？

### 11.3 兼容性项

- [ ] example（minimal/SPMe_Thermal/czm_cycle/testexample/jellyroll_stress_displacement）跑通且结果在 1% 内？
- [ ] `JuBat.jl` 的 `export` 表未删减？

---

## 12. 风险登记与回滚

### 12.1 风险等级

| 风险 | 文件 | 缓解 |
|---|---|---|
| 数值核心改动引入回归 | czm.jl / CzmSolve.jl / ThermalDistributed.jl / SPMe.jl | 标 High-risk-leave-alone；如必须动，先写 characterization test |
| 公共 API 误删 | SPM.jl / P2D.jl / dead code 嫌疑 | grep 验证 + changelog 标注 breaking change |
| 兼容入口删除破坏旧脚本 | CouplingState.jl:611 等 | grep 验证；保留的迁移到 deprecation 路径 |
| 字符串键集中化引入 typo | Variables.jl / CallModel.jl / CouplingState.jl | characterization test：dump 所有键快照对比 |

### 12.2 回滚策略

- 每个文件一个独立 commit，便于 `git revert`
- `Simplify/baseline.md` 记录每步行数与测试状态
- 任何一步测试不通过立即停止，回滚到上一个绿色 commit
- 数值核心（High-risk-leave-alone）改动必须独立 PR，不与其它清理混合

---

## 附录 A：执行顺序建议（供 executing-plans 参考）

按 Strategy B 优先级 + 风险递增：

1. **低风险预热**：install.jl / SetCase.jl 兼容入口 / parameters/* dead 验证
2. **第 1 层完成**：VariableKeys.jl 新建 → Variables.jl / CallModel.jl / CouplingState.jl 键集中化
3. **第 5 层主菜**：CycleData.jl D1 删除 + CycleSolver callback 化
4. **第 3 层重复消除**：ThermalDistributed.jl D3/D4 合并
5. **兜底清理**：类别 A/B 跨文件处理
6. **Dead code 批量**：SPM.jl / P2D.jl / ThermalLumped / ring 验证后删除
7. **D9 triage**：CouplingState overload 评估
8. **High-risk 仅在测试到位时**：czm.jl / CzmSolve.jl 内部 diff 确认项

---

## 附录 B：基线数据（2026-08-04 快照）

```
src/*.jl：36 文件 / 11,184 行
src/parameters/*.jl：5 文件 / 711 行（单独统计）
最长文件：czm.jl (873)
variables["..."] 使用：477 跨 17 文件
try/catch：~13 处
@warn：~25 处
向后兼容入口/注释：~15 处
```

执行 executing-plans 时每完成一项更新本表。

---

**本 spec 完成。下一步：spec-document-reviewer 评审循环 → 用户审核 → writing-plans 生成实现计划。**
