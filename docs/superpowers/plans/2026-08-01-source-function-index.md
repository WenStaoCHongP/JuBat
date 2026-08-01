# 源码函数索引文档 实施计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 `src/` 下 35 个根目录源文件（不含 `install.jl`）生成函数级注释文档（含 [DEBUG]/[PLACEHOLDER]/[COMPLEX-CHECK] 三类标注），输出到 `md/源码函数索引/`。

**Architecture:** 按物理模块 9 批次生成；每个源文件 → 一份 md；最后用 `_索引.md` 汇总 TOC 和全局三类标注统计。所有引用使用「文件名:行号」格式；不修改任何源码。

**Tech Stack:** 纯文档（Markdown）；依赖 spec: `docs/superpowers/specs/2026-08-01-source-function-index-design.md`

**基准 commit:** `dffa4c5`（spec 评审通过时）

---

## 文件结构

新建 `md/源码函数索引/` 目录，含 36 份文件：

| 输出文件 | 对应源文件 | 批次 |
|----------|------------|------|
| `_索引.md` | （汇总） | 9 |
| `JuBat.md` | `src/JuBat.jl` | 1 |
| `Option.md` | `src/Option.jl` | 1 |
| `SetParams.md` | `src/SetParams.jl` | 1 |
| `SetCase.md` | `src/SetCase.jl` | 1 |
| `SPMe.md` | `src/SPMe.jl` | 2 |
| `SPM.md` | `src/SPM.jl` | 2 |
| `P2D.md` | `src/P2D.jl` | 2 |
| `ElectrodeDiffusion.md` | `src/ElectrodeDiffusion.jl` | 2 |
| `ElectrodePotential.md` | `src/ElectrodePotential.jl` | 2 |
| `ElectrolyteDiffusion.md` | `src/ElectrolyteDiffusion.jl` | 2 |
| `ElectrolytePotential.md` | `src/ElectrolytePotential.jl` | 2 |
| `Thermal.md` | `src/Thermal.jl` | 3 |
| `ThermalDistributed.md` | `src/ThermalDistributed.jl` | 3 |
| `ThermalPolar2D.md` | `src/ThermalPolar2D.jl` | 3 |
| `czm.md` | `src/czm.jl` | 4 |
| `CzmSolve.md` | `src/CzmSolve.jl` | 4 |
| `CzmPostProcess.md` | `src/CzmPostProcess.jl` | 4 |
| `CzmUnitMesh.md` | `src/CzmUnitMesh.jl` | 4 |
| `Mechanical.md` | `src/Mechanical.jl` | 4 |
| `Materialmatrix.md` | `src/Materialmatrix.jl` | 4 |
| `SetMesh.md` | `src/SetMesh.jl` | 5 |
| `Jellyrollmodel.md` | `src/Jellyrollmodel.jl` | 5 |
| `ring.md` | `src/ring.jl` | 5 |
| `CouplingState.md` | `src/CouplingState.jl` | 6 |
| `Variables.md` | `src/Variables.jl` | 6 |
| `Initialisation.md` | `src/Initialisation.jl` | 6 |
| `PostProcessing.md` | `src/PostProcessing.jl` | 6 |
| `Solve.md` | `src/Solve.jl` | 7 |
| `CallModel.md` | `src/CallModel.jl` | 7 |
| `CycleSolver.md` | `src/CycleSolver.jl` | 7 |
| `CycleData.md` | `src/CycleData.jl` | 7 |
| `Parallelsolution.md` | `src/Parallelsolution.jl` | 7 |
| `Assemble.md` | `src/Assemble.jl` | 7 |
| `CsvExport.md` | `src/CsvExport.jl` | 8 |
| `Tools.md` | `src/Tools.jl` | 8 |

**不修改的现有文件**：任何 `src/*.jl` 源码、`md/00-15` 现有技术文档、CLAUDE.md。

---

## 单文件生成模板（适用于 Task 1.x ~ 8.x）

每个源文件任务统一遵循以下步骤。**`<SRC>` 是源文件相对路径**（如 `src/CouplingState.jl`），**`<NAME>` 是去 `.jl` 后缀的基名**（如 `CouplingState`），**`<OUT>` 是输出路径**（`md/源码函数索引/<NAME>.md`）。

### 步骤

- [ ] **S1: 完整读取源文件 `<SRC>`**

  Run: 用 Read 工具读取 `<SRC>`（务必读完整文件，不限于前 2000 行）

- [ ] **S2: 提取所有顶层定义**

  识别所有 `function`、`mutable struct`、`struct`、`const`，记录每个的：
  - 名称、签名、起止行号
  - 是否为 trivial（一行 getter、纯 re-export、字面常量）→ 归入「省略」组

- [ ] **S3: 识别三类标注**

  按 spec §4 启发式规则全文扫描 `<SRC>`，记录命中行号 + 内容片段 + 推测说明：
  - **[DEBUG]**: `println(` / `print(` / `@show` / `@info`（**结构化日志除外**：结构化日志 = `@info` 出现在 `finally` 块、函数末尾、或 phase/cycle/初始化的固定位置且消息含这些语义；其余 `@info` 一律标 [DEBUG]）/ 注释含 `debug/调试/临时/temp`
  - **[PLACEHOLDER]**: `TODO/FIXME/HACK/XXX` 注释；注释含 `placeholder/fallback/hardcoded/占位/兜底/临时/防止/避免`；裸 `NaN/Inf`；try-catch 静默吞错（catch 块只 return 默认值，**无** `@warn`/`@error`/`println`）；魔数 + 兜底注释组合（**以相邻注释为主要信号**：注释含上述关键词 → 标；无此类注释 → 默认不标）
  - **[COMPLEX-CHECK]**: `&&` 链 ≥3 条件；嵌套 `if` ≥3 层；连续 ≥2 个 `hasproperty/isdefined/!== nothing` 检查；单条件表达式 >100 字符

  每条标注需独立判断是否真符合（如：物理合理的零初值不算 PLACEHOLDER）。
  **遇到歧义**：在「用途推测/风险」列写明歧义点并标 `[UNCERTAIN]` 前缀，不强行归类。

- [ ] **S4: 生成 md**

  按 spec §3.1 模板写 `<OUT>`，包含以下小节（无则显式写"无"）：
  1. **文件头**：源文件路径、规模（用 `wc -l` 数据）、函数/struct 计数、职责一句话、相关技术文档链接
  2. **数据结构**：每个 struct 一条目（字段、用途）
  3. **函数清单**：每个非 trivial 函数一独立条目（签名+行号、职责、关键逻辑 bullets、跨文件依赖）
  4. **省略项**：合并 trivial getter/常量为一行说明
  5. **[DEBUG] 表**：表头固定为 `### [DEBUG]`，列 = 行号 | 内容 | 用途推测
  6. **[PLACEHOLDER] 表**：表头固定为 `### [PLACEHOLDER]`，列 = 行号 | 内容 | 风险
  7. **[COMPLEX-CHECK] 表**：表头固定为 `### [COMPLEX-CHECK]`，列 = 行号 | 内容 | 简化建议

  **数据行格式严格统一**：以 `| L<数字> |` 开头（如 `| L142 | println(...) | 临时检查 |`），便于 Task 8.3 awk 脚本聚合。

  **边缘情况**：若文件无 function/struct 定义（纯 include/export，如 `JuBat.jl`），第 2-3 节写「本文件无独立函数/struct 定义」，第 1 节「职责」字段说明实际内容（include 顺序、export 列表）；若某标注表为空，写「无」并保留表头。

  写完用 Write 工具保存到 `<OUT>`。

- [ ] **S5: 验证行号引用**

  抽查 3 处行号引用：随机挑 md 里的 3 个 `L<数字>` 引用，回到 `<SRC>` 看该行确实是对应内容。如有偏移，修正 md。

- [ ] **S6: 提交**

  ```bash
  git add <OUT>
  git commit -m "docs(src-index): add <NAME>.md"
  ```

---

## Chunk 1: 批次 1 — 入口与参数

### Task 1.0: 创建目录骨架

**Files:**
- Create: `md/源码函数索引/.gitkeep`（占位，确保目录存在）

- [ ] **Step 1: 创建目录**

  Run:
  ```bash
  mkdir -p "md/源码函数索引"
  touch "md/源码函数索引/.gitkeep"
  ```

- [ ] **Step 2: 提交**

  ```bash
  git add "md/源码函数索引/.gitkeep"
  git commit -m "docs(src-index): scaffold directory"
  ```

### Task 1.1: JuBat.md

套用【单文件生成模板】，`<SRC>` = `src/JuBat.jl`，`<NAME>` = `JuBat`，`<OUT>` = `md/源码函数索引/JuBat.md`。

**注意**：这是模块入口文件，主要是 `include` 和 `export`，函数定义极少。重点说明模块加载顺序、对外 export 的 API 清单（可分组列出）。

### Task 1.2: Option.md

套用模板，`<SRC>` = `src/Option.jl`，`<NAME>` = `Option`，`<OUT>` = `md/源码函数索引/Option.md`。

**注意**：主要是 `Option` struct + 字段说明。每个字段已部分在 CLAUDE.md §5 列出，**不要照抄**，要按字段在源码的真实顺序、补上 CLAUDE.md 没列的字段（如 `czm_*` 系列），含默认值与单位。

### Task 1.3: SetParams.md

`<SRC>` = `src/SetParams.jl`。包含 `Electrode/Cell/Cohesive/Scale/Params` struct。重点列各 struct 字段及其归一化含义（参考 `md/01_参数定义与归一化.md`），交叉链接该文档。

### Task 1.4: SetCase.md

`<SRC>` = `src/SetCase.jl`。包含 `SetCase`、`ChooseCell`、`NormaliseParam` 等。重点说明 `case` 对象构建流程、归一化入口。

---

## Chunk 2: 批次 2 — 电化学

### Task 2.1: SPMe.md
`<SRC>` = `src/SPMe.jl`。SPMe 主体。交叉链接 `md/04_电化学模型_SPMe.md`。

### Task 2.2: SPM.md
`<SRC>` = `src/SPM.jl`。

### Task 2.3: P2D.md
`<SRC>` = `src/P2D.jl`。伪二维模型。

### Task 2.4: ElectrodeDiffusion.md
`<SRC>` = `src/ElectrodeDiffusion.jl`。注意此文件极短（~30 行），可能只有 1 个函数 — 函数清单可能只有 1 条，标注可能为"无"。

### Task 2.5: ElectrodePotential.md
`<SRC>` = `src/ElectrodePotential.jl`。同上，极短。

### Task 2.6: ElectrolyteDiffusion.md
`<SRC>` = `src/ElectrolyteDiffusion.jl`。

### Task 2.7: ElectrolytePotential.md
`<SRC>` = `src/ElectrolytePotential.jl`。

---

## Chunk 3: 批次 3 — 热模型

### Task 3.1: Thermal.md
`<SRC>` = `src/Thermal.jl`。交叉链接 `md/05_热模型_二维分布式.md`。

### Task 3.2: ThermalDistributed.md
`<SRC>` = `src/ThermalDistributed.jl`。含 10 个函数，包括边界条件、各向异性导热、分层热源。重点说明分层热源计算路径。

### Task 3.3: ThermalPolar2D.md
`<SRC>` = `src/ThermalPolar2D.jl`。极坐标 FVM。

---

## Chunk 4: 批次 4 — 力学 / CZM

### Task 4.1: czm.md
`<SRC>` = `src/czm.jl`。CZM 本构 + 本构参数 + `create_czm_mesh`。交叉链接 `md/06_内聚力模型_CZM.md`。13 个函数/struct，文档体量较大，预计 ~150 行 md。

### Task 4.2: CzmSolve.md
`<SRC>` = `src/CzmSolve.jl`。Newton-Raphson 求解器，13 个函数。重点说明迭代收敛准则、面积缩减处理（`czm_area_loss_*`）。

### Task 4.3: CzmPostProcess.md
`<SRC>` = `src/CzmPostProcess.jl`。

### Task 4.4: CzmUnitMesh.md
`<SRC>` = `src/CzmUnitMesh.jl`。单元测试用 strip mesh（参考 `test/unit_czm_*.jl`）。

### Task 4.5: Mechanical.md
`<SRC>` = `src/Mechanical.jl`。应力计算。注意区分颗粒 vs 极片模量（CLAUDE.md §9.4 / `md/15_颗粒与极片模量区分.md`）。

### Task 4.6: Materialmatrix.md
`<SRC>` = `src/Materialmatrix.jl`。本构矩阵构造，12 个函数。

---

## Chunk 5: 批次 5 — 网格与几何

### Task 5.1: SetMesh.md
`<SRC>` = `src/SetMesh.jl`。项目最长文件之一（739 行、16 个函数）。重点说明各 `setup_*_mesh` 函数与 `CohesiveMesh` 结构。交叉链接 `md/02_几何与网格.md`。

### Task 5.2: Jellyrollmodel.md
`<SRC>` = `src/Jellyrollmodel.jl`。阿基米德螺旋线、collector-seeded 网格。8 个函数。

### Task 5.3: ring.md
`<SRC>` = `src/ring.jl`。

---

## Chunk 6: 批次 6 — 状态与变量

### Task 6.1: CouplingState.md
`<SRC>` = `src/CouplingState.jl`。18 个函数/struct，跨周期状态管理。重点说明损伤累积、SOH 更新、面积缩减触发逻辑。

### Task 6.2: Variables.md
`<SRC>` = `src/Variables.jl`。`StandardVariables` / `create_element_workspace` / `copy_element_results`。重点列出 hardcoded 键表（已知技术债，参考 MEMORY.md）。

### Task 6.3: Initialisation.md
`<SRC>` = `src/Initialisation.jl`。状态初始化。

### Task 6.4: PostProcessing.md
`<SRC>` = `src/PostProcessing.jl`。结果提取与归一化还原。

---

## Chunk 7: 批次 7 — 求解器

### Task 7.1: Solve.md
`<SRC>` = `src/Solve.jl`。主求解器 + `CallModel_MultiSPMe`。交叉链接 `md/09_分流求解器.md`。

### Task 7.2: CallModel.md
`<SRC>` = `src/CallModel.jl`。多 SPMe 调用框架。

### Task 7.3: CycleSolver.md
`<SRC>` = `src/CycleSolver.jl`。循环求解器、`PhaseResult/CycleResult`。

### Task 7.4: CycleData.md
`<SRC>` = `src/CycleData.jl`。

### Task 7.5: Parallelsolution.md
`<SRC>` = `src/Parallelsolution.jl`。分流求解器、`solve_branch_currents_newton`。交叉链接 `md/09_分流求解器.md`。

### Task 7.6: Assemble.md
`<SRC>` = `src/Assemble.jl`。短文件，~40 行。

---

## Chunk 8: 批次 8 — IO 与工具 + 索引

### Task 8.1: CsvExport.md
`<SRC>` = `src/CsvExport.jl`。25 KB / 16 个函数。重点列出 `variables[...]` 字符串键用法（已知技术债）。

### Task 8.2: Tools.md
`<SRC>` = `src/Tools.jl`。通用工具函数。

### Task 8.3: 生成 `_索引.md`

**Files:**
- Create: `md/源码函数索引/_索引.md`

**关键约定**：所有源文件 md 的三类标注表必须使用固定表头 `### [DEBUG]`、`### [PLACEHOLDER]`、`### [COMPLEX-CHECK]`（已在 Task 模板 S4 规定），数据行以 `^| L` 开头。本任务基于此确定性结构进行聚合。

- [ ] **Step 1: 用脚本确定性聚合各 md 的三类标注计数**

  Run（bash，在仓库根目录）：
  ```bash
  for f in "md/源码函数索引"/*.md; do
    [[ "$f" == *_索引.md ]] && continue
    name=$(basename "$f" .md)
    # 抽取每个 ### 表节的数据行数：awk 进入/退出模式
    counts=$(awk '
      /^### \[DEBUG\]/         {sec="D"; next}
      /^### \[PLACEHOLDER\]/   {sec="P"; next}
      /^### \[COMPLEX-CHECK\]/ {sec="C"; next}
      /^### /                  {sec=""; next}
      /^---/                   {sec=""; next}
      /^\| L/ && sec=="D" {d++}
      /^\| L/ && sec=="P" {p++}
      /^\| L/ && sec=="C" {c++}
      END {printf "D=%d P=%d C=%d", d+0, p+0, c+0}
    ' "$f")
    echo "$name: $counts"
  done
  ```
  预期输出（每行一份）：
  ```
  JuBat: D=0 P=0 C=0
  CouplingState: D=3 P=2 C=5
  ...
  ```
  将输出原样记录到下一步的统计表里。**若某 md 的 D/P/C 计数与肉眼抽查明显不符**（抽查 2 份），先回到该 md 修正表头/数据行格式，重跑脚本。

- [ ] **Step 2: 按物理模块组织 TOC**

  按 spec §3.2 模板，分 7 个模块组列出所有 35 个源文件 md 链接：
  - 电化学（SPMe/SPM/P2D/Electrode×2/Electrolyte×2）
  - 热模型（Thermal/ThermalDistributed/ThermalPolar2D）
  - 力学·CZM（czm/CzmSolve/CzmPostProcess/CzmUnitMesh/Mechanical/Materialmatrix）
  - 网格·几何（SetMesh/Jellyrollmodel/ring）
  - 参数·状态（Option/SetParams/SetCase/CouplingState/Variables/Initialisation/PostProcessing）
  - 求解器（Solve/CallModel/CycleSolver/CycleData/Parallelsolution/Assemble）
  - 入口·IO·工具（JuBat/CsvExport/Tools）

- [ ] **Step 3: 写三类标注全局统计表 + 优先处理建议**

  **3a 总数表**：三类总数合计（基于 Step 1 脚本输出求和）。

  **3b Top 5 热点文件表**：按 `(D + P + C)` 降序取前 5，列：文件名 | D | P | C | 合计。

  **3c 优先处理建议**（确定性规则）：
  - **[PLACEHOLDER] 高风险项**：列出所有 [PLACEHOLDER] 条目中，**「风险」列含以下任一关键词**的条目：`掩盖 / coupling / 耦合 bug / nothing / 未就绪 / 默认值 / NaN`。每条带「文件名:行号 + 风险描述」。
  - **[COMPLEX-CHECK] Top 5 简化候选**：按表达式**字符长度降序**取前 5（字符长度 = md 中 `^| L<n> | <内容> |` 一行的「内容」字段长度）。每条带「文件名:行号 + 简化建议」。

  文件顶部注明「以 commit `<HASH>` 为准」（运行 `git rev-parse --short HEAD` 取当前 hash）。

- [ ] **Step 4: 写入文件并提交**

  ```bash
  git add "md/源码函数索引/_索引.md"
  # 删除 Task 1.0 创建的占位文件
  git rm "md/源码函数索引/.gitkeep"
  git commit -m "docs(src-index): add _索引.md and finalize directory"
  ```
  注：`.gitkeep` 由 Task 1.0 创建，必须在此处删除以满足验收「无 `.gitkeep`」。

---

## 验收（Plan-level）

执行完所有任务后：

- [ ] `ls "md/源码函数索引/"` 应见 35 份 `<NAME>.md` + 1 份 `_索引.md`，无 `.gitkeep`
- [ ] `_索引.md` 含 commit hash 标注 + Top 5 表 + 优先处理建议
- [ ] 随机抽 5 份 md 验证行号引用准确（每份抽 3 处；其中至少 1 份必须是源码 >400 行的文件，如 SetMesh/CouplingState/czm，因长文件行号偏移风险高）
- [ ] 随机抽 5 处三类标注核对源码（无误报重大问题）

Task 8.3 的最终 commit 即为完成标记，无需额外空 commit。

---

## 风险与边界

| 风险 | 应对 |
|------|------|
| 单文件 md 体量大（如 SetMesh、CzmSolve）→ 单 commit 大 | 允许；不拆分单个 md |
| 行号引用偏移 | 每任务 S5 强制抽查 3 处 |
| 三类标注误报 | 文档统一声明"仅供参考、需人工复核" |
| 部分源文件极短（如 Electrode*.jl）→ 函数清单可能为空 | 在文档里显式说明"该文件仅含 X，无独立函数" |
| 工作量 30 个文件 | 按 chunk 提交，每个 chunk 可独立交付 |

---

## 后续工作（不在本计划范围）

- 在 CLAUDE.md「文档索引」处添加新目录链接（待执行完毕后单独 commit）
- 维护策略：新增函数/源文件时如何同步更新（待落地后讨论）
