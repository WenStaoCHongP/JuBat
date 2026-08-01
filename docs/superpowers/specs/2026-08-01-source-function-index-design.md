# 源码函数索引文档（src/ 函数级注释 + 三类标注）

**日期**: 2026-08-01
**作者**: 与 Claude Code 协作
**状态**: 设计稿（待评审）

---

## 1. 背景与目标

### 1.1 背景

JuBat 项目的 `src/` 目录共 36 个 `.jl` 文件（含 `parameters/` 子目录中 5 个参数集文件）、约 203 个函数/struct、约 11,184 行代码。现有 `md/` 目录的 17 篇技术文档以**物理模型/算法**为主线组织（电化学、热、CZM、耦合等），并不按源文件/函数粒度展开。

这导致两个问题：

1. **代码定位困难** — 看到 md 里的概念，找不到对应源码函数；看到源码函数，没有"它做什么、为什么这么写"的高层说明。
2. **技术债隐蔽** — 代码里散落着为调试临时加的 `println`、为绕过缺失数据而硬写的兜底赋值、深层嵌套的属性判断。这些是隐患但从未被汇总。

### 1.2 目标

为 `src/` 每个 `.jl` 文件（`parameters/` 子目录除外）生成一份 markdown，包含：

- **函数级职责注释** — 每个非 trivial 函数的用途、关键逻辑、依赖关系
- **三类标注** — 显式标记并解释：
  - `[DEBUG]` 为调试临时添加的语句
  - `[PLACEHOLDER]` 为保证程序运行而添加的兜底赋值
  - `[COMPLEX-CHECK]` 复杂的属性/条件判断

**非目标**：
- 不修改任何源代码（注释只写进 md，源码保持原样）
- 不覆盖 `src/parameters/`（参数集文件，函数极少）
- 不替换现有 `md/00-15` 技术文档，仅作为代码级索引补充

---

## 2. 总体结构

### 2.1 目录布局

新建 `md/源码函数索引/` 子目录：

```
md/
├── 00_文档索引与命名规范.md          (现有)
├── 01_参数定义与归一化.md            (现有)
├── ...
├── 15_颗粒与极片模量区分.md          (现有)
└── 源码函数索引/                     ← 新建
    ├── _索引.md                      ← 总入口（TOC + 三类标注全局统计）
    ├── JuBat.md                      ← 对应 src/JuBat.jl
    ├── CouplingState.md              ← 对应 src/CouplingState.jl
    ├── czm.md                        ← 对应 src/czm.jl
    ├── CzmSolve.md
    ├── ... (共 35 个文件 + _索引.md)
```

**命名约定**：md 文件名 = 对应 src 文件去 `.jl` 后缀（如 `CouplingState.jl` → `CouplingState.md`），便于双向跳转。

### 2.2 覆盖范围

| 路径 | 是否覆盖 | 说明 |
|------|----------|------|
| `src/*.jl`（35 个根目录核心文件，不含 `install.jl`） | ✓ | 每个文件单独一份 md |
| `src/parameters/*.jl`（5 个） | ✗ | 参数集文件，函数极少，跳过 |
| `src/install.jl` | ✗ | 仅依赖声明（107 字节） |
| `src/.DS_Store` | ✗ | 系统文件 |

合计覆盖：35 个根目录源文件 + 1 个 `_索引.md` = 36 份新文档。

### 2.3 与现有文档的关系

- 现有 `md/00-15`：按**物理模型/算法**主线（"SPMe 怎么实现"、"CZM 本构是什么"）
- 新 `md/源码函数索引/`：按**源文件/函数**主线（"这个文件里有哪些函数、每个做什么"）
- 两者互补，不冲突。在 `_索引.md` 中会显式交叉引用现有文档（例如 `czm.md` 顶部链接到 `md/06_内聚力模型_CZM.md`）。

---

## 3. 文档结构

### 3.1 单文件 md 模板

```markdown
# CouplingState.jl

**源文件**: `src/CouplingState.jl` · **~33 KB / 800 行** · **18 个函数/struct**
**职责**: 维护电-热-CZM 耦合过程中的跨周期状态（损伤累积、SOC 传递、温度记忆等）。
**相关技术文档**: `md/10_参数传递与模块架构.md`

---

## 数据结构

### `CouplingState`  (mutable struct, L12-L48)
跨周期耦合状态的容器。
- 字段 `D_history`: 各 cohesive 单元的累积损伤向量
- 字段 `T_init`: 每相初始温度
- 字段 `soh`: 当前健康状态（0-1）

---

## 函数清单

### `init_coupling_state(case)`  L52-L78
**职责**: 从 `case` 初始化耦合状态，分配各单元损伤/SOC 数组。
**关键逻辑**:
- 调用 `setup_damage_arrays` 生成初始损伤向量
- 若 `opt.czm_enabled == false`，跳过损伤分配
**依赖**: `SetParams.jl::get_cell_geometry`, `Variables.jl::create_workspace`

### `accumulate_damage!(state, czm_result)`  L120-L155
**职责**: 把单步 CZM 损伤增量累加到 `state.D_history`，并更新 SOH。
**关键逻辑**:
- D > 0.83 时触发面积缩减逻辑（见 `czm_area_loss_enabled`）
- 累加后断言 `all(0 .<= state.D_history .<= 1)`

> *省略: `get_D(state)` 等 4 个一行 getter，见源文件 L210-L225*

---

## 三类标注

### [DEBUG] 调试语句
| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L142 | `println("D_step = ", D_step)` | 临时损伤增量检查 |
| L287 | `@show "T_e check", T_e` | 单元温度调试 |

### [PLACEHOLDER] 运行兜底
| 行号 | 内容 | 风险 |
|------|------|------|
| L95 | `q_e fallback = 0.0` | SPMe 未就绪时热源置零，可能掩盖耦合 bug |

### [COMPLEX-CHECK] 复杂属性判断
| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L168 | `hasproperty(case, :opt) && case.opt !== nothing && hasproperty(case.opt, :czm_enabled) && case.opt.czm_enabled` | 抽取 `is_czm_enabled(case)` 助手 |
```

### 3.2 函数条目颗粒度规则

**单独列条目**：
- struct（mutable/immutable struct、`mutable struct`、`struct`）
- ≥3 行的函数
- 任何包含分支、循环、状态修改的函数

**合并到「省略」说明**：
- 一行 getter（如 `get_D(state) = state.D`）
- 纯 re-export（`export foo, bar`）
- 字面常量赋值（`const X = 1.0`）

每个独立条目包含：
- **签名 + 行号范围**（如 `L52-L78`）
- **职责**：1-2 句话
- **关键逻辑**：bullet 列出非平凡步骤（条件分支、副作用、关键调用）
- **依赖**：跨文件调用（仅跨文件的，文件内部不重复列）

### 3.3 `_索引.md` 结构

```markdown
# src/ 源码函数索引

JuBat 项目 `src/` 目录函数级注释与三类标注汇总。

## 按物理模块浏览

### 电化学
- [SPMe.md](SPMe.md) — SPMe 模型主体
- [SPM.md](SPM.md) — 单颗粒模型
- [P2D.md](P2D.md) — 伪二维模型
- ...

### 热模型
- [ThermalDistributed.md](ThermalDistributed.md)
- [ThermalPolar2D.md](ThermalPolar2D.md)
- ...

### 力学 / CZM
- [czm.md](czm.md) — CZM 本构与本构参数
- [CzmSolve.md](CzmSolve.md) — Newton-Raphson 求解
- [Mechanical.md](Mechanical.md) — 应力计算
- ...

### 网格 / 几何
### 参数与状态
### 求解器
### IO / 后处理

---

## 三类标注全局统计

共发现 **N1** 处 [DEBUG]、**N2** 处 [PLACEHOLDER]、**N3** 处 [COMPLEX-CHECK]。

### Top 5 热点文件（按三类标注总数）
| 文件 | DEBUG | PLACEHOLDER | COMPLEX-CHECK | 合计 |
|------|-------|-------------|---------------|------|
| CouplingState.jl | 3 | 2 | 5 | 10 |
| ... |

### 建议优先处理
- **PLACEHOLDER 高风险项**: ...(若存在)
- **COMPLEX-CHECK 简化候选**: ...
```

---

## 4. 三类标注识别规则（启发式）

### 4.1 [DEBUG] 调试语句

| 触发模式 | 示例 |
|----------|------|
| `println(` / `print(` / `@show` / `@info` 调用 | `println("T_e = ", T_e)` |
| 注释含关键词 | `# debug`、`# 调试`、`# 临时`、`# temp`、`# for debug` |
| 明显临时的中间变量打印 | `@show "checkpoint 1", x` |

**排除**：
- **结构化日志**：`@info` 出现在 `finally` 块、函数末尾、phase/cycle/初始化的固定位置且消息含这些语义 → 不算 [DEBUG]。**判定准则**：若该输出语句删除后程序逻辑不变，且非 phase/循环/初始化的固定结构化日志，则判为 [DEBUG]；否则排除。
- **用户配置字段**：名为 `debug_*`/`verbose_*`/`log_*` 的 struct 字段（如 `Option.debug_coupling`）是用户主动开启的诊断开关，**不属于** [DEBUG]。它们是合法 API。
- **参数校验警告**：`@warn` 配合参数完整性检查（缺失字段、范围越界）属于结构性校验，不算 [DEBUG]。但若 `@warn` 后跟硬编码兜底值，则该**兜底赋值**可能属于 [PLACEHOLDER]。

### 4.2 [PLACEHOLDER] 运行兜底

| 触发模式 | 示例 |
|----------|------|
| `TODO`/`FIXME`/`HACK`/`XXX` 注释 | `# TODO: 用真实 E_coat 替换` |
| 注释含 `placeholder`/`fallback`/`hardcoded`/`hardcode`/`占位`/`兜底`/`临时` | |
| 裸 `NaN`/`Inf` 兜底赋值 | `q_e = NaN  # SPMe 未就绪` |
| try-catch 静默吞错（catch 块只 return 默认值，无日志） | `catch; return 0.0 end` |
| 配合魔数 0.0/1.0 的"避免崩溃"赋值（**以相邻注释为主要信号**：注释含 `fallback/防止/避免/保证运行/兜底` → 判为 [PLACEHOLDER]；无此类注释 → 默认不标） | `E_coat = 5e8  # 防止 nothing` |
| try-catch 静默吞错（catch 块只 return 默认值，**无** `@warn`/`@error`/`println`） | `catch; return 0.0 end` |

**排除**：物理上确实合理的零初值（如 `D = zeros(n)` 初始损伤为 0）——这是物理初值，不是兜底。

### 4.3 [COMPLEX-CHECK] 复杂属性判断

| 触发模式 | 示例 |
|----------|------|
| `&&` 链 ≥3 个条件 | `a !== nothing && hasproperty(a, :b) && a.b > 0 && isfinite(a.b)` |
| 嵌套 `if` ≥3 层 | |
| 连续 ≥2 个 `hasproperty`/`isdefined`/`!== nothing`/`isa` 检查 | |
| 单表达式长度 >100 字符的条件判断 | |

**简化建议生成规则**：
- `&&` 链检查同一对象的多个属性 → 建议抽取助手函数（如 `is_czm_enabled(case)`）
- 深层 `hasproperty`/`isdefined` 链 → 建议用 `try`-`catch` 或重构数据结构
- 嵌套 `if` → 建议早返回（guard clause）

### 4.4 标注可信度

每个标注条目含「用途推测/风险/简化建议」字段，由生成者**基于上下文推断**填写，不保证 100% 准确。文档顶部统一声明：

> 标注基于启发式规则自动识别 + 上下文推测，**仅供定位参考**，需人工复核后再决定是否修改源码。

---

## 5. 工作流程

### 5.1 阅读顺序（按物理模块分批）

按以下批次生成文档，每批完成后可独立交付：

| 批次 | 文件 | 主线 |
|------|------|------|
| 1 | `JuBat.jl`, `Option.jl`, `SetParams.jl`, `SetCase.jl` | 入口与参数 |
| 2 | `SPMe.jl`, `SPM.jl`, `P2D.jl`, `ElectrodeDiffusion.jl`, `ElectrodePotential.jl`, `ElectrolyteDiffusion.jl`, `ElectrolytePotential.jl` | 电化学 |
| 3 | `Thermal.jl`, `ThermalDistributed.jl`, `ThermalPolar2D.jl` | 热模型 |
| 4 | `czm.jl`, `CzmSolve.jl`, `CzmPostProcess.jl`, `CzmUnitMesh.jl`, `Mechanical.jl`, `Materialmatrix.jl` | 力学 / CZM |
| 5 | `SetMesh.jl`, `Jellyrollmodel.jl`, `ring.jl` | 网格与几何 |
| 6 | `CouplingState.jl`, `Variables.jl`, `Initialisation.jl`, `PostProcessing.jl` | 状态与变量 |
| 7 | `Solve.jl`, `CallModel.jl`, `CycleSolver.jl`, `CycleData.jl`, `Parallelsolution.jl`, `Assemble.jl` | 求解器 |
| 8 | `CsvExport.jl`, `Tools.jl` | IO 与工具 |
| 9 | `_索引.md` | 汇总 |

### 5.2 单文件工作流

1. **完整读取**该 `.jl` 文件
2. **提取顶层定义**：所有 `function`/`struct`/`mutable struct`/`const` + 行号
3. **分类**：独立条目 vs 合并省略
4. **识别三类标注**：按第 4 节规则扫描全文
5. **撰写对应 md**：按第 3.1 节模板
6. **交叉验证**：行号引用回查源码（避免行号偏移）

### 5.3 验收标准

每个生成的 md 文件：
- [ ] 顶部含元信息（源文件路径、规模、职责、相关技术文档链接）
- [ ] 所有非 trivial 函数/struct 有独立条目
- [ ] trivial getter/常量合并到省略说明
- [ ] 三类标注表格完整（无则显式写"无"）
- [ ] 所有行号引用准确（抽查 ≥3 处与源码一致）

---

## 6. 风险与权衡

### 6.1 已识别风险

| 风险 | 缓解 |
|------|------|
| 行号引用易过期（源码改动后失效） | 每个 md 顶部标注「以 `<commit-hash>` 为准」；建议后续在 CLAUDE.md 提醒改源码时同步更新 |
| 三类标注误报（启发式过于激进） | 文档显式声明"仅供参考、需人工复核"；[PLACEHOLDER] 排除物理合理的初值 |
| 函数职责概括不准确 | 每条目含「关键逻辑」bullet，便于用户快速校验 |
| 工作量大（35 个文件） | 按 5.1 分批交付，每批独立可用 |

### 6.2 显式放弃的方案

- **回写源码标记**（`# [DEBUG]` 等）：用户已确认不修改源码
- **覆盖 `parameters/`**：函数极少，价值低
- **自动生成（AST 解析）**：能提取签名但读不出"职责"和"三类标注"，必须人工阅读

---

## 7. 后续工作（不在本 spec 范围）

- 在 CLAUDE.md「文档索引」处添加新目录链接
- 维护策略：新增函数/源文件时如何同步更新（待文档落地后讨论）
- 是否把三类标注做成可机器检查的 lint 规则（长期方向）

---

## 8. 验收

完成时交付：
- [ ] `md/源码函数索引/` 目录及 35 份源文件 md + 1 份 `_索引.md`
- [ ] `_索引.md` 含 TOC、三类标注全局统计、记录生成时所基于的 git commit hash
- [ ] 抽查 3 个文件，所有行号引用准确
- [ ] 抽查 5 处三类标注，识别合理（无误报/漏报重大问题）
