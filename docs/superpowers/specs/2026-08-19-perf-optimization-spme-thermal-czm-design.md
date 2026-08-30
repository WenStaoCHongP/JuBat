# SPMe-分布式热-CZM 支线仿真提速优化设计规格

> 日期: 2026-08-19
> 状态: 已实施（2026-08-19）。结论：批次 0 裁剪后仅 Task 7 开门并已合并（四判据 bit 级一致，CZM -1.7%）；profile 再定位真正瓶颈为 CzmSolve.jl:222 每迭代稀疏 LU 分解，增补方案（K_bc 因子内容判据复用）见 docs/planning-with-files/29_仿真提速-SPMe热CZM/findings.md
> 范围: 仅 SPMe (`per_element_spme=true`) - 分布式热 (`distributed2D`) - CZM 耦合支线
> 关联既有工作: `2026-04-20-czm-vectorized-solver-design.md`（实施中，本计划 D 组为其延续）、`2026-04-02-initialisation-optimization-design.md`（已落地，`MultiSPMeLayout`/`MeshGeometry` 为本计划缓存挂载点）

---

## 1. 背景与目标

对 SPMe-分布式热-CZM 耦合支线执行**纯工程性能优化**：在不改变任何浮点运算数值与顺序的前提下，通过缓存、预分配、消除重复计算降低仿真墙钟时间。

**目标场景**: `example/testexample.jl`（nθ=80、60 s 全耦合），自带 `timing_totals` 分段计时（spme/branch/thermal/czm）作为各模块提速效果的验收依据。

**不做的事**（YAGNI 裁剪）:
- 不改时间步自适应逻辑来减少步数（基线网格/步数必须不变）
- 不换线性求解器类型（直接法 → Krylov 等被数值约束排除）
- 不做多线程化（9.6 基线约定单线程运行，优化必须在单线程下见效）
- 不做循环仿真场景专项优化（可作为后续独立计划）

## 2. 约束与已确认决策

| 决策点 | 内容 | 来源 |
|--------|------|------|
| 数值等价性 | **基线 bit 级一致**：只做缓存/预分配/消除重复计算；浮点运算数值与顺序不变；任何会重排浮点的手段单独报批 | 用户 2026-08-19 确认 |
| 验收判据 | 退出码、网格/步数、`metrics.toml` 科学结果（按记录精度）、PNG SHA-256（**本机前后对比**，不用跨机档案 SHA） | AGENTS.md 9.6 + 记忆"基线机器陷阱" |
| 运行约定 | Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no` | AGENTS.md 9.6 |
| 实施策略 | **方案一：测量驱动 + 收益门槛**：批次 0 实测占比后，只为占比 ≥10% 的模块开批次 | 用户 2026-08-19 确认 |
| 占比口径 | `debug_coupling=true` 的 `timing_totals` 为主，Profile 标准库火焰图交叉确认（debug 计时自身有开销） | — |
| 门槛值 | 10%（占基线墙钟）；低于门槛的嫌疑点记入 findings 不开批次 | 用户可调 |

## 3. 静态热点排查结论（批次 0 之前的工作假设）

按嫌疑排序（来源：2026-08-19 代码探索，file:line 以当前分支为准）:

1. **CZM Newton 求解**: 每电化学步一次，每次 Newton 迭代全量重组装 + 稀疏直接求解，线搜索内再组装；`K_total = K_bulk + K_coh` 每迭代分配新矩阵（`src/CzmSolve.jl:194-235`）。K_bulk 已由 `CZMAssemblyCache` 缓存，但 cohesive 部分、K_total、f_int_bulk 每迭代重建。
2. **热矩阵每步 5 次 `Assemble` 全量重建**（`src/ThermalDistributed.jl:19-37, 160-196`）: MT/KT 稀疏模式与大部分系数为不变量；BC 每步先 `copy`；`apply_cool_method` 逐高斯点 `K[ni,nj] -=`。
3. **逐单元 SPMe 循环每步重分配**（`src/CallModel.jl:120-131`、`src/SPMe.jl:118-121`）: `ws_pool` 每步重建、`copy_element_results` 每单元新建 18 键 Dict、`case.index` keys 收集 + filter 每单元每步执行。
4. **两级 blockdiag 拼装 + 每步稀疏 `\` 直接法**（`src/CallModel.jl:136-186`、`src/Solve.jl:214`）: 结构恒定但每步 splat 重建、无分解复用。
5. **分流求解器**（`src/Parallelsolution.jl:200-294`）: 已预分配、较轻，低嫌疑。

该排序仅为工作假设，由批次 0 实测裁决。

## 4. 批次 0：测量与本机基线锚点

**批次 0 做四件事:**

1. **本机基线锚点**: 当前 HEAD 按 9.6 约定运行 `example/testexample.jl`，归档墙钟时间、`metrics.toml`、PNG SHA-256、`timing_totals`。**重跑两次验证 PNG SHA 本机稳定性**——不稳定则验收改用 metrics + 网格/步数 + 退出码三判据（并记录）。
2. **分段占比**: `debug_coupling=true` 跑一次取 `timing_totals` 四项占比；Profile 标准库对非 debug 运行采样火焰图交叉确认。
3. **两个关键事实核查**（决定收益上限）:
   - **3a 热系统矩阵系数是否随步变化**: 重点查界面热阻 `k_eff(D)`（`md/07_界面热阻模型.md`）——若损伤每步更新且热阻随之变化，KT 非常量，A2/C2 从"缓存分解因子"降级为"`lu!(F, A)` 复用符号结构、每步重算数值分解"。方法：运行中每步打印 `norm(KT_t - KT_{t-1})`。
   - **3b CZM basic 法实际迭代数**: Newton 迭代 × 线搜索重组装次数，决定 D 组值不值得开批。
4. **CZM 既有计划残留核对**: 对照 `docs/superpowers/plans/2026-04-20-czm-vectorized-solver-plan.md` 与 `docs/planning-with-files/09_向量化CZM/`、`docs/planning-with-files/01_CZM瓶颈/`，列出已完成/未完成项，D 组只做残留项之外的部分，避免重复或冲突。

**输出**: `docs/planning-with-files/29_仿真提速-SPMe热CZM/`（按 AGENTS.md 9.5 约定，三件套由 planning-with-files 技能管理）下写 findings.md：实测占比表、3a/3b 答案、CZM 残留核对结论、据此裁剪的批次清单；同步更新 `docs/planning-with-files/index.md`。

## 5. 候选优化手段目录（批次 0 从中裁剪）

每项附 bit 一致性论证。**全部手段不改变任何浮点运算的数值与顺序**。

### A. 热模块（`src/ThermalDistributed.jl`）

- **A1 矩阵不变量缓存**: MT/KT 系数不随步变化时只算一次缓存复用；`ThermalDistributed2D_BC` 每步 `copy(KT)`/`copy(FT)` 改预分配缓冲 `copyto!`。*一致性*: 元素数值与求和顺序不变，只是不再重复计算。前置: 依赖事实 3a 结论——KT 若随 D/T 变化，则只缓存不变部分（如 MT），KT 保留重建但改为预分配 COO 缓冲。
- **A2 分解复用**: 系统矩阵 (M/dt + K) 每步不变 → 缓存一次分解因子，每步只 `ldiv!`。*一致性*: 同一因子 + 同一 rhs → 解 bit 一致。KT 随损伤变化时降级为 `lu!(F, A)` 符号结构复用（数值分解确定性，同输入同输出）。前置: 事实 3a。
- **A3 对流 BC 并入缓存**: `apply_cool_method` 对流项系数（h、边界几何）不变，一次性并入缓存矩阵。前置: A1。

### B. 逐单元 SPMe 循环（`src/CallModel.jl`、`src/SPMe.jl`）

- **B1 `ws_pool` 每步重建 → case 层缓存**: 按 objectid 失效，复用 `ensure_czm_cache` 既有模式，挂载点用 `MultiSPMeLayout`/`MeshGeometry` 所在的 `CouplingState.jl` 体系。
- **B2 `copy_element_results` 18 键 Dict → 预分配复用**: 对外接口不变。
- **B3 `case.index` keys 收集 + filter → 预计算一次缓存 var_list**。

### C. 全局拼装与求解（`src/CallModel.jl`、`src/Solve.jl`）

- **C1 两级 blockdiag → 预分配 CSC 缓冲 + `copyto!`**: *一致性*: 拼出的矩阵数值完全相同。前置: 无。
- **C2 求解复用**: 依赖 A2 同一结论（M、K 不变 → 全局系统矩阵不变）。

### D. CZM（`src/CzmSolve.jl`）——定位为 2026-04-20 向量化计划的延续

- **D1 `K_total = K_bulk + K_coh` 每迭代分配 → 预分配缓冲 + 同序逐元素加**。*前置小实验*: 随机矩阵上手写同序加法与 SparseArrays `+` 结果 `==` 比对（同 pattern 下数值路径一致才做，否则放弃）。
- **D2 线搜索内最多 8 次全量重组装 → 只更新随损伤变化部分**。深度视事实 3b 与 CZM 实测占比决定；**开批前必须先完成批次 0 第 4 项残留核对**，与既有计划重叠的部分并入既有计划或明确接手。
- **D3 `clone_czm_mesh_with_damage` 每步克隆 → 视测量决定**。

### E. 分流求解器

静态分析已预分配、较轻。除非实测占比 ≥10% 才考虑，否则不动。

## 6. 批次排序规则与验收流程

**排序**: 按批次 0 实测占比从高到低开批；同模块内按依赖序（A1→A2→A3、C1→C2、D1 在 D2 前）。每批可包含同一模块内依赖链上的多个字母项，但一个批次不跨模块。

**每批固定流程**:
1. 改动单独一个 commit；
2. 前后各跑一次本机基线，对比退出码、网格/步数、`metrics.toml` 记录精度、PNG SHA-256（若批次 0 验证 SHA 稳定）——**任一不一致即回退该批**，定位或放弃该手段；
3. 该批模块的 `timing_totals`（+ 墙钟）前后对比作为提速证据；
4. findings.md 记录本批结果与下一批决定。

## 7. 风险清单

| 风险 | 应对 |
|------|------|
| `k_eff(D)` 损伤-热阻耦合使 KT 每步变 | A2/C2 自动降级为符号分解复用；A1/A3 缓存不变部分不受影响 |
| 手写同序加法与 SparseArrays `+` 数值路径不一致 | D1 前置等价性小实验，不过则放弃 D1 |
| PNG SHA 本机不稳定 | 批次 0 重跑两次验证；不稳定则改三判据 |
| debug 计时自身开销扭曲占比 | Profile 火焰图交叉确认 |
| 分配消除改变 GC 时机影响 timing 口径 | 提速证据同时看墙钟与分段两项 |
| 与 2026-04-20 CZM 向量化计划冲突/重复 | 批次 0 第 4 项残留核对；D 组开批前完成 |
| `case.index` 等 `Dict{String,Any}` 深层类型不稳定 | 本计划只做 B3 预计算缓存，不做类型重构（归 CouplingState 后续工作） |

## 8. 成功标准

- 全部批次完成后: 基线场景（testexample.jl）墙钟时间相对批次 0 锚点**可测量下降**（具体目标值在批次 0 拿到占比数据后，按各批次理论收益上限折算并经用户确认后写入 findings.md）；
- 所有已实施批次保持 bit 级基线一致（判据见第 2 节）；
- findings.md 完整记录: 占比表、事实核查结论、每批前后 timing 对比、未做项及原因。
