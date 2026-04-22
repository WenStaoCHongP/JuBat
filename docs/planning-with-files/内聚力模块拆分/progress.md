# Progress Log

## Session: 2026-04-22

### Phase 1: 责任边界盘点
- **Status:** complete
- **Started:** 2026-04-22
- Actions taken:
  - 读取 planning-with-files 相关模板，确认计划文件结构。
  - 读取 `src/CzmSolve.jl`、`src/czm.jl` 和 `src/JuBat.jl`，确认当前职责分布。
  - 参考仓库记忆，确认归一化、后处理和收敛语义上的历史约束。
  - 量化基线统计：CzmSolve.jl 1001 行/19 函数，czm.jl 659 行/11 函数。
  - 重复代码统计：~314 行重复（BC 提取 60 行、clone 44 行、Newton 循环 180 行、末尾重组 30 行）。
  - 产出完整函数清单表（含行数、职责、调用者、嵌套深度）。
  - 产出调用依赖关系图和分层架构图。
- Files created/modified:
  - `docs/planning-with-files/内聚力模块拆分/task_plan.md` (created)
  - `docs/planning-with-files/内聚力模块拆分/findings.md` (created → 大幅更新)
  - `docs/planning-with-files/内聚力模块拆分/progress.md` (created)

### Phase 2: 模块拆分设计
- **Status:** complete
- **Started:** 2026-04-22
- Actions taken:
  - 回答了 3 个 Key Questions：
    - Q1: 保留 CzmSolve.jl 作为调度入口，拆出 CzmPostProcess.jl，并把耦合 helpers 并入现有 CouplingState.jl
    - Q2: 列出所有必须保持原名的符号（已导出 + 外部脚本调用）
    - Q3: 先拆后处理（最低风险）→ 参数计算 → helper 提取 → 求解器
  - 设计了目标文件布局和 include 顺序
  - 制定了 Protected Constraints（5 条不可破坏的语义约束）
  - 制定了回滚策略（feature branch + per-phase commit + revert）
  - 创建了 CZM 行为基线探针脚本 `tools/czm_baseline_probe.jl`
  - 根据 reviewer 反馈收窄 Task 6：`backtrack_line_search!` 只服务 `solve_czm_basic_step`，`newton_raphson_czm` 因惩罚式 BC 残差语义差异而排除在外
- Files created/modified:
  - `docs/planning-with-files/内聚力模块拆分/task_plan.md` (大幅更新)
  - `docs/planning-with-files/内聚力模块拆分/findings.md` (更新依赖图)
  - `docs/planning-with-files/内聚力模块拆分/progress.md` (更新)
  - `tools/czm_baseline_probe.jl` (新建)

### Phase 3-5: 待执行（helper 提取 + CouplingState.jl 扩展）
- **Status:** pending
- 待用户确认基线并批准后开始实施

## Test Results
| Test | Input | Expected | Actual | Status |
|------|-------|----------|--------|--------|
| CZM Baseline Probe | nθ=40, dT=5K/T0, Δsoc=0.1 | 待运行记录 | 待运行 | pending |

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
| None | - | - | - |

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Phase 2 已完成，Phase 3 待启动 |
| Where am I going? | Phase 3 代码结构简化（helper 提取 + 文件拆分） |
| What's the goal? | 把内聚力模块拆成边界清晰、易维护的结构 |
| What have I learned? | 314 行重复代码，19+11 函数，5 条不可破坏约束 |
| What have I done? | Phase 1+2 完成，3 个 Key Questions 闭合，基线探针就绪 |
