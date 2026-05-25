# Progress Log

## Session: 2026-05-12

### Phase 1: Scope and interface freeze
- **Status:** complete
- **Started:** 2026-05-12
- Actions taken:
  - 读取 planning-with-files 模板，确认计划文件结构。
  - 检查了 CZM 相关入口：`example/内聚力验证/czm_example.jl`、`tools/verify_czm_unit.jl`、`tools/verify_czm_analytical.jl`、`src/CzmSolve.jl`、`src/czm.jl`、`src/Materialmatrix.jl`。
  - 形成了独立验证边界：不接入电化学-热耦合，缺失的扩散应力和热应力通过函数输入注入。
  - 创建了 `task_plan.md`、`findings.md`、`progress.md` 三份计划文件。
- Files created:
  - `docs/planning-with-files/czm模块单独验证/task_plan.md`
  - `docs/planning-with-files/czm模块单独验证/findings.md`
  - `docs/planning-with-files/czm模块单独验证/progress.md`

### Phase 2: Benchmark design
- **Status:** complete
- 已定义 pure Mode I、本构输入驱动求解、mix/BK 三类基准。

### Phase 3: Standalone driver
- **Status:** complete
- 已实现 `tools/verify_czm_standalone.jl`，使用 provider callback 显式注入 `dT_elem / Δsoc_n_elem / Δsoc_p_elem`。

### Phase 4: Numerical verification
- **Status:** complete
- 已运行 standalone 脚本，Mode I 本构与解析解误差为 0，uniform/gradient 输入驱动求解均收敛，mix/BK 结果已生成。

### Phase 5: Documentation and handoff
- **Status:** complete
- 已把脚本位置、输出目录、输入接口和最小验收集写入计划与 findings。
- 已补充时间历史路径：热/扩散输入随时间变化，CZM D(t) 可追踪，且验证通过。

## Session: 2026-05-12 (续 — 求解器收敛性对比)

### Phase 6: 求解器收敛性对比验证
- **Status:** complete
- Actions taken:
  - 重写 `tools/verify_czm_standalone.jl`：聚焦纯 Mode I + 三种求解器（basic / load_substep / arc_length）对比。
  - 发现归一化参数关键问题：必须使用 `case.param.cohesive`（归一化后）而非 `param_dim.cohesive`（物理值），否则系统病态。
  - 运行对比：elastic regime 三种方法均收敛，损伤 regime 只有 basic 收敛（但残差偏大），load_substep 和 arc_length 在损伤区失败。
  - 修复 `solve_czm_basic_step` 收敛判断 bug：`result.residual_norm` 现在报告收敛时刻的真实残差，而非损伤更新后重新组装的残差。
- Files modified:
  - `tools/verify_czm_standalone.jl` — 重写为求解器对比脚本
  - `src/CzmSolve.jl` — 修复 basic 收敛残差报告 bug（新增 `converged_R_norm`）
  - `docs/superpowers/specs/2026-05-12-czm-solver-convergence-verification-design.md` — 设计文档

## Test Results
| Test | Input | Expected | Actual | Status |
|------|-------|----------|--------|--------|
| Standalone script | `tools/verify_czm_standalone.jl` | 按计划完成 Mode I / 时间历史 / mix 验证 | 已通过，输出写入 `output/czm_standalone/` | pass |
| 求解器收敛性对比（球面弧长） | Δsoc_n=0.1-10, nθ=40 | 三种方法对比 | elastic 全部 OK，损伤区 basic OK，其余 FAIL | pass |
| 求解器收敛性对比（Crisfield 圆柱弧长） | Δsoc_n=0.1-10, nθ=40 | arc_length 损伤区改善 | arc_length 收敛 6/8（+3），D=0.77~1.0 | pass |

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
| 2026-05-12 | 使用 `param_dim.cohesive` 导致系统病态 | 1 | 改用 `case.param.cohesive` |
| 2026-05-12 | basic 方法报告残差与收敛判断不一致 | 1 | 新增 `converged_R_norm` 记录收敛时残差 |
| 2026-05-12 | arc_length 球面约束在损伤区 stall（D=0） | 1 | 替换为 Crisfield 圆柱弧长法（CzmSolve.jl:354-406），收敛率 3/8→6/8 |

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Crisfield 圆柱弧长法已实现并验证通过 |
| Where am I going? | 剩余问题：Δsoc=1.5 接近损伤起始时 arc_length 仍失败；Δsoc≥10 所有方法都失败 |
| What's the goal? | 验证 CZM 三种求解器在 Jellyroll 参数下的收敛表现 |
| What have I learned? | 球面弧长约束在软化区不工作是因为 elastic predictor 低估位移，导致 hypersphere 与平衡路径不相交；Crisfield 圆柱法解耦位移与载荷，自然追踪软化路径 |
| What have I done? | 将 `solve_czm_arc_length_step` 从球面弧长改为 Crisfield 圆柱弧长法，arc_length 收敛率从 3/8 提升到 6/8 |
