# Findings & Decisions

## Requirements
- 分析文档中提到的 CZM 性能瓶颈。
- 重点判断“向量化 CZM 求解器”能否作为主要缓解手段。
- 输出一份符合 superpowers 风格的优化计划。

## Research Findings
- 文档明确指出 `solve_czm_step` 每步约 12s，占仿真时间的绝大部分，瓶颈在 Newton-Raphson + 20 载荷子步 × 801 个 CZM 单元。
- 现有结论认为这类瓶颈可通过增大 `czm_update_interval`、减少 `czm_load_steps` 或向量化 CZM 求解器缓解。
- `Solve` 主循环已把 CZM 损伤更新并入每步流程，因此优化计划需要同时考虑主循环调用频率和单步求解代价。
- 现有 superpowers 计划文件采用“Goal / Architecture / File Structure / Chunk / Task”结构，适合直接套用。
- `newton_raphson_czm` / `assemble_coupled_system` 才是实际热区，热点集中在每个 cohesive 单元的循环、矩阵组装和 Gauss 积分，而不是单独一个函数调用。
- `assemble_coupled_system` 每次都会重建 `K_bulk`，但 bulk 刚度对一次 CZM 更新来说通常是常量，存在直接缓存的空间。
- `compute_czm_strain_inputs` 仍在逐单元遍历节点求平均，属于容易做预分配和批量化处理的输入准备步骤。
- `CycleSolver.solve_phase` 在 `Solve` 已更新损伤后还会再次更新一次，计划里需要明确删除或重构这条冗余路径。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 采用“先分析、后计划”的两阶段方式 | 先确认瓶颈来源，再把优化拆成可验证步骤 |
| 将向量化作为核心方向之一 | 与文档中的性能建议一致，且对大量 CZM 单元更可能产生收益 |
| 计划文档放在 docs/superpowers/plans | 与仓库现有 superpowers 体系一致 |
| 把“向量化”和“缓存”一起列为方案 | 仅向量化单元循环不一定够，bulk stiffness 和几何不变量也应一起收敛 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| None | N/A |

## Resources
- `docs/代码优化/15_CZM耦合修复与单步损伤演化.md`
- `docs/superpowers/plans/2026-04-07-simulation-speedup.md`
- `docs/superpowers/plans/2026-04-02-initialisation-optimization.md`

## Visual/Browser Findings
- None yet.
