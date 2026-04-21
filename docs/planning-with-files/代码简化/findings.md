# 代码简化研究记录

## Requirements
- 制定长期的 src 代码简化计划。
- 目标是减少行数，提高可读性，并尽量降低回归风险。
- 先完成 src 臃肿代码识别和可简化性排序。
- 计划文件放在 docs/planning-with-files/代码简化/。

## Research Findings

### 体量最大的 src 文件
| 文件 | 行数 | 初步判断 |
|------|------|----------|
| [src/SetMesh.jl](src/SetMesh.jl) | 703 | 体量最大，重复网格操作多，但几何正确性风险较高 |
| [src/CycleData.jl](src/CycleData.jl) | 605 | 导出与格式化密集，适合第一波清理 |
| [src/CzmSolve.jl](src/CzmSolve.jl) | 546 | 求解框架冗长，Newton 循环重复明显 |
| [src/Jellyrollmodel.jl](src/Jellyrollmodel.jl) | 492 | 几何与拓扑逻辑耦合较深 |
| [src/CycleSolver.jl](src/CycleSolver.jl) | 484 | 状态流控制较重，属于高耦合核心层 |
| [src/czm.jl](src/czm.jl) | 459 | CZM 核心逻辑加网格创建，重构需谨慎 |
| [src/Solve.jl](src/Solve.jl) | 449 | 主求解器编排，属于系统主入口 |
| [src/ThermalDistributed.jl](src/ThermalDistributed.jl) | 434 | 热源与边界处理重复明显 |
| [src/SetParams.jl](src/SetParams.jl) | 418 | 核心参数与归一化入口，适合局部抽取 |
| [src/Parallelsolution.jl](src/Parallelsolution.jl) | 368 | 分流牛顿逻辑可抽象 |

### 结构性候选文件
| 文件 | 行数 | 说明 |
|------|------|------|
| [src/Variables.jl](src/Variables.jl) | 255 | 体量不大但 Dict-heavy，适合结构重整 |
| [src/CallModel.jl](src/CallModel.jl) | 219 | 体量不大但分发逻辑集中 |
| [src/SPMe.jl](src/SPMe.jl) | 271 | 核心电化学模块，后续通过接口降低重复 |
| [src/P2D.jl](src/P2D.jl) | 271 | 核心电化学模块，后续通过接口降低重复 |

### 可简化性排序
1. [src/CycleData.jl](src/CycleData.jl) - 导出和格式化重复多，低风险，适合先收缩。
2. [src/PostProcessing.jl](src/PostProcessing.jl) - 单位还原和变量映射重复明显，适合抽表驱动逻辑。
3. [src/Materialmatrix.jl](src/Materialmatrix.jl) - 插值和材料矩阵装配较规整，可做局部抽象。
4. [src/SetMesh.jl](src/SetMesh.jl) - 体量最大，但要兼顾网格正确性，适合在上面三项之后推进。
5. [src/Parallelsolution.jl](src/Parallelsolution.jl) - 牛顿求解和截止检测混在一起，可提取通用求解器骨架。
6. [src/ThermalDistributed.jl](src/ThermalDistributed.jl) - 热源和边界处理可继续合并，但需要保留数值行为。
7. [src/CzmSolve.jl](src/CzmSolve.jl) - 重复迭代骨架较多，收益高，但属于核心求解路径。
8. [src/Jellyrollmodel.jl](src/Jellyrollmodel.jl) - 几何与拓扑耦合深，适合后置重构。
9. [src/Variables.jl](src/Variables.jl) - 体量不大但中央数据结构很重，适合做结构层简化。
10. [src/CallModel.jl](src/CallModel.jl) - 分发入口集中，适合统一接口收口。
11. [src/SetParams.jl](src/SetParams.jl) - 归一化和核心参数层，宜做局部抽取，不宜大拆。
12. [src/Solve.jl](src/Solve.jl) - 主求解器编排，适合在前几波稳定后再拆分。
13. [src/CycleSolver.jl](src/CycleSolver.jl) - 状态机复杂，回归风险高。
14. [src/czm.jl](src/czm.jl) - CZM 和网格创建耦合，重构需很谨慎。
15. [src/SPMe.jl](src/SPMe.jl) / [src/P2D.jl](src/P2D.jl) - 核心电化学模块，需在接口层稳定后再处理。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 先做低耦合重复逻辑，再做核心求解器 | 能更快减少行数，同时降低回归风险 |
| 用“体量 + 可简化性 + 风险”联合排序 | 只看行数会把高风险核心文件排得过早 |
| 把核心模块分到后续波次 | Solve、CycleSolver、czm、SPMe、P2D 改动面大，适合等基础层稳定后再动 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| PowerShell 在统计行数时出现空管道元素语法错误 | 改成数组 + ForEach-Object 的写法后重新执行 |

## Resources
- 统计来源：src 目录 Julia 文件行数盘点
- 只读探查来源：src 简化优先级探索结果

## Visual/Browser Findings
- 无