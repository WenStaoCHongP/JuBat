# src 代码简化长期计划

## Goal
建立一套长期、分波次的 src 简化路线图，优先清理重复、易抽取的代码，持续缩减行数并提升可读性，同时避免破坏电化学、热学和 CZM 的核心行为。

## Current Phase
Phase 1: 基线盘点与排序 complete

## Phases

### Phase 1: 基线盘点与排序
- [x] 统计 src 主要文件体量
- [x] 按可简化性和回归风险划分优先级
- [x] 记录第一批候选文件和理由
- **Status:** complete

### Phase 2: 第一波快速收缩
- [ ] 优先处理 PostProcessing.jl、CycleData.jl、Materialmatrix.jl
- [ ] 抽取重复的单位还原、导出格式化和插值辅助函数
- [ ] 用小范围回归验证保持输出一致
- **Status:** pending

### Phase 3: 中等风险高收益重构
- [ ] 拆分 SetMesh.jl、ThermalDistributed.jl、Parallelsolution.jl
- [ ] 收敛重复的网格、边界和求解器控制逻辑
- [ ] 降低工具函数散落和长条件分支
- **Status:** pending

### Phase 4: 核心流程重整
- [ ] 重构 CzmSolve.jl、Jellyrollmodel.jl、Solve.jl、CycleSolver.jl
- [ ] 把 Newton 循环、状态流和几何/网格职责拆开
- [ ] 为关键路径补回归基线
- **Status:** pending

### Phase 5: 结构层收口
- [ ] 评估 Variables.jl、CallModel.jl、SetParams.jl、czm.jl、SPMe.jl、P2D.jl
- [ ] 优先做接口收口和数据结构简化，再考虑大改
- [ ] 只在前几波稳定后推进
- **Status:** pending

## Key Questions
1. 哪些文件在不改行为的前提下最容易减少 LOC？
2. 哪些核心文件需要先补回归测试再动？
3. 是否需要新增统一的 helper/interface 层来吸收重复逻辑？

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 采用“可简化性 + 回归风险”排序，而不是只看行数 | 避免把体量最大的文件误判为第一优先级 |
| 第一波只做局部、重复、低耦合代码 | 先拿到可见收益，再推进核心模块 |
| 核心求解器、状态机和几何网格放到后续波次 | 这些区域收益高，但回归风险也最高 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| PowerShell 空管道元素语法错误 | 1 | 改成数组 + ForEach-Object 的写法后重新统计行数 |

## Notes
- 先把简化顺序定下来，再按波次推进实现。
- 每一波都要同步更新 findings.md 和 progress.md。