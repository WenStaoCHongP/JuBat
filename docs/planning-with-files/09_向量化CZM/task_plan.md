# Task Plan: CZM 向量化优化分析

## Goal
分析 CZM 求解器的性能瓶颈，判断向量化优化的可行方向，并形成一份可执行的 superpowers 计划文档。

## Current Phase
Complete

## Phases

### Phase 1: Requirements & Discovery
- [ ] 明确用户要解决的问题范围
- [ ] 梳理现有 CZM 求解流程与性能瓶颈
- [ ] 收集与向量化相关的代码和文档线索
- **Status:** complete

### Phase 2: Direction Analysis
- [ ] 归纳可行的优化方向
- [ ] 识别收益、风险和依赖条件
- [ ] 形成优先级排序
- **Status:** complete

### Phase 3: Superpowers Plan Drafting
- [ ] 产出 superpowers 风格计划文档
- [ ] 明确分阶段实施路径与验证项
- [ ] 对齐现有文档结构和命名
- **Status:** complete

### Phase 4: Review & Delivery
- [ ] 检查计划是否自洽
- [ ] 记录残余风险与后续建议
- [ ] 向用户交付结论
- **Status:** complete

## Key Questions
1. CZM 的主要耗时来自哪里：Newton 迭代、子步循环，还是单元级装配？
2. 向量化应优先作用于哪些数据结构：单元循环、载荷子步，还是局部线性代数？
3. 需要把结果落到新的 superpowers 计划文件，还是仅在现有计划上补充？

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 使用 superpowers 计划格式输出 | 与仓库现有优化文档一致，便于后续执行 |
| 先分析方向再落计划 | 避免在未确认瓶颈结构前过早固化实现方案 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None | 1 | N/A |

## Notes
- 优先对齐 docs/superpowers/plans 的既有写法。
- 需要把发现及时记录到 findings.md，避免上下文丢失。