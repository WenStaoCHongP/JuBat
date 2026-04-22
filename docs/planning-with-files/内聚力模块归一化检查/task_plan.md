# Task Plan: 内聚力模块参数归一化检查

## Goal
核对 JuBat 内聚力模块的参数归一化链路，确认强度、位移、刚度和尺度量在 NormaliseParam 中的定义、传递与使用一致。

## Current Phase
Complete

## Phases

### Phase 1: Requirements & Discovery
- [x] 明确要检查的内聚力参数范围
- [x] 找到归一化入口与相关尺度字段
- [x] 记录初步发现到 findings.md
- **Status:** complete

### Phase 2: Code Path Tracing
- [x] 追踪 cohesive 参数从原始参数到归一化参数的赋值链路
- [x] 检查 CZM、本构和后处理中的使用位置
- [x] 标记任何依赖 scale 的计算
- **Status:** complete

### Phase 3: Targeted Verification
- [x] 运行最小化 Julia 检查脚本或现有例子
- [x] 对比归一化前后关键字段
- [x] 记录是否存在量纲不一致或重复缩放
- **Status:** complete

### Phase 4: Findings Summary
- [x] 汇总结论与风险点
- [x] 给出是否需要代码修正的判断
- [x] 更新 progress.md
- **Status:** complete

## Key Questions
1. cohesive 的强度、临界位移和刚度是否都在同一归一化体系下处理？
2. scale.σ_czm、scale.δ_czm、scale.K_czm 的来源是否与文档和代码一致？
3. CZM 相关后续计算是否再次缩放，导致重复归一化或漏归一化？

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 以 NormaliseParam 和 CZM 相关调用链为主线 | 这是最直接的归一化入口，能最快验证尺度是否一致 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None yet | 1 | N/A |

## Notes
- 先确认参数定义，再确认使用位置，最后做运行时验证。
- 如果发现问题，优先定位到归一化入口而不是后处理输出。