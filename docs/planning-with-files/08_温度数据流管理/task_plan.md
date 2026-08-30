# Task Plan: 温度模块数据流管理

## Goal
清理温度模块在 Variables.jl 中的变量定义，消除死内存分配、统一输出路径、修复已知 bug，使其对齐"Variables.jl 定义、PostProcessing.jl 还原"的设计规范。

## Current Phase
Complete

## Phases

### Phase 1: 审计与问题定位
- [x] 审计全部 thermal2D 变量键的定义、写入、还原三端
- [x] 识别死内存分配、缺失还原、键名 bug、输出路径不一致
- **Status:** complete

### Phase 2: 清理死内存分配
- [x] 2.1 删除 Variables.jl 中 9 个从未写入的变量键（节省 ~5-15 MB/仿真） — `6346299`
- [x] 2.2 修复 `"thermal2D temperature  at nodes history"` 双空格 bug — `6346299`
- **Status:** complete

### Phase 3: 统一输出路径
- [x] 3.1 将 Solve.jl:392-401 中的热变量直接输出迁移到 PostProcessing.jl — `dffc0cd`
- [x] 3.2 补齐 PostProcessing.jl 中缺失的有价值输出（element current, SOC, voltages 等） — `dffc0cd`
- [x] 3.3 统一 PostProcessing 的输出键命名风格（带物理单位） — `dffc0cd`
- **Status:** complete

### Phase 4: 验证
- [x] 模块加载验证通过
- [x] 导出符号验证通过
- **Status:** complete

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 删除死分配而非修复写入 | 这些变量从未被任何代码路径使用，是历史残留 |
| 将 Solve.jl 热输出迁移到 PostProcessing.jl | 统一还原入口，与电化学和 CZM 模式一致 |
| 保留中间变量在 StandardVariables 中 | soc_n/soc_p 等虽用于内部耦合，但历史数据对分析有价值 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None | - | - |
