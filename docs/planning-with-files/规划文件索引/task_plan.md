# Task Plan: planning-with-files 总索引

## Goal
为 `docs/planning-with-files/` 建立可维护的总索引，说明每个任务目录和文件的用途、建立时间、最近修改时间、Git 修改次数、状态及跟踪口径。

## Current Phase
Phase 3

## Phases

### Phase 1: 目录与元数据盘点
- [x] 盘点所有任务目录和文件
- [x] 提取任务目标与状态
- [x] 获取 Git 首次/最近提交时间和提交次数
- **Status:** complete

### Phase 2: 索引设计与生成
- [x] 定义目录汇总表与逐文件明细表
- [x] 对未跟踪文件采用文件系统时间回退
- [x] 创建 `docs/planning-with-files/index.md`
- **Status:** complete

### Phase 3: 校验与交付
- [x] 核对索引覆盖所有现有文件
- [x] 核对链接、时间和计数口径
- [x] 更新项目约定并交付
- **Status:** complete

## Key Questions
1. “建立时间”和“修改次数”应如何保证可复现？
2. 如何同时覆盖标准三文件任务目录与含探针/独立设计文档的目录？
3. 索引如何避免在新增任务后迅速过时？

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 索引放在 `docs/planning-with-files/index.md` | 作为该目录统一入口，路径直观 |
| 优先使用 Git 历史，未跟踪文件回退文件系统元数据 | Git 数据可复现，同时不遗漏新增文件 |
| 同时提供任务目录汇总和逐文件明细 | 兼顾快速浏览与完整审计 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
