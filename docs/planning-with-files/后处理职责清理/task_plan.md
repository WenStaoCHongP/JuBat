# 后处理职责清理计划

## Goal

重命名四个后处理文件中以下划线开头的函数，消除职责交叉和重复实现，并将 `PostProcessing.jl` 第 117 行后的函数迁移到更合适的现有文件，使 `PostProcessing.jl` 仅保留通用结果还原职责且数值/API 语义清晰。

## Current Phase

Complete

## Phases

### Phase 1：定义、调用与职责审计
- [x] 核对四份源码索引与实际定义
- [x] 收集所有下划线前缀函数及调用者
- [x] 分析第 117 行后函数的职责和重复实现
- **Status:** complete

### Phase 2：迁移与重命名设计
- [x] 确定每个函数的新名称和目标文件
- [x] 明确公开 API、include 顺序和兼容边界
- [x] 记录不迁移/不合并项及理由
- **Status:** complete

### Phase 3：实施源码重构
- [x] 移动 `PostProcessing.jl` 第 117 行后的函数
- [x] 重命名下划线前缀函数并更新全部调用点
- [x] 删除确认重复的实现或合并共享逻辑
- **Status:** complete

### Phase 4：索引与文档同步
- [x] 更新四份源码函数索引及总索引
- [x] 更新受影响技术文档和规划总索引
- [x] 新增职责/命名架构测试
- **Status:** complete

### Phase 5：验证
- [x] 运行后处理、CSV、循环和 CZM 聚焦测试
- [x] 运行冻结的 `example/testexample.jl`
- [x] 检查活动源码/索引无旧名称和职责残留
- **Status:** complete

## Key Questions

1. “下划线前缀函数”范围是否完全落在指定四个后处理文件及其调用链？
2. `PostProcessing.jl` 第 117 行后每个函数属于 CSV、循环数据还是 CZM 后处理？
3. 重复功能应合并到哪个单一权威实现？

## Decisions Made

| Decision | Rationale |
|---|---|
| 保持结果字典键、CSV 列名和数值语义不变 | 本任务是职责与命名清理，不应改变外部结果契约 |
| 不保留下划线旧名 façade | 用户明确要求重命名，应同步全部调用者并清除旧定义 |
| 规划文件存放在本任务独立目录 | 遵守项目 planning-with-files 约定 |

## Errors Encountered

| Error | Attempt | Resolution |
|---|---:|---|
| 暂无 | — | — |

## Constraints

- 不回退或覆盖当前脏工作树中的其他修改。
- 不改变后处理结果键、单位、CSV 列名和冻结科学结果。
- `PostProcessing.jl` 最终只保留通用结果后处理，不放循环/CZM/CSV 专属逻辑。
- 修改后必须通过冻结 `testexample` 数值和 PNG 哈希门禁。
