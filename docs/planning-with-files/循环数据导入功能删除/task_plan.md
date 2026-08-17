# CycleData 外部循环数据导入求解功能删除计划

## Goal

删除 `src/CycleData.jl` 中读取外部循环数据并据此继续/驱动求解的相关功能，同时保留正常循环求解、仿真过程数据导出和 CSV 输出语义。

## Current Phase

Complete

## Phases

### Phase 1：边界与调用点审计
- [x] 对照源码函数索引与实际定义
- [x] 查清导出、测试、示例和下游调用
- [x] 明确保留与删除的 API 边界
- **Status:** complete

### Phase 2：实现删除
- [x] 删除外部循环数据导入/恢复求解函数
- [x] 删除对应顶层导出、示例或测试引用
- [x] 更新 `md/源码函数索引/CycleData.md`
- **Status:** complete

### Phase 3：验证与索引
- [x] 运行聚焦测试
- [x] 运行受影响路径的导出烟雾测试和冻结 `testexample`
- [x] 更新 `docs/planning-with-files/index.md`
- **Status:** complete

## Constraints

- 不删除正常 `solve_phase_with_export` / `solve_cycle_with_export` 数据采集与导出能力，除非调用审计证明其属于外部导入求解链。
- 不改变 `CycleSolver.jl` 的普通循环求解结果语义。
- 删除后不得残留无效 export、索引条目或调用点。

## Errors Encountered

| Timestamp | Error | Attempt | Resolution |
|---|---|---:|---|
