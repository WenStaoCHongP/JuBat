# CZM 分离输出修复计划

## Goal

保留现有法向分离结果契约，新增有效分离历史；`mix` 模式打印全时空最大有效分离，`model1` 打印全时空最大法向分离，并消除 CZM 历史的一步滞后与逐单元分离未写回问题。

## Current Phase

已完成

## Phases

### Phase 1：输出契约 RED 测试

- [x] 建立不复用生产归约逻辑的受控 CZM 结果契约
- [x] 证明当前缺少 `czm δ_max_eff [m]`
- [x] 证明求解后历史必须采用本步收敛结果
- **Status:** completed

### Phase 2：最小生产修复

- [x] 初始化、归一化还原并保留有效分离历史
- [x] 在收敛的 CZM 更新后、历史写入前刷新当前场与峰值摘要
- [x] 按 `czm.model` 选择示例打印量和标签
- **Status:** completed

### Phase 3：验证与交付

- [x] RED 测试转 GREEN
- [x] 运行受影响 CZM 测试与 60 s `testexample.jl`
- [x] 用快照/节点复算确认打印值不滞后
- **Status:** completed

## Decisions Made

| Decision | Rationale |
|---|---|
| 保留 `czm δ_max_n [m]` | 不破坏现有结果键和网格敏感性脚本 |
| 新增 `czm δ_max_eff [m]` | `mix` 损伤驱动量是有效分离，不应冒充法向分离 |
| 收敛后刷新所有 CZM 输出 | 解决峰值滞后和逐单元历史全零的共同根因 |
| 不删除原有输出 | 遵守科学/API 契约保留原则 |

## Errors Encountered

| Error | Attempt | Resolution |
|---|---:|---|
| 暂无 | 1 | — |
