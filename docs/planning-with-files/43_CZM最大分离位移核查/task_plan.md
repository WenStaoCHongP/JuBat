# CZM 最大分离位移核查计划

## Goal

判断当前 `testexample.jl` 报告的最大法向分离位移是否由正确的节点位移、拓扑定向、高斯积分、时间历史归约和尺度还原得到，并给出可复现的数值证据与适用边界。

## Current Phase

已完成

## Phases

### Phase 1：数据流与现状复现

- [x] 追踪 `CZMResult.separation_n` 到最终输出的完整调用链
- [x] 核对归一化尺度与“最大值”的时间/空间归约语义
- [x] 记录当前工作树与近期边界变化
- **Status:** completed

### Phase 2：单元级与手算交叉验证

- [x] 运行已有 cohesive 分离/法向单元测试
- [x] 构造受控刚体、纯法向、纯切向和节点顺序测试
- [x] 用独立公式重算高斯点及单元平均分离
- **Status:** completed

### Phase 3：端到端定位与敏感性核查

- [x] 定位 `testexample` 最大值对应的时间层、cohesive 单元和节点
- [x] 对比结果字典、原始 `CZMResult` 与位移场独立重算值
- [x] 核对边界修改前后差异是否来自物理解而非后处理错误
- **Status:** completed

### Phase 4：结论与交付

- [x] 汇总正确性判据、发现的问题和证据
- [x] 记录测试结果与已知限制
- [x] 向用户交付分析结论；未经请求不修改生产代码
- **Status:** completed

## Key Questions

1. 输出的 `czm δ_max_n [m]` 是每步空间最大值历史，还是全时空单一标量？
2. `compute_separation` 的法向、节点顺序、积分与单位还原是否彼此一致？
3. 当前 60 s 快速基线的最大分离能否由位移场独立逐位复算？
4. `fix_inner=false` 新增 23 个固定点后，结果字典是否完整记录每次力学更新？

## Decisions Made

| Decision | Rationale |
|---|---|
| 先诊断不修复 | 用户要求运行单元级脚本并分析，尚未授权生产代码修复 |
| 同时检查单元级与端到端链路 | 单元公式正确不等于时间历史和有量纲输出正确 |
| 不用零值或近似掩盖异常 | 保持 JuBat 工程状态与科学输出的严格失败语义 |
| `fix_inner=false` 也叠加分层端点 | 用户更正边界契约：相对原外圈净增 23 点 |
| 不新写失败测试 | 用户显式要求复用现有 `true` 模式测试；额外以直接计数核对 `false` |

## Errors Encountered

| Error | Attempt | Resolution |
|---|---:|---|
| 暂无 | 1 | — |

## Notes

- 当前工作树包含用户已有的 `src/CouplingState.jl`、`src/CzmSolve.jl` 等未提交修改；调查不得覆盖。
- 当前 `testexample.jl` 为 `fix_inner=false`；修改后该模式为外圈 + 分层端点，相对原外圈净增 23 个固定节点。
- 修改前后最大分离由 `7.4202e-13 m` 变为 `1.0667e-12 m`，该变化本身不证明正确或错误。
