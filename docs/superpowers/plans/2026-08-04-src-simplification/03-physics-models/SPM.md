# SPM.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。Option 默认模型（`opt.model="SPM"`），活跃 API。

## 现状（85 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `SPM` | 1 | 入口（颗粒模型，无电解液） |

**调用者**：`CallModel.jl:219`、`Initialisation.jl:3`、`example/change_time_step.jl:10`

## Audit

- [ ] 通读 85 行，确认无重复（与 SPMe.jl 的差异是有意设计）
- [ ] baseline 记录

## Result
无修改。
