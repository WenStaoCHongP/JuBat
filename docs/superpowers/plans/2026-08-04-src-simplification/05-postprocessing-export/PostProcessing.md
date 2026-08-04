# PostProcessing.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 5 后处理/导出 | **桶:** Leave alone

**Goal:** 仅审查。后处理主入口（按 model 分支提取结果）。

## 现状（350 行）

| 函数 | 用途 |
|---|---|
| `Postprocessing` | 主入口（SPM/SPMe/P2D 分支） |

**调用者**：export `JuBat.jl:40`

## Audit

- [ ] 通读 350 行
- [ ] 确认 model 分支仅做结果提取，无重复计算
- [ ] baseline 记录

## Result
无修改。
