# Initialisation.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 4 求解器 | **桶:** Leave alone

**Goal:** 仅审查。模型初始化分支（SPM/SPMe/P2D/y0 注入）。

## 现状（152 行）

| 函数 | 用途 |
|---|---|
| `ModelInitialisation` | 按 `opt.model` 分支初始化 y0 |

**调用者**：`Solve.jl`、`SetCase.jl`

## Audit

- [ ] 通读 152 行
- [ ] 确认 3 个 model 分支（SPM/SPMe/P2D）无重复（共同部分已抽公共）
- [ ] baseline 记录

## Result
无修改。
