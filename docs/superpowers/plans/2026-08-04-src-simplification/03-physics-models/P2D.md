# P2D.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。Pseudo-two-dimensional 模型，主线 API。

## 现状（287 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `P2D` | 1 | 入口 |

**调用者**：`CallModel.jl:225`、`Initialisation.jl:17,41`、`example/minimal_example.jl:9`、`example/change_current.jl:13`

## Audit

- [ ] 通读 287 行
- [ ] 确认 `ElectrodeDiffusion`/`ElectrodePotential`/`ElectrolyteDiffusion`/`ElectrolytePotential` 调用参数顺序
- [ ] baseline 记录

## Result
无修改。
