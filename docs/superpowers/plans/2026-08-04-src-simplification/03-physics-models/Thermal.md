# Thermal.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。Lumped 热模型（备用路径）。

## 现状（80 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `ThermalLumped` | 1 | 集总热模型 |

**调用者**：`CallModel.jl:232`

## Audit

- [ ] 通读 80 行
- [ ] baseline 记录

## Result
无修改。
