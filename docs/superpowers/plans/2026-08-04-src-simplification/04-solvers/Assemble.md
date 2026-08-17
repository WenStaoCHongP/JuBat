# Assemble.jl Audit Plan

**Status:** ✅ Completed（审计保留） | **Layer:** 4 求解器 | **桶:** High-risk-leave-alone

**Goal:** 仅审查。FEM 组装基础工具（`Assemble1D`/`Assemble`），所有物理模型依赖。

## 现状（40 行）

| 函数 | 用途 |
|---|---|
| `Assemble` | 稀疏矩阵组装（i,j,coeff 三元组） |
| `Assemble1D` | 1D 简化版 |

**调用者**：`ThermalDistributed.jl`、`ElectrodeDiffusion.jl` 等

## Audit

- [ ] 通读 40 行
- [ ] baseline 记录

## Result
无修改。数值基础设施不动。
