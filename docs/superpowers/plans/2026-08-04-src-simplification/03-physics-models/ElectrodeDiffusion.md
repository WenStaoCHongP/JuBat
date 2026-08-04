# ElectrodeDiffusion.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。FEM 矩阵组装单元（M du/dt = Ku + F），D-cluster 未命中。

## 现状（18 行，最小单元）

**调用者**：`SPM.jl:22-23`、`SPMe.jl`、`P2D.jl:29-30`、export `JuBat.jl:40`

## Audit

- [ ] 通读 18 行
- [ ] baseline 记录

## Result
无修改。
