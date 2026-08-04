# SetMesh.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 2 几何/网格 | **桶:** High-risk-leave-alone

**Goal:** 仅审查，不动。基础设施影响所有模型。

## Audit
- [ ] 通读 767 行，确认无 try/catch 兜底、无向后兼容入口
- [ ] grep `Assemble1D`（如有计划清查）— 实际归属 Assemble.jl，本文件无关
- [ ] baseline 记录：`echo "$(date +%F): SetMesh.jl 已审查，桶=High-risk-leave-alone" >> Simplify/baseline.md`

## Result
无修改。数值基础设施不动。
