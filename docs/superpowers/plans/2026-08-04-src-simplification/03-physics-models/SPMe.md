# SPMe.jl Audit Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Leave alone

**Goal:** 仅审查。主线核心电化学模型，D-cluster 未命中。

## 现状（290 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `SPMe` | 1 | 主入口 |
| `SPMe_element` | — | per-element 子模型（多 SPMe 架构） |

**调用者**：`CallModel.jl:222`（`opt.model=="SPMe"` 主分支）

## Audit

- [ ] 通读 290 行，确认无 try/catch、无重复矩阵组装
- [ ] grep `ElectrolyteDiffusion` 调用点（line 28, 79），确认参数传递一致
- [ ] baseline 记录：

```bash
echo "$(date +%F): SPMe.jl 已审查，桶=Leave alone" >> Simplify/baseline.md
```

## Result
无修改。
