# JuBat.jl Module Entry Audit Plan

**Status:** ⬜ Pending | **Layer:** 5 后处理/导出 | **桶:** Leave alone + export 清理

**Goal:** 审查 36 个 include + ~50 个 export。export 清理在相关文件被删/合并时同步更新。

## 现状（89 行）

| 区块 | 行号 | 用途 |
|---|---|---|
| includes | 4-37 | 36 个 src 文件装载 |
| exports | 40-89 | 公开 API |

## Tasks

### Task 1: include 顺序审查

**Files:** 仅审查 `src/JuBat.jl:4-37`

- [ ] 确认依赖顺序合理（如 `CouplingState.jl` 先于 `SetCase.jl`，`czm.jl` 先于 `CzmSolve.jl`）
- [ ] 确认无循环 include

### Task 2: export 清理（联动其他 plan）

**Files:** Modify `src/JuBat.jl:40-89`

各 plan 完成后同步：

- [ ] `czm.md` Task 1 → 移除 `assemble_coupled_system_full`（line 69）
- [ ] `ring.md` 若删 ring.jl → 移除 `ring_mesh`（line 51）
- [ ] parameters/Enertech.md 等若删参数文件 → 检查 ChooseCell 字符串分支

### Task 3: include 顺序优化（可选，不强制）

**Files:** 仅审查

- [ ] 若发现 include 顺序混乱（如 `Variables.jl` 在模型之后），记录但不强制重排（依赖关系可能脆弱）

### Task 4: baseline

```bash
echo "$(date +%F): JuBat.jl 已审查，桶=Leave alone（export 联动其他 plan）" >> Simplify/baseline.md
```

## Result
本文件自身基本无修改；export 清理随相关 plan 同步进行。
