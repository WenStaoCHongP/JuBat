# ring.jl Dead Code Verification Plan

**Status:** ⬜ Pending | **Layer:** 2 几何/网格 | **桶:** Delete 候选

**Goal:** grep 验证 `ring_mesh` 调用者；无则删除整个文件。联动 `parameters/Ring.jl`、`ThermalPolar2D.jl`。

## Tasks

- [ ] **Step 1: grep 验证**

```bash
grep -rn "ring_mesh\b" --include="*.jl" --include="*.md" .
```

- [ ] **Step 2: 决策**

| 调用点情况 | 决策 |
|---|---|
| 仅 `ring.jl` 自身 + `JuBat.jl` export | 删 `ring.jl` + 移除 export + 联动删 `ThermalPolar2D.jl`（如仅被 ring 路径用） |
| example/ 有调用 | 保留，标 Leave alone |
| test/ 有调用 | 保留 |

- [ ] **Step 3: 联动决策**

参考 `parameters-Ring.md` 与 `ThermalPolar2D.md`。三者一同通过 PR。

- [ ] **Step 4: Commit + baseline**

```bash
git rm src/ring.jl  # 或保留
git commit -m "chore(ring): 决策=[删除/保留]"
```

## Risk
低；grep 是可靠判据。
