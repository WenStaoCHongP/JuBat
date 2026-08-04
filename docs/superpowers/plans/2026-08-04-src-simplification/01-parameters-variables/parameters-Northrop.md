# parameters/Northrop.jl Dead Code Verification Plan

**Status:** ⬜ Pending
**Layer:** 1 - 参数与变量
**桶:** Delete 候选

**Goal:** grep 验证 `ChooseCell("Northrop")` 是否有调用者；无则删除整个文件。

---

## Tasks

- [ ] **Step 1: grep 验证**

```bash
grep -rn 'ChooseCell("Northrop"\|"Northrop"' --include="*.jl" --include="*.md" --include="*.toml" .
```

- [ ] **Step 2: 决策与执行（同 Enertech 模板）**

- 若仅 `SetParams.jl:ChooseCell` 自身匹配 → 删整文件 + 移除 ChooseCell 分支 + 移除 JuBat.jl export（如有）
- 若有外部调用者 → 改桶为 Leave alone

- [ ] **Step 3: 跑 `example/minimal_example.jl` 验证**

- [ ] **Step 4: Commit + baseline 更新**

```bash
git add -A
git commit -m "chore(parameters): 删除未使用的 Northrop 参数集（dead code）"
echo "$(date +%F): parameters/Northrop.jl 决策=[删除/保留]" >> Simplify/baseline.md
```

---

## Risk

低；先 grep 再删。
