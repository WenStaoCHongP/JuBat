# parameters/LGM50.jl Dead Code Verification Plan

**Status:** ✅ Completed（保留）
**Layer:** 1 - 参数与变量
**桶:** Delete 候选

**Goal:** grep 验证 `ChooseCell("LG M50")` / `"LGM50"` 是否有调用者；无则删除整个文件。

---

## Tasks

- [ ] **Step 1: grep 验证**

```bash
grep -rn 'ChooseCell("LG M50"\|ChooseCell("LGM50"\|"LGM50"' --include="*.jl" --include="*.md" --include="*.toml" .
```

注意 `LG M50` 与 `LGM50` 两种拼写。

- [ ] **Step 2: 决策与执行（同 Enertech 模板）**

- 若仅 `SetParams.jl:ChooseCell` 自身匹配 → 删整文件 + 移除 ChooseCell 分支 + 移除 JuBat.jl export（如有）
- 若有外部调用者 → 改桶为 Leave alone

- [ ] **Step 3: 跑 `example/minimal_example.jl` 验证**

- [ ] **Step 4: Commit + baseline 更新**

```bash
git add -A
git commit -m "chore(parameters): 删除未使用的 LGM50 参数集（dead code）"
echo "$(date +%F): parameters/LGM50.jl 决策=[删除/保留]" >> Simplify/baseline.md
```

---

## Risk

低；先 grep 再删。

## Execution Result (2026-08-05)

- `example/minimal_example.jl`、`thermal_example.jl`、`mechanical_example.jl` 等多个活跃示例调用 `ChooseCell("LG M50")`。
- 明确为活跃参数集，保留。
