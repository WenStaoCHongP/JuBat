# SetCase.jl Compatibility Ctor Removal Plan

**Status:** ✅ Completed（保留）
**Layer:** 1 - 参数与变量
**桶:** Leave alone + 兼容入口验证（spec §7.3）

**Goal:** grep 验证 5 参数兼容构造器 `Case` (line 114) 是否仍被调用；无则删除。

---

## Files

- Modify: `src/SetCase.jl:114-115`（5 参数兼容构造器）

---

## Tasks

- [ ] **Step 1: grep 5 参数构造器调用点**

```bash
grep -rn "Case(" --include="*.jl" src/ example/ test/ | grep -v "SetCase\|CycleResult\|PhaseResult\|MultiSPMeLayout"
```

5 参数签名为 `Case(param_dim, param, opt, mesh, index)`；3 参数为 `SetCase(param_dim, opt, y0)`。

- [ ] **Step 2: 决策**

- 若 5 参数版仅自身定义、无外部调用者 → 删除 line 114-115
- 若有调用者 → 保留并加注释 "5 参数兼容入口；调用者：[列出]"

- [ ] **Step 3: 跑主线 example 验证**

Run: `julia example/minimal_example.jl`
Expected: PASS

- [ ] **Step 4: Commit + baseline**

```bash
git add src/SetCase.jl
git commit -m "chore(SetCase): 删除无调用者的 5 参数兼容构造器"
echo "$(date +%F): SetCase.jl 决策=[删除/保留 5 参数构造器]" >> Simplify/baseline.md
```

---

## Risk

低；grep 是可靠判据。

## Execution Result (2026-08-05)

- `src/SetCase.jl:95` 直接调用五参数 `Case(param_dim, param, opt, mesh, index)`，因此不是 dead compatibility entry。
- 保留构造器；未修改源码，无需重复运行行为基线。
