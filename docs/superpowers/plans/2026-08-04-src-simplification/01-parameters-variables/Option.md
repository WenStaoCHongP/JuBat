# Option.jl Audit Plan

**Status:** ⬜ Pending
**Layer:** 1 - 参数与变量
**桶:** Leave alone

**Goal:** 审查 `Option.jl`；验证 `solveType` 的 `forward`/`backward` 分支是否被实际使用，未使用可清理。

---

## Audit Checklist

- [ ] **Step 1: grep `solveType` 实际取值**

```bash
grep -rn "solveType\s*==" src/ example/ test/ | sort | uniq -c
```

Expected: 主要看到 `"Crank-Nicolson"`；如 `forward`/`backward` 0 次匹配，可清分支。

- [ ] **Step 2: 通读 `Option.jl` (88 行)**，确认：
  - 所有字段默认值合理
  - `CycleOption` + `PhaseType` 仍被 CycleSolver 使用
  - 无 try/catch 或兼容入口

- [ ] **Step 3: 决策**

- 若 `forward`/`backward` 完全未用 → 在 spec §7.5 补录"分支可清"；本 plan 不动（避免破坏公共 API）
- 若有使用 → 无动作

- [ ] **Step 4: baseline 记录**

```bash
echo "$(date +%F): Option.jl 已审查，桶=Leave alone" >> Simplify/baseline.md
```

---

## Result

无修改（保守起见保留所有 Option 字段，公共 API 不变）。
