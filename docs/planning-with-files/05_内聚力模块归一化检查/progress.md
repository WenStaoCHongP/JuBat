# Progress Log

## Session: 2026-04-22

### Phase 1: Requirements & Discovery
- **Status:** complete
- **Started:** 2026-04-22
- Actions taken:
  - Created task_plan.md, findings.md, and progress.md for this check
  - Scoped the investigation to cohesive normalization entrypoints and usage paths
- Files created/modified:
  - docs/planning-with-files/05_内聚力模块归一化检查/task_plan.md
  - docs/planning-with-files/05_内聚力模块归一化检查/findings.md
  - docs/planning-with-files/05_内聚力模块归一化检查/progress.md

### Phase 2: Code Path Tracing
- **Status:** complete
- Actions taken:
  - Traced NormaliseParam in src/SetParams.jl
  - Traced cohesive parameter definitions in src/parameters/Jellyroll.jl
  - Traced consumption in src/Materialmatrix.jl
  - Confirmed scale definition and normalization formulas for cohesive fields
  - Inspected src/PostProcessing.jl and confirmed actual CZM result restoration scales
- Files created/modified:
  - docs/planning-with-files/05_内聚力模块归一化检查/findings.md
  - docs/planning-with-files/05_内聚力模块归一化检查/progress.md

## Test Results
| Test | Input | Expected | Actual | Status |
|------|-------|----------|--------|--------|
| Planning files created | create_file calls | Files created successfully | Success | ✓ |
| Cohesive normalization on Jellyroll | ChooseCell("Jellyroll") + NormaliseParam | Finite normalized cohesive values | scale=(8.2e7, 6.17e-7, 50.6, 1.33e14); coh_norm=(1.0, 1.0, 18.06) | ✓ |
| Cohesive normalization on LG M50 | ChooseCell("LG M50") + NormaliseParam | Finite normalized cohesive values | zero scale fields and NaN normalized cohesive values | ✗ |

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
| 2026-04-22 | type Scale has no field ? during direct Unicode field access in Julia check | 1 | Switch to getfield(Symbol(...)) for σ_czm / δ_czm / K_czm |
| 2026-04-22 | NormaliseParam on ChooseCell("LG M50") produced NaN cohesive values because scale fields were zero | 2 | Treat LG M50 as an invalid validation preset for CZM and switch to a CZM-enabled cell preset |
| 2026-04-22 | type Cohesive has no field Γ_c_t during shear-side runtime check | 3 | Use getfield(..., :G_c_t) instead of a Greek-letter field name |

## Notes
- LG M50 is not a valid cohesive normalization preset in the current codebase.
- Next runtime check should use ChooseCell("Jellyroll").
- Jellyroll passes the normalization check, but its shear parameters make δ_0_t > δ_c_t, so the CZM shear branch degenerates to an abrupt damage jump.
- CZM traction postprocessing is restored with E_n/E_p, not σ_czm; this should be treated as the current implementation contract until the doc is updated.

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Complete |
| Where am I going? | Deliver summary to user |
| What's the goal? | Check cohesive parameter normalization consistency |
| What have I learned? | NormaliseParam scaling is consistent on Jellyroll; LG M50 is not a CZM validation preset; shear parameters make δ_0_t > δ_c_t |
| What have I done? | Created planning files, traced the normalization path, and validated with Julia runtime checks |
