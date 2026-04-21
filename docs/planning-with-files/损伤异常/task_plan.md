# Task Plan: 定位损伤为 0 的原因

## Goal
定位 `example/testexample.jl` 运行后损伤相关输出始终为 0 的根因，并给出可验证的修复建议。

## Current Phase
Phase 7

## Phases

### Phase 1: Requirements & Discovery
- [x] Understand user intent
- [x] Identify constraints and requirements
- [x] Document findings in findings.md
- **Status:** complete

### Phase 2: Code Path Tracing
- [x] Trace damage variables from solver to postprocessing
- [x] Identify where values are initialized or overwritten
- [x] Compare with repository conventions and known issues
- **Status:** complete

### Phase 3: Root Cause Verification
- [x] Reproduce the zero-damage path in code
- [x] Check variable naming and key guards
- [x] Check model activation conditions and output wiring
- **Status:** complete

### Phase 4: Fix or Report
- [x] Apply fix if root cause is clear
- [x] Validate with targeted checks
- [x] Summarize evidence and next steps
- **Status:** complete

## Key Questions
1. Is the damage variable produced by the CZM solver, postprocessing, or both?
2. Are the expected output keys in this branch named differently from the script?
3. Is the damage state actually evolving, or is only the reported output zero?

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Use planning-with-files for this investigation | Needed to keep tracing evidence and avoid losing context across many reads |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| `DamageState(::Float64, ...)` MethodError inside `clone_damage_states` | 1 | Confirmed as root cause by direct CZM probe; the caller currently swallows the failure |

## Conclusion
The zero-damage result comes from a silent CZM update failure, not from a lack of thermal/SOC loading. `Solve.jl` catches errors from `update_czm_damage!` at debug level only, and `clone_damage_states` in `src/CzmSolve.jl` calls a non-existent `DamageState(s.D, ...)` constructor. That throws during every CZM update, so damage states never advance and the saved `D_max`, `D_mean`, and `δ_max_*` histories remain zero.

Follow-up probes show a separate CZM convergence problem in the first coupled time step: the basic solver reaches a residual floor near `1.8e-4`, and the later integrated run still stalls even when `czm_tol` is relaxed to `3e-4`. That is a distinct issue from the original zero-damage bug.

### Phase 5: Newton Branch Diagnosis
- [x] Probe branch classification across Newton iterations
- [x] Compare separations with damage thresholds
- [x] Rule out active-set switching as the stall cause
- **Status:** complete

### Phase 6: Conditioning Check
- [x] Measure K_bc condition number and smallest singular value
- [x] Compare K_bc with free-submatrix conditioning
- [x] Confirm whether the stall is driven by numerical ill-conditioning
- **Status:** complete

### Phase 7: Post-fix Validation
- [x] Compare physical vs normalized cohesive inputs
- [x] Confirm normalized CZM solve converges
- [x] Confirm the full solve now accumulates damage
- [x] Check whether `load_substep` with 5 steps actually improves residuals
- [x] Probe `bilinear_tangent` continuity near `δ_n = 0`
- [ ] Trace the later-step CZM warnings separately
- [ ] Decide whether to smooth the unilateral opening/closing switch or keep the current non-smooth tangent
- **Status:** in-progress

## Follow-up Note
- The third-Newton stall is not caused by cohesive softening or fracture activation.
- The branch probe shows all 320 cohesive elements remain elastic through iterations 1-3, with `softening = 0`, `failed = 0`, and `Dmax = 0`.
- Next suspect is solver conditioning or a residual/BC mismatch, not active-set switching in `bilinear_tangent`.

## Conditioning Note
- `K_bc` is extremely ill-conditioned in the current probe: `σ_max = 1.764222e19`, `σ_min = 1.143284e0`, `cond2 = 1.543118e19`.
- The free-DOF submatrix is similarly scaled: `σ_max = 1.764222e19`, `σ_min = 9.907785e-1`, `cond2 = 1.780642e19`.
- The matrix is not singular, but the scaling is severe enough to make the Newton direction very sensitive.

## Post-fix Note
- The new unit-scale probe shows the physical cohesive parameter set is the source of the `1e19` conditioning: `K_n = 2.4e17`, `K_t = 4.3e12`.
- The normalized cohesive parameters are orders of magnitude smaller (`K_n ≈ 1.81e3`, `K_t ≈ 3.23e-2`), and the same first CZM step then converges in 2 iterations.
- The full solve now reaches `D_max ≈ 0.9959` and `D_mean ≈ 0.0322`, so the original zero-damage symptom is no longer present.
- The load-substep probe shows `czm_load_steps = 5` is not enough to guarantee an accurate solve: on the first coupled state it can stop at residual `8.7e-4`, while `czm_load_steps = 20` is needed to recover a residual near machine precision.
- The tangent continuity probe shows a hard jump at `δ_n = 0` in damaged-but-not-fractured states (`K_nn` drops from `1.806068e3` in compression to `1.010661e-2` at zero opening), so the remaining instability is now traced to the unilateral constitutive kink.

## Notes
- Start from `example/testexample.jl`, then trace CZM outputs in `src/czm.jl`, `src/CzmSolve.jl`, `src/Mechanical.jl`, and postprocessing.
- Re-read this plan before any major branch in the investigation.
