# Findings & Decisions

## Requirements
- Locate the root cause of damage output being 0 in the test example.
- Distinguish whether the issue is in the solver, state update, naming, or postprocessing layer.
- Provide evidence-based explanation and, if possible, a fix.

## Research Findings
- The investigation starts from `example/testexample.jl`, which enables `czm_enabled = true` and reads `czm D_max`, `czm D_mean`, `czm δ_max_n [m]`, `czm δ_mean_n [m]`, and `czm n_fractured`.
- Repository memory notes that CZM key names in this branch may be lowercase and that postprocessing already uses guards for missing CZM keys.
- Initial code search showed several mechanical-output paths are gated by `case.opt.mechanicalmodel == "full"`, but the CZM update path itself is separate and still runs when `czm_enabled = true`.
- The test example sets `opt.mechanicalmodel = "none"`, so the remaining suspects are: weak thermal/chemical loading, CZM convergence rollback, or a postprocessing mismatch rather than a hard mechanical gate.
- The cohesive law in `src/Materialmatrix.jl` uses `δ_n_pos = max(0.0, δ_n)` and, for the default `czm_model = "model1"`, damage only increases from positive normal opening. Tangential separation is ignored in that mode.
- The plotted `δ_max_n` / `δ_mean_n` values are maxima of the positive normal opening history, so a flat zero plot means the interface never accumulated positive opening above the threshold, not necessarily that the solver never ran.
- A direct probe of the full Jellyroll case showed the thermal and SOC fields do evolve: `soc_n delta max = 0.8870214911440929`, `soc_p delta max = 0.5475957242778637`, and node temperature spans `298.15 K` to `321.58944549539024 K`.
- The same probe showed the final CZM damage state is still identically zero: `D max/mean = 0`, `δ_max_n max/mean = 0`, `δ_max_t max/mean = 0`, `δ_max_eff max/mean = 0`, `accumulated max/mean = 0`, `fractured count = 0`.
- Running `solve_czm_step` directly on the final state raises `MethodError: no method matching DamageState(::Float64, ::Float64, ::Float64, ::Float64, ::Bool, ::Float64)` from `clone_damage_states` in `src/CzmSolve.jl`.
- `Solve.jl` wraps `update_czm_damage!` in `try/catch` and logs the failure only at `@debug`, so this constructor error is silently swallowed during normal runs.

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Trace the full CZM output path before changing code | Avoids patching symptoms without confirming the source of the zero values |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| None yet | Pending investigation |

## Resources
- `example/testexample.jl`
- `src/czm.jl`
- `src/CzmSolve.jl`
- `src/Mechanical.jl`
- `src/PostProcessing.jl`

## Visual/Browser Findings

## Conclusion
## Conclusion
- The zero-damage result is caused by a silent CZM update failure, not by missing thermal/SOC loading.
- The failure originates in `clone_damage_states` in `src/CzmSolve.jl`, which calls a non-existent `DamageState(s.D, ...)` constructor.
- Because `Solve.jl` catches that exception at debug level only, the damage update never commits and the reported `D_max`, `D_mean`, and `δ_max_*` histories stay at zero.
- The next step, if you want to fix rather than just localize, is to replace the invalid clone constructor call with an explicit field-by-field copy or add the matching inner constructor.

## Follow-up CZM Diagnosis
- The first actual coupled CZM step is numerically much harder than the standalone `t=0` probe: with the default `czm_tol = 1e-4`, `solve_czm_step(..., iter_method="basic")` stalls after the residual falls to about `1.86e-4`.
- A manual Newton trace on that first-step state shows the third iteration has no acceptable descent step even when the line search is reduced to `α = 0.00025`; the residual bottoms out around `1.8e-4`.
- The stall is not caused by the BC penalty method or by a visible active-set jump at the cohesive interface: exact BC elimination gives a similar residual floor, and the interface separations remain near zero (`δ_n` on the order of `1e-21`).
- Relaxing the CZM tolerance to `3e-4` lets the first coupled step converge, but the full `example/testexample.jl` run still emits repeated `CZM solve issue` warnings and ends with `D_max = 0`, so the remaining problem is solver robustness across later steps, not the original clone bug.

## Follow-up Newton Branch Probe
- The dedicated branch probe `tangent_branch_probe.jl` shows iterations 1-3 remain entirely in the elastic branch.
- Per-iteration branch counts were stable: `elastic = 320`, `softening = 0`, `failed = 0`, `Dmax = 0`.
- The largest observed separation was only `max δeff = 2.293358e-22`, while the thresholds are `δ_0_n = 3.416667e-10 m` and `δ_c_n = 6.170732e-7 m` (`δ_0_n* = 5.536891e-4`, `δ_c_n* = 1`).
- This rules out constitutive softening or active-set switching inside `bilinear_tangent` as the reason the third Newton direction loses descent.

## Unit Scale Confirmation
- A side-by-side probe separated the two cohesive input scales.
- Physical cohesive inputs are `K_n = 2.4e17`, `K_t = 4.285714e12`, `δ_0_n = 3.416667e-10`, `δ_c_n = 6.170732e-7`, and they produce `K_bc` with `cond2 = 1.543118e19` and a failed first-step solve.
- Normalized cohesive inputs are `K_n ≈ 1.806068e3`, `K_t ≈ 3.225121e-2`, `δ_0_n ≈ 5.536891e-4`, `δ_c_n = 1`, and they reduce `K_bc` to `cond2 = 9.729879e11`.
- With normalized cohesive inputs, the same first-step CZM solve converges in 2 iterations with residual `2.183138e-12`.

## Current Conclusion
- The unit mismatch was the real cause of the first-step Newton stall: feeding physical cohesive parameters into the CZM solve makes the stiffness matrix badly scaled and the line search fail.
- Once the cohesive parameters are normalized, the first CZM step converges and the full coupled solve now accumulates damage (`D_max ≈ 0.9959`, `D_mean ≈ 0.0322`).
- The remaining issue is separate: later in the full solve, `CZM solve issue` warnings still appear, so there is a follow-on robustness problem after the original unit mismatch is removed.

## Current Branch Note
- In the current code path, `SetParams.NormaliseParam` already normalizes the cohesive parameters, so `testexample.jl` is not failing because raw physical CZM units are still being passed into `solve_czm_step`.
- The remaining convergence bottleneck is the non-smooth unilateral cohesive law at `δ_n = 0` in `bilinear_tangent`, amplified by the very stiff cohesive scale and the fact that the real coupled state stays extremely close to the opening/compression switch.

## Conditioning Probe
- Re-running `tangent_branch_probe.jl` after the user's fix still shows the third Newton step failing the line search.
- The global BC stiffness matrix is severely ill-conditioned: `σ_max = 1.764222e19`, `σ_min = 1.143284e0`, `cond2 = 1.543118e19`.
- The free submatrix is equally problematic: `σ_max = 1.764222e19`, `σ_min = 9.907785e-1`, `cond2 = 1.780642e19`.
- So the immediate blocker is numerical conditioning, not a singular matrix and not an active-set jump.

## Load Substep Probe
- A dedicated compare script against the first coupled CZM state shows `load_substep` is not automatically better than `basic`.
- With the current tolerance, `basic` converges in 2 iterations with residual `2.183138e-12`.
- `load_substep` with 2 substeps converges in 4 iterations with residual `1.091569e-12`.
- `load_substep` with 5 substeps converges in 7 iterations but stops at residual `8.713331e-04`, which is accepted because the implementation uses a looser final threshold than `basic`.
- `load_substep` with 10 substeps behaves similarly and stops at residual `4.356665e-04`.
- `load_substep` with 20 substeps reaches residual `5.457843e-13`, so the method only becomes genuinely accurate when the substep count is much larger.
- Conclusion: the current `czm_load_steps = 5` setting is too small to show a real convergence benefit.

## Tangent Continuity Probe
- A dedicated sweep of `bilinear_tangent` around the opening/closing switch shows a hard jump in the normal tangent at `δ_n = 0`.
- In a damaged-but-not-fractured state with `D ≈ 0.999944`, the tangent changes from `K_nn = 1.806068e3` for `δ_n = -1e-9` to `K_nn = 1.010661e-2` at `δ_n = 0` and `δ_n = 1e-9`.
- That jump comes from the current unilateral formulation in `src/Materialmatrix.jl`: `δ_n_pos = max(0.0, δ_n)` in the traction law and the explicit `if δ_n >= 0` branch in the tangent.
- The fracture branch also drops to a small residual floor (`1e-10 * K_n`), but the zero-crossing kink is the more important source of Newton sensitivity.
- Conclusion: `bilinear_tangent` is not continuous by construction at the compression/opening transition, so the remaining Newton instability is a constitutive non-smoothness issue rather than a load-step control issue.
