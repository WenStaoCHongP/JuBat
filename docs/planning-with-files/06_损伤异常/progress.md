# Progress Log

## Session: 2026-04-21

### Phase 1: Requirements & Discovery
- **Status:** complete
- **Started:** 2026-04-21
- Actions taken:
  - Read planning-with-files skill instructions and templates.
  - Created planning files for the damage-zero investigation.
  - Captured initial scope from `example/testexample.jl` and repository memory.
  - Identified an initial hypothesis: CZM damage update paths may be gated by `mechanicalmodel == "full"` while the test example leaves it at `"none"`.
  - Refined the hypothesis after code tracing: the default CZM law only accumulates damage from positive normal opening, and the current run appears to produce zero positive opening on the interface.
  - Wrote the final conclusion into the planning files.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/task_plan.md`
  - `docs/planning-with-files/06_损伤异常/findings.md`
  - `docs/planning-with-files/06_损伤异常/progress.md`

### Phase 2: [Title]
- **Status:** complete
- Actions taken:
  - Traced CZM damage flow from solver update to cohesive law and postprocessing.
- Files created/modified:
  -

### Phase 3: [Title]
- **Status:** complete
- Actions taken:
  - Verified that the cohesive law only increases damage for positive normal opening.
  - Confirmed the zero plots match a non-opening interface response.
- Files created/modified:
  -

### Phase 4: [Title]
- **Status:** complete
- Actions taken:
  - Captured the conclusion in `task_plan.md` and `findings.md`.
  - Closed out the investigation as a report-only outcome; no code fix was applied.
  - Ran a direct CZM probe showing thermal/SOC loading changes but zero damage-state growth.
  - Confirmed the direct `solve_czm_step` path fails in `clone_damage_states` with a `DamageState` constructor error.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/task_plan.md`
  - `docs/planning-with-files/06_损伤异常/findings.md`
  - `docs/planning-with-files/06_损伤异常/tmp_czm_probe.jl`

## Test Results
| Test | Input | Expected | Actual | Status |
|------|-------|----------|--------|--------|
| Full Jellyroll probe | `tmp_czm_probe.jl` | Nonzero damage or at least nonzero separation history | SOC and temperature change, but `D`, `δ_max_n`, `δ_max_t`, `δ_max_eff`, and `accumulated_damage` all remain 0 | ✓ |
| Direct CZM step | `solve_czm_step(...)` on final state | Return separation/damage history | Fails with `MethodError: no method matching DamageState(::Float64, ...)` from `clone_damage_states` | ✓ |

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
| None yet | - | 1 | - |
| 2026-04-21 | Invalid array length from accidental `memory` call | 1 | Ignored and continued with the investigation |
| 2026-04-21 | `LoadError: type JellyrollMesh has no field Jellyroll_czm` in `tools/verify_czm_system.jl` | 1 | Switched to a direct probe against `mesh_data.thermal2D` |
| 2026-04-21 | `KeyError: key "czm separation normal [m]" not found` | 1 | Confirmed that this key is not emitted by current `PostProcessing` |
| 2026-04-21 | `MethodError: no method matching DamageState(::Float64, ...)` | 1 | Identified as the root cause inside `clone_damage_states` |
| 2026-04-21 | `syntax: more than one semicolon in argument list` | 1 | Fixed the direct CZM probe call syntax |

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Phase 4 complete in task_plan.md |
| Where am I going? | No further action required unless the user wants the clone fix patched into the code |
| What's the goal? | Locate why damage is reported as 0 |
| What have I learned? | See findings.md |
| What have I done? | Created planning files, traced the CZM path, ran direct probes, and recorded the final conclusion |

## Session: 2026-04-21 Follow-up CZM Convergence Probe

### Phase 5: First-step solver diagnosis
- **Status:** complete
- Actions taken:
  - Built `first_czm_step_probe.jl` to compare the initial `t=0` CZM solve, the first coupled step, and multiple tolerances.
  - Confirmed the first coupled step stalls at the default tolerance: `tol = 1e-4` fails, while `tol = 3e-4` and `tol = 1e-3` accept the same state after 2 iterations.
  - Ran a manual Newton trace and found the third iteration cannot find a descending line-search step even after shrinking to `α = 0.00025`.
  - Checked exact BC elimination and found a similar residual floor, so the BC penalty was not the cause.
  - Verified the first-step interface separations stay essentially zero (`δ_n` around `1e-21`), so the stall is not driven by a visible opening jump.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/first_czm_step_probe.jl`

### Phase 6: Full-solve tolerance probe
- **Status:** complete
- Actions taken:
  - Ran a full `Solve(case)` probe with `opt.czm_tol = 3e-4`.
  - Confirmed the run still emits repeated `CZM solve issue` warnings later in the simulation.
  - Confirmed the final damage outputs remain zero (`D_max = 0`, `D_mean = 0`, `n_fractured = 0`).
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/solve_tol_3e-4_probe.jl`

### Session: 2026-04-21 Follow-up Newton Branch Probe
- **Status:** complete
- Actions taken:
  - Added `tangent_branch_probe.jl` to print branch counts and threshold distances during the first coupled CZM step.
  - Ran the probe through iterations 1-3.
  - Verified all 320 cohesive elements remain in the elastic branch with no softening or fracture activation.
  - Confirmed the third Newton direction loss is not caused by active-set switching in the cohesive law.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/tangent_branch_probe.jl`

### Session: 2026-04-21 Conditioning Check
- **Status:** complete
- Actions taken:
  - Extended `tangent_branch_probe.jl` to compute `svdvals` for `K_bc` and the free submatrix.
  - Re-ran the probe after the user's unit-mismatch fix.
  - Recorded that both `K_bc` and `K_free` remain extremely ill-conditioned, with condition numbers around `1e19`.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/tangent_branch_probe.jl`

### Session: 2026-04-21 Unit Scale Comparison
- **Status:** complete
- Actions taken:
  - Added `unit_scale_probe.jl` to compare physical and normalized cohesive parameters side by side.
  - Verified that physical cohesive inputs produce `K_bc` with `cond2 ≈ 1.54e19` and a stalled solve.
  - Verified that normalized cohesive inputs reduce `K_bc` to `cond2 ≈ 9.73e11` and the same first CZM step converges in 2 iterations.
  - Re-ran the full `Solve(case)` probe and confirmed the final damage values are now nonzero (`D_max ≈ 0.9959`, `D_mean ≈ 0.0322`).
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/unit_scale_probe.jl`

### Session: 2026-04-21 Load Substep Comparison
- **Status:** complete
- Actions taken:
  - Added `load_substep_compare.jl` to compare `basic` and `load_substep` on the same first coupled CZM state.
  - Confirmed `load_substep` with 5 and 10 substeps can terminate with residuals around `8.7e-4` and `4.4e-4`, respectively, because the current final acceptance threshold is looser than the basic solver.
  - Confirmed `load_substep` only becomes as accurate as `basic` when the substep count is raised to 20.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/load_substep_compare.jl`

### Session: 2026-04-21 Tangent Continuity Probe
- **Status:** complete
- Actions taken:
  - Added `tangent_continuity_probe.jl` to sweep `bilinear_tangent` across loading and unloading states.
  - Confirmed a hard tangent jump at `δ_n = 0` in a damaged-but-not-fractured state: `K_nn` drops from `1.806068e3` for negative opening to `1.010661e-2` at zero opening and positive infinitesimal opening.
  - Confirmed the fracture branch also uses a small residual floor, but the dominant non-smoothness is the unilateral opening/compression switch.
- Files created/modified:
  - `docs/planning-with-files/06_损伤异常/tangent_continuity_probe.jl`

### Final Diagnosis
- **Status:** complete
- Outcome:
  - Current `testexample.jl` does not suffer from a remaining physical/normalized cohesive-unit mismatch in the CZM entry path.
  - The active convergence blocker is the non-smooth cohesive tangent at `δ_n = 0`, which is made much harder to solve by the very stiff cohesive scale and the near-zero separation state produced by the real coupled loading.
