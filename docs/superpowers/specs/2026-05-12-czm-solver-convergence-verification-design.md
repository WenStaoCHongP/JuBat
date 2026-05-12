# CZM Solver Convergence Verification Design

## Goal

Rewrite `tools/verify_czm_standalone.jl` to verify the correctness and convergence
performance of the three CZM iteration methods (`basic`, `load_substep`, `arc_length`)
in `src/CzmSolve.jl` under **pure Mode I (normal)** loading using **Jellyroll.jl
production parameters** on a **production-grade mesh**.

## Scope

- Pure Mode I (normal direction only). No mixed-mode verification.
- Single uniform thermal-chemical load scenario.
- Compare three solver methods on identical load path.
- Output: console summary table only (no plots).

## Parameters (from Jellyroll.jl)

All parameters obtained via `ChooseCell("Jellyroll")` + `SetCase` + `SetParams`,
which handles normalization. Raw physical values for reference:

| Parameter | Value |
|-----------|-------|
| σ_max_n | 92e6 Pa |
| K_n | 1.2e17 Pa/m |
| δ_0_n | σ_max_n / K_n ≈ 7.67e-10 m |
| G_c_n | 6.2 J/m² |
| δ_c_n | 2 * G_c_n / σ_max_n ≈ 1.35e-7 m |

## Effective Material Parameters

These are computed from `compute_czm_effective_params(case)` which averages
across NE and PE layers. For reference:

| Parameter | Approximate Value |
|-----------|-------------------|
| E_eff | ~2.56e10 Pa |
| ν_eff | ~0.24 |
| α_eff (thermal expansion) | ~6.55e-6 /K |
| β_n (NE partial molar vol / 3) | ~1.033e-6 m³/mol |
| β_p (PE partial molar vol / 3) | ~-2.43e-7 m³/mol |

**Implementation note**: Do NOT hardcode these. Use `compute_czm_effective_params`
or extract from `case.param` after `SetParams` normalization. This ensures the
standalone test uses the same code path as production.

## Mesh

- `nθ = 80`, `gsorder = 2` (production-grade)
- Created via `jellyroll_collector_seed_mesh` + `create_czm_mesh`

## Load Scenario

**Single-shot call**: Call `solve_czm_step` once per method with the maximum load
level that should drive elements through complete fracture. No external load ramp.

The three methods handle internal load progression differently:
- **basic**: Applies full load in one shot. Converges or fails. No internal substeps.
- **load_substep**: Internally splits load into `n_load_steps` substeps with adaptive
  step size control.
- **arc_length**: Internally splits load into substeps with arc-length constraint for
  tracking snap-through behavior.

Each method starts from the same initial state (u=0, D=0) with identical load inputs.
This makes the comparison fair and unambiguous.

### Load magnitude

To drive fracture via thermal expansion alone (Δsoc = 0):
- Required dT ≈ 415 K to reach σ_max_n (damage initiation)
- Use dT values in range [0, 500] K across test levels
- `dT_elem` is uniform (same value for all elements)

With combined thermal + SOC loading:
- Δsoc_n ~ 0.1-0.5 contributes ~1e-7 to ~5e-7 strain via β_n
- Use both dT and Δsoc_n to drive loading; set Δsoc_p = 0 for simplicity

**Recommended approach**: Use `compute_czm_effective_params` to get α_eff, then
calculate the dT range needed: `dT_init = δ_0_n * K_n / (E_eff * α_eff)` for damage
initiation, `dT_fracture ≈ 4× dT_init` for complete fracture. Set Δsoc_n_elem = 0,
Δsoc_p_elem = 0 to isolate thermal driving.

## Comparison Metrics

| Metric | Description |
|--------|-------------|
| converged | Whether the solver converged (bool) |
| iterations | Total Newton iterations consumed |
| D_max | Maximum damage across all cohesive elements |
| D_mean | Mean damage across all cohesive elements |
| residual_norm | Final residual norm |

Test multiple load levels to characterize convergence behavior across the
traction-separation curve: elastic (no damage), damage initiation, softening,
near-complete fracture.

## Output Format

Console-printed summary table comparing three methods across load levels.

Example:

```
CZM Solver Convergence Comparison (Mode I, nθ=80)
═══════════════════════════════════════════════════════════════════
dT (K)  | basic               | load_substep        | arc_length
────────┼─────────────────────┼─────────────────────┼───────────
  100   | ✓ 5 it D=0.00 r=1e-12| ✓ 5 it D=0.00 r=1e-12| ✓ 5 it D=0.00
  200   | ✓ 12 it D=0.15 r=3e-5| ✓ 8 it D=0.15 r=2e-5| ✓ 8 it D=0.15
  300   | ✓ 18 it D=0.52 r=8e-5| ✓ 10 it D=0.52 r=4e-5| ✓ 9 it D=0.52
  400   | FAIL (20 it)         | ✓ 14 it D=0.89 r=6e-5| ✓ 12 it D=0.88
  500   | FAIL (20 it)         | ✓ 18 it D=0.97 r=9e-5| ✓ 15 it D=0.96
═══════════════════════════════════════════════════════════════════
```

Optional: write the same table to `output/czm_standalone/solver_comparison.txt`.

Note: When basic fails, D_max shown is from the state at failure (non-zero if
damage had begun in previous converged iterations). The script runs each load
level independently (fresh state), so failure at one level does not affect others.

## Implementation Structure

The rewritten `tools/verify_czm_standalone.jl` will contain:

1. **Setup**: Build param via `ChooseCell("Jellyroll")`, build mesh via
   `jellyroll_collector_seed_mesh` + `create_czm_mesh`, compute effective params
   via `compute_czm_effective_params` or equivalent extraction from normalized params.
2. **Cache builder**: Call `build_czm_cache(czm_mesh, E_eff, ν_eff, param)` once,
   pass to each solver call.
3. **Load level loop**: For each dT level, call `solve_czm_step` with each of the
   three methods, starting from u=0 with fresh damage states each time.
4. **Record & compare**: Collect metrics into a table.
5. **Report**: Print summary table to console.

### Key implementation details

- Each method call gets a fresh `czm_mesh` clone (with fresh `damage_states`)
  so load levels are independent.
- Use `build_czm_cache` directly (not `ensure_czm_cache` which requires a `Case`).
- `F_ext = zeros(ndof)` (no external mechanical force; purely thermo-chemical driven).
- Solver tolerance: `tol = 1e-4` (matching `opt.czm_tol` default).
- `max_iter = 200`, `n_load_steps = 50` for load_substep and arc_length.
- Viscous regularization disabled: `visc_beta = 1.0`.
- **basic** has no internal load substeps; it applies full load in one step. This
  is a known limitation — basic may fail on large load increments where
  load_substep/arc_length succeed.

### Parameter construction

The standalone script should use the standard code path to ensure normalization
consistency:

```julia
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.czm_enabled = true
case = JuBat.SetCase(param_dim, opt)
# After SetCase, case.param contains normalized parameters
# E_eff, ν_eff etc. can be computed from case.param
```

## Acceptance Criteria

1. Script runs without errors: `julia --project=. tools/verify_czm_standalone.jl`
2. Summary table is printed to console
3. At least one of the three methods converges through complete fracture
4. Methods that converge show consistent D_max values (within 5% of each other)
5. No dependency on electrochemical-thermal main chain (standalone)
