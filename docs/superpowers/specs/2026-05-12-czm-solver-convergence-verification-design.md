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

| Parameter | Value |
|-----------|-------|
| σ_max_n | 92e6 Pa |
| K_n | 1.2e17 Pa/m |
| δ_0_n | σ_max_n / K_n ≈ 7.67e-10 m |
| G_c_n | 6.2 J/m² |
| δ_c_n | 2 * G_c_n / σ_max_n ≈ 1.35e-7 m |

## Mesh

- `nθ = 80`, `gsorder = 2` (production-grade)
- Created via `jellyroll_collector_seed_mesh` + `create_czm_mesh`

## Load Scenario

Single uniform load: all elements receive identical `dT_elem`, `Δsoc_n_elem`,
`Δsoc_p_elem`. Load is ramped through multiple levels (e.g. 5-10 steps) from
zero to a level sufficient to drive all cohesive elements from elastic through
damage initiation to complete fracture (δ > δ_c_n).

At each load level, call `solve_czm_step` with each of the three methods and
record convergence metrics.

## Comparison Metrics

| Metric | Description |
|--------|-------------|
| converged | Whether the solver converged (bool) |
| iterations | Total Newton iterations consumed |
| D_max | Maximum damage across all cohesive elements |
| D_mean | Mean damage across all cohesive elements |
| residual_norm | Final residual norm |

## Output Format

Console-printed summary table comparing three methods across load levels.

Example:

```
CZM Solver Convergence Comparison (Mode I, nθ=80)
═══════════════════════════════════════════════════
Load level | basic        | load_substep | arc_length
───────────┼──────────────┼──────────────┼───────────
  1/8      | conv ✓ 5 it  | conv ✓ 5 it  | conv ✓ 5 it
  2/8      | conv ✓ 8 it  | conv ✓ 6 it  | conv ✓ 6 it
  3/8      | conv ✓ 12 it | conv ✓ 7 it  | conv ✓ 7 it
  4/8      | FAIL         | conv ✓ 10 it | conv ✓ 9 it
  5/8      | FAIL         | conv ✓ 14 it | conv ✓ 11 it
  ...
═══════════════════════════════════════════════════
Final: D_max(basic)=0.00, D_max(load_substep)=0.97, D_max(arc_length)=0.95
```

Optional: write the same table to `output/czm_standalone/solver_comparison.txt`.

## Implementation Structure

The rewritten `tools/verify_czm_standalone.jl` will contain:

1. **Setup**: Build mesh and cohesive params from Jellyroll.jl
2. **Cache builder**: Build `CZMAssemblyCache` once, clone for each method
3. **Load ramp loop**: For each load level, call `solve_czm_step` with
   `iter_method="basic"`, `"load_substep"`, `"arc_length"`
4. **Record & compare**: Collect metrics into a table
5. **Report**: Print summary table to console

Key implementation details:
- Each method gets an independent copy of `czm_mesh` (with fresh `damage_states`)
  so methods don't interfere with each other.
- `E_eff` and `ν_eff` computed from existing parameter helpers.
- `α_eff`, `β_n`, `β_p` taken from param or set to reasonable thermal expansion
  coefficients.
- `F_ext = zeros(ndof)` (no external mechanical force; purely thermo-chemical driven).
- Solver tolerance: `tol = 1e-4` (matching `opt.czm_tol` default).
- `max_iter = 200`, `n_load_steps = 50` for load_substep and arc_length.

## Acceptance Criteria

1. Script runs without errors: `julia --project=. tools/verify_czm_standalone.jl`
2. Summary table is printed to console
3. At least one of the three methods converges through complete fracture
4. Methods that converge show consistent D_max values (within 5% of each other)
5. No dependency on electrochemical-thermal main chain (standalone)
