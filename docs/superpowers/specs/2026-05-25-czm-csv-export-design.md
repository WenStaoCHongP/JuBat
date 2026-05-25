# CZM Post-Processing CSV Export Design

**Date:** 2026-05-25
**Status:** Reviewed
**Scope:** Export cycling simulation results to CSV for Python/Matplotlib post-processing

---

## 1. Goal

Provide an independent post-processing function `export_cycling_csv` that writes `solve_cycling` results to a set of CSV files (tidy/long format), enabling Python analysis of:
- Element-level electrochemical currents and thermal sources
- Node-level temperature and displacement fields
- Cohesive element damage, separation, traction evolution
- Driving force decomposition (thermal vs diffusion strain)
- CZM solver convergence diagnostics

This supports root-cause analysis of cycle-count vs. CZM damage relationships.

## 2. File Organization (7 files)

All files use long/tidy format with a common set of index columns: `time_s`, `cycle`, `phase`.

### 2.1 `cycle_summary.csv`

Per-phase cycle-level aggregates.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | discharge/charge/rest1/rest2 |
| `capacity_ah` | Float | Ah | Phase capacity |
| `soh` | Float | [0,1] | State of health |
| `D_max` | Float | [0,1] | Maximum damage |
| `D_mean` | Float | [0,1] | Mean damage |
| `n_fractured` | Int | - | Fully fractured cohesive elements |
| `T_max_K` | Float | K | Maximum temperature |
| `T_mean_end_K` | Float | K | Mean temperature at phase end |
| `V_start` | Float | V | Starting voltage |
| `V_end` | Float | V | Ending voltage |

**Source:** `CyclingResult.cycle_results[i].charge/discharge/etc.` fields + `result.soh[i]`.

### 2.2 `element_currents.csv`

Per-element electrochemical and thermal data at each time step.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `elem_id` | Int | - | Thermal element index |
| `I_e` | Float | A | Element branch current (× `scale.I_typ`) |
| `area` | Float | m^2 | Element area (× `scale.L^2`) |
| `T_e` | Float | K | Element mean temperature (× `scale.T_ref`) |
| `soc_n` | Float | [0,1] | Negative electrode SOC (intrinsic dimensionless) |
| `soc_p` | Float | [0,1] | Positive electrode SOC (intrinsic dimensionless) |
| `eta_n` | Float | V | Negative overpotential (× `scale.phi`) |
| `eta_p` | Float | V | Positive overpotential (× `scale.phi`) |
| `q_rxn_ne` | Float | W/m^3 | Negative reaction heat (× `scale.q`) |
| `q_sp` | Float | W/m^3 | Separator ohmic heat (× `scale.q`) |
| `q_rxn_pe` | Float | W/m^3 | Positive reaction heat (× `scale.q`) |

**Source:** `variables["thermal2D element current"]`, `variables["thermal2D element soc_n"]`, etc. **All values are normalized** in the `variables` dictionary; must be denormalized using `case.param.scale` factors before writing to CSV. See Section 5.3 for the complete conversion table.

### 2.3 `node_temperature.csv`

Per-node temperature field at each time step.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `node_id` | Int | - | Node index |
| `x` | Float | m | Node x coordinate |
| `y` | Float | m | Node y coordinate |
| `T_K` | Float | K | Temperature |

**Source:** `variables["thermal2D temperature at nodes"]`, node coordinates from `case.mesh["thermal2D"]`.

### 2.4 `cohesive_damage.csv`

Per-cohesive-element damage state at each CZM update step.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `coh_id` | Int | - | Cohesive element index |
| `length` | Float | m | Element length |
| `D` | Float | [0,1] | Damage value |
| `delta_n` | Float | m | Normal separation |
| `delta_t` | Float | m | Tangential separation |
| `T_n` | Float | Pa | Normal traction |
| `T_t` | Float | Pa | Tangential traction |
| `fractured` | Bool | - | Whether fully fractured |
| `theta_deg` | Float | deg | Angular position |

**Source:** CZM snapshots saved during `solve_cycling` (requires enhancement, see Section 4). **theta_deg** is computed from the element midpoint coordinates using `atan(y, x) * 180/π`.

**Unit conversion:**
- `delta_n`, `delta_t`: multiply by `scale.r0` (not `scale.L`) — confirmed in `PostProcessing.jl:110`
- `T_n`, `T_t`: multiply by `scale.σ_czm`
- `length`: multiply by `scale.L`
- `D`, `fractured`: already physical units, no scaling needed

### 2.5 `node_displacement.csv`

Per-node displacement field at each CZM update step.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `node_id` | Int | - | Node index |
| `x` | Float | m | Node x coordinate |
| `y` | Float | m | Node y coordinate |
| `ux` | Float | m | x-displacement |
| `uy` | Float | m | y-displacement |

**Source:** CZM displacement snapshots (requires enhancement, see Section 4).

### 2.6 `cohesive_driving_force.csv`

Per-cohesive-element driving force decomposition.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `coh_id` | Int | - | Cohesive element index |
| `dT_neighbor` | Float | K | Temperature rise of adjacent elements |
| `dsoc_n_neighbor` | Float | [0,1] | Negative electrode SOC change |
| `dsoc_p_neighbor` | Float | [0,1] | Positive electrode SOC change |
| `eps_0_thermal` | Float | - | Thermal strain component alpha*DeltaT |
| `eps_0_diffusion` | Float | - | Diffusion strain component beta*DeltaSOC |
| `eps_0_total` | Float | - | Total initial strain |

**Source:** Computed from `variables` thermal/SOC data + CZM element-to-thermal-element mapping.

### 2.7 `czm_solver_diagnostics.csv`

Per-step CZM solver convergence information.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `time_s` | Float | s | Physical time |
| `cycle` | Int | - | Cycle number |
| `phase` | String | - | Phase name |
| `converged` | Bool | - | Whether solver converged |
| `iterations` | Int | - | Newton iterations used |
| `residual_norm` | Float | - | Final residual norm |
| `method` | String | - | Solver method (basic/load_substep/arc_length) |

**Source:** CZM solver return info (requires enhancement, see Section 4).

## 3. Function Interface

```julia
"""
    export_cycling_csv(result::CyclingResult, case::Case, czm_mesh::CohesiveMesh;
                       output_dir="output/csv", overwrite=false)

Export solve_cycling results to CSV files for post-processing.

# Arguments
- `result`: CyclingResult from solve_cycling
- `case`: Case object with mesh and parameters
- `czm_mesh`: CohesiveMesh with damage history

# Keyword Arguments
- `output_dir`: Output directory (default "output/csv")
- `overwrite`: Overwrite existing files (default false)

# Generated Files
7 CSV files: cycle_summary, element_currents, node_temperature,
cohesive_damage, node_displacement, cohesive_driving_force, czm_solver_diagnostics
"""
function export_cycling_csv(result, case, czm_mesh;
                            output_dir="output/csv", overwrite=false)
```

Usage in `coupled_czm_thermal_example.jl`:
```julia
result = JuBat.solve_cycling(case, cycle_opt, czm_mesh; verbose=true, save_detailed=true)
JuBat.export_cycling_csv(result, case, czm_mesh; output_dir="output/csv/czm_study_1")
```

## 4. Required Data Enhancement

The current `solve_cycling` with `save_detailed=true` does not save per-step CZM snapshots. The following enhancements are needed in `src/CycleSolver.jl`:

### 4.1 CZM Snapshot Structure

```julia
mutable struct CZMSnapshot
    time_s::Float64
    cycle::Int
    phase::String
    displacement::Vector{Float64}        # ndof-length vector (normalized)
    damage::Vector{Float64}              # n_coh-length vector
    separation_n::Vector{Float64}        # n_coh-length vector (normalized)
    separation_t::Vector{Float64}        # n_coh-length vector (normalized)
    traction_n::Vector{Float64}          # n_coh-length vector (normalized)
    traction_t::Vector{Float64}          # n_coh-length vector (normalized)
    converged::Bool
    iterations::Int
    residual_norm::Float64
    method::String                       # "basic", "load_substep", or "arc_length"
end
```

### 4.2 Storage in CyclingResult

Add a field to `CyclingResult`:
```julia
czm_snapshots::Vector{CZMSnapshot}    # filled when save_detailed=true
```

### 4.3 Collection Point

The CZM solve happens inside `CallModel_MultiSPMe` (in `src/Solve.jl`), not directly in `solve_cycling`. The snapshot collection must be placed in `CallModel_MultiSPMe` after the CZM solve call, inside the CZM update block.

**Implementation approach:** Add an optional `czm_snapshots` vector parameter to `CallModel_MultiSPMe`. When provided (non-nothing), after each CZM solve, push a snapshot. `solve_cycling` creates this vector and passes it through.

Alternatively, the CZM result data is already stored in the `variables` dictionary via `CzmPostProcess.czmo_output_to_variables`. We can extract CZM per-element data from `variables` keys like `"czm D_max"`, `"czm separation normal"`, etc. However, these are aggregate statistics, not per-element data. **The per-element CZM snapshot is required for the cohesive_damage.csv and node_displacement.csv files.**

## 5. Implementation Notes

### 5.1 File Location

New file: `src/CsvExport.jl`, included in `src/JuBat.jl`.

### 5.2 Dependencies

Only standard library: `DelimitedFiles` (for `writedlm`) or `CSV.jl` if available. To minimize dependencies, use `open(..., "w")` with `println` for CSV writing.

### 5.3 Unit Conversion

All values in the `variables` dictionary are **normalized** (dimensionless). The following table specifies the conversion for each column:

| Variable | Raw Source | Scale Factor | Physical Unit |
|----------|-----------|-------------|---------------|
| Temperature | `T_nodes` | × `scale.T_ref` | K |
| Element current | `I_e` | × `scale.I_typ` | A |
| Area | `A_elem` (from mesh) | × `scale.L^2` | m^2 |
| Node coordinates | `mesh.node` | × `scale.L` | m |
| Time | `t` (normalized) | × `scale.t0` | s |
| Voltage | `V` (normalized) | × `scale.phi` | V |
| Overpotential | `eta_n/p` | × `scale.phi` | V |
| Heat source | `q_*` | × `scale.q` | W/m^3 |
| SOC | `soc_n/p` | — | dimensionless [0,1] |
| CZM displacement | `u` | × `scale.L` | m |
| CZM separation | `δ_n/t` | × `scale.r0` | m |
| CZM traction | `T_n/t` | × `scale.σ_czm` | Pa |
| CZM damage | `D` | — | dimensionless [0,1] |
| CZM length | `L_elem` | × `scale.L` | m |
| Thermal strain | `α·ΔT` | — | dimensionless |
| Diffusion strain | `β·ΔSOC` | — | dimensionless |

**Key note:** CZM separation uses `scale.r0` (particle radius scale), NOT `scale.L` (electrode thickness scale). This is confirmed by `PostProcessing.jl:110` which uses `case.param.scale.r0` for separation denormalization.

### 5.4 Memory Considerations

For long simulations (many cycles, many time steps), the snapshot vector can grow large. The function should:
- Write files incrementally (stream, not batch) if possible
- Or accept a `max_snapshots` parameter with LRU-style eviction
- For this initial design, we accumulate all snapshots and write at the end (simpler, acceptable for typical 1-50 cycle runs)

### 5.5 Node Coordinate Denormalization

Node coordinates in `case.mesh["thermal2D"]` are normalized. Must multiply by `scale.L` before writing to CSV.

## 6. Prerequisites and Backward Compatibility

### 6.1 Usage Prerequisites

- `solve_cycling` must be called with `save_detailed=true` for per-step data (element_currents, node_temperature, cohesive_damage, node_displacement)
- CZM must be enabled (`opt.czm_enabled = true`) for CZM-related CSV files
- Thermal must be enabled (`opt.thermal_enabled = true`) for temperature data
- `cohesive_driving_force.csv` requires `α_eff`, `β_n`, `β_p` parameters to be non-zero

### 6.2 Backward Compatibility

- `CZMSnapshot` and the `czm_snapshots` field in `CyclingResult` are only populated when `save_detailed=true`
- Existing code that does not call `export_cycling_csv` is unaffected
- The `CyclingResult` struct change (adding `czm_snapshots` field with default empty vector) is backward compatible

### 6.3 Graceful Degradation

When certain data is unavailable, the function should:
- Skip writing the corresponding CSV file (e.g., if CZM is disabled, skip cohesive_damage.csv)
- Print a warning to stdout listing which files were skipped and why
- Always write `cycle_summary.csv` if any cycle data is available

## 7. Testing Plan

1. Run `coupled_czm_thermal_example.jl` with `n_cycles=1`
2. Call `export_cycling_csv` on the result
3. Verify all 7 files are created with correct dimensions
4. Load each CSV in Python and verify:
   - No NaN/Inf values
   - Units are physical (T ~ 298-310 K, not 0-1)
   - Row counts match expected steps x entities
   - Damage values in [0,1] range
