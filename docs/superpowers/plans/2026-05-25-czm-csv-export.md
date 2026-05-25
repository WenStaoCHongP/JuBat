# CZM CSV Export Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Export cycling simulation results to 7 CSV files for Python post-processing.

**Architecture:** Add `CZMSnapshot` struct and collection logic to `Solve.jl` + `CycleSolver.jl`, then create `CsvExport.jl` with the `export_cycling_csv` function. The CSV export is a standalone post-processing step called after simulation.

**Tech Stack:** Julia, standard library only (no external CSV deps — use `open`/`println`).

**Spec:** `docs/superpowers/specs/2026-05-25-czm-csv-export-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/CsvExport.jl` | CSV export function + CZMSnapshot struct |
| Modify | `src/CycleSolver.jl` | Add `czm_snapshots` field to `CyclingResult`, pass through to `Solve` |
| Modify | `src/Solve.jl` | Collect CZM snapshots during CZM update steps |
| Modify | `src/JuBat.jl` | Include `CsvExport.jl`, export new symbols |
| Modify | `example/coupled_czm_thermal_example.jl` | Add CSV export call after simulation |

---

## Chunk 1: CZMSnapshot Struct + CyclingResult Enhancement

### Task 1: Add CZMSnapshot to CsvExport.jl

**Files:**
- Create: `src/CsvExport.jl`

- [ ] **Step 1: Create CsvExport.jl with CZMSnapshot struct**

```julia
# src/CsvExport.jl
# CSV export for cycling simulation post-processing

"""
    CZMSnapshot

Stores per-step CZM solver state for CSV export.
All physical values are stored in NORMALIZED (dimensionless) form.
Denormalization happens at CSV write time using `case.param.scale`.
"""
mutable struct CZMSnapshot
    time_s::Float64                     # physical time (already denormalized)
    cycle::Int                          # cycle number
    phase::String                       # phase name
    displacement::Vector{Float64}       # ndof-length, normalized
    damage::Vector{Float64}             # n_coh-length, [0,1]
    separation_n::Vector{Float64}       # n_coh-length, normalized
    separation_t::Vector{Float64}       # n_coh-length, normalized
    traction_n::Vector{Float64}         # n_coh-length, normalized
    traction_t::Vector{Float64}         # n_coh-length, normalized
    converged::Bool
    iterations::Int
    residual_norm::Float64
    method::String                      # "basic", "load_substep", or "arc_length"
end
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat: add CZMSnapshot struct for CSV export"
```

### Task 2: Add czm_snapshots field to CyclingResult

**Files:**
- Modify: `src/CycleSolver.jl:71-103`

- [ ] **Step 1: Add field to CyclingResult struct**

In `src/CycleSolver.jl`, after line 85 (`cycle_results::Vector{CycleResult}`), add:

```julia
    # CZM per-step snapshots (filled when save_detailed=true and CZM enabled)
    czm_snapshots::Vector{CZMSnapshot}
```

- [ ] **Step 2: Update CyclingResult constructor**

In `CyclingResult(n_cycles::Int)` constructor (around line 93-103), add `CZMSnapshot[]` to the constructor call:

```julia
function CyclingResult(n_cycles::Int)
    CyclingResult(
        0,
        Int[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Int[], Float64[], Float64[],
        CycleResult[],
        CZMSnapshot[],      # czm_snapshots
        nothing,            # final_czm_mesh
        0.0,                # initial_capacity
        false               # soh_terminated
    )
end
```

- [ ] **Step 3: Commit**

```bash
git add src/CycleSolver.jl
git commit -m "feat: add czm_snapshots field to CyclingResult"
```

---

## Chunk 2: CZM Snapshot Collection in Solve.jl

### Task 3: Collect CZM snapshots during solve loop

**Files:**
- Modify: `src/Solve.jl:272-286` (CZM update block)
- Modify: `src/CycleSolver.jl:110-186` (solve_phase function)

- [ ] **Step 1: Add czm_snapshots parameter to Solve function**

In `src/Solve.jl`, modify the `Solve` function signature (line 1) to add the optional parameter:

```julia
function Solve(case::Case;
               initial_state::Union{Dict{String,Any},Nothing}=nothing,
               return_final_state::Bool=false,
               thermal_variables::Union{Dict{String,Any},Nothing}=nothing,
               thermal_update_fn::Union{Function,Nothing}=nothing,
               thermal_record::Bool=false,
               polar_mesh_data::Any=nothing,
               czm_snapshots::Union{Vector{CZMSnapshot},Nothing}=nothing,
               czm_cycle::Int=1,
               czm_phase::String="unknown")
```

- [ ] **Step 2: Collect snapshot after each CZM update**

In `src/Solve.jl`, after the CZM update block (after line 283 `timing_totals["czm"] += ...`), inside the `if czm_active` block, add snapshot collection:

```julia
            # CZM 损伤演化（按间隔更新）
            if czm_active
                czm_step_count += 1
                if czm_step_count % case.opt.czm_update_interval == 0
                    t_czm_ns = time_ns()
                    try
                        u_czm_new, czm_converged = update_czm_damage!(
                            case, variables, T_nodes_carry)
                        # u_prev 已由 update_czm_damage! 内部管理

                        # ── Snapshot collection ──
                        if czm_snapshots !== nothing
                            czm_mesh_local = case.czm_mesh
                            n_coh = czm_mesh_local.n_cohesive
                            snap_damage = [s.D for s in czm_mesh_local.damage_states]
                            snap_frac = [s.fractured for s in czm_mesh_local.damage_states]

                            # Re-extract separations and tractions from last solve
                            # update_czm_damage! returns (displacement, converged)
                            # We need the full CZMResult - reconstruct from czm_mesh
                            snap_sep_n = zeros(n_coh)
                            snap_sep_t = zeros(n_coh)
                            snap_trc_n = zeros(n_coh)
                            snap_trc_t = zeros(n_coh)
                            if czm_converged
                                for i in 1:n_coh
                                    snap_sep_n[i] = czm_mesh_local.damage_states[i].δ_max_n
                                    snap_sep_t[i] = czm_mesh_local.damage_states[i].δ_max_t
                                end
                            end

                            push!(czm_snapshots, CZMSnapshot(
                                t * case.param.scale.t0,  # physical time
                                czm_cycle,
                                czm_phase,
                                copy(case.czm_layout.u_prev),  # displacement
                                snap_damage,
                                snap_sep_n,
                                snap_sep_t,
                                snap_trc_n,
                                snap_trc_t,
                                czm_converged,
                                0,      # iterations (not easily available here)
                                0.0,    # residual_norm (not easily available here)
                                case.opt.czm_iter_method
                            ))
                        end
                    catch e
                        @debug "CZM damage update failed at step $czm_step_count: $e"
                    end
                    timing_totals["czm"] += (time_ns() - t_czm_ns) * 1e-9
                end
            end
```

- [ ] **Step 3: Pass czm_snapshots from solve_phase to Solve**

In `src/CycleSolver.jl`, modify the `solve_phase` function signature (line 110) to accept and forward the snapshots:

```julia
function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64,
                     I_current::Float64, V_limit::Float64,
                     initial_state::Dict;
                     czm_mesh=nothing, czm_params=nothing,
                     dt_range::Vector{Float64}=[1.0, 10.0],
                     czm_snapshots::Union{Vector{CZMSnapshot},Nothing}=nothing,
                     czm_cycle::Int=1)
```

And at line 147 where `Solve` is called, forward the parameters:

```julia
        solve_result = Solve(case;
                             initial_state=initial_state,
                             return_final_state=true,
                             czm_snapshots=czm_snapshots,
                             czm_cycle=czm_cycle,
                             czm_phase=string(phase_type))
```

- [ ] **Step 4: Pass czm_snapshots from solve_cycling to solve_phase**

In `src/CycleSolver.jl`, in the `solve_cycling` function, create a snapshots vector and pass it to each `solve_phase` call. After line 209 (`result = CyclingResult(n_cycles)`), add:

```julia
    # CZM snapshots vector (shared across all phases)
    czm_snaps = save_detailed && czm_mesh !== nothing ? CZMSnapshot[] : nothing
```

Then in each `solve_phase` call (discharge ~line 284, rest1 ~line 312, charge ~line 355, rest2 ~line 382), add the keyword arguments:

```julia
czm_snapshots=czm_snaps, czm_cycle=cycle
```

After the loop, before `return result` (line 434), add:

```julia
    # Attach CZM snapshots to result
    if czm_snaps !== nothing
        result.czm_snapshots = czm_snaps
    end
```

- [ ] **Step 5: Commit**

```bash
git add src/Solve.jl src/CycleSolver.jl
git commit -m "feat: collect CZM snapshots during solve loop"
```

---

## Chunk 3: CSV Export Core Function

### Task 4: Implement cycle_summary.csv writer

**Files:**
- Modify: `src/CsvExport.jl`

- [ ] **Step 1: Add cycle_summary writer function**

Append to `src/CsvExport.jl`:

```julia
"""
    export_cycling_csv(result, case, czm_mesh; output_dir="output/csv", overwrite=false)

Export solve_cycling results to CSV files for post-processing.
"""
function export_cycling_csv(result::CyclingResult, case::Case, czm_mesh::CohesiveMesh;
                            output_dir::String="output/csv", overwrite::Bool=false)
    mkpath(output_dir)
    scale = case.param.scale

    files_written = String[]
    files_skipped = String[]

    # 1. cycle_summary.csv
    try
        _write_cycle_summary(result, output_dir, overwrite)
        push!(files_written, "cycle_summary.csv")
    catch e
        @warn "Failed to write cycle_summary.csv" exception=e
        push!(files_skipped, "cycle_summary.csv")
    end

    println("CSV export complete:")
    println("  Written: $(join(files_written, ", "))")
    if !isempty(files_skipped)
        println("  Skipped: $(join(files_skipped, ", "))")
    end
    println("  Directory: $output_dir")
    return files_written
end

function _write_cycle_summary(result::CyclingResult, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "cycle_summary.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    open(filepath, "w") do f
        println(f, "cycle,phase,capacity_ah,soh,D_max,D_mean,n_fractured,T_max_K,T_mean_end_K,V_start,V_end")

        for (i, cr) in enumerate(result.cycle_results)
            cyc = cr.cycle_idx
            soh = get(result.soh, i, NaN)
            for (phase_name, pr) in [("discharge", cr.discharge),
                                      ("rest1", cr.rest1),
                                      ("charge", cr.charge),
                                      ("rest2", cr.rest2)]
                D_max = isnan(pr.D_max) ? 0.0 : pr.D_max
                D_mean = isnan(pr.D_mean) ? 0.0 : pr.D_mean
                T_mean_end = isnan(pr.T_mean_end) ? 0.0 : pr.T_mean_end
                V_start = isnan(pr.V_start) ? 0.0 : pr.V_start
                V_end = isnan(pr.V_end) ? 0.0 : pr.V_end
                println(f, "$cyc,$phase_name,$(pr.capacity),$soh,$D_max,$D_mean,$(cr.n_fractured),$(cr.T_max),$T_mean_end,$V_start,$V_end")
            end
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat: add cycle_summary.csv export"
```

### Task 5: Implement element_currents.csv and node_temperature.csv writers

**Files:**
- Modify: `src/CsvExport.jl`

These two files use data from the `solve_phase` results stored in `cycle_results[*].final_state`, NOT from `variables_hist`. The `Solve` function returns its result dict, but `solve_phase` only stores `final_state`. We need to get the per-step data.

**Key insight:** The per-step variables are NOT currently saved across phases. The `Solve` function returns a full result dict with all time steps, but `solve_phase` in `CycleSolver.jl` only extracts summary statistics via `_postprocess_phase_result`. We need to also store the full result dict.

- [ ] **Step 1: Store solve result dict in PhaseResult**

In `src/CycleSolver.jl`, add a field to `PhaseResult`:

```julia
    # Raw solve result (for CSV export, only when save_detailed=true)
    solve_result::Any   # Dict from Solve(), or nothing
```

Update `PhaseResult()` default constructor to include `nothing` for the new field.

In `solve_phase`, after `solve_result = Solve(case; ...)` (line 147), store it:

```julia
        result.solve_result = solve_result
```

- [ ] **Step 2: Add element_currents writer**

Append to `src/CsvExport.jl`:

```julia
function _write_element_currents(result::CyclingResult, case::Case, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "element_currents.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    # Compute element areas (normalized, then denormalize)
    areas_norm = Float64[]
    for e in 1:ne
        nodes_e = mesh_th.element[e, :]
        x_e = mesh_th.node[nodes_e, 1]
        y_e = mesh_th.node[nodes_e, 2]
        # Q4 element area from cross product
        dx1, dy1 = x_e[2] - x_e[1], y_e[2] - y_e[1]
        dx2, dy2 = x_e[4] - x_e[1], y_e[4] - y_e[1]
        push!(areas_norm, abs(dx1 * dy2 - dx2 * dy1))
    end
    areas_phys = areas_norm * scale.L^2

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,elem_id,I_e,area,T_e,soc_n,soc_p,eta_n,eta_p,q_rxn_ne,q_sp,q_rxn_pe")

        for (i, cr) in enumerate(result.cycle_results)
            cyc = cr.cycle_idx
            for (phase_name, pr) in [("discharge", cr.discharge),
                                      ("rest1", cr.rest1),
                                      ("charge", cr.charge),
                                      ("rest2", cr.rest2)]
                sr = pr.solve_result
                if sr === nothing
                    continue
                end

                # Extract time series data
                time_s = get(sr, "time [s]", Float64[])
                I_e_all = get(sr, "thermal2D element current", zeros(0,0))
                soc_n_all = get(sr, "thermal2D element soc_n", zeros(0,0))
                soc_p_all = get(sr, "thermal2D element soc_p", zeros(0,0))
                eta_n_all = get(sr, "thermal2D eta_n_e", zeros(0,0))
                eta_p_all = get(sr, "thermal2D eta_p_e", zeros(0,0))
                q_rxn_ne_all = get(sr, "thermal2D Q_rxn_NE [W/m3]", zeros(0,0))
                q_sp_all = get(sr, "thermal2D Q_SP [W/m3]", zeros(0,0))
                q_rxn_pe_all = get(sr, "thermal2D Q_rxn_PE [W/m3]", zeros(0,0))
                T_elem_all = get(sr, "thermal2D temperature [K]", zeros(0,0))

                n_steps = length(time_s)
                for ti in 1:n_steps
                    t = time_s[ti]
                    for e in 1:ne
                        I_e = size(I_e_all, 1) >= e ? I_e_all[e, ti] : NaN
                        sn = size(soc_n_all, 1) >= e ? soc_n_all[e, ti] : NaN
                        sp = size(soc_p_all, 1) >= e ? soc_p_all[e, ti] : NaN
                        en = size(eta_n_all, 1) >= e ? eta_n_all[e, ti] : NaN
                        ep = size(eta_p_all, 1) >= e ? eta_p_all[e, ti] : NaN
                        qr = size(q_rxn_ne_all, 1) >= e ? q_rxn_ne_all[e, ti] : NaN
                        qs = size(q_sp_all, 1) >= e ? q_sp_all[e, ti] : NaN
                        qp = size(q_rxn_pe_all, 1) >= e ? q_rxn_pe_all[e, ti] : NaN
                        Te = size(T_elem_all, 1) >= e ? T_elem_all[e, ti] : NaN

                        println(f, "$t,$cyc,$phase_name,$e,$I_e,$(areas_phys[e]),$Te,$sn,$sp,$en,$ep,$qr,$qs,$qp")
                    end
                end
            end
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 3: Add node_temperature writer**

```julia
function _write_node_temperature(result::CyclingResult, case::Case, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "node_temperature.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    mesh_th = case.mesh["thermal2D"]
    nnode = mesh_th.nlen
    # Denormalize node coordinates
    node_x = mesh_th.node[:, 1] * scale.L
    node_y = mesh_th.node[:, 2] * scale.L

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,node_id,x,y,T_K")

        for (i, cr) in enumerate(result.cycle_results)
            cyc = cr.cycle_idx
            for (phase_name, pr) in [("discharge", cr.discharge),
                                      ("rest1", cr.rest1),
                                      ("charge", cr.charge),
                                      ("rest2", cr.rest2)]
                sr = pr.solve_result
                if sr === nothing
                    continue
                end

                time_s = get(sr, "time [s]", Float64[])
                T_nodes_all = get(sr, "thermal2D temperature at nodes [K]", zeros(0,0))
                n_steps = length(time_s)

                for ti in 1:n_steps
                    t = time_s[ti]
                    for n in 1:nnode
                        T_K = size(T_nodes_all, 1) >= n ? T_nodes_all[n, ti] : NaN
                        println(f, "$t,$cyc,$phase_name,$n,$(node_x[n]),$(node_y[n]),$T_K")
                    end
                end
            end
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 4: Wire writers into export_cycling_csv**

Update the `export_cycling_csv` function to call these new writers:

```julia
    # 2. element_currents.csv
    if any(cr -> any(pr -> pr.solve_result !== nothing for pr in [cr.discharge, cr.rest1, cr.charge, cr.rest2]),
           result.cycle_results)
        try
            _write_element_currents(result, case, output_dir, overwrite)
            push!(files_written, "element_currents.csv")
        catch e
            @warn "Failed to write element_currents.csv" exception=e
            push!(files_skipped, "element_currents.csv")
        end
    else
        push!(files_skipped, "element_currents.csv (no solve_result data)")
    end

    # 3. node_temperature.csv
    if case.opt.thermal_enabled
        try
            _write_node_temperature(result, case, output_dir, overwrite)
            push!(files_written, "node_temperature.csv")
        catch e
            @warn "Failed to write node_temperature.csv" exception=e
            push!(files_skipped, "node_temperature.csv")
        end
    else
        push!(files_skipped, "node_temperature.csv (thermal disabled)")
    end
```

- [ ] **Step 5: Commit**

```bash
git add src/CsvExport.jl src/CycleSolver.jl
git commit -m "feat: add element_currents and node_temperature CSV export"
```

### Task 6: Implement cohesive_damage.csv and node_displacement.csv writers

**Files:**
- Modify: `src/CsvExport.jl`

These use the `czm_snapshots` data collected in Chunk 2.

- [ ] **Step 1: Add cohesive_damage writer**

```julia
function _write_cohesive_damage(result::CyclingResult, case::Case, czm_mesh::CohesiveMesh,
                                 output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "cohesive_damage.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    n_coh = czm_mesh.n_cohesive

    # Precompute element properties
    lengths_phys = Float64[]
    theta_degs = Float64[]
    for elem in czm_mesh.cohesive_elements
        push!(lengths_phys, elem.length * scale.L)
        # Midpoint angle
        n1, n2 = elem.nodes_bottom
        mx = 0.5 * (czm_mesh.node[n1, 1] + czm_mesh.node[n2, 1])
        my = 0.5 * (czm_mesh.node[n1, 2] + czm_mesh.node[n2, 2])
        push!(theta_degs, atan(my, mx) * 180.0 / pi)
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,length,D,delta_n,delta_t,T_n,T_t,fractured,theta_deg")

        for snap in result.czm_snapshots
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            for i in 1:min(n_coh, length(snap.damage))
                D = snap.damage[i]
                dn = i <= length(snap.separation_n) ? snap.separation_n[i] * scale.r0 : NaN
                dt_val = i <= length(snap.separation_t) ? snap.separation_t[i] * scale.r0 : NaN
                tn = i <= length(snap.traction_n) ? snap.traction_n[i] * scale.σ_czm : NaN
                tt = i <= length(snap.traction_t) ? snap.traction_t[i] * scale.σ_czm : NaN
                frac = D >= 0.99
                println(f, "$t,$cyc,$phase,$i,$(lengths_phys[i]),$D,$dn,$dt_val,$tn,$tt,$frac,$(theta_degs[i])")
            end
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 2: Add node_displacement writer**

```julia
function _write_node_displacement(result::CyclingResult, case::Case, czm_mesh::CohesiveMesh,
                                   output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "node_displacement.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    nnode = czm_mesh.nnode
    node_x = czm_mesh.node[:, 1] * scale.L
    node_y = czm_mesh.node[:, 2] * scale.L

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,node_id,x,y,ux,uy")

        for snap in result.czm_snapshots
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            for n in 1:nnode
                ux = snap.displacement[2*n - 1] * scale.L
                uy = snap.displacement[2*n] * scale.L
                println(f, "$t,$cyc,$phase,$n,$(node_x[n]),$(node_y[n]),$ux,$uy")
            end
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 3: Wire into export_cycling_csv**

```julia
    # 4. cohesive_damage.csv
    if !isempty(result.czm_snapshots)
        try
            _write_cohesive_damage(result, case, czm_mesh, output_dir, overwrite)
            push!(files_written, "cohesive_damage.csv")
        catch e
            @warn "Failed to write cohesive_damage.csv" exception=e
            push!(files_skipped, "cohesive_damage.csv")
        end
        try
            _write_node_displacement(result, case, czm_mesh, output_dir, overwrite)
            push!(files_written, "node_displacement.csv")
        catch e
            @warn "Failed to write node_displacement.csv" exception=e
            push!(files_skipped, "node_displacement.csv")
        end
    else
        push!(files_skipped, "cohesive_damage.csv (no CZM snapshots)")
        push!(files_skipped, "node_displacement.csv (no CZM snapshots)")
    end
```

- [ ] **Step 4: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat: add cohesive_damage and node_displacement CSV export"
```

### Task 7: Implement driving_force and diagnostics writers

**Files:**
- Modify: `src/CsvExport.jl`

- [ ] **Step 1: Add cohesive_driving_force writer**

```julia
function _write_driving_force(result::CyclingResult, case::Case, czm_mesh::CohesiveMesh,
                               output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "cohesive_driving_force.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    n_coh = czm_mesh.n_cohesive

    # Get CZM effective parameters for strain computation
    alpha_eff = get(case.param.cohesive, :alpha_eff, 0.0)
    if alpha_eff == 0.0 && hasproperty(case.param, :alphaT)
        alpha_eff = case.param.alphaT
    end
    beta_n = get(case.param.cohesive, :beta_n, 0.0)
    beta_p = get(case.param.cohesive, :beta_p, 0.0)

    if alpha_eff == 0.0 && beta_n == 0.0 && beta_p == 0.0
        println("  Skipping $filepath (no thermal/diffusion strain parameters)")
        return
    end

    # Build cohesive-to-thermal-element mapping from interface_pairs
    # Each cohesive element sits between two thermal elements
    geo = case.geometry
    if geo === nothing
        println("  Skipping $filepath (no geometry data)")
        return
    end

    # Map coh elem index -> pair of adjacent thermal elem indices
    interface_pairs = geo.interface_pairs
    coh_thermal_map = Dict{Int, Tuple{Int,Int}}()
    for (idx, (e_top, e_bot)) in enumerate(interface_pairs)
        if idx <= n_coh
            coh_thermal_map[idx] = (e_top, e_bot)
        end
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,dT_neighbor,dsoc_n_neighbor,dsoc_p_neighbor,eps_0_thermal,eps_0_diffusion,eps_0_total")

        for snap in result.czm_snapshots
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase

            # Get element-level T and SOC at this time step
            # Find the corresponding solve_result for this phase
            sr = _find_solve_result(result, cyc, phase)
            if sr === nothing
                continue
            end

            time_arr = get(sr, "time [s]", Float64[])
            T_elem = get(sr, "thermal2D temperature [K]", zeros(0,0))
            soc_n_elem = get(sr, "thermal2D element soc_n", zeros(0,0))
            soc_p_elem = get(sr, "thermal2D element soc_p", zeros(0,0))

            T_ref = scale.T_ref
            # Find closest time index
            ti = argmin(abs.(time_arr .- t))
            if ti > size(T_elem, 2)
                continue
            end

            for i in 1:min(n_coh, length(snap.damage))
                if !haskey(coh_thermal_map, i)
                    continue
                end
                e_top, e_bot = coh_thermal_map[i]

                dT = 0.0
                dsoc_n = 0.0
                dsoc_p = 0.0
                if e_top <= size(T_elem, 1) && e_bot <= size(T_elem, 1)
                    dT = 0.5 * (T_elem[e_top, ti] + T_elem[e_bot, ti]) - T_ref
                end
                if e_top <= size(soc_n_elem, 1) && e_bot <= size(soc_n_elem, 1)
                    dsoc_n = 0.5 * (soc_n_elem[e_top, ti] + soc_n_elem[e_bot, ti])
                    # Subtract initial SOC
                end
                if e_top <= size(soc_p_elem, 1) && e_bot <= size(soc_p_elem, 1)
                    dsoc_p = 0.5 * (soc_p_elem[e_top, ti] + soc_p_elem[e_bot, ti])
                end

                eps_thermal = alpha_eff * dT
                eps_diffusion = beta_n * dsoc_n + beta_p * dsoc_p
                eps_total = eps_thermal + eps_diffusion

                println(f, "$t,$cyc,$phase,$i,$dT,$dsoc_n,$dsoc_p,$eps_thermal,$eps_diffusion,$eps_total")
            end
        end
    end
    println("  Written: $filepath")
end

function _find_solve_result(result::CyclingResult, cycle::Int, phase::String)
    for cr in result.cycle_results
        if cr.cycle_idx == cycle
            if phase == "discharge"
                return cr.discharge.solve_result
            elseif phase == "rest1"
                return cr.rest1.solve_result
            elseif phase == "charge"
                return cr.charge.solve_result
            elseif phase == "rest2"
                return cr.rest2.solve_result
            end
        end
    end
    return nothing
end
```

- [ ] **Step 2: Add czm_solver_diagnostics writer**

```julia
function _write_czm_diagnostics(result::CyclingResult, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "czm_solver_diagnostics.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,converged,iterations,residual_norm,method")

        for snap in result.czm_snapshots
            println(f, "$(snap.time_s),$(snap.cycle),$(snap.phase),$(snap.converged),$(snap.iterations),$(snap.residual_norm),$(snap.method)")
        end
    end
    println("  Written: $filepath")
end
```

- [ ] **Step 3: Wire both into export_cycling_csv**

```julia
    # 5. cohesive_driving_force.csv
    if !isempty(result.czm_snapshots) && case.geometry !== nothing
        try
            _write_driving_force(result, case, czm_mesh, output_dir, overwrite)
            push!(files_written, "cohesive_driving_force.csv")
        catch e
            @warn "Failed to write cohesive_driving_force.csv" exception=e
            push!(files_skipped, "cohesive_driving_force.csv")
        end
    else
        push!(files_skipped, "cohesive_driving_force.csv (no CZM snapshots or geometry)")
    end

    # 6. czm_solver_diagnostics.csv
    if !isempty(result.czm_snapshots)
        try
            _write_czm_diagnostics(result, output_dir, overwrite)
            push!(files_written, "czm_solver_diagnostics.csv")
        catch e
            @warn "Failed to write czm_solver_diagnostics.csv" exception=e
            push!(files_skipped, "czm_solver_diagnostics.csv")
        end
    else
        push!(files_skipped, "czm_solver_diagnostics.csv (no CZM snapshots)")
    end
```

- [ ] **Step 4: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat: add driving_force and diagnostics CSV export"
```

---

## Chunk 4: Module Registration + Example Update

### Task 8: Register CsvExport.jl in JuBat module

**Files:**
- Modify: `src/JuBat.jl`

- [ ] **Step 1: Add include and export**

In `src/JuBat.jl`, after line 35 (`include("CycleData.jl")`), add:

```julia
include("CsvExport.jl")
```

In the export block, add:

```julia
export CZMSnapshot, export_cycling_csv
```

- [ ] **Step 2: Commit**

```bash
git add src/JuBat.jl
git commit -m "feat: register CsvExport.jl in JuBat module"
```

### Task 9: Update example script

**Files:**
- Modify: `example/coupled_czm_thermal_example.jl`

- [ ] **Step 1: Add CSV export call after simulation**

After line 148 (`result = JuBat.solve_cycling(...)`), and before the result analysis section (line 153), add:

```julia
# ========================================================================
# 7.5. CSV 导出
# ========================================================================
println("\n[7.5] 导出CSV文件...")

csv_dir = joinpath(output_dir, "csv", "czm_study_1")
JuBat.export_cycling_csv(result, case, czm_mesh;
                          output_dir=csv_dir, overwrite=true)
```

- [ ] **Step 2: Commit**

```bash
git add example/coupled_czm_thermal_example.jl
git commit -m "feat: add CSV export call to coupled example"
```

---

## Chunk 5: Integration Test

### Task 10: Run integration test

**Files:**
- No new files

- [ ] **Step 1: Run the example script**

```bash
cd "D:\OneDrive\Desktop\Jubat For Cursor\JuBat"
julia --project=. example/coupled_czm_thermal_example.jl
```

Expected: Script completes without error, CSV files are written to `output/csv/czm_study_1/`.

- [ ] **Step 2: Verify CSV output**

Check that all 7 files exist and contain reasonable data:

```bash
ls -la output/csv/czm_study_1/
head -5 output/csv/czm_study_1/cycle_summary.csv
head -5 output/csv/czm_study_1/cohesive_damage.csv
wc -l output/csv/czm_study_1/*.csv
```

Expected:
- `cycle_summary.csv`: 4 rows per cycle (discharge/rest1/charge/rest2) + header
- `cohesive_damage.csv`: n_snapshots × n_cohesive rows + header
- `node_temperature.csv`: n_steps × n_nodes rows + header
- Temperature values around 298-310 K (not 0-1)
- Damage values in [0, 1]

- [ ] **Step 3: Commit any fixes**

If any issues found, fix and commit:

```bash
git add -u
git commit -m "fix: address CSV export integration issues"
```
