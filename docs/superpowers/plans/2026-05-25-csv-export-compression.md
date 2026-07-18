# CSV Export Compression Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add configurable sampling to `export_cycling_csv` so users can reduce CSV output from ~1.86 GB to ~20-30 MB.

**Architecture:** Introduce a `CsvExportOptions` struct that controls which time steps and CZM snapshots are written. Existing write functions gain a `csv_opt` parameter and call helper predicates to skip non-essential rows. Default mode is `:full` for backward compatibility.

**Tech Stack:** Julia, existing `src/CsvExport.jl` module.

**Spec:** `docs/superpowers/specs/2026-05-25-csv-export-compression-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `src/CsvExport.jl` | Modify | Add `CsvExportOptions` struct, helpers, modify all write functions |
| `src/JuBat.jl` | Modify | Add `CsvExportOptions` to export list (line 83) |
| `example/coupled_czm_thermal_example.jl` | Modify | Add `csv_opt` usage demonstration |

No new files needed.

---

## Chunk 1: Struct, Helpers, and Export

### Task 1: Add `CsvExportOptions` struct, helpers, and module export

**Files:**
- Modify: `src/CsvExport.jl` (insert after line 3, before the main export function)
- Modify: `src/JuBat.jl` line 83 (add export)

- [ ] **Step 1: Add `CsvExportOptions` struct and constructors to `src/CsvExport.jl`**

Insert after line 3 (after the module comment block):

```julia
# ========================================================================
# Export options (configurable sampling)
# ========================================================================

"""
    CsvExportOptions

Controls which time steps and snapshots are written to CSV files.

# Fields
- `mode::Symbol`: `:full` (all steps), `:phase_ends` (first+last per phase), `:custom` (every N steps)
- `save_every::Int`: Step interval for `:custom` mode (ignored for other modes)
- `full_output_cycles::Vector{Int}`: Cycles that always get full output regardless of mode
- `skip_files::Vector{String}`: Filenames to skip entirely (e.g. `["node_displacement.csv"]`)

# Examples
```julia
# Default: full output (backward compatible)
CsvExportOptions()

# Only first+last step per phase, with full output for cycles 1 and 10
CsvExportOptions(:phase_ends, 1, [1, 10], String[])

# Every 50th step, skip displacement file
CsvExportOptions(:custom, 50, Int[], ["node_displacement.csv"])
```
"""
struct CsvExportOptions
    mode::Symbol
    save_every::Int
    full_output_cycles::Vector{Int}
    skip_files::Vector{String}
end

CsvExportOptions() = CsvExportOptions(:full, 1, [1], String[])
CsvExportOptions(mode::Symbol) = CsvExportOptions(mode, 1, [1], String[])
```

- [ ] **Step 2: Add validation function**

Insert after the constructors:

```julia
function _validate_csv_opt(csv_opt::CsvExportOptions)
    if csv_opt.mode ∉ (:full, :phase_ends, :custom)
        error("CsvExportOptions: invalid mode $(csv_opt.mode), must be :full, :phase_ends, or :custom")
    end
    if csv_opt.mode == :custom && csv_opt.save_every < 1
        error("CsvExportOptions: save_every must be ≥ 1 for :custom mode")
    end
end
```

- [ ] **Step 3: Add `_should_output_step` helper**

```julia
"""Whether time step `ti` (of `n_steps` total) in `cycle` should be written."""
function _should_output_step(csv_opt::CsvExportOptions, cycle::Int, ti::Int, n_steps::Int)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return (ti == 1 || ti == n_steps)
    csv_opt.mode == :custom && return (ti == 1 || ti == n_steps || ti % csv_opt.save_every == 1)
    return true
end
```

- [ ] **Step 4: Add `_compute_last_snap_keys` and `_should_output_snapshot` helpers**

```julia
"""Pre-compute Set of (cycle, phase) keys that correspond to the last snapshot in each phase group.
Returns empty Set for empty snapshot arrays."""
function _compute_last_snap_keys(czm_snapshots)
    last_keys = Set{Tuple{Int,String}}()
    isempty(czm_snapshots) && return last_keys
    prev_key = (czm_snapshots[1].cycle, czm_snapshots[1].phase)
    for i in 2:length(czm_snapshots)
        key = (czm_snapshots[i].cycle, czm_snapshots[i].phase)
        if key != prev_key
            push!(last_keys, prev_key)
        end
        prev_key = key
    end
    push!(last_keys, prev_key)  # always include the final group
    return last_keys
end

"""Whether a CZM snapshot should be written, given it is (or isn't) the last in its phase group."""
function _should_output_snapshot(csv_opt::CsvExportOptions, cycle::Int, is_last_in_phase::Bool)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return is_last_in_phase
    csv_opt.mode == :custom && return true
    return true
end
```

- [ ] **Step 5: Add export to `src/JuBat.jl`**

On line 83, change:
```julia
export CZMSnapshot, export_cycling_csv
```
to:
```julia
export CZMSnapshot, export_cycling_csv, CsvExportOptions
```

- [ ] **Step 6: Commit**

```bash
git add src/CsvExport.jl src/JuBat.jl
git commit -m "feat(csv): add CsvExportOptions struct, helpers, and module export"
```

---

## Chunk 2: Modify Main Export Function and Skip-Files

### Task 2: Update `export_cycling_csv` signature and add skip-files logic

**Files:**
- Modify: `src/CsvExport.jl` — `export_cycling_csv` function

- [ ] **Step 1: Update function signature and docstring**

Replace the existing docstring and function signature (lines 9-16) with:

```julia
"""
    export_cycling_csv(result, case, czm_mesh;
                       output_dir="output/csv", overwrite=false, csv_opt=CsvExportOptions())

Export solve_cycling results to CSV files for post-processing.

# Arguments
- `csv_opt::CsvExportOptions`: Controls sampling mode and which files to write.
  Default is `CsvExportOptions()` (`:full` mode, backward compatible).
"""
function export_cycling_csv(result, case, czm_mesh;
                            output_dir::String="output/csv", overwrite::Bool=false,
                            csv_opt::CsvExportOptions=CsvExportOptions())
    _validate_csv_opt(csv_opt)
    mkpath(output_dir)

    files_written = String[]
    files_skipped = String[]
```

- [ ] **Step 2: Add skip-files guards and pass `csv_opt` to all write calls**

Replace the body (lines 22-104) with the following structure. Each file block gets:
1. A `skip_files` guard
2. `csv_opt` passed to the `_write_*` call

```julia
    # 1. cycle_summary.csv — always written, no sampling
    try
        _write_cycle_summary(result, output_dir, overwrite)
        push!(files_written, "cycle_summary.csv")
    catch e
        @warn "Failed to write cycle_summary.csv" exception=e
        push!(files_skipped, "cycle_summary.csv")
    end

    # 2. element_currents.csv
    if "element_currents.csv" in csv_opt.skip_files
        push!(files_skipped, "element_currents.csv (skipped by option)")
    else
        has_data = !isempty(result.cycle_results) &&
                   any(cr -> cr.discharge.solve_result !== nothing, result.cycle_results)
        if has_data
            try
                _write_element_currents(result, case, output_dir, overwrite, csv_opt)
                push!(files_written, "element_currents.csv")
            catch e
                @warn "Failed to write element_currents.csv" exception=e
                push!(files_skipped, "element_currents.csv")
            end
        else
            push!(files_skipped, "element_currents.csv (no solve_result data)")
        end
    end

    # 3. node_temperature.csv
    if "node_temperature.csv" in csv_opt.skip_files
        push!(files_skipped, "node_temperature.csv (skipped by option)")
    else
        if case.opt.thermal_enabled
            try
                _write_node_temperature(result, case, output_dir, overwrite, csv_opt)
                push!(files_written, "node_temperature.csv")
            catch e
                @warn "Failed to write node_temperature.csv" exception=e
                push!(files_skipped, "node_temperature.csv")
            end
        else
            push!(files_skipped, "node_temperature.csv (thermal disabled)")
        end
    end

    # 4. cohesive_damage.csv
    if "cohesive_damage.csv" in csv_opt.skip_files
        push!(files_skipped, "cohesive_damage.csv (skipped by option)")
    elseif !isempty(result.czm_snapshots)
        try
            _write_cohesive_damage(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            push!(files_written, "cohesive_damage.csv")
        catch e
            @warn "Failed to write cohesive_damage.csv" exception=e
            push!(files_skipped, "cohesive_damage.csv")
        end
    else
        push!(files_skipped, "cohesive_damage.csv (no CZM snapshots)")
    end

    # 5. node_displacement.csv
    if "node_displacement.csv" in csv_opt.skip_files
        push!(files_skipped, "node_displacement.csv (skipped by option)")
    elseif !isempty(result.czm_snapshots)
        try
            _write_node_displacement(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            push!(files_written, "node_displacement.csv")
        catch e
            @warn "Failed to write node_displacement.csv" exception=e
            push!(files_skipped, "node_displacement.csv")
        end
    else
        push!(files_skipped, "node_displacement.csv (no CZM snapshots)")
    end

    # 6. cohesive_driving_force.csv
    if "cohesive_driving_force.csv" in csv_opt.skip_files
        push!(files_skipped, "cohesive_driving_force.csv (skipped by option)")
    elseif !isempty(result.czm_snapshots) && case.geometry !== nothing
        try
            _write_driving_force(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            push!(files_written, "cohesive_driving_force.csv")
        catch e
            @warn "Failed to write cohesive_driving_force.csv" exception=e
            push!(files_skipped, "cohesive_driving_force.csv")
        end
    else
        push!(files_skipped, "cohesive_driving_force.csv (no CZM snapshots or geometry)")
    end

    # 7. czm_solver_diagnostics.csv — always written (small), no sampling
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

The println/reporting block at the end remains unchanged.

- [ ] **Step 3: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add csv_opt parameter and skip_files support to export_cycling_csv"
```

---

## Chunk 3: Modify Solve-Result-Driven Write Functions

### Task 3: Update `_write_element_currents` with step filtering

**Files:**
- Modify: `src/CsvExport.jl` — `_write_element_currents` function

- [ ] **Step 1: Add `csv_opt` parameter and step filter**

Change the signature to:
```julia
function _write_element_currents(result, case, output_dir::String, overwrite::Bool,
                                  csv_opt::CsvExportOptions)
```

Inside the `for ti in 1:n_steps` loop, add the filter as the first check:

```julia
                n_steps = length(time_s)
                for ti in 1:n_steps
                    if !_should_output_step(csv_opt, cyc, ti, n_steps)
                        continue
                    end
                    t = time_s[ti]
                    for e in 1:ne
                        # ... (existing element loop body unchanged)
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add step sampling to element_currents export"
```

### Task 4: Update `_write_node_temperature` with step filtering

**Files:**
- Modify: `src/CsvExport.jl` — `_write_node_temperature` function

- [ ] **Step 1: Add `csv_opt` parameter and step filter**

Change the signature to:
```julia
function _write_node_temperature(result, case, output_dir::String, overwrite::Bool,
                                  csv_opt::CsvExportOptions)
```

Inside the `for ti in 1:n_steps` loop, add the filter:

```julia
                n_steps = length(time_s)

                for ti in 1:n_steps
                    if !_should_output_step(csv_opt, cyc, ti, n_steps)
                        continue
                    end
                    t = time_s[ti]
                    for n in 1:nnode
                        # ... (existing node loop body unchanged)
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add step sampling to node_temperature export"
```

---

## Chunk 4: Modify Snapshot-Driven Write Functions

### Task 5: Update `_write_cohesive_damage` with snapshot filtering

**Files:**
- Modify: `src/CsvExport.jl` — `_write_cohesive_damage` function

- [ ] **Step 1: Add `csv_opt` parameter and snapshot filter**

Change the signature to:
```julia
function _write_cohesive_damage(result, case, czm_mesh,
                                 output_dir::String, overwrite::Bool,
                                 csv_opt::CsvExportOptions)
```

Before the snapshot loop, compute phase boundaries:
```julia
    last_snap_keys = _compute_last_snap_keys(result.czm_snapshots)
```

Inside the snapshot loop, add the filter:
```julia
        for snap in result.czm_snapshots
            key = (snap.cycle, snap.phase)
            is_last = key in last_snap_keys
            if !_should_output_snapshot(csv_opt, snap.cycle, is_last)
                continue
            end
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            for i in 1:min(n_coh, length(snap.damage))
                # ... (existing cohesive element loop body unchanged)
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add snapshot sampling to cohesive_damage export"
```

### Task 6: Update `_write_node_displacement` with snapshot filtering

**Files:**
- Modify: `src/CsvExport.jl` — `_write_node_displacement` function

- [ ] **Step 1: Add `csv_opt` parameter and snapshot filter**

Change the signature to:
```julia
function _write_node_displacement(result, case, czm_mesh,
                                   output_dir::String, overwrite::Bool,
                                   csv_opt::CsvExportOptions)
```

Before the snapshot loop:
```julia
    last_snap_keys = _compute_last_snap_keys(result.czm_snapshots)
```

Inside the snapshot loop:
```julia
        for snap in result.czm_snapshots
            key = (snap.cycle, snap.phase)
            is_last = key in last_snap_keys
            if !_should_output_snapshot(csv_opt, snap.cycle, is_last)
                continue
            end
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            nd = div(length(snap.displacement), 2)
            for n in 1:min(nnode, nd)
                # ... (existing node loop body unchanged)
```

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add snapshot sampling to node_displacement export"
```

### Task 7: Update `_write_driving_force` with snapshot filtering

**Files:**
- Modify: `src/CsvExport.jl` — `_write_driving_force` function

- [ ] **Step 1: Add `csv_opt` parameter and snapshot filter**

Change the signature to:
```julia
function _write_driving_force(result, case, czm_mesh,
                               output_dir::String, overwrite::Bool,
                               csv_opt::CsvExportOptions)
```

Same pattern as Tasks 5-6: compute `last_snap_keys` before the loop, add `_should_output_snapshot` check inside the loop.

Note: This function is currently a no-op (early return because `alpha_eff == 0.0`), so the filter won't have visible effect until thermal/diffusion strain parameters are integrated. Still add it for consistency.

- [ ] **Step 2: Commit**

```bash
git add src/CsvExport.jl
git commit -m "feat(csv): add snapshot sampling to driving_force export"
```

---

## Chunk 5: Update Example and Verify

### Task 8: Update example script to demonstrate `csv_opt`

**Files:**
- Modify: `example/coupled_czm_thermal_example.jl` lines 155-161

- [ ] **Step 1: Replace the CSV export call**

Replace lines 155-161 with:

```julia
# ========================================================================
# 7.5. CSV 导出
# ========================================================================
println("\n[7.5] 导出CSV文件...")

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)
csv_dir = joinpath(output_dir, "csv", "czm_study_1")

# 配置CSV导出选项：仅输出每个阶段首尾步，指定循环完整输出
csv_opt = JuBat.CsvExportOptions(
    :phase_ends,                # 仅阶段首尾
    1,                          # save_every (仅 :custom 模式使用)
    [1, result.n_cycles],       # 第1个和最后一个循环完整输出
    String[]                    # 不跳过任何文件
)

JuBat.export_cycling_csv(result, case, czm_mesh;
                          output_dir=csv_dir, overwrite=true,
                          csv_opt=csv_opt)
```

- [ ] **Step 2: Commit**

```bash
git add example/coupled_czm_thermal_example.jl
git commit -m "feat(csv): demonstrate CsvExportOptions in example script"
```

### Task 9: Verify backward compatibility

**Files:** None (manual verification)

- [ ] **Step 1: Verify default behavior is unchanged**

Confirm that calling `export_cycling_csv(result, case, czm_mesh)` without `csv_opt` still produces the same full output. This is guaranteed because:
- Default `CsvExportOptions()` uses `:full` mode
- All `_should_output_step` calls return `true` in `:full` mode
- All `_should_output_snapshot` calls return `true` in `:full` mode

- [ ] **Step 2: Verify compression mode works**

Run the example script and check that:
- CSV files are created
- `element_currents.csv` has fewer rows than before (only first+last step per phase, plus full output for cycles 1 and N)
- `cycle_summary.csv` is unchanged
- `czm_solver_diagnostics.csv` is unchanged

---

## Summary

| Task | Description | Est. Lines Changed |
|------|-------------|--------------------|
| 1 | Add struct + helpers + export | +70 lines |
| 2 | Update main function + skip_files | ~90 lines modified |
| 3 | Filter `_write_element_currents` | ~5 lines modified |
| 4 | Filter `_write_node_temperature` | ~5 lines modified |
| 5 | Filter `_write_cohesive_damage` | ~7 lines modified |
| 6 | Filter `_write_node_displacement` | ~7 lines modified |
| 7 | Filter `_write_driving_force` | ~7 lines modified |
| 8 | Update example | ~10 lines modified |
| 9 | Verification | 0 lines |
