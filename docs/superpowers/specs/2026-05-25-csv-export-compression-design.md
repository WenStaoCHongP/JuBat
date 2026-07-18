# CSV Export Compression Design

**Date**: 2026-05-25
**Status**: Approved
**Scope**: `src/CsvExport.jl` + `example/coupled_czm_thermal_example.jl`

## Problem

Running `coupled_czm_thermal_example.jl` produces ~1.86 GB of CSV files:

| File | Size |
|------|------|
| node_displacement.csv | 594 MB |
| element_currents.csv | 541 MB |
| node_temperature.csv | 447 MB |
| cohesive_damage.csv | 279 MB |
| czm_solver_diagnostics.csv | 280 KB |
| cycle_summary.csv | 5 KB |

Root cause: every time step × every element/node is written in full. With ~500 steps/cycle × 10 cycles and ~80+ elements/nodes, this produces millions of rows.

## Design

### 1. New `CsvExportOptions` struct

```julia
struct CsvExportOptions
    mode::Symbol                        # :full, :phase_ends, :custom
    save_every::Int                     # step interval for :custom mode
    full_output_cycles::Vector{Int}     # cycles with complete output
    skip_files::Vector{String}          # files to skip entirely
end

# Convenience constructors
# Default is :full for backward compatibility
CsvExportOptions() = CsvExportOptions(:full, 1, [1], String[])
CsvExportOptions(mode::Symbol) = CsvExportOptions(mode, 1, [1], String[])
```

### 2. Modified function signature

```julia
function export_cycling_csv(result, case, czm_mesh;
                            output_dir="output/csv", overwrite=false,
                            csv_opt=CsvExportOptions())
```

Backward-compatible: default `CsvExportOptions()` uses `:full` mode (original behavior). Users opt into compression by explicitly passing `mode=:phase_ends` or `mode=:custom`.

### 3. Sampling modes

| Mode | Behavior |
|------|----------|
| `:full` | All time steps output (default, original behavior) |
| `:phase_ends` | Only first and last step per phase |
| `:custom` | Every `save_every`-th step, plus always first and last step |

For cycles in `full_output_cycles`, the mode is overridden to `:full`.

### 4. Per-file logic

**element_currents.csv, node_temperature.csv** (solve_result driven):
- Filter time indices via `_should_output_step(csv_opt, cycle, ti, n_steps)`
- No change to column layout
- Caller still checks `solve_result === nothing` before applying filter

**cohesive_damage.csv, node_displacement.csv** (snapshot driven):
- `:phase_ends` mode: keep only the last snapshot per (cycle, phase) combination
- Phase boundary detection: pre-compute a `Set{Tuple{Int,String}}` of last-snapshot keys from the sorted `czm_snapshots` array
- Phase name note: CZM snapshots use `string(phase_type)` which produces `"PHASE_DISCHARGE"`, `"PHASE_CHARGE"`, `"PHASE_REST"` — the existing code already handles this correctly in the CSV output column

**czm_solver_diagnostics.csv, cycle_summary.csv**:
- Always output in full (already small)

### 5. Helper functions

```julia
function _should_output_step(csv_opt::CsvExportOptions, cycle::Int, ti::Int, n_steps::Int)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return (ti == 1 || ti == n_steps)
    csv_opt.mode == :custom && return (ti == 1 || ti == n_steps || ti % csv_opt.save_every == 1)
    return true
end

"""Pre-compute set of (cycle, phase) keys for the last snapshot in each group."""
function _compute_last_snap_keys(czm_snapshots)
    last_keys = Set{Tuple{Int,String}}()
    prev_key = nothing
    for snap in czm_snapshots
        key = (snap.cycle, snap.phase)
        if prev_key !== nothing && key != prev_key
            push!(last_keys, prev_key)
        end
        prev_key = key
    end
    if prev_key !== nothing
        push!(last_keys, prev_key)
    end
    return last_keys
end

function _should_output_snapshot(csv_opt::CsvExportOptions, cycle::Int, phase::String,
                                  is_last_in_phase::Bool)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return is_last_in_phase
    csv_opt.mode == :custom && return true
    return true
end
```

Phase boundary detection for snapshot filtering:
```julia
# In _write_cohesive_damage and _write_node_displacement:
last_snap_keys = _compute_last_snap_keys(result.czm_snapshots)
for snap in result.czm_snapshots
    key = (snap.cycle, snap.phase)
    is_last = key in last_snap_keys
    if !_should_output_snapshot(csv_opt, snap.cycle, snap.phase, is_last)
        continue
    end
    # ... write data
end
```

### 6. skip_files support

In `export_cycling_csv`, check `csv_opt.skip_files` before writing each file. Skipped files are added to `files_skipped` list with a descriptive reason.

### 7. Input validation

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

### 8. Estimated compression

With 10 cycles × ~500 steps/cycle × ~80 elements:

- **`:phase_ends`**: ~80 rows vs ~5000 rows per file → **~60x compression**
- **`:phase_ends` + `full_output_cycles=[1]`**: ~580 rows → **~8x compression**
- Total reduction: ~1.86 GB → **~20-30 MB** (with `:phase_ends` mode)

## Files to modify

1. `src/CsvExport.jl` — Add struct, helpers, modify write functions
2. `example/coupled_czm_thermal_example.jl` — Demonstrate new options

## Testing

- Run with default `:full` mode, verify output matches original exactly
- Run with `:phase_ends` mode, verify CSV files are correct and reduced in size
- Run with `full_output_cycles=[1,5,10]`, verify those cycles have full output
- Run with `skip_files=["node_displacement.csv"]`, verify the file is not created
- Run with `:custom` mode and `save_every=50`, verify sampling interval
