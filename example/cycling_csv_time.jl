# CSV time coordinates for the two cycling examples. Solver snapshots remain local.
function cycling_csv_phases(result)
    Dict((cr.cycle_idx, name) => ph for cr in result.cycle_results
        for (name, ph) in (("discharge", cr.discharge), ("rest1", cr.rest1),
                           ("charge", cr.charge), ("rest2", cr.rest2)) if ph !== nothing)
end

function cycling_csv_snapshots(result)
    phases = cycling_csv_phases(result)
    rows = map(result.czm_snapshots) do snapshot
        phase = phases[(snapshot.cycle, snapshot.phase)]
        t = snapshot.time_s
        isfinite(t) && 0.0 <= t <= phase.duration + 16eps(max(1.0, phase.duration)) ||
            error("snapshot time outside its recorded phase: $(snapshot.cycle)/$(snapshot.phase)/$t")
        (; snapshot, phase=snapshot.phase, time_s=phase.t_start + t)
    end
    sort!(rows; by=row -> row.time_s)
    return rows
end

function cycling_csv_times(phase, ncols)
    local_time = vec(phase.solve_result["time [s]"])
    length(local_time) == ncols || throw(DimensionMismatch("CSV field columns differ from recorded times"))
    all(isfinite, local_time) && issorted(local_time) || error("invalid recorded CSV times")
    return phase.t_start .+ local_time
end
