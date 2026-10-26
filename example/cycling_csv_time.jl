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
        JuBat.check_czm_snapshot(snapshot)
        t = snapshot.time_s
        isfinite(t) && 0.0 <= t <= phase.duration + 16eps(max(1.0, phase.duration)) ||
            error("snapshot time outside its recorded phase: $(snapshot.cycle)/$(snapshot.phase)/$t")
        (; snapshot, phase=snapshot.phase, time_s=phase.t_start + t)
    end
    sort!(rows; by=row -> row.time_s)
    return rows
end

# The two cycling examples share the physical snapshot tables and their units.
# Validate all rows before opening files, so a failed result never becomes a physical CSV.
function write_czm_snapshot_csv(result, case, outdir; cohesive_long=false)
    rows = cycling_csv_snapshots(result)
    isempty(rows) && return rows
    mesh = case.mesh["czm"]
    n_coh = length(mesh.cohesive_elements)
    for row in rows
        s = row.snapshot
        length(s.damage) == n_coh && length(s.displacement) == 2size(mesh.node,1) ||
            throw(DimensionMismatch("snapshot fields differ from the exported mechanical mesh"))
        s.gp_points == JuBat.NCweight(mesh.czm_submesh.mesh_bonded.gs.order)[2] ||
            throw(DimensionMismatch("snapshot GP point order differs from the exported mechanical mesh"))
        if s.plasticity
            size(s.bulk_kappa,1) == size(mesh.bulk_element,1) ||
                throw(DimensionMismatch("snapshot bulk kappa differs from the exported bulk mesh"))
        end
    end
    scale = case.param.scale
    mkpath(outdir)
    open(joinpath(outdir,"element_map.csv"),"w") do io
        println(io,"elem_id,interface,x_m,y_m")
        for (j,elem) in enumerate(mesh.cohesive_elements)
            ns = vcat(elem.nodes_bottom,elem.nodes_top)
            x = sum(mesh.node[ns,1]) / length(ns) * scale.L
            y = sum(mesh.node[ns,2]) / length(ns) * scale.L
            println(io,"$j,$(string(elem.interface_type)),$x,$y")
        end
    end
    header = "t_s," * join(("e$j" for j in 1:n_coh),",")
    for cycle in sort(unique(row.snapshot.cycle for row in rows))
        cycle_rows = filter(row -> row.snapshot.cycle == cycle,rows)
        for (name,values) in (
            ("damage_history",s -> s.damage),
            ("sep_n_history",s -> s.separation_n .* scale.δ_czm),
            ("sep_t_history",s -> s.separation_t .* scale.δ_czm),
            ("traction_n_history",s -> s.traction_n .* scale.σ_czm),
            ("traction_t_history",s -> s.traction_t .* scale.σ_czm),
            ("damage_visc_history",s -> s.damage_visc),
            ("fractured_history",s -> s.fractured))
            open(joinpath(outdir,"$(name)_cyc$cycle.csv"),"w") do io
                println(io,header)
                for row in cycle_rows
                    println(io,string(row.time_s) * "," * join(values(row.snapshot),","))
                end
            end
        end
    end
    if cohesive_long
        open(joinpath(outdir,"cohesive_damage.csv"),"w") do io
            println(io,"cycle,phase,time_s,coh_id,D,sep_n_m,sep_t_m,traction_n_Pa,traction_t_Pa,D_visc,fractured")
            for row in rows
                s = row.snapshot
                for j in eachindex(s.damage)
                    println(io,join((s.cycle,row.phase,row.time_s,j,s.damage[j],
                        s.separation_n[j]*scale.δ_czm,s.separation_t[j]*scale.δ_czm,
                        s.traction_n[j]*scale.σ_czm,s.traction_t[j]*scale.σ_czm,
                        s.damage_visc[j],s.fractured[j]),","))
                end
            end
        end
    end
    open(joinpath(outdir,"cohesive_gp_history.csv"),"w") do io
        println(io,"cycle,phase,time_s,mechanical_dt_s,coh_id,gp_id,xi,D,D_visc,fractured,sep_n_m,sep_t_m,traction_n_Pa,traction_t_Pa,delta_max_n_m,delta_max_t_m,delta_max_eff_m,czm_model,plasticity_enabled")
        for row in rows
            s = row.snapshot
            for e in axes(s.gp_damage_states,1), g in axes(s.gp_damage_states,2)
                state = s.gp_damage_states[e,g]
                println(io,join((s.cycle,row.phase,row.time_s,s.mechanical_dt_s,e,g,s.gp_points[g],
                    state.D,state.D_visc,state.fractured,s.gp_separation_n[e,g]*scale.δ_czm,
                    s.gp_separation_t[e,g]*scale.δ_czm,s.gp_traction_n[e,g]*scale.σ_czm,
                    s.gp_traction_t[e,g]*scale.σ_czm,state.δ_max_n*scale.δ_czm,
                    state.δ_max_t*scale.δ_czm,state.δ_max_eff*scale.δ_czm,s.czm_model,s.plasticity),","))
            end
        end
    end
    # No rows for J2-off snapshots: absence is not a measured zero kappa.
    open(joinpath(outdir,"bulk_plastic_history.csv"),"w") do io
        println(io,"cycle,phase,time_s,mechanical_dt_s,bulk_elem_id,gp_id,kappa")
        for row in rows
            s = row.snapshot
            s.plasticity || continue
            for e in axes(s.bulk_kappa,1), g in axes(s.bulk_kappa,2)
                println(io,join((s.cycle,row.phase,row.time_s,s.mechanical_dt_s,e,g,s.bulk_kappa[e,g]),","))
            end
        end
    end
    return rows
end

function cycling_csv_times(phase, ncols)
    local_time = vec(phase.solve_result["time [s]"])
    length(local_time) == ncols || throw(DimensionMismatch("CSV field columns differ from recorded times"))
    all(isfinite, local_time) && issorted(local_time) || error("invalid recorded CSV times")
    return phase.t_start .+ local_time
end
