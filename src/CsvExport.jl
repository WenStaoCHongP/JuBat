# src/CsvExport.jl
# CSV export for cycling simulation post-processing
# Note: CZMSnapshot struct is defined in CouplingState.jl (included before Solve.jl)

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
"""
struct CsvExportOptions
    mode::Symbol
    save_every::Int
    full_output_cycles::Vector{Int}
    skip_files::Vector{String}
end

CsvExportOptions() = CsvExportOptions(:full, 1, [1], String[])
CsvExportOptions(mode::Symbol) = CsvExportOptions(mode, 1, [1], String[])

function validate_csv_options(csv_opt::CsvExportOptions)
    if csv_opt.mode ∉ (:full, :phase_ends, :custom)
        error("CsvExportOptions: invalid mode $(csv_opt.mode), must be :full, :phase_ends, or :custom")
    end
    if csv_opt.mode == :custom && csv_opt.save_every < 1
        error("CsvExportOptions: save_every must be ≥ 1 for :custom mode")
    end
end

"""Whether time step `ti` (of `n_steps` total) in `cycle` should be written."""
function should_export_step(csv_opt::CsvExportOptions, cycle::Int, ti::Int, n_steps::Int)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return (ti == 1 || ti == n_steps)
    csv_opt.mode == :custom && return (ti == 1 || ti == n_steps || ti % csv_opt.save_every == 1)
    return true
end

"""Pre-compute Set of snapshot indices that correspond to the last snapshot in each phase group.
Returns empty Set for empty snapshot arrays."""
function compute_last_snapshot_indices(czm_snapshots)
    last_indices = Set{Int}()
    isempty(czm_snapshots) && return last_indices
    prev_key = (czm_snapshots[1].cycle, czm_snapshots[1].phase)
    prev_idx = 1
    for i in 2:length(czm_snapshots)
        key = (czm_snapshots[i].cycle, czm_snapshots[i].phase)
        if key != prev_key
            push!(last_indices, prev_idx)
        end
        prev_key = key
        prev_idx = i
    end
    push!(last_indices, prev_idx)
    return last_indices
end

"""Whether a CZM snapshot at index `snap_idx` should be written."""
function should_export_snapshot(csv_opt::CsvExportOptions, cycle::Int, snap_idx::Int,
                                  last_indices::Set{Int})
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return snap_idx in last_indices
    csv_opt.mode == :custom && return true
    return true
end

# ========================================================================
# Main export function
# ========================================================================

function write_csv_guarded!(write_fn::Function, filename::String,
                             files_written::Vector{String}, files_skipped::Vector{String})
    try
        write_fn()
        push!(files_written, filename)
    catch e
        @warn "Failed to write $filename" exception=e
        push!(files_skipped, filename)
    end
    return nothing
end

"""
    export_cycling_csv(result, case, czm_mesh;
                       output_dir="output/csv", overwrite=false,
                       csv_opt=CsvExportOptions())

Export solve_cycling results to CSV files for post-processing.
"""
function export_cycling_csv(result, case, czm_mesh;
                            output_dir::String="output/csv", overwrite::Bool=false,
                            csv_opt::CsvExportOptions=CsvExportOptions())
    validate_csv_options(csv_opt)
    mkpath(output_dir)

    files_written = String[]
    files_skipped = String[]

    # 1. cycle_summary.csv (always full output)
    write_csv_guarded!("cycle_summary.csv", files_written, files_skipped) do
        write_cycle_summary_csv(result, output_dir, overwrite)
    end

    # 2. element_currents.csv
    has_data = !isempty(result.cycle_results) &&
               any(cr -> cr.discharge.solve_result !== nothing, result.cycle_results)
    if has_data && !("element_currents.csv" in csv_opt.skip_files)
        write_csv_guarded!("element_currents.csv", files_written, files_skipped) do
            write_element_currents_csv(result, case, output_dir, overwrite, csv_opt)
        end
    else
        if "element_currents.csv" in csv_opt.skip_files
            push!(files_skipped, "element_currents.csv (skipped by csv_opt)")
        else
            push!(files_skipped, "element_currents.csv (no solve_result data)")
        end
    end

    # 3. node_temperature.csv
    if case.opt.thermal_enabled && !("node_temperature.csv" in csv_opt.skip_files)
        write_csv_guarded!("node_temperature.csv", files_written, files_skipped) do
            write_node_temperature_csv(result, case, output_dir, overwrite, csv_opt)
        end
    else
        if "node_temperature.csv" in csv_opt.skip_files
            push!(files_skipped, "node_temperature.csv (skipped by csv_opt)")
        else
            push!(files_skipped, "node_temperature.csv (thermal disabled)")
        end
    end

    # 4. cohesive_damage.csv + 5. node_displacement.csv
    if !isempty(result.czm_snapshots)
        if !("cohesive_damage.csv" in csv_opt.skip_files)
            write_csv_guarded!("cohesive_damage.csv", files_written, files_skipped) do
                write_cohesive_damage_csv(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            end
        else
            push!(files_skipped, "cohesive_damage.csv (skipped by csv_opt)")
        end
        if !("node_displacement.csv" in csv_opt.skip_files)
            write_csv_guarded!("node_displacement.csv", files_written, files_skipped) do
                write_node_displacement_csv(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            end
        else
            push!(files_skipped, "node_displacement.csv (skipped by csv_opt)")
        end
    else
        push!(files_skipped, "cohesive_damage.csv (no CZM snapshots)")
        push!(files_skipped, "node_displacement.csv (no CZM snapshots)")
    end

    # 6. cohesive_driving_force.csv
    if !isempty(result.czm_snapshots) && case.geometry !== nothing
        if !("cohesive_driving_force.csv" in csv_opt.skip_files)
            write_csv_guarded!("cohesive_driving_force.csv", files_written, files_skipped) do
                write_cohesive_driving_force_csv(result, case, czm_mesh, output_dir, overwrite, csv_opt)
            end
        else
            push!(files_skipped, "cohesive_driving_force.csv (skipped by csv_opt)")
        end
    else
        push!(files_skipped, "cohesive_driving_force.csv (no CZM snapshots or geometry)")
    end

    # 7. czm_solver_diagnostics.csv (always full output)
    if !isempty(result.czm_snapshots)
        write_csv_guarded!("czm_solver_diagnostics.csv", files_written, files_skipped) do
            write_czm_diagnostics_csv(result, output_dir, overwrite)
        end
    else
        push!(files_skipped, "czm_solver_diagnostics.csv (no CZM snapshots)")
    end

    println("CSV export complete:")
    println("  Written: $(join(files_written, ", "))")
    if !isempty(files_skipped)
        println("  Skipped: $(join(files_skipped, ", "))")
    end
    println("  Directory: $output_dir")
    return files_written
end

# ========================================================================
# 1. cycle_summary.csv
# ========================================================================

function write_cycle_summary_csv(result, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "cycle_summary.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    n_cycles = length(result.cycle_results)
    length(result.soh) == n_cycles || throw(DimensionMismatch(
        "cycle summary SOH length $(length(result.soh)) does not match cycle count $n_cycles"))

    open(filepath, "w") do f
        println(f, "cycle,phase,capacity_ah,soh,D_max,D_mean,n_fractured,T_max_K,T_mean_end_K,V_start,V_end")

        for (i, cr) in enumerate(result.cycle_results)
            cyc = cr.cycle_idx
            soh_val = result.soh[i]
            for (phase_name, pr) in existing_cycle_phases(cr)
                println(f, "$cyc,$phase_name,$(pr.capacity),$soh_val,$(pr.D_max),$(pr.D_mean),$(cr.n_fractured),$(cr.T_max),$(pr.T_mean_end),$(pr.V_start),$(pr.V_end)")
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 2. element_currents.csv
# ========================================================================

function write_element_currents_csv(result, case, output_dir::String, overwrite::Bool,
                                  csv_opt::CsvExportOptions)
    filepath = joinpath(output_dir, "element_currents.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    # Compute element areas from mesh geometry (normalized)
    elem_areas = compute_csv_element_areas(mesh_th)
    length(elem_areas) == ne || throw(DimensionMismatch(
        "computed element area count $(length(elem_areas)) does not match thermal element count $ne"))

    phase_exports = Any[]
    for cr in result.cycle_results
        cyc = cr.cycle_idx
        for (phase_name, pr) in existing_cycle_phases(cr)
            context = "cycle $cyc phase $phase_name"
            sr = require_csv_solve_result(pr, context)
            time_s = require_csv_vector(sr, "time [s]", context)
            n_steps = length(time_s)
            fields = (
                I_e = require_csv_matrix(sr, "thermal2D element current", ne, n_steps, context),
                soc_n = require_csv_matrix(sr, "thermal2D element soc_n", ne, n_steps, context),
                soc_p = require_csv_matrix(sr, "thermal2D element soc_p", ne, n_steps, context),
                eta_n = require_csv_matrix(sr, "thermal2D eta_n_e", ne, n_steps, context),
                eta_p = require_csv_matrix(sr, "thermal2D eta_p_e", ne, n_steps, context),
                q_rxn_ne = require_csv_matrix(sr, "thermal2D Q_rxn_NE [W/m3]", ne, n_steps, context),
                q_sp = require_csv_matrix(sr, "thermal2D Q_SP [W/m3]", ne, n_steps, context),
                q_rxn_pe = require_csv_matrix(sr, "thermal2D Q_rxn_PE [W/m3]", ne, n_steps, context),
                temperature = require_csv_matrix(sr, "thermal2D temperature [K]", ne, n_steps, context),
            )
            push!(phase_exports, (cyc=cyc, phase_name=phase_name, time_s=time_s, fields=fields))
        end
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,elem_id,I_e,area,T_e,soc_n,soc_p,eta_n,eta_p,q_rxn_ne,q_sp,q_rxn_pe")

        for data in phase_exports
            n_steps = length(data.time_s)
            for ti in 1:n_steps
                should_export_step(csv_opt, data.cyc, ti, n_steps) || continue
                t = data.time_s[ti]
                for e in 1:ne
                    area_phys = elem_areas[e] * scale.L^2
                    println(f, "$t,$(data.cyc),$(data.phase_name),$e,$(data.fields.I_e[e, ti]),$area_phys,$(data.fields.temperature[e, ti]),$(data.fields.soc_n[e, ti]),$(data.fields.soc_p[e, ti]),$(data.fields.eta_n[e, ti]),$(data.fields.eta_p[e, ti]),$(data.fields.q_rxn_ne[e, ti]),$(data.fields.q_sp[e, ti]),$(data.fields.q_rxn_pe[e, ti])")
                end
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 3. node_temperature.csv
# ========================================================================

function write_node_temperature_csv(result, case, output_dir::String, overwrite::Bool,
                                  csv_opt::CsvExportOptions)
    filepath = joinpath(output_dir, "node_temperature.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    mesh_th = case.mesh["thermal2D"]
    nnode = mesh_th.nlen
    node_x = mesh_th.node[:, 1] * scale.L
    node_y = mesh_th.node[:, 2] * scale.L

    phase_exports = Any[]
    for cr in result.cycle_results
        cyc = cr.cycle_idx
        for (phase_name, pr) in existing_cycle_phases(cr)
            context = "cycle $cyc phase $phase_name"
            sr = require_csv_solve_result(pr, context)
            time_s = require_csv_vector(sr, "time [s]", context)
            temperatures = require_csv_matrix(
                sr, "thermal2D temperature at nodes [K]", nnode, length(time_s), context)
            push!(phase_exports, (cyc=cyc, phase_name=phase_name,
                                  time_s=time_s, temperatures=temperatures))
        end
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,node_id,x,y,T_K")

        for data in phase_exports
            n_steps = length(data.time_s)
            for ti in 1:n_steps
                should_export_step(csv_opt, data.cyc, ti, n_steps) || continue
                t = data.time_s[ti]
                for n in 1:nnode
                    println(f, "$t,$(data.cyc),$(data.phase_name),$n,$(node_x[n]),$(node_y[n]),$(data.temperatures[n, ti])")
                end
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 4. cohesive_damage.csv
# ========================================================================

function write_cohesive_damage_csv(result, case, czm_mesh,
                                 output_dir::String, overwrite::Bool,
                                 csv_opt::CsvExportOptions)
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
        n1, n2 = elem.nodes_bottom
        mx = 0.5 * (czm_mesh.node[n1, 1] + czm_mesh.node[n2, 1])
        my = 0.5 * (czm_mesh.node[n1, 2] + czm_mesh.node[n2, 2])
        push!(theta_degs, atan(my, mx) * 180.0 / pi)
    end
    length(czm_mesh.cohesive_elements) == n_coh || throw(DimensionMismatch(
        "cohesive element count $(length(czm_mesh.cohesive_elements)) does not match n_cohesive $n_coh"))

    last_snap_indices = compute_last_snapshot_indices(result.czm_snapshots)
    selected_snapshots = Any[]
    for (si, snap) in enumerate(result.czm_snapshots)
        should_export_snapshot(csv_opt, snap.cycle, si, last_snap_indices) || continue
        context = "CZM snapshot $si (cycle $(snap.cycle), phase $(snap.phase))"
        require_csv_length(snap.damage, n_coh, "damage", context)
        require_csv_length(snap.separation_n, n_coh, "separation_n", context)
        require_csv_length(snap.separation_t, n_coh, "separation_t", context)
        require_csv_length(snap.traction_n, n_coh, "traction_n", context)
        require_csv_length(snap.traction_t, n_coh, "traction_t", context)
        push!(selected_snapshots, snap)
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,length,D,delta_n,delta_t,T_n,T_t,fractured,theta_deg")

        for snap in selected_snapshots
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            for i in 1:n_coh
                D = snap.damage[i]
                # 分离位移以 scale.δ_czm 归一（重设计 v2；修正原误用 scale.r0 颗粒半径尺度）
                dn = snap.separation_n[i] * scale.δ_czm
                dt_val = snap.separation_t[i] * scale.δ_czm
                tn = snap.traction_n[i] * scale.σ_czm
                tt = snap.traction_t[i] * scale.σ_czm
                frac = D >= 0.95
                println(f, "$t,$cyc,$phase,$i,$(lengths_phys[i]),$D,$dn,$dt_val,$tn,$tt,$frac,$(theta_degs[i])")
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 5. node_displacement.csv
# ========================================================================

function write_node_displacement_csv(result, case, czm_mesh,
                                   output_dir::String, overwrite::Bool,
                                   csv_opt::CsvExportOptions)
    filepath = joinpath(output_dir, "node_displacement.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    nnode = czm_mesh.nnode
    node_x = czm_mesh.node[:, 1] * scale.L
    node_y = czm_mesh.node[:, 2] * scale.L

    last_snap_indices = compute_last_snapshot_indices(result.czm_snapshots)
    selected_snapshots = Any[]
    for (si, snap) in enumerate(result.czm_snapshots)
        should_export_snapshot(csv_opt, snap.cycle, si, last_snap_indices) || continue
        require_csv_length(snap.displacement, 2 * nnode, "displacement",
            "CZM snapshot $si (cycle $(snap.cycle), phase $(snap.phase))")
        push!(selected_snapshots, snap)
    end

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,node_id,x,y,ux,uy")

        for snap in selected_snapshots
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

# ========================================================================
# 6. cohesive_driving_force.csv
# ========================================================================

function write_cohesive_driving_force_csv(result, case, czm_mesh,
                               output_dir::String, overwrite::Bool,
                               csv_opt::CsvExportOptions)
    filepath = joinpath(output_dir, "cohesive_driving_force.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    scale = case.param.scale
    n_coh = czm_mesh.n_cohesive
    T_ref = scale.T_ref

    # CZM effective params (thermal/diffusion strain coupling)
    # These are NOT fields on Cohesive struct; default to 0.0 until
    # compute_czm_effective_params(case) is integrated.
    alpha_eff = 0.0
    beta_n = 0.0
    beta_p = 0.0

    if alpha_eff == 0.0 && beta_n == 0.0 && beta_p == 0.0
        println("  Skipping $filepath (no thermal/diffusion strain parameters)")
        return
    end

    geo = case.geometry
    if geo === nothing
        println("  Skipping $filepath (no geometry data)")
        return
    end

    # Build cohesive-to-thermal-element mapping
    # interface_pairs is Vector{Tuple{Int,Int}}
    interface_pairs = geo.interface_pairs
    coh_thermal_map = Dict{Int, Tuple{Int,Int}}()
    for (idx, pair) in enumerate(interface_pairs)
        if idx <= n_coh
            coh_thermal_map[idx] = pair
        end
    end

    last_snap_indices = compute_last_snapshot_indices(result.czm_snapshots)

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,dT_neighbor,dsoc_n_neighbor,dsoc_p_neighbor,eps_0_thermal,eps_0_diffusion,eps_0_total")

        for (si, snap) in enumerate(result.czm_snapshots)
            if !should_export_snapshot(csv_opt, snap.cycle, si, last_snap_indices)
                continue
            end

            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase

            sr = find_cycle_solve_result(result, cyc, phase)
            if sr === nothing
                continue
            end

            time_arr = get(sr, "time [s]", Float64[])
            isempty(time_arr) && continue
            T_elem = get(sr, "thermal2D temperature [K]", zeros(0,0))
            soc_n_elem = get(sr, "thermal2D element soc_n", zeros(0,0))
            soc_p_elem = get(sr, "thermal2D element soc_p", zeros(0,0))

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

# ========================================================================
# 7. czm_solver_diagnostics.csv
# ========================================================================

function write_czm_diagnostics_csv(result, output_dir::String, overwrite::Bool)
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

# ========================================================================
# Helpers
# ========================================================================

function require_csv_solve_result(phase_result, context::String)
    phase_result.solve_result === nothing && error("$context has no solve_result")
    phase_result.solve_result isa AbstractDict || throw(ArgumentError(
        "$context solve_result must be a dictionary, got $(typeof(phase_result.solve_result))"))
    return phase_result.solve_result
end

function require_csv_vector(solve_result::AbstractDict, key::String, context::String)
    haskey(solve_result, key) || throw(KeyError("$context: $key"))
    value = solve_result[key]
    value isa AbstractVector || throw(DimensionMismatch(
        "$context field \"$key\" must be a vector, got size $(size(value))"))
    isempty(value) && error("$context field \"$key\" must not be empty")
    return value
end

function require_csv_matrix(solve_result::AbstractDict, key::String,
                            nrows::Int, ncols::Int, context::String)
    haskey(solve_result, key) || throw(KeyError("$context: $key"))
    value = solve_result[key]
    value isa AbstractMatrix || throw(DimensionMismatch(
        "$context field \"$key\" must be a matrix, got size $(size(value))"))
    size(value) == (nrows, ncols) || throw(DimensionMismatch(
        "$context field \"$key\" has size $(size(value)); expected ($nrows, $ncols)"))
    return value
end

function require_csv_length(value, expected::Int, name::String, context::String)
    value isa AbstractVector || throw(DimensionMismatch(
        "$context field \"$name\" must be a vector, got $(typeof(value))"))
    length(value) == expected || throw(DimensionMismatch(
        "$context field \"$name\" has length $(length(value)); expected $expected"))
    return value
end

function existing_cycle_phases(cycle_result)
    phases = Tuple{String,PhaseResult}[("discharge", cycle_result.discharge)]
    cycle_result.rest1 === nothing || push!(phases, ("rest1", cycle_result.rest1))
    push!(phases, ("charge", cycle_result.charge))
    cycle_result.rest2 === nothing || push!(phases, ("rest2", cycle_result.rest2))
    return phases
end

"""Find solve_result for a specific cycle+phase"""
function find_cycle_solve_result(result, cycle::Int, phase::String)
    for cr in result.cycle_results
        if cr.cycle_idx == cycle
            if phase == "discharge"
                return cr.discharge.solve_result
            elseif phase == "rest1"
                return cr.rest1 === nothing ? nothing : cr.rest1.solve_result
            elseif phase == "charge"
                return cr.charge.solve_result
            elseif phase == "rest2"
                return cr.rest2 === nothing ? nothing : cr.rest2.solve_result
            end
        end
    end
    return nothing
end

"""Compute normalized element areas from thermal mesh (Q4 elements)"""
function compute_csv_element_areas(mesh::Mesh)
    ne = size(mesh.element, 1)
    areas = zeros(ne)
    for e in 1:ne
        n1, n2, n3, n4 = mesh.element[e, 1], mesh.element[e, 2], mesh.element[e, 3], mesh.element[e, 4]
        x1, y1 = mesh.node[n1, 1], mesh.node[n1, 2]
        x2, y2 = mesh.node[n2, 1], mesh.node[n2, 2]
        x3, y3 = mesh.node[n3, 1], mesh.node[n3, 2]
        x4, y4 = mesh.node[n4, 1], mesh.node[n4, 2]
        # Shoelace formula for quadrilateral area
        areas[e] = 0.5 * abs((x1*y2 - x2*y1) + (x2*y3 - x3*y2) +
                             (x3*y4 - x4*y3) + (x4*y1 - x1*y4))
    end
    return areas
end

# ========================================================================
# Raw single-cycle data export
# ========================================================================

"""
    export_cycle_data_to_csv(export_data, output_dir; prefix="cycle")

将循环数据导出为CSV文件。

# 输出文件
- `{prefix}_timesteps.csv`: 时间步汇总数据
- `{prefix}_T_nodes.csv`: 节点温度场（每行一个时间步，每列一个节点）
- `{prefix}_soc_n.csv`: 负极SOC场（每行一个时间步，每列一个单元）
- `{prefix}_soc_p.csv`: 正极SOC场（每行一个时间步，每列一个单元）
- `{prefix}_mesh_nodes.csv`: 网格节点坐标
- `{prefix}_mesh_elements.csv`: 网格单元连接
"""
function export_cycle_data_to_csv(export_data::CycleExportData, output_dir::String;
                                   prefix::String="cycle")
    isdir(output_dir) || mkpath(output_dir)

    n_steps = length(export_data.timesteps)
    nT = export_data.nT
    ne = export_data.ne

    # 1. 时间步汇总数据
    timesteps_file = joinpath(output_dir, "$(prefix)_timesteps.csv")
    open(timesteps_file, "w") do io
        println(io, "step,time_s,phase,V,I_A,T_max_K,T_mean_K,soc_mean")
        for (i, ts) in enumerate(export_data.timesteps)
            phase_str = ts.phase == PHASE_DISCHARGE ? "discharge" :
                        ts.phase == PHASE_CHARGE ? "charge" : "rest"
            @printf(io, "%d,%.6f,%s,%.6f,%.6f,%.4f,%.4f,%.6f\n",
                    i, ts.time, phase_str, ts.V, ts.I, ts.T_max, ts.T_mean, ts.soc_mean)
        end
    end
    println("  ✓ 保存: $timesteps_file")

    # 2. 节点温度场
    T_nodes_file = joinpath(output_dir, "$(prefix)_T_nodes.csv")
    open(T_nodes_file, "w") do io
        header = join(["node_$(i)" for i in 1:nT], ",")
        println(io, "step,time_s,$header")
        for (i, ts) in enumerate(export_data.timesteps)
            T_str = join([@sprintf("%.4f", T) for T in ts.T_nodes], ",")
            @printf(io, "%d,%.6f,%s\n", i, ts.time, T_str)
        end
    end
    println("  ✓ 保存: $T_nodes_file")

    # 3. 负极SOC场
    soc_n_file = joinpath(output_dir, "$(prefix)_soc_n.csv")
    open(soc_n_file, "w") do io
        header = join(["elem_$(i)" for i in 1:ne], ",")
        println(io, "step,time_s,$header")
        for (i, ts) in enumerate(export_data.timesteps)
            soc_str = join([@sprintf("%.6f", s) for s in ts.soc_n], ",")
            @printf(io, "%d,%.6f,%s\n", i, ts.time, soc_str)
        end
    end
    println("  ✓ 保存: $soc_n_file")

    # 4. 正极SOC场
    soc_p_file = joinpath(output_dir, "$(prefix)_soc_p.csv")
    open(soc_p_file, "w") do io
        header = join(["elem_$(i)" for i in 1:ne], ",")
        println(io, "step,time_s,$header")
        for (i, ts) in enumerate(export_data.timesteps)
            soc_str = join([@sprintf("%.6f", s) for s in ts.soc_p], ",")
            @printf(io, "%d,%.6f,%s\n", i, ts.time, soc_str)
        end
    end
    println("  ✓ 保存: $soc_p_file")

    # 5. 网格节点坐标
    nodes_file = joinpath(output_dir, "$(prefix)_mesh_nodes.csv")
    open(nodes_file, "w") do io
        println(io, "node_id,x,y")
        for i in 1:nT
            @printf(io, "%d,%.8f,%.8f\n", i, export_data.node_coords[i, 1], export_data.node_coords[i, 2])
        end
    end
    println("  ✓ 保存: $nodes_file")

    # 6. 网格单元连接
    elements_file = joinpath(output_dir, "$(prefix)_mesh_elements.csv")
    open(elements_file, "w") do io
        println(io, "elem_id,n1,n2,n3,n4")
        for e in 1:ne
            nodes = export_data.element_connectivity[e, :]
            println(io, "$e,$(nodes[1]),$(nodes[2]),$(nodes[3]),$(nodes[4])")
        end
    end
    println("  ✓ 保存: $elements_file")

    return (timesteps_file, T_nodes_file, soc_n_file, soc_p_file, nodes_file, elements_file)
end
