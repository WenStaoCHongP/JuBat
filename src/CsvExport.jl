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

function _validate_csv_opt(csv_opt::CsvExportOptions)
    if csv_opt.mode ∉ (:full, :phase_ends, :custom)
        error("CsvExportOptions: invalid mode $(csv_opt.mode), must be :full, :phase_ends, or :custom")
    end
    if csv_opt.mode == :custom && csv_opt.save_every < 1
        error("CsvExportOptions: save_every must be ≥ 1 for :custom mode")
    end
end

"""Whether time step `ti` (of `n_steps` total) in `cycle` should be written."""
function _should_output_step(csv_opt::CsvExportOptions, cycle::Int, ti::Int, n_steps::Int)
    cycle in csv_opt.full_output_cycles && return true
    csv_opt.mode == :full && return true
    csv_opt.mode == :phase_ends && return (ti == 1 || ti == n_steps)
    csv_opt.mode == :custom && return (ti == 1 || ti == n_steps || ti % csv_opt.save_every == 1)
    return true
end

"""Pre-compute Set of snapshot indices that correspond to the last snapshot in each phase group.
Returns empty Set for empty snapshot arrays."""
function _compute_last_snap_indices(czm_snapshots)
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
function _should_output_snapshot(csv_opt::CsvExportOptions, cycle::Int, snap_idx::Int,
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

"""
    export_cycling_csv(result, case, czm_mesh;
                       output_dir="output/csv", overwrite=false,
                       csv_opt=CsvExportOptions())

Export solve_cycling results to CSV files for post-processing.
"""
function export_cycling_csv(result, case, czm_mesh;
                            output_dir::String="output/csv", overwrite::Bool=false,
                            csv_opt::CsvExportOptions=CsvExportOptions())
    _validate_csv_opt(csv_opt)
    mkpath(output_dir)

    files_written = String[]
    files_skipped = String[]

    # 1. cycle_summary.csv (always full output)
    try
        _write_cycle_summary(result, output_dir, overwrite)
        push!(files_written, "cycle_summary.csv")
    catch e
        @warn "Failed to write cycle_summary.csv" exception=e
        push!(files_skipped, "cycle_summary.csv")
    end

    # 2. element_currents.csv
    has_data = !isempty(result.cycle_results) &&
               any(cr -> cr.discharge.solve_result !== nothing, result.cycle_results)
    if has_data && !("element_currents.csv" in csv_opt.skip_files)
        try
            _write_element_currents(result, case, output_dir, overwrite, csv_opt)
            push!(files_written, "element_currents.csv")
        catch e
            @warn "Failed to write element_currents.csv" exception=e
            push!(files_skipped, "element_currents.csv")
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
        try
            _write_node_temperature(result, case, output_dir, overwrite, csv_opt)
            push!(files_written, "node_temperature.csv")
        catch e
            @warn "Failed to write node_temperature.csv" exception=e
            push!(files_skipped, "node_temperature.csv")
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
            try
                _write_cohesive_damage(result, case, czm_mesh, output_dir, overwrite, csv_opt)
                push!(files_written, "cohesive_damage.csv")
            catch e
                @warn "Failed to write cohesive_damage.csv" exception=e
                push!(files_skipped, "cohesive_damage.csv")
            end
        else
            push!(files_skipped, "cohesive_damage.csv (skipped by csv_opt)")
        end
        if !("node_displacement.csv" in csv_opt.skip_files)
            try
                _write_node_displacement(result, case, czm_mesh, output_dir, overwrite, csv_opt)
                push!(files_written, "node_displacement.csv")
            catch e
                @warn "Failed to write node_displacement.csv" exception=e
                push!(files_skipped, "node_displacement.csv")
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
            try
                _write_driving_force(result, case, czm_mesh, output_dir, overwrite, csv_opt)
                push!(files_written, "cohesive_driving_force.csv")
            catch e
                @warn "Failed to write cohesive_driving_force.csv" exception=e
                push!(files_skipped, "cohesive_driving_force.csv")
            end
        else
            push!(files_skipped, "cohesive_driving_force.csv (skipped by csv_opt)")
        end
    else
        push!(files_skipped, "cohesive_driving_force.csv (no CZM snapshots or geometry)")
    end

    # 7. czm_solver_diagnostics.csv (always full output)
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

function _write_cycle_summary(result, output_dir::String, overwrite::Bool)
    filepath = joinpath(output_dir, "cycle_summary.csv")
    if isfile(filepath) && !overwrite
        println("  Skipping $filepath (already exists)")
        return
    end

    open(filepath, "w") do f
        println(f, "cycle,phase,capacity_ah,soh,D_max,D_mean,n_fractured,T_max_K,T_mean_end_K,V_start,V_end")

        for (i, cr) in enumerate(result.cycle_results)
            cyc = cr.cycle_idx
            soh_val = i <= length(result.soh) ? result.soh[i] : NaN
            for (phase_name, pr) in [("discharge", cr.discharge),
                                      ("rest1", cr.rest1),
                                      ("charge", cr.charge),
                                      ("rest2", cr.rest2)]
                D_max = isnan(pr.D_max) ? 0.0 : pr.D_max
                D_mean = isnan(pr.D_mean) ? 0.0 : pr.D_mean
                T_mean_end = isnan(pr.T_mean_end) ? 0.0 : pr.T_mean_end
                V_start = isnan(pr.V_start) ? 0.0 : pr.V_start
                V_end = isnan(pr.V_end) ? 0.0 : pr.V_end
                println(f, "$cyc,$phase_name,$(pr.capacity),$soh_val,$D_max,$D_mean,$(cr.n_fractured),$(cr.T_max),$T_mean_end,$V_start,$V_end")
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 2. element_currents.csv
# ========================================================================

function _write_element_currents(result, case, output_dir::String, overwrite::Bool,
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
    elem_areas = _compute_element_areas(mesh_th)

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

                time_s = get(sr, "time [s]", Float64[])
                isempty(time_s) && continue

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
                    if !_should_output_step(csv_opt, cyc, ti, n_steps)
                        continue
                    end
                    t = time_s[ti]
                    for e in 1:ne
                        I_e = _safe_get(I_e_all, e, ti)
                        sn = _safe_get(soc_n_all, e, ti)
                        sp = _safe_get(soc_p_all, e, ti)
                        en = _safe_get(eta_n_all, e, ti)
                        ep = _safe_get(eta_p_all, e, ti)
                        qr = _safe_get(q_rxn_ne_all, e, ti)
                        qs = _safe_get(q_sp_all, e, ti)
                        qp = _safe_get(q_rxn_pe_all, e, ti)
                        Te = _safe_get(T_elem_all, e, ti)
                        area_phys = e <= length(elem_areas) ? elem_areas[e] * scale.L^2 : 0.0

                        println(f, "$t,$cyc,$phase_name,$e,$I_e,$area_phys,$Te,$sn,$sp,$en,$ep,$qr,$qs,$qp")
                    end
                end
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 3. node_temperature.csv
# ========================================================================

function _write_node_temperature(result, case, output_dir::String, overwrite::Bool,
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
                isempty(time_s) && continue
                T_nodes_all = get(sr, "thermal2D temperature at nodes [K]", zeros(0,0))
                n_steps = length(time_s)

                for ti in 1:n_steps
                    if !_should_output_step(csv_opt, cyc, ti, n_steps)
                        continue
                    end
                    t = time_s[ti]
                    for n in 1:nnode
                        T_K = _safe_get(T_nodes_all, n, ti)
                        println(f, "$t,$cyc,$phase_name,$n,$(node_x[n]),$(node_y[n]),$T_K")
                    end
                end
            end
        end
    end
    println("  Written: $filepath")
end

# ========================================================================
# 4. cohesive_damage.csv
# ========================================================================

function _write_cohesive_damage(result, case, czm_mesh,
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

    last_snap_indices = _compute_last_snap_indices(result.czm_snapshots)

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,length,D,delta_n,delta_t,T_n,T_t,fractured,theta_deg")

        for (si, snap) in enumerate(result.czm_snapshots)
            if !_should_output_snapshot(csv_opt, snap.cycle, si, last_snap_indices)
                continue
            end
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            for i in 1:min(n_coh, length(snap.damage))
                D = snap.damage[i]
                dn = i <= length(snap.separation_n) ? snap.separation_n[i] * scale.r0 : 0.0
                dt_val = i <= length(snap.separation_t) ? snap.separation_t[i] * scale.r0 : 0.0
                tn = i <= length(snap.traction_n) ? snap.traction_n[i] * scale.σ_czm : 0.0
                tt = i <= length(snap.traction_t) ? snap.traction_t[i] * scale.σ_czm : 0.0
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

function _write_node_displacement(result, case, czm_mesh,
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

    last_snap_indices = _compute_last_snap_indices(result.czm_snapshots)

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,node_id,x,y,ux,uy")

        for (si, snap) in enumerate(result.czm_snapshots)
            if !_should_output_snapshot(csv_opt, snap.cycle, si, last_snap_indices)
                continue
            end
            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase
            nd = div(length(snap.displacement), 2)
            for n in 1:min(nnode, nd)
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

function _write_driving_force(result, case, czm_mesh,
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

    last_snap_indices = _compute_last_snap_indices(result.czm_snapshots)

    open(filepath, "w") do f
        println(f, "time_s,cycle,phase,coh_id,dT_neighbor,dsoc_n_neighbor,dsoc_p_neighbor,eps_0_thermal,eps_0_diffusion,eps_0_total")

        for (si, snap) in enumerate(result.czm_snapshots)
            if !_should_output_snapshot(csv_opt, snap.cycle, si, last_snap_indices)
                continue
            end

            t = snap.time_s
            cyc = snap.cycle
            phase = snap.phase

            sr = _find_solve_result(result, cyc, phase)
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

function _write_czm_diagnostics(result, output_dir::String, overwrite::Bool)
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

"""Safe 2D array access with NaN fallback"""
function _safe_get(arr::AbstractArray{<:Real,2}, row::Int, col::Int)
    if row <= size(arr, 1) && col <= size(arr, 2)
        return arr[row, col]
    end
    return NaN
end
_safe_get(arr, row, col) = NaN

"""Find solve_result for a specific cycle+phase"""
function _find_solve_result(result, cycle::Int, phase::String)
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

"""Compute normalized element areas from thermal mesh (Q4 elements)"""
function _compute_element_areas(mesh::Mesh)
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
