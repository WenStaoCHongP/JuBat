using CSV, DataFrames
include("../src/JuBat.jl")
using Plots
# Use GR backend to ensure file-based rendering works in headless environments
gr()
using Dates

# Jellyroll coupled SPMe <-> 2D thermal example
# Environment overrides (optional):
#   JUBAT_MODE = strong|weak
#   JUBAT_UNITS = SI|nd
#   JUBAT_PER_ELEM_SPMe = 1|0
#   JUBAT_PARALLEL_V = 1|0
#   JUBAT_COLLECTOR_SEEDED = 1|0
#   JUBAT_NTHETA = integer (segments per turn for collector-seeded)
#   JUBAT_TIME = total seconds (e.g., 60)
#   JUBAT_DT0 = min dt seconds; JUBAT_DT1 = max dt seconds
#   JUBAT_QUICK = 1 to use a faster small run
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
# enable thermal solve / mesh creation
opt.thermal_enabled = true
# explicitly use SPMe model for electrochemistry
opt.model = "SPMe"
# Use collector-seeded mesh and distributed 2D thermal model
opt.collector_seeded = true
opt.per_element_spme = true
opt.parallel_solve_V = true
opt.coupling_mode = "strong"
opt.thermalmodel = "distributed2D"
opt.units_thermal = "SI"
opt.debug_coupling = true
opt.debug_sample_elems = true

# Apply ENV overrides
mode_env = get(ENV, "JUBAT_MODE", "")
if !isempty(mode_env)
    opt.coupling_mode = lowercase(mode_env) in ("strong","weak") ? lowercase(mode_env) : opt.coupling_mode
end
units_env = get(ENV, "JUBAT_UNITS", "")
if !isempty(units_env)
    opt.units_thermal = (uppercase(units_env) == "SI") ? "SI" : "nd"
end
opt.per_element_spme = parse(Bool, get(ENV, "JUBAT_PER_ELEM_SPMe", string(opt.per_element_spme)))
opt.parallel_solve_V = parse(Bool, get(ENV, "JUBAT_PARALLEL_V", string(opt.parallel_solve_V)))
opt.collector_seeded = parse(Bool, get(ENV, "JUBAT_COLLECTOR_SEEDED", string(opt.collector_seeded)))

# time and current settings (with optional quick mode)
quick = parse(Bool, get(ENV, "JUBAT_QUICK", "false"))
# default run
opt.dt = [1.0, 10.0]
opt.time = [0.0, 60.0]
if quick
    opt.dt = [0.5, 5.0]
    opt.time = [0.0, 20.0]
end
if haskey(ENV, "JUBAT_DT0")
    try opt.dt[1] = parse(Float64, ENV["JUBAT_DT0"]) catch end
end
if haskey(ENV, "JUBAT_DT1")
    try opt.dt[2] = parse(Float64, ENV["JUBAT_DT1"]) catch end
end
if haskey(ENV, "JUBAT_TIME")
    try opt.time[2] = parse(Float64, ENV["JUBAT_TIME"]) catch end
end

# modest discharge current for example
opt.Current = t -> 1.0 * param_dim.cell.I1C

case = JuBat.SetCase(param_dim, opt)
println("Config: mode=$(opt.coupling_mode), units=$(opt.units_thermal), per-element=$(opt.per_element_spme), V-solve=$(opt.parallel_solve_V), collector-seeded=$(opt.collector_seeded)")
println("Starting jellyroll coupled example")
# create a thermal mesh (collector-seeded) and attach to case so thermal assembly runs
try
    nθ = try parse(Int, get(ENV, "JUBAT_NTHETA", quick ? "120" : "180")) catch; 180 end
    mesh_th = JuBat.jellyroll_Q4_mesh(param_dim; nx=nθ, gsorder=2, crop_mode=:collector_seeded)
    case.mesh["thermal2D"] = mesh_th
    println("Created thermal2D mesh with ", size(mesh_th.element,1), " elements; nθ=", nθ)
catch e
    println("Warning creating thermal mesh: ", e)
end

res = nothing
try
    global res = JuBat.Solve(case)
catch e
    # Strong coupling / V-solve path may fail in some setups; fallback to conservative options
    println("Warning: initial Solve failed with error: ", e)
    println("Falling back: disabling per-element SPMe and parallel V-solve, switching to weak coupling and re-running")
    opt.per_element_spme = false
    opt.parallel_solve_V = false
    opt.coupling_mode = "weak"
    # recreate case and thermal mesh
    case = JuBat.SetCase(param_dim, opt)
    try
    mesh_th = JuBat.jellyroll_Q4_mesh(param_dim; nx=120, gsorder=2, crop_mode=:collector_seeded)
        case.mesh["thermal2D"] = mesh_th
    println("(Fallback) Created thermal2D mesh with ", size(mesh_th.element,1), " elements")
    catch e2
        println("Warning creating thermal mesh in fallback: ", e2)
    end
    try
        global res = JuBat.Solve(case)
    catch e3
        println("Fallback Solve also failed: ", e3)
        rethrow(e3)
    end
end
println("Finished")
println("Debug: res is ", res === nothing ? "nothing" : "not nothing")
if res !== nothing
    println("Debug: res keys: ", collect(keys(res)))
end

# --- collect initial / mid / final temperature fields ---
# nodes (x,y) available from final result if thermal mesh exists
if res === nothing
    println("No valid solution result; skipping temperature field analysis")
    println("Example script finished.")
    exit(0)
end

nodes = haskey(res, "thermal2D nodes xy [m]") ? res["thermal2D nodes xy [m]"] : (haskey(case.mesh, "thermal2D") ? case.mesh["thermal2D"].node : Float64[])
tf = opt.time[2]
tmid = 0.5 * tf

println("Diagnostic: nodes size = ", size(nodes))

# initial temperature field (K)
T0 = fill(param_dim.cell.T0, size(nodes,1))

# mid-run: run a shorter solve to get mid temperature snapshot
opt_mid = deepcopy(opt)
opt_mid.time = [0.0, tmid]
case_mid = JuBat.SetCase(param_dim, opt_mid)
println("Running mid-point solve to t=", tmid, " s")
res_mid = JuBat.Solve(case_mid)
T_mid = haskey(res_mid, "thermal2D T_nodes [K]") ? res_mid["thermal2D T_nodes [K]"] : T0

# final temperature (from earlier run)
T_final = haskey(res, "thermal2D T_nodes [K]") ? res["thermal2D T_nodes [K]"] : T0

# plot three snapshots side-by-side
if size(nodes,1) > 0
    x = nodes[:,1]; y = nodes[:,2]
    p = plot(layout=(1,3), size=(1200,400))
    # compute unified color limits so all three panels use same scale
    all_T = vcat(T0, T_mid, T_final)
    tmin = minimum(all_T)
    tmax = maximum(all_T)
    clim = (tmin, tmax)

    scatter!(p[1], x, y, marker_z=T0, ms=6, mc=:inferno, markerstrokewidth=0, colorbar=true, title = "Initial T (t=0s)", clim=clim)
    scatter!(p[2], x, y, marker_z=T_mid, ms=6, mc=:inferno, markerstrokewidth=0, colorbar=true, title = "Mid T (t=$(Int(round(tmid)))s)", clim=clim)
    scatter!(p[3], x, y, marker_z=T_final, ms=6, mc=:inferno, markerstrokewidth=0, colorbar=true, title = "Final T (t=$(Int(round(tf)))s)", clim=clim)

    outdir = joinpath(@__DIR__, "..", "output")
    # ensure output directory exists
    try
        if !isdir(outdir)
            mkpath(outdir)
            println("Created output directory: ", outdir)
        else
            println("Output directory exists: ", outdir)
        end
    catch e
        println("Warning: could not create output directory: ", e)
    end
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
    outpath = joinpath(outdir, "jellyroll_T_snapshots_$(timestamp).png")
    savefig(p, outpath)
    println("Saved $(outpath)")
else
    # if no nodes, write a diagnostic file so user can see output dir and variables
    outdir = joinpath(@__DIR__, "..", "output")
    try
        if !isdir(outdir)
            mkpath(outdir)
        end
        diagpath = joinpath(outdir, "jellyroll_diag_$(Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")).txt")
        open(diagpath, "w") do io
            println(io, "No thermal nodes found; nodes size: ", size(nodes))
            println(io, "res keys: ", collect(keys(res)))
            println(io, "case.mesh keys: ", collect(keys(case.mesh)))
            if haskey(res, "variables")
                println(io, "variables keys: ", collect(keys(res["variables"])))
            end
        end
        println("Wrote diagnostic file: ", diagpath)
    catch e
        println("Warning writing diagnostic file: ", e)
    end
end

# Try to print diagnostics if available
variables = haskey(res, "variables") ? res["variables"] : nothing
if variables !== nothing
    # scales for unit conversion
    I1C = param_dim.cell.I1C
    phi = param_dim.scale.phi
    q_ref = param_dim.scale.q_th
    # current diagnostics (note: currents are nondimensional vs I1C)
    if haskey(variables, "thermal2D element current")
        Ivec_nd = variables["thermal2D element current"]
        Ivec_A = Ivec_nd .* I1C
        println("element current count: ", length(Ivec_nd))
        println("I_e (nd) min/mean/max: ", minimum(Ivec_nd), ", ", mean(Ivec_nd), ", ", maximum(Ivec_nd))
        println("I_e [A]  min/mean/max: ", minimum(Ivec_A),  ", ", mean(Ivec_A),  ", ", maximum(Ivec_A))
        # conservation check vs total current at final time
        I_tot_nd = 0.0
        if haskey(variables, "cell current")
            I_tot_nd = variables["cell current"]
        end
        println("sum(Ie) nd = ", sum(Ivec_nd), ", I_total nd = ", I_tot_nd, ", abs diff = ", abs(sum(Ivec_nd) - I_tot_nd))
    end
    # heat source diagnostics (convert to SI for printing)
    if haskey(variables, "heat_source_fields")
        q = variables["heat_source_fields"]
        is_SI = false
        if haskey(variables, "heat_source_units_code")
            code = variables["heat_source_units_code"]
            is_SI = (isa(code, Float64) && code > 0.5) || (isa(code, AbstractVector) && length(code)>0 && code[1] > 0.5)
        end
        q_SI = is_SI ? q : (q .* q_ref)
        println("heat_source_fields shape: ", size(q))
        println("q [W/m^3] min/mean/max: ", minimum(q_SI), ", ", mean(q_SI), ", ", maximum(q_SI))
    end
    if haskey(variables, "thermal2D common voltage")
        V_nd = variables["thermal2D common voltage"]
        println("common voltage: ", V_nd, " (nd), ", V_nd*phi, " [V]")
    end
    if haskey(variables, "thermal2D energy residual")
        println("energy residual (nd): ", variables["thermal2D energy residual"])
    end
end

# Save a small CSV for inspection
try
    if variables !== nothing && haskey(variables, "thermal2D element current")
        df = DataFrame(Ie_nd = variables["thermal2D element current"], Ie_A = variables["thermal2D element current"] .* param_dim.cell.I1C)
        csvpath = joinpath(outdir, "jellyroll_Ie_$(timestamp).csv")
        CSV.write(csvpath, df)
        println("Saved $(csvpath)")
    end
    if variables !== nothing && haskey(variables, "heat_source_fields")
        q = variables["heat_source_fields"]
        is_SI = false
        if haskey(variables, "heat_source_units_code")
            code = variables["heat_source_units_code"]
            is_SI = (isa(code, Float64) && code > 0.5) || (isa(code, AbstractVector) && length(code)>0 && code[1] > 0.5)
        end
        if is_SI
            df2 = DataFrame(q_Wm3 = q)
            csvpath2 = joinpath(outdir, "jellyroll_q_SI_$(timestamp).csv")
            CSV.write(csvpath2, df2)
            println("Saved $(csvpath2)")
        else
            df2 = DataFrame(q_nd = q, q_Wm3 = q .* param_dim.scale.q_th)
            csvpath2 = joinpath(outdir, "jellyroll_q_nd_$(timestamp).csv")
            CSV.write(csvpath2, df2)
            println("Saved $(csvpath2)")
        end
    end
catch e
    println("Warning writing CSV: ", e)
end

println("Example script finished.")
