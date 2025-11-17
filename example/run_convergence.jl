using LinearAlgebra, SparseArrays, Statistics, Plots
include("../src/JuBat.jl")

function run_case(nx; ny=nx, total_time_s=20.0, dt_s=0.5)
    # load params
    param_dim = include("../src/parameters/ThermalUniform.jl")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.units_thermal = "nd"
    opt.time = [0.0, max(total_time_s, 1.0)]
    opt.dt = [dt_s, dt_s]
    opt.solveType = "backward"
    case = JuBat.SetCase(param_dim, opt)

    L = param_dim.cell.length
    mesh_sq = JuBat.SetMesh([0.0, L, 0.0, L], [nx, ny], "Q4", 2)
    case.mesh["thermal2D"] = mesh_sq

    # analytic initial sine mode
    xy = mesh_sq.node
    k = param_dim.PE.lambda
    rho = param_dim.cell.rho
    cp = param_dim.cell.heat_Q
    alpha_SI = k / (rho * cp)

    # nondim scaling
    t_th = max(case.param_dim.scale.t_th, 1e-16)
    dt_th = dt_s / t_th
    nsteps = Int(clamp(floor(total_time_s / dt_s), 1, 10^6))

    # initial condition in nd
    T_amp_SI = 1.0
    T_amp_nd = T_amp_SI / case.param_dim.scale.T_ref
    nN = size(xy,1)
    Tnum = zeros(Float64, nN)
    for i in 1:nN
        x = xy[i,1]; y = xy[i,2]
        Tnum[i] = T_amp_nd * sin(pi * x / L) * sin(pi * y / L)
    end

    vars = JuBat.StandardVariables(case, 1)
    vars["T_nodes"] = copy(Tnum)
    vars["heat_source_fields"] = zeros(Float64, size(mesh_sq.element,1))
    vars["heat_source_units_code"] = 0.0

    for step in 1:nsteps
    MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
    JuBat.ThermalDistributed2D_BC(KT, FT, case, 0.0)
        T_prev = Tnum
        A = (1.0/dt_th) .* MT + KT
        rhsT = (1.0/dt_th) .* (MT * T_prev) + FT
        T_new = A \ rhsT
        Tnum = T_new
        vars["T_nodes"] = Tnum
    end

    # compute final L2 error between numeric and analytic (SI)
    Tnum_SI = Tnum .* case.param_dim.scale.T_ref
    t_final = nsteps * dt_s
    decay = exp(-alpha_SI * 2*pi^2 / (L^2) * t_final)
    Tana_SI = similar(Tnum_SI)
    for i in 1:nN
        x = xy[i,1]; y = xy[i,2]
        Tana_SI[i] = T_amp_SI * decay * sin(pi * x / L) * sin(pi * y / L)
    end
    err = abs.(Tnum_SI .- Tana_SI)
    L2 = sqrt(mean(err.^2))
    return L2, L / nx
end

function main()
    ns = [32, 64, 128]
    errs = Float64[]; dxs = Float64[]; labels = String[]
    for nx in ns
        println("Running nx=$nx ...")
        L2, dx = run_case(nx)
        push!(errs, L2); push!(dxs, dx); push!(labels, string(nx))
        println("nx=$nx dx=$(round(dx, sigdigits=4)) L2_err=$L2")
    end

    # log-log fit via linear least squares on log data
    logdx = log.(dxs); logerr = log.(errs)
    A = hcat(logdx, ones(length(logdx)))
    coeffs = A \ logerr
    slope = coeffs[1]; intercept = coeffs[2]
    println("Fitted slope (convergence order) = ", slope)

    plt = scatter(dxs, errs, xscale=:log10, yscale=:log10, xlabel="Δx [m]", ylabel="L2 error [K]",
        title = "Convergence: L2 error vs Δx", label="data", ms=6)
    xs = range(minimum(dxs), stop=maximum(dxs), length=50)
    fitline = exp.(intercept) .* xs .^ (slope)
    plot!(plt, xs, fitline, lw=2, label = "fit slope=$(round(slope, digits=3))")
    savefig(plt, "convergence_L2_vs_dx.pdf")
    println("Saved convergence plot: convergence_L2_vs_dx.pdf")
end

main()
