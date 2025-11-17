using LinearAlgebra, SparseArrays, Statistics, Plots
include("../src/JuBat.jl")

function main(; nx::Int=64, ny::Int=64, total_time_s::Float64=20.0, dt_s::Float64=0.5, q_SI::Float64=1e5)
    # Use our new uniform parameter set by loading directly
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

    # mesh square [0,L]x[0,L]
    L = param_dim.cell.length
    mesh_sq = JuBat.SetMesh([0.0, L, 0.0, L], [nx, ny], "Q4", 2)
    case.mesh["thermal2D"] = mesh_sq

    vars = JuBat.StandardVariables(case, 1)
    # initialize to ambient
    vars["T_nodes"] = fill(param_dim.cell.T0 / case.param_dim.scale.T_ref, mesh_sq.nlen)
    # heat source: convert SI volumetric heat [W/m^3] to non-dimensional using case.param_dim.scale if needed
    # for simplicity we place q in SI and mark units code = 1 so ThermalDistributed2D assembly will treat it consistently
    q_elem = fill(q_SI, size(mesh_sq.element,1))
    vars["heat_source_fields"] = q_elem
    vars["heat_source_units_code"] = 1.0

    # time scaling
    t_th = max(case.param_dim.scale.t_th, 1e-16)
    dt_th = dt_s / t_th
    nsteps = Int(clamp(floor(total_time_s / dt_s), 1, 10^6))

    T_hist_mean = Float64[]; t_hist = Float64[]

    for step in 1:nsteps
        MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
    JuBat.ThermalDistributed2D_BC(KT, FT, case, 0.0)
        T_prev = vars["T_nodes"]
        A = (1.0/dt_th) .* MT + KT
        rhsT = (1.0/dt_th) .* (MT * T_prev) + FT
        T_new = A \ rhsT
        vars["T_nodes"] = T_new
        push!(T_hist_mean, mean(T_new) * case.param_dim.scale.T_ref)
        push!(t_hist, step * dt_s)
    end

    # analytic validation: set initial condition to sine mode and compare decay
    println("Running analytic (1,1) sine-mode comparison")
    xy = mesh_sq.node
    # alpha from uniform k, rho, cp
    k = param_dim.PE.lambda
    rho = param_dim.cell.rho
    cp = param_dim.cell.heat_Q
    alpha_SI = k / (rho * cp)

    # set initial T = sin(pi x/L) sin(pi y/L) amplitude 1 K
    T_amp_SI = 1.0
    T_amp_nd = T_amp_SI / case.param_dim.scale.T_ref
    nN = size(xy,1)
    T0 = zeros(Float64, nN)
    for i in 1:nN
        x = xy[i,1]; y = xy[i,2]
        T0[i] = T_amp_nd * sin(pi * x / L) * sin(pi * y / L)
    end

    vars["heat_source_fields"] .= 0.0
    vars["heat_source_units_code"] = 0.0
    vars["T_nodes"] = T0

    analytic_err_L2 = Float64[]; analytic_err_max = Float64[]; analytic_center_num = Float64[]; analytic_center_ana = Float64[]
    distances = sum((xy .- [L/2 L/2]).^2, dims=2); center_idx = argmin(vec(distances))

    Tnum = copy(T0)
    for step in 1:nsteps
        MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
    JuBat.ThermalDistributed2D_BC(KT, FT, case, 0.0)
        T_prev = Tnum
        A = (1.0/dt_th) .* MT + KT
        rhsT = (1.0/dt_th) .* (MT * T_prev) + FT
        T_new = A \ rhsT
        Tnum = T_new
        vars["T_nodes"] = Tnum

        Tnum_SI = Tnum .* case.param_dim.scale.T_ref
        tnow = step * dt_s
        decay = exp(-alpha_SI * 2*pi^2 / (L^2) * tnow)
        Tana_SI = similar(Tnum_SI)
        for i in 1:nN
            x = xy[i,1]; y = xy[i,2]
            Tana_SI[i] = T_amp_SI * decay * sin(pi * x / L) * sin(pi * y / L)
        end
        err = abs.(Tnum_SI .- Tana_SI)
        push!(analytic_err_L2, sqrt(mean(err.^2)))
        push!(analytic_err_max, maximum(err))
        push!(analytic_center_num, Tnum_SI[center_idx])
        push!(analytic_center_ana, Tana_SI[center_idx])
    end

    println("alpha_SI = ", alpha_SI)
    println("Final L2 error [K]: ", analytic_err_L2[end])
    println("Final max abs error [K]: ", analytic_err_max[end])

    # save plots
    Tfield = vars["T_nodes"] .* case.param_dim.scale.T_ref
    Tmin, Tmax = minimum(Tfield), maximum(Tfield)
    cmap = cgrad(:viridis, 256)
    plt = plot(xlabel="x [m]", ylabel="y [m]", aspect_ratio=1, title="Uniform square final T field", legend=false)
    ne = size(mesh_sq.element, 1)
    for e in 1:ne
        n1,n2,n3,n4 = mesh_sq.element[e,1], mesh_sq.element[e,2], mesh_sq.element[e,3], mesh_sq.element[e,4]
        xs = [xy[n1,1], xy[n2,1], xy[n3,1], xy[n4,1], xy[n1,1]]
        ys = [xy[n1,2], xy[n2,2], xy[n3,2], xy[n4,2], xy[n1,2]]
        Te = (Tfield[n1] + Tfield[n2] + Tfield[n3] + Tfield[n4]) / 4
        α = (Te - Tmin) / max(Tmax - Tmin, eps())
        idx = clamp(1 + floor(Int, α * 255), 1, 256)
        col = cmap.colors[idx]
        plot!(plt, xs, ys, seriestype=:shape, c=col, linecolor=:black, lw=0.2, label=false)
    end
    savefig(plt, "uniform_T_field.png")

    p1 = plot(collect(1:nsteps) .* dt_s, analytic_err_L2, xlabel="time [s]", ylabel="L2 error [K]",
              title="L2 error vs time (analytic mode)", lw=2)
    savefig(p1, "uniform_analytic_L2_error.pdf")

    p2 = plot(collect(1:nsteps) .* dt_s, analytic_center_num, lw=2, label="numeric")
    plot!(p2, collect(1:nsteps) .* dt_s, analytic_center_ana, lw=2, ls=:dash, label="analytic")
    xlabel!(p2, "time [s]"); ylabel!(p2, "T at center [K]"); title!(p2, "Center temperature: numeric vs analytic")
    savefig(p2, "uniform_analytic_center_compare.pdf")

    println("Outputs: uniform_T_field.png, uniform_analytic_L2_error.pdf, uniform_analytic_center_compare.pdf")
end

main()
