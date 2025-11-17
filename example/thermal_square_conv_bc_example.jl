using LinearAlgebra, SparseArrays, Statistics, Plots
include("../src/JuBat.jl")

"""
Square 2D thermal diffusion with convection on all four sides

Problem setup:
- Domain: square [0, L] × [0, L], L = 100 × (t_pos + t_sep + t_neg) from cell params (dimensional)
- Mesh: structured Q4, nx = ny = 64 (configurable)
- PDE: (ρ c) ∂T/∂t = ∇·(k ∇T) + q, with q=0 by default
- Boundary: convection at all sides: -k ∂T/∂n = h (T - T_amb)
- Units: use Scheme-B dimensionless thermal internally; final plots in SI [K]

This example validates the ThermalDistributed2D assembly + ThermalDistributed2D_BC on a simple square.
"""
function main(; nx::Int=64, ny::Int=64, total_time_s::Float64=60.0, dt_s::Float64=0.5, q0_nd::Float64=1.0,
                 validate_analytic::Bool=false, T_amp_SI::Float64=1.0)
    # 1) Options and params (no electrochemistry, only thermal)
    opt = JuBat.Option()
    opt.model = "SPMe"            # electrochemical model present but unused (q=0); keep defaults
    opt.Nn = 8; opt.Ns = 4; opt.Np = 8
    opt.Nrn = 6; opt.Nrp = 6
    opt.dimension = 1
    opt.gsorder = 2
    # time in electrochem scale; we only care thermal → pick a simple window
    opt.time = [0.0, max(total_time_s, 1.0)]
    opt.dt = [dt_s, dt_s]
    opt.solveType = "backward"
    # Thermal 2D toggles
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.units_thermal = "nd"      # keep dimensionless internally
    opt.coupling_mode = "loose"   # no strong coupling needed here
    opt.per_element_spme = false
    opt.parallel_solve_V = false
    opt.debug_coupling = false

    # Use a standard cell to get dimensional parameters
    param_dim = JuBat.ChooseCell("Jellyroll")
    # Compute square length L = 100 × (t_pos + t_sep + t_neg)
    t_pos = hasproperty(param_dim.PE, :thickness) ? param_dim.PE.thickness : 0.0
    t_sep = hasproperty(param_dim.SP, :thickness) ? param_dim.SP.thickness : 0.0
    t_neg = hasproperty(param_dim.NE, :thickness) ? param_dim.NE.thickness : 0.0
    L = 100.0 * (t_pos + t_sep + t_neg)
    L = max(L, 1e-3)  # ensure non-zero size

    # Convection coefficient and ambient; keep the defaults already placed into scaling via param_dim.scale.h_th
    # If needed, you can tweak param_dim.cell.h and param_dim.cell.T_amb before SetCase

    # 2) Case
    case = JuBat.SetCase(param_dim, opt)

    # 3) Build a square Q4 mesh in SI meters
    mesh_sq = JuBat.SetMesh([0.0, L, 0.0, L], [nx, ny], "Q4", opt.gsorder)
    case.mesh["thermal2D"] = mesh_sq

    # 4) Initial variables and constant heat source (nd)
    #    T_nodes in dimensionless units expected by ThermalDistributed2D path
    vars = JuBat.StandardVariables(case, 1)
    # default: uniform initial temperature (dimensionless).
    vars["T_nodes"] = fill(case.param.cell.T0, mesh_sq.nlen)
    # constant internal source (non-dimensional)
    vars["heat_source_fields"] = fill(q0_nd, size(mesh_sq.element,1))
    vars["heat_source_units_code"] = 0.0  # 0=nd, 1=SI(W/m^3)

    # 5) Assemble and step in time: BE scheme on thermal alone
    # Use thermal scaling only: t_nd = t_s / t_th
    t_th = max(case.param_dim.scale.t_th, 1e-16)
    dt_th = dt_s / t_th
    nsteps = Int(clamp(floor(total_time_s / dt_s), 1, 10^7))

    # storage for plots
    T_hist_mean = Float64[]
    t_hist = Float64[]

    for step in 1:nsteps
    # keep constant heat source each step (can be modified to time-dependent if needed)
    vars["heat_source_fields"] .= q0_nd
        # If analytic validation requested, compute analytic solution at current time
        # (we compute analytic in SI and compare to numeric converted to SI below)
    MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
    JuBat.ThermalDistributed2D_BC(KT, FT, case, step * dt_s)
        T_prev = vars["T_nodes"]
        A = (1.0/dt_th) .* MT + KT
        rhsT = (1.0/dt_th) .* (MT * T_prev) + FT
        # mild anchoring to ambient temp for robustness
        alpha = 1e-12
        if alpha > 0
            nT = size(A,1)
            @inbounds for i in 1:nT
                A[i,i] += alpha
            end
            T_amb_nd = case.param_dim.cell.T_amb / case.param_dim.scale.T_ref
            rhsT .+= alpha .* T_amb_nd
        end
        T_new = A \ rhsT
        vars["T_nodes"] = T_new
        push!(T_hist_mean, mean(T_new) * case.param_dim.scale.T_ref)
        push!(t_hist, step * dt_s)
    end

    # If analytic validation enabled, re-run a short analytic comparison run (we computed nothing yet)
    if validate_analytic
        # Recompute analytic history by stepping again but now capturing analytic error each step.
        println("Running analytic validation: initial sine mode, comparing numeric field to analytic decay.")
        # Prepare analytic initial condition and overwrite vars and forcing
        xy = mesh_sq.node
        # Ensure ambient = 0 and very large h to approximate Dirichlet BC
        case.param_dim.cell.T_amb = 0.0
        case.param_dim.cell.h = 1e9
        # pick k, rho, cp from params for alpha
        # choose a reference conductivity from available layers
        k_ref = hasproperty(case.param_dim.PE, :lambda) ? case.param_dim.PE.lambda : 0.0
        if !(k_ref > 0)
            k_ref = hasproperty(case.param_dim.NE, :lambda) ? case.param_dim.NE.lambda : (hasproperty(case.param_dim.SP, :lambda) ? case.param_dim.SP.lambda : 1.0)
        end
        rho_eff = case.param_dim.cell.rho
        cp_eff = case.param_dim.cell.heat_Q
        alpha_SI = k_ref / (rho_eff * cp_eff)

        # set initial numeric field to the dimensionless version of T_amp*sin(pi x/L) sin(pi y/L)
        T_amp_nd = T_amp_SI / case.param_dim.scale.T_ref
        nN = size(xy,1)
        for i in 1:nN
            x = xy[i,1]; y = xy[i,2]
            vars["T_nodes"][i] = T_amp_nd * sin(pi * x / L) * sin(pi * y / L)
        end

        # zero internal source
        vars["heat_source_fields"] .= 0.0

        # histories for analytic comparison
        analytic_err_L2 = Float64[]
        analytic_err_max = Float64[]
        analytic_center_num = Float64[]
        analytic_center_ana = Float64[]

        # find index closest to center (L/2,L/2)
        distances = sum((xy .- [L/2 L/2]).^2, dims=2)
        center_idx = argmin(vec(distances))

        # time stepping again from initial condition
        Tnum = copy(vars["T_nodes"]) # dimensionless
        for step in 1:nsteps
            MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
            JuBat.ThermalDistributed2D_BC(KT, FT, case, step * dt_s)
            T_prev = Tnum
            A = (1.0/dt_th) .* MT + KT
            rhsT = (1.0/dt_th) .* (MT * T_prev) + FT
            alpha_reg = 1e-12
            if alpha_reg > 0
                nT = size(A,1)
                @inbounds for i in 1:nT
                    A[i,i] += alpha_reg
                end
                T_amb_nd = case.param_dim.cell.T_amb / case.param_dim.scale.T_ref
                rhsT .+= alpha_reg .* T_amb_nd
            end
            T_new = A \ rhsT
            Tnum = T_new
            vars["T_nodes"] = Tnum

            # numeric in SI
            Tnum_SI = Tnum .* case.param_dim.scale.T_ref

            # analytic in SI for mode (1,1): T = T_amp * exp(-alpha * 2*pi^2 / L^2 * t) * sin(pi x/L) sin(pi y/L)
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

        println("Analytic validation finished.")
        println("alpha_SI = ", alpha_SI)
        println("Final L2 error [K]: ", analytic_err_L2[end])
        println("Final max abs error [K]: ", analytic_err_max[end])

        # quick plots for comparison
        p1 = plot(collect(1:nsteps) .* dt_s, analytic_err_L2, xlabel="time [s]", ylabel="L2 error [K]",
                  title="L2 error vs time (analytic mode)", lw=2)
        savefig(p1, "analytic_validation_L2_error.pdf")

        p2 = plot(collect(1:nsteps) .* dt_s, analytic_center_num, lw=2, label="numeric")
        plot!(p2, collect(1:nsteps) .* dt_s, analytic_center_ana, lw=2, ls=:dash, label="analytic")
        xlabel!(p2, "time [s]"); ylabel!(p2, "T at center [K]"); title!(p2, "Center temperature: numeric vs analytic")
        savefig(p2, "analytic_validation_center_compare.pdf")
    end

    # 6) Post plots: field and mean temperature history
    xy = mesh_sq.node
    Tn = vars["T_nodes"] .* case.param_dim.scale.T_ref
    Tmin, Tmax = minimum(Tn), maximum(Tn)
    cmap = cgrad(:turbo, 256, categorical=false)
    plt = plot(xlabel="x [m]", ylabel="y [m]", aspect_ratio=1,
               title="Square 2D thermal: final T field", legend=false)
    ne = size(mesh_sq.element, 1)
    for e in 1:ne
        n1,n2,n3,n4 = mesh_sq.element[e,1], mesh_sq.element[e,2], mesh_sq.element[e,3], mesh_sq.element[e,4]
        xs = [xy[n1,1], xy[n2,1], xy[n3,1], xy[n4,1], xy[n1,1]]
        ys = [xy[n1,2], xy[n2,2], xy[n3,2], xy[n4,2], xy[n1,2]]
        Te = (Tn[n1] + Tn[n2] + Tn[n3] + Tn[n4]) / 4
        α = (Te - Tmin) / max(Tmax - Tmin, eps())
        idx = clamp(1 + floor(Int, α * 255), 1, 256)
        col = cmap.colors[idx]
        plot!(plt, xs, ys, seriestype=:shape, c=col, linecolor=:gray70, lw=0.2, label=false)
    end
    scatter!(plt, [0.0, 0.0], [0.0, 0.0], marker_z=[Tmin, Tmax], ms=0,
             c=:turbo, colorbar_title="T [K]", clims=(Tmin, Tmax))
    savefig(plt, "thermal_square_conv_field.png")

    pmean = plot(t_hist, T_hist_mean, xlabel="time [s]", ylabel="mean T [K]",
                 lw=2, grid=true, title="Mean temperature vs time (square conv)")
    savefig(pmean, "thermal_square_conv_Tmean.pdf")

    # 6.1) Temperature distribution along centerline y=L/2 at final time
    cy = L/2
    dy = L / max(ny, 1)
    tol = dy/4
    # collect nodes near the horizontal centerline
    idxs = findall(abs.(xy[:,2] .- cy) .<= tol)
    if !isempty(idxs)
        xs = xy[idxs, 1]
        Ts = Tn[idxs]
        # sort by x
        ord = sortperm(xs)
        xs_sorted = xs[ord]
        Ts_sorted = Ts[ord]
        pprof = plot(xs_sorted, Ts_sorted, lw=2, xlabel="x [m]", ylabel="T [K]", grid=true,
                     title="Centerline temperature profile (y = L/2)")
        savefig(pprof, "thermal_square_conv_T_profile_centerline.png")
        savefig(pprof, "thermal_square_conv_T_profile_centerline.pdf")
    end

    println("Square 2D convection example finished. Outputs: thermal_square_conv_field.png, thermal_square_conv_Tmean.pdf, thermal_square_conv_T_profile_centerline.png/pdf")
end

main(validate_analytic=true)
