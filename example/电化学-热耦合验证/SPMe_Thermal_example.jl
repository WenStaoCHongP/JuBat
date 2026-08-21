"""
SPMe + 2D distributed thermal coupling verification for a jellyroll cell.

Scope:
- Electrochemical SPMe + distributed 2D thermal model only
- No mechanics, no cohesive zone model
- Mesh convergence test
- Conservation checks: current (charge), lithium inventory (mass), energy balance
- Outputs:
  - Mesh outline image
  - 60s temperature cloud
  - Mean temperature vs capacity
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat

const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "output", "SPMe_Thermal_example")
mkpath(OUTPUT_DIR)

function element_areas(mesh)
    A = zeros(Float64, size(mesh.element, 1))
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    return A
end

function plot_mesh_outline(mesh, out_path)
    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    plt = plot(size=(700, 700), aspect_ratio=1, legend=false,
        xlabel="x [m]", ylabel="y [m]", title="Jellyroll Mesh Outline")
    ne = size(mesh.element, 1)
    for e in 1:ne
        n = mesh.element[e, :]
        xs = [x[n[1]], x[n[2]], x[n[3]], x[n[4]], x[n[1]]]
        ys = [y[n[1]], y[n[2]], y[n[3]], y[n[4]], y[n[1]]]
        plot!(plt, xs, ys, lw=0.4, color=:black, alpha=0.7)
    end
    savefig(plt, out_path)
end

function effective_cp_mass(param_dim)
    cell = param_dim.cell
    area = cell.area

    m_NE = 2.0 * param_dim.NE.rho * param_dim.NE.thickness * area
    m_PE = 2.0 * param_dim.PE.rho * param_dim.PE.thickness * area
    m_SP = 2.0 * param_dim.SP.rho * param_dim.SP.thickness * area
    m_PCC = param_dim.PCC.rho * param_dim.PCC.thickness * area
    m_NCC = param_dim.NCC.rho * param_dim.NCC.thickness * area

    m_total = m_NE + m_PE + m_SP + m_PCC + m_NCC
    cp_eff = (m_NE * param_dim.NE.heat_Q +
              m_PE * param_dim.PE.heat_Q +
              m_SP * param_dim.SP.heat_Q +
              m_PCC * param_dim.PCC.heat_Q +
              m_NCC * param_dim.NCC.heat_Q) / max(m_total, 1e-12)

    return m_total, cp_eff
end

function plot_temperature_cloud(mesh, T_nodes, out_path; title="Temperature Cloud")
    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    plt = scatter(x, y, marker_z=T_nodes, markersize=3,
        xlabel="x [m]", ylabel="y [m]", title=title,
        color=:inferno, colorbar=true, aspect_ratio=1, legend=false)
    savefig(plt, out_path)
end

function mean_temperature(result, mesh)
    return result["temperature [K]"]
end

function conservation_checks(result, mesh, areas, param_dim)
    checks = Dict{String, Float64}()

    # Charge conservation: area-weighted element current vs total current
    I_e_hist = result["thermal2D element current"]
    I_total = result["cell current [A]"]
    w = areas ./ sum(areas)
    I_e_weighted = vec(sum(I_e_hist .* w; dims=1))
    err = maximum(abs.(I_e_weighted .- I_total) ./ max.(abs.(I_total), 1e-12))
    checks["charge_conservation_rel"] = err

    # Lithium inventory conservation: area-weighted soc_n + soc_p should be stable
    soc_n = result["thermal2D element soc_n"]
    soc_p = result["thermal2D element soc_p"]
    inv = vec(sum(soc_n .* w; dims=1)) + vec(sum(soc_p .* w; dims=1))
    drift = maximum(abs.(inv .- inv[1]))
    checks["lithium_inventory_drift"] = drift

    # Energy balance (approx): heat generation - convective loss vs internal energy
    q_hist = result["heat_source_fields"]
    t = result["time [s]"]
    cell = param_dim.cell
    depth = cell.width

        # total heat generation [W]
        q_total = vec(sum(q_hist .* areas; dims=1)) .* depth

        # convective loss using mean temperature
        T_mean = mean_temperature(result, mesh)
        q_loss = cell.h .* cell.cooling_surface .* (T_mean .- cell.T_amb)

        # integrate over time (trapezoid)
    function trapz(t, y)
        s = 0.0
        @inbounds for k in 2:length(t)
            s += 0.5 * (y[k] + y[k-1]) * (t[k] - t[k-1])
        end
        return s
    end

    E_gen = trapz(t, q_total)
    E_loss = trapz(t, q_loss)

    m_total, cp_eff = effective_cp_mass(param_dim)
    T0 = T_mean[1]
    dU = m_total * cp_eff * (T_mean[end] - T0)

    denom = max(abs(E_gen), 1e-12)
    checks["energy_balance_rel"] = abs((E_gen - E_loss) - dU) / denom

    return checks
end

function run_case(param_dim, crate, t_end; ntheta=80)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, t_end]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"

    opt.per_element_spme = true
    opt.mechanicalmodel = "none"
    opt.czm_enabled = false

    I1C = param_dim.cell.I1C
    I = I1C * crate
    opt.Current = x -> I

    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=ntheta, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    result = JuBat.Solve(case)
    areas = element_areas(mesh_th)
    return (; case, result, mesh=mesh_th, areas)
end

function mesh_convergence(param_dim, t_end)
    ntheta_list = [40, 60, 80, 120]
    @printf("\n[Mesh convergence] 1C, t_end=%.1f s\n", t_end)
    @printf("  ntheta | T_mean_end [K] | V_end [V]\n")
    for nθ in ntheta_list
        run = run_case(param_dim, 1.0, t_end; ntheta=nθ)
        T_mean = mean_temperature(run.result, run.mesh)
        V = run.result["cell voltage [V]"]
        @printf("  %6d | %13.4f | %9.4f\n", nθ, T_mean[end], V[end])
    end
end

function main()
    param_dim = JuBat.ChooseCell("Jellyroll")

    # Requested C-rates and durations
    cases = [
        (crate=0.3, t_end=12000.0),
        (crate=0.7, t_end=4800.0),
        (crate=1.0, t_end=3600.0)
    ]

    # Base mesh for visualization
    ntheta_base = 80

    # Run convergence test (1C)
    mesh_convergence(param_dim, 3600.0)

    # Run each case
    T_mean_curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for c in cases
        run = run_case(param_dim, c.crate, c.t_end; ntheta=ntheta_base)
        result = run.result
        mesh = run.mesh
        areas = run.areas

        # Mesh outline (once)
        if c == cases[1]
            plot_mesh_outline(mesh, joinpath(OUTPUT_DIR, "spme_thermal_mesh.png"))
        end

        # 60 s temperature cloud
        T_nodes_hist = result["thermal2D T_nodes [K]"]
        t = result["time [s]"]
        _, idx = findmin(abs.(t .- 60.0))
        T_nodes_60 = ndims(T_nodes_hist) == 2 ? T_nodes_hist[:, idx] : T_nodes_hist
        out_path = joinpath(OUTPUT_DIR, @sprintf("spme_thermal_T60_%.1fC.png", c.crate))
        plot_temperature_cloud(mesh, T_nodes_60, out_path;
            title=@sprintf("T at 60 s (%.1fC)", c.crate))

        # Mean temperature vs capacity
        T_mean = mean_temperature(result, mesh)
        t = result["time [s]"]
        I = abs(result["cell current [A]"][1])
        capacity_Ah = I .* t ./ 3600.0
        T_mean_curves[c.crate] = (capacity_Ah, T_mean)

        # Conservation checks
        checks = conservation_checks(result, mesh, areas, param_dim)
        println("\n[Checks] ", @sprintf("%.1fC", c.crate))
        @printf("  Charge conservation rel error: %.3e\n", checks["charge_conservation_rel"])
        @printf("  Lithium inventory drift: %.3e\n", checks["lithium_inventory_drift"])
        @printf("  Energy balance rel error: %.3e\n", checks["energy_balance_rel"])
    end

    # Plot mean temperature vs capacity
    p = plot(xlabel="Capacity [Ah]", ylabel="Mean temperature [K]",
        title="Mean Temperature vs Capacity")
    for (crate, (cap, Tmean)) in sort(collect(T_mean_curves))
        plot!(p, cap, Tmean, lw=2, label=@sprintf("%.1fC", crate))
    end
    savefig(p, joinpath(OUTPUT_DIR, "spme_thermal_Tmean_vs_capacity.png"))

    println("\nOutputs:")
    println("  - ", joinpath(OUTPUT_DIR, "spme_thermal_mesh.png"))
    println("  - ", joinpath(OUTPUT_DIR, "spme_thermal_T60_0.3C.png"))
    println("  - ", joinpath(OUTPUT_DIR, "spme_thermal_T60_0.7C.png"))
    println("  - ", joinpath(OUTPUT_DIR, "spme_thermal_T60_1.0C.png"))
    println("  - ", joinpath(OUTPUT_DIR, "spme_thermal_Tmean_vs_capacity.png"))
end

main()
