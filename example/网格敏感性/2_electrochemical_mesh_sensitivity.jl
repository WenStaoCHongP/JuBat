"""
Script 2: 电化学网格敏感性分析

4 组 (Nn, Ns, Np) 下运行 SPMe + lumped thermal：
  - (40, 20, 40)  ← 参考解
  - (20, 10, 20)
  - (20,  5, 20)
  - (10,  5, 10)

指标：V(t) RMSPE、T(t) RMSPE、dT/dt(t) RMSPE
输出：收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_rmspe_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const MESH_CONFIGS = [
    (40, 20, 40),
    (20, 10, 20),
    (20,  5, 20),
    (10,  5, 10),
]

function run_echem_case(param_dim, Nn, Ns, Np)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "lumped"
    opt.dimension = 1
    opt.Nn = Nn; opt.Ns = Ns; opt.Np = Np
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.per_element_spme = false
    opt.mechanicalmodel = "none"
    opt.czm.enabled = false

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)

    result = JuBat.Solve(case)
    return result
end

function main()
    println("=" ^ 70)
    println("Script 2: 电化学网格敏感性分析")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    results = Dict{Tuple{Int,Int,Int}, Dict}()
    wall_times = Float64[]

    for (Nn, Ns, Np) in MESH_CONFIGS
        label = "($Nn, $Ns, $Np)"
        @printf("\n--- 运行 %s ---\n", label)
        t0 = time_ns()
        r = run_echem_case(param_dim, Nn, Ns, Np)
        dt_wall = (time_ns() - t0) * 1e-9
        push!(wall_times, dt_wall)
        results[(Nn, Ns, Np)] = r
        @printf("  耗时 %.2f s,  %d 步\n", dt_wall, length(r["time [s]"]))
    end

    # ── 参考解 ──
    ref = results[MESH_CONFIGS[1]]
    t_ref = ref["time [s]"]
    V_ref = ref["cell voltage [V]"]
    T_ref = ref["temperature [K]"]

    # ── 汇总表 (RMSPE) ──
    println("\n" * "=" ^ 70)
    println("电化学网格收敛性汇总 (RMSPE)")
    println("=" ^ 70)
    @printf("%-16s  %12s  %12s  %12s  %10s  %10s  %10s  %10s\n",
            "Mesh", "V(t) RMSPE%", "T(t) RMSPE%", "dT/dt RMSPE%",
            "V_skip%", "T_skip%", "dT_skip%", "Wall [s]")
    println("-" ^ 100)

    V_errors = Float64[]
    T_errors = Float64[]
    dTdt_errors = Float64[]

    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        t_c = r["time [s]"]
        V_c = r["cell voltage [V]"]
        T_c = r["temperature [K]"]

        # 对齐到参考时间网格
        V_aligned = align_to_ref(t_c, V_c, t_ref)
        T_aligned = align_to_ref(t_c, T_c, t_ref)

        # 电压曲线 RMSPE
        err_V, skip_V = i == 1 ? (0.0, 0.0) : rmspe(V_aligned, V_ref)

        # 温度曲线 RMSPE
        err_T, skip_T = i == 1 ? (0.0, 0.0) : rmspe(T_aligned, T_ref)

        # dT/dt 曲线 RMSPE
        if i == 1
            err_dTdt, skip_dTdt = (0.0, 0.0)
        else
            dTdt_ref_vals = diff(T_ref) ./ diff(t_ref)
            dTdt_cand_vals = diff(T_aligned) ./ diff(t_ref)
            err_dTdt, skip_dTdt = rmspe(dTdt_cand_vals, dTdt_ref_vals)
        end

        push!(V_errors, err_V)
        push!(T_errors, err_T)
        push!(dTdt_errors, err_dTdt)

        @printf("(%2d,%2d,%2d)       %12.4f  %12.4f  %12.4f  %10.1f  %10.1f  %10.1f  %10.2f\n",
                Nn, Ns, Np, err_V, err_T, err_dTdt,
                skip_V*100, skip_T*100, skip_dTdt*100, wall_times[i])
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "2_electrochemical_mesh_sensitivity")
    mkpath(out_dir)

    colors = [:blue, :orange, :green, :red]

    # 图1: 电压曲线对比
    p1 = plot(xlabel="Time [s]", ylabel="Cell Voltage [V]",
              title="Electrochemical Mesh: Voltage", legend=:topright)
    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        label = i == 1 ? "ref ($Nn,$Ns,$Np)" : "($Nn,$Ns,$Np)"
        ls = i == 1 ? :solid : :dash
        plot!(p1, r["time [s]"], r["cell voltage [V]"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p1, joinpath(out_dir, "echem_voltage_convergence.png"))

    # 图2: 温度曲线对比
    p2 = plot(xlabel="Time [s]", ylabel="Temperature [K]",
              title="Electrochemical Mesh: Temperature", legend=:bottomright)
    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        label = i == 1 ? "ref ($Nn,$Ns,$Np)" : "($Nn,$Ns,$Np)"
        ls = i == 1 ? :solid : :dash
        plot!(p2, r["time [s]"], r["temperature [K]"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p2, joinpath(out_dir, "echem_temperature_convergence.png"))

    # 图3: 收敛误差折线图 (RMSPE)
    xp = collect(1:4)
    p3 = plot(xlabel="Mesh No.", ylabel="RMSPE [%]",
              title="Electrochemical Mesh: Convergence Error (RMSPE)",
              xticks=xp, xlims=(0.5, 4.5),
              legend=:topright)
    plot!(p3, xp, V_errors, marker=:o, lw=2, label="V(t) RMSPE", color=:blue)
    plot!(p3, xp, T_errors, marker=:s, lw=2, label="T(t) RMSPE", color=:orange)
    plot!(p3, xp, dTdt_errors, marker=:d, lw=2, label="dT/dt RMSPE", color=:green)
    hline!([5.0], label="5% threshold", color=:black, ls=:dash, lw=2)
    hline!([1.0], label="1%", color=:gray, ls=:dot, lw=1)
    mesh_txt = "1→(40,20,40) ref\n2→(20,10,20)\n3→(20,5,20)\n4→(10,5,10)"
    annotate!(p3, 4.4, maximum([V_errors; T_errors; dTdt_errors]) * 0.95,
              text(mesh_txt, 8, :right, :top, :Courier))
    savefig(p3, joinpath(out_dir, "echem_convergence_error.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()
