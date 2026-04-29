"""
Script 3: 热网格敏感性分析

4 组 nθ = {20, 40, 80, 160} 下运行纯热模型（均匀体积热源 + 表面冷却），
基于 Jellyroll spiral mesh，使用 opt.model = "thermal" 直接求解热方程。

指标：T_max、T_min、T_range、角度变化收敛
输出：收敛图
"""

using Printf, Plots, Statistics

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const THERMAL_Nθ = [20, 40, 80, 160]
const Q0 = 2.0e5  # 均匀体积热源 [W/m³]

function run_thermal_case(param_dim, nθ)
    # 先用 SPMe 通过 SetCase（thermal 模式不支持 SetCase）
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "distributed2D"
    opt.solveType = "Crank-Nicolson"
    opt.time = [0.0, 3600]
    opt.dt = [5, 5]
    opt.gsorder = 2
    opt.Nn = 5; opt.Ns = 3; opt.Np = 5
    opt.Nrn = 5; opt.Nrp = 5

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    # 切换为纯热模式
    case.opt.model = "thermal"

    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    nnode = mesh.nlen
    scale = param_dim.scale

    # 归一化均匀热源
    q0_nd = Q0 / scale.q

    # 初始化热变量
    variables = Dict{String,Any}()
    variables["thermal2D temperature at nodes"] = fill(param_dim.cell.T0 / scale.T_ref, nnode)
    variables["heat_source_fields"] = fill(q0_nd, ne)

    # 更新函数（热源恒定）
    update_fn = (t, vars) -> begin
        vars["heat_source_fields"] = fill(q0_nd, ne)
    end

    solve_t0 = time_ns()
    result = JuBat.Solve(case;
        thermal_variables=variables,
        thermal_update_fn=update_fn,
        thermal_record=true)
    solve_wall = (time_ns() - solve_t0) * 1e-9

    T_nodes = result.T_nodes .* scale.T_ref
    times   = result.time .* scale.t0
    T_hist  = result.T_hist .* scale.T_ref

    return (T_nodes=T_nodes, times=times, T_hist=T_hist, mesh=mesh, solve_wall=solve_wall)
end

function main()
    println("=" ^ 70)
    println("Script 3: 热网格敏感性分析 (纯热模型, q0=$(Q0) W/m³)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")

    all_Tmax   = Float64[]
    all_Tmin   = Float64[]
    all_Trange = Float64[]
    all_wall   = Float64[]
    all_results = Dict{Int, Any}()

    for nθ in THERMAL_Nθ
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        r = run_thermal_case(param_dim, nθ)
        dt_solve = r.solve_wall

        all_results[nθ] = r

        T_final = r.T_nodes
        Tmax = maximum(T_final)
        Tmin = minimum(T_final)
        Trange = Tmax - Tmin

        push!(all_Tmax, Tmax)
        push!(all_Tmin, Tmin)
        push!(all_Trange, Trange)
        push!(all_wall, dt_solve)

        @printf("  耗时 %.2f s,  T_max = %.3f K,  T_min = %.3f K,  ΔT = %.4f K\n",
            dt_solve, Tmax, Tmin, Trange)
    end

    # ── 参考解 (nθ=160) ──
    ref_idx = length(THERMAL_Nθ)

    println("\n" * "=" ^ 70)
    println("热网格收敛性汇总")
    println("=" ^ 70)
    @printf("%-10s  %10s  %10s  %10s  %12s  %10s\n",
            "nθ", "T_max [K]", "T_min [K]", "ΔT [K]", "T_max err [%]", "Solve [s]")
    println("-" ^ 70)

    T_errors = Float64[]
    for (i, nθ) in enumerate(THERMAL_Nθ)
        err = abs(all_Tmax[i] - all_Tmax[ref_idx]) / all_Tmax[ref_idx] * 100
        push!(T_errors, err)
        @printf("%-10d  %10.3f  %10.3f  %10.4f  %12.4f  %10.2f\n",
                nθ, all_Tmax[i], all_Tmin[i], all_Trange[i], err, all_wall[i])
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_sensitivity")
    mkpath(out_dir)

    colors = [:red, :orange, :green, :blue]

    # 图1: 温度场时间演化（取体积平均温度）
    p1 = plot(xlabel="Time [s]", ylabel="Mean Temperature [K]",
              title="Thermal Mesh: Temperature (pure thermal, q0=$(Q0) W/m³)",
              legend=:bottomright)
    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        if !isempty(T_hist)
            nnode = size(T_hist, 1)
            nt = size(T_hist, 2)
            T_mean = [mean(T_hist[:, k]) for k in 1:nt]
            label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
            ls = i == ref_idx ? :solid : :dash
            plot!(p1, r.times[1:nt], T_mean, label=label, lw=1.5, color=colors[i], ls=ls)
        end
    end
    savefig(p1, joinpath(out_dir, "thermal_temperature_convergence.png"))

    # 图2: T_max / T_min 收敛
    p2 = plot(THERMAL_Nθ, all_Tmax,
              xlabel="nθ (per revolution)", ylabel="Temperature [K]",
              title="Thermal Mesh: Peak Temperature Convergence",
              marker=:o, lw=2, color=:blue, label="T_max")
    plot!(p2, THERMAL_Nθ, all_Tmin,
          marker=:s, lw=2, color=:red, ls=:dash, label="T_min")
    savefig(p2, joinpath(out_dir, "thermal_peak_convergence.png"))

    # 图3: 相对误差折线图
    xp = collect(1:length(THERMAL_Nθ))
    nθ_labels = ["$n" for n in THERMAL_Nθ]
    p3 = plot(xlabel="Mesh No.", ylabel="T_max Relative Error [%]",
              title="Thermal Mesh: Convergence Error",
              xticks=(xp, nθ_labels), legend=:topright)
    plot!(p3, xp, T_errors, marker=:o, lw=2, color=:blue, label="|T_max - T_ref| / T_ref")
    hline!([5.0], label="5%", color=:red, ls=:dash)
    hline!([1.0], label="1%", color=:green, ls=:dash)
    # 网格说明文本框
    mesh_txt = "1→nθ=20\n2→nθ=40\n3→nθ=80\n4→nθ=160 ref"
    annotate!(p3, maximum(xp) - 0.1, maximum(T_errors) * 0.9,
              text(mesh_txt, 8, :right, :top, :Courier))
    savefig(p3, joinpath(out_dir, "thermal_convergence_error.png"))

    # 图4: ΔT 收敛（空间均匀性指标）
    p4 = plot(THERMAL_Nθ, all_Trange,
              xlabel="nθ (per revolution)", ylabel="ΔT = T_max - T_min [K]",
              title="Thermal Mesh: Spatial Uniformity",
              marker=:diamond, lw=2, color=:purple, label="ΔT")
    savefig(p4, joinpath(out_dir, "thermal_spatial_uniformity.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
