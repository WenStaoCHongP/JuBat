"""
示例：Jellyroll_very
- 仅保留电化学-热耦合（多SPMe并行）
- 绘制不同倍率下电压/温度曲线，并与 PyBaMM 对比
- 优化最终温度场云图（低温浅蓝，高温红）
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf, CSV, DataFrames
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function run_case(Crate; nθ=80)
    # =====================
    # 1. 参数设置
    # =====================
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2

    # 电流与时间
    I1C = 5.0
    Iapp = I1C * Crate
    opt.Current = x -> Iapp
    opt.time = [0.0, 3600 / Crate]
    opt.dt = [0.5, 10.0]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    # 热模型
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"

    # 多SPMe并行
    opt.per_element_spme = true

    case = JuBat.SetCase(param_dim, opt)

    # Jellyroll collector-seeded 网格
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th

    result = JuBat.Solve(case)
    return case, mesh_th, result
end

function plot_temperature_field(case, mesh_th, result; outfile="Jellyroll_very_Tfield.png")
    if !haskey(result, "thermal2D T_nodes [K]")
        @warn "未找到温度场数据"
        return
    end

    T_nodes_final = result["thermal2D T_nodes [K]"]
    xnod = mesh_th.node[:, 1]
    ynod = mesh_th.node[:, 2]

    # 插值网格
    nx, ny = 400, 400
    xs = range(minimum(xnod), stop=maximum(xnod), length=nx)
    ys = range(minimum(ynod), stop=maximum(ynod), length=ny)

    dx = step(xs); dy = step(ys)
    sigma = 1.0 * max(dx, dy)
    two_sigma2 = 2.0 * sigma^2

    Z = fill(NaN, ny, nx)

    Rin = getfield(case.param_dim.cell, :Rin)
    Rout = getfield(case.param_dim.cell, :Rout)

    @inbounds for j in 1:ny
        yv = ys[j]
        for i in 1:nx
            xv = xs[i]
            r = sqrt(xv^2 + yv^2)
            if r < Rin || r > Rout
                continue
            end

            dxv = xnod .- xv
            dyv = ynod .- yv
            d2 = dxv .* dxv .+ dyv .* dyv
            w = exp.(-d2 ./ two_sigma2)
            s = sum(w)
            if s > 0
                Z[j, i] = sum(w .* T_nodes_final) / s
            end
        end
    end

    valid = .!isnan.(Z)
    if !any(valid)
        @warn "无有效温度数据可视化"
        return
    end

    Zvals = Z[valid]
    vmin = minimum(Zvals)
    vmax = maximum(Zvals)

    # 低温浅蓝，高温红
    cmap = cgrad([:lightblue, :red])

    p = plot(size=(800, 800), title="Final Temperature Field")
    heatmap!(p, xs, ys, Z;
        aspect_ratio=1,
        color=cmap,
        colorbar=true,
        xlabel="x (m)",
        ylabel="y (m)",
        clims=(vmin, vmax))
    contour!(p, xs, ys, Z; levels=10, linewidth=1, linecolor=:black, alpha=0.4)
    savefig(p, outfile)
end

function main()
    println("="^80)
    println("Jellyroll_very：电化学-热耦合仿真")
    println("="^80)

    # 倍率设置
    Crates = [0.5, 1.0, 2.0]
    colors = [:black, :blue, :red]

    pV = plot(xlabel="Output capacity [Ah]", ylabel="Cell voltage [V]")
    ylims!(pV, 2.5, 4.3)
    pT = plot(xlabel="Output capacity [Ah]", ylabel="Temperature [K]", legend=:bottomright)

    path = joinpath(pwd(), "src", "data")
    baseline_csv = joinpath(path, "pybamm_SPMe_1C.csv")
    baseline_data = nothing
    try
        baseline_data = CSV.read(baseline_csv, DataFrame, header=1)
        baseline_data = Matrix(baseline_data)
        println("✓ 已加载基线: pybamm_SPMe_1C.csv")
    catch e
        @warn "基线数据读取失败" exception=(e, catch_backtrace())
    end

    for i in eachindex(Crates)
        Crate = Crates[i]
        case, mesh_th, result = run_case(Crate)

        # 仿真曲线用散点
        t = result["time [s]"]
        V = result["cell voltage [V]"]
        T = result["temperature [K]"]

        cap = t ./ 3600 .* Crate * 5

        scatter!(pV, cap, V; label="$(Crate)C (JuBat)", color=colors[i], markersize=3)
        scatter!(pT, cap, T; label="$(Crate)C (JuBat)", color=colors[i], markersize=3)

        # PyBaMM 基线对比曲线（统一 1C 参考）
        if baseline_data !== nothing
            plot!(pV, baseline_data[:, 1] ./ 3600 .* 1.0 * 5, baseline_data[:, 2],
                label="PyBaMM SPMe 1C (baseline)", linestyle=:dot, linecolor=:gray, lw=2)
            # 基线文件不含温度列，仅绘制电压
        end

        # 仅对最后一组倍率绘制温度场
        if i == lastindex(Crates)
            plot_temperature_field(case, mesh_th, result; outfile="Jellyroll_very_Tfield.png")
        end
    end

    savefig(pV, "Jellyroll_very_voltage.pdf")
    savefig(pT, "Jellyroll_very_temperature.pdf")

    println("✓ 完成：Jellyroll_very_voltage.pdf")
    println("✓ 完成：Jellyroll_very_temperature.pdf")
    println("✓ 完成：Jellyroll_very_Tfield.png")
end

main()
