"""
Script 3: 热学网格收敛性分析 (GCI 框架)

4 组 nθ = {20, 40, 80, 160} 下运行纯热模型（均匀体积热源 + 表面冷却），
基于 Jellyroll spiral mesh，使用 opt.model = "thermal" 直接求解热方程。

指标：空间场 L2_rel, T_max(t) L2_rel, T_range(t) L2
GCI 标量：T_max(t_end)
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const THERMAL_Nθ = [20, 40, 80, 160]
const Q0 = 2.0e5  # 均匀体积热源 [W/m³]

function run_thermal_case(param_dim, nθ)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "distributed2D"
    opt.solveType = "Crank-Nicolson"
    opt.time = [0.0, 3600]
    opt.dt = [5, 5]
    opt.gsorder = 2
    opt.thermal_enabled = true
    opt.per_element_spme = false  # 纯热模型，禁用 SPMe

    # 禁用 CZM 和电化学
    opt.czm_enabled = false
    opt.mechanicalmodel = "none"

    # 创建案例
    case = JuBat.SetCase(param_dim, opt)

    # 创建 Jellyroll 热网格
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    # 添加均匀体积热源（纯热模型）
    ne = length(case.mesh["thermal2D"].element)
    q_source = fill(Q0, ne)
    case.q_source = q_source

    # 求解
    result = JuBat.Solve(case)

    # 提取结果
    t = result.time
    T = result["thermal2D temperature [K]"]
    T_max = maximum(T, dims=2)[:]
    T_range = maximum(T, dims=2)[:] .- minimum(T, dims=2)[:]

    return t, T, T_max, T_range
end

function main()
    println("=" ^ 70)
    println("Script 3: 热学网格收敛性分析 (GCI 框架)")
    println("=" ^ 70)

    # 选择 Jellyroll 电池参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    # 运行所有网格配置
    results = Dict()
    for nθ in THERMAL_Nθ
        @printf("运行 nθ = %d...\n", nθ)
        t, T, T_max, T_range = run_thermal_case(param_dim, nθ)
        results[nθ] = (t, T, T_max, T_range)
    end

    # 参考解：最细网格 (nθ=160)
    ref_nθ = maximum(THERMAL_Nθ)
    t_ref, T_ref, T_max_ref, T_range_ref = results[ref_nθ]

    # 创建输出目录
    out_dir = joinpath(@__DIR__, "..", "..", "output", "mesh_convergence")
    mkpath(out_dir)

    # 计算收敛误差
    error_table = []
    for nθ in THERMAL_Nθ
        t, T, T_max, T_range = results[nθ]

        # 时间对齐
        T_aligned = [align_to_ref(t, T[:,i], t_ref) for i in 1:size(T,2)]
        T_max_aligned = align_to_ref(t, T_max, t_ref)
        T_range_aligned = align_to_ref(t, T_range, t_ref)

        # 空间场 L2 误差（选择 t=1800s 时刻）
        t_idx = findfirst(x -> x >= 1800, t_ref)
        T_1800_ref = T_ref[t_idx, :]
        T_1800 = T_aligned[t_idx, :]

        # IDW 插值到参考网格
        mesh = results[ref_nθ][2].mesh["thermal2D"]
        x_ref = mesh.node[:,1]
        y_ref = mesh.node[:,2]

        T_1800_interp = [interpolate_to_ref_field(T_1800, t, t_ref, x_ref, y_ref; k=4)
                       for (T_1800, t) in zip(eachcol(T_1800), t_ref)]
        T_1800_interp = hcat(T_1800_interp...)

        l2_err = l2_rel_norm(T_1800_interp, T_1800_ref)
        max_err = max_norm(T_1800_interp, T_1800_ref)

        # 时间序列 L2 误差
        l2_err_Tmax = l2_rel_norm(T_max_aligned, T_max_ref)
        l2_err_Trange = l2_norm(T_range_aligned, T_range_ref)

        push!(error_table, (nθ, l2_err, max_err, l2_err_Tmax, l2_err_Trange))
    end

    # 打印误差表
    @printf("\n热学收敛误差表 (参考解: nθ=%d)\n", ref_nθ)
    @printf("%-8s %-12s %-12s %-12s %-12s\n", "nθ", "L2_rel(场)", "L∞(场)", "L2_rel(Tmax)", "L2(T_range)")
    @printf("%-8s %-12s %-12s %-12s %-12s\n", "-"^8, "-"^12, "-"^12, "-"^12, "-"^12)
    for (nθ, l2_err, max_err, l2_err_Tmax, l2_err_Trange) in error_table
        @printf("%-8d %-12.4e %-12.4e %-12.4e %-12.4e\n",
                nθ, l2_err, max_err, l2_err_Tmax, l2_err_Trange)
    end

    # GCI 分析
    # 计算等效 h（基于单元面积）
    h_values = []
    for nθ in THERMAL_Nθ
        mesh = results[nθ][2].mesh["thermal2D"]
        h = effective_h(mesh)
        push!(h_values, h)
    end

    # GCI 标量：T_max(t_end)
    T_max_end = [results[nθ][3][end] for nθ in THERMAL_Nθ]

    # 计算收敛阶和 GCI
    p, gci = compute_gci(T_max_end, h_values)

    @printf("\nGCI 汇总表 (T_max(t_end))\n")
    @printf("%-8s %-12s %-12s %-12s %-12s\n", "nθ", "h", "T_max(end)", "p", "GCI(95%)")
    @printf("%-8s %-12s %-12s %-12s %-12s\n", "-"^8, "-"^12, "-"^12, "-"^12, "-"^12)
    for i in 1:length(THERMAL_Nθ)
        @printf("%-8d %-12.4e %-12.4f %-12.4f %-12.4f\n",
                THERMAL_Nθ[i], h_values[i], T_max_end[i], p[i], gci[i])
    end

    # 渐近检查
    asymptotic_check(T_max_end, h_values)

    # 绘制收敛图
    p1 = plot(h_values, T_max_end,
             xlabel="h (等效单元尺寸)", ylabel="T_max(t_end) [K]",
             title="热学 GCI 收敛性",
             xscale=:log10, yscale=:log10,
             marker=:circle, label="数据点")

    # 添加理论收敛线
    if !isnan(p[end])
        h_ref = h_values[end]
        T_ref = T_max_end[end]
        p_val = p[end]
        h_range = [h_ref * 0.5, h_ref * 2.0]
        T_theory = T_ref .* (h_range ./ h_ref).^p_val
        plot!(h_range, T_theory, label="理论 p=$(round(p_val, digits=3))",
              linestyle=:dash, color=:red)
    end

    savefig(p1, joinpath(out_dir, "thermal_gci_convergence.png"))

    # 绘制时间序列误差
    p2 = plot(t_ref, T_max_ref, label="参考解 (nθ=$(ref_nθ))",
              xlabel="时间 [s]", ylabel="T_max [K]",
              title="最大温度时间序列")

    for (nθ, t, T_max) in zip(THERMAL_Nθ, [results[nθ][1] for nθ in THERMAL_Nθ],
                            [results[nθ][3] for nθ in THERMAL_Nθ])
        plot!(t, T_max, label="nθ=$(nθ)")
    end

    savefig(p2, joinpath(out_dir, "thermal_Tmax_time_series.png"))

    # 绘制空间场误差
    t_plot = 1800  # 选择 t=1800s 时刻
    t_idx = findfirst(x -> x >= t_plot, t_ref)

    p3 = plot(t_ref, T_max_ref, label="参考解 (nθ=$(ref_nθ))",
              xlabel="时间 [s]", ylabel="T_max [K]",
              title="最大温度时间序列")

    for (nθ, t, T_max) in zip(THERMAL_Nθ, [results[nθ][1] for nθ in THERMAL_Nθ],
                            [results[nθ][3] for nθ in THERMAL_Nθ])
        plot!(t, T_max, label="nθ=$(nθ)")
    end

    savefig(p3, joinpath(out_dir, "thermal_Tmax_time_series.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end