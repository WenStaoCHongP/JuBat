"""
测试案例：Jellyroll电池多SPMe并行电化学-热耦合仿真

功能：
- 多SPMe并行架构 + 二维分布式热模型 + CZM内聚力模型
- 输出各模块仿真时间占比及总耗时
- 温度及电压随容量变化曲线
- 最大/平均分离位移随时间变化曲线
- 损伤D及容量随时间变化曲线

日期：2025-12-31
"""

using Printf
using Plots
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("Jellyroll电池多SPMe并行电化学-热耦合仿真")
    println("="^80)

    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/4] 参数设置...")

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = JuBat.Option()
    Crates = 1.0
    i = 5 * Crates
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"

    opt.time = [0.0, 60]
    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.debug_coupling = true
    opt.debug_log_path = joinpath(@__DIR__, "..", "output", "testexample", "simple_coupling_debug.log")
    opt.czm_enabled = true
    opt.czm_fix_inner = false
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    println("OK: 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 多SPMe并行\n")
    @printf("  CZM迭代法: %s\n", opt.czm_iter_method)
    @printf("  CZM载荷子步数: %d\n", opt.czm_load_steps)

    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/4] 创建案例和Jellyroll网格...")

    case = JuBat.SetCase(param_dim, opt)

    n_theta = 80
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # 创建 CZM 网格（czm_enabled = true 时必须）
    if opt.czm_enabled
        case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    end

    ne = size(mesh_th.element, 1)
    println("OK: Jellyroll网格创建完成")
    @printf("  周向单元数 n_theta: %d\n", n_theta)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", mesh_th.nlen)

    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)

    # ========================================================================
    # 3. 运行求解器
    # ========================================================================
    println("\n[3/4] 运行多SPMe并行求解器...")

    t_wall_start = time_ns()
    result = JuBat.Solve(case)
    t_wall_s = (time_ns() - t_wall_start) * 1e-9

    println("OK: 求解成功完成")

    # ========================================================================
    # 4. 输出各模块仿真时间占比
    # ========================================================================
    println("\n[4/4] 仿真耗时统计")
    println("-"^60)

    num_steps = length(result["time [s]"])
    V = result["cell voltage [V]"]

    @printf("  总时间步数: %d\n", num_steps)
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    @printf("  电压降: %.4f V\n", V[1] - V[end])

    @printf("\n  总仿真耗时 (wall-clock): %.3f s\n", t_wall_s)

    if haskey(result, "timing SPMe solve total [s]")
        println("\n  各模块耗时分解（CallModel累计）:")
        @printf("    %-25s %10s %8s %15s\n", "模块", "累计 [s]", "占比 [%]", "平均 [ms/step]")
        println("    " * "-"^62)

        modules = [
            ("SPMe 求解",             "timing SPMe solve total [s]",           "timing SPMe solve ratio [%]",           "timing SPMe solve avg [ms]"),
            ("分流求解器",             "timing branch solver total [s]",        "timing branch solver ratio [%]",        "timing branch solver avg [ms]"),
            ("热分布式模型",           "timing thermal distributed total [s]",  "timing thermal distributed ratio [%]",  "timing thermal distributed avg [ms]"),
            ("CZM 模型",              "timing CZM model total [s]",            "timing CZM model ratio [%]",            "timing CZM model avg [ms]"),
        ]

        total_cm = 0.0
        for (name, tot_key, ratio_key, avg_key) in modules
            tot = get(result, tot_key, 0.0)
            ratio = get(result, ratio_key, 0.0)
            avg = get(result, avg_key, 0.0)
            total_cm += tot
            @printf("    %-25s %10.3f %8.2f %15.3f\n", name, tot, ratio, avg)
        end
        println("    " * "-"^62)
        @printf("    %-25s %10.3f\n", "CallModel 合计", total_cm)
    end

    println("\n" * "="^80)
    println("仿真完成")
    println("="^80)

    # ========================================================================
    # 5. 结果输出
    # ========================================================================
    println("\n[结果输出]")

    t_s = result["time [s]"]
    V_V = result["cell voltage [V]"]
    I_A = result["cell current [A]"]
    T_K = result["temperature [K]"]

    # 计算累积放电容量 (Ah)
    capacity_Ah = zeros(Float64, length(t_s))
    for k in 2:length(t_s)
        dt_k = t_s[k] - t_s[k-1]
        capacity_Ah[k] = capacity_Ah[k-1] + abs(I_A[k]) * dt_k / 3600.0
    end

    @printf("  最终容量: %.4f Ah\n", capacity_Ah[end])
    @printf("  温度范围: %.2f K ~ %.2f K\n", minimum(T_K), maximum(T_K))

    # CZM 结果输出
    has_czm = haskey(result, "czm D_max")
    if has_czm
        D_max = result["czm D_max"]
        D_mean = result["czm D_mean"]
        δ_max = result["czm δ_max_n [m]"]
        δ_mean = result["czm δ_mean_n [m]"]
        n_frac = result["czm n_fractured"]
        @printf("  CZM 最终 D_max: %.4f%%\n", D_max[end] * 100)
        @printf("  CZM 最终 D_mean: %.4f%%\n", D_mean[end] * 100)
        @printf("  CZM 最大法向分离位移: %.4e m\n", maximum(δ_max))
        @printf("  CZM 断裂单元数: %d\n", Int(n_frac[end]))
    end

    # ========================================================================
    # 6. 绘图
    # ========================================================================
    println("\n[绘图]")

    output_dir = joinpath(@__DIR__, "..", "output", "testexample")
    mkpath(output_dir)

    # ── 图1: 温度及电压随容量变化 ──
    p1_T = plot(capacity_Ah, T_K,
        xlabel="Capacity (Ah)", ylabel="Temperature (K)",
        label="T", lw=2, color=:red, legend=:topleft,
        title="Temperature vs Capacity")
    p1_V = twinx(p1_T)
    plot!(p1_V, capacity_Ah, V_V,
        ylabel="Voltage (V)", label="V", lw=2, color=:blue,
        legend=:topright)

    # ── 图2: 最大/平均分离位移随时间变化 ──
    if has_czm
        p2 = plot(t_s, δ_max .* 1e9,
            xlabel="Time (s)", ylabel="Separation Displacement (nm)",
            label="δ_max", lw=2, color=:red, marker=:none)
        plot!(p2, t_s, δ_mean .* 1e9,
            label="δ_mean", lw=2, color=:blue, ls=:dash)
        title!(p2, "Normal Separation vs Time")
    else
        p2 = plot(title="CZM not enabled", grid=false)
    end

    # ── 图3: 损伤D及容量随时间变化 ──
    if has_czm
        p3 = plot(t_s, D_max .* 100,
            xlabel="Time (s)", ylabel="Damage (%)",
            label="D_max", lw=2, color=:red, legend=:topleft,
            title="Damage & Capacity vs Time")
        plot!(p3, t_s, D_mean .* 100,
            label="D_mean", lw=2, color=:blue, ls=:dash)
        p3_C = twinx(p3)
        plot!(p3_C, t_s, capacity_Ah,
            ylabel="Capacity (Ah)", label="Capacity", lw=2, color=:green,
            legend=:topright)
    else
        p3 = plot(title="CZM not enabled", grid=false)
    end

    # ── 组合输出 ──
    p_combined = plot(p1_T, p2, p3, layout=(3, 1), size=(800, 900))
    savefig(p_combined, joinpath(output_dir, "testexample_results.png"))
    @printf("  结果图已保存: %s\n", joinpath(output_dir, "testexample_results.png"))

    println("\n" * "="^80)
    println("全部完成")
    println("="^80)
end

main()
