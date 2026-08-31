"""
快速文字基线：Jellyroll电池多SPMe并行电化学-热-CZM仿真（仅文字结果，60 秒）

功能：
- 多SPMe并行架构 + 二维分布式热模型 + CZM内聚力模型
- 输出各模块仿真时间占比及总耗时
- 输出全部求解文字指标（电压/容量/温度/CZM/层分辨应力范围），作为修改后的快速行为基线

说明：本文件不含绘图（AGENTS.md §9.6 快速门）；全套绘图验证入口见
`example/couple_example.jl`（输出 `output/couple_example/`）。

日期：2025-12-31（2026-08-29 拆分为纯文字基线）
"""

using Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function rotate_stress_to_polar(node, element, sigma_xx, sigma_yy, sigma_xy)
    ne = size(element, 1)
    length(sigma_xx) == length(sigma_yy) == length(sigma_xy) == ne ||
        error("stress component lengths must match the mechanical bulk element count")

    sigma_theta_theta = Vector{Float64}(undef, ne)
    tau_r_theta = Vector{Float64}(undef, ne)
    @inbounds for e in 1:ne
        nodes = element[e, :]
        x = sum(node[nodes, 1]) / 4
        y = sum(node[nodes, 2]) / 4
        r2 = x^2 + y^2
        r2 > 0 || error("mechanical bulk element $e has its centroid at the polar origin")

        sigma_theta_theta[e] = (
            y^2 * sigma_xx[e] - 2x * y * sigma_xy[e] + x^2 * sigma_yy[e]
        ) / r2
        tau_r_theta[e] = (
            x * y * (sigma_yy[e] - sigma_xx[e]) + (x^2 - y^2) * sigma_xy[e]
        ) / r2
    end
    return sigma_theta_theta, tau_r_theta
end

function main()
    println("="^80)
    println("Jellyroll电池多SPMe并行电化学-热耦合仿真（文字基线）")
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
    opt.czm.enabled = true
    opt.czm.model = "mix"
    opt.czm.fix_inner = false
    opt.czm.iter_method = "basic"
    opt.czm.load_steps = 10
    opt.czm.tol = 1e-3

    println("OK: 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 多SPMe并行\n")
    @printf("  CZM迭代法: %s\n", opt.czm.iter_method)
    @printf("  CZM载荷子步数: %d\n", opt.czm.load_steps)

    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/4] 创建案例和Jellyroll网格...")

    case = JuBat.SetCase(param_dim, opt)

    n_theta = 80
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # 创建 CZM 网格（czm_enabled = true 时必须）与演化状态
    if opt.czm.enabled
        case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
        case.mech = JuBat.MechState(case.czm_mesh)
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
    # 5. 结果输出（仅文字）
    # ========================================================================
    println("\n[结果输出]")

    t_s = result["time [s]"]
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
        n_frac = result["czm n_fractured"]
        @printf("  CZM 最终 D_max: %.4f%%\n", D_max[end] * 100)
        @printf("  CZM 最终 D_mean: %.4f%%\n", D_mean[end] * 100)
        @printf("  CZM 最大法向分离位移: %.4e m\n", maximum(δ_max))
        @printf("  CZM 断裂单元数: %d\n", Int(n_frac[end]))
    end

    # 层分辨应力文字指标（求解中在线导出）
    sigma_xx_MPa = result["diffusion stress xx [Pa]"][:, end] .* 1e-6
    sigma_yy_MPa = result["diffusion stress yy [Pa]"][:, end] .* 1e-6
    sigma_xy_MPa = result["diffusion stress xy [Pa]"][:, end] .* 1e-6
    sigma_theta_theta_MPa, tau_r_theta_MPa = rotate_stress_to_polar(
        case.czm_mesh.node,
        case.czm_mesh.bulk_element,
        sigma_xx_MPa, sigma_yy_MPa, sigma_xy_MPa)

    @printf("  最终环向应力范围: %.4e ~ %.4e MPa\n",
        minimum(sigma_theta_theta_MPa), maximum(sigma_theta_theta_MPa))
    @printf("  最终切向剪应力范围: %.4e ~ %.4e MPa\n",
        minimum(tau_r_theta_MPa), maximum(tau_r_theta_MPa))

    println("\n" * "="^80)
    println("全部完成")
    println("="^80)
end

main()
