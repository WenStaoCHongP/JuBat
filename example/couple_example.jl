"""
完整耦合示例：Jellyroll电池多SPMe并行电化学-热-CZM仿真（求解 + 全套绘图）

功能：
- 多SPMe并行架构 + 二维分布式热模型 + CZM内聚力模型
- 输出各模块仿真时间占比及总耗时
- 输出最终时刻温度场
- 输出最终时刻层分辨环向应力场和切向剪应力场（求解中在线导出）

说明：本文件为全量验证入口（含绘图产物）；日常修改的快速文字基线见
`example/testexample.jl`（60 秒，仅文字结果）。输出按 AGENTS.md §9.9 写入
`output/couple_example/`。

日期：2025-12-31（2026-08-29 自 testexample.jl 拆分）
"""

using Printf
using Statistics
using Plots
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function q4_nodal_to_element_mean(mesh, nodal_field)
    length(nodal_field) == mesh.nlen ||
        error("nodal field length must match the thermal node count")
    ne = size(mesh.element, 1)
    field_elem = Vector{Float64}(undef, ne)
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        field_elem[e] = sum(nodal_field[nodes]) / 4
    end
    return field_elem
end

function plot_q4_cloud(node, element, field_elem, coordinate_scale;
                       title, colorbar_title, cmap, clims=nothing)
    ne = size(element, 1)
    length(field_elem) == ne || error("field length must match the plotted Q4 element count")
    all(isfinite, field_elem) || error("spatial field values contain NaN or Inf")

    plt = plot(
        clims=clims,
        colorbar=true,
        colorbar_title=colorbar_title,
        xlabel="x [m]",
        ylabel="y [m]",
        title=title,
        aspect_ratio=:equal,
        legend=false,
        size=(800, 700))

    @inbounds for e in 1:ne
        nodes = element[e, :]
        x = node[nodes, 1] .* coordinate_scale
        y = node[nodes, 2] .* coordinate_scale
        all(isfinite, x) || error("Q4 element $e has NaN or Inf x coordinates")
        all(isfinite, y) || error("Q4 element $e has NaN or Inf y coordinates")
        plot!(plt, Plots.Shape(x, y);
            fill_z=[field_elem[e]],
            color=cmap,
            clims=clims,
            linealpha=0.0,
            linewidth=0.0,
            colorbar_entry=true,
            label=false)
    end
    return plt
end

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
    opt.debug_log_path = joinpath(@__DIR__, "..", "output", "couple_example", "simple_coupling_debug.log")
    opt.czm.enabled = true
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

    # 创建 CZM 网格与演化状态（czm_enabled = true 时必须）
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
    # 5. 结果输出
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
        δ_key = JuBat.czm_max_separation_key(opt.czm.model)
        δ_max = result[δ_key]
        n_frac = result["czm n_fractured"]
        @printf("  CZM 最终 D_max: %.4f%%\n", D_max[end] * 100)
        @printf("  CZM 最终 D_mean: %.4f%%\n", D_mean[end] * 100)
        if opt.czm.model == "mix"
            @printf("  CZM 最大有效分离位移: %.4e m\n", maximum(δ_max))
        else
            @printf("  CZM 最大法向分离位移: %.4e m\n", maximum(δ_max))
        end
        @printf("  CZM 断裂单元数: %d\n", Int(n_frac[end]))
    end

    # ========================================================================
    # 6. 绘图
    # ========================================================================
    println("\n[绘图]")

    output_dir = joinpath(@__DIR__, "..", "output", "couple_example")
    mkpath(output_dir)

    T_nodes_final_K = result["thermal2D final temperature at nodes [K]"]

    # 层分辨应力直接取自求解过程中在线导出的结果键（Pa → MPa）
    sigma_xx_MPa = result["diffusion stress xx [Pa]"][:, end] .* 1e-6
    sigma_yy_MPa = result["diffusion stress yy [Pa]"][:, end] .* 1e-6
    sigma_xy_MPa = result["diffusion stress xy [Pa]"][:, end] .* 1e-6

    mechanical_mesh = case.czm_mesh
    sigma_theta_theta_MPa, tau_r_theta_MPa = rotate_stress_to_polar(
        mechanical_mesh.node,
        mechanical_mesh.bulk_element,
        sigma_xx_MPa,
        sigma_yy_MPa,
        sigma_xy_MPa)

    coordinate_scale = case.param.scale.L
    T_elem_final_K = q4_nodal_to_element_mean(mesh_th, T_nodes_final_K)
    final_time_s = t_s[end]

    # 对称色标取 |场| 的 99% 分位：端部应力集中不压扁涂层尺度的层间对比
    symmetric_p99_clims(field) = begin
        q = quantile(abs.(field), 0.99)
        q > 0 ? (-q, q) : nothing
    end
    hoop_stress_clims = symmetric_p99_clims(sigma_theta_theta_MPa)
    tangential_stress_clims = symmetric_p99_clims(tau_r_theta_MPa)

    temperature_plot = plot_q4_cloud(mesh_th.node, mesh_th.element, T_elem_final_K, coordinate_scale;
        title=@sprintf("Final Temperature Field, t = %.1f s", final_time_s),
        colorbar_title="T [K]",
        cmap=:thermal)
    hoop_stress_plot = plot_q4_cloud(
        mechanical_mesh.node, mechanical_mesh.bulk_element,
        sigma_theta_theta_MPa, coordinate_scale;
        title=@sprintf("Final Hoop Stress Field on Mechanical Submesh, t = %.1f s", final_time_s),
        colorbar_title="sigma_theta_theta [MPa]",
        cmap=:RdBu,
        clims=hoop_stress_clims)
    tangential_stress_plot = plot_q4_cloud(
        mechanical_mesh.node, mechanical_mesh.bulk_element,
        tau_r_theta_MPa, coordinate_scale;
        title=@sprintf("Final Tangential Shear Stress Field on Mechanical Submesh, t = %.1f s", final_time_s),
        colorbar_title="tau_r_theta [MPa]",
        cmap=:RdBu,
        clims=tangential_stress_clims)

    temperature_path = joinpath(output_dir, "final_temperature_field.png")
    hoop_stress_path = joinpath(output_dir, "final_hoop_stress_field.png")
    tangential_stress_path = joinpath(output_dir, "final_tangential_shear_stress_field.png")
    savefig(temperature_plot, temperature_path)
    savefig(hoop_stress_plot, hoop_stress_path)
    savefig(tangential_stress_plot, tangential_stress_path)

    @printf("  最终温度场已保存: %s\n", temperature_path)
    @printf("  最终环向应力场已保存: %s\n", hoop_stress_path)
    @printf("  最终切向剪应力场已保存: %s\n", tangential_stress_path)
    @printf("  最终环向应力范围: %.4e ~ %.4e MPa\n",
        minimum(sigma_theta_theta_MPa), maximum(sigma_theta_theta_MPa))
    @printf("  最终切向剪应力范围: %.4e ~ %.4e MPa\n",
        minimum(tau_r_theta_MPa), maximum(tau_r_theta_MPa))

    println("\n" * "="^80)
    println("全部完成")
    println("="^80)
end

main()
