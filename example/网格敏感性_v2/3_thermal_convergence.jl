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

    q0_nd = Q0 / scale.q

    variables = Dict{String,Any}()
    variables["thermal2D temperature at nodes"] = fill(param_dim.cell.T0 / scale.T_ref, nnode)
    variables["heat_source_fields"] = fill(q0_nd, ne)

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

    return (T_nodes=T_nodes, times=times, T_hist=T_hist,
            mesh=mesh, mesh_data=mesh_data, solve_wall=solve_wall)
end

function main()
    println("=" ^ 70)
    println("Script 3: 热学网格收敛性分析 (GCI, q0=$(Q0) W/m³)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")

    all_results = Dict{Int, Any}()
    h_vals = Float64[]
    wall_times = Float64[]

    for nθ in THERMAL_Nθ
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        r = run_thermal_case(param_dim, nθ)
        all_results[nθ] = r

        h = effective_h(r.mesh)
        push!(h_vals, h)
        push!(wall_times, r.solve_wall)

        T_final = r.T_nodes
        @printf("  耗时 %.2f s,  T_max = %.3f K,  T_min = %.3f K,  h = %.6f\n",
            r.solve_wall, maximum(T_final), minimum(T_final), h)
    end

    # ── 参考解 (nθ=160) ──
    ref_idx = length(THERMAL_Nθ)
    ref_result = all_results[THERMAL_Nθ[ref_idx]]
    T_hist_ref = ref_result.T_hist
    t_ref = ref_result.times
    nt_ref = length(t_ref)
    mesh_ref = ref_result.mesh

    # 参考节点笛卡尔坐标
    ref_x = mesh_ref.node[:, 1]
    ref_y = mesh_ref.node[:, 2]

    Tmax_ref_curve = [maximum(T_hist_ref[:, k]) for k in 1:nt_ref]
    Trange_ref_curve = [maximum(T_hist_ref[:, k]) - minimum(T_hist_ref[:, k]) for k in 1:nt_ref]

    # ── 收敛误差表 ──
    Tmax_l2 = Float64[]; Trange_l2 = Float64[]; Spatial_l2 = Float64[]
    Tmax_linf = Float64[]

    println("\n" * "=" ^ 70)
    println("热学网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %12s  %12s  %12s  %12s  %10s\n",
            "nθ", "h", "Tmax L2_rel%", "Trange L2", "Spatial L2%", "Tmax Linf", "Wall [s]")
    println("-" ^ 100)

    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        t_c = r.times
        nt_c = length(t_c)
        mesh_c = r.mesh

        if i == ref_idx
            push!(Tmax_l2, 0.0); push!(Trange_l2, 0.0)
            push!(Spatial_l2, 0.0); push!(Tmax_linf, 0.0)
            @printf("%-10d  %10.6f  %12s  %12s  %12s  %12s  %10.2f\n",
                    nθ, h_vals[i], "ref", "ref", "ref", "ref", wall_times[i])
            continue
        end

        # T_max(t) 曲线 L2_rel
        Tmax_cand = [maximum(T_hist[:, k]) for k in 1:nt_c]
        Tmax_cand_aligned = align_to_ref(t_c, Tmax_cand, t_ref)
        err_Tmax = l2_rel_norm(Tmax_cand_aligned, Tmax_ref_curve) * 100
        err_Tmax_inf = max_norm(Tmax_cand_aligned, Tmax_ref_curve)

        # T_range(t) L2（绝对值）
        Trange_cand = [maximum(T_hist[:, k]) - minimum(T_hist[:, k]) for k in 1:nt_c]
        Trange_cand_aligned = align_to_ref(t_c, Trange_cand, t_ref)
        err_Trange = l2_norm(Trange_cand_aligned, Trange_ref_curve)

        # 空间场 L2_rel（需要 IDW 插值到参考网格）
        err_spatial = NaN
        nt_overlap = min(nt_c, nt_ref)
        if size(T_hist, 1) != size(T_hist_ref, 1)
            # 节点数不同，使用 IDW 插值
            coarse_x = mesh_c.node[:, 1]
            coarse_y = mesh_c.node[:, 2]
            spatial_errs = Float64[]
            for k in 1:nt_overlap
                T_cand_interp = interpolate_to_ref_field(
                    T_hist[:, k], coarse_x, coarse_y, ref_x, ref_y)
                val = l2_rel_norm(T_cand_interp, T_hist_ref[:, k])
                isnan(val) || push!(spatial_errs, val)
            end
            err_spatial = isempty(spatial_errs) ? NaN : mean(spatial_errs) * 100
        else
            # 节点数相同，直接计算
            spatial_errs = Float64[]
            for k in 1:nt_overlap
                val = l2_rel_norm(T_hist[:, k], T_hist_ref[:, k])
                isnan(val) || push!(spatial_errs, val)
            end
            err_spatial = isempty(spatial_errs) ? NaN : mean(spatial_errs) * 100
        end

        push!(Tmax_l2, err_Tmax); push!(Trange_l2, err_Trange)
        push!(Spatial_l2, err_spatial); push!(Tmax_linf, err_Tmax_inf)

        @printf("%-10d  %10.6f  %12.4f  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, h_vals[i], err_Tmax, err_Trange, err_spatial,
                err_Tmax_inf, wall_times[i])
    end

    # ── GCI 计算 ──
    Tmax_end_vals = [maximum(all_results[nθ].T_nodes) for nθ in THERMAL_Nθ]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %12s  %8s  %12s\n",
            "Pair", "r", "epsilon", "p_obs", "GCI(T_max)")
    println("-" ^ 60)

    gci_results = []
    for i in 1:length(THERMAL_Nθ)-1
        # r > 1: 粗网格 h / 细网格 h（THERMAL_Nθ 从粗到细排列）
        r21 = h_vals[i] / h_vals[i+1]
        r32 = i+1 < length(THERMAL_Nθ) ? h_vals[i+1] / h_vals[i+2] : r21

        if i + 2 <= length(THERMAL_Nθ)
            p_T = observed_order(Tmax_end_vals[i], Tmax_end_vals[i+1],
                                 Tmax_end_vals[i+2], r21, r32)
        else
            p_T = 2.0
        end

        # compute_gci(f_fine, f_coarse, r): fine=细网格(i+1), coarse=粗网格(i)
        gci_T = compute_gci(Tmax_end_vals[i+1], Tmax_end_vals[i], r21;
                            p=isnan(p_T) ? 2.0 : p_T)
        push!(gci_results, (r21, p_T, gci_T))

        eps_T = abs(Tmax_end_vals[i+1] - Tmax_end_vals[i]) / abs(Tmax_end_vals[i])
        @printf("L%d→L%d       %8.3f  %12.6f  %8.2f  %12.4f\n",
                i, i+1, r21, eps_T, isnan(p_T) ? -1.0 : p_T, gci_T)
    end

    # ── 渐近收敛检查 ──
    if length(gci_results) >= 2
        println("\n渐近收敛检查:")
        for (name, gcis) in [("T_max", [r[3] for r in gci_results])]
            ratio = asymptotic_check(gcis[1], gcis[2],
                                     gci_results[1][1], gci_results[1][2])
            @printf("  %-8s: GCI_12/GCI_23 ratio = %.4f %s\n",
                    name, isnan(ratio) ? 0.0 : ratio,
                    isnan(ratio) ? "(insufficient data)" :
                    0.8 <= ratio <= 1.2 ? "✓ asymptotic" : "✗ not asymptotic")
        end
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    colors = [:red, :orange, :green, :blue]

    # 图1: 温度场时间演化
    p1 = plot(xlabel="Time [s]", ylabel="Mean Temperature [K]",
              title="Thermal Mesh: Temperature (q0=$(Q0) W/m³)",
              legend=:bottomright)
    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        if !isempty(T_hist)
            nt = size(T_hist, 2)
            T_mean = [mean(T_hist[:, k]) for k in 1:nt]
            label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
            ls = i == ref_idx ? :solid : :dash
            plot!(p1, r.times[1:nt], T_mean, label=label, lw=1.5, color=colors[i], ls=ls)
        end
    end
    savefig(p1, joinpath(out_dir, "thermal_temperature_convergence.png"))

    # 图2: log-log 收敛误差图（排除参考解的 0 误差）
    h_plot = h_vals[1:end-1]
    Tmax_plot = Tmax_l2[1:end-1]
    Spatial_plot = Spatial_l2[1:end-1]

    p3 = plot(xlabel="h (element size)", ylabel="L2_rel error [%]",
              title="Thermal Mesh: Convergence",
              xscale=:log10, yscale=:log10, legend=:topright)

    plot!(p3, h_plot, Tmax_plot, marker=:o, lw=2, label="T_max(t) L2_rel", color=:blue)
    plot!(p3, h_plot, Spatial_plot, marker=:d, lw=2, label="Spatial L2_rel", color=:green)

    # 理论斜率线 O(h^2) (Q4, k=1)
    if length(h_plot) >= 2 && Tmax_plot[1] > 0
        h_ref_line = [minimum(h_plot), maximum(h_plot)]
        e_ref = Tmax_plot[1]
        h_ref_val = h_plot[1]
        theory_line = e_ref .* (h_ref_line ./ h_ref_val).^2
        plot!(p3, h_ref_line, theory_line, lw=1, ls=:dash, color=:blue, label="O(h²) theory")
    end

    savefig(p3, joinpath(out_dir, "thermal_spatial_error.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
