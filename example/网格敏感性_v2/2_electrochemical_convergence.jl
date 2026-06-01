"""
Script 2: 电化学网格收敛性分析 (GCI 框架)

4 组 (Nn, Ns, Np) 下运行 SPMe + lumped thermal：
  - (40, 20, 40)  ← 参考解
  - (20, 10, 20)
  - (20,  5, 20)
  - (10,  5, 10)

指标：V(t) L2_rel, T(t) L2_rel, dT/dt L2
GCI 标量：V_end, T_peak, max|dT/dt|
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const MESH_CONFIGS = [
    (40, 20, 40),  # Level 1 (参考解)
    (20, 10, 20),  # Level 2
    (20,  5, 20),  # Level 3
    (10,  5, 10),  # Level 4 (最粗)
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
    opt.czm_enabled = false

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    return result
end

"""电化学 Track 的近似 h（无热网格时）"""
function echem_h(Nn, Ns, Np)
    return 1.0 / sqrt(Nn + Ns + Np)
end

function main()
    println("=" ^ 70)
    println("Script 2: 电化学网格收敛性分析 (GCI)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    results = Dict{Tuple{Int,Int,Int}, Dict}()
    wall_times = Float64[]
    h_vals = Float64[]

    for (Nn, Ns, Np) in MESH_CONFIGS
        label = "($Nn, $Ns, $Np)"
        @printf("\n--- 运行 %s ---\n", label)
        t0 = time_ns()
        r = run_echem_case(param_dim, Nn, Ns, Np)
        dt_wall = (time_ns() - t0) * 1e-9
        push!(wall_times, dt_wall)
        results[(Nn, Ns, Np)] = r
        push!(h_vals, echem_h(Nn, Ns, Np))
        @printf("  耗时 %.2f s,  %d 步\n", dt_wall, length(r["time [s]"]))
    end

    # ── 参考解 ──
    ref = results[MESH_CONFIGS[1]]
    t_ref = ref["time [s]"]
    V_ref = ref["cell voltage [V]"]
    T_ref = ref["temperature [K]"]

    # ── 收敛误差表 ──
    V_l2 = Float64[]; T_l2 = Float64[]; dTdt_l2 = Float64[]
    V_linf = Float64[]; T_linf = Float64[]

    println("\n" * "=" ^ 70)
    println("电化学网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-16s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n",
            "Mesh", "h",
            "V L2_rel%", "T L2_rel%", "dT/dt L2",
            "V Linf", "T Linf", "Wall [s]")
    println("-" ^ 110)

    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        t_c = r["time [s]"]
        V_c = r["cell voltage [V]"]
        T_c = r["temperature [K]"]

        V_aligned = align_to_ref(t_c, V_c, t_ref)
        T_aligned = align_to_ref(t_c, T_c, t_ref)

        if i == 1
            push!(V_l2, 0.0); push!(T_l2, 0.0); push!(dTdt_l2, 0.0)
            push!(V_linf, 0.0); push!(T_linf, 0.0)
            @printf("(%2d,%2d,%2d)       %10.4f  %12s  %12s  %12s  %12s  %12s  %10.2f\n",
                    Nn, Ns, Np, h_vals[i], "ref", "ref", "ref", "ref", "ref", wall_times[i])
        else
            err_V_l2 = l2_rel_norm(V_aligned, V_ref) * 100
            err_T_l2 = l2_rel_norm(T_aligned, T_ref) * 100
            err_V_inf = max_norm(V_aligned, V_ref)
            err_T_inf = max_norm(T_aligned, T_ref)

            dTdt_ref_vals = diff(T_ref) ./ diff(t_ref)
            dTdt_cand_vals = diff(T_aligned) ./ diff(t_ref)
            err_dTdt = l2_norm(dTdt_cand_vals, dTdt_ref_vals)

            push!(V_l2, err_V_l2); push!(T_l2, err_T_l2); push!(dTdt_l2, err_dTdt)
            push!(V_linf, err_V_inf); push!(T_linf, err_T_inf)

            @printf("(%2d,%2d,%2d)       %10.4f  %12.4f  %12.4f  %12.4f  %12.6f  %12.6f  %10.2f\n",
                    Nn, Ns, Np, h_vals[i], err_V_l2, err_T_l2, err_dTdt,
                    err_V_inf, err_T_inf, wall_times[i])
        end
    end

    # ── GCI 计算 ──
    # GCI 标量：V_end, T_peak
    V_end_vals = [results[c]["cell voltage [V]"][end] for c in MESH_CONFIGS]
    T_peak_vals = [maximum(results[c]["temperature [K]"]) for c in MESH_CONFIGS]
    dTdt_max_vals = [let
        r = results[c]; t = r["time [s]"]; T = r["temperature [K]"]
        maximum(abs.(diff(T) ./ diff(t)))
    end for c in MESH_CONFIGS]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %12s  %8s  %12s  %12s  %12s\n",
            "Pair", "r", "epsilon", "p_obs", "GCI(V_end)", "GCI(T_peak)", "GCI(dTdt)")
    println("-" ^ 90)

    gci_results = []
    p_V = NaN; p_T = NaN; p_dTdt = NaN  # 在循环外初始化，避免作用域问题
    for i in 1:length(MESH_CONFIGS)-1
        r21 = h_vals[i+1] / h_vals[i]
        r32 = i+1 < length(MESH_CONFIGS) ? h_vals[i+2] / h_vals[i+1] : r21

        # 三组值计算 p（如果可用）
        if i + 2 <= length(MESH_CONFIGS)
            p_V = observed_order(V_end_vals[i], V_end_vals[i+1], V_end_vals[i+2], r21, r32)
            p_T = observed_order(T_peak_vals[i], T_peak_vals[i+1], T_peak_vals[i+2], r21, r32)
            p_dTdt = observed_order(dTdt_max_vals[i], dTdt_max_vals[i+1], dTdt_max_vals[i+2], r21, r32)
        else
            p_V = 2.0; p_T = 2.0; p_dTdt = 2.0  # 回退到理论值
        end

        gci_V = compute_gci(V_end_vals[i], V_end_vals[i+1], r21; p=isnan(p_V) ? 2.0 : p_V)
        gci_T = compute_gci(T_peak_vals[i], T_peak_vals[i+1], r21; p=isnan(p_T) ? 2.0 : p_T)
        gci_dT = compute_gci(dTdt_max_vals[i], dTdt_max_vals[i+1], r21; p=isnan(p_dTdt) ? 2.0 : p_dTdt)

        push!(gci_results, (r21, p_V, gci_V, gci_T, gci_dT))

        eps_V = abs(V_end_vals[i+1] - V_end_vals[i]) / abs(V_end_vals[i])
        @printf("L%d→L%d       %8.3f  %12.6f  %8.2f  %12.4f  %12.4f  %12.4f\n",
                i, i+1, r21, eps_V,
                isnan(p_V) ? -1.0 : p_V,
                gci_V, gci_T, gci_dT)
    end

    # ── 渐近收敛检查 ──
    if length(gci_results) >= 2
        println("\n渐近收敛检查:")
        for (name, gcis) in [("V_end", [r[3] for r in gci_results]),
                             ("T_peak", [r[4] for r in gci_results]),
                             ("dTdt", [r[5] for r in gci_results])]
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

    # 图3: log-log 收敛误差图（含理论斜率线）
    h_plot = h_vals[2:end]
    V_plot = V_l2[2:end]; T_plot = T_l2[2:end]

    p3 = plot(xlabel="h (element size)", ylabel="L2_rel error [%]",
              title="Electrochemical Mesh: Convergence",
              xscale=:log10, yscale=:log10, legend=:topright)

    plot!(p3, h_plot, V_plot, marker=:o, lw=2, label="V(t) L2_rel", color=:blue)
    plot!(p3, h_plot, T_plot, marker=:s, lw=2, label="T(t) L2_rel", color=:orange)

    # 理论斜率线 O(h^2)
    if length(h_plot) >= 2
        h_ref_line = [minimum(h_plot), maximum(h_plot)]
        # 基于最细网格点的误差值，画理论斜率
        e_ref_V = V_plot[1]
        e_ref_T = T_plot[1]
        h_ref_val = h_plot[1]
        theory_V = e_ref_V .* (h_ref_line ./ h_ref_val).^2
        theory_T = e_ref_T .* (h_ref_line ./ h_ref_val).^2
        plot!(p3, h_ref_line, theory_V, lw=1, ls=:dash, color=:blue, label="O(h²) V ref")
        plot!(p3, h_ref_line, theory_T, lw=1, ls=:dash, color=:orange, label="O(h²) T ref")
    end

    # 标注观测收敛阶
    if !isnan(p_V) || !isnan(p_T)
        annotate_text = ""
        !isnan(p_V) && (annotate_text *= "p_obs(V)=$(round(p_V, digits=2))\n")
        !isnan(p_T) && (annotate_text *= "p_obs(T)=$(round(p_T, digits=2))\n")
        annotate_text *= "p_theory=2.0"
        annotate!(p3, maximum(h_plot), minimum([V_plot; T_plot]) * 1.5,
                  text(annotate_text, 8, :left, :bottom, :Courier))
    end

    savefig(p3, joinpath(out_dir, "echem_error_convergence.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()