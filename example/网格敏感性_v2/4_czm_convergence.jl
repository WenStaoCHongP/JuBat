"""
Script 4: CZM 网格收敛性分析 (GCI 框架)

先运行 Script 1 确定 nθ 区间，再用 4 组 nθ 运行全耦合模型。

主指标：断裂能耗散 E_frac(t) 的 GCI
辅指标：D_max(t) L2_rel, n_fractured(t) L2, δ_max_n(t) L2_rel, 牵引-分离面积偏差
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

"""从参数计算 CZM nθ 候选值"""
function get_czm_nθ(param_dim)
    E_NE  = param_dim.NE.E;  t_NE = param_dim.NE.thickness
    E_PE  = param_dim.PE.E;  t_PE = param_dim.PE.thickness
    E_eff = (E_NE * t_NE + E_PE * t_PE) / (t_NE + t_PE)
    # G_c   = param_dim.cohesive.G_c_n                                      # TODO Chunk 2 Task 2.1
    # σ_max = param_dim.cohesive.σ_max_n                                    # TODO Chunk 2 Task 2.1
    G_c   = NaN  # TODO Chunk 2 Task 2.1
    σ_max = NaN  # TODO Chunk 2 Task 2.1
    l_c   = G_c * E_eff / σ_max^2

    R_in  = param_dim.cell.Rin
    R_out = param_dim.cell.Rout

    nθ_inner = ceil(Int, 2π * R_in  / l_c)
    nθ_outer = ceil(Int, 2π * R_out / l_c)

    # 确保最小 nθ 足够大（全耦合模型需要足够细的网格才能收敛）
    # 参考 testexample.jl 使用 nθ=360，CZM 收敛分析使用 80-320 范围
    nθ_list = [80, 160, 240, 320]

    return nθ_list, l_c, E_eff
end

function run_czm_case(param_dim, nθ)
    opt = JuBat.Option()
    Crates = 3.0
    i = 5 * Crates
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 600]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
    case.czm_mesh = czm_mesh

    result = JuBat.Solve(case)
    return result, czm_mesh, mesh_data
end

"""计算断裂能耗散 E_frac(t)"""
function compute_fracture_energy(result, czm_mesh, param_dim)
    # G_c = param_dim.cohesive.G_c_n                                        # TODO Chunk 2 Task 2.1
    G_c = NaN  # TODO Chunk 2 Task 2.1
    scale = param_dim.scale
    nt = length(result["time [s]"])

    coh_lengths = [elem.length * scale.L for elem in czm_mesh.cohesive_elements]

    E_frac = zeros(Float64, nt)
    if haskey(result, "czm damage [0-1]")
        D_all = result["czm damage [0-1]"]
        if ndims(D_all) == 2
            n_coh = min(length(coh_lengths), size(D_all, 1))
            for k in 1:min(nt, size(D_all, 2))
                for e in 1:n_coh
                    E_frac[k] += G_c * coh_lengths[e] * D_all[e, k]
                end
            end
        end
    end
    return E_frac
end

function main()
    println("=" ^ 70)
    println("Script 4: CZM 网格收敛性分析 (GCI)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    nθ_list, l_c, E_eff = get_czm_nθ(param_dim)
    @printf("\nl_c = %.1f μm,  E_eff = %.2f GPa\n", l_c * 1e6, E_eff * 1e-9)
    @printf("CZM nθ 候选值: %s\n", string(nθ_list))

    all_results  = Dict{Int, Dict}()
    all_czm      = Dict{Int, Any}()
    all_wall     = Float64[]
    h_vals       = Float64[]
    E_frac_curves = Dict{Int, Vector{Float64}}()

    for nθ in nθ_list
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        t0 = time_ns()
        r, czm, mesh_data = run_czm_case(param_dim, nθ)
        dt_wall = (time_ns() - t0) * 1e-9

        all_results[nθ] = r
        all_czm[nθ] = czm

        h = effective_h(mesh_data.thermal2D)
        push!(h_vals, h)
        push!(all_wall, dt_wall)

        E_frac = compute_fracture_energy(r, czm, param_dim)
        E_frac_curves[nθ] = E_frac

        Dmax_end = haskey(r, "czm D_max") ? r["czm D_max"][end] : NaN
        nfrac_end = haskey(r, "czm n_fractured") ? r["czm n_fractured"][end] : NaN

        nfrac_display = isnan(nfrac_end) ? 0 : Int(nfrac_end)
        @printf("  耗时 %.2f s,  D_max = %.4f,  n_frac = %d,  E_frac = %.4e,  h = %.6f\n",
                dt_wall, Dmax_end, nfrac_display, E_frac[end], h)
    end

    # ── 参考解 ──
    ref_idx = length(nθ_list)
    ref_result = all_results[nθ_list[ref_idx]]
    t_ref = ref_result["time [s]"]

    # ── 收敛误差表 ──
    D_l2 = Float64[]; nfrac_l2 = Float64[]; delta_l2 = Float64[]
    area_errors = Float64[]; Efrac_l2 = Float64[]

    println("\n" * "=" ^ 70)
    println("CZM 网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n",
            "nθ", "n_coh", "h",
            "E_frac L2", "D_max L2%", "n_frac L2", "δ L2%", "Area err%", "Wall [s]")
    println("-" ^ 115)

    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        t_c = r["time [s]"]
        n_coh = all_czm[nθ].n_cohesive

        if i == ref_idx
            push!(D_l2, 0.0); push!(nfrac_l2, 0.0); push!(delta_l2, 0.0)
            push!(area_errors, 0.0); push!(Efrac_l2, 0.0)
            @printf("%-10d  %10d  %10.6f  %12s  %12s  %12s  %12s  %12s  %10.2f\n",
                    nθ, n_coh, h_vals[i], "ref", "ref", "ref", "ref", "ref", all_wall[i])
            continue
        end

        # E_frac(t) L2（主指标）
        E_aligned = align_to_ref(t_c, E_frac_curves[nθ], t_ref)
        err_Efrac = l2_norm(E_aligned, E_frac_curves[nθ_list[ref_idx]])
        push!(Efrac_l2, err_Efrac)

        # D_max(t) L2_rel
        err_D = NaN
        if haskey(r, "czm D_max") && haskey(ref_result, "czm D_max")
            D_aligned = align_to_ref(t_c, r["czm D_max"], t_ref)
            err_D = l2_rel_norm(D_aligned, ref_result["czm D_max"]) * 100
        end
        push!(D_l2, err_D)

        # n_fractured(t) L2
        err_nf = NaN
        if haskey(r, "czm n_fractured") && haskey(ref_result, "czm n_fractured")
            nf_aligned = align_to_ref(t_c, r["czm n_fractured"], t_ref)
            err_nf = l2_norm(nf_aligned, ref_result["czm n_fractured"])
        end
        push!(nfrac_l2, err_nf)

        # δ_max_n(t) L2_rel
        err_d = NaN
        if haskey(r, "czm δ_max_n [m]") && haskey(ref_result, "czm δ_max_n [m]")
            d_aligned = align_to_ref(t_c, r["czm δ_max_n [m]"], t_ref)
            err_d = l2_rel_norm(d_aligned, ref_result["czm δ_max_n [m]"]) * 100
        end
        push!(delta_l2, err_d)

        # 牵引-分离面积偏差
        err_area = NaN
        if haskey(r, "czm traction normal [Pa]") && haskey(ref_result, "czm traction normal [Pa]")
            traction_c = r["czm traction normal [Pa]"]
            sep_c = r["czm separation normal [m]"]
            traction_ref = ref_result["czm traction normal [Pa]"]
            sep_ref = ref_result["czm separation normal [m]"]

            if ndims(traction_c) == 2 && ndims(traction_ref) == 2
                D_c = haskey(r, "czm damage [0-1]") ? r["czm damage [0-1]"] : nothing
                D_r = haskey(ref_result, "czm damage [0-1]") ? ref_result["czm damage [0-1]"] : nothing

                peak_c = D_c !== nothing && ndims(D_c) == 2 ? argmax(D_c[:, end]) : 1
                peak_r = D_r !== nothing && ndims(D_r) == 2 ? argmax(D_r[:, end]) : 1

                err_area = area_error(
                    sep_c[peak_c, :], traction_c[peak_c, :],
                    sep_ref[peak_r, :], traction_ref[peak_r, :])
            end
        end
        push!(area_errors, err_area)

        @printf("%-10d  %10d  %10.6f  %12.4e  %12.4f  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, n_coh, h_vals[i], err_Efrac, err_D, err_nf, err_d, err_area, all_wall[i])
    end

    # ── GCI 计算 ──
    Efrac_end_vals = [E_frac_curves[nθ][end] for nθ in nθ_list]
    Dmax_end_vals = [haskey(all_results[nθ], "czm D_max") ? all_results[nθ]["czm D_max"][end] : NaN for nθ in nθ_list]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %8s  %12s  %12s\n",
            "Pair", "r", "p_obs", "GCI(E_frac)", "GCI(D_max)")
    println("-" ^ 65)

    gci_results = []
    for i in 1:length(nθ_list)-1
        # r > 1: 粗网格 h / 细网格 h（nθ_list 从粗到细排列）
        r21 = h_vals[i] / h_vals[i+1]
        r32 = i+1 < length(nθ_list) ? h_vals[i+1] / h_vals[i+2] : r21

        if i + 2 <= length(nθ_list)
            p_E = observed_order(Efrac_end_vals[i], Efrac_end_vals[i+1],
                                 Efrac_end_vals[i+2], r21, r32)
            p_D = observed_order(Dmax_end_vals[i], Dmax_end_vals[i+1],
                                 Dmax_end_vals[i+2], r21, r32)
        else
            p_E = NaN; p_D = NaN
        end

        # compute_gci(f_fine, f_coarse, r): fine=细网格(i+1), coarse=粗网格(i)
        gci_E = compute_gci(Efrac_end_vals[i+1], Efrac_end_vals[i], r21;
                            p=isnan(p_E) ? 2.0 : p_E)
        gci_D = isnan(Dmax_end_vals[i]) ? NaN :
                compute_gci(Dmax_end_vals[i+1], Dmax_end_vals[i], r21;
                            p=isnan(p_D) ? 2.0 : p_D)

        push!(gci_results, (r21, p_E, gci_E, gci_D))

        @printf("L%d→L%d       %8.3f  %8.2f  %12.4f  %12.4f\n",
                i, i+1, r21,
                isnan(p_E) ? -1.0 : p_E,
                gci_E, isnan(gci_D) ? NaN : gci_D)
    end

    # ── 渐近收敛检查 ──
    if length(gci_results) >= 2
        println("\n渐近收敛检查:")
        for (name, gcis) in [("E_frac", [r[3] for r in gci_results]),
                             ("D_max", [r[4] for r in gci_results])]
            gcis_valid = filter(!isnan, gcis)
            if length(gcis_valid) >= 2
                ratio = asymptotic_check(gcis_valid[1], gcis_valid[2],
                                         gci_results[1][1], gci_results[1][2])
                @printf("  %-8s: GCI_12/GCI_23 ratio = %.4f %s\n",
                        name, isnan(ratio) ? 0.0 : ratio,
                        isnan(ratio) ? "(insufficient data)" :
                        0.8 <= ratio <= 1.2 ? "✓ asymptotic" : "✗ not asymptotic")
            end
        end
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    colors = [:red, :orange, :green, :blue]

    # 图1: D_max 演化对比
    p1 = plot(xlabel="Time [s]", ylabel="D_max",
              title="CZM Mesh: Damage Evolution", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        haskey(r, "czm D_max") || continue
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p1, r["time [s]"], r["czm D_max"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p1, joinpath(out_dir, "czm_damage_evolution.png"))

    # 图2: 断裂能耗散对比（新增主指标图）
    p2 = plot(xlabel="Time [s]", ylabel="E_frac [J/m]",
              title="CZM Mesh: Fracture Energy Dissipation", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p2, r["time [s]"], E_frac_curves[nθ],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p2, joinpath(out_dir, "czm_fracture_energy.png"))

    # 图3: log-log 收敛误差图
    if ref_idx > 1
        h_plot = h_vals[1:end-1]
        Efrac_plot = Efrac_l2[1:end-1]
        D_plot = D_l2[1:end-1]

        p4 = plot(xlabel="h (element size)", ylabel="L2 error",
                  title="CZM Mesh: Convergence",
                  xscale=:log10, yscale=:log10, legend=:topright)
        plot!(p4, h_plot, Efrac_plot, marker=:o, lw=2,
              label="E_frac(t) L2", color=:blue)
        plot!(p4, h_plot, D_plot, marker=:s, lw=2,
              label="D_max(t) L2_rel%", color=:orange)

        # 理论斜率线
        if length(h_plot) >= 2 && Efrac_plot[1] > 0
            h_ref_line = [minimum(h_plot), maximum(h_plot)]
            e_ref = Efrac_plot[1]
            h_ref_val = h_plot[1]
            theory_line = e_ref .* (h_ref_line ./ h_ref_val).^2
            plot!(p4, h_ref_line, theory_line, lw=1, ls=:dash,
                  color=:blue, label="O(h²) theory")
        end

        savefig(p4, joinpath(out_dir, "czm_error_convergence.png"))
    end

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()