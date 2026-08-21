"""
Script 4: CZM 网格敏感性分析

先运行 Script 1 确定 nθ 区间，再用 4 组 nθ 运行全耦合模型
（热-化学载荷驱动 CZM 损伤），对比 D_max(t)、n_fractured(t)、
δ_max_n(t) 的 RMSPE 和牵引-分离面积偏差。

指标：D_max(t) RMSPE、n_frac(t) RMSPE、δ_max_n(t) RMSPE、Traction-Sep 面积偏差
输出：收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_rmspe_utils.jl"))

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

    # 4 等分区间，取整
    step = max(1, round(Int, (nθ_outer - nθ_inner) / 3))
    nθ_list = [nθ_inner + i * step for i in 0:3]
    nθ_list[4] = nθ_outer

    return nθ_list, l_c, E_eff
end

function run_czm_case(param_dim, nθ; E_eff::Float64)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

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
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    case.czm_mesh = czm_mesh

    result = JuBat.Solve(case)
    return result, czm_mesh
end

function main()
    println("=" ^ 70)
    println("Script 4: CZM 网格敏感性分析")
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
    all_Dmax_end = Float64[]
    all_Dmean_end = Float64[]
    all_nfrac_end = Float64[]

    for nθ in nθ_list
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        t0 = time_ns()
        r, czm = run_czm_case(param_dim, nθ; E_eff=E_eff)
        dt_wall = (time_ns() - t0) * 1e-9

        all_results[nθ] = r
        all_czm[nθ] = czm

        Dmax_end  = haskey(r, "czm D_max")  ? r["czm D_max"][end]  : NaN
        Dmean_end = haskey(r, "czm D_mean") ? r["czm D_mean"][end] : NaN
        nfrac_end = haskey(r, "czm n_fractured") ? r["czm n_fractured"][end] : NaN

        push!(all_wall, dt_wall)
        push!(all_Dmax_end, Dmax_end)
        push!(all_Dmean_end, Dmean_end)
        push!(all_nfrac_end, nfrac_end)

        @printf("  耗时 %.2f s,  D_max = %.4f,  D_mean = %.6f,  n_frac = %d\n",
                dt_wall, Dmax_end, Dmean_end, Int(nfrac_end))
    end

    # ── 参考解 (最细) ──
    ref_idx = length(nθ_list)
    ref_result = all_results[nθ_list[ref_idx]]
    t_ref = ref_result["time [s]"]

    println("\n" * "=" ^ 70)
    println("CZM 网格收敛性汇总 (RMSPE)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %12s  %12s  %12s  %12s  %10s  %10s  %10s\n",
            "nθ", "n_coh",
            "D_max RMSPE%", "n_frac RMSPE%", "δ_max_n RMSPE%", "Area err%",
            "D_skip%", "nf_skip%", "Wall [s]")
    println("-" ^ 105)

    D_rmspe = Float64[]
    nfrac_rmspe = Float64[]
    delta_rmspe = Float64[]
    area_errors = Float64[]

    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        t_c = r["time [s]"]
        n_coh = all_czm[nθ].n_cohesive

        if i == ref_idx
            push!(D_rmspe, 0.0)
            push!(nfrac_rmspe, 0.0)
            push!(delta_rmspe, 0.0)
            push!(area_errors, 0.0)
            @printf("%-10d  %10d  %12s  %12s  %12s  %12s  %10s  %10s  %10.2f\n",
                    nθ, n_coh, "ref", "ref", "ref", "ref", "-", "-", all_wall[i])
            continue
        end

        # D_max(t) RMSPE
        skip_D = NaN
        if haskey(r, "czm D_max") && haskey(ref_result, "czm D_max")
            D_aligned = align_to_ref(t_c, r["czm D_max"], t_ref)
            err_D, skip_D = rmspe(D_aligned, ref_result["czm D_max"])
        else
            err_D = NaN
        end

        # n_fractured(t) RMSPE
        skip_nf = NaN
        if haskey(r, "czm n_fractured") && haskey(ref_result, "czm n_fractured")
            nf_aligned = align_to_ref(t_c, r["czm n_fractured"], t_ref)
            err_nf, skip_nf = rmspe(nf_aligned, ref_result["czm n_fractured"])
        else
            err_nf = NaN
        end

        # δ_max_n(t) RMSPE
        if haskey(r, "czm δ_max_n [m]") && haskey(ref_result, "czm δ_max_n [m]")
            d_aligned = align_to_ref(t_c, r["czm δ_max_n [m]"], t_ref)
            err_d, _ = rmspe(d_aligned, ref_result["czm δ_max_n [m]"])
        else
            err_d = NaN
        end

        # 牵引-分离面积偏差（峰值损伤单元）
        err_area = NaN
        if haskey(r, "czm traction normal [Pa]") && haskey(ref_result, "czm traction normal [Pa]")
            traction_c = r["czm traction normal [Pa]"]
            sep_c = r["czm separation normal [m]"]
            traction_ref = ref_result["czm traction normal [Pa]"]
            sep_ref = ref_result["czm separation normal [m]"]

            if ndims(traction_c) == 2 && ndims(traction_ref) == 2
                # 选取最终时刻 D 值最大的单元
                D_c = haskey(r, "czm damage [0-1]") ? r["czm damage [0-1]"] : nothing
                D_r = haskey(ref_result, "czm damage [0-1]") ? ref_result["czm damage [0-1]"] : nothing

                peak_c = D_c !== nothing && ndims(D_c) == 2 ? argmax(D_c[:, end]) : 1
                peak_r = D_r !== nothing && ndims(D_r) == 2 ? argmax(D_r[:, end]) : 1

                err_area = area_error(
                    sep_c[peak_c, :], traction_c[peak_c, :],
                    sep_ref[peak_r, :], traction_ref[peak_r, :])
            end
        end

        push!(D_rmspe, err_D)
        push!(nfrac_rmspe, err_nf)
        push!(delta_rmspe, err_d)
        push!(area_errors, err_area)

        @printf("%-10d  %10d  %12.4f  %12.4f  %12.4f  %12.4f  %10.1f  %10.1f  %10.2f\n",
                nθ, n_coh, err_D, err_nf, err_d, err_area,
                skip_D*100, skip_nf*100, all_wall[i])
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "4_czm_mesh_sensitivity")
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

    # 图2: n_fractured 演化
    p2 = plot(xlabel="Time [s]", ylabel="n_fractured",
              title="CZM Mesh: Fractured Elements", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        haskey(r, "czm n_fractured") || continue
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p2, r["time [s]"], r["czm n_fractured"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p2, joinpath(out_dir, "czm_fractured_evolution.png"))

    # 图3: δ_max_n 演化
    p3 = plot(xlabel="Time [s]", ylabel="δ_max_n [nm]",
              title="CZM Mesh: Max Normal Separation", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        haskey(r, "czm δ_max_n [m]") || continue
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p3, r["time [s]"], r["czm δ_max_n [m]"] .* 1e9,
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p3, joinpath(out_dir, "czm_separation_evolution.png"))

    # 图4: RMSPE 收敛误差
    if ref_idx > 1
        p4 = plot(xlabel="nθ (per revolution)", ylabel="RMSPE [%]",
                  title="CZM Mesh: Convergence Error (RMSPE)",
                  legend=:topright, yscale=:log10)
        plot!(p4, nθ_list[1:end-1], D_rmspe[1:end-1],
              marker=:o, lw=2, label="D_max(t)", color=:blue)
        plot!(p4, nθ_list[1:end-1], nfrac_rmspe[1:end-1],
              marker=:s, lw=2, label="n_frac(t)", color=:orange)
        plot!(p4, nθ_list[1:end-1], delta_rmspe[1:end-1],
              marker=:d, lw=2, label="δ_max_n(t)", color=:green)
        plot!(p4, nθ_list[1:end-1], area_errors[1:end-1],
              marker=:p, lw=2, label="Traction-Sep area", color=:purple)
        hline!(p4, [5.0], label="5%", color=:red, ls=:dash)
        savefig(p4, joinpath(out_dir, "czm_convergence_error.png"))
    end

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()
