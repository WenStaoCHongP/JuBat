
# verify_czm_standalone.jl
#
# CZM 求解器收敛性对比验证
# - 纯 Mode I 法向加载
# - 使用 Jellyroll.jl 生产级参数
# - 对比 basic / load_substep / arc_length 三种求解器
# - 生产级网格 (nθ=80)
#
# 运行方式:
#   julia --project=. tools/verify_czm_standalone.jl

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "czm_standalone")

# ========================================================================
# 结果记录
# ========================================================================

struct MethodResult
    converged::Bool
    iterations::Int
    D_max::Float64
    D_mean::Float64
    residual_norm::Float64
    n_fractured::Int
end

# ========================================================================
# 主函数
# ========================================================================

function main()
    mkpath(OUTPUT_DIR)

    println("=" ^ 70)
    println("CZM Solver Convergence Comparison (Mode I)")
    println("=" ^ 70)

    # ── 1. 参数与网格 ──────────────────────────────────────────────
    param_dim = JuBat.ChooseCell("Jellyroll")

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.gsorder = 2
    opt.thermal_enabled = false
    opt.thermalmodel = "none"
    opt.per_element_spme = false
    opt.czm_enabled = true
    opt.czm_model = "model1"
    opt.czm_iter_method = "load_substep"
    opt.czm_max_iter = 200
    opt.czm_tol = 1e-4
    opt.czm_load_steps = 50
    opt.czm_viscous_enabled = false

    case = JuBat.SetCase(param_dim, opt)

    nθ = 40
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, nθ_czm=40, gsorder=2)
    czm_mesh_template = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_data.thermal2D, case.param)

    println("  Mesh: nθ=$nθ, gsorder=2")
    println("  Nodes: $(czm_mesh_template.nnode)")
    println("  Bulk elements: $(size(czm_mesh_template.bulk_element, 1))")
    println("  Cohesive elements: $(czm_mesh_template.n_cohesive)")

    # ── 2. 有效参数 ────────────────────────────────────────────────
    # TODO Chunk 4: compute_czm_effective_params 已被 compute_czm_params_per_interface 替换
    # czm_param_cache = JuBat.compute_czm_params_per_interface(case)
    # pe = czm_param_cache.by_interface[:PE_PCC]
    # E_eff, ν_eff, α_eff = pe.E_eff, pe.ν, pe.α
    # β_n = case.param.NE.Omega / 3.0
    # β_p = case.param.PE.Omega / 3.0
    E_eff = ν_eff = α_eff = β_n = β_p = NaN  # placeholder
    @printf("  E_eff = %.6e, ν_eff = %.4f\n", E_eff, ν_eff)
    @printf("  α_eff = %.6e, β_n = %.6e, β_p = %.6e\n", α_eff, β_n, β_p)

    param = case.param
    czm_params = case.param.cohesive  # 归一化后的参数（不是 param_dim.cohesive）
    ndof = 2 * czm_mesh_template.nnode

    # @printf("  Cohesive: σ_max_n=%.2e, K_n=%.2e, δ_0_n=%.2e, G_c_n=%.2e, δ_c_n=%.2e\n",  # TODO Chunk 2 Task 2.1
    #     czm_params.σ_max_n, czm_params.K_n, czm_params.δ_0_n,                            # TODO Chunk 2 Task 2.1
    #     czm_params.G_c_n, czm_params.δ_c_n)                                              # TODO Chunk 2 Task 2.1

    # ── 3. 构建缓存 ───────────────────────────────────────────────
    cache = JuBat.build_czm_cache(czm_mesh_template, E_eff, ν_eff, param)

    # ── 4. 载荷水平 ────────────────────────────────────────────────
    # 使用 Δsoc_n 驱动（生产中损伤主要由 SOC 变化引起）
    # dT_elem 设为 0，Δsoc_p 设为 0
    # Δsoc_n 的量级需要覆盖弹性 → 损伤起始 → 断裂
    # baseline_probe 使用 Δsoc=0.0002（弹性区），生产中 Δsoc 可达 0.3-0.5
    ne = size(czm_mesh_template.bulk_element, 1)
    soc_n_levels = [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]

    methods = ["basic", "load_substep", "arc_length"]

    println("\n--- Solver parameters ---")
    println("  tol = 1e-4, max_iter = 200, n_load_steps = 50")
    println("  visc_beta = 1.0 (no viscous regularization)")
    println("  F_ext = zeros (no external mechanical force)")
    println("  dT_elem = 0, Δsoc_p_elem = 0 (SOC driving only)")

    # ── 5. 运行对比 ───────────────────────────────────────────────
    println("\n" * "=" ^ 70)
    println("Convergence Comparison Table")
    println("=" ^ 70)

    header = @sprintf("Δsoc_n | %-22s | %-22s | %-22s", "basic", "load_substep", "arc_length")
    println(header)
    println("-" ^ length(header))

    results_all = Dict{String, Vector{MethodResult}}()
    for m in methods
        results_all[m] = MethodResult[]
    end

    for soc_n_val in soc_n_levels
        dT_elem = fill(0.0, ne)
        Δsoc_n_elem = fill(soc_n_val, ne)
        Δsoc_p_elem = fill(0.0, ne)
        F_ext = zeros(Float64, ndof)

        row_parts = String[]

        for method in methods
            # 每个方法使用独立的 czm_mesh（fresh damage_states）
            czm_mesh_fresh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_data.thermal2D, param_dim)
            u_prev = zeros(Float64, ndof)

            result, updated_mesh = JuBat.solve_czm_step(
                czm_mesh_fresh, F_ext, E_eff, ν_eff, czm_params, param, u_prev;
                α_eff=α_eff, β_n=β_n, β_p=β_p,
                dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
                max_iter=200, tol=1e-4, n_load_steps=50,
                iter_method=method, cache=cache
            )

            stats = JuBat.get_damage_statistics(updated_mesh)

            mr = MethodResult(
                result.converged,
                result.iterations,
                stats.max_D,
                stats.mean_D,
                result.residual_norm,
                stats.n_fractured
            )
            push!(results_all[method], mr)

            if mr.converged
                part = @sprintf("OK %3dit D=%.4f r=%.1e", mr.iterations, mr.D_max, mr.residual_norm)
            else
                part = @sprintf("FAIL %3dit D=%.4f r=%.1e", mr.iterations, mr.D_max, mr.residual_norm)
            end
            push!(row_parts, part)
        end

        @printf("%6.3f | %s | %s | %s\n",
            soc_n_val, row_parts[1], row_parts[2], row_parts[3])
    end

    # ── 6. 汇总 ───────────────────────────────────────────────────
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)

    for method in methods
        results = results_all[method]
        n_converged = count(r -> r.converged, results)
        max_D_final = isempty(results) ? 0.0 : results[end].D_max
        max_D_overall = isempty(results) ? 0.0 : maximum(r -> r.D_max, results)
        total_iter = sum(r -> r.iterations, results)

        @printf("  %-15s: converged %d/%d levels, D_max=%.4f (final) / %.4f (peak), total_iter=%d\n",
            method, n_converged, length(soc_n_levels), max_D_final, max_D_overall, total_iter)
    end

    # ── 7. 写入报告 ───────────────────────────────────────────────
    report_path = joinpath(OUTPUT_DIR, "solver_comparison.txt")
    open(report_path, "w") do io
        println(io, "CZM Solver Convergence Comparison (Mode I)")
        println(io, "=" ^ 60)
        println(io, "")
        println(io, "Mesh: nθ=$nθ, gsorder=2")
        println(io, "Cohesive elements: $(czm_mesh_template.n_cohesive)")
        @printf(io, "E_eff=%.6e, ν_eff=%.4f, α_eff=%.6e\n", E_eff, ν_eff, α_eff)
        # @printf(io, "σ_max_n=%.2e, K_n=%.2e, δ_0_n=%.2e, G_c_n=%.2e, δ_c_n=%.2e\n",  # TODO Chunk 2 Task 2.1
        #     czm_params.σ_max_n, czm_params.K_n, czm_params.δ_0_n,                      # TODO Chunk 2 Task 2.1
        #     czm_params.G_c_n, czm_params.δ_c_n)                                        # TODO Chunk 2 Task 2.1
        println(io, "tol=1e-4, max_iter=200, n_load_steps=50, visc_beta=1.0")
        println(io, "")
        println(io, @sprintf("%-8s | %-22s | %-22s | %-22s", "Δsoc_n", "basic", "load_substep", "arc_length"))
        println(io, "-" ^ 80)

        for (i, soc_n_val) in enumerate(soc_n_levels)
            parts = String[]
            for method in methods
                mr = results_all[method][i]
                if mr.converged
                    push!(parts, @sprintf("OK %3dit D=%.4f r=%.1e", mr.iterations, mr.D_max, mr.residual_norm))
                else
                    push!(parts, @sprintf("FAIL %3dit D=%.4f r=%.1e", mr.iterations, mr.D_max, mr.residual_norm))
                end
            end
            @printf(io, "%-8.3f | %s | %s | %s\n", soc_n_val, parts[1], parts[2], parts[3])
        end

        println(io, "")
        println(io, "Summary:")
        for method in methods
            results = results_all[method]
            n_converged = count(r -> r.converged, results)
            max_D_final = isempty(results) ? 0.0 : results[end].D_max
            total_iter = sum(r -> r.iterations, results)
            @printf(io, "  %-15s: converged %d/%d, D_max(final)=%.4f, total_iter=%d\n",
                method, n_converged, length(soc_n_levels), max_D_final, total_iter)
        end
    end

    println("\nReport written to: $(report_path)")
    println("=" ^ 70)
    println("Verification complete.")
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
