# tools/czm_baseline_probe.jl
#
# CZM 模块拆分前的行为基线探针
# 用途：记录当前 CZM 求解器的数值输出，作为回归检测基准
#
# 运行方式: julia tools/czm_baseline_probe.jl
#
# 预期输出格式:
#   BASELINE_START
#   method=basic: converged=X, iterations=N, residual_norm=R, D_max=D, D_mean=M
#   method=load_substep: converged=X, iterations=N, residual_norm=R, D_max=D, D_mean=M
#   method=arc_length: converged=X, iterations=N, residual_norm=R, D_max=D, D_mean=M
#   BASELINE_END
#
# 设计要点:
#   - 使用小扰动 (dT=1K, Δsoc=0.005) 确保 elastic regime 求解器收敛
#   - 使用 cache 提高效率和 BC 一致性
#   - 使用与 Option 默认值匹配的 tol=1e-4, max_iter=200

using LinearAlgebra
using SparseArrays
using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function run_baseline()
    println("============================================================")
    println("CZM Baseline Probe - Pre-Refactoring Snapshot")
    println("============================================================")

    # 1. 准备参数和网格
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.gsorder = 2
    opt.czm_model = "model1"          # Mode I only

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=40, gsorder=2)
    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)

    println("  Nodes: $(czm_mesh.nnode)")
    println("  Bulk Elements: $(size(czm_mesh.bulk_element, 1))")
    println("  Cohesive Elements: $(czm_mesh.n_cohesive)")

    # 2. 计算有效参数
    E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
    @printf("  E_eff = %.6e, ν_eff = %.6f\n", E_eff, ν_eff)
    @printf("  α_eff = %.6e, β_n = %.6e, β_p = %.6e\n", α_eff, β_n, β_p)

    # 3. 构造热化学载荷（小扰动确保 elastic regime 求解器收敛）
    ne = size(czm_mesh.bulk_element, 1)
    dT_elem = fill(0.05 / param_dim.cell.T0, ne)   # 归一化温升 0.05K
    Δsoc_n_elem = fill(0.0002, ne)                  # 极小 SOC 变化
    Δsoc_p_elem = fill(0.0002, ne)

    # 4. 构造 CZM 参数（归一化后的参数，与 update_czm_damage! 一致）
    czm_params = case.param.cohesive
    param = case.param
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)

    # 5. 构建 cache（与实际 Solve.jl 工作流一致）
    cache = JuBat.build_czm_cache(czm_mesh, E_eff, ν_eff, param)

    # 6. 测试三种求解方法
    println("\nBASELINE_START")

    for method in ["basic", "load_substep", "arc_length"]
        # 重建 czm_mesh 以确保每个方法从相同初始状态开始
        czm_mesh_fresh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
        u_prev = zeros(Float64, ndof)

        result, updated_mesh = JuBat.solve_czm_step(
            czm_mesh_fresh, F_ext, E_eff, ν_eff, czm_params, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=200, tol=1e-4, n_load_steps=50,
            iter_method=method, cache=cache
        )

        stats = JuBat.get_damage_statistics(updated_mesh)

        @printf("method=%s: converged=%s, iterations=%d, residual_norm=%.6e, D_max=%.6e, D_mean=%.6e, n_fractured=%d\n",
            method, result.converged, result.iterations, result.residual_norm,
            stats.max_D, stats.mean_D, stats.n_fractured)
    end

    println("BASELINE_END")

    # 7. 本构层测试
    println("\n--- Constitutive Baseline ---")
    state = JuBat.DamageState()
    δ_n, δ_t = 1e-3, 0.5e-3  # 归一化分离位移
    T_n, T_t, D, new_state = JuBat.bilinear_traction_state(δ_n, δ_t, state, czm_params)
    dT = JuBat.bilinear_tangent(δ_n, δ_t, new_state, czm_params)
    @printf("constitutive: T_n=%.6e, T_t=%.6e, D=%.6f\n", T_n, T_t, D)
    @printf("tangent: dT11=%.6e, dT12=%.6e, dT21=%.6e, dT22=%.6e\n", dT[1,1], dT[1,2], dT[2,1], dT[2,2])

    println("\n============================================================")
    println("Baseline probe complete.")
    println("Please record the BASELINE_START ... BASELINE_END output")
    println("for regression comparison after refactoring.")
    println("============================================================")
end

run_baseline()
