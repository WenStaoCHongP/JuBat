include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat
using Printf

"""
网格收敛：nθ_czm ∈ [40, 80, 160] 下 δ_max 稳定性
按 spec v2 §8.3
"""
param_dim = JuBat.ChooseCell("Jellyroll")

results = []
for nθ_czm in [40, 80, 160]
    @printf("Running nθ_czm=%d...\n", nθ_czm)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_model = "model1"
    opt.solveType = "Crank-Nicolson"
    opt.time = [0, 100.0]
    opt.dt = [1e-6, 1.0]

    case = JuBat.SetCase(param_dim, opt)
    nθ_thermal = 80
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ_thermal, nθ_czm=nθ_czm, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    czm_submesh = mesh_data.czm_submesh
    case.czm_mesh = JuBat.create_czm_mesh(czm_submesh, mesh_data.thermal2D, case.param)
    case.czm_param_cache = JuBat.compute_czm_params_per_interface(case)

    t_elapsed = @elapsed result = JuBat.Solve(case)

    sep_n = get(result, "czm separation normal [m]", zeros(1, 1))
    δ_max = maximum(abs.(sep_n))
    D_max = maximum(get(result, "czm D_max", [0.0]))
    n_fractured = maximum(get(result, "czm n_fractured", [0]))

    push!(results, (nθ_czm, δ_max, D_max, n_fractured, t_elapsed))
    @printf("  nθ_czm=%d  δ_max=%.4e  D_max=%.4f  n_fractured=%d  time=%.2fs\n",
            nθ_czm, δ_max, D_max, n_fractured, t_elapsed)
end

println("\n=== 网格收敛摘要 ===")
for r in results
    @printf("nθ_czm=%d: δ_max=%.4e, D_max=%.4f, time=%.2fs\n", r...)
end

δ_40 = results[1][2]
δ_160 = results[3][2]
rel_diff = abs(δ_160 - δ_40) / max(abs(δ_160), 1e-30)
@printf("δ_max(40→160) 相对变化: %.2f%%\n", rel_diff * 100)
if rel_diff < 0.05
    println("✓ 网格收敛满足 < 5% 准则")
else
    println("⚠ 网格未收敛，建议进一步加密")
end

println("\n=== 绝对计时（仅参考，不与旧版本对比）===")
for r in results
    @printf("nθ_czm=%d: %.2fs\n", r[1], r[5])
end
