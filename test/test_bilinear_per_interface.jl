using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "bilinear_* with CzmInterfaceParams" begin
    # 构造 PE-PCC 接面参数
    params = JuBat.CzmInterfaceParams(
        E_eff = 1.0, ν = 0.3, α = 1.5e-5,
        σ_max = 82e6 / 1e10,
        K_n = 2.4e17 / (1e10 / 1e-6),
        δ_0_n = 82e6 / 2.4e17,
        δ_c_n = 2 * 25.3 / 82e6,
        G_c = 25.3,
        τ_max = 82e6 / 1e10,
        K_t = 2.4e17 / (1e10 / 1e-6),
        δ_0_t = 82e6 / 2.4e17,
        δ_c_t = 2 * 25.3 / 82e6,
        G_c_t = 25.3,
        η = 1.45,
        czm_model = "model1",
        h_c0 = 1e7, k_air = 0.026, lambda_m = 70e-9,
        beta = 1.0, threshold = 70e-9,
    )

    D = JuBat.DamageState()

    # 弹性段：δ_n < δ_0_n
    δ_n_small = params.δ_0_n / 2
    T_n, T_t, _, D_new = JuBat.bilinear_traction_state(δ_n_small, 0.0, D, params)
    @test T_n ≈ params.K_n * δ_n_small
    @test T_t ≈ 0.0
    @test D_new.D ≈ 0.0

    # 软化段：δ_0_n < δ_n < δ_c_n
    δ_n_mid = 0.5 * (params.δ_0_n + params.δ_c_n)
    T_n2, _, _, D_new2 = JuBat.bilinear_traction_state(δ_n_mid, 0.0, D, params)
    @test 0 < T_n2 < params.σ_max
    @test 0 < D_new2.D < 1

    # 完全失效：δ_n > δ_c_n
    δ_n_big = 2 * params.δ_c_n
    T_n3, _, _, D_new3 = JuBat.bilinear_traction_state(δ_n_big, 0.0, D, params)
    @test T_n3 ≈ 0.0
    @test D_new3.D ≈ 1.0

    # 切向（Mode II）
    _, T_t2, _, _ = JuBat.bilinear_traction_state(0.0, params.δ_0_t / 2, D, params)
    @test T_t2 ≈ params.K_t * params.δ_0_t / 2

    # NE-NCC 接面参数（不同 K_n）应给出不同结果
    params_ne = JuBat.CzmInterfaceParams(;
        E_eff = params.E_eff, ν = params.ν, α = params.α,
        σ_max = 2 * params.σ_max,
        K_n = 2 * params.K_n,
        δ_0_n = params.δ_0_n,
        δ_c_n = params.δ_c_n / 2,
        G_c = params.G_c,
        τ_max = params.τ_max,
        K_t = params.K_t,
        δ_0_t = params.δ_0_t,
        δ_c_t = params.δ_c_t,
        G_c_t = params.G_c_t,
        η = params.η,
        czm_model = params.czm_model,
        h_c0 = params.h_c0, k_air = params.k_air, lambda_m = params.lambda_m,
        beta = params.beta, threshold = params.threshold,
    )
    T_n_ne, _, _, _ = JuBat.bilinear_traction_state(params_ne.δ_0_n / 2, 0.0, D, params_ne)
    @test T_n_ne ≈ params_ne.K_n * params_ne.δ_0_n / 2
    @test T_n_ne ≠ T_n
end
