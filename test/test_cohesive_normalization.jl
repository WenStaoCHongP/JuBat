using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 2026-08-30 重构：界面字段挂 CurrentCollector，归一化逐字段显式
@testset "Cohesive per-interface normalization on collectors" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    # NormaliseParam 由 SetCase 内部调用，这里通过 SetCase 触发
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    pcc = case.param.PCC
    ncc = case.param.NCC
    pcc_dim = param_dim.PCC
    ncc_dim = param_dim.NCC
    scale = param_dim.scale

    # ========== 归一化一致性：归一化值 = 量纲值 / 正确除子 ==========
    field_divisor = [
        # σ 类 → scale.σ_czm
        ("σ_max", :σ_czm), ("τ_max", :σ_czm),
        # δ 类 → scale.δ_czm
        ("δ_0", :δ_czm), ("δ_c", :δ_czm), ("δ_0_t", :δ_czm), ("δ_c_t", :δ_czm),
        # G_c 类 → scale.G_czm
        ("G_c", :G_czm), ("G_c_t", :G_czm),
        # K 类 → scale.K_czm
        ("K_n", :K_czm), ("K_t", :K_czm),
    ]

    for (field, divisor_key) in field_divisor
        divisor = getfield(scale, divisor_key)
        @test getfield(pcc, Symbol(field)) ≈ getfield(pcc_dim, Symbol(field)) / divisor rtol=1e-6
        @test getfield(ncc, Symbol(field)) ≈ getfield(ncc_dim, Symbol(field)) / divisor rtol=1e-6
    end

    # ========== 负向回归测试：错 divisor 必须不成立 ==========
    wrong_σ = pcc_dim.σ_max / scale.δ_czm
    @test !isapprox(pcc.σ_max, wrong_σ; rtol=1e-6)
    wrong_K = pcc_dim.K_n / scale.δ_czm
    @test !isapprox(pcc.K_n, wrong_K; rtol=1e-6)
end
