using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "Cohesive per-interface normalization" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    # NormaliseParam 由 SetCase 内部调用，这里通过 SetCase 触发
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    coh = case.param.cohesive
    coh_dim = param_dim.cohesive
    scale = param_dim.scale

    # ========== 全 20 字段归一化一致性测试 ==========
    # (field_name, divisor_scale_key) 映射；每条断言验证归一化值 = 量纲值 / 正确除子
    field_divisor = [
        # σ 类 → scale.σ_czm
        ("σ_max_pe_pcc", :σ_czm),
        ("σ_max_ne_ncc", :σ_czm),
        ("τ_max_pe_pcc", :σ_czm),
        ("τ_max_ne_ncc", :σ_czm),
        # δ 类 → scale.δ_czm
        ("δ_0_pe_pcc",   :δ_czm),
        ("δ_c_pe_pcc",   :δ_czm),
        ("δ_0_pe_pcc_t", :δ_czm),
        ("δ_c_pe_pcc_t", :δ_czm),
        ("δ_0_ne_ncc",   :δ_czm),
        ("δ_c_ne_ncc",   :δ_czm),
        ("δ_0_ne_ncc_t", :δ_czm),
        ("δ_c_ne_ncc_t", :δ_czm),
        # G_c 类 → scale.G_czm
        ("G_c_pe_pcc",   :G_czm),
        ("G_c_pe_pcc_t", :G_czm),
        ("G_c_ne_ncc",   :G_czm),
        ("G_c_ne_ncc_t", :G_czm),
        # K 类 → scale.K_czm
        ("K_n_pe_pcc",   :K_czm),
        ("K_t_pe_pcc",   :K_czm),
        ("K_n_ne_ncc",   :K_czm),
        ("K_t_ne_ncc",   :K_czm),
    ]

    for (field, divisor_key) in field_divisor
        divisor = getfield(scale, divisor_key)
        norm_val = getfield(coh, Symbol(field))
        dim_val  = getfield(coh_dim, Symbol(field))
        @test norm_val ≈ dim_val / divisor rtol=1e-6
    end

    # ========== 负向回归测试：错 divisor 必须不成立 ==========
    # 若有人误写 σ_max_pe_pcc / scale.δ_czm（错 divisor），归一化值不应当等于该错误结果
    wrong_σ = coh_dim.σ_max_pe_pcc / scale.δ_czm
    @test !isapprox(coh.σ_max_pe_pcc, wrong_σ; rtol=1e-6)
    # K_n_pe_pcc 不应当等于 K_n_dim / scale.δ_czm（错 divisor）
    wrong_K = coh_dim.K_n_pe_pcc / scale.δ_czm
    @test !isapprox(coh.K_n_pe_pcc, wrong_K; rtol=1e-6)
end
