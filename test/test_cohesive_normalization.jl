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
    scale = param_dim.scale

    # 归一化一致性：σ_max_pe_pcc_norm ≈ σ_max_pe_pcc_dim / scale.σ_czm
    # 注意：scale.σ_czm 由 σ_max_pe_pcc 派生（=σ_max_pe_pcc），所以比值应精确 = 1
    expected = param_dim.cohesive.σ_max_pe_pcc / scale.σ_czm
    @test coh.σ_max_pe_pcc ≈ expected rtol=1e-6

    # 同样验证 G_c / K / δ
    @test coh.G_c_pe_pcc ≈ param_dim.cohesive.G_c_pe_pcc / scale.G_czm rtol=1e-6
    @test coh.K_n_pe_pcc ≈ param_dim.cohesive.K_n_pe_pcc / scale.K_czm rtol=1e-6
    @test coh.δ_c_pe_pcc ≈ param_dim.cohesive.δ_c_pe_pcc / scale.δ_czm rtol=1e-6

    # NE-NCC 同样验证
    @test coh.σ_max_ne_ncc ≈ param_dim.cohesive.σ_max_ne_ncc / scale.σ_czm rtol=1e-6
    @test coh.G_c_ne_ncc ≈ param_dim.cohesive.G_c_ne_ncc / scale.G_czm rtol=1e-6
end
