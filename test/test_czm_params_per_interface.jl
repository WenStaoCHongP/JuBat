using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "compute_czm_params_per_interface" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)

    cache = JuBat.compute_czm_params_per_interface(case)

    @test cache isa JuBat.CzmParamCache
    @test haskey(cache.by_interface, :PE_PCC)
    @test haskey(cache.by_interface, :NE_NCC)

    pe = cache.by_interface[:PE_PCC]
    ne = cache.by_interface[:NE_NCC]

    # PE-PCC 用 PE.E_coat（非全栈均一化）
    scale = param_dim.scale
    expected_E_pe = case.param.PE.E_coat * scale.E_coat / scale.σ_czm
    @test pe.E_eff ≈ expected_E_pe rtol=1e-6

    # NE-NCC 用 NE.E_coat
    expected_E_ne = case.param.NE.E_coat * scale.E_coat / scale.σ_czm
    @test ne.E_eff ≈ expected_E_ne rtol=1e-6

    # σ_max 与参数集一致（归一化值）
    @test pe.σ_max ≈ case.param.cohesive.σ_max_pe_pcc rtol=1e-6
    @test ne.σ_max ≈ case.param.cohesive.σ_max_ne_ncc rtol=1e-6
end
