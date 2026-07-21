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

    # ν（泊松比）— 防 PE/NE 复制粘贴错误
    @test pe.ν ≈ case.param.PE.nu_coat rtol=1e-6
    @test ne.ν ≈ case.param.NE.nu_coat rtol=1e-6

    # G_c — 防后缀拼写错误
    @test pe.G_c ≈ case.param.cohesive.G_c_pe_pcc rtol=1e-6
    @test ne.G_c ≈ case.param.cohesive.G_c_ne_ncc rtol=1e-6

    # K_n — 法向刚度
    @test pe.K_n ≈ case.param.cohesive.K_n_pe_pcc rtol=1e-6
    @test ne.K_n ≈ case.param.cohesive.K_n_ne_ncc rtol=1e-6

    # δ_0_n — 法向初始强度对应的分离
    @test pe.δ_0_n ≈ case.param.cohesive.δ_0_pe_pcc rtol=1e-6
    @test ne.δ_0_n ≈ case.param.cohesive.δ_0_ne_ncc rtol=1e-6

    # τ_max — 切向强度
    @test pe.τ_max ≈ case.param.cohesive.τ_max_pe_pcc rtol=1e-6
    @test ne.τ_max ≈ case.param.cohesive.τ_max_ne_ncc rtol=1e-6

    # CzmParamCache 契约：id 与 param_ref
    @test cache.id == objectid(case.param)
    @test cache.param_ref === case.param
end

@testset "compute_czm_params_per_interface defensive asserts" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    # Tamper with E_coat to trigger assert
    case.param.PE.E_coat = 0.0
    @test_throws AssertionError JuBat.compute_czm_params_per_interface(case)
end
