using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "Jellyroll cohesive per-interface params" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive

    # PE-PCC 必须有非零值
    @test coh.σ_max_pe_pcc > 0
    @test coh.G_c_pe_pcc > 0
    @test coh.K_n_pe_pcc > 0
    @test coh.δ_0_pe_pcc > 0
    @test coh.δ_c_pe_pcc > 0
    @test coh.τ_max_pe_pcc > 0
    @test coh.K_t_pe_pcc > 0

    # δ_c = 2G_c/σ_max 一致性
    @test coh.δ_c_pe_pcc ≈ 2 * coh.G_c_pe_pcc / coh.σ_max_pe_pcc rtol=1e-6
    @test coh.δ_0_pe_pcc ≈ coh.σ_max_pe_pcc / coh.K_n_pe_pcc rtol=1e-6

    # NE-NCC 同上
    @test coh.σ_max_ne_ncc > 0
    @test coh.G_c_ne_ncc > 0
    @test coh.K_n_ne_ncc > 0
    @test coh.δ_c_ne_ncc ≈ 2 * coh.G_c_ne_ncc / coh.σ_max_ne_ncc rtol=1e-6
    @test coh.δ_0_ne_ncc ≈ coh.σ_max_ne_ncc / coh.K_n_ne_ncc rtol=1e-6
end
