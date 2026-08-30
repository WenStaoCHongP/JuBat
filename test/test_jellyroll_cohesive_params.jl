using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 2026-08-30 重构：界面参数挂 CurrentCollector（PCC=PE-PCC 界面，NCC=NE-NCC 界面）

@testset "Jellyroll per-interface params on collectors" begin
    param_dim = JuBat.ChooseCell("Jellyroll")

    # PE-PCC 界面（PCC）必须有非零值
    @test param_dim.PCC.σ_max > 0
    @test param_dim.PCC.G_c > 0
    @test param_dim.PCC.K_n > 0
    @test param_dim.PCC.δ_0 > 0
    @test param_dim.PCC.δ_c > 0
    @test param_dim.PCC.τ_max > 0
    @test param_dim.PCC.K_t > 0

    # δ_c = 2G_c/σ_max 一致性
    @test param_dim.PCC.δ_c ≈ 2 * param_dim.PCC.G_c / param_dim.PCC.σ_max rtol=1e-6
    @test param_dim.PCC.δ_0 ≈ param_dim.PCC.σ_max / param_dim.PCC.K_n rtol=1e-6

    # NE-NCC 界面（NCC）同上
    @test param_dim.NCC.σ_max > 0
    @test param_dim.NCC.G_c > 0
    @test param_dim.NCC.K_n > 0
    @test param_dim.NCC.δ_c ≈ 2 * param_dim.NCC.G_c / param_dim.NCC.σ_max rtol=1e-6
    @test param_dim.NCC.δ_0 ≈ param_dim.NCC.σ_max / param_dim.NCC.K_n rtol=1e-6
end
