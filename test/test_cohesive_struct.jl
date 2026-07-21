using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "Cohesive struct per-interface fields" begin
    coh = JuBat.Cohesive()

    # PE-PCC Mode I
    @test hasproperty(coh, :σ_max_pe_pcc)
    @test hasproperty(coh, :K_n_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc)
    @test hasproperty(coh, :G_c_pe_pcc)
    @test hasproperty(coh, :δ_c_pe_pcc)

    # PE-PCC Mode II
    @test hasproperty(coh, :τ_max_pe_pcc)
    @test hasproperty(coh, :K_t_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc_t)
    @test hasproperty(coh, :G_c_pe_pcc_t)
    @test hasproperty(coh, :δ_c_pe_pcc_t)

    # NE-NCC Mode I
    @test hasproperty(coh, :σ_max_ne_ncc)
    @test hasproperty(coh, :K_n_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc)
    @test hasproperty(coh, :G_c_ne_ncc)
    @test hasproperty(coh, :δ_c_ne_ncc)

    # NE-NCC Mode II
    @test hasproperty(coh, :τ_max_ne_ncc)
    @test hasproperty(coh, :K_t_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc_t)
    @test hasproperty(coh, :G_c_ne_ncc_t)
    @test hasproperty(coh, :δ_c_ne_ncc_t)

    # 旧字段已移除
    @test !hasproperty(coh, :σ_max_n)
    @test !hasproperty(coh, :K_n)
    @test !hasproperty(coh, :δ_0_n)
    @test !hasproperty(coh, :G_c_n)
    @test !hasproperty(coh, :δ_c_n)
    @test !hasproperty(coh, :τ_max_t)
    @test !hasproperty(coh, :K_t)
    @test !hasproperty(coh, :δ_0_t)
    @test !hasproperty(coh, :G_c_t)
    @test !hasproperty(coh, :δ_c_t)

    # 保留字段（界面热阻、粘性、BK eta、czm_model）
    @test hasproperty(coh, :eta)
    @test hasproperty(coh, :czm_model)
    @test hasproperty(coh, :h_c0)
    @test hasproperty(coh, :k_air)
    @test hasproperty(coh, :lambda_m)
    @test hasproperty(coh, :beta)
    @test hasproperty(coh, :threshold)
    @test hasproperty(coh, :tau_visc)
end

@testset "CohesiveElement interface_type + host elems" begin
    elem = JuBat.CohesiveElement(
        1,                          # id
        [1, 2, 3, 4],               # nodes
        [1, 2],                     # nodes_bottom
        [4, 3],                     # nodes_top
        1.0,                        # length
        :PE_PCC,                    # interface_type
        10,                         # host_outer_elem
        7                           # host_inner_elem
    )
    @test elem.interface_type == :PE_PCC
    @test elem.host_outer_elem == 10
    @test elem.host_inner_elem == 7

    # 旧 layer_idx 字段已删除
    @test !hasproperty(elem, :layer_idx)
end
