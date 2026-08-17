include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test

@testset "create_unit_czm_strip topology" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    param = case.param

    czm_mesh, meta = JuBat.create_unit_czm_strip(param; y0=1.0, gsorder=2)

    @test size(czm_mesh.bulk_element, 1) == 8
    @test czm_mesh.n_cohesive == 4
    types = [e.interface_type for e in czm_mesh.cohesive_elements]
    @test count(==(:PE_PCC), types) == 2
    @test count(==(:NE_NCC), types) == 2
    @test length(meta.layer_materials) == 8
    @test meta.layer_materials == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    @test length(meta.bottom_nodes) == 2
    @test length(meta.top_nodes_after_czm) == 2
    @test length(meta.pcc_nodes) == 4
    @test length(meta.ncc_nodes) == 4
    @test Set(meta.pcc_nodes) == Set(Int.(czm_mesh.bulk_element[2, :]))
    @test Set(meta.ncc_nodes) == Set(Int.(czm_mesh.bulk_element[6, :]))
    # 底边 y 最小、顶边 y 最大
    @test all(czm_mesh.node[n, 2] ≈ meta.y_interfaces[1] for n in meta.bottom_nodes)
    @test maximum(czm_mesh.node[meta.top_nodes_after_czm, 2]) ≈ meta.y_interfaces[end] atol=1e-12
end
