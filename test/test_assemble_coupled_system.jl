using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "assemble_coupled_system with CzmParamCache" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, nθ_czm=20, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = mesh_data.czm_submesh
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    param_cache = JuBat.compute_czm_params_per_interface(case)
    @test param_cache isa JuBat.CzmParamCache
    @test haskey(param_cache.by_interface, :PE_PCC)
    @test haskey(param_cache.by_interface, :NE_NCC)

    pe = param_cache.by_interface[:PE_PCC]
    ne = param_cache.by_interface[:NE_NCC]
    @test pe.E_eff ≈ case.param.PE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm
    @test ne.E_eff ≈ case.param.NE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm

    u = zeros(2 * case.czm_mesh.nnode)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    K, f, seps, tracts = JuBat.assemble_coupled_system(case.czm_mesh, u, param_cache;
                                                       damage_states=case.czm_mesh.damage_states,
                                                       K_bulk_cached=cache.K_bulk,
                                                       geom_cache=cache.cohesive_geom,
                                                       ws=cache.ws)
    @test size(K, 1) == length(f) == 2 * case.czm_mesh.nnode
    @test !any(isnan, K)
    @test !any(isnan, f)
end
