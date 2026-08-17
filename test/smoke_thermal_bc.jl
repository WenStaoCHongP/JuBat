using Test
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function thermal_bc_fixture()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    sim_case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
    sim_case = JuBat.setup_thermal2D_mesh(sim_case, mesh_data)

    mesh = sim_case.mesh["thermal2D"]
    nnode = size(mesh.node, 1)
    return sim_case, mesh, spzeros(nnode, nnode), zeros(nnode)
end

@testset "Thermal boundary in-place behavior" begin
    sim_case, mesh, K0, F0 = thermal_bc_fixture()

    K = copy(K0)
    F = copy(F0)
    K_result, F_result = JuBat.apply_convection_bc(K, F, mesh, nothing, sim_case)
    @test K_result === K
    @test F_result === F
    @test all(isfinite, nonzeros(K))
    @test all(isfinite, F)

    for cool_method in ("none", "surface", "tab")
        sim_case.opt.cool_method = cool_method
        K = copy(K0)
        F = copy(F0)
        K_result, F_result = JuBat.apply_cool_method(K, F, mesh, sim_case)
        @test K_result === K
        @test F_result === F
        @test all(isfinite, nonzeros(K))
        @test all(isfinite, F)
    end

    @test !isdefined(JuBat, Symbol("apply_convection_bc!"))
    @test !isdefined(JuBat, Symbol("apply_cool_method!"))
end
