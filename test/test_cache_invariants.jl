using Test

include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

@testset "distributed thermal geometry is precomputed" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    sim_case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
    sim_case = JuBat.setup_thermal2D_mesh(sim_case, mesh_data)

    @test sim_case.geometry !== nothing
    @test !isempty(sim_case.geometry.layer_weights)
    @test sim_case.geometry.boundary_edges !== nothing
    @test length(sim_case.geometry.boundary_edges.edges) == length(sim_case.geometry.boundary_edges.L_edge)
    @test fieldtype(JuBat.MeshGeometry, :boundary_edges) === JuBat.BoundaryEdgeCache
end

@testset "redundant optional-cache probes stay removed" begin
    thermal_source = read(joinpath(@__DIR__, "..", "src", "ThermalDistributed.jl"), String)
    coupling_source = read(joinpath(@__DIR__, "..", "src", "CallModel.jl"), String)

    @test !occursin("case.geometry !== nothing ? case.geometry.layer_weights", thermal_source)
    @test !occursin("edge_cache = case.geometry !== nothing", thermal_source)
    @test !occursin("hasproperty(case, :I_e_cache)", coupling_source)
    @test !occursin("I_e_prev", coupling_source)
    @test !occursin("hasfield(typeof(geom), :czm_element_map)", thermal_source)
    @test !occursin("hasfield(typeof(geom), :czm_element_map)", coupling_source)
end
