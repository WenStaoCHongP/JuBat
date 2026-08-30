using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "compute_czm_strain_inputs 按材料类型分发" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = mesh_data.czm_submesh
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    case.czm_param_cache = JuBat.compute_czm_params_per_interface(case)

    ne_czm = size(submesh.mesh.element, 1)
    ne_thermal = size(case.mesh["thermal2D"].element, 1)

    T0 = case.param.cell.T0
    n_thermal_node = case.mesh["thermal2D"].nlen
    T_nodes = T0 .+ 5.0 * sin.(collect(1.0:n_thermal_node))

    cs0_p = case.param.PE.cs0
    cs0_n = case.param.NE.cs0
    variables = Dict{String, Any}(
        "thermal2D element soc_p" => cs0_p .+ 100.0 * cos.(collect(1.0:ne_thermal)),
        "thermal2D element soc_n" => cs0_n .+ 50.0  * sin.(collect(1.0:ne_thermal)),
    )

    out = JuBat.compute_czm_strain_inputs(case, variables, T_nodes)
    @test haskey(out, :dT_czm)
    @test haskey(out, :Δsoc_p_czm)
    @test haskey(out, :Δsoc_n_czm)
    @test !haskey(out, :T_czm_nodes)

    @test length(out.dT_czm) == ne_czm
    @test length(out.Δsoc_p_czm) == ne_czm
    @test length(out.Δsoc_n_czm) == ne_czm

    @test minimum(out.dT_czm) >= -5.5
    @test maximum(out.dT_czm) <=  5.5
    @test !any(isnan, out.dT_czm)

    for e in 1:ne_czm
        mt = submesh.material_type[e]
        e_thermal = submesh.thermal_elem_map[e]
        thermal_nodes = case.mesh["thermal2D"].element[e_thermal, :]
        expected_dT = sum(T_nodes[thermal_nodes]) / 4 - case.param.cell.T0
        @test out.dT_czm[e] ≈ expected_dT atol=1e-12
        if mt == :PE
            expected = variables["thermal2D element soc_p"][e_thermal] - case.param.PE.cs0
            @test out.Δsoc_p_czm[e] ≈ expected atol=1e-6
            @test out.Δsoc_n_czm[e] == 0.0
        elseif mt == :NE
            expected = variables["thermal2D element soc_n"][e_thermal] - case.param.NE.cs0
            @test out.Δsoc_n_czm[e] ≈ expected atol=1e-6
            @test out.Δsoc_p_czm[e] == 0.0
        else
            @test out.Δsoc_p_czm[e] == 0.0
            @test out.Δsoc_n_czm[e] == 0.0
        end
    end
end
