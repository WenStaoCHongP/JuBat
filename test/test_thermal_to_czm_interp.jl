using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using SparseArrays

@testset "build_thermal_to_czm_interp" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = mesh_data.czm_submesh

    M = JuBat.build_thermal_to_czm_interp(case.mesh["thermal2D"], submesh)

    n_czm_node = submesh.mesh.nlen
    n_thermal_node = case.mesh["thermal2D"].nlen
    @test size(M) == (n_czm_node, n_thermal_node)

    # 每行 ≤4 个非零元，行和 = 1（双线性插值 partition of unity）
    for i in 1:n_czm_node
        row = M[i, :]
        @test nnz(row) <= 4
        @test sum(row) ≈ 1.0 atol=1e-10
        @test all(v >= 0 for v in nonzeros(row))
    end

    # 温度场正向插值验证
    T_thermal = collect(1.0:n_thermal_node)
    T_czm = M * T_thermal
    @test length(T_czm) == n_czm_node
    @test !any(isnan, T_czm)

    # 边界节点插值约束
    @test minimum(T_czm) >= minimum(T_thermal) - 1e-10
    @test maximum(T_czm) <= maximum(T_thermal) + 1e-10

    # 拓扑父热单元不包含节点时必须报错，不得回退到最近热节点
    bad_submesh = deepcopy(submesh)
    bad_submesh.mesh.node[1, :] .= 10 .* maximum(abs, case.mesh["thermal2D"].node)
    @test_throws ErrorException JuBat.build_thermal_to_czm_interp(case.mesh["thermal2D"], bad_submesh)
end

@testset "create_czm_mesh 填充 thermal_to_czm" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = mesh_data.czm_submesh
    czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    @test czm_mesh.thermal_to_czm !== nothing
    @test size(czm_mesh.thermal_to_czm, 1) == submesh.mesh.nlen
end
