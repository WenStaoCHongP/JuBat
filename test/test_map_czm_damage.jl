using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "map_czm_damage_to_thermal max 归约" begin
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
    ne_thermal = size(case.mesh["thermal2D"].element, 1)

    # 手动设置部分 damage_states 验证 max 归约
    n_coh = case.czm_mesh.n_cohesive
    for i in 1:min(n_coh, length(case.czm_mesh.damage_states))
        case.czm_mesh.damage_states[i].D = (i % 5) / 4.0   # 0, 0.25, 0.5, 0.75, 1.0 循环
    end

    D_per_thermal = JuBat.map_czm_damage_to_thermal(case.czm_mesh, ne_thermal)
    @test length(D_per_thermal) == ne_thermal
    @test all(0 .<= D_per_thermal .<= 1.0)

    # 构造合成场景：手动指定同一 e_thermal 对应多个 cohesive 单元，验证 max
    czm_mesh_synthetic = case.czm_mesh
    # 把前 3 个 cohesive 单元强制映射到粗热单元 1
    czm_mesh_synthetic.cohesive_to_thermal[1:3] .= 1
    czm_mesh_synthetic.damage_states[1].D = 0.3
    czm_mesh_synthetic.damage_states[2].D = 0.7
    czm_mesh_synthetic.damage_states[3].D = 0.1

    D2 = JuBat.map_czm_damage_to_thermal(czm_mesh_synthetic, ne_thermal)
    @test D2[1] ≈ 0.7   # max(0.3, 0.7, 0.1)

    # 完全断裂值
    czm_mesh_synthetic.damage_states[1].D = 1.0
    D3 = JuBat.map_czm_damage_to_thermal(czm_mesh_synthetic, ne_thermal)
    @test D3[1] ≈ 1.0

    # 未覆盖 CZM 的粗热单元应返回 0
    uncovered = findfirst(e -> !in(e, czm_mesh_synthetic.cohesive_to_thermal), 1:ne_thermal)
    if uncovered !== nothing
        @test D3[uncovered] == 0.0
    end
end
