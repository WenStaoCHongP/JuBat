using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "create_czm_mesh from CzmSubmesh" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true                   # 选择未合并热网格
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, nθ_czm=20, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    submesh = mesh_data.czm_submesh
    @test submesh !== nothing
    czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    @test czm_mesh isa JuBat.CohesiveMesh
    @test czm_mesh.n_cohesive > 0

    # interface_type 取值合法
    for elem in czm_mesh.cohesive_elements
        @test elem.interface_type in (:PE_PCC, :NE_NCC)
        @test elem.host_outer_elem >= 1
        @test elem.host_inner_elem >= 1
    end

    # 数量预期（v3）：每卷绕圈 4 个界面（PE-PCC / PCC-PE / NE-NCC / NCC-NE），每界面 n_segments 个 cohesive 单元
    n_segments_per_turn = size(submesh.mesh.element, 1) ÷ 8 ÷ maximum(submesh.winding_turn)
    n_turns_active = length(unique(submesh.winding_turn))
    n_expected = 4 * n_segments_per_turn * n_turns_active
    # 容差：边界裁剪允许 ±n_segments_per_turn
    @test abs(czm_mesh.n_cohesive - n_expected) <= n_segments_per_turn

    # czm_submesh 字段已设置
    @test czm_mesh.czm_submesh === submesh

    # thermal_to_czm 字段已设置（在 Chunk 5 中真正填充，此处先 nothing）
    @test hasfield(JuBat.CohesiveMesh, :thermal_to_czm)

    # cohesive_to_thermal 长度 = n_cohesive，值合法
    @test length(czm_mesh.cohesive_to_thermal) == czm_mesh.n_cohesive
    n_thermal = size(case.mesh["thermal2D"].element, 1)
    @test all(1 <= e <= n_thermal for e in czm_mesh.cohesive_to_thermal)

    # 节点复制 + 重写外层 bulk 正确性自检（spec §4.3）
    for coh in czm_mesh.cohesive_elements
        n_a, n_b, n_b_copy, n_a_copy = coh.nodes
        # 副本节点坐标与原节点一致
        @test czm_mesh.node[n_a, :] ≈ czm_mesh.node[n_a_copy, :] atol=1e-12
        @test czm_mesh.node[n_b, :] ≈ czm_mesh.node[n_b_copy, :] atol=1e-12
        # 4 节点不重复
        @test length(unique(coh.nodes)) == 4
    end

    # 外层 bulk 单元的共边位置必须是副本节点
    for coh in czm_mesh.cohesive_elements
        outer_nodes = czm_mesh.bulk_element[coh.host_outer_elem, :]
        n_a, n_b, n_b_copy, n_a_copy = coh.nodes
        # 副本 n_a_copy/n_b_copy 必须在外层单元中
        @test n_a_copy in outer_nodes
        @test n_b_copy in outer_nodes
        # 原节点 n_a/n_b 不应在外层单元中（已被副本替换）
        @test !(n_a in outer_nodes)
        @test !(n_b in outer_nodes)
    end
end
