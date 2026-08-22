using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "jellyroll_collector_seed_mesh 扩展 czm_submesh (v3)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")

    # 默认不构造力学子网格
    mesh_data_default = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    @test mesh_data_default.czm_submesh === nothing

    # 启用分层 Q4 子网格；周向分段直接继承热网格
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, czm_enabled=true, gsorder=2)
    submesh = mesh_data.czm_submesh
    @test submesh isa JuBat.CzmSubmesh
    @test submesh.mesh.type == "Q4"
    ne = size(submesh.mesh.element, 1)

    n_thermal = size(mesh_data.thermal2D.element, 1)
    @test ne == 8 * n_thermal
    @test submesh.mesh.nlen == 9 * (n_thermal + 1)

    # material_type 取值合法
    @test all(m in (:PE, :PCC, :SP, :NE, :NCC) for m in submesh.material_type)

    # 8 层周期序列：通过实际长度反推 n_segments，独立于实现的下界策略
    ne_actual = length(submesh.material_type)
    n_segments_actual = ne_actual ÷ 8
    layer_first_indices = [(l - 1) * n_segments_actual + 1 for l in 1:8]
    @test submesh.material_type[layer_first_indices] == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]

    # 每层力学单元严格继承热单元父拓扑
    for layer in 1:8
        ids = ((layer - 1) * n_thermal + 1):(layer * n_thermal)
        @test submesh.thermal_elem_map[ids] == collect(1:n_thermal)
    end

    # winding_turn 独立检查（不重新实现 floor 公式）
    # 沿周向（第 1 层，最内）winding_turn 应从 1 单调非递减
    n_segments_total = length(submesh.winding_turn) ÷ 8
    layer1_winding = submesh.winding_turn[1:n_segments_total]
    @test layer1_winding[1] == 1
    @test layer1_winding[end] >= layer1_winding[1]
    @test issorted(layer1_winding)
    # 沿径向（同一 segment 的 8 层）winding_turn 应非递减，且相邻层至多差 1
    mid_seg = n_segments_total ÷ 2
    radial_profile = [submesh.winding_turn[(l - 1) * n_segments_total + mid_seg] for l in 1:8]
    @test issorted(radial_profile)
    @test all(diff(radial_profile) .<= 1)

    # 力学 Φ 配对继承热网格配对计数；v1.5 合并后为 bonded 索引下的 (i,i)（完美粘结记录）
    @test length(submesh.phi_pairs) == length(mesh_data.interface_pairs)
    for (outer_node, inner_node) in submesh.phi_pairs
        @test outer_node == inner_node
        @test outer_node <= submesh.mesh_bonded.nlen
    end
    @test submesh.mesh_bonded.nlen == submesh.mesh.nlen - length(submesh.phi_pairs)
end
