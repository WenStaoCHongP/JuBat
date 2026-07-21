using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "jellyroll_collector_seed_mesh 扩展 czm_submesh (v3)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")

    # 默认：nθ_czm=nothing，czm_submesh 应为 nothing（向后兼容）
    mesh_data_default = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    @test mesh_data_default.czm_submesh === nothing

    # 启用分层 Q4 子网格
    nθ_czm = 40
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, nθ_czm=nθ_czm, gsorder=2)
    submesh = mesh_data.czm_submesh
    @test submesh isa JuBat.CzmSubmesh
    @test submesh.mesh.type == "Q4"
    ne = size(submesh.mesh.element, 1)

    # 总单元数 = 8 层 × n_turns × nθ_czm
    # n_turns 从螺旋几何直接计算（与 thermal mesh 解耦）
    # b = s_total / (2π), theta1 = (Rout - Rin - s_total) * (2π) / s_total
    # n_turns = max(1, round(Int, theta1 / (2π))) = max(1, round(Int, (Rout - Rin - s_total) / s_total))
    n_turns_expected = max(1, round(Int, (param_dim.cell.Rout - param_dim.cell.Rin - param_dim.cell.layer) / param_dim.cell.layer))
    @test ne == 8 * n_turns_expected * nθ_czm

    # material_type 取值合法
    @test all(m in (:PE, :PCC, :SP, :NE, :NCC) for m in submesh.material_type)

    # 8 层周期序列：build_czm_submesh 采用 layer-outer / segment-inner 遍历，
    # 故同一 segment 索引上跨 8 层的材料序列应为 [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    # 每层段数 = n_segments，seg=1 在每层中的索引为 [1, n_segments+1, 2*n_segments+1, ...]
    n_turns_expected = max(1, round(Int, (param_dim.cell.Rout - param_dim.cell.Rin - param_dim.cell.layer) / param_dim.cell.layer))
    n_segments_per_turn_expected = max(3, nθ_czm)
    n_segments_expected = n_turns_expected * n_segments_per_turn_expected
    indices_seg1_across_layers = [(l - 1) * n_segments_expected + 1 for l in 1:8]
    @test submesh.material_type[indices_seg1_across_layers] == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]

    # 每个 CZM 体单元映射到一个合法粗热单元
    n_thermal = size(mesh_data.thermal2D.element, 1)
    @test all(1 <= e <= n_thermal for e in submesh.thermal_elem_map)

    # winding_turn 与径向位置一致
    s_total = param_dim.cell.layer
    a = param_dim.cell.Rin
    for e in 1:length(submesh.winding_turn)
        n1 = submesh.mesh.element[e, 1]
        n3 = submesh.mesh.element[e, 3]
        r = 0.5 * (hypot(submesh.mesh.node[n1, 1], submesh.mesh.node[n1, 2]) +
                   hypot(submesh.mesh.node[n3, 1], submesh.mesh.node[n3, 2]))
        expected_turn = max(1, Int(floor((r - a) / s_total)) + 1)
        @test submesh.winding_turn[e] == expected_turn
    end
end
