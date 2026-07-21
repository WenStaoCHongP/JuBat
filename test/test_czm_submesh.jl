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
    # n_turns = max(1, ceil(Int, theta1 / (2π))) = max(1, ceil(Int, (Rout - Rin - s_total) / s_total))
    n_turns_expected = max(1, ceil(Int, (param_dim.cell.Rout - param_dim.cell.Rin - param_dim.cell.layer) / param_dim.cell.layer))
    @test ne == 8 * n_turns_expected * nθ_czm

    # material_type 取值合法
    @test all(m in (:PE, :PCC, :SP, :NE, :NCC) for m in submesh.material_type)

    # 8 层周期序列：通过实际长度反推 n_segments，独立于实现的下界策略
    ne_actual = length(submesh.material_type)
    n_segments_actual = ne_actual ÷ 8
    layer_first_indices = [(l - 1) * n_segments_actual + 1 for l in 1:8]
    @test submesh.material_type[layer_first_indices] == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]

    # 每个 CZM 体单元映射到一个合法粗热单元
    n_thermal = size(mesh_data.thermal2D.element, 1)
    @test all(1 <= e <= n_thermal for e in submesh.thermal_elem_map)

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
end

@testset "build_czm_submesh 最小 nθ_czm 边界" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    # nθ_czm = 3 是 max(3, nθ_czm) 下界
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, nθ_czm=3, gsorder=2)
    submesh = mesh_data.czm_submesh
    @test submesh isa JuBat.CzmSubmesh
    @test size(submesh.mesh.element, 1) >= 8 * 3  # 至少 1 turn × 3 segments × 8 layers
    @test all(m in (:PE, :PCC, :SP, :NE, :NCC) for m in submesh.material_type)
end
