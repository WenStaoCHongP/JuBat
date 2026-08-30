using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "CurrentCollector interface fields" begin
    # 2026-08-30 重构：Cohesive 已删除，界面字段挂 CurrentCollector（PCC/NCC 实例）
    cc = JuBat.CurrentCollector()

    # Mode I
    @test hasproperty(cc, :σ_max)
    @test hasproperty(cc, :K_n)
    @test hasproperty(cc, :δ_0)
    @test hasproperty(cc, :G_c)
    @test hasproperty(cc, :δ_c)

    # Mode II
    @test hasproperty(cc, :τ_max)
    @test hasproperty(cc, :K_t)
    @test hasproperty(cc, :δ_0_t)
    @test hasproperty(cc, :G_c_t)
    @test hasproperty(cc, :δ_c_t)

    # BK 指数与界面热阻
    @test hasproperty(cc, :eta)
    @test hasproperty(cc, :h_c0)
    @test hasproperty(cc, :k_air)
    @test hasproperty(cc, :lambda_m)
    @test hasproperty(cc, :beta)
    @test hasproperty(cc, :threshold)

    # Params 不再有 cohesive 字段；Cohesive 类型已删除
    @test !isdefined(JuBat, :Cohesive)
    @test isdefined(JuBat, :CurrentCollector)
end

@testset "CohesiveElement interface_type + host elems" begin
    elem = JuBat.CohesiveElement(
        1,                          # id
        [1, 2, 3, 4],               # nodes
        [1, 2],                     # nodes_bottom
        [4, 3],                     # nodes_top
        1.0,                        # length
        :PE_PCC,                    # interface_type
        10,                         # host_outer_elem
        7                           # host_inner_elem
    )
    @test elem.interface_type == :PE_PCC
    @test elem.host_outer_elem == 10
    @test elem.host_inner_elem == 7

    # 第二种 interface_type 也应工作
    elem2 = JuBat.CohesiveElement(
        2,                          # id
        [5, 6, 7, 8],               # nodes
        [5, 6],                     # nodes_bottom
        [8, 7],                     # nodes_top
        2.0,                        # length
        :NE_NCC,                    # interface_type
        20,                         # host_outer_elem
        15                          # host_inner_elem
    )
    @test elem2.interface_type == :NE_NCC
    @test elem2.host_outer_elem == 20
    @test elem2.host_inner_elem == 15

    # 旧 layer_idx 字段已删除
    @test !hasproperty(elem, :layer_idx)
end

@testset "MechState fields" begin
    # 2026-08-30 重构：CzmInterfaceParams/CzmParamCache 已删除，演化状态聚合 MechState
    @test !isdefined(JuBat, :CzmInterfaceParams)
    @test !isdefined(JuBat, :CzmParamCache)
    mesh = JuBat.CohesiveMesh()
    @test !hasproperty(mesh, :damage_states)   # 损伤态已迁至 MechState
    # MechState 需要真实网格构造（字段长度依网格）——字段检查在 delta_core/c4lite 覆盖
    @test isdefined(JuBat, :MechState)
end

@testset "CzmSubmesh struct" begin
    # 仅验证字段存在，构造完整子网格在 Chunk 3 测试
    @test JuBat.CzmSubmesh === JuBat.CzmSubmesh  # 类型存在性检查
end

@testset "CohesiveMesh czm_submesh field" begin
    mesh = JuBat.CohesiveMesh()
    @test hasproperty(mesh, :czm_submesh)
    @test !hasproperty(mesh, :thermal_to_czm)
    @test hasproperty(mesh, :cohesive_to_thermal)   # v5 新增：反向映射
    @test isnothing(mesh.czm_submesh)
    @test isnothing(mesh.cohesive_to_thermal)
end
