using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "Cohesive struct per-interface fields" begin
    coh = JuBat.Cohesive()

    # PE-PCC Mode I
    @test hasproperty(coh, :σ_max_pe_pcc)
    @test hasproperty(coh, :K_n_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc)
    @test hasproperty(coh, :G_c_pe_pcc)
    @test hasproperty(coh, :δ_c_pe_pcc)

    # PE-PCC Mode II
    @test hasproperty(coh, :τ_max_pe_pcc)
    @test hasproperty(coh, :K_t_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc_t)
    @test hasproperty(coh, :G_c_pe_pcc_t)
    @test hasproperty(coh, :δ_c_pe_pcc_t)

    # NE-NCC Mode I
    @test hasproperty(coh, :σ_max_ne_ncc)
    @test hasproperty(coh, :K_n_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc)
    @test hasproperty(coh, :G_c_ne_ncc)
    @test hasproperty(coh, :δ_c_ne_ncc)

    # NE-NCC Mode II
    @test hasproperty(coh, :τ_max_ne_ncc)
    @test hasproperty(coh, :K_t_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc_t)
    @test hasproperty(coh, :G_c_ne_ncc_t)
    @test hasproperty(coh, :δ_c_ne_ncc_t)

    # 旧字段已移除
    @test !hasproperty(coh, :σ_max_n)
    @test !hasproperty(coh, :K_n)
    @test !hasproperty(coh, :δ_0_n)
    @test !hasproperty(coh, :G_c_n)
    @test !hasproperty(coh, :δ_c_n)
    @test !hasproperty(coh, :τ_max_t)
    @test !hasproperty(coh, :K_t)
    @test !hasproperty(coh, :δ_0_t)
    @test !hasproperty(coh, :G_c_t)
    @test !hasproperty(coh, :δ_c_t)

    # 保留字段（界面热阻、粘性、BK eta、czm_model）
    @test hasproperty(coh, :eta)
    @test hasproperty(coh, :czm_model)
    @test hasproperty(coh, :h_c0)
    @test hasproperty(coh, :k_air)
    @test hasproperty(coh, :lambda_m)
    @test hasproperty(coh, :beta)
    @test hasproperty(coh, :threshold)
    @test hasproperty(coh, :tau_visc)
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

@testset "CzmInterfaceParams and CzmParamCache" begin
    params_pe_pcc = JuBat.CzmInterfaceParams(;
        E_eff = 1.0e3,
        ν = 0.3,
        α = 2.0e-6,
        σ_max = 50e6,
        K_n = 1.0e17,
        δ_0_n = 1.0e-9,
        δ_c_n = 5.0e-7,
        G_c = 25.0,
        τ_max = 50e6,
        K_t = 1.0e17,
        δ_0_t = 1.0e-9,
        δ_c_t = 5.0e-7,
        G_c_t = 25.0,
        η = 1.45,
        czm_model = "model1",
        h_c0 = 1e7,
        k_air = 0.026,
        lambda_m = 70e-9,
        beta = 1.0,
        threshold = 70e-9,
    )
    @test params_pe_pcc.σ_max == 50e6
    @test params_pe_pcc.czm_model == "model1"
    @test params_pe_pcc.threshold == 70e-9

    # CzmParamCache 含 param_ref 与 id 字段（spec §3.5.2）
    # Task 4.4 fix：id 现为内容哈希（hash(CzmInterfaceParams)）而非 objectid(param)。
    # 此处手工构造，仅检查字段存取；id 取 hash(params_pe_pcc) 以反映新语义。
    param_dim = JuBat.ChooseCell("Jellyroll")
    param = JuBat.NormaliseParam(param_dim)
    content_id = hash(params_pe_pcc)
    cache = JuBat.CzmParamCache(Dict(:PE_PCC => params_pe_pcc), param, content_id)
    @test haskey(cache.by_interface, :PE_PCC)
    @test cache.by_interface[:PE_PCC].E_eff == 1.0e3
    @test cache.id == content_id
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
