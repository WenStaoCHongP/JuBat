using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "ensure_czm_cache 失效判据" begin
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
    param_cache = JuBat.compute_czm_params_per_interface(case)

    # 首次调用：cache 为 nothing，触发构建
    cache1 = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    @test cache1 isa JuBat.CZMAssemblyCache
    @test cache1.czm_mesh_id == objectid(case.czm_mesh)
    @test cache1.param_cache_id == param_cache.id

    # 第二次调用：cache 未失效，应返回同一对象
    cache2 = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    @test cache2 === cache1

    # 网格变化时失效
    mesh_data2 = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, czm_enabled=true, gsorder=2)
    submesh2 = mesh_data2.czm_submesh
    czm_mesh2 = JuBat.create_czm_mesh(submesh2, case.mesh["thermal2D"], case.param)
    cache3 = JuBat.ensure_czm_cache(case, czm_mesh2, param_cache)
    @test cache3 !== cache1
    @test cache3.czm_mesh_id == objectid(czm_mesh2)

    # v5: param_cache 变化时也应失效
    # Task 4.4 fix：CzmParamCache.id 现为内容哈希 hash((hash(pe_pcc), hash(ne_ncc)))，
    # 因此只要 CzmInterfaceParams 的字段值变化，id 就会变化。
    # 注意：compute_czm_params_per_interface 读取 case.param（已归一化），
    # 而非 param_dim。所以修改 param_dim.cohesive.σ_max_pe_pcc 后，必须用
    # NormaliseParam(param_dim) 把改动同步进 case.param，否则下一次
    # compute_czm_params_per_interface 仍读旧值——这是 Case (A) 路径下必要的
    # 同步步骤（不是 band-aid）。
    param_dim.cohesive.σ_max_pe_pcc = param_dim.cohesive.σ_max_pe_pcc * 1.1
    case.param = JuBat.NormaliseParam(param_dim)
    param_cache2 = JuBat.compute_czm_params_per_interface(case)
    @test param_cache2.id ≠ param_cache.id
    cache4 = JuBat.ensure_czm_cache(case, czm_mesh2, param_cache2)
    @test cache4 !== cache3
    @test cache4.param_cache_id == param_cache2.id
end
