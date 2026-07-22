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

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, nθ_czm=20, gsorder=2)
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
    mesh_data2 = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, nθ_czm=30, gsorder=2)
    submesh2 = mesh_data2.czm_submesh
    czm_mesh2 = JuBat.create_czm_mesh(submesh2, case.mesh["thermal2D"], case.param)
    cache3 = JuBat.ensure_czm_cache(case, czm_mesh2, param_cache)
    @test cache3 !== cache1
    @test cache3.czm_mesh_id == objectid(czm_mesh2)

    # v5: param_cache 变化时也应失效
    # 注意：param_cache.id = objectid(param)，修改 param_dim.cohesive 字段不会
    # 改变 objectid。为了让 id 真正变化，需要重新构造 param 对象。
    # 这里通过修改字段 + NormaliseParam 重建 param 来获得新的 param_cache。
    param_dim.cohesive.σ_max_pe_pcc = param_dim.cohesive.σ_max_pe_pcc * 1.1
    # 重算 param（NormaliseParam 会创建新 Params 实例，objectid 不同）
    case.param = JuBat.NormaliseParam(param_dim)
    param_cache2 = JuBat.compute_czm_params_per_interface(case)
    @test param_cache2.id ≠ param_cache.id
    cache4 = JuBat.ensure_czm_cache(case, czm_mesh2, param_cache2)
    @test cache4 !== cache3
    @test cache4.param_cache_id == param_cache2.id
end
