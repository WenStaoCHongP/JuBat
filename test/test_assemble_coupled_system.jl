using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 2026-08-30 重构：界面参数直读 param.PCC/param.NCC；装配缓存为 czm_mesh 惰性字段。
@testset "assemble_coupled_system with direct param" begin
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
    case.mech = JuBat.MechState(case.czm_mesh)

    # 归一化后界面参数就位（σ*、K*、δ* 与锚点一致）
    @test case.param.PCC.σ_max ≈ 1.0 atol=1e-12          # PE-PCC 为锚定界面（σ_czm = PCC.σ_max_dim）
    @test case.param.PCC.δ_c ≈ 1.0 atol=1e-12            # 锚定界面 δ_c* ≡ 1
    @test case.param.NCC.σ_max > 0
    @test case.param.NCC.δ_0 > 0

    u = zeros(2 * case.czm_mesh.nnode)
    K, f, seps, tracts = JuBat.assemble_coupled_system(case.czm_mesh, u, case.param;
                                                       damage_states=case.mech.damage_states,
                                                       K_bulk_cached=JuBat.bulk_stiffness(case.czm_mesh, case.param))
    @test size(K, 1) == length(f) == 2 * case.czm_mesh.nnode
    @test !any(isnan, K)
    @test !any(isnan, f)
end

# 原 test_thermal_resistance_disabled.jl 的唯一实质断言（spec v2 §2.4）：
# czm_enabled=true 时热求解仍使用合并网格（连续径向导热路径）
@testset "czm_enabled uses merged thermal mesh" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.czm.enabled = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    @test case.mesh["thermal2D"] === mesh_data.thermal2D_merged
end
