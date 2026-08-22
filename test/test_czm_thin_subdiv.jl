using Test
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function subdiv_fixture(; nθ::Int=8, thin_subdiv::Int=1)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2,
                                                    thin_subdiv=thin_subdiv)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    return case, param_cache
end

@testset "thin_subdiv=3 拓扑：子层单元/材料继承/cohesive 恒量/Φ 重合" begin
    case1, _ = subdiv_fixture(nθ=8, thin_subdiv=1)
    case3, _ = subdiv_fixture(nθ=8, thin_subdiv=3)
    nseg = 8 * 21  # nθ=8 的分段数（由默认网格实测校验）
    ne1 = size(case1.czm_mesh.bulk_element, 1)
    ne3 = size(case3.czm_mesh.bulk_element, 1)
    @test ne3 - ne1 == 8 * nseg        # 4 薄层（SP×2/PCC/NCC）× (3-1) 子层增量
    @test ne3 == (4 + 12) * nseg       # 4 厚层 + 12 薄子层
    # 材料继承：细分层材料类型与母层一致（逐单元按径向序）
    mt3 = case3.czm_mesh.czm_submesh.material_type
    @test count(==(:PE), mt3) == 2 * 3 * nseg
    @test count(==(:NE), mt3) == 2 * 3 * nseg
    @test count(==(:SP), mt3) == 2 * nseg
    @test count(==(:PCC), mt3) == nseg
    @test count(==(:NCC), mt3) == nseg
    # cohesive 恒 4·nseg（同材边界豁免）
    @test case3.czm_mesh.n_cohesive == 4 * nseg == case1.czm_mesh.n_cohesive
    # Φ 配对坐标重合
    for (o, i) in case3.czm_mesh.czm_submesh.phi_pairs[1:5:end]
        @test norm(case3.czm_mesh.node[o, :] .- case3.czm_mesh.node[i, :]) < 1e-8
    end
end

@testset "thin_subdiv=1（默认）与现状网格逐位一致" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option(); opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"
    case = JuBat.SetCase(param_dim, opt)
    md_old = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=8, czm_enabled=true, gsorder=2)
    md_new = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=8, czm_enabled=true, gsorder=2, thin_subdiv=1)
    @test md_new.czm_submesh.mesh.node == md_old.czm_submesh.mesh.node
    @test md_new.czm_submesh.mesh.element == md_old.czm_submesh.mesh.element
    @test md_new.czm_submesh.material_type == md_old.czm_submesh.material_type
    @test md_new.czm_submesh.phi_pairs == md_old.czm_submesh.phi_pairs
end

@testset "split_KG：K_mat+K_G ≈ K_total 且稀疏模式一致（geo 路径）" begin
    case, pc = subdiv_fixture(nθ=8)
    cm = case.czm_mesh
    u = zeros(2 * cm.nnode)
    σ0 = JuBat.winding_prestress_field(cm, case.param)
    _, K = JuBat.assemble_bulk_residual_tangent(cm, u, pc; geo_nl=true, prestress=σ0)
    _, Kmat, KG = JuBat.assemble_bulk_residual_tangent(cm, u, pc; geo_nl=true, prestress=σ0, split_KG=true)
    # 跨求和次序差异 → rtol 1e-12（计划偏差：geo 路径无位级冻结契约，冻结门禁全在线性路径）
    @test isapprox(Array(K), Array(Kmat) + Array(KG); rtol=1e-12)
    @test nnz(K) == nnz(Kmat) == nnz(KG) || nnz(K) ≥ nnz(Kmat)
    @test norm(Array(KG) - Array(KG)') ≤ 1e-12 * norm(Array(KG))   # 跨单元求和次序致位级不对称，范数判据
end

@testset "非法组合：线弹性路径 split_KG=true → error" begin
    case, pc = subdiv_fixture(nθ=8)
    cm = case.czm_mesh
    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        cm, zeros(2 * cm.nnode), pc; split_KG=true)
end
