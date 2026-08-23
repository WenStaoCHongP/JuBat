using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec v1.5 §3.4：Φ 缝默认完美粘结——phi_pairs 节点合并入网格构造（非 opt-in）。

function merge_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    return case, md
end

@testset "Φ 合并拓扑：真实配对、cohesive 恒量、末端豁免" begin
    case, md = merge_fixture(nθ=8)
    sub = md.czm_submesh
    # v1.5 双网格契约：phi_pairs 保留未合并 .mesh 上的真实节点对；mesh_bonded 承担拓扑合并
    pp = md.czm_submesh.phi_pairs
    @test !isempty(pp)
    nseg = size(case.czm_mesh.bulk_element, 1) ÷ 8
    @test case.czm_mesh.n_cohesive == 4 * nseg          # 合并不产生新 cohesive
    @test all(o != i for (o, i) in pp)
    @test all(isapprox(sub.mesh.node[o, :], sub.mesh.node[i, :]; atol=1e-12) for (o, i) in pp)
    @test all(!(o in sub.phi_keep) && i in sub.phi_keep for (o, i) in pp)
    # nnode 减少量 == 合并对数（重映射后没有节点消失两次）
    # mesh_bonded 节点数 = mesh 原始数 − 合并对数；单元最大索引 = bonded nlen（重映射连续）
    @test md.czm_submesh.mesh_bonded.nlen == md.czm_submesh.mesh.nlen - length(pp)
    @test maximum(md.czm_submesh.mesh_bonded.element) == md.czm_submesh.mesh_bonded.nlen
    @test length(md.czm_submesh.phi_keep) == md.czm_submesh.mesh_bonded.nlen
    @test size(case.czm_mesh.thermal_to_czm, 1) == md.czm_submesh.mesh.nlen  # 插值行数契约（未合并布局）
end

@testset "绑定对位移恒等（装配级）" begin
    case, _ = merge_fixture(nθ=8)
    cm = case.czm_mesh
    pc = JuBat.compute_czm_params_per_interface(case)
    ndof = 2 * cm.nnode
    u = zeros(ndof)
    for n in 1:cm.nnode
        u[2*n-1] = 1e-4 * sin(0.3 * cm.node[n, 1] + 1.0)
        u[2*n] = 1e-4 * cos(0.2 * cm.node[n, 2])
    end
    K, f = JuBat.assemble_coupled_system(cm, u, pc)
    # 合并节点在 K 中表现为同一自由度行/列（无重复索引）——用 f 有限性与 K 对称近似性作烟雾断言
    @test all(isfinite, f) && all(isfinite, K)
    @test norm(Array(K) - Array(K)') ≤ 1e-9 * norm(Array(K))
end

@testset "坐标不重合即拒绝合并（AGENTS 9.7）" begin
    _, md = merge_fixture(nθ=8)
    sub = md.czm_submesh
    bad_mesh = deepcopy(sub.mesh)
    outer, _ = first(sub.phi_pairs)
    bad_mesh.node[outer, 1] += 1e-4
    @test_throws ErrorException JuBat.merge_phi_pairs(bad_mesh, sub.phi_pairs, 2)
end
