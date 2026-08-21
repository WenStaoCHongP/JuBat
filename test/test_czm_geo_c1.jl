using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# C1（Theory/08 §7.5.1）：塑性关 + K_G→0 + 接触关 ⟹ 严格退回工况 R 线性理论。
# 数值判据：ε*→0 极限下 geo_nl 解趋近线性解；geo_nl=false 与 Batch 1 冻结路径逐位一致。

function build_c1_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    return case, param_cache, cache
end

@testset "C1：K_G→0 极限下 geo_nl 解趋近线性解" begin
    case, param_cache, cache = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(ndof)
    u0 = zeros(ndof)
    εtiny = 1e-6
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0,
           dT=fill(εtiny, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    # 线性参照：F_tc 显式（Batch 1 冻结路径）
    r_lin, _ = JuBat.solve_czm_step(
        czm_mesh, F_ext, param_cache, case.param, u0;
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache)
    # GL 路径：eigenstrain 内嵌（D-B2-1）
    r_geo, _ = JuBat.solve_czm_step(
        czm_mesh, F_ext, param_cache, case.param, u0;
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache,
        geo_nl=true, eigenstrain=eig)
    @test r_geo.converged && r_lin.converged
    @test isapprox(r_geo.displacement, r_lin.displacement; rtol=1e-4, atol=1e-12)
end

@testset "geo_nl=false 与 Batch 1 冻结解逐位一致（回归锚）" begin
    case, param_cache, cache = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0, dT=fill(1e-4, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    r1, _ = JuBat.solve_czm_step(
        czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache)
    r2, _ = JuBat.solve_czm_step(
        czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache,
        geo_nl=false)
    @test r1.displacement == r2.displacement
end
