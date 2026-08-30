using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# C1（Theory/08 §7.5.1）：塑性关 + K_G→0 + 接触关 ⟹ 严格退回工况 R 线性理论。
# 数值判据：ε*→0 极限下 geo_nl 解趋近线性解；geo_nl=false 两次求解逐位一致。
# 2026-08-30 重构适配：solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; 载荷/状态)。

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
    case.mech = JuBat.MechState(case.czm_mesh)
    return case
end

@testset "C1：K_G→0 极限下 geo_nl 解趋近线性解" begin
    case = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(ndof)
    εtiny = 1e-6
    eig = (dT=fill(εtiny, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.max_iter = 100
    case.opt.czm.tol = 1e-10
    # 线性参照：F_tc 显式（Batch 1 冻结路径）
    ms_lin = JuBat.MechState(czm_mesh)
    r_lin = JuBat.solve_czm_step(
        czm_mesh, ms_lin, case.param, F_ext, case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp)
    # GL 路径：eigenstrain 内嵌（D-B2-1）；ms 复位零初态保证同起点
    ms_geo = JuBat.MechState(czm_mesh)
    case.opt.czm.geo_nonlinear = true
    r_geo = JuBat.solve_czm_step(
        czm_mesh, ms_geo, case.param, F_ext, case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp, eigenstrain=eig)
    @test r_geo.converged && r_lin.converged
    @test isapprox(r_geo.displacement, r_lin.displacement; rtol=1e-4, atol=1e-12)
end

@testset "geo_nl=false 两次求解逐位一致（回归锚）" begin
    case = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (dT=fill(1e-4, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.max_iter = 100
    case.opt.czm.tol = 1e-10
    ms1 = JuBat.MechState(czm_mesh)
    ms2 = JuBat.MechState(czm_mesh)
    r1 = JuBat.solve_czm_step(
        czm_mesh, ms1, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp)
    r2 = JuBat.solve_czm_step(
        czm_mesh, ms2, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp)
    @test r1.displacement == r2.displacement
end
