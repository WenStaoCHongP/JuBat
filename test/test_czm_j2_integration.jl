using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# C2-lite 集成（spec §7 Batch 3）：强本征应变下 PCC/NCC 屈服、塑性关退回 C1、
# 状态累计与不可逆性、缺参显式报错。夹具同 C1（归一化网格）。
# 2026-08-30 重构适配：solve_czm_step(czm_mesh, ms, param, F_ext, czm_opt; 载荷/状态)；
# 复用同一 ms 时第二次求解自动从上次收敛位移出发（收敛提交语义）。

function build_j2_fixture(; nθ::Int=8)
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

@testset "塑性关（geo 关/开）与既有路径逐位一致（回归锚）" begin
    case = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (dT=fill(1e-6, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.max_iter = 100
    case.opt.czm.tol = 1e-10
    ms1 = JuBat.MechState(czm_mesh)
    r1 = JuBat.solve_czm_step(czm_mesh, ms1, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp)
    ms2 = JuBat.MechState(czm_mesh)
    case.opt.czm.geo_nonlinear = true
    r2 = JuBat.solve_czm_step(czm_mesh, ms2, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp, eigenstrain=eig)
    @test r1.converged && r2.converged
    @test isapprox(r2.displacement, r1.displacement; rtol=1e-4, atol=1e-12)  # K_G→0 极限（同 C1 判据）
end

@testset "C2-lite：强压缩本征应变下 PCC/NCC 屈服、其余层 κ=0" begin
    case = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    # 层间失配驱动（与 unit_czm_eigenstrain 同机理）。α/β 分层化（2026-08-29）后箔
    # 本征应变归零，±0.3 失配传到硬箔的应力低于屈服，需 Δsoc=−1.0 才超 Cu 箔屈服应变
    # （合成驱动，非物理 soc 区间；探测：−0.3 → κ=0，−1.0 → κ≈9.4e-3）
    eig = (dT=zeros(ne), Δsn=fill(-1.0, ne), Δsp=fill(-1.0, ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    case.opt.czm.geo_nonlinear = true
    case.opt.czm.j2_plasticity = true
    ms = JuBat.MechState(czm_mesh)
    r = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig, mech_state=states)
    @test r.converged
    # 提交（复刻生产 D-B3-2 路径）
    JuBat.assemble_coupled_system(czm_mesh, r.displacement, case.param;
        damage_states=ms.damage_states,
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    mt = czm_mesh.czm_submesh.material_type
    κ_foil = [states[e, g].kappa for e in 1:ne, g in 1:4 if mt[e] in (:PCC, :NCC)]
    κ_soft = [states[e, g].kappa for e in 1:ne, g in 1:4 if !(mt[e] in (:PCC, :NCC))]
    @test maximum(κ_foil) > 0.0
    @test all(iszero, κ_soft)
end

@testset "状态不可逆：同载荷重解后提交的 κ 不再增长" begin
    case = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    # 层间失配驱动：Δsn=−1.0/Δsp=+1.0 反号放电方向，箔零本征应变 ⟹ 失配超箔屈服
    # （β 分层化重标定：±0.3 → κ=0，±1.0 → κ≈7e-4）
    eig = (dT=zeros(ne), Δsn=fill(-1.0, ne), Δsp=fill(1.0, ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    case.opt.czm.geo_nonlinear = true
    case.opt.czm.j2_plasticity = true
    ms = JuBat.MechState(czm_mesh)
    r1 = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig, mech_state=states)
    JuBat.assemble_coupled_system(czm_mesh, r1.displacement, case.param;
        damage_states=ms.damage_states,
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    κ_after1 = maximum(s.kappa for s in states)
    # 同载荷以 ms.u_prev 重解（生产时序），再提交：理想塑性下 κ 不应显著增长（已屈服面停留）
    r2 = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm;
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig, mech_state=states)
    JuBat.assemble_coupled_system(czm_mesh, r2.displacement, case.param;
        damage_states=ms.damage_states,
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    κ_after2 = maximum(s.kappa for s in states)
    @test κ_after1 > 0.0
    @test κ_after2 ≥ κ_after1 * (1 - 1e-8)   # 数值级持平，不得回退
    @test isapprox(κ_after2, κ_after1; rtol=1e-6, atol=1e-12)
end

@testset "缺参/组合非法必须报错（AGENTS 9.7 / D-B3-1）" begin
    case = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (dT=fill(-1e-3, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    # plasticity 需 geo_nl
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), case.param; plasticity=true, mech_state=states)
    # plasticity 无 mech_state
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), case.param; geo_nl=true, eigenstrain=eig, plasticity=true)
    # Batch 5 起 arc_length + geo_nl 可路由（Crisfield geo 弧长，不再显式拒绝）
    case.opt.czm.iter_method = "arc_length"
    case.opt.czm.geo_nonlinear = true
    ms = JuBat.MechState(czm_mesh)
    r_arc = JuBat.solve_czm_step(
        czm_mesh, ms, case.param, zeros(ndof), case.opt.czm;
        eigenstrain=eig)
    @test r_arc isa JuBat.CZMResult
end
