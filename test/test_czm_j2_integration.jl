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
    case.mesh["czm"] = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    return case
end

@testset "material point response keeps S local and returns trial Cauchy state" begin
    committed = JuBat.PlasticState()
    S, C_ep, trial = JuBat.material_point_response(
        (0.1, 0.2, 0.3), (1.0, 0.0, 0.0, 1.0),
        Matrix{Float64}(I, 3, 3), 0.0, committed;
        prestress=(1.0, 2.0, 3.0))
    @test S == (1.1, 2.2, 3.3)
    @test C_ep == Matrix{Float64}(I, 3, 3)
    @test trial.eps_p == committed.eps_p
    @test trial.kappa == committed.kappa
    @test trial.sigma == S
end

@testset "PlasticState carries committed Cauchy stress" begin
    state = JuBat.PlasticState()
    @test state.eps_p == (0.0, 0.0, 0.0)
    @test state.kappa == 0.0
    @test all(isnan, state.sigma)

    S = (2.0, 3.0, 0.5)
    F = [1.1 0.2; -0.1 0.9]
    Smat = [S[1] S[3]; S[3] S[2]]
    expected = F * Smat * transpose(F) / det(F)
    sigma = JuBat.cauchy_from_second_piola(
        S, F[1, 1], F[1, 2], F[2, 1], F[2, 2])
    @test collect(sigma) ≈
          [expected[1, 1], expected[2, 2], expected[1, 2]] rtol=1e-14
    @test_throws ErrorException JuBat.cauchy_from_second_piola(
        S, 1.0, 0.0, 0.0, 0.0)
end

@testset "bulk assembly writes trial states without mutating committed states" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    trial_u = zeros(ndof)
    for n in 1:size(czm_mesh.node, 1)
        trial_u[2n - 1] = 0.02 * czm_mesh.node[n, 1]
    end
    eig = (dT=zeros(ne), Δsn=zeros(ne), Δsp=zeros(ne))
    committed = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    committed_before = JuBat.clone_plastic_states(committed)
    trial = JuBat.clone_plastic_states(committed)
    ms = JuBat.MechState(czm_mesh)
    JuBat.assemble_coupled_system(czm_mesh, trial_u, case.param; czm_workspace=czm_workspace,
        damage_states=ms.damage_states,
        geo_nl=true, eigenstrain=eig, plasticity=true,
        committed_plastic_states=committed,
        trial_plastic_states=trial)
    @test all(committed[e, g].eps_p == committed_before[e, g].eps_p &&
              committed[e, g].kappa == committed_before[e, g].kappa &&
              isequal(committed[e, g].sigma, committed_before[e, g].sigma)
              for e in 1:ne, g in 1:4)
    @test any(s -> s.kappa > 0.0, trial)
    @test all(s -> all(isfinite, s.sigma), trial)
end

@testset "塑性关（geo 关/开）与既有路径逐位一致（回归锚）" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    eig = (dT=fill(1e-6, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.max_iter = 100
    case.opt.czm.tol = 1e-10
    ms1 = JuBat.MechState(czm_mesh)
    r1 = JuBat.solve_czm_step(czm_mesh, ms1, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp)
    ms2 = JuBat.MechState(czm_mesh)
    case.opt.czm.geo_nonlinear = true
    r2 = JuBat.solve_czm_step(czm_mesh, ms2, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp, eigenstrain=eig)
    @test r1.converged && r2.converged
    @test isapprox(r2.displacement, r1.displacement; rtol=1e-4, atol=1e-12)  # K_G→0 极限（同 C1 判据）
end

@testset "C2-lite：强压缩本征应变下 PCC/NCC 屈服、其余层 κ=0" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    # 层间失配驱动（与 unit_czm_eigenstrain 同机理）。α/β 分层化（2026-08-29）后箔
    # 本征应变归零，±0.3 失配传到硬箔的应力低于屈服，需 Δsoc=−1.0 才超 Cu 箔屈服应变
    # （合成驱动，非物理 soc 区间；探测：−0.3 → κ=0，−1.0 → κ≈9.4e-3）
    eig = (dT=zeros(ne), Δsn=fill(-1.0, ne), Δsp=fill(-1.0, ne))
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    case.opt.czm.geo_nonlinear = true
    case.opt.czm.j2_plasticity = true
    ms = JuBat.MechState(czm_mesh)
    r = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig)
    @test r.converged
    @test ms.plastic_states !== nothing
    @test all(s -> all(isfinite, s.sigma), ms.plastic_states)
    mt = czm_mesh.czm_submesh.material_type
    κ_foil = [ms.plastic_states[e, g].kappa for e in 1:ne, g in 1:4 if mt[e] in (:PCC, :NCC)]
    κ_soft = [ms.plastic_states[e, g].kappa for e in 1:ne, g in 1:4 if !(mt[e] in (:PCC, :NCC))]
    @test maximum(κ_foil) > 0.0
    @test all(iszero, κ_soft)
end

@testset "状态不可逆：同载荷重解后提交的 κ 不再增长" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    # 层间失配驱动：Δsn=−1.0/Δsp=+1.0 反号放电方向，箔零本征应变 ⟹ 失配超箔屈服
    # （β 分层化重标定：±0.3 → κ=0，±1.0 → κ≈7e-4）
    eig = (dT=zeros(ne), Δsn=fill(-1.0, ne), Δsp=fill(1.0, ne))
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    case.opt.czm.geo_nonlinear = true
    case.opt.czm.j2_plasticity = true
    ms = JuBat.MechState(czm_mesh)
    r1 = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig)
    @test r1.converged
    @test ms.plastic_states !== nothing
    κ_after1 = maximum(s.kappa for s in ms.plastic_states)
    # 同载荷以 ms.u_prev 重解（生产时序），再提交：理想塑性下 κ 不应显著增长（已屈服面停留）
    r2 = JuBat.solve_czm_step(czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        eigenstrain=eig)
    @test r2.converged
    κ_after2 = maximum(s.kappa for s in ms.plastic_states)
    @test κ_after1 > 0.0
    @test κ_after2 ≥ κ_after1 * (1 - 1e-8)   # 数值级持平，不得回退
    @test isapprox(κ_after2, κ_after1; rtol=1e-6, atol=1e-12)
end

@testset "失败求解不创建或覆盖 committed plastic state" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    eig = (dT=zeros(ne), Δsn=fill(-1.0, ne), Δsp=fill(-1.0, ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.geo_nonlinear = true
    case.opt.czm.j2_plasticity = true
    case.opt.czm.max_iter = 1
    case.opt.czm.tol = 1e-16
    ms_fail = JuBat.MechState(czm_mesh)
    r_fail = JuBat.solve_czm_step(
        czm_mesh, ms_fail, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        eigenstrain=eig)
    @test !r_fail.converged
    @test ms_fail.plastic_states === nothing
    @test all(iszero, ms_fail.u_prev)

    ms = JuBat.MechState(czm_mesh)
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    r_ok = JuBat.solve_czm_step(
        czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        eigenstrain=eig)
    @test r_ok.converged
    committed_before = JuBat.clone_plastic_states(ms.plastic_states)
    eig_retry = (dT=zeros(ne), Δsn=fill(1.0, ne), Δsp=fill(-1.0, ne))
    case.opt.czm.iter_method = "basic"
    case.opt.czm.max_iter = 1
    case.opt.czm.tol = 1e-16
    r_retry = JuBat.solve_czm_step(
        czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        eigenstrain=eig_retry)
    @test !r_retry.converged
    @test all(ms.plastic_states[e, g].eps_p == committed_before[e, g].eps_p &&
              ms.plastic_states[e, g].kappa == committed_before[e, g].kappa &&
              ms.plastic_states[e, g].sigma == committed_before[e, g].sigma
              for e in axes(ms.plastic_states, 1), g in axes(ms.plastic_states, 2))
end

@testset "basic/load-substep/geo arc-length commit the same state contract" begin
    for method in ("basic", "load_substep", "arc_length")
        case = build_j2_fixture()
        czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
        ne = size(czm_mesh.bulk_element, 1)
        eig = (dT=fill(1e-6, ne), Δsn=zeros(ne), Δsp=zeros(ne))
        case.opt.czm.iter_method = method
        case.opt.czm.geo_nonlinear = true
        case.opt.czm.j2_plasticity = true
        case.opt.czm.max_iter = 100
        case.opt.czm.tol = method == "basic" ? 1e-10 : 1e-8
        case.opt.czm.load_steps = 10
        ms = JuBat.MechState(czm_mesh)
        r = JuBat.solve_czm_step(
            czm_mesh, ms, case.param, zeros(2 * size(czm_mesh.node, 1)), case.opt.czm; czm_workspace=czm_workspace,
            dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
            eigenstrain=eig)
        @test r.converged
        @test ms.plastic_states !== nothing
        @test all(s -> all(isfinite, s.sigma), ms.plastic_states)
        mt = czm_mesh.czm_submesh.material_type
        @test all(ms.plastic_states[e, g].kappa >= 0.0
                  for e in 1:ne, g in 1:4 if mt[e] in (:PCC, :NCC))
        @test all(iszero(ms.plastic_states[e, g].kappa)
                  for e in 1:ne, g in 1:4 if !(mt[e] in (:PCC, :NCC)))
    end
end

@testset "缺参/组合非法必须报错（AGENTS 9.7 / D-B3-1）" begin
    case = build_j2_fixture()
    czm_mesh = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * size(czm_mesh.node, 1)
    eig = (dT=fill(-1e-3, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    # plasticity 需 geo_nl
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), case.param; czm_workspace=czm_workspace, plasticity=true)
    # plasticity 无 committed/trial states
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), case.param; czm_workspace=czm_workspace, geo_nl=true, eigenstrain=eig, plasticity=true)
    # Batch 5 起 arc_length + geo_nl 可路由（Crisfield geo 弧长，不再显式拒绝）
    case.opt.czm.iter_method = "arc_length"
    case.opt.czm.geo_nonlinear = true
    ms = JuBat.MechState(czm_mesh)
    r_arc = JuBat.solve_czm_step(
        czm_mesh, ms, case.param, zeros(ndof), case.opt.czm; czm_workspace=czm_workspace,
        eigenstrain=eig)
    @test r_arc isa JuBat.CZMResult
end
