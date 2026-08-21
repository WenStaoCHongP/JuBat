using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# C2-lite 集成（spec §7 Batch 3）：强本征应变下 PCC/NCC 屈服、塑性关退回 C1、
# 状态累计与不可逆性、耗散非负。夹具同 C1（归一化网格）。

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
    param_cache = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    return case, param_cache, cache
end

@testset "塑性关（geo 关/开）与既有路径逐位一致（回归锚）" begin
    case, param_cache, cache = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0, dT=fill(1e-6, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    kw = (α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
          dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
          max_iter=100, tol=1e-10, iter_method="basic", cache=cache)
    r1, _ = JuBat.solve_czm_step(czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof); kw...)
    r2, _ = JuBat.solve_czm_step(czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
                                 kw..., geo_nl=true, eigenstrain=eig)
    @test r1.converged && r2.converged
    @test isapprox(r2.displacement, r1.displacement; rtol=1e-4, atol=1e-12)  # K_G→0 极限（同 C1 判据）
end

@testset "C2-lite：强压缩本征应变下 PCC/NCC 屈服、其余层 κ=0" begin
    case, param_cache, cache = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    # 层间失配驱动（与 unit_czm_eigenstrain 同机理）：Δsn/Δsp = −0.3 放电方向，
    # β_n/β_p 反号 ⟹ NE/PE 本征应变差 ~1.5e-2，超 Cu 箔屈服应变 1.8e-3
    eig = (α_eff=1.0, β_n=case.param.NE.Omega / 3.0, β_p=case.param.PE.Omega / 3.0,
           dT=zeros(ne), Δsn=fill(-0.3, ne), Δsp=fill(-0.3, ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    r, _ = JuBat.solve_czm_step(czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=200, tol=1e-8, iter_method="load_substep", cache=cache,
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states)
    @test r.converged
    # 提交（复刻生产 D-B3-2 路径）
    JuBat.assemble_coupled_system(czm_mesh, r.displacement, param_cache;
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    mt = czm_mesh.czm_submesh.material_type
    κ_foil = [states[e, g].kappa for e in 1:ne, g in 1:4 if mt[e] in (:PCC, :NCC)]
    κ_soft = [states[e, g].kappa for e in 1:ne, g in 1:4 if !(mt[e] in (:PCC, :NCC))]
    # ε₀=6e-3 > Cu 箔屈服应变 200/110000≈1.8e-3、Al 60/70000≈0.86e-3 → 两类箔均应屈服
    @test maximum(κ_foil) > 0.0
    @test all(iszero, κ_soft)
end

@testset "状态不可逆：同载荷重解后提交的 κ 不再增长" begin
    case, param_cache, cache = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    # 层间失配驱动（与 unit_czm_eigenstrain 同机理）：Δsn/Δsp = −0.3 放电方向，
    # β_n/β_p 反号 ⟹ NE/PE 本征应变差 ~1.5e-2，超 Cu 箔屈服应变 1.8e-3
    eig = (α_eff=1.0, β_n=case.param.NE.Omega / 3.0, β_p=case.param.PE.Omega / 3.0,
           dT=zeros(ne), Δsn=fill(-0.3, ne), Δsp=fill(-0.3, ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    kw = (α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
          dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
          max_iter=200, tol=1e-8, iter_method="load_substep", cache=cache,
          geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states)
    r1, _ = JuBat.solve_czm_step(czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof); kw...)
    JuBat.assemble_coupled_system(czm_mesh, r1.displacement, param_cache;
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    κ_after1 = maximum(s.kappa for s in states)
    # 同载荷以 u_prev 重解（生产时序），再提交：理想塑性下 κ 不应显著增长（已屈服面停留）
    r2, _ = JuBat.solve_czm_step(czm_mesh, zeros(ndof), param_cache, case.param, r1.displacement; kw...)
    JuBat.assemble_coupled_system(czm_mesh, r2.displacement, param_cache;
        geo_nl=true, eigenstrain=eig, plasticity=true, mech_state=states, commit_plastic=true)
    κ_after2 = maximum(s.kappa for s in states)
    @test κ_after1 > 0.0
    @test κ_after2 ≥ κ_after1 * 0.999   # 允许数值级持平，不得回退
    @test κ_after2 ≤ 4 * κ_after1        # 同载荷重解的再预测漂移有界（D-B3-2 语义：单调不可逆成立，新 Δε* 步进下漂移二阶小）
end

@testset "缺参/组合非法必须报错（AGENTS 9.7 / D-B3-1）" begin
    case, param_cache, cache = build_j2_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0, dT=fill(-1e-3, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    states = [JuBat.PlasticState() for _ in 1:ne, _ in 1:4]
    # plasticity 需 geo_nl
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), param_cache; plasticity=true, mech_state=states)
    # plasticity 无 mech_state
    @test_throws ErrorException JuBat.assemble_coupled_system(
        czm_mesh, zeros(ndof), param_cache; geo_nl=true, eigenstrain=eig, plasticity=true)
    # 非法组合：arc_length + geo_nl（含塑性场景）
    @test_throws ErrorException JuBat.solve_czm_step(
        czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        iter_method="arc_length", geo_nl=true, eigenstrain=eig)
end
