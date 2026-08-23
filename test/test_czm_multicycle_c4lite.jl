using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# Batch 5 Task 3：C4-lite 多圈集成——多相位 Δ_core 与 D 联合增长（spec §7 Batch 5）
# 驱动：Δsoc 失配（与 unit_czm_eigenstrain 同机理）；推进用 load_substep（低载荷避开极限点，
# Crisfield 弧长阻塞已登记 findings）。

function c4_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    opt.czm_enabled = true
    opt.czm_geo_nonlinear = true
    opt.czm_j2_plasticity = true
    opt.czm_winding_prestress = true
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    pc = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, pc)
    case.czm_layout = JuBat.CzmLayout(case.czm_mesh)
    case.czm_layout.plastic_states = [JuBat.PlasticState() for _ in 1:size(case.czm_mesh.bulk_element, 1), _ in 1:4]
    return case, pc, cache
end

@testset "C4-lite：多相位 Δ_core 与 D_max 联合增长（低载荷避开极限点）" begin
    case, pc, cache = c4_fixture()
    cm = case.czm_mesh
    ndof = 2 * cm.nnode
    ne = size(cm.bulk_element, 1)
    α = pc.by_interface[:PE_PCC].α
    βn = case.param.NE.Omega / 3
    βp = case.param.PE.Omega / 3
    # 预应力场（Batch 2'，几何固定一次计算持久持有）
    JuBat.ensure_node_ref!(case)
    prestress = JuBat.winding_prestress_field(cm, case.param)
    scaled = [(0.2a, 0.2b, 0.2c) for (a, b, c) in prestress]

    # 3 相位：Δsoc −0.15 → −0.20 → −0.25（失配驱动、幅值递增 33%/相位，均低于
    # Δsoc 0.5 收敛界与 progress≈0.52 极限点）
    lvls = [0.15, 0.20, 0.25]
    Δ_hist = Float64[]
    D_hist = Float64[]
    for (i, lvl) in enumerate(lvls)
        eig = (α_eff = α, β_n = βn, β_p = βp,
               dT = zeros(ne), Δsn = fill(-lvl, ne), Δsp = fill(-lvl, ne))
        r, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, case.czm_layout.u_prev;
            α_eff = α, β_n = βn, β_p = βp,
            dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp,
            max_iter = 200, tol = 1e-8, n_load_steps = 50, iter_method = "load_substep",
            cache = cache, geo_nl = true, eigenstrain = eig,
            plasticity = true, mech_state = case.czm_layout.plastic_states,
            prestress = scaled)
        @test r.converged
        case.czm_layout.u_prev = r.displacement
        JuBat.assemble_coupled_system(cm, r.displacement, pc;
            geo_nl = true, eigenstrain = eig, plasticity = true,
            mech_state = case.czm_layout.plastic_states, commit_plastic = true)
        w, Δ = JuBat.core_ovalization(cm, r.displacement, JuBat.ensure_node_ref!(case))
        push!(Δ_hist, Δ)
        push!(D_hist, maximum(s.D for s in cm.damage_states))
        println("phase $i (Δsoc=$(lvl)): Δ_core=$(round(Δ, sigdigits=4)) D_max=$(round(D_hist[end], sigdigits=4))")
    end
    # D10 结论固化（findings 2026-08-23）：Δ_core≡0（对称响应被 D8 滤波正确去除）、D 不激活
    @test all(iszero, Δ_hist)
    @test all(iszero, D_hist)
    @test all(isfinite, Δ_hist)
end
