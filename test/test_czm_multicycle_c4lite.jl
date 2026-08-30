using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# Batch 5 Task 3：C4-lite 多圈集成——多相位状态持久化与对称 D10 结论
# 驱动：Δsoc 失配（与 unit_czm_eigenstrain 同机理）；推进用 load_substep（低载荷避开极限点，
# Crisfield 弧长阻塞已登记 findings）。
# 2026-08-30 重构适配：状态在 case.mech（MechState）上收敛提交；
# 界面参数直读 param.PCC/param.NCC（collector_params）。

function c4_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    opt.czm.enabled = true
    opt.czm.geo_nonlinear = true
    opt.czm.j2_plasticity = true
    opt.czm.winding_prestress = true
    opt.czm.fix_inner = false
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    case.mech.plastic_states = [JuBat.PlasticState() for _ in 1:size(case.czm_mesh.bulk_element, 1), _ in 1:4]
    return case
end

@testset "C4-lite：Γ_in,free 不得被位移边界固定" begin
    case = c4_fixture()
    wt = case.czm_mesh.czm_submesh.winding_turn
    nθn = count(==(1), wt[1:count(==(wt[1]), wt)])
    core_dofs = reduce(vcat, ([2n - 1, 2n] for n in 1:nθn))
    bc_dofs, _ = JuBat.extract_bc_dofs(case.czm_mesh, case.param; fix_inner=case.opt.czm.fix_inner)
    @test isempty(intersect(bc_dofs, core_dofs))
end

@testset "C4-lite：自由芯部多相位状态产生不圆度；尚未触发损伤" begin
    case = c4_fixture()
    cm = case.czm_mesh
    ndof = 2 * cm.nnode
    ne = size(cm.bulk_element, 1)
    ms = case.mech
    # 预应力场（Batch 2'，几何固定一次计算持久持有）
    JuBat.ensure_node_ref!(case)
    prestress = JuBat.winding_prestress_field(cm, case.param)
    scaled = [(0.2a, 0.2b, 0.2c) for (a, b, c) in prestress]
    case.opt.czm.iter_method = "load_substep"
    case.opt.czm.max_iter = 200
    case.opt.czm.tol = 1e-8
    case.opt.czm.load_steps = 50

    # 3 相位：Δsoc −0.15 → −0.20 → −0.25（失配驱动、幅值递增 33%/相位，均低于
    # Δsoc 0.5 收敛界与 progress≈0.52 极限点）
    lvls = [0.15, 0.20, 0.25]
    Δ_hist = Float64[]
    D_hist = Float64[]
    initiation_ratio_hist = Float64[]
    for (i, lvl) in enumerate(lvls)
        eig = (dT = zeros(ne), Δsn = fill(-lvl, ne), Δsp = fill(-lvl, ne))
        r = JuBat.solve_czm_step(cm, ms, case.param, zeros(ndof), case.opt.czm;
            dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp,
            eigenstrain = eig, mech_state = ms.plastic_states, prestress = scaled)
        @test r.converged
        # 损伤/位移已由求解器收敛提交到 ms
        JuBat.assemble_coupled_system(cm, r.displacement, case.param;
            damage_states = ms.damage_states,
            geo_nl = true, eigenstrain = eig, plasticity = true,
            mech_state = ms.plastic_states, commit_plastic = true)
        w, Δ = JuBat.core_ovalization(cm, r.displacement, JuBat.ensure_node_ref!(case))
        push!(Δ_hist, Δ)
        push!(D_hist, maximum(s.D for s in ms.damage_states))
        initiation_ratio = maximum(
            max(r.separation_n[j], 0.0) /
            JuBat.collector_params(case.param, cm.cohesive_elements[j].interface_type).δ_0
            for j in eachindex(r.separation_n))
        push!(initiation_ratio_hist, initiation_ratio)
        println("phase $i (Δsoc=$(lvl)): Δ_core=$(round(Δ, sigdigits=4)) " *
                "D_max=$(round(D_hist[end], sigdigits=4)) " *
                "δn/δ0=$(round(initiation_ratio, sigdigits=4))")
    end
    # 自由芯部使螺旋固有不对称可发展为 n≥2 位移；当前低载三相位仍未激活 CZM 损伤。
    @test all(>(0.0), Δ_hist)
    @test all(diff(Δ_hist) .> 0.0)
    @test all(0.0 .< initiation_ratio_hist .< 1.0)
    @test all(iszero, D_hist)
    @test all(isfinite, Δ_hist)
end

@testset "受控分离产生非零损伤并可跨步提交" begin
    case = c4_fixture()
    cm = case.czm_mesh
    ms = case.mech
    separations = [begin
        p = JuBat.collector_params(case.param, coh.interface_type)
        (0.5 * (p.δ_0 + p.δ_c), 0.0)
    end for coh in cm.cohesive_elements]
    trial = JuBat.update_damage_per_interface(cm, JuBat.clone_damage_states(ms.damage_states),
                                               separations, case.param, case.opt.czm.model)
    @test maximum(s.D for s in trial) > 0.0
    ms.damage_states = trial   # 跨步提交
    @test maximum(s.D for s in ms.damage_states) > 0.0
end
