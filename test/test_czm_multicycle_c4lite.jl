using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# Batch 5 Task 3：C4-lite 多圈集成——多相位状态持久化与对称 D10 结论
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
    opt.czm_fix_inner = false
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    pc = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, pc; fix_inner=case.opt.czm_fix_inner)
    case.czm_layout = JuBat.CzmLayout(case.czm_mesh)
    case.czm_layout.plastic_states = [JuBat.PlasticState() for _ in 1:size(case.czm_mesh.bulk_element, 1), _ in 1:4]
    return case, pc, cache
end

@testset "C4-lite：Γ_in,free 不得被位移边界固定" begin
    case, _, cache = c4_fixture()
    wt = case.czm_mesh.czm_submesh.winding_turn
    nθn = count(==(1), wt[1:count(==(wt[1]), wt)])
    core_dofs = reduce(vcat, ([2n - 1, 2n] for n in 1:nθn))
    @test isempty(intersect(cache.bc_dofs, core_dofs))
end

@testset "C4-lite：自由芯部多相位状态产生不圆度；尚未触发损伤" begin
    case, pc, cache = c4_fixture()
    cm = case.czm_mesh
    ndof = 2 * cm.nnode
    ne = size(cm.bulk_element, 1)
    # 预应力场（Batch 2'，几何固定一次计算持久持有）
    JuBat.ensure_node_ref!(case)
    prestress = JuBat.winding_prestress_field(cm, case.param)
    scaled = [(0.2a, 0.2b, 0.2c) for (a, b, c) in prestress]

    # 3 相位：Δsoc −0.15 → −0.20 → −0.25（失配驱动、幅值递增 33%/相位，均低于
    # Δsoc 0.5 收敛界与 progress≈0.52 极限点）
    lvls = [0.15, 0.20, 0.25]
    Δ_hist = Float64[]
    D_hist = Float64[]
    initiation_ratio_hist = Float64[]
    for (i, lvl) in enumerate(lvls)
        eig = (dT = zeros(ne), Δsn = fill(-lvl, ne), Δsp = fill(-lvl, ne))
        r, updated = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, case.czm_layout.u_prev;
            dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp,
            max_iter = 200, tol = 1e-8, n_load_steps = 50, iter_method = "load_substep",
            cache = cache, geo_nl = true, eigenstrain = eig,
            plasticity = true, mech_state = case.czm_layout.plastic_states,
            prestress = scaled)
        @test r.converged
        case.czm_layout.u_prev = r.displacement
        cm.damage_states = updated.damage_states
        JuBat.assemble_coupled_system(cm, r.displacement, pc;
            geo_nl = true, eigenstrain = eig, plasticity = true,
            mech_state = case.czm_layout.plastic_states, commit_plastic = true)
        w, Δ = JuBat.core_ovalization(cm, r.displacement, JuBat.ensure_node_ref!(case))
        push!(Δ_hist, Δ)
        push!(D_hist, maximum(s.D for s in cm.damage_states))
        initiation_ratio = maximum(
            max(r.separation_n[j], 0.0) /
            pc.by_interface[cm.cohesive_elements[j].interface_type].δ_0_n
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
    case, pc, _ = c4_fixture()
    cm = case.czm_mesh
    separations = [begin
        p = pc.by_interface[coh.interface_type]
        (0.5 * (p.δ_0_n + p.δ_c_n), 0.0)
    end for coh in cm.cohesive_elements]
    trial = JuBat.update_damage_per_interface(cm, JuBat.clone_damage_states(cm.damage_states),
                                               separations, pc)
    updated = JuBat.clone_czm_mesh_with_damage(cm, trial)
    @test maximum(s.D for s in updated.damage_states) > 0.0
    @test updated.czm_submesh === cm.czm_submesh
    @test updated.cohesive_to_thermal === cm.cohesive_to_thermal
    cm.damage_states = updated.damage_states
    @test maximum(s.D for s in cm.damage_states) > 0.0
end
