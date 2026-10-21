using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function material_tangent_error(δn, δt, state, ip; model="mix")
    h = 1e-6 * ip.δ_0
    J = JuBat.bilinear_tangent(δn, δt, state, ip, model)
    J_fd = zeros(2, 2)
    for j in 1:2
        plus = j == 1 ? (δn + h, δt) : (δn, δt + h)
        minus = j == 1 ? (δn - h, δt) : (δn, δt - h)
        T_plus = JuBat.bilinear_traction_state(plus..., state, ip, model)
        T_minus = JuBat.bilinear_traction_state(minus..., state, ip, model)
        J_fd[:, j] .= ([T_plus[1], T_plus[2]] .- [T_minus[1], T_minus[2]]) ./ (2h)
    end
    return norm(J - J_fd) / norm(J_fd)
end

@testset "CZM mix tangent follows the active smooth branch" begin
    case = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), JuBat.Option())
    ip = deepcopy(case.param.PCC)
    fresh = JuBat.DamageState()

    unequal = deepcopy(ip)
    unequal.δ_0_t = 2ip.δ_0
    unequal.δ_c_t = 3ip.δ_c
    unequal.τ_max = unequal.K_t * unequal.δ_0_t
    unequal.G_c_t = unequal.τ_max * unequal.δ_c_t / 2
    @test material_tangent_error(2ip.δ_0, 3ip.δ_0, fresh, unequal) < 1e-4

    @test material_tangent_error(-ip.δ_0, 3ip.δ_0, fresh, ip) < 1e-4
    @test material_tangent_error(-ip.δ_0, -3ip.δ_0, fresh, ip) < 1e-4

    _, _, _, normal_peak = JuBat.bilinear_traction_state(
        1.5ip.δ_0, 0.0, fresh, unequal, "mix")
    _, _, _, changed_mode = JuBat.bilinear_traction_state(
        0.0, 1.7ip.δ_0, normal_peak, unequal, "mix")
    @test changed_mode.D >= normal_peak.D
end

function arc_gp_fixture()
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param;
        nθ=8, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.mesh["czm"] = JuBat.create_czm_mesh(md.czm_submesh,
        case.mesh["thermal2D"], case.param)
    return case, JuBat.MechState(case.mesh["czm"])
end

@testset "Arc GP trial keeps each integration point's peak through unloading" begin
    case, ms = arc_gp_fixture()
    cm = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    bc_dofs, bc_vals = JuBat.extract_bc_dofs(cm, case.param)
    selected = findfirst(elem -> all(node ->
        !(2node - 1 in bc_dofs) && !(2node in bc_dofs), elem.nodes_top),
        cm.cohesive_elements)
    @test selected !== nothing
    elem = cm.cohesive_elements[selected]
    _, normal, _, _ = JuBat.cohesive_local_frame(cm, elem)
    ip = JuBat.collector_params(case.param, elem.interface_type)
    opening = ip.δ_0 * case.param.scale.δ_czm / case.param.scale.L
    u_peak = zeros(2 * size(cm.node, 1))
    for (node, factor) in zip(elem.nodes_top, (10.0, 50.0))
        u_peak[2node - 1] = factor * opening * normal[1]
        u_peak[2node] = factor * opening * normal[2]
    end
    gp_start = JuBat.gp_history_for_step(ms, cm)
    zero_force = zeros(2 * size(cm.node, 1))
    peak = JuBat.arc_equilibrium_trial(cm, u_peak, case.param,
        gp_start, zero_force, bc_dofs, bc_vals; czm_workspace=czm_workspace, czm_model="mix")
    unloaded = JuBat.arc_equilibrium_trial(cm, 0.9u_peak, case.param,
        peak.gp_damage_states, zero_force, bc_dofs, bc_vals; czm_workspace=czm_workspace,
        czm_model="mix")
    peak_D = [s.D for s in peak.gp_damage_states[selected, :]]
    @test maximum(peak_D) > 0.9
    @test maximum(peak_D) - minimum(peak_D) > 1e-3
    @test all(unloaded.gp_damage_states[selected, g].D == peak_D[g]
        for g in eachindex(peak_D))
    @test all(unloaded.gp_damage_states[selected, g].D_visc ==
              peak.gp_damage_states[selected, g].D_visc
        for g in eachindex(peak_D))
    @test all(s.D == 0.0 for s in gp_start)

    h = 1e-6
    plus = JuBat.arc_equilibrium_trial(cm, u_peak + h * u_peak,
        case.param, gp_start, zero_force, bc_dofs, bc_vals; czm_workspace=czm_workspace,
        czm_model="mix")
    minus = JuBat.arc_equilibrium_trial(cm, u_peak - h * u_peak,
        case.param, gp_start, zero_force, bc_dofs, bc_vals; czm_workspace=czm_workspace,
        czm_model="mix")
    free = setdiff(1:length(u_peak), bc_dofs)
    finite_difference = -((plus.R - minus.R) ./ (2h))[free]
    assembled = (peak.K_total * u_peak)[free]
    @test norm(assembled - finite_difference) / norm(finite_difference) < 1e-4
end

@testset "Arc step owns GP history and rejects an unrecoverable average" begin
    case, ms = arc_gp_fixture()
    cm = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
    case.opt.czm.iter_method = "arc_length"
    case.opt.czm.geo_nonlinear = false
    case.opt.czm.tol = 1e-10
    ne = size(cm.bulk_element, 1)
    zero_elem = zeros(ne)
    result = JuBat.solve_czm_step(cm, ms, case.param,
        zeros(2 * size(cm.node, 1)), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=zero_elem, Δsoc_n_elem=zero_elem, Δsoc_p_elem=zero_elem)
    @test result.converged
    @test ms.gp_damage_states !== nothing
    @test size(ms.gp_damage_states, 1) == length(cm.cohesive_elements)

    invalid_ms = JuBat.MechState(cm)
    invalid_ms.damage_states[1].D = 0.5
    @test_throws ArgumentError JuBat.solve_czm_step(cm, invalid_ms,
        case.param, zeros(2 * size(cm.node, 1)), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=zero_elem, Δsoc_n_elem=zero_elem, Δsoc_p_elem=zero_elem)

    lagged_ms = JuBat.MechState(cm)
    lagged_ms.gp_damage_states = JuBat.gp_history_for_step(lagged_ms, cm)
    lagged_ms.gp_damage_states[1, 1].D = 0.5
    lagged_ms.gp_damage_states[1, 1].D_visc = 0.25
    @test_throws ArgumentError JuBat.solve_czm_step(cm, lagged_ms,
        case.param, zeros(2 * size(cm.node, 1)), case.opt.czm; czm_workspace=czm_workspace,
        dT_elem=zero_elem, Δsoc_n_elem=zero_elem, Δsoc_p_elem=zero_elem)
end

@testset "Failed arc step saves complete mechanical input for replay" begin
    case, ms = arc_gp_fixture()
    @test hasproperty(case.opt.czm, :snapshot_path)
    if hasproperty(case.opt.czm, :snapshot_path)
        cm = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
        outdir = joinpath(@__DIR__, "..", "output", "test_czm_arc_gp_history")
        mkpath(outdir)
        snapshot_path = joinpath(outdir,
            "failed_step_$(getpid())_$(time_ns()).bin")
        case.opt.czm.snapshot_path = snapshot_path
        case.opt.czm.iter_method = "arc_length"
        case.opt.czm.geo_nonlinear = false
        case.opt.czm.max_iter = 0
        case.opt.czm.tol = 1e-30
        F_ext = zeros(2 * size(cm.node, 1))
        bc_dofs, _ = JuBat.extract_bc_dofs(cm, case.param)
        free_dof = first(setdiff(1:length(F_ext), bc_dofs))
        F_ext[free_dof] = 1.0
        zero_elem = zeros(size(cm.bulk_element, 1))
        result = JuBat.solve_czm_step(cm, ms, case.param, F_ext,
            case.opt.czm; czm_workspace=czm_workspace, dT_elem=zero_elem,
            Δsoc_n_elem=zero_elem, Δsoc_p_elem=zero_elem)
        @test !result.converged
        @test isfile(snapshot_path)
        @test all(iszero, ms.u_prev)
        @test ms.gp_damage_states === nothing
        replay_result, replay_ms = JuBat.replay_czm_step_snapshot(snapshot_path)
        @test replay_result.converged == result.converged
        @test replay_result.residual_norm == result.residual_norm
        @test all(iszero, replay_ms.u_prev)
    end
end

@testset "Arc wall-time limit exits before unbounded correction" begin
    case, ms = arc_gp_fixture()
    @test hasproperty(case.opt.czm, :wall_time_limit_seconds)
    if hasproperty(case.opt.czm, :wall_time_limit_seconds)
        cm = case.mesh["czm"]
    czm_workspace = JuBat.CZMWorkspace(case.mesh["czm"], case.param)
        case.opt.czm.iter_method = "arc_length"
        case.opt.czm.geo_nonlinear = false
        case.opt.czm.wall_time_limit_seconds = 1e-9
        F_ext = zeros(2 * size(cm.node, 1))
        bc_dofs, _ = JuBat.extract_bc_dofs(cm, case.param)
        F_ext[first(setdiff(1:length(F_ext), bc_dofs))] = 1.0
        zero_elem = zeros(size(cm.bulk_element, 1))
        result = JuBat.solve_czm_step(cm, ms, case.param, F_ext,
            case.opt.czm; czm_workspace=czm_workspace, dT_elem=zero_elem,
            Δsoc_n_elem=zero_elem, Δsoc_p_elem=zero_elem)
        @test !result.converged
        @test result.iterations == 0
        @test all(iszero, ms.u_prev)
        @test ms.gp_damage_states === nothing
    end
end
