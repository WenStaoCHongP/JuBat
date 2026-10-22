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

# A small real layered Q4/COH mesh with only the selected interface's top
# nodes free. Balanced nodal forces prescribe an exact initial equilibrium;
# no Newton tolerance or correction-budget change is needed for these tests.
function history73_fixture(; order=2)
    case = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), JuBat.Option())
    mesh, _ = JuBat.create_unit_czm_strip(case.param; gsorder=order)
    selected = 1
    free_nodes = mesh.cohesive_elements[selected].nodes_top
    pins = Dict{Int64,Symbol}(n=>:fixed_xy for n in axes(mesh.node,1) if !(n in free_nodes))
    boundary = JuBat.MechanicalBoundary(pins,mesh.node,Dict{Int64,Vector{Symbol}}(),
        Set{Int64}(),Set{Int64}(),Set{Int64}(),Set{Int64}(),0,0,0,0)
    return case.param,mesh,boundary,selected
end

history73_keys(ms) = (copy(ms.u_prev),
    [Tuple(getfield(s,f) for f in fieldnames(JuBat.DamageState)) for s in ms.damage_states],
    ms.gp_damage_states === nothing ? nothing :
        [Tuple(getfield(s,f) for f in fieldnames(JuBat.DamageState)) for s in ms.gp_damage_states],
    deepcopy(ms.plastic_states),deepcopy(ms.winding_prestress),deepcopy(ms.node_ref),ms.mechanics_initialized)

function history73_displacement(param,mesh,selected,factors)
    elem=mesh.cohesive_elements[selected]
    _,normal,_,_=JuBat.cohesive_local_frame(mesh,elem)
    ip=JuBat.collector_params(param,elem.interface_type)
    opening=ip.δ_0*param.scale.δ_czm/param.scale.L
    u=zeros(2size(mesh.node,1))
    for (node,factor) in zip(elem.nodes_top,factors)
        u[2node-1]=factor*opening*normal[1]
        u[2node]=factor*opening*normal[2]
    end
    return u
end

function history73_force(param,mesh,boundary,u,gp;geo=false,beta=1.0)
    zeros_e=zeros(size(mesh.bulk_element,1))
    eig=(dT=zeros_e,Δsn=zeros_e,Δsp=zeros_e)
    weights=JuBat.NCweight(mesh.czm_submesh.mesh_bonded.gs.order)[1]
    _,force,_,_=JuBat.assemble_coupled_system(mesh,u,param;
        czm_workspace=JuBat.CZMWorkspace(mesh,param),
        damage_states=JuBat.aggregate_gp_damage_states(gp,weights),gp_damage_states=gp,
        geo_nl=geo,eigenstrain=geo ? eig : nothing,visc_beta=beta)
    force=copy(force)
    dofs,_=JuBat.mechanical_boundary_dofs(boundary)
    force[dofs].=0.0
    return force
end

function history73_solve(param,mesh,boundary,ms,force,method;geo=false,max_iter=50)
    z=zeros(size(mesh.bulk_element,1))
    opt=JuBat.CzmOptions(iter_method=method,geo_nonlinear=geo,viscous_enabled=false,
        max_iter=max_iter,tol=1e-8,load_steps=2)
    return JuBat.solve_czm_step(mesh,ms,param,force,opt;
        czm_workspace=JuBat.CZMWorkspace(mesh,param),boundary=boundary,
        dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,
        eigenstrain=geo ? (dT=z,Δsn=z,Δsp=z) : nothing)
end

@testset "73 actual solver entries retain nonuniform GP history" begin
    for (method,geo) in (("basic",false),("gp_basic",false),("load_substep",false),("arc_length",true))
        @testset "$method geo=$geo" begin
            p,m,bc,e=history73_fixture()
            ms=JuBat.MechState(m)
            u=history73_displacement(p,m,e,(10.0,50.0))
            gp=JuBat.gp_history_for_step(ms,m)
            force=history73_force(p,m,bc,u,gp;geo=geo)
            ms.u_prev=copy(u)
            result=history73_solve(p,m,bc,ms,force,method;geo=geo)
            @test result.converged
            @test ms.gp_damage_states !== nothing
            if ms.gp_damage_states !== nothing
                ip=JuBat.collector_params(p,m.cohesive_elements[e].interface_type)
                for (g,factor) in enumerate((10.0,50.0))
                    gap=factor*ip.δ_0
                    expected=min(1.0,ip.δ_c*(gap-ip.δ_0)/(gap*(ip.δ_c-ip.δ_0)))
                    @test ms.gp_damage_states[e,g].D ≈ expected
                    @test ms.gp_damage_states[e,g].δ_max_n ≈ gap
                end
                @test ms.gp_damage_states[e,1].D != ms.gp_damage_states[e,2].D
                @test ms.damage_states[e].D ≈ sum(s.D for s in ms.gp_damage_states[e,:])/2

                # Export/statistics and frozen output assembly must not advance
                # the accepted GP history, even when repeated.
                before=history73_keys(ms)
                vars=Dict{String,Union{Array{Float64},Float64}}()
                JuBat.czm_output_to_variables(result,vars;committed_states=ms.damage_states)
                JuBat.czm_output_to_variables(result,vars;committed_states=ms.damage_states)
                JuBat.get_damage_statistics(ms.damage_states)
                JuBat.assemble_czm_system(m,ms.u_prev,p;damage_states=ms.damage_states,
                    gp_damage_states=ms.gp_damage_states,visc_beta=0.0)
                @test history73_keys(ms)==before

                # Resume from the true accepted state; only the target load
                # changes. No mean-to-GP seed and no replacement displacement.
                peak=deepcopy(ms)
                for factor in (0.9,1.1)
                    target=factor*u
                    force=history73_force(p,m,bc,target,ms.gp_damage_states;geo=geo)
                    resumed=history73_solve(p,m,bc,ms,force,method;geo=geo)
                    @test resumed.converged
                    @test ms.u_prev ≈ target rtol=1e-5
                    @test all(s.D >= old.D for (s,old) in zip(ms.gp_damage_states,peak.gp_damage_states))
                    factor<1 && @test [s.D for s in ms.gp_damage_states]==[s.D for s in peak.gp_damage_states]
                end
                # A failed attempt leaves every external history unchanged,
                # and the same accepted state can be retried successfully.
                impossible=copy(force)
                free=first(setdiff(1:length(force),JuBat.mechanical_boundary_dofs(bc)[1]))
                impossible[free]+=1e6
                before=history73_keys(ms)
                failed=try
                    history73_solve(p,m,bc,ms,impossible,method;geo=geo,max_iter=geo ? 1 : 0)
                catch err
                    @test err isa ErrorException
                    nothing
                end
                @test failed===nothing || !failed.converged
                @test history73_keys(ms)==before
                @test history73_solve(p,m,bc,ms,force,method;geo=geo).converged

                bad=deepcopy(ms);bad.gp_damage_states=nothing;before=history73_keys(bad)
                @test_throws ArgumentError history73_solve(p,m,bc,bad,force,method;geo=geo)
                @test history73_keys(bad)==before
                bad=JuBat.MechState(m);bad.damage_states[e].D=0.2;before=history73_keys(bad)
                @test_throws ArgumentError history73_solve(p,m,bc,bad,force,method;geo=geo)
                @test history73_keys(bad)==before
            end
        end
    end
end

@testset "73 signed mean cancellation cannot erase GP opening" begin
    for (method,geo) in (("basic",false),("gp_basic",false),("load_substep",false),("arc_length",true))
        p,m,bc,e=history73_fixture();ms=JuBat.MechState(m)
        u=history73_displacement(p,m,e,(-10.0,10.0))
        force=history73_force(p,m,bc,u,JuBat.gp_history_for_step(ms,m);geo=geo)
        ms.u_prev=copy(u)
        r=history73_solve(p,m,bc,ms,force,method;geo=geo)
        @test r.converged
        @test abs(r.separation_n[e])<1e-12
        @test ms.gp_damage_states[e,1].D==0.0
        @test ms.gp_damage_states[e,2].D>0.0
        @test r.damage[e]>0.0
    end
end

@testset "73 exported basic API defaults to real GP history" begin
    p,m,bc,_=history73_fixture()
    a=JuBat.MechState(m);b=deepcopy(a);z=zeros(size(m.bulk_element,1));force=zeros(length(a.u_prev))
    kwargs=(boundary=bc,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z)
    r1=JuBat.solve_czm_basic_step(m,force,p,a;czm_workspace=JuBat.CZMWorkspace(m,p),kwargs...)
    r2=JuBat.solve_czm_basic_step(m,force,p,b;czm_workspace=JuBat.CZMWorkspace(m,p),gp_history=true,kwargs...)
    @test r1.converged && r2.converged
    @test history73_keys(a)==history73_keys(b)
    @test r1.displacement==r2.displacement && r1.damage==r2.damage && r1.residual_norm==r2.residual_norm
    before=history73_keys(a)
    @test_throws ArgumentError JuBat.solve_czm_basic_step(m,force,p,a;
        czm_workspace=JuBat.CZMWorkspace(m,p),gp_history=false,kwargs...)
    @test history73_keys(a)==before
end

@testset "73 cycle API preserves real weighted GP history" begin
    p,m,bc,e=history73_fixture(order=3)
    ms=JuBat.MechState(m)
    ms.gp_damage_states=JuBat.gp_history_for_step(ms,m)
    for g in 1:3
        s=ms.gp_damage_states[1,g]
        s.D=0.1g;s.D_visc=0.05g;s.accumulated_damage=0.1g
    end
    ms.damage_states=JuBat.aggregate_gp_damage_states(ms.gp_damage_states,JuBat.NCweight(3)[1])
    ms.mechanics_initialized=true
    before=history73_keys(ms);old_gp=ms.gp_damage_states;old_elements=ms.damage_states
    @test JuBat.accumulate_cycle_damage!(ms,0.0) === ms
    @test history73_keys(ms)==before
    @test ms.gp_damage_states === old_gp && ms.damage_states === old_elements
    JuBat.accumulate_cycle_damage!(ms,0.85)
    @test ms.gp_damage_states !== nothing
    if ms.gp_damage_states !== nothing
        @test [s.D for s in ms.gp_damage_states[1,:]]==[0.1,1.0,1.0]
        @test [s.D_visc for s in ms.gp_damage_states[1,:]] ≈ [0.05,0.1,0.15]
        @test !any(s.fractured for s in ms.gp_damage_states)
        @test ms.damage_states[1].D ≈ (0.1+4+1)/6
        @test ms.damage_states[1].accumulated_damage ≈ (0.95+4*1.05+1.15)/6
        # Physical-time gp_basic can continue this lagged helper state;
        # the helper itself did not invent an elapsed time or relax D_visc.
        z=zeros(size(m.bulk_element,1))
        old=deepcopy(ms.gp_damage_states)
        r=JuBat.solve_czm_step(m,ms,p,zeros(length(ms.u_prev)),
            JuBat.CzmOptions(iter_method="gp_basic",viscous_enabled=true,viscous_tau=1.0);
            czm_workspace=JuBat.CZMWorkspace(m,p),dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,dt_seconds=1.0)
        @test r.converged
        @test all(s.D_visc ≈ o.D_visc+0.5*(o.D-o.D_visc) for (s,o) in zip(ms.gp_damage_states,old))
    end
    JuBat.reset_damage_states!(ms)
    @test ms.gp_damage_states !== nothing
    @test ms.mechanics_initialized
    if ms.gp_damage_states !== nothing
        @test all(s.D==s.D_visc==s.accumulated_damage==0.0 for s in ms.gp_damage_states)
    end
end

@testset "73 auxiliary invalid input and saturation publish atomically" begin
    p,m,bc,_=history73_fixture(order=3)
    seed=JuBat.MechState(m)
    seed.gp_damage_states=JuBat.gp_history_for_step(seed,m)
    seed.mechanics_initialized=true
    for increment in (-0.1,NaN,Inf)
        ms=deepcopy(seed);before=history73_keys(ms)
        @test_throws ArgumentError JuBat.accumulate_cycle_damage!(ms,increment)
        @test isequal(history73_keys(ms),before)
    end
    for corrupt in (s->(s.D_visc=0.2),s->(s.fractured=true),
                    s->(s.δ_max_eff=Inf),s->(s.accumulated_damage=-0.1),
                    s->(s.D=1.0;s.D_visc=1-5e-11))
        ms=deepcopy(seed);corrupt(ms.gp_damage_states[end]);before=history73_keys(ms)
        @test_throws ArgumentError JuBat.accumulate_cycle_damage!(ms,0.1)
        @test isequal(history73_keys(ms),before)
    end
    for initialized in (false,true)
        ms=JuBat.MechState(m);ms.mechanics_initialized=initialized;before=history73_keys(ms)
        @test JuBat.accumulate_cycle_damage!(ms,0.0) === ms
        @test_throws ArgumentError JuBat.accumulate_cycle_damage!(ms,0.1)
        @test history73_keys(ms)==before
        if initialized
            @test_throws ArgumentError JuBat.reset_damage_states!(ms)
            @test history73_keys(ms)==before
        else
            @test JuBat.reset_damage_states!(ms) === ms
            @test ms.gp_damage_states === nothing
        end
    end
    # New saturation follows 72's canonical threshold; old noncanonical
    # histories above were rejected rather than silently migrated.
    ms=deepcopy(seed);s=ms.gp_damage_states[1,1]
    s.D=1-1e-11;s.D_visc=1-5e-11;s.accumulated_damage=0.9
    JuBat.accumulate_cycle_damage!(ms,0.2)
    @test ms.gp_damage_states[1,1].D==ms.gp_damage_states[1,1].D_visc==1.0
    @test ms.gp_damage_states[1,1].fractured
    before=deepcopy(ms.gp_damage_states[1,1])
    JuBat.accumulate_cycle_damage!(ms,0.1)
    @test ms.gp_damage_states[1,1].accumulated_damage==before.accumulated_damage
    ms=deepcopy(seed);s=ms.gp_damage_states[1,1]
    s.D=1-1e-11;s.D_visc=1-5e-11;s.accumulated_damage=s.D
    oldD=s.D;oldDv=s.D_visc
    JuBat.accumulate_cycle_damage!(ms,1e-13)
    @test ms.gp_damage_states[1,1].D==oldD
    @test ms.gp_damage_states[1,1].D_visc==oldDv
    @test !ms.gp_damage_states[1,1].fractured
    # Both a GP sum overflow and weighted aggregate overflow must stop before publication.
    for amount in (floatmax(Float64),floatmax(Float64)/2)
        ms=deepcopy(seed);foreach(s->(s.accumulated_damage=amount),ms.gp_damage_states)
        before=history73_keys(ms)
        @test_throws ArgumentError JuBat.accumulate_cycle_damage!(ms,amount)
        @test isequal(history73_keys(ms),before)
    end
    ms=deepcopy(seed);ms.u_prev=history73_displacement(p,m,1,(0.2,0.3));ms.node_ref=copy(m.node)
    old_u=ms.u_prev;old_ref=ms.node_ref;shape=size(ms.gp_damage_states)
    JuBat.reset_damage_states!(ms)
    @test size(ms.gp_damage_states)==shape
    @test ms.u_prev===old_u && ms.node_ref===old_ref && ms.mechanics_initialized
    force=history73_force(p,m,bc,ms.u_prev,ms.gp_damage_states)
    @test history73_solve(p,m,bc,ms,force,"basic").converged
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
