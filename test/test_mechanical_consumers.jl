using Test
include(joinpath(@__DIR__,"..","tools","validate_case_batch_c_trace.jl"))
const OUT=normpath(joinpath(@__DIR__,"..","output","test_mechanical_consumers"))
mkpath(OUT)
@testset "Physical consumers reject invalid values without committing" begin
    c=fixture_c("gp");m=c.mesh["czm"];w=JuBat.CZMWorkspace(m,c.param);z=zeros(size(m.bulk_element,1));F=zeros(2size(m.node,1))
    for (corrupt,expected) in ((s->pop!(s.u_prev),DimensionMismatch),(s->(s.u_prev[1]=NaN),ArgumentError),(s->pop!(s.damage_states),DimensionMismatch),(s->(s.damage_states[1].D_visc=1.1),ArgumentError),(s->(s.damage_states[1].δ_max_eff=NaN),ArgumentError))
        ms=JuBat.MechState(m);corrupt(ms);before=mech_numbers(ms)
        @test_throws expected JuBat.solve_czm_step(m,ms,c.param,F,c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,dt_seconds=0.02)
        @test isequal(mech_numbers(ms),before)
    end
    ms=JuBat.MechState(m);r=JuBat.solve_czm_step(m,ms,c.param,F,c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,dt_seconds=0.02)
    @test r.converged && ms.mechanics_initialized
    for (corrupt,expected) in ((s->(s.gp_damage_states=nothing),ArgumentError),(s->(s.gp_damage_states=s.gp_damage_states[:,1:1]),DimensionMismatch),(s->(s.gp_damage_states[1].D_visc=-0.1),ArgumentError),(s->(s.gp_damage_states[1].δ_max_n=Inf),ArgumentError))
        bad=deepcopy(ms);corrupt(bad);before=mech_numbers(bad)
        @test_throws expected JuBat.gp_history_for_step(bad,m)
        @test isequal(mech_numbers(bad),before)
    end
    @test_throws DimensionMismatch JuBat.recover_bulk_stress(m.node,m.bulk_element,m.czm_submesh.material_type,F[1:end-1],z,c.param)
    badu=copy(F);badu[1]=Inf
    @test_throws ArgumentError JuBat.recover_bulk_stress(m.node,m.bulk_element,m.czm_submesh.material_type,badu,z,c.param)
    states=[JuBat.PlasticState() for _ in 1:length(z),_ in 1:4]
    bad=deepcopy(ms);bad.plastic_states=states;bad.mechanics_initialized=true
    bad.plastic_states[1].eps_p=(NaN,0.0,0.0)
    @test_throws ArgumentError JuBat.prepare_plastic_state_buffers(bad,m,true)
    bad.plastic_states[1].eps_p=(0.0,0.0,0.0)
    @test_throws ErrorException JuBat.project_bulk_stress(m,bad.plastic_states)
    foreach(s->(s.sigma=(0.0,0.0,0.0)),bad.plastic_states)
    bad.plastic_states[1].kappa=Inf
    @test_throws ErrorException JuBat.project_bulk_stress(m,bad.plastic_states)
    ms.node_ref=copy(m.node);ms.node_ref[1]=Inf
    @test_throws ArgumentError JuBat.ensure_node_ref!(c,ms)
    ms.node_ref=zeros(1,2)
    @test_throws DimensionMismatch JuBat.ensure_node_ref!(c,ms)
    fresh=JuBat.MechState(m);oldnode=m.node[1];m.node[1]=NaN
    try
        @test_throws ArgumentError JuBat.ensure_node_ref!(c,fresh)
        @test fresh.node_ref===nothing
    finally
        m.node[1]=oldnode
    end
    # The actual GL kernel consumes S0; no complete-state scan is added to Newton iterations.
    @test_throws ArgumentError JuBat.gl_element_residual_tangent([0.0,1.0,1.0,0.0],[0.0,0.0,1.0,1.0],zeros(8),Matrix{Float64}(I,3,3),0.0,2;sigma0=(NaN,0.0,0.0))
end
@testset "Initial consumers fail before replacing an existing session" begin
    for kind in ("gp","geo")
        c=fixture_c(kind);a=JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.03,nothing;dt_range=[0.01,0.01]).final_state
        bad=deepcopy(a)
        if kind=="gp";bad.mech.u_prev[1]=NaN
        else;bad.mech.plastic_states=nothing;end
        old=(y=[-99.0],mech=nothing,t_global=-1.0);h=Dict("sentinel"=>[42.0]);cache=Ref(7)
        rt=Dict{String,Any}("accepted"=>old,"history"=>h,"czm_workspace"=>cache);shots=JuBat.CZMSnapshot[]
        expected=kind=="gp" ? ArgumentError : ErrorException
        @test_throws expected JuBat.Solve(c;initial_state=bad,variables=rt,czm_snapshots=shots)
        @test rt["accepted"]===old && rt["history"]===h && rt["czm_workspace"]===cache
        @test h==Dict("sentinel"=>[42.0]) && isempty(shots)
    end
    c=fixture_c("export");y=JuBat.ModelInitialisation_MultiSPMe(c);y[1]=NaN
    old=(y=[-99.0],mech=nothing,t_global=-1.0);h=Any[:sentinel];rt=Dict{String,Any}("accepted"=>old,"history"=>h)
    @test_throws ArgumentError JuBat.solve_phase(c,JuBat.PHASE_REST,0.01,(y=y,mech=nothing,t_global=0.0);dt_range=[0.01,0.01],variables=rt)
    @test rt["accepted"]===old && rt["history"]===h && h==Any[:sentinel]
end
@testset "Unused histories are deferred until their real consumer" begin
    c=fixture_c("gp");a=JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.03,nothing;dt_range=[0.01,0.01]).final_state
    a.mech.gp_damage_states=nothing;a.mech.node_ref=zeros(1,2)
    short=JuBat.solve_phase(c,JuBat.PHASE_REST,0.01,a;dt_range=[0.01,0.01])
    @test short.final_state.mech.gp_damage_states===nothing
    @test size(short.final_state.mech.node_ref)==(1,2)
    rt=Dict{String,Any}();shots=JuBat.CZMSnapshot[];prior=Ref{Any}();history=Ref{Any}();history_values=Ref{Any}();history_aliases=Ref{Any}();numbers=Ref{Any}()
    observe=ev->if ev.kind==:accepted_commit
        prior[]=rt["accepted"];numbers[]=mech_numbers(prior[].mech);history[]=rt["history"];history_values[]=deepcopy(history[]);history_aliases[]=copy(history[])
    end
    @test_throws ArgumentError JuBat.solve_phase(c,JuBat.PHASE_REST,0.02,a;dt_range=[0.01,0.01],variables=rt,czm_snapshots=shots,chemical_step_observer=observe)
    @test isassigned(prior) && rt["accepted"]===prior[] && isequal(mech_numbers(rt["accepted"].mech),numbers[])
    @test rt["history"]===history[] && isequal(rt["history"],history_values[])
    @test all(rt["history"][k]===history_aliases[][k] for k in keys(history_aliases[])) && isempty(shots)
    @test_throws DimensionMismatch JuBat.ensure_node_ref!(c,short.final_state.mech)
end
