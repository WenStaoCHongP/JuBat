using Test, Serialization
include(joinpath(@__DIR__,"..","tools","validate_case_batch_c_trace.jl"))
const OUT_C=normpath(joinpath(@__DIR__,"..","output","test_case_batch_c"))
mkpath(OUT_C)
state_numbers(a)=(;y=copy(a.y),mech=mech_numbers(a.mech),t=a.t_global)

if !("--entry-only" in ARGS)
@testset "C real ownership and restoration contracts" begin
    @test fieldnames(JuBat.Case)==(:param_dim,:param,:opt,:mesh,:index)
    @test !(:contact in fieldnames(JuBat.MechState))
    for kind in ("gp","geo")
        c=fixture_c(kind); rt=Dict{String,Any}(); shots=JuBat.CZMSnapshot[]
        mechanical_dt=Float64[]
        JuBat.mechanical_step_observer[]=ev->(ev.phase==:before && push!(mechanical_dt,ev.dt_seconds))
        r=try
            JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.03,nothing;dt_range=[0.01,0.01],variables=rt,czm_snapshots=shots)
        finally
            JuBat.mechanical_step_observer[]=nothing
        end
        @test mechanical_dt==[0.02]
        a=r.final_state
        @test state_numbers(a)==state_numbers(rt["accepted"])
        @test a !== rt["accepted"] && a.y !== rt["accepted"].y && a.mech !== rt["accepted"].mech
        @test length(shots)==1
        @test shots[1].displacement !== a.mech.u_prev
        @test !isdefined(JuBat,:validate_mechanical_state)
        @test a.mech.mechanics_initialized
        @test_throws ArgumentError JuBat.Solve(c;initial_state=(;a...,y=reshape(copy(a.y),:,1)))
        bad=deepcopy(a); bad.y[1]=NaN
        @test_throws ArgumentError JuBat.Solve(c;initial_state=bad)
        bad=(;a...,t_global=NaN)
        @test_throws ArgumentError JuBat.Solve(c;initial_state=bad)
        bad=deepcopy(a)
        if kind=="gp"
            @test all(iszero(s.D) for s in a.mech.gp_damage_states)
            bad.mech.gp_damage_states=nothing
        else
            bad.mech.plastic_states=nothing
        end
        history_error=kind=="gp" ? ArgumentError : ErrorException
        @test_throws history_error JuBat.Solve(c;initial_state=bad)
        old=state_numbers(a)
        # Same input physical checkpoint on an independently built Case, with a new workspace.
        fresh=fixture_c(kind); same_rt=Dict{String,Any}(); new_rt=Dict{String,Any}()
        same=JuBat.solve_phase(c,JuBat.PHASE_REST,0.03,a;dt_range=[0.01,0.01],variables=same_rt)
        other=JuBat.solve_phase(fresh,JuBat.PHASE_REST,0.03,a;dt_range=[0.01,0.01],variables=new_rt)
        aa=numerical_result(same.solve_result);bb=numerical_result(other.solve_result)
        for key in keys(aa)
            aa[key]==bb[key] || println("C_RESULT_EQUALITY key=",key," isequal=",isequal(aa[key],bb[key])," value=",aa[key])
        end
        @test isequal(aa,bb)
        @test state_numbers(same.final_state)==state_numbers(other.final_state)
        @test state_numbers(same.final_state)==state_numbers(same_rt["accepted"])
        @test state_numbers(other.final_state)==state_numbers(new_rt["accepted"])
        @test same_rt["czm_workspace"] !== new_rt["czm_workspace"]
        @test state_numbers(a)==old
    end
end

@testset "C production preparation failure leaves prior acceptance" begin
    c=fixture_c("gp");c.opt.debug_coupling=true;original_current=c.opt.Current
    rt=Dict{String,Any}(); shots=JuBat.CZMSnapshot[]; before=Ref{Any}(); hist=Ref{Any}(); observed=Ref(false)
    JuBat.mechanical_step_observer[]=ev->begin
        if ev.phase==:after && ev.result.converged
            observed[]=true
            # Real history conversion failure after successful candidate mechanics.
            rt["history"]["cell voltage"]=zeros(2,size(rt["history"]["cell voltage"],2))
            before[]=state_numbers(rt["accepted"]);hist[]=deepcopy(rt["history"])
        end
    end
    try
        @test_throws DimensionMismatch JuBat.solve_phase(c,JuBat.PHASE_CHARGE,0.03,nothing;dt_range=[0.01,0.01],variables=rt,czm_snapshots=shots)
    finally
        JuBat.mechanical_step_observer[]=nothing
    end
    @test observed[]
    @test c.opt.Current===original_current
    @test state_numbers(rt["accepted"])==before[]
    @test rt["history"]==hist[]
    @test isempty(shots)
    # A user callback after publication sees a committed state even if it throws.
    c=fixture_c("gp");rt=Dict{String,Any}(); seen=Ref{Any}()
    cb=ev->if ev.kind==:accepted_commit
        seen[]=state_numbers(rt["accepted"]);error("post-commit observer fixture")
    end
    @test_throws ErrorException JuBat.Solve(c;variables=rt,chemical_step_observer=cb)
    @test state_numbers(rt["accepted"])==seen[]
end

@testset "C actual chemical rejection is transactional" begin
    p=JuBat.ChooseCell("Jellyroll");o=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,thermalmodel="none",per_element_spme=false)
    o.Current=t->t<0.006 ? 0.0 : 0.5;o.time=[0.0,0.025];o.dt=[0.001,0.008];o.dtType="auto";o.dtThreshold=1e-10
    c=JuBat.SetCase(p,o);rt=Dict{String,Any}();shots=JuBat.CZMSnapshot[];saved=Ref{Any}();rejected=Ref(0)
    cb=ev->begin
        if ev.kind==:attempt_begin
            saved[]=(state_numbers(rt["accepted"]),deepcopy(rt["history"]),deepcopy(shots))
        elseif ev.kind==:rejected
            rejected[]+=1
            @test state_numbers(rt["accepted"])==saved[][1]
            @test isequal(rt["history"],saved[][2])
            @test shots==saved[][3]
        end
    end
    JuBat.Solve(c;variables=rt,czm_snapshots=shots,chemical_step_observer=cb)
    @test rejected[]==2
end

@testset "C real mechanical solver failure is transactional" begin
    c=fixture_c("geo");c.opt.czm.max_iter=1;c.opt.czm.update_interval=1
    a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy((haskey(c.mesh,"czm") ? JuBat.MechState(c.mesh["czm"]) : nothing)),t_global=0.0);a.y[c.index["thermal2D temperature"]].=1.2
    rt=Dict{String,Any}();shots=JuBat.CZMSnapshot[];saved=Ref{Any}();failed=Ref(false)
    JuBat.mechanical_step_observer[]=ev->begin
        if ev.phase==:before
            saved[]=(state_numbers(rt["accepted"]),deepcopy(rt["history"]))
        elseif ev.phase==:after
            failed[]=!ev.result.converged
        elseif ev.phase==:error
            failed[]=true
        end
    end
    try
        @test_throws ErrorException JuBat.Solve(c;initial_state=a,variables=rt,czm_snapshots=shots)
    finally
        JuBat.mechanical_step_observer[]=nothing
    end
    @test failed[]
    @test state_numbers(rt["accepted"])==saved[][1]
    @test isequal(rt["history"],saved[][2])
    @test isempty(shots)
end

@testset "C cycling uses a copied seed and inherited temperature" begin
    c=fixture_c("gp");rt=Dict{String,Any}();seed=JuBat.MechState(c.mesh["czm"]);before=mech_numbers(seed)
    opts=JuBat.CycleOption(n_cycles=1,t_charge=0.03,t_rest1=0.03,t_discharge=0.03,t_rest2=0.03,
        SOC_init=0.65,dt_cycle=[0.01,0.01])
    old_current=c.opt.Current;old_time=copy(c.opt.time);old_dt=copy(c.opt.dt)
    r=JuBat.solve_cycling(c,opts,seed;verbose=false,save_detailed=true,variables=rt)
    @test r.n_cycles==1
    @test mech_numbers(seed)==before
    @test r.final_mech !== rt["accepted"].mech
    @test mech_numbers(r.final_mech)==mech_numbers(rt["accepted"].mech)
    @test c.opt.Current===old_current && c.opt.time==old_time && c.opt.dt==old_dt
end

@testset "C active plastic checkpoint carries and unloads" begin
    # Reuse the B default-parameter synthetic opposite-sign SOC load checkpoint,
    # captured on the same ntheta=8 collector topology; no temporary yield factor.
    evidence=deserialize(joinpath(@__DIR__,"..","output","validate_case_batch_b_j2","original_repeat_candidate","repeat.bin"))
    c=fixture_c("geo");c.opt.czm.winding_prestress=false
    c.opt.czm.iter_method="load_substep";c.opt.czm.max_iter=200;c.opt.czm.tol=1e-8
    m=c.mesh["czm"];ms=JuBat.MechState(m)
    @test length(evidence.displacement)==length(ms.u_prev)
    @test size(evidence.states)==(size(m.bulk_element,1),4)
    ms.u_prev=copy(evidence.displacement)
    ms.plastic_states=map(s->JuBat.PlasticState(s.eps_p,s.kappa,s.sigma),evidence.states)
    for (s,values) in zip(ms.damage_states,evidence.damage)
        for (f,v) in zip(fieldnames(JuBat.DamageState),values);setfield!(s,f,v);end
    end
    JuBat.ensure_node_ref!(c,ms)
    println("C_ACTIVE_RESTORE max_kappa=",maximum(s.kappa for s in ms.plastic_states)," max_D=",maximum(s.D for s in ms.damage_states))
    @test maximum(s.kappa for s in ms.plastic_states)>0
    a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy(ms),t_global=0.0)
    ms.mechanics_initialized=true
    a=(;a...,mech=deepcopy(ms))
    a=(;a...,t_global=0.1)
    before=state_numbers(a);branchD=Any[]; rt=Dict{String,Any}()
    JuBat.branch_consumption_observer[]=ev->push!(branchD,copy(ev.D_elem))
    try
        r=JuBat.solve_phase(c,JuBat.PHASE_REST,0.01,a;dt_range=[0.01,0.01],variables=rt)
        @test mech_numbers(r.final_state.mech)==mech_numbers(ms)
        @test r.final_state.mech.node_ref !== ms.node_ref
        @test maximum(s.kappa for s in r.final_state.mech.plastic_states)==maximum(s.kappa for s in ms.plastic_states)
    finally
        JuBat.branch_consumption_observer[]=nothing
    end
    @test state_numbers(a)==before
    fresh=fixture_c("geo");fresh.opt.czm.winding_prestress=false
    fresh.opt.czm.iter_method="load_substep";fresh.opt.czm.max_iter=200;fresh.opt.czm.tol=1e-8
    same=JuBat.solve_phase(c,JuBat.PHASE_REST,0.02,a;dt_range=[0.01,0.01])
    other=JuBat.solve_phase(fresh,JuBat.PHASE_REST,0.02,a;dt_range=[0.01,0.01])
    @test same.final_state.mech.mechanics_initialized
    @test maximum(s.kappa for s in same.final_state.mech.plastic_states)>=maximum(s.kappa for s in ms.plastic_states)
    @test state_numbers(same.final_state)==state_numbers(other.final_state)
    @test isequal(numerical_result(same.solve_result),numerical_result(other.solve_result))
end

@testset "C constitutive nonzero damage seed feeds actual area ON branches" begin
    c=fixture_c("gp");m=c.mesh["czm"];ms=JuBat.MechState(m)
    G=JuBat.gp_history_for_step(ms,m)
    e=1;ip=JuBat.collector_params(c.param,m.cohesive_elements[e].interface_type)
    _,_,_,peak=JuBat.bilinear_traction_state(1.5ip.δ_0,0.0,JuBat.DamageState(),ip,"model1")
    @test 0<peak.D<1
    for g in axes(G,2);G[e,g]=JuBat.clone_damage_state(peak);end
    ms.gp_damage_states=G
    ms.damage_states=JuBat.aggregate_gp_damage_states(G,JuBat.NCweight(m.czm_submesh.mesh_bonded.gs.order)[1])
    a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy(ms),t_global=0.0);before=state_numbers(a);Dseen=Any[]
    JuBat.branch_consumption_observer[]=ev->push!(Dseen,copy(ev.D_elem))
    try
        r=JuBat.solve_phase(c,JuBat.PHASE_REST,0.02,a;dt_range=[0.01,0.01])
        @test maximum(first(Dseen))>0
        @test all(n.D>=old.D for (n,old) in zip(r.final_state.mech.gp_damage_states,a.mech.gp_damage_states))
        @test r.final_state.mech.mechanics_initialized
        @test state_numbers(a)==before
        fresh=fixture_c("gp")
        rr=JuBat.solve_phase(fresh,JuBat.PHASE_REST,0.02,a;dt_range=[0.01,0.01])
        @test state_numbers(r.final_state)==state_numbers(rr.final_state)
        @test isequal(numerical_result(r.solve_result),numerical_result(rr.solve_result))
    finally
        JuBat.branch_consumption_observer[]=nothing
    end
end
end

@testset "Plain state input copying and voltage fail-fast" begin
    @test all(!isdefined(JuBat,n) for n in (:AcceptedState,:AcceptedStamp,:PhaseInputs,:initial_accepted_state,:phase_inputs,:replace_accepted,:accepted_from_step))
    p=JuBat.ChooseCell("Jellyroll");o=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,thermalmodel="none",per_element_spme=false)
    o.Current=_ ->0.5;o.time=[0.0,0.02];o.dt=[0.01,0.01];o.dtType="constant"
    c=JuBat.SetCase(p,o);c.opt.y0=copy(JuBat.ModelInitialisation(c));seed=copy(c.opt.y0);rt=Dict{String,Any}();seen=Ref(false)
    cb=ev->if ev.kind==:init_BE
        seen[]=true
        @test !haskey(rt,"accepted") && !haskey(rt,"history")
        error("copy boundary fixture")
    end
    @test_throws ErrorException JuBat.Solve(c;variables=rt,chemical_step_observer=cb)
    @test seen[] && c.opt.y0==seed && isempty(rt)
    valid=JuBat.Solve(c;return_final_state=true)
    @test valid["final_state"].y !== c.opt.y0 && c.opt.y0==seed
    # A real electrochemical module can produce non-finite voltage with finite y;
    # the numerical-entry guard must stop before history/accepted publication.
    c.param.NE.U=x->NaN;rt=Dict{String,Any}()
    @test_throws ErrorException JuBat.Solve(c;variables=rt)
    @test isempty(rt)
    c=fixture_c("gp");m=c.mesh["czm"];ms=JuBat.MechState(m);w=JuBat.CZMWorkspace(m,c.param);ne=size(m.bulk_element,1)
    @test !ms.mechanics_initialized
    r=JuBat.solve_czm_step(m,ms,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,
        dT_elem=zeros(ne),Δsoc_n_elem=zeros(ne),Δsoc_p_elem=zeros(ne),dt_seconds=0.02)
    @test r.converged && ms.mechanics_initialized
    @test all(iszero(s.D) for s in ms.gp_damage_states)
end
@testset "Initialized mechanics never reconstructs missing histories" begin
    c=fixture_c("gp");m=c.mesh["czm"];w=JuBat.CZMWorkspace(m,c.param);n=size(m.bulk_element,1);z=zeros(n)
    ms=JuBat.MechState(m)
    r=JuBat.solve_czm_step(m,ms,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,dt_seconds=0.02)
    @test r.converged && ms.mechanics_initialized && ms.gp_damage_states!==nothing
    @test all(iszero(s.D) for s in ms.gp_damage_states)
    bad=deepcopy(ms);bad.gp_damage_states=nothing;before=mech_numbers(bad)
    @test_throws ArgumentError JuBat.solve_czm_step(m,bad,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,dt_seconds=0.02)
    @test mech_numbers(bad)==before
    c=fixture_c("geo");c.opt.czm.winding_prestress=false;m=c.mesh["czm"];w=JuBat.CZMWorkspace(m,c.param);n=size(m.bulk_element,1);z=zeros(n);ms=JuBat.MechState(m)
    r=JuBat.solve_czm_step(m,ms,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,eigenstrain=(dT=z,Δsn=z,Δsp=z),dt_seconds=0.02)
    @test r.converged && ms.mechanics_initialized && ms.plastic_states!==nothing
    bad=deepcopy(ms);bad.plastic_states=nothing;before=mech_numbers(bad)
    @test_throws ArgumentError JuBat.solve_czm_step(m,bad,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,eigenstrain=(dT=z,Δsn=z,Δsp=z),dt_seconds=0.02)
    @test mech_numbers(bad)==before
    c=fixture_c("geo");m=c.mesh["czm"];w=JuBat.CZMWorkspace(m,c.param);n=size(m.bulk_element,1);z=zeros(n);ms=JuBat.MechState(m);before=mech_numbers(ms)
    @test_throws ArgumentError JuBat.solve_czm_step(m,ms,c.param,zeros(2size(m.node,1)),c.opt.czm;czm_workspace=w,dT_elem=z,Δsoc_n_elem=z,Δsoc_p_elem=z,eigenstrain=(dT=z,Δsn=z,Δsp=z),dt_seconds=0.02)
    @test mech_numbers(ms)==before
    y=JuBat.ModelInitialisation_MultiSPMe(c)
    vars=JuBat.CallModel(c,y,0.0;jacobi="update",mech=ms)[4]
    r=JuBat.update_czm_damage!(c,ms,vars,JuBat.get_thermal_dofs(y,c);dt_seconds=0.02,czm_workspace=w)
    @test r.converged && ms.mechanics_initialized && ms.winding_prestress!==nothing
    @test any(any(!iszero(v) for v in s) for s in ms.winding_prestress)
    bad=deepcopy(ms);bad.winding_prestress=nothing;before=mech_numbers(bad)
    @test_throws ArgumentError JuBat.update_czm_damage!(c,bad,vars,JuBat.get_thermal_dofs(y,c);dt_seconds=0.02,czm_workspace=w)
    @test mech_numbers(bad)==before
end
