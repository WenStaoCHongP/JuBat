using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","tools","validate_case_batch_c_trace.jl"))
const TC_OUT=normpath(joinpath(@__DIR__,"..","output","test_temperature_current"))
mkpath(TC_OUT)

if !("--dt-scope-only" in ARGS)
@testset "Fixed temperature reference and phase waveform" begin
    @test !isdefined(JuBat,:AcceptedState)
    @test !isdefined(JuBat,:PhaseInputs)
    @test !isdefined(JuBat,:reset_cycle_temperature!)
    @test !(:reset_T_each_cycle in fieldnames(JuBat.CycleOption))
    @test !(:reset_T_before_charge in fieldnames(JuBat.CycleOption))
    @test !(:I_charge in fieldnames(JuBat.CycleOption))
    @test !(:I_discharge in fieldnames(JuBat.CycleOption))
    c=fixture_c("export");wave=t->0.5+2t;c.opt.Current=wave;T0=c.param.cell.T0
    a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy((haskey(c.mesh,"czm") ? JuBat.MechState(c.mesh["czm"]) : nothing)),t_global=0.0);idx=c.index["thermal2D temperature"]
    a.y[idx].=T0+0.02
    for (phase,sgn) in ((JuBat.PHASE_DISCHARGE,1),(JuBat.PHASE_REST,0),(JuBat.PHASE_CHARGE,-1))
        before=copy(a.y);runtime=Dict{String,Any}()
        r=JuBat.solve_phase(c,phase,0.03,a;dt_range=[0.01,0.01],variables=runtime)
        @test c.opt.Current===wave
        sr=r.solve_result;times=sr["time [s]"];current=sr["cell current [A]"]
        @test sr["thermal2D temperature at nodes [K]"][:,1]==before[idx]*c.param.scale.T_ref
        @test all(isapprox.(current,sgn.*wave.(times);rtol=8eps(Float64),atol=0.0))
        expected=sgn==0 ? 0.0 : (0.5*times[end]+times[end]^2)/3600
        @test isapprox(r.capacity,expected;rtol=8eps(Float64),atol=0.0)
        @test c.param.cell.T0==T0 && a.y==before
        @test maximum(JuBat.get_thermal_dofs(r.final_state.y,c))>T0+0.01
        a=r.final_state
    end
    # Both public entries restore the identical function when initialization fails.
    bad=(;a...,y=a.y[1:end-1])
    @test_throws DimensionMismatch JuBat.solve_phase(c,JuBat.PHASE_CHARGE,0.03,bad;dt_range=[0.01,0.01])
    @test c.opt.Current===wave
    @test_throws DimensionMismatch JuBat.solve_phase_with_export(c,JuBat.PHASE_REST,0.03,bad;dt_range=[0.01,0.01])
    @test c.opt.Current===wave
end

@testset "Export uses physical amperes and inherited accepted temperature" begin
    c=fixture_c("export");wave=t->0.5+2t;c.opt.Current=wave
    @test c.param.scale.I_typ!=1
    a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy((haskey(c.mesh,"czm") ? JuBat.MechState(c.mesh["czm"]) : nothing)),t_global=0.0);idx=c.index["thermal2D temperature"];a.y[idx].+=0.02
    oldT=copy(a.y[idx]);T0=c.param.cell.T0
    for (phase,sgn) in ((JuBat.PHASE_DISCHARGE,1),(JuBat.PHASE_REST,0),(JuBat.PHASE_CHARGE,-1))
        rt=Dict{String,Any}();r,samples=JuBat.solve_phase_with_export(c,phase,0.03,a;dt_range=[0.01,0.01],variables=rt)
        @test c.opt.Current===wave
        local_times=[s.time-a.t_global for s in samples]
        @test all(isapprox.([s.I for s in samples],sgn.*wave.(local_times);rtol=8eps(Float64),atol=0.0))
        expected=sgn==0 ? 0.0 : (0.5*0.03+0.03^2)/3600
        @test isapprox(r.capacity,expected;rtol=8eps(Float64),atol=0.0)
        @test c.param.cell.T0==T0
        @test a.y[idx]==oldT
        @test minimum(r.final_state.y[idx])>T0+0.01
        @test rt["accepted"].y==r.final_state.y
        a=r.final_state;oldT=copy(a.y[idx])
    end
end

@testset "Two real microcycles inherit every phase and cycle seam" begin
    c=fixture_c("export");wave=t->0.5+2t;c.opt.Current=wave
    opts=JuBat.CycleOption(n_cycles=2,t_discharge=0.03,t_rest1=0.03,t_charge=0.03,t_rest2=0.03,SOC_init=0.65,dt_cycle=[0.01,0.01])
    rt=Dict{String,Any}();r=JuBat.solve_cycling(c,opts;verbose=false,save_detailed=true,variables=rt)
    @test r.n_cycles==2 && c.opt.Current===wave
    previous=nothing;changed=false
    for cycle in r.cycle_results
        for phase in (cycle.discharge,cycle.rest1,cycle.charge,cycle.rest2)
            idx=c.index["thermal2D temperature"]
            firstT=phase.solve_result["thermal2D temperature at nodes [K]"][:,1]
            if previous!==nothing
                @test firstT==previous*c.param.scale.T_ref
                changed |= any(previous .!= c.param.cell.T0)
            end
            previous=copy(phase.final_state.y[idx])
        end
    end
    @test changed
    @test rt["accepted"].y[c.index["thermal2D temperature"]]==previous
end

@testset "Lumped temperature inherits; no-thermal model uses fixed T0" begin
    for model in ("SPM","SPMe","P2D")
        o=JuBat.Option(model=model,Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,thermalmodel="none",per_element_spme=false)
        c=JuBat.SetCase(JuBat.ChooseCell("Jellyroll"),o);y=JuBat.ModelInitialisation(c)
        vars=model=="SPM" ? JuBat.SPM_variables(c,y,0.0) : model=="SPMe" ? JuBat.SPMe_variables(c,y,0.0) : JuBat.P2D_variables(c,y,0.0)
        @test vars["temperature"]==c.param.cell.T0
    end
    o=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,thermalmodel="lumped",per_element_spme=false)
    o.Current=_ ->0.5;o.thermal_enabled=true
    c=JuBat.SetCase(JuBat.ChooseCell("Jellyroll"),o);a=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy((haskey(c.mesh,"czm") ? JuBat.MechState(c.mesh["czm"]) : nothing)),t_global=0.0);a.y[c.index["temperature"]].+=0.02
    for phase in (JuBat.PHASE_DISCHARGE,JuBat.PHASE_REST,JuBat.PHASE_CHARGE)
        before=copy(a.y);r=JuBat.solve_phase(c,phase,0.03,a;dt_range=[0.01,0.01])
        @test r.solve_result["temperature [K]"][1]==only(before[c.index["temperature"]])*c.param.scale.T_ref
        a=r.final_state
    end
end

end

@testset "Phase time and dt scope restore original configuration references" begin
    for exported in (false,true)
        c=fixture_c("export");original_dt=[0.01 0.02];c.opt.dt=original_dt
        original_time=[2.0 3.0];c.opt.time=original_time;original_time_contents=copy(original_time)
        original_contents=copy(original_dt);a=nothing
        for phase in (JuBat.PHASE_DISCHARGE,JuBat.PHASE_CHARGE)
            phase_dt=[0.005,0.005];seen=Bool[]
            c.opt.Current=t->begin
                push!(seen,c.opt.dt===phase_dt && size(c.opt.time)==(1,2) && c.opt.time==[0.0 0.02])
                0.5
            end
            original_current=c.opt.Current
            result=exported ? JuBat.solve_phase_with_export(c,phase,0.02,a;dt_range=phase_dt) : JuBat.solve_phase(c,phase,0.02,a;dt_range=phase_dt)
            @test !isempty(seen) && all(seen)
            @test c.opt.dt===original_dt && c.opt.dt==original_contents
            @test c.opt.time===original_time && c.opt.time==original_time_contents
            @test phase_dt==[0.005,0.005] && c.opt.Current===original_current
            times=exported ? [x.time for x in result[2]] : vec(result.solve_result["time [s]"])
            @test length(times)>=3 && all(isapprox.(diff(times),0.005;rtol=8eps(Float64),atol=0.0))
            phase_result=exported ? result[1] : result
            @test isapprox(phase_result.duration,0.02;rtol=8eps(Float64),atol=0.0)
            a=phase_result.final_state
        end
        bad=(;a...,y=a.y[1:end-1]);phase_dt=[0.005,0.005]
        if exported
            @test_throws DimensionMismatch JuBat.solve_phase_with_export(c,JuBat.PHASE_CHARGE,0.02,bad;dt_range=phase_dt)
        else
            @test_throws DimensionMismatch JuBat.solve_phase(c,JuBat.PHASE_CHARGE,0.02,bad;dt_range=phase_dt)
        end
        @test c.opt.dt===original_dt && c.opt.dt==original_contents
        @test c.opt.time===original_time && c.opt.time==original_time_contents
        times_seen=Float64[]
        wave=t->begin
            push!(times_seen,t)
            t>0.0 && error("dt scope mid-solve fixture")
            0.5
        end
        c.opt.Current=wave
        if exported
            @test_throws ErrorException JuBat.solve_phase_with_export(c,JuBat.PHASE_DISCHARGE,0.02,a;dt_range=phase_dt)
        else
            @test_throws ErrorException JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.02,a;dt_range=phase_dt)
        end
        @test any(>(0.0),times_seen)
        @test c.opt.dt===original_dt && c.opt.dt==original_contents
        @test c.opt.time===original_time && c.opt.time==original_time_contents
        @test c.opt.Current===wave && phase_dt==[0.005,0.005]
    end
end
