using Test, Logging
include(joinpath(@__DIR__,"..","tools","validate_case_batch_c_trace.jl"))
using .JuBat
const H=Dict{String,Union{Array{Float64},Float64}}
const OUT=normpath(joinpath(@__DIR__,"..","output","test_history_update"))
mkpath(OUT)
if !("--production-only" in ARGS)
@testset "Complete history columns fail without mutation" begin
    for debug in (false,true), grow in (false,true), failure in (:missing,:length,:scalar)
        hist=H("x"=>[1.0 2.0;3.0 4.0],"y"=>[5.0 6.0])
        before=deepcopy(hist);aliases=copy(hist);v=grow ? 3 : 2
        values=failure==:missing ? H("x"=>[7.0,8.0]) : failure==:length ? H("x"=>[7.0,8.0,9.0],"y"=>10.0) : H("x"=>7.0,"y"=>10.0)
        expected=failure==:missing ? KeyError : DimensionMismatch
        @test_throws expected JuBat.Variable_update!(hist,values,v;debug_coupling=debug)
        @test hist==before
        @test all(hist[k]===aliases[k] for k in keys(hist))
    end
end
end
@testset "Production missing required column after mechanics is transactional" begin
    c=fixture_c("gp");c.opt.debug_coupling=false
    rt=Dict{String,Any}();shots=JuBat.CZMSnapshot[];observed=Ref(false)
    accepted=Ref{Any}();history=Ref{Any}();contents=Ref{Any}();aliases=Ref{Any}();y_before=Ref{Any}();mech_before=Ref{Any}()
    JuBat.mechanical_step_observer[]=ev->begin
        if ev.phase==:after && ev.result.converged
            observed[]=true
            # Controlled output-schema injection: existing fields remain intact,
            # while the candidate producer deliberately has no value for this key.
            rt["history"]["fixture required missing output"]=zeros(1,size(rt["history"]["time"],2))
            accepted[]=rt["accepted"];history[]=rt["history"];contents[]=deepcopy(history[]);aliases[]=copy(history[])
            y_before[]=copy(accepted[].y);mech_before[]=mech_numbers(accepted[].mech)
        end
    end
    try
        @test_throws KeyError JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.03,nothing;dt_range=[0.01,0.01],variables=rt,czm_snapshots=shots)
    finally
        JuBat.mechanical_step_observer[]=nothing
    end
    @test observed[]
    @test rt["accepted"]===accepted[] && rt["accepted"].y==y_before[]
    @test mech_numbers(rt["accepted"].mech)==mech_before[]
    @test rt["history"]===history[] && isequal(rt["history"],contents[])
    @test all(rt["history"][k]===aliases[][k] for k in keys(aliases[]))
    @test isempty(shots)
end
