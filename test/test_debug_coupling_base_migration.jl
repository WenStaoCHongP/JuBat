using Test, Logging
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

const H = Dict{String,Union{Array{Float64},Float64}}
const OUT = joinpath(@__DIR__, "..", "output", "test_debug_coupling_base_migration")
mkpath(OUT)

function captured(f)
    try
        return (value=f(), error=nothing)
    catch err
        return (value=nothing, error=err)
    end
end

function history_calls!(found, node)
    node isa Expr || return found
    if node.head == :call && node.args[1] == :Variable_update!
        push!(found, node)
    end
    for child in node.args
        history_calls!(found, child)
    end
    return found
end

function simple_case(debug; invalid_window=false)
    p = JuBat.ChooseCell("Jellyroll")
    if invalid_window
        p.cell.v_l = 5.0
        p.cell.v_h = 4.2
    end
    opt = JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4)
    opt.Current = t -> 0.5
    opt.time = [0.0,0.03]; opt.dt = [0.01,0.01]; opt.dtType = "constant"
    opt.thermalmodel = "none"; opt.mechanicalmodel = "none"
    opt.debug_coupling = debug
    return JuBat.SetCase(p,opt)
end

@testset "P1 on native 3df history updates" begin
    @test !JuBat.Option().debug_coupling
    @testset "both native Solve calls forward the real option" begin
        source = read(joinpath(@__DIR__,"..","src","Solve.jl"),String)
        calls = history_calls!(Expr[], Meta.parseall(source))
        @test length(calls) == 2
        for call in calls
            parameters = filter(x -> x isa Expr && x.head == :parameters, call.args)
            keywords = Any[]
            for parameter in parameters
                append!(keywords,filter(x -> x isa Expr && x.head == :kw && x.args[1] == :debug_coupling,parameter.args))
            end
            @test length(keywords) == 1
            @test all(x -> x.args[2] == :(case.opt.debug_coupling), keywords)
        end
    end

    @testset "default bypasses diagnostics while true preserves original errors" begin
        hist = H("x"=>zeros(2,1))
        result = captured(() -> JuBat.Variable_update!(hist,H("x"=>[1.0 2.0]),1))
        @test result.error === nothing
        @test hist["x"][:,1] == [1.0,2.0]
        bad_cases = [
            (H("x"=>[0.0,0.0]), H("x"=>[1.0,2.0]), DimensionMismatch,
                "history variable \"x\" must be a matrix"),
            (H("x"=>zeros(2,1)), H("x"=>[1.0 2.0]), DimensionMismatch,
                "variable \"x\" must have exactly one column, got size (1, 2)"),
            (H("x"=>zeros(2,1)), H("x"=>reshape([1.0,2.0],2,1,1)), DimensionMismatch,
                "variable \"x\" must be a vector or single-column matrix, got 3 dimensions"),
            (H("x"=>zeros(2,1)), H("x"=>[1.0]), DimensionMismatch,
                "variable \"x\" length 1 does not match history row count 2"),
            (H("x"=>zeros(2,1)), H("x"=>1.0), DimensionMismatch,
                "scalar variable \"x\" cannot populate history with 2 rows"),
            (H("x"=>1.0), H("x"=>[2.0]), ArgumentError,
                "scalar history variable \"x\" requires a Float64 value, got Vector{Float64}")]
        for (hist,vars,kind,message) in bad_cases
            result = captured(() -> JuBat.Variable_update!(hist,vars,1;debug_coupling=true))
            @test result.error isa kind
            @test result.error isa kind && result.error.msg == message
        end
        for value in ([1.0 2.0], reshape([1.0,2.0],2,1,1))
            hist = H("x"=>zeros(2,1))
            result = captured(() -> JuBat.Variable_update!(hist,H("x"=>value),1;debug_coupling=false))
            @test result.error === nothing
            @test hist["x"][:,1] == [1.0,2.0]
        end
    end

    @testset "valid copies aliases growth and partial updates retain native contracts" begin
        finals = H[]
        for debug in (false,true), matrix_input in (false,true)
            hist=H("x"=>[1.0 3.0;2.0 4.0],"scalar"=>[5.0 6.0],"latest"=>7.0,"untouched"=>[20.0 21.0])
            before=deepcopy(hist); old_alias=hist["x"]
            data=matrix_input ? reshape([8.0,9.0],2,1) : [8.0,9.0]
            vars=H("x"=>data,"scalar"=>10.0,"latest"=>11.0,"extra"=>[1.0])
            result=captured(() -> JuBat.Variable_update!(hist,vars,2;debug_coupling=debug))
            @test result.error === nothing && result.value === hist
            @test hist["x"]==[1.0 8.0;2.0 9.0] && hist["x"]===old_alias
            @test hist["scalar"]==[5.0 10.0] && hist["latest"]==11.0
            @test hist["untouched"]==before["untouched"] && !haskey(hist,"extra")
            data .= -100.0
            @test hist["x"][:,2]==[8.0,9.0]
            data .= 12.0
            result=captured(() -> with_logger(NullLogger()) do
                JuBat.Variable_update!(hist,vars,3;debug_coupling=debug)
            end)
            @test result.error === nothing && result.value === hist
            @test size(hist["x"])==(2,1002) && hist["x"]!==old_alias
            @test hist["x"][:,1:2]==[1.0 8.0;2.0 9.0] && old_alias==[1.0 8.0;2.0 9.0]
            @test size(hist["x"],2)>=3 && hist["x"][:,3]==[12.0,12.0] && all(iszero,hist["x"][:,4:end])
            push!(finals,hist)
        end
        @test all(h -> h==first(finals),finals)
        for debug in (false,true)
            empty=H()
            result=captured(() -> JuBat.Variable_update!(empty,H(),1;debug_coupling=debug))
            @test result.error === nothing && result.value===empty
            hist=H("x"=>zeros(1,1),"y"=>zeros(1,2))
            result=captured(() -> JuBat.Variable_update!(hist,H("x"=>[1.0]),1;debug_coupling=debug))
            @test result.error === nothing && hist["x"][1]==1.0 && hist["y"]==zeros(1,2)
            for bad in (NaN,Inf)
                hist=H("x"=>zeros(1,1))
                result=captured(() -> JuBat.Variable_update!(hist,H("x"=>[bad]),1;debug_coupling=debug))
                @test result.error === nothing && isequal(hist["x"][1],bad)
            end
        end
        result=captured(() -> JuBat.Variable_update!(H("x"=>zeros(2,1)),H("x"=>[1.0,2.0]),0))
        @test result.error isa BoundsError
        result=captured(() -> JuBat.Variable_update!(H("x"=>zeros(2,1)),H("x"=>[1.0,2.0,3.0]),1))
        @test result.error isa DimensionMismatch
    end

    @testset "source and model legality gates remain unconditional" begin
        for debug in (false,true)
            @test_throws ArgumentError JuBat.Solve(simple_case(debug;invalid_window=true))
            @test_throws ArgumentError JuBat.Solve(simple_case(debug);initial_state=Dict{String,Any}())
            @test_throws ArgumentError JuBat.Solve(simple_case(debug);initial_state=Dict{String,Any}("y"=>1.0))
            case=simple_case(debug)
            mesh=JuBat.jellyroll_collector_seed_mesh(case.param;nθ=8,czm_enabled=false,gsorder=2)
            case=JuBat.setup_thermal2D_mesh(case,mesh)
            case.opt.model="thermal"
            @test_throws KeyError JuBat.Solve(case;thermal_variables=Dict{String,Any}())
        end
    end
end

@testset "native short Solve is identical with diagnostics off and on" begin
    outputs=Any[]; current_calls=Vector{Float64}[]; observer_events=Any[]
    for debug in (false,true)
        calls=Float64[]; events=Any[]
        case=simple_case(debug)
        case.opt.Current=t->(push!(calls,t);0.5)
        result=JuBat.Solve(case;return_final_state=true,
            chemical_step_observer=e->push!(events,deepcopy(e)))
        push!(outputs,Dict(k=>v for (k,v) in result if !startswith(k,"timing ") || k=="timing CallModel calls"))
        push!(current_calls,calls);push!(observer_events,events)
        @test all(isfinite,result["cell voltage [V]"])
        @test result["termination_reason"]=="time_limit"
    end
    @test isequal(outputs[1],outputs[2])
    @test current_calls[1]==current_calls[2]
    @test isequal(observer_events[1],observer_events[2])
end
