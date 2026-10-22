using Test
isdefined(Main,:fixture_c) || include(joinpath(@__DIR__,"..","tools","validate_case_batch_c_trace.jl"))
@testset "Enabled CZM cannot silently lose its topology" begin
    c=fixture_c("export")
    c.opt.czm.enabled=true
    c.opt.czm.area_loss_enabled=false
    old=(y=[123.0],mech=nothing,t_global=7.0)
    history=Dict("sentinel"=>[9.0])
    runtime=Dict{String,Any}("accepted"=>old,"history"=>history)
    @test_throws KeyError JuBat.Solve(c;variables=runtime)
    @test runtime["accepted"]===old
    @test runtime["history"]===history
    time=c.opt.time;dt=c.opt.dt;current=c.opt.Current
    @test_throws KeyError JuBat.solve_phase(c,JuBat.PHASE_DISCHARGE,0.03,nothing;variables=runtime)
    @test c.opt.time===time && c.opt.dt===dt && c.opt.Current===current
    @test runtime["accepted"]===old && runtime["history"]===history
    @test_throws KeyError JuBat.solve_cycling(c,JuBat.CycleOption();variables=runtime,verbose=false)
    @test runtime["accepted"]===old && runtime["history"]===history
end
