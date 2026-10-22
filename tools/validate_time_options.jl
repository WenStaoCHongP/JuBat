using Test
include(joinpath(@__DIR__,"..","src","JuBat.jl"))
using .JuBat
const OUT=normpath(joinpath(@__DIR__,"..","output","validate_time_options"))
mkpath(OUT)
function fixture(thermal, row)
    o=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,thermalmodel=thermal,per_element_spme=false)
    o.time=row ? [0.0 0.03] : [0.0,0.03]
    o.dt=row ? [0.01 0.01] : [0.01,0.01]
    o.Current=_ ->0.5;o.thermal_enabled=thermal!="none"
    JuBat.SetCase(JuBat.ChooseCell("Jellyroll"),o)
end
if !("--thermal-only" in ARGS)
@testset "Direct Solve preserves Option array shapes" begin
    @test size(JuBat.Option().time)==(1,2)
    for thermal in ("lumped","none")
        results=[]
        for row in (false,true)
            c=fixture(thermal,row);t=copy(c.opt.time);d=copy(c.opt.dt)
            r=JuBat.Solve(c;return_final_state=true)
            @test c.opt.time==t && size(c.opt.time)==size(t)
            @test c.opt.dt==d && size(c.opt.dt)==size(d)
            @test all(isfinite,r["final_state"].y)
            push!(results,r)
        end
        @test results[1]["final_state"].y==results[2]["final_state"].y
        @test isequal(Dict(k=>v for (k,v) in results[1] if !startswith(k,"timing ") && k!="final_state"),Dict(k=>v for (k,v) in results[2] if !startswith(k,"timing ") && k!="final_state"))
        c=fixture(thermal,true);t=copy(c.opt.time);d=copy(c.opt.dt)
        r=JuBat.Solve(c;return_final_state=true)
        @test r["final_state"].y==results[1]["final_state"].y
        @test c.opt.time==t && c.opt.dt==d
    end
end
end
@testset "Pure thermal default row-array entry" begin
    c=fixture("distributed2D",true)
    md=JuBat.jellyroll_collector_seed_mesh(c.param;nθ=20,gsorder=2)
    c=JuBat.setup_thermal2D_mesh(c,md);c.opt.model="thermal"
    mesh=c.mesh["thermal2D"]
    vars=Dict{String,Any}("thermal2D temperature at nodes"=>fill(c.param.cell.T0,mesh.nlen),"heat_source_fields"=>zeros(size(mesh.element,1)))
    t=copy(c.opt.time);d=copy(c.opt.dt)
    r=JuBat.Solve(c;thermal_variables=vars,thermal_record=true)
    @test length(r.time)>=3 && all(isfinite,r.T_nodes)
    @test size(r.T_hist,1)==mesh.nlen
    @test c.opt.time==t && c.opt.dt==d
end
