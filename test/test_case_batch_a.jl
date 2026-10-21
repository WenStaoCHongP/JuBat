using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

const OUT_A = normpath(joinpath(@__DIR__, "..", "output", "test_case_batch_a"))
mkpath(OUT_A)

function case_a(; model="SPMe", thermal="distributed2D", multi=true, y0=Float64[])
    opt = JuBat.Option(model=model, Nn=3, Ns=2, Np=3, Nrn=4, Nrp=4,
        thermalmodel=thermal, thermal_enabled=thermal != "none", per_element_spme=multi)
    opt.y0 = y0
    opt.Current = t -> 0.5
    opt.time = [0.0, 0.02]; opt.dt = [0.01, 0.01]; opt.dtType = "constant"
    opt.mechanicalmodel = "none"; opt.czm.enabled = false
    c = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), opt)
    md = JuBat.jellyroll_collector_seed_mesh(c.param; nθ=12, gsorder=2)
    return JuBat.setup_thermal2D_mesh(c, md), md
end

@testset "A fixed mesh and global state contracts" begin
    c, md = case_a()
    m = c.mesh["thermal2D"]; ne = size(m.element, 1); nT = m.nlen
    nc = c.mesh["negative particle"].nlen + c.mesh["positive particle"].nlen + c.mesh["electrolyte"].nlen
    @test c.mesh isa Dict{String,Any}
    @test !hasfield(JuBat.Case, :layout) && !hasfield(JuBat.Case, :geometry)
    @test c.index["SPMe chemical block length"] == nc
    @test c.index["thermal2D temperature"] == collect(ne*nc+1:ne*nc+nT)
    # Independent sum uses the former layout integration order.
    original_areas = zeros(ne)
    for g in eachindex(m.gs.detJ)
        original_areas[m.gs.ele[g]] += m.gs.weight[g] * m.gs.detJ[g]
    end
    @test c.mesh["thermal2D element area"] == original_areas
    @test c.mesh["thermal2D layer weights"] == last(JuBat.jellyroll_element_properties(m,c.param))
    @test isempty(c.mesh["thermal2D interface pairs"])
    unmerged = JuBat.setup_thermal2D_mesh(c, md; use_merged=false)
    @test unmerged.mesh["thermal2D interface pairs"] == md.interface_pairs
    @test !isempty(md.interface_pairs)

    y = JuBat.ModelInitialisation_MultiSPMe(c)
    chem = vcat(fill(c.param.NE.cs0,c.mesh["negative particle"].nlen),
        fill(c.param.PE.cs0,c.mesh["positive particle"].nlen),fill(c.param.EL.ce0,c.mesh["electrolyte"].nlen))
    @test y == vcat(repeat(chem,ne),fill(c.param.cell.T0,nT))
    @test JuBat.extract_element_state(y, ne, c) == chem
    @test JuBat.get_thermal_dofs(y, c) == fill(c.param.cell.T0,nT)
    yn = JuBat.update_state(y,c; element_index=ne,element_state=chem.+0.001,thermal_nodes=fill(1.01,nT))
    @test JuBat.extract_element_state(yn,ne,c) == chem.+0.001
    @test JuBat.get_thermal_dofs(yn,c) == fill(1.01,nT)
    @test y == vcat(repeat(chem,ne),fill(c.param.cell.T0,nT))
    @test_throws AssertionError JuBat.update_state(y,c; element_index=ne+1,element_state=chem)
    @test_throws AssertionError JuBat.update_state(y,c; element_index=1,element_state=chem[1:end-1])
    @test_throws AssertionError JuBat.update_state(y,c; thermal_nodes=ones(nT-1))
    @test_throws BoundsError JuBat.get_thermal_dofs(y[1:end-1],c)
    missing = deepcopy(c); delete!(missing.index,"SPMe chemical block length")
    @test_throws KeyError JuBat.CallModel_MultiSPMe(missing,y,0.0;jacobi="update")
    @test_throws DimensionMismatch JuBat.Solve(c; initial_state=Dict{String,Any}("y"=>y[1:end-1]))
    @test_throws DimensionMismatch JuBat.solve_phase_with_export(c,JuBat.PHASE_REST,0.02,0.0,Dict{String,Any}("y"=>y[1:end-1]))

    custom = chem .* 0.99
    cc, _ = case_a(y0=reshape(custom,:,1))
    yc = JuBat.ModelInitialisation_MultiSPMe(cc)
    @test cc.opt.y0 == reshape(custom,:,1)
    @test yc == vcat(repeat(custom,ne),fill(cc.param.cell.T0,nT))
    long, _ = case_a(y0=vcat(custom,0.123))
    @test long.index["SPMe chemical block length"] == nc+1
    @test JuBat.ModelInitialisation_MultiSPMe(long) == vcat(repeat(vcat(custom,0.123),ne),fill(long.param.cell.T0,nT))

    # Local SPMe must consume only chemistry despite global indices beyond its length.
    for T in (c.param.cell.T0, c.param.cell.T0+0.01)
        a = JuBat.SPMe_element(c,chem,0.0,1; I_e=0.1,T_e=T)
        b = JuBat.SPMe_element(c,chem,0.0,1; I_e=0.1,T_e=T,workspace=JuBat.create_element_workspace(c))
        @test a[1:3] == b[1:3]
        @test a[4]["temperature"] == b[4]["temperature"] == T
        for k in ("cell voltage","negative particle lithium concentration","positive particle lithium concentration",
                  "electrolyte lithium concentration","negative electrode overpotential","positive electrode overpotential")
            @test a[4][k] == b[4][k]
        end
        @test !haskey(a[4],"SPMe chemical block length")
        @test !haskey(a[4],"thermal2D temperature")
    end
    # Normal/external Solve paths exercise global assembly and temperature extraction.
    a = JuBat.Solve(deepcopy(c);return_final_state=true)
    b = JuBat.Solve(deepcopy(c);initial_state=Dict{String,Any}("y"=>copy(y)),return_final_state=true)
    @test a["final_state"]["y"] == b["final_state"]["y"]
    @test a["cell voltage [V]"] == b["cell voltage [V]"]
    phase = JuBat.solve_phase_with_export(deepcopy(c),JuBat.PHASE_REST,0.02,0.0,
        Dict{String,Any}("y"=>copy(y),"t_global"=>0.0);dt_range=[0.01,0.01])
    @test length(phase[1].final_state["y"]) == length(y)
    state = Dict{String,Any}("y"=>copy(yn))
    JuBat.reset_cycle_temperature!(c,state)
    @test JuBat.get_thermal_dofs(state["y"],c) == fill(c.param.cell.T0,nT)
    @test all(isfinite,JuBat.state_concentration_variance(c,y))

    # Current production CSV is intentionally skipped for zero effective parameters.
    # Its unreachable row body is unchanged except for the mesh-owned pair lookup.
    for (tag, obj) in (("merged",c),("unmerged",unmerged))
        dir = joinpath(OUT_A,tag); mkpath(dir)
        JuBat.write_cohesive_driving_force_csv(nothing,obj,(cohesive_elements=[nothing],),dir,true,JuBat.CsvExportOptions())
        @test !isfile(joinpath(dir,"cohesive_driving_force.csv"))
    end
end

@testset "Single-model indices and temperature remain local" begin
    for model in ("SPM","SPMe","P2D"), thermal in ("none","lumped","distributed2D")
        c, _ = case_a(model=model,thermal=thermal,multi=false)
        @test !haskey(c.index,"SPMe chemical block length")
        @test !haskey(c.index,"thermal2D temperature")
        y = JuBat.ModelInitialisation(c)
        vars = model == "SPM" ? JuBat.SPM_variables(c,y,0.0) :
            model == "SPMe" ? JuBat.SPMe_variables(c,y,0.0) : JuBat.P2D_variables(c,y,0.0)
        @test vars["temperature"] == c.param.cell.T0
        @test all(isfinite,vars["negative particle lithium concentration"])
        if thermal != "none"
            @test only(c.index["temperature"]) <= length(y)
        end
    end
end
