using Test

include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

@testset "temperature result has a fixed shape" begin
    source = read(joinpath(@__DIR__, "..", "src", "Variables.jl"), String)
    @test occursin("variables[\"temperature\"] = zeros(Float64, 1, num)", source)
    @test !occursin("n_temperature = haskey(case.index, \"temperature\")", source)
end

@testset "temperature source follows configuration" begin
    for source_file in ("SPM.jl", "SPMe.jl", "P2D.jl")
        source = read(joinpath(@__DIR__, "..", "src", source_file), String)
        @test !occursin("haskey(case.index, \"temperature\")", source)
        @test !occursin("\"temperature\" in var_list", source)
    end

    param_dim = JuBat.ChooseCell("LG M50")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "none"
    case = JuBat.SetCase(param_dim, opt)
    state = JuBat.ModelInitialisation(case)
    variables_none = JuBat.SPMe_variables(case, state, 0.0)
    @test variables_none["temperature"] == case.param.cell.T0

    opt_lumped = JuBat.Option()
    opt_lumped.model = "SPMe"
    opt_lumped.thermalmodel = "lumped"
    case_lumped = JuBat.SetCase(param_dim, opt_lumped)
    state_lumped = JuBat.ModelInitialisation(case_lumped)
    expected = case_lumped.param.cell.T0 + 0.5
    state_lumped[case_lumped.index["temperature"]] .= expected
    variables_lumped = JuBat.SPMe_variables(case_lumped, state_lumped, 0.0)
    @test variables_lumped["temperature"] == expected

    # 单元温度注入优先于状态温度
    variables_supplied = JuBat.SPMe_variables(case_lumped, state_lumped, 0.0; T_e=expected + 0.25)
    @test variables_supplied["temperature"] == expected + 0.25
end

@testset "models without thermal state" begin
    for model in ("SPM", "SPMe", "P2D")
        param_dim = JuBat.ChooseCell("LG M50")
        opt = JuBat.Option()
        opt.model = model
        opt.thermalmodel = "none"
        opt.Current = _ -> 0.0

        case = JuBat.SetCase(param_dim, opt)
        @test !haskey(case.index, "temperature")

        variables = JuBat.StandardVariables(case, 2)
        @test size(variables["temperature"]) == (1, 2)
    end
end

@testset "standalone SPMe uses ambient temperature" begin
    param_dim = JuBat.ChooseCell("LG M50")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "none"
    opt.Current = _ -> 0.0

    case = JuBat.SetCase(param_dim, opt)
    state = JuBat.ModelInitialisation(case)
    variables = JuBat.SPMe_variables(case, state, 0.0)

    @test variables["temperature"] == case.param.cell.T0
end

@testset "mechanical example uses repository data" begin
    script = read(joinpath(@__DIR__, "..", "example", "mechanical_example.jl"), String)
    @test occursin("joinpath(@__DIR__", script)
    @test !occursin("C:\\Users\\", script)
end
