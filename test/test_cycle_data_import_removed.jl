using Test

include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

@testset "external cycle-data solve path is removed" begin
    cycle_source = read(joinpath(@__DIR__, "..", "src", "CycleData.jl"), String)
    entry_source = read(joinpath(@__DIR__, "..", "src", "JuBat.jl"), String)
    export_example = read(joinpath(@__DIR__, "..", "example", "循环验证", "export_cycle_data_example.jl"), String)

    @test !isdefined(JuBat, :load_cycle_data_from_csv)
    @test !occursin("load_cycle_data_from_csv", cycle_source)
    @test !occursin("load_cycle_data_from_csv", entry_source)
    @test !isfile(joinpath(@__DIR__, "..", "example", "循环验证", "czm_from_precomputed_example.jl"))
    @test !occursin("czm_from_precomputed_example.jl", export_example)

@test isdefined(JuBat, :solve_phase_with_export)
@test isdefined(JuBat, :solve_cycle_with_export)
@test isdefined(JuBat, :export_cycle_data_to_csv)

@testset "cycle data remains export-only" begin
    step = JuBat.TimeStepData(
        0.0,
        JuBat.PHASE_REST,
        3.7,
        0.0,
        [298.15, 298.15],
        298.15,
        298.15,
        [0.5],
        [0.5],
        0.5,
    )
    data = JuBat.CycleExportData(
        1,
        [step],
        [0.0 0.0; 1.0 0.0],
        reshape([1, 2, 2, 1], 1, 4),
        1,
        2,
    )

    mktempdir() do output_dir
        files = JuBat.export_cycle_data_to_csv(data, output_dir; prefix="smoke")
        @test length(files) == 6
        @test all(isfile, files)
    end
end
end
