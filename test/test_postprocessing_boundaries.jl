using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_DIR = joinpath(PROJECT_ROOT, "src")
const INDEX_DIR = joinpath(PROJECT_ROOT, "md", "源码函数索引")

include(joinpath(SRC_DIR, "JuBat.jl"))
using .JuBat

@testset "post-processing responsibilities and naming" begin
    post_source = read(joinpath(SRC_DIR, "PostProcessing.jl"), String)
    cycle_post_source = read(joinpath(SRC_DIR, "CyclePostProcess.jl"), String)
    cycle_data_source = read(joinpath(SRC_DIR, "CycleData.jl"), String)
    csv_source = read(joinpath(SRC_DIR, "CsvExport.jl"), String)
    czm_post_source = read(joinpath(SRC_DIR, "CzmPostProcess.jl"), String)
    entry_source = read(joinpath(SRC_DIR, "JuBat.jl"), String)

    @test length(readlines(joinpath(SRC_DIR, "PostProcessing.jl"))) == 117
    @test count("function ", post_source) == 1
    @test occursin("function PostProcessing", post_source)
    @test !occursin("plot_cycling_results", post_source)

    @test occursin("function postprocess_phase_result", cycle_post_source)
    @test occursin("function postprocess_cycle_result!", cycle_post_source)
    @test occursin("function plot_cycling_results", cycle_post_source)
    @test !occursin(r"(?m)^function\s+_", cycle_post_source)
    @test JuBat.phase_termination_symbol(JuBat.PHASE_REST, "voltage_cutoff_low") == :time
    @test JuBat.phase_termination_symbol(JuBat.PHASE_CHARGE, "voltage_cutoff_high") == :voltage
    @test JuBat.phase_termination_symbol(JuBat.PHASE_DISCHARGE, "time_limit") == :time

    @test !occursin("function export_cycle_data_to_csv", cycle_data_source)
    @test occursin("function export_cycle_data_to_csv", csv_source)
    @test occursin("function export_cycling_csv", csv_source)
    @test occursin("function get_damage_statistics", czm_post_source)

    for source in (post_source, cycle_post_source, cycle_data_source, csv_source, czm_post_source)
        @test !occursin(r"(?m)^function\s+_", source)
        @test !occursin(r"(?m)^_[A-Za-z][A-Za-z0-9_!]*\(", source)
    end

    solver_pos = findfirst("include(\"CycleSolver.jl\")", entry_source)
    post_pos = findfirst("include(\"CyclePostProcess.jl\")", entry_source)
    data_pos = findfirst("include(\"CycleData.jl\")", entry_source)
    csv_pos = findfirst("include(\"CsvExport.jl\")", entry_source)
    @test solver_pos !== nothing
    @test post_pos !== nothing
    @test data_pos !== nothing
    @test csv_pos !== nothing
    @test first(solver_pos) < first(post_pos) < first(data_pos) < first(csv_pos)

    @test isfile(joinpath(INDEX_DIR, "CyclePostProcess.md"))
end
