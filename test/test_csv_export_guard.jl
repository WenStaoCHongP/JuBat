using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))

@testset "CSV export guarded write" begin
    files_written = String[]
    files_skipped = String[]
    calls = Symbol[]

    result = JuBat.write_csv_guarded!("success.csv", files_written, files_skipped) do
        push!(calls, :success)
    end
    @test result === nothing
    @test calls == [:success]
    @test files_written == ["success.csv"]
    @test isempty(files_skipped)

    @test_logs (:warn, r"Failed to write failure.csv") JuBat.write_csv_guarded!(
        "failure.csv", files_written, files_skipped) do
        error("expected test failure")
    end
    @test files_written == ["success.csv"]
    @test files_skipped == ["failure.csv"]
end

@testset "CSV export minimal public path" begin
    result = (cycle_results=Any[], soh=Float64[], czm_snapshots=Any[])
    sim_case = (opt=(thermal_enabled=false,), geometry=nothing)

    mktempdir() do output_dir
        files_written = JuBat.export_cycling_csv(
            result, sim_case, nothing; output_dir=output_dir, overwrite=true)
        @test files_written == ["cycle_summary.csv"]
        @test isfile(joinpath(output_dir, "cycle_summary.csv"))
        @test startswith(read(joinpath(output_dir, "cycle_summary.csv"), String), "cycle,phase")
    end
end
