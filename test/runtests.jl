using Test
using JuBat

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const TEST_FILES = sort(filter(
    path -> endswith(path, ".jl") && basename(path) != "runtests.jl",
    readdir(@__DIR__; join=true),
))

@testset "package entry point" begin
    @test JuBat isa Module
end

@testset "isolated test files" begin
    for path in TEST_FILES
        @testset "$(basename(path))" begin
            command = `$(Base.julia_cmd()) --startup-file=no --project=$(PROJECT_ROOT) $(path)`
            command = addenv(command, "GKSwstype" => "100", "JULIA_NUM_THREADS" => "1")
            process = run(ignorestatus(command))
            @test process.exitcode == 0
        end
    end
end
