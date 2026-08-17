using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_DIR = joinpath(PROJECT_ROOT, "src")
const INDEX_DIR = joinpath(PROJECT_ROOT, "md", "源码函数索引")

@testset "CZM source is split by responsibility" begin
    src_names = readdir(SRC_DIR)
    @test "CzmMesh.jl" in src_names
    @test "Czm.jl" in src_names
    @test "CzmBC.jl" in src_names
    @test !("czm.jl" in src_names)

    mesh_source = read(joinpath(SRC_DIR, "CzmMesh.jl"), String)
    core_source = read(joinpath(SRC_DIR, "Czm.jl"), String)
    bc_source = read(joinpath(SRC_DIR, "CzmBC.jl"), String)
    entry_source = read(joinpath(SRC_DIR, "JuBat.jl"), String)

    @test occursin("mutable struct CohesiveElement", mesh_source)
    @test occursin("function create_czm_mesh", mesh_source)
    @test !occursin("function assemble_czm_system", mesh_source)
    @test !occursin("function apply_bc_czm", mesh_source)

    @test occursin("mutable struct DamageState", core_source)
    @test occursin("function assemble_czm_system", core_source)
    @test occursin("function build_czm_cache", core_source)
    @test !occursin("function create_czm_mesh", core_source)
    @test !occursin("function apply_bc_czm", core_source)

    @test occursin("function apply_bc_czm", bc_source)
    @test occursin("function identify_bc_nodes_czm", bc_source)
    @test !occursin("function assemble_czm_system", bc_source)

    bc_pos = findfirst("include(\"CzmBC.jl\")", entry_source)
    core_pos = findfirst("include(\"Czm.jl\")", entry_source)
    mesh_pos = findfirst("include(\"CzmMesh.jl\")", entry_source)
    @test bc_pos !== nothing
    @test core_pos !== nothing
    @test mesh_pos !== nothing
    @test first(bc_pos) < first(core_pos) < first(mesh_pos)
    @test !occursin("include(\"czm.jl\")", entry_source)

    index_names = readdir(INDEX_DIR)
    @test "CzmMesh.md" in index_names
    @test "Czm.md" in index_names
    @test "CzmBC.md" in index_names
    @test !("czm.md" in index_names)
end
