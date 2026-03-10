SCRIPT_DIR = @__DIR__
project_root = joinpath(SCRIPT_DIR, "..")
include(joinpath(project_root, "src", "JuBat.jl"))
println("JuBat loaded ok")
