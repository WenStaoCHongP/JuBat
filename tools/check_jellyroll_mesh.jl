using Random
include("../src/Jellyrollmodel.jl")
include("../src/SetMesh.jl")

function summarize_mesh(mesh)
    println("mesh type: ", mesh.type)
    println("nnode: ", mesh.nlen, ", nelem: ", size(mesh.element,1))
    # element center stats
    try
        centers = jellyroll_element_centers(mesh)
        println("element centers sample (first 3): ", centers[1:min(3,size(centers,1)), :])
    catch err
    end
end

# Build a minimal param_dim-like object with the needed fields for jellyroll_spiral_params
mutable struct LayerStub
    thickness::Float64
    lambda::Float64
end
mutable struct CellStub
    Rin::Float64
    Rout::Float64
    height::Float64
end
mutable struct ParamDimStub
    PCC::LayerStub
    PE::LayerStub
    SP::LayerStub
    NE::LayerStub
    NCC::LayerStub
    cell::CellStub
end

param_dim = ParamDimStub(
    LayerStub(1e-4, 200.0),
    LayerStub(2e-5, 0.2),
    LayerStub(1e-5, 0.1),
    LayerStub(2e-5, 0.2),
    LayerStub(1e-4, 400.0),
    CellStub(0.01, 0.05, 0.065)
)

println("=== collector_seeded ===")
mesh_col = jellyroll_collector_seed_mesh(param_dim; nθ=120, gsorder=2)
summarize_mesh(mesh_col)
fw = jellyroll_get_layer_weights(mesh_col)
println("layer_weights present? ", fw !== nothing)
if fw !== nothing
    println("fw shape: ", size(fw))
    println("unique rows in fw (first 5): ")
    uniq = unique(fw, dims=1)
    println(uniq[1:min(5,size(uniq,1)), :])
end

println("\ndone")
