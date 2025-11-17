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


println("=== inscribed ===")
mesh_ins = jellyroll_Q4_mesh(param_dim; nx=40, ny=40, gsorder=2, crop_to_annulus=true, crop_mode=:inscribed)
summarize_mesh(mesh_ins)
fw = jellyroll_get_layer_weights(mesh_ins)
println("layer_weights present? ", fw !== nothing)
if fw !== nothing
    println("fw shape: ", size(fw))
    println("fw row sample: ", fw[1:min(5,size(fw,1)), :])
end

println("\n=== center ===")
mesh_cen = jellyroll_Q4_mesh(param_dim; nx=40, ny=40, gsorder=2, crop_to_annulus=true, crop_mode=:center)
summarize_mesh(mesh_cen)
fw2 = jellyroll_get_layer_weights(mesh_cen)
println("layer_weights present? ", fw2 !== nothing)
if fw2 !== nothing
    println("fw shape: ", size(fw2))
    println("fw row sample: ", fw2[1:min(5,size(fw2,1)), :])
else
    # compute sample weights
    fks = jellyroll_element_layer_weights(mesh_cen, param_dim; nsamples_per_dim=4, logic=:spiral)
    println("computed fks shape: ", size(fks))
    println("fks sample: ", fks[1:min(5,size(fks,1)), :])
end

println("\n=== collector_seeded ===")
mesh_col = jellyroll_Q4_mesh(param_dim; nx=120, ny=40, gsorder=2, crop_to_annulus=true, crop_mode=:collector_seeded)
summarize_mesh(mesh_col)
fw3 = jellyroll_get_layer_weights(mesh_col)
println("layer_weights present? ", fw3 !== nothing)
if fw3 !== nothing
    println("fw shape: ", size(fw3))
    println("unique rows in fw (first 5): ")
    uniq = unique(fw3, dims=1)
    println(uniq[1:min(5,size(uniq,1)), :])
end

println("done")
