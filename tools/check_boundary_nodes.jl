using Plots, LinearAlgebra, Statistics
include("../src/JuBat.jl")

# Simple diagnostic: mark tab nodes (red) and outer boundary nodes (green)
function main()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    case = JuBat.SetCase(param_dim, opt)

    # build mesh like examples
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)
    case.mesh["thermal2D"] = mesh_th

    # obtain tab node indices (use JuBat module-qualified name)
    pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
    pos_nodes = unique(pos_idx)
    neg_nodes = unique(neg_idx)
    tab_nodes = unique(vcat(pos_nodes, neg_nodes))

    # identify outer boundary nodes using new edge_boundary logic
    # which=:outer 表示外螺旋终圈 (θ_cum ∈ [2π(N-1),2πN])
    outer_nodes = Int[]
    for i in 1:mesh_th.nlen
        if JuBat.edge_boundary(:node_on, mesh_th, i, param_dim; which=:outer)
            push!(outer_nodes, i)
        end
    end
    x = mesh_th.node[:,1]; y = mesh_th.node[:,2]

    println("Found tab nodes: ", length(tab_nodes), ", outer boundary nodes: ", length(outer_nodes))

    # plot
    plt = plot(size=(900,900), title="Boundary nodes (edge_boundary): red=tab, green=outer boundary")
    scatter!(plt, x, y; ms=1.0, color=:gray, alpha=0.3, label=false)
    if !isempty(outer_nodes)
        scatter!(plt, x[outer_nodes], y[outer_nodes]; ms=3.5, color=:green, label="outer boundary")
    end
    if !isempty(pos_nodes)
        scatter!(plt, x[pos_nodes], y[pos_nodes]; ms=4.0, color=:red, label="pos tabs")
    end
    if !isempty(neg_nodes)
        scatter!(plt, x[neg_nodes], y[neg_nodes]; ms=4.0, color=:yellow, label="neg tabs")
    end
    # save regular PNG
    savefig(plt, "boundary_nodes.png")
    println("Saved: boundary_nodes.png")
    # also save vector SVG and a high-resolution bitmap to avoid pixelation when zooming
    try
        savefig(plt, "boundary_nodes.svg")
        println("Saved: boundary_nodes.svg")
    catch e
        @warn "Failed to save SVG" exception=(e, catch_backtrace())
    end
    try
        # create a high-resolution raster version (e.g., 3600x3600)
    plt_hr = plot(size=(3600,3600), title="Boundary nodes (edge_boundary): red=tab, green=outer boundary")
        scatter!(plt_hr, x, y; ms=1.0, color=:gray, alpha=0.3, label=false)
        if !isempty(outer_nodes)
            scatter!(plt_hr, x[outer_nodes], y[outer_nodes]; ms=6.0, color=:green, label="outer boundary")
        end
        if !isempty(pos_nodes)
            scatter!(plt_hr, x[pos_nodes], y[pos_nodes]; ms=8.0, color=:red, label="pos tabs")
        end
        if !isempty(neg_nodes)
            scatter!(plt_hr, x[neg_nodes], y[neg_nodes]; ms=8.0, color=:yellow, label="neg tabs")
        end
        savefig(plt_hr, "boundary_nodes_hr.png")
        println("Saved: boundary_nodes_hr.png")
    catch e
        @warn "Failed to save high-res PNG" exception=(e, catch_backtrace())
    end
end

main()
