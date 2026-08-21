using Plots, LinearAlgebra, Statistics
include("../src/JuBat.jl")

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "check_boundary_nodes")
mkpath(OUTPUT_DIR)

# Simple diagnostic: mark tab nodes (red) and outer boundary nodes (green)
function main()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.czm_enabled = false
    case = JuBat.SetCase(param_dim, opt)

    # build mesh like examples
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=180, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # obtain tab node indices (use JuBat module-qualified name)
    pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
    pos_nodes = unique(pos_idx)
    neg_nodes = unique(neg_idx)
    tab_nodes = unique(vcat(pos_nodes, neg_nodes))

    # identify outer boundary nodes using edge_boundary logic
    # which=:outer 表示外螺旋
    # 计算 theta_range 参数
    Rin, Rout = param_dim.cell.Rin, param_dim.cell.Rout
    bval = max(param_dim.cell.layer / (2 * pi), 1e-12)
    theta0_mesh = max(0.0, (Rin - param_dim.cell.Rin) / bval)
    theta1_mesh = min((Rout - param_dim.cell.Rin - param_dim.cell.layer) / bval, (Rout - param_dim.cell.Rin) / bval)
    theta_in_range = (theta0_mesh, min(theta0_mesh + 2.0 * pi, theta1_mesh))
    theta_out_range = (max(theta1_mesh - 2.0 * pi, theta0_mesh), theta1_mesh)

    inner_nodes = Int[]
    for i in 1:mesh_th.nlen
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:inner, theta_range=theta_in_range)
            push!(inner_nodes, i)
        end
    end
    outer_nodes = Int[]
    for i in 1:mesh_th.nlen
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:outer, theta_range=theta_out_range)
            push!(outer_nodes, i)
        end
    end
    x = mesh_th.node[:,1]; y = mesh_th.node[:,2]

    println("Found tab nodes: ", length(tab_nodes), ", inner boundary nodes: ", length(inner_nodes), ", outer boundary nodes: ", length(outer_nodes))

    # plot
    plt = plot(size=(900,900), title="Boundary nodes: red=tab, green=outer, blue=inner")
    scatter!(plt, x, y; ms=1.0, color=:gray, alpha=0.3, label=false)
    if !isempty(inner_nodes)
        scatter!(plt, x[inner_nodes], y[inner_nodes]; ms=3.5, color=:blue, label="inner boundary")
    end
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
    savefig(plt, joinpath(OUTPUT_DIR, "boundary_nodes.png"))
    println("Saved: ", joinpath(OUTPUT_DIR, "boundary_nodes.png"))
    # also save vector SVG and a high-resolution bitmap to avoid pixelation when zooming
    try
        savefig(plt, joinpath(OUTPUT_DIR, "boundary_nodes.svg"))
        println("Saved: ", joinpath(OUTPUT_DIR, "boundary_nodes.svg"))
    catch e
        @warn "Failed to save SVG" exception=(e, catch_backtrace())
    end
    try
        # create a high-resolution raster version (e.g., 3600x3600)
    plt_hr = plot(size=(3600,3600), title="Boundary nodes: red=tab, green=outer, blue=inner")
        scatter!(plt_hr, x, y; ms=1.0, color=:gray, alpha=0.3, label=false)
        if !isempty(inner_nodes)
            scatter!(plt_hr, x[inner_nodes], y[inner_nodes]; ms=6.0, color=:blue, label="inner boundary")
        end
        if !isempty(outer_nodes)
            scatter!(plt_hr, x[outer_nodes], y[outer_nodes]; ms=6.0, color=:green, label="outer boundary")
        end
        if !isempty(pos_nodes)
            scatter!(plt_hr, x[pos_nodes], y[pos_nodes]; ms=8.0, color=:red, label="pos tabs")
        end
        if !isempty(neg_nodes)
            scatter!(plt_hr, x[neg_nodes], y[neg_nodes]; ms=8.0, color=:yellow, label="neg tabs")
        end
        savefig(plt_hr, joinpath(OUTPUT_DIR, "boundary_nodes_hr.png"))
        println("Saved: ", joinpath(OUTPUT_DIR, "boundary_nodes_hr.png"))
    catch e
        @warn "Failed to save high-res PNG" exception=(e, catch_backtrace())
    end
end

main()
