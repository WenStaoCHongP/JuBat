using LinearAlgebra, Statistics
using Printf
include("../src/JuBat.jl")

function main()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    case = JuBat.SetCase(param_dim, opt)

    # build mesh like examples
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)
    case.mesh["thermal2D"] = mesh_th

    # obtain tab node indices
    pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))

    # identify boundary edges/nodes using same logic as check_boundary_nodes.jl
    ne = size(mesh_th.element, 1)
    counts = Dict{Tuple{Int,Int}, Int}()
    edge_meta = Dict{Tuple{Int,Int}, Tuple{Int,Int,Int,Int}}()
    for e in 1:ne
        n1, n2, n3, n4 = mesh_th.element[e,1], mesh_th.element[e,2], mesh_th.element[e,3], mesh_th.element[e,4]
        local_edges = ((1,2,n1,n2), (2,3,n2,n3), (3,4,n3,n4), (4,1,n4,n1))
        for (i1,i2,a,b) in local_edges
            key = a < b ? (a,b) : (b,a)
            counts[key] = get(counts, key, 0) + 1
            if !haskey(edge_meta, key)
                edge_meta[key] = (e, (i1==1 && i2==2) ? 1 : (i1==2 && i2==3) ? 2 : (i1==3 && i2==4) ? 3 : 4, a, b)
            end
        end
    end

    x = mesh_th.node[:,1]; y = mesh_th.node[:,2]
    boundary_nodes = Int[]
    rmeans = Float64[]
    edge_lengths = Float64[]
    for (key, c) in counts
        c == 1 || continue
        (_, _, a, b) = edge_meta[key]
        push!(boundary_nodes, a); push!(boundary_nodes, b)
        ra = hypot(x[a], y[a]); rb = hypot(x[b], y[b])
        push!(rmeans, 0.5*(ra + rb))
        push!(edge_lengths, hypot(x[b]-x[a], y[b]-y[a]))
    end
    boundary_nodes = unique(boundary_nodes)

    # Heuristic: pick outer boundary nodes by radius proximity to max r
    rvals = hypot.(x, y)
    unique_radii = Float64[]
    if !isempty(rmeans)
        rmax = maximum(rmeans)
        L_med = isempty(edge_lengths) ? 0.0 : median(edge_lengths)
        tol = max(1e-6, 0.6 * L_med)
        outer_nodes = Int[]
        for n in boundary_nodes
            if (rmax - rvals[n]) <= tol
                push!(outer_nodes, n)
            end
        end
        outer_nodes = unique(outer_nodes)
        unique_radii = sort(unique(rvals[outer_nodes]))
    else
        outer_nodes = Int[]
    end

    println("Found tab nodes: ", length(tab_nodes), ", outer boundary nodes: ", length(outer_nodes))
    println("Unique outer radii (count=$(length(unique_radii))):")
    for rr in unique_radii
        @printf("  %0.12f\n", rr)
    end

    # print some spiral parameters
    try
        sp = JuBat.jellyroll_spiral_params(param_dim)
        println("Spiral params:")
        for (k,v) in sp
            println("  ", k, " => ", v)
        end
    catch
        println("jellyroll_spiral_params not available in this build")
    end
end

main()
