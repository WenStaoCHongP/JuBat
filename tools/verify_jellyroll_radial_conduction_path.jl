"""
验证 Jellyroll 热模型边界节点温度异常是否由缺失径向导热路径导致。
运行: julia --project=. example/热模块验证/verify_jellyroll_radial_conduction_path.jl
"""
using LinearAlgebra, Statistics, Printf
include(joinpath(@__DIR__, "..", "..", "src", "JuBat.jl"))
using .JuBat

const SEP = "="^80

function build_jellyroll_option()
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn, opt.Ns, opt.Np = 10, 5, 10
    opt.Nrn, opt.Nrp = 10, 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"
    opt.Current = _ -> 5.0
    opt.time = [0.0, 1.0]
    opt.dt = [1e-3, 1.0]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "none"
    opt.per_element_spme = true
    opt.czm_enabled = false
    return opt
end

function build_ring_option()
    opt = JuBat.Option()
    opt.model = "thermal"
    opt.thermal_enabled = true
    opt.thermalmodel = "ring2D"
    opt.cool_method = "none"
    return opt
end

function unique_edges(mesh)
    seen = Set{Tuple{Int, Int}}()
    for e in 1:size(mesh.element, 1)
        n1, n2, n3, n4 = mesh.element[e, :]
        for (a, b) in ((n1, n2), (n2, n3), (n3, n4), (n4, n1))
            key = a < b ? (a, b) : (b, a)
            push!(seen, key)
        end
    end
    return collect(seen)
end

function build_adjacency(nnode, edges)
    adj = [Int[] for _ in 1:nnode]
    for (a, b) in edges
        push!(adj[a], b)
        push!(adj[b], a)
    end
    return adj
end

function count_components(adj)
    n = length(adj)
    visited = falses(n)
    ncomp = 0
    for s in 1:n
        visited[s] && continue
        ncomp += 1
        q = [s]
        visited[s] = true
        head = 1
        while head <= length(q)
            u = q[head]
            head += 1
            for v in adj[u]
                if !visited[v]
                    visited[v] = true
                    push!(q, v)
                end
            end
        end
    end
    return ncomp
end

function degree_stats(adj)
    deg = [length(v) for v in adj]
    return (; min=minimum(deg), max=maximum(deg), mean=mean(deg), median=median(deg))
end

function quantized_key(x, y, tol)
    return (round(Int, x / tol), round(Int, y / tol))
end

function find_coincident_node_pairs(mesh; tol=1e-12)
    bins = Dict{Tuple{Int, Int}, Vector{Int}}()
    for i in 1:mesh.nlen
        x, y = mesh.node[i, 1], mesh.node[i, 2]
        key = quantized_key(x, y, tol)
        if !haskey(bins, key)
            bins[key] = Int[]
        end
        push!(bins[key], i)
    end

    pairs = Tuple{Int, Int}[]
    for idxs in values(bins)
        if length(idxs) > 1
            for i in 1:(length(idxs) - 1)
                for j in (i + 1):length(idxs)
                    a, b = idxs[i], idxs[j]
                    push!(pairs, a < b ? (a, b) : (b, a))
                end
            end
        end
    end
    return unique(pairs)
end

function bfs_distance(adj, s::Int, t::Int)
    s == t && return 0
    n = length(adj)
    dist = fill(-1, n)
    q = [s]
    dist[s] = 0
    head = 1
    while head <= length(q)
        u = q[head]
        head += 1
        du = dist[u]
        for v in adj[u]
            if dist[v] < 0
                dist[v] = du + 1
                v == t && return dist[v]
                push!(q, v)
            end
        end
    end
    return typemax(Int)
end

function sample_path_distances(adj, pairs; nsample=40)
    isempty(pairs) && return Int[]
    step = max(1, fld(length(pairs), nsample))
    sampled = pairs[1:step:end]
    out = Int[]
    for (a, b) in sampled
        d = bfs_distance(adj, a, b)
        if d != typemax(Int)
            push!(out, d)
        end
    end
    return out
end

function element_dr_dtheta_stats(mesh)
    ne = size(mesh.element, 1)
    dr = zeros(Float64, ne)
    dtheta = zeros(Float64, ne)
    for e in 1:ne
        nodes = mesh.element[e, :]
        rs = [hypot(mesh.node[n, 1], mesh.node[n, 2]) for n in nodes]
        ts = [atan(mesh.node[n, 2], mesh.node[n, 1]) for n in nodes]
        dr[e] = maximum(rs) - minimum(rs)

        ts_sorted = sort(ts)
        gaps = [ts_sorted[i + 1] - ts_sorted[i] for i in 1:(length(ts_sorted) - 1)]
        push!(gaps, (ts_sorted[1] + 2pi) - ts_sorted[end])
        dtheta[e] = 2pi - maximum(gaps)
    end
    return (; dr_min=minimum(dr), dr_max=maximum(dr), dr_mean=mean(dr),
            dtheta_min=minimum(dtheta), dtheta_max=maximum(dtheta), dtheta_mean=mean(dtheta))
end

function analyze_jellyroll(param_dim; ntheta=24, gsorder=2)
    opt = build_jellyroll_option()
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=ntheta, gsorder=gsorder)

    case_un = JuBat.setup_thermal2D_mesh(case, mesh_data; use_merged=false)
    mesh_un = case_un.mesh["thermal2D"]

    case_mg = JuBat.setup_thermal2D_mesh(case, mesh_data; use_merged=true)
    mesh_mg = case_mg.mesh["thermal2D"]

    edges_un = unique_edges(mesh_un)
    adj_un = build_adjacency(mesh_un.nlen, edges_un)

    coincident_pairs = if isempty(mesh_data.interface_pairs)
        find_coincident_node_pairs(mesh_un)
    else
        unique([(min(a, b), max(a, b)) for (a, b) in mesh_data.interface_pairs])
    end

    direct_connected = count(((a, b),) -> (b in adj_un[a]), coincident_pairs)

    vars_un = Dict{String, Union{Array{Float64}, Float64}}(
        "heat_source_fields" => zeros(Float64, size(mesh_un.element, 1))
    )
    _, KT_un, _ = JuBat.ThermalDistributed2D(case_un, vars_un)

    pair_coupling_un = Float64[]
    for (a, b) in coincident_pairs
        push!(pair_coupling_un, abs(KT_un[a, b]) + abs(KT_un[b, a]))
    end
    nonzero_pair_coupling_un = count(>(1e-14), pair_coupling_un)

    merge_map = mesh_data.merge_map
    merged_same_node = count(((a, b),) -> merge_map[a] == merge_map[b], coincident_pairs)

    vars_mg = Dict{String, Union{Array{Float64}, Float64}}(
        "heat_source_fields" => zeros(Float64, size(mesh_mg.element, 1))
    )
    _, KT_mg, _ = JuBat.ThermalDistributed2D(case_mg, vars_mg)

    mapped_distinct_pairs = Tuple{Int, Int}[]
    for (a, b) in coincident_pairs
        na, nb = merge_map[a], merge_map[b]
        na == nb && continue
        push!(mapped_distinct_pairs, na < nb ? (na, nb) : (nb, na))
    end
    unique!(mapped_distinct_pairs)

    pair_coupling_mg = Float64[]
    for (a, b) in mapped_distinct_pairs
        push!(pair_coupling_mg, abs(KT_mg[a, b]) + abs(KT_mg[b, a]))
    end

    path_dists = sample_path_distances(adj_un, coincident_pairs; nsample=40)

    return (
        mesh_un=mesh_un,
        mesh_mg=mesh_mg,
        n_comp_un=count_components(adj_un),
        deg_un=degree_stats(adj_un),
        coincident_pairs=coincident_pairs,
        direct_connected=direct_connected,
        nonzero_pair_coupling_un=nonzero_pair_coupling_un,
        merged_same_node=merged_same_node,
        mapped_distinct_pairs=mapped_distinct_pairs,
        pair_coupling_mg=pair_coupling_mg,
        path_dists=path_dists,
        drdtheta_un=element_dr_dtheta_stats(mesh_un),
        drdtheta_mg=element_dr_dtheta_stats(mesh_mg),
    )
end

function analyze_ring(param_dim; ntheta=24)
    param_ring = JuBat.ChooseCell("Ring")

    # Keep geometry and anisotropic conductivity consistent with Jellyroll base cell.
    for name in (:Rin, :Rout, :lambda_r, :lambda_t, :h, :T0, :T_amb)
        setproperty!(param_ring.cell, name, getproperty(param_dim.cell, name))
    end

    opt = build_ring_option()
    case = JuBat.SetCase(param_ring, opt)

    # Approximate number of radial layers from Jellyroll repeat thickness.
    nr_est = round(Int, (param_ring.cell.Rout - param_ring.cell.Rin) / max(param_ring.cell.layer, 1e-12))
    nr = clamp(nr_est, 4, 40)
    mesh_data = JuBat.ring_mesh(case.param; ntheta=ntheta, nr=nr, gsorder=2)
    mesh = mesh_data.mesh
    case.mesh["thermal2D"] = mesh

    edges = unique_edges(mesh)
    adj = build_adjacency(mesh.nlen, edges)

    vars = Dict{String, Any}("heat_source_fields" => zeros(Float64, size(mesh.element, 1)))
    _, KT, _ = JuBat.ThermalDistributed2D_Ring(case, vars)

    # Radial neighbor pairs: same theta index, adjacent radial index.
    idx(ir, it) = (ir - 1) * mesh_data.ntheta + it
    radial_pairs = Tuple{Int, Int}[]
    for ir in 1:mesh_data.nr
        for it in 1:mesh_data.ntheta
            a = idx(ir, it)
            b = idx(ir + 1, it)
            push!(radial_pairs, a < b ? (a, b) : (b, a))
        end
    end

    direct_connected = count(((a, b),) -> (b in adj[a]), radial_pairs)
    radial_coupling = [abs(KT[a, b]) + abs(KT[b, a]) for (a, b) in radial_pairs]
    nonzero_radial_coupling = count(>(1e-14), radial_coupling)

    return (
        mesh=mesh,
        n_comp=count_components(adj),
        deg=degree_stats(adj),
        radial_pairs=radial_pairs,
        direct_connected=direct_connected,
        nonzero_radial_coupling=nonzero_radial_coupling,
        drdtheta=element_dr_dtheta_stats(mesh),
        nr=mesh_data.nr,
        ntheta=mesh_data.ntheta,
    )
end

function print_summary(jelly, ring)
    println(SEP)
    println("Jellyroll vs Ring: radial conduction path diagnostics")
    println(SEP)

    np = length(jelly.coincident_pairs)
    frac_direct = np == 0 ? 0.0 : jelly.direct_connected / np
    frac_nz_un = np == 0 ? 0.0 : jelly.nonzero_pair_coupling_un / np
    frac_merged_same = np == 0 ? 0.0 : jelly.merged_same_node / np

    println("[Jellyroll default (unmerged)]")
    @printf("  nodes=%d, elements=%d, connected_components=%d\n", jelly.mesh_un.nlen, size(jelly.mesh_un.element, 1), jelly.n_comp_un)
    @printf("  degree(min/median/mean/max)=%.0f / %.1f / %.2f / %.0f\n", jelly.deg_un.min, jelly.deg_un.median, jelly.deg_un.mean, jelly.deg_un.max)
    @printf("  coincident inter-layer node pairs=%d\n", np)
    @printf("  directly connected coincident pairs=%d (%.2f%%)\n", jelly.direct_connected, 100 * frac_direct)
    @printf("  nonzero K_ij on coincident pairs=%d (%.2f%%)\n", jelly.nonzero_pair_coupling_un, 100 * frac_nz_un)
    if !isempty(jelly.path_dists)
        @printf("  graph distance between coincident pairs (sample): min=%d, median=%.1f, max=%d\n",
                minimum(jelly.path_dists), median(jelly.path_dists), maximum(jelly.path_dists))
    end

    println("\n[Jellyroll merged connectivity reference]")
    @printf("  nodes=%d, elements=%d\n", jelly.mesh_mg.nlen, size(jelly.mesh_mg.element, 1))
    @printf("  coincident pairs collapsed to same node=%d (%.2f%%)\n", jelly.merged_same_node, 100 * frac_merged_same)

    println("\n[Ring mesh reference]")
    nrp = length(ring.radial_pairs)
    frac_ring_direct = nrp == 0 ? 0.0 : ring.direct_connected / nrp
    frac_ring_nz = nrp == 0 ? 0.0 : ring.nonzero_radial_coupling / nrp
    @printf("  nodes=%d, elements=%d, nr=%d, ntheta=%d, connected_components=%d\n",
            ring.mesh.nlen, size(ring.mesh.element, 1), ring.nr, ring.ntheta, ring.n_comp)
    @printf("  radial neighbor pairs=%d\n", nrp)
    @printf("  directly connected radial pairs=%d (%.2f%%)\n", ring.direct_connected, 100 * frac_ring_direct)
    @printf("  nonzero K_ij on radial pairs=%d (%.2f%%)\n", ring.nonzero_radial_coupling, 100 * frac_ring_nz)

    println("\n[Mesh partition comparison]")
    j = jelly.drdtheta_un
    r = ring.drdtheta
    @printf("  Jellyroll element dr(mean)=%.4e, dtheta(mean)=%.4e\n", j.dr_mean, j.dtheta_mean)
    @printf("  Ring      element dr(mean)=%.4e, dtheta(mean)=%.4e\n", r.dr_mean, r.dtheta_mean)

    println("\n[Diagnosis]")
    likely_missing_radial_path = (frac_nz_un < 0.05) && (frac_merged_same > 0.95) && (frac_ring_nz > 0.95)
    if likely_missing_radial_path
        println("  RESULT: Evidence strongly supports that direct radial conduction paths are missing in default Jellyroll unmerged thermal mesh.")
        println("          This can produce boundary-node temperature anomalies compared with ring-like radial connectivity.")
    else
        println("  RESULT: Current evidence is not sufficient to conclude missing radial conduction is the dominant cause.")
        println("          Please increase ntheta and inspect pair-coupling/path-distance outputs.")
    end

    println(SEP)
end

function main()
    param_dim = JuBat.ChooseCell("Jellyroll")

    # Use a moderate ntheta so this diagnostics script stays lightweight.
    jelly = analyze_jellyroll(param_dim; ntheta=24, gsorder=2)
    ring = analyze_ring(param_dim; ntheta=24)

    print_summary(jelly, ring)
end

main()
