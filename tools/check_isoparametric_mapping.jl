using Statistics
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))

param = JuBat.ChooseCell("Jellyroll")
mesh = JuBat.jellyroll_collector_seed_mesh(param; nθ=40, gsorder=2)

gs = mesh.gs
# gs.xi: (ngs, dim), gs.x: (ngs, dim), gs.ele: length ngs
ngs = size(gs.x,1)
ne = size(mesh.element,1)
max_abs_err = 0.0
rms_err = 0.0
count = 0

# For speed, test a random subset if too many gauss points
indices = collect(1:ngs)
if ngs > 4000
    using Random
    Random.seed!(42)
    # use randperm to get a random subset without replacement
    indices = sort(Random.randperm(ngs)[1:4000])
end

errs = Float64[]
neg_detj = Int[]
for idx in indices
    xi = gs.xi[idx, :]
    ele = gs.ele[idx]
    nds = mesh.element[ele, :]
    nodes = mesh.node[nds, :]
    # compute N at xi using Q4 shape functions (same formula as in ShapeFunction2D / LagrangeBasis)
    ξ = xi[1]; η = xi[2]
    N = 0.25 .* [(1-ξ)*(1-η), (1+ξ)*(1-η), (1+ξ)*(1+η), (1-ξ)*(1+η)]
    x_ref = N[1]*nodes[1,:] + N[2]*nodes[2,:] + N[3]*nodes[3,:] + N[4]*nodes[4,:]
    x_gs = gs.x[idx, :]
    err = maximum(abs.(x_ref .- x_gs))
    push!(errs, err)
    # check detJ at this gauss point as provided by gs.detJ
    detj = gs.detJ[idx]
    if detj <= 0.0
        push!(neg_detj, idx)
    end
end

println("Checked $(length(indices)) gauss points (out of $ngs)")
println("max mapping abs error = ", maximum(errs))
println("median mapping abs error = ", median(errs))
println("mean mapping abs error = ", mean(errs))
println("gauss points with non-positive detJ: ", length(neg_detj))
if length(neg_detj) > 0
    println("first few indices with non-positive detJ: ", neg_detj[1:min(end,20)])
end

# report per-element detJ stats (min of gauss points per element)
ne = size(mesh.element,1)
# try to reshape detJ into (ngauss_per_ele, ne) if possible
ng_per_ele = Int(round(length(gs.detJ)/ne))
if ng_per_ele * ne == length(gs.detJ)
    detJ_mat = reshape(gs.detJ, ng_per_ele, ne)
    per_ele_min = mapslices(minimum, detJ_mat; dims=1)'
    println("per-element detJ: min=$(minimum(per_ele_min)), median=$(median(per_ele_min)), mean=$(mean(per_ele_min))")
    # list elements with very small detJ relative to median
    med = median(per_ele_min)
    small = findall(per_ele_min .< med*1e-3)
    println("elements with detJ < 1e-3 * median: $(length(small))")
    if length(small) > 0
        println("first few elements: ", small[1:min(end,20)])
    end
else
    println("Could not reshape gs.detJ into neat per-element matrix: ngs=$(length(gs.detJ)), ne=$ne")
end
