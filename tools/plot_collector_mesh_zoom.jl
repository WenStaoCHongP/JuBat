using Plots, Statistics
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "plot_collector_mesh_zoom")
mkpath(OUTPUT_DIR)
param=JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.czm_enabled = false
case = JuBat.SetCase(param, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=40, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
mesh = case.mesh["thermal2D"]

# compute element centers (centroid of nodes)
ne = size(mesh.element,1)
centers = zeros(ne,2)
r = zeros(ne)
for e in 1:ne
    nds = mesh.element[e, :]
    xy = mesh.node[nds, :]
    centers[e,1] = mean(xy[:,1])
    centers[e,2] = mean(xy[:,2])
    r[e] = sqrt(centers[e,1]^2 + centers[e,2]^2)
end

# select inner ring: either by radius threshold (e.g., 20th percentile) or fixed count
pct = 0.20
rth = quantile(r, pct)
inner_idx = findall(r .<= rth)
# if too many or too few, clamp
n_selected = clamp(length(inner_idx), 40, 800)
# sort inner by radius and take closest n_selected
sorted_inner = sort(inner_idx, by=i->r[i])
sel = sorted_inner[1:n_selected]

# compute bounding box for selected elements and expand margin
sel_nodes = unique(vcat([mesh.element[e,:] for e in sel]...))
xs = mesh.node[sel_nodes,1]; ys = mesh.node[sel_nodes,2]
xmin, xmax = minimum(xs), maximum(xs)
ymin, ymax = minimum(ys), maximum(ys)
padx = 0.05*(xmax-xmin); pady = 0.05*(ymax-ymin)

plt = plot(; aspect_ratio=1, legend=false, title="collector inner ring zoom (nθ=40)")
# draw selected element filled polygons with light gray border
for e in sel
    nds = mesh.element[e, :]
    xy = mesh.node[nds, :]
    xs_poly = [xy[1,1], xy[2,1], xy[3,1], xy[4,1]]
    ys_poly = [xy[1,2], xy[2,2], xy[3,2], xy[4,2]]
    plot!(plt, xs_poly, ys_poly, seriestype=:shape, linecolor=:black, fillcolor=:antiquewhite, alpha=0.9)
end

# draw node markers for nodes in selection
scatter!(plt, mesh.node[sel_nodes,1], mesh.node[sel_nodes,2], ms=4, color=:blue)

# annotate first 40 nodes
for (i,n) in enumerate(sel_nodes[1:min(40,end)])
    annotate!(plt, mesh.node[n,1], mesh.node[n,2], text(string(n), 8, :green))
end

xlims!(plt, xmin-padx, xmax+padx)
ylims!(plt, ymin-pady, ymax+pady)

png_path = joinpath(OUTPUT_DIR, "collector_mesh_zoom.png")
svg_path = joinpath(OUTPUT_DIR, "collector_mesh_zoom.svg")

savefig(plt, png_path)
println("Saved: $png_path")
# also save svg
savefig(plt, svg_path)
println("Saved: $svg_path")
