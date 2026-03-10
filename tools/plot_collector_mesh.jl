using Plots, Statistics
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
param=JuBat.ChooseCell("Jellyroll")
mesh=JuBat.jellyroll_collector_seed_mesh(param; nθ=8, gsorder=2)
# plot nodes and element edges
x = mesh.node[:,1]; y = mesh.node[:,2]
plt = plot(; aspect_ratio=1, legend=false, title="collector-seeded mesh (nθ=20)")
# draw all element edges
ne = size(mesh.element,1)
for e in 1:ne
    nds = mesh.element[e, :]
    xy = mesh.node[nds, :]
    # close polygon
    xs = [xy[1,1], xy[2,1], xy[3,1], xy[4,1], xy[1,1]]
    ys = [xy[1,2], xy[2,2], xy[3,2], xy[4,2], xy[1,2]]
    plot!(plt, xs, ys, lw=0.6, color=:black)
end
# scatter nodes small
scatter!(plt, x, y, ms=3, color=:blue)
# annotate first 20 elements with their id and node indices
for e in 1:min(20, ne)
    nds = mesh.element[e, :]
    xy = mesh.node[nds, :]
    # center
    cx = mean(xy[:,1]); cy = mean(xy[:,2])
    annotate!(plt, cx, cy, text(string(e), 8, :red))
    # annotate nodes near element
    for (i, n) in enumerate(nds)
        nx = mesh.node[n,1]; ny = mesh.node[n,2]
        annotate!(plt, nx, ny, text(string(n), 6, :green))
    end
end
savefig(plt, joinpath(@__DIR__, "collector_mesh.png"))
println("Saved: tools/collector_mesh.png")
