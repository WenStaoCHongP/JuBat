using Plots, Statistics
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "plot_collector_mesh")
mkpath(OUTPUT_DIR)
param=JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.czm_enabled = false
case = JuBat.SetCase(param, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=8, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
mesh = case.mesh["thermal2D"]
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
savefig(plt, joinpath(OUTPUT_DIR, "collector_mesh.png"))
println("Saved: ", joinpath(OUTPUT_DIR, "collector_mesh.png"))
