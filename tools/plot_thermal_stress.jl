using Plots
using Statistics
include("../src/JuBat.jl")

function main()
    # setup
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.gsorder = 2
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.czm_enabled = false
    case = JuBat.SetCase(param_dim, opt)

    # build thermal mesh
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=120, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # prepare variables and set a non-uniform T_nodes for visualization
    variables = Dict{String, Union{Array{Float64}, Float64}}()
    # create a radial temperature gradient for demonstration
    x = mesh_th.node[:,1]; y = mesh_th.node[:,2]
    r = hypot.(x, y)
    rmin = minimum(r); rmax = maximum(r)
    # nondimensional temperature nodes, small radial ramp from T0 to T0+0.1
    T0 = case.param.cell.T0
    T_nodes = fill(T0, length(r)) .+ 0.1 .* ((r .- rmin) ./ (rmax - rmin))
    variables["T_nodes"] = T_nodes
    # also provide scalar "temperature" used by thermal_stress (1D avg T)
    variables["temperature"] = mean(T_nodes)

    # compute thermal stress (per-element)
    variables = JuBat.thermal_stress(case, variables)
    σ_elem = variables["thermal2D element thermal stress"]

    # element centers
    centers = JuBat.jellyroll_element_centers(mesh_th)
    cx = centers[:,1]; cy = centers[:,2]

    # scatter plot of element stress
    plt = scatter(cx, cy, marker_z = σ_elem, ms=6, markershape=:rect, colorbar=true,
        title = "Thermal element stress (nondim)", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "thermal_element_stress.png")
    println("Saved: thermal_element_stress.png")
    try
        savefig(plt, "thermal_element_stress.svg")
        println("Saved: thermal_element_stress.svg")
    catch e
        @warn "Failed to save SVG" exception=(e, catch_backtrace())
    end
end

main()
