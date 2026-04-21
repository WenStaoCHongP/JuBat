using Printf
include(joinpath(@__DIR__, "..", "..", "..", "src", "JuBat.jl"))
using .JuBat

param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
param_dim.cell.v_h = 4.2

opt = JuBat.Option()
opt.Current = x -> 5.0
opt.model = "SPMe"
opt.Nn = 10
opt.Ns = 5
opt.Np = 10
opt.Nrn = 10
opt.Nrp = 10
opt.gsorder = 2
opt.dimension = 1
opt.mechanicalmodel = "none"
opt.time = [0.0, 3600.0]
opt.dt = [0.5, 10.0]
opt.dtType = "auto"
opt.jacobi = "update"
opt.solveType = "Crank-Nicolson"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.cool_method = "surface"
opt.per_element_spme = true
opt.debug_coupling = false
opt.czm_enabled = true
opt.czm_tol = 3e-4

case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=16, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)

result = JuBat.Solve(case)
println("solve finished")
println("final D_max = ", result["czm D_max"][end])
println("final D_mean = ", result["czm D_mean"][end])
println("final n_fractured = ", result["czm n_fractured"][end])
