include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Statistics

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
opt.time = [0.0, 30.0]
opt.dt = [0.01, 0.5]
opt.dtType = "auto"
opt.jacobi = "update"
opt.solveType = "Crank-Nicolson"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.per_element_spme = true
opt.collector_seeded = true

case = JuBat.SetCase(param_dim, opt)
mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=16, gsorder=2)
case.mesh["thermal2D"] = mesh_th

y0 = JuBat.ModelInitialisation_MultiSPMe(case)

println("layout", case.multi_spme_layout)
println("initial T0 =", case.param.cell.T0)
println("t0 scale =", case.param.scale.t0)
println("dt_min nd =", case.opt.dt[1] / case.param.scale.t0)
println("dt_max nd =", case.opt.dt[2] / case.param.scale.t0)
println("Biot number =", case.param_dim.scale.h_th)
println("ambient T_nd =", case.param_dim.cell.T_amb / case.param_dim.scale.T_ref)

M, K, F, vars, y_phi = JuBat.CallModel_MultiSPMe(case, y0, 0.0, jacobi="update")

println("mean voltage nd =", vars["thermal2D common voltage"])
if haskey(vars, "T_nodes")
	T_nodes = vars["T_nodes"]
	println("CallModel_MultiSPMe returned T_nodes stats: min=$(minimum(T_nodes)), max=$(maximum(T_nodes)), mean=$(mean(T_nodes)))")
end

dt_init = 1e-8
vc = 1:size(M, 1)
y_c = (M - K * dt_init) \ (M * y0[vc] + F * dt_init)
y_old = vcat(y_c, y_phi)

M2, K2, F2, vars2, y_phi2 = JuBat.CallModel(case, y_old, dt_init, jacobi="update")

println("post-step voltage nd =", vars2["thermal2D common voltage"])
if haskey(vars2, "T_nodes")
	T_nodes2 = vars2["T_nodes"]
	println("CallModel post-step T_nodes stats: min=$(minimum(T_nodes2)), max=$(maximum(T_nodes2)), mean=$(mean(T_nodes2)))")
end
if haskey(vars2, "heat_source_fields")
	q2 = vars2["heat_source_fields"]
	println("post-step heat_source min=$(minimum(q2)), max=$(maximum(q2)), mean=$(mean(q2)))")
else
	println("post-step heat_source not available")
end

t_step = case.opt.dt[1] / case.param.scale.t0
M3, K3, F3, vars3, y_phi3 = JuBat.CallModel(case, y_old, t_step, jacobi="update")
println("CallModel at t=dt_min voltage nd =", vars3["thermal2D common voltage"])
if haskey(vars3, "heat_source_fields")
	q3 = vars3["heat_source_fields"]
	println("dt_min heat_source min=$(minimum(q3)), max=$(maximum(q3)), mean=$(mean(q3)))")
else
	println("dt_min heat_source not available")
end