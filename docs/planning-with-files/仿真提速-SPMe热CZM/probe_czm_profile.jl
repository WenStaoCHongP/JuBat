# CZM 单次更新 profile：定位 10.4s/次的真实热点
include(joinpath(@__DIR__, "../../../src/JuBat.jl"))
using .JuBat
using Profile

param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
param_dim.cell.v_h = 4.2
opt = JuBat.Option()
opt.Current = x -> 5.0
opt.model = "SPMe"
opt.Nn = 10; opt.Ns = 5; opt.Np = 10
opt.Nrn = 10; opt.Nrp = 10
opt.gsorder = 2
opt.dimension = 1
opt.mechanicalmodel = "none"
opt.time = [0.0, 60]
opt.dt = [0.5, 10]
opt.dtType = "auto"
opt.jacobi = "update"
opt.solveType = "Crank-Nicolson"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.cool_method = "surface"
opt.per_element_spme = true
opt.czm_enabled = true
opt.czm_fix_inner = false
opt.czm_iter_method = "basic"
opt.czm_load_steps = 10
opt.czm_tol = 1e-3

case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, czm_enabled=true, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
case.czm_layout = JuBat.CzmLayout(case.czm_mesh)

t = 0.0
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
_, _, _, vars1, _ = JuBat.CallModel(case, y0, t; jacobi="update")
T_nodes = JuBat.get_thermal_dofs(y0, case.layout)

# 预热（编译+缓存首建），然后计时+profile 第二次
JuBat.update_czm_damage!(case, vars1, T_nodes)
t_ns = time_ns()
Profile.clear()
@profile JuBat.update_czm_damage!(case, vars1, T_nodes)
elapsed = (time_ns() - t_ns) * 1e-9
println("单次 update_czm_damage! 耗时 = ", round(elapsed, digits=2), " s")
println("czm ndof = ", 2 * case.czm_mesh.nnode, ", n_coh = ", case.czm_mesh.n_cohesive)
println("="^70)
Profile.print(fmt = :flat, noisefromat = :none, mincount = 200)
