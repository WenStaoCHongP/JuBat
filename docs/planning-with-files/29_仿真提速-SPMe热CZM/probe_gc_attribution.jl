# 短仿真实验：区分 CZM 真实计算 vs GC 归因偏差
include(joinpath(@__DIR__, "../../../src/JuBat.jl"))
using .JuBat
using Printf

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
opt.time = [0.0, 20]
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
opt.czm_tol = 1e-3

case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, czm_enabled=true, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

GC.gc()
println("=== 短仿真（20s 仿真时间）===")
result = @time JuBat.Solve(case)
for k in ("timing SPMe solve total [s]", "timing branch solver total [s]",
          "timing thermal distributed total [s]", "timing CZM model total [s]",
          "timing CZM model avg [ms]", "timing CallModel calls")
    @printf("%-42s %s\n", k, result[k])
end
