# 端到端冒烟：SPMe + 2D 热 + CZM 全耦合（缩短版），验证重设计 v2 下求解收敛
include(joinpath(@__DIR__, "../src/JuBat.jl"))
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
opt.time = [0.0, 300.0]
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
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=60, czm_enabled=true, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

result = JuBat.Solve(case)

t = result["time [s]"]
V = result["cell voltage [V]"]
@printf("steps=%d, V: %.4f -> %.4f\n", length(t), V[1], V[end])
@assert !any(isnan, V) "电压出现 NaN"
if haskey(result, "czm D_max")
    D_max = result["czm D_max"]
    δ_max = result["czm δ_max_n [m]"]
    @printf("D_max(end)=%.6f, δ_max_n(end)=%.4e m\n", D_max[end], maximum(δ_max))
    @assert !any(isnan, D_max) "D_max 出现 NaN"
    @assert all(0.0 .<= D_max .<= 1.0) "D_max 超出 [0,1]"
end
println("SMOKE OK")
