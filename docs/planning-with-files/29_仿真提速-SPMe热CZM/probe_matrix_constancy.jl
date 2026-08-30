# 探针：验证 CallModel 返回的全局 M/K 是否跨状态常量；测量 CZM 单次更新迭代数
# case 构造块逐行复制自 example/testexample.jl:29-86，参数完全一致
include(joinpath(@__DIR__, "../../../src/JuBat.jl"))
using .JuBat

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

t = 0.0
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M1, K1, F1, vars1, _ = JuBat.CallModel(case, y0, t; jacobi="update")

# 状态扰动：整体 ±5%（涵盖浓度与温度 DOF）
y0b = copy(y0)
y0b .*= 1.05
M2, K2, F2, vars2, _ = JuBat.CallModel(case, y0b, t; jacobi="update")

println("=== 矩阵常量性 ===")
println("M identical: ", M1 == M2)
println("K identical: ", K1 == K2)
println("F identical (预期 false): ", F1 == F2)

# 若不等，定位变化块（化学块 chem_range vs 热块 thermal_range）
if M1 != M2
    lay = case.layout
    D = abs.(M1 - M2)
    nz = findall(>(0), D)
    println("M 差异非零数 = ", length(nz))
    chem_cnt = count(idx -> idx[1] <= lay.n_chem * lay.ne && idx[2] <= lay.n_chem * lay.ne, nz)
    println("M 化学块差异数 = ", chem_cnt, ", 其余（热块/耦合）= ", length(nz) - chem_cnt)
end
if K1 != K2
    lay = case.layout
    D = abs.(K1 - K2)
    nz = findall(>(0), D)
    println("K 差异非零数 = ", length(nz))
    chem_cnt = count(idx -> idx[1] <= lay.n_chem * lay.ne && idx[2] <= lay.n_chem * lay.ne, nz)
    println("K 化学块差异数 = ", chem_cnt, ", 其余（热块/耦合）= ", length(nz) - chem_cnt)
end

println("=== CZM 单次更新 ===")
case.czm_layout = JuBat.CzmLayout(case.czm_mesh)   # Solve.jl:200-202 主循环中才建，探针补齐
T_nodes = JuBat.get_thermal_dofs(y0, case.layout)
res = JuBat.update_czm_damage!(case, vars1, T_nodes)
println("iterations = ", res.iterations, ", converged = ", res.converged)
