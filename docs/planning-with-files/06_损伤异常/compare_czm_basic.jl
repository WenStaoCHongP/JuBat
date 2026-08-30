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

case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=16, nθ_czm=40, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M, K, F, variables, _ = JuBat.CallModel_MultiSPMe(case, y0, 0.0; jacobi="update")

E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
T_nodes = variables["thermal2D temperature at nodes"]
dT_elem, Δsoc_n_elem, Δsoc_p_elem = JuBat.compute_czm_strain_inputs(case, variables, case.czm_mesh, T_nodes)

u0 = zeros(Float64, 2 * case.czm_mesh.nnode)
result_basic, mesh_basic = JuBat.solve_czm_step(
    deepcopy(case.czm_mesh), zeros(Float64, 2 * case.czm_mesh.nnode),
    E_eff, ν_eff, case.param_dim.cohesive, case.param, u0;
    α_eff=α_eff, β_n=β_n, β_p=β_p,
    dT_elem=dT_elem,
    Δsoc_n_elem=Δsoc_n_elem,
    Δsoc_p_elem=Δsoc_p_elem,
    max_iter=30, tol=1e-6, n_load_steps=10, arc_length_alpha=1.0, iter_method="basic")

println("basic converged = ", result_basic.converged)
println("basic iterations = ", result_basic.iterations)
println("basic residual = ", result_basic.residual_norm)
println("basic damage max = ", maximum(result_basic.damage))
println("mesh damage max = ", maximum(s.D for s in mesh_basic.damage_states))
