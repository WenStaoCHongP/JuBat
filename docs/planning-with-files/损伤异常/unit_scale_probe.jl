include(joinpath(@__DIR__, "..", "..", "..", "src", "JuBat.jl"))
using .JuBat
using LinearAlgebra
using Printf

function build_case()
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
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_data.thermal2D, case.param)
    return case
end

function branch_summary(separations, damage_states, cohesive_params)
    δn = [s[1] for s in separations]
    δeff = max.(0.0, δn)
    n_elastic = count(x -> x <= cohesive_params.δ_0_n, δeff)
    n_soft = count(x -> x > cohesive_params.δ_0_n && x < cohesive_params.δ_c_n, δeff)
    n_failed = count(x -> x >= cohesive_params.δ_c_n, δeff)
    Dmax = isempty(damage_states) ? 0.0 : maximum(s.D for s in damage_states)
    return maximum(abs.(δn)), maximum(δeff), n_elastic, n_soft, n_failed, Dmax
end

function spectrum_summary(K::AbstractMatrix)
    sv = svdvals(Matrix(K))
    return maximum(sv), minimum(sv), maximum(sv) / minimum(sv)
end

case = build_case()
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M0, K0, F0, vars0, yphi0 = JuBat.CallModel(case, y0, 0.0; jacobi="update")
vc = 1:size(M0, 1)
dt_init = 1e-8
y_c0 = (M0 - K0 * dt_init) \ (M0 * y0[vc] + F0 * dt_init)
y_old = vcat(y_c0, yphi0)

dt = case.opt.dt[1] / case.param.scale.t0
theta = 0.5
E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
M1, K1, F1, vars1, yphi1 = JuBat.CallModel(case, y_old, dt; jacobi="update")
Mt = M1 - theta * K1 * dt
Kt = (1.0 - theta) * K0 * dt + M1
Ft = theta * F1 * dt + (1.0 - theta) * F0 * dt
y_c1 = Mt \ (Kt * y_old[vc] + Ft)
T1 = y_c1[(end - case.layout.nT + 1):end]
dT1, dsn1, dsp1 = JuBat.compute_czm_strain_inputs(case, vars1, case.czm_mesh, T1)
Fth1 = JuBat.assemble_thermal_chemical_load(case.czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT1, dsn1, dsp1)
cache = JuBat.ensure_czm_cache(case, case.czm_mesh, E_eff, ν_eff)

coh_phys = deepcopy(case.param_dim.cohesive)
coh_norm = deepcopy(case.param.cohesive)

println("=== Cohesive scale comparison ===")
@printf("physical K_n = %.6e, K_t = %.6e, δ0_n = %.6e, δc_n = %.6e\n", coh_phys.K_n, coh_phys.K_t, coh_phys.δ_0_n, coh_phys.δ_c_n)
@printf("normalized K_n = %.6e, K_t = %.6e, δ0_n = %.6e, δc_n = %.6e\n", coh_norm.K_n, coh_norm.K_t, coh_norm.δ_0_n, coh_norm.δ_c_n)
@printf("scale K_czm = %.6e, σ_czm = %.6e, δ_czm = %.6e\n", case.param_dim.scale.K_czm, case.param_dim.scale.σ_czm, case.param_dim.scale.δ_czm)
@printf("Fth1 norm = %.6e\n", norm(Fth1))

for (label, coh) in [("physical", coh_phys), ("normalized", coh_norm)]
    K_total, f_int_total, separations, _ = JuBat.assemble_coupled_system(
        case.czm_mesh, zeros(Float64, 2 * case.czm_mesh.nnode), E_eff, ν_eff, coh;
        damage_states=case.czm_mesh.damage_states,
        K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom,
        ws=cache.ws,
    )
    R = Fth1 - f_int_total
    JuBat.apply_czm_dirichlet!(R, cache.bc_dofs, zeros(length(cache.bc_dofs)))
    K_bc, R_bc = JuBat.apply_bc_czm(K_total, R; bc_dofs=cache.bc_dofs, bc_vals=zeros(length(cache.bc_dofs)))
    σ_max, σ_min, κ = spectrum_summary(K_bc)
    δn_max, δeff_max, n_elastic, n_soft, n_failed, Dmax = branch_summary(separations, case.czm_mesh.damage_states, coh)
    @printf("%s K_bc: σ_max = %.6e, σ_min = %.6e, cond2 = %.6e\n", label, σ_max, σ_min, κ)
    @printf("%s branch: max|δn| = %.6e, max δeff = %.6e, elastic = %d, softening = %d, failed = %d, Dmax = %.6e\n", label, δn_max, δeff_max, n_elastic, n_soft, n_failed, Dmax)
    result, _ = JuBat.solve_czm_step(
        deepcopy(case.czm_mesh), zeros(Float64, 2 * case.czm_mesh.nnode),
        E_eff, ν_eff, coh, case.param, zeros(Float64, 2 * case.czm_mesh.nnode);
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT1, Δsoc_n_elem=dsn1, Δsoc_p_elem=dsp1,
        max_iter=30, tol=1e-6, n_load_steps=10, arc_length_alpha=1.0,
        iter_method="basic", cache=cache,
    )
    @printf("%s solve: converged = %s, iterations = %d, residual = %.6e\n", label, string(result.converged), result.iterations, result.residual_norm)
end
