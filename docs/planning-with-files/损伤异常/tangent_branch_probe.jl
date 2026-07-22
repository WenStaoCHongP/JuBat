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

function summarise_branch(label::AbstractString, separations, damage_states, cohesive_params)
    δn = [s[1] for s in separations]
    δeff = max.(0.0, δn)
    D = [s.D for s in damage_states]
    n_elastic = count(x -> x <= cohesive_params.δ_0_n, δeff)
    n_soft = count(x -> x > cohesive_params.δ_0_n && x < cohesive_params.δ_c_n, δeff)
    n_failed = count(x -> x >= cohesive_params.δ_c_n, δeff)
    n_loading = count(x -> x > -1e-15, δeff)
    @printf("%s: max|δn| = %.6e, max δeff = %.6e, δ0 = %.6e, δc = %.6e, elastic = %d, softening = %d, failed = %d, loading = %d, Dmax = %.6e\n",
        label, maximum(abs.(δn)), maximum(δeff), cohesive_params.δ_0_n, cohesive_params.δ_c_n, n_elastic, n_soft, n_failed, n_loading, isempty(D) ? 0.0 : maximum(D))
end

function matrix_spectrum_summary(K::AbstractMatrix)
    sv = svdvals(Matrix(K))
    σ_max = maximum(sv)
    σ_min = minimum(sv)
    κ = σ_min == 0.0 ? Inf : σ_max / σ_min
    return σ_max, σ_min, κ
end
case = build_case()
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M0, K0, F0, vars0, yphi0 = JuBat.CallModel(case, y0, 0.0; jacobi="update")
vc = 1:size(M0, 1)
dt_init = 1e-8
y_c0 = (M0 - K0 * dt_init) \ (M0 * y0[vc] + F0 * dt_init)
y_old = vcat(y_c0, yphi0)

nT = case.layout.nT
dt = case.opt.dt[1] / case.param.scale.t0
theta = 0.5
E_eff0, ν_eff0, α_eff0, β_n0, β_p0 = JuBat.compute_czm_effective_params(case)
M1, K1, F1, vars1, yphi1 = JuBat.CallModel(case, y_old, dt; jacobi="update")
Mt = M1 - theta * K1 * dt
Kt = (1.0 - theta) * K0 * dt + M1
Ft = theta * F1 * dt + (1.0 - theta) * F0 * dt
y_c1 = Mt \ (Kt * y_old[vc] + Ft)
T1 = y_c1[(end - nT + 1):end]
dT1, dsn1, dsp1 = JuBat.compute_czm_strain_inputs(case, vars1, case.czm_mesh, T1)
Fth1 = JuBat.assemble_thermal_chemical_load(case.czm_mesh, E_eff0, ν_eff0, α_eff0, β_n0, β_p0, dT1, dsn1, dsp1)
cache = JuBat.ensure_czm_cache(case, case.czm_mesh, E_eff0, ν_eff0)

println("=== Thresholds ===")
@printf("cohesive δ0_n = %.6e, δc_n = %.6e, K_n = %.6e, K_t = %.6e\n", case.param_dim.cohesive.δ_0_n, case.param_dim.cohesive.δ_c_n, case.param_dim.cohesive.K_n, case.param_dim.cohesive.K_t)
@printf("normalized δ0_n = %.6e, δc_n = %.6e\n", case.param.cohesive.δ_0_n, case.param.cohesive.δ_c_n)
@printf("Fth1 norm = %.6e\n", norm(Fth1))

let u_trace = zeros(Float64, 2 * case.czm_mesh.nnode)
println("=== Newton trace with branch summary ===")
for iter in 1:4
    K_trace, f_int_trace, separations_trace, _ = JuBat.assemble_coupled_system(
        case.czm_mesh, u_trace, E_eff0, ν_eff0, case.param_dim.cohesive;
        damage_states=case.czm_mesh.damage_states,
        K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom,
        ws=cache.ws
    )
    R_trace = Fth1 - f_int_trace
    JuBat.apply_czm_dirichlet!(R_trace, cache.bc_dofs, zeros(length(cache.bc_dofs)))
    R_trace_norm = norm(R_trace)
    @printf("iter %d: residual = %.6e, max|u| = %.6e\n", iter, R_trace_norm, maximum(abs.(u_trace)))
    summarise_branch("  branch", separations_trace, case.czm_mesh.damage_states, case.param_dim.cohesive)

    K_trace_bc, R_trace_bc = JuBat.apply_bc_czm(K_trace, R_trace; bc_dofs=cache.bc_dofs, bc_vals=zeros(length(cache.bc_dofs)))
    σ_max_bc, σ_min_bc, κ_bc = matrix_spectrum_summary(K_trace_bc)
    free_dofs = setdiff(collect(1:length(R_trace)), cache.bc_dofs)
    K_free = K_trace[free_dofs, free_dofs]
    σ_max_free, σ_min_free, κ_free = matrix_spectrum_summary(K_free)
    @printf("  K_bc spectrum: σ_max = %.6e, σ_min = %.6e, cond2 = %.6e\n", σ_max_bc, σ_min_bc, κ_bc)
    @printf("  K_free spectrum: σ_max = %.6e, σ_min = %.6e, cond2 = %.6e\n", σ_max_free, σ_min_free, κ_free)
    Δu_trace = K_trace_bc \ R_trace_bc
    @printf("  |Δu| = %.6e\n", norm(Δu_trace))

    accepted = false
    α = 1.0
    for ls in 1:8
        u_trial = u_trace .+ α .* Δu_trace
        JuBat.apply_czm_dirichlet!(u_trial, cache.bc_dofs, cache.bc_vals)
        _, f_int_trial, _, _ = JuBat.assemble_coupled_system(
            case.czm_mesh, u_trial, E_eff0, ν_eff0, case.param_dim.cohesive;
            damage_states=case.czm_mesh.damage_states,
            K_bulk_cached=cache.K_bulk,
            geom_cache=cache.cohesive_geom,
            ws=cache.ws
        )
        R_trial = Fth1 - f_int_trial
        JuBat.apply_czm_dirichlet!(R_trial, cache.bc_dofs, zeros(length(cache.bc_dofs)))
        @printf("    ls %d: α = %.6f, residual = %.6e\n", ls, α, norm(R_trial))
        if norm(R_trial) < R_trace_norm
            u_trace = u_trial
            accepted = true
            break
        end
        α *= 0.5
    end

    if !accepted
        println("  line search failed; stopping trace")
        break
    end
end
end
