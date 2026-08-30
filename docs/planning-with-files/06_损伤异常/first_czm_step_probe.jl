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
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    return case
end

function solve_czm_state(case, czm_mesh, dT_elem, dsn_elem, dsp_elem, u_prev; tol=1e-6)
    E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
    result, new_mesh = JuBat.solve_czm_step(
        deepcopy(czm_mesh),
        zeros(Float64, 2 * czm_mesh.nnode),
        E_eff, ν_eff, case.param_dim.cohesive, case.param, u_prev;
        α_eff=α_eff,
        β_n=β_n,
        β_p=β_p,
        dT_elem=dT_elem,
        Δsoc_n_elem=dsn_elem,
        Δsoc_p_elem=dsp_elem,
        max_iter=30,
        tol=tol,
        n_load_steps=10,
        arc_length_alpha=1.0,
        iter_method="basic"
    )
    return result, new_mesh
end

case = build_case()
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M0, K0, F0, vars0, yphi0 = JuBat.CallModel(case, y0, 0.0; jacobi="update")
vc = 1:size(M0, 1)
dt_init = 1e-8
y_c0 = (M0 - K0 * dt_init) \ (M0 * y0[vc] + F0 * dt_init)
y_old = vcat(y_c0, yphi0)

nT = case.layout.nT
T0 = y_c0[(end - nT + 1):end]
dT0, dsn0, dsp0 = JuBat.compute_czm_strain_inputs(case, vars0, case.czm_mesh, T0)
E_eff0, ν_eff0, α_eff0, β_n0, β_p0 = JuBat.compute_czm_effective_params(case)
Fth0 = JuBat.assemble_thermal_chemical_load(case.czm_mesh, E_eff0, ν_eff0, α_eff0, β_n0, β_p0, dT0, dsn0, dsp0)

println("=== Initial state probe (t = 0) ===")
@printf("T0 range: %.6f .. %.6f\n", minimum(T0), maximum(T0))
@printf("dT0 range: %.6f .. %.6f, norm = %.6e\n", minimum(dT0), maximum(dT0), norm(dT0))
@printf("dsn0 range: %.6f .. %.6f, norm = %.6e\n", minimum(dsn0), maximum(dsn0), norm(dsn0))
@printf("dsp0 range: %.6f .. %.6f, norm = %.6e\n", minimum(dsp0), maximum(dsp0), norm(dsp0))
@printf("Fth0 norm = %.6e, max|Fth0| = %.6e\n", norm(Fth0), maximum(abs.(Fth0)))
res0, mesh0 = solve_czm_state(case, case.czm_mesh, dT0, dsn0, dsp0, zeros(Float64, 2 * case.czm_mesh.nnode))
@printf("t=0 basic converged = %s, iterations = %d, residual = %.6e\n", string(res0.converged), res0.iterations, res0.residual_norm)
@printf("t=0 displacement norm = %.6e, max|u| = %.6e\n", norm(res0.displacement), maximum(abs.(res0.displacement)))

println("=== First step probe (after one coupled solve step) ===")
dt = case.opt.dt[1] / case.param.scale.t0
theta = 0.5
M1, K1, F1, vars1, yphi1 = JuBat.CallModel(case, y_old, dt; jacobi="update")
Mt = M1 - theta * K1 * dt
Kt = (1.0 - theta) * K0 * dt + M1
Ft = theta * F1 * dt + (1.0 - theta) * F0 * dt
y_c1 = Mt \ (Kt * y_old[vc] + Ft)
T1 = y_c1[(end - nT + 1):end]
dT1, dsn1, dsp1 = JuBat.compute_czm_strain_inputs(case, vars1, case.czm_mesh, T1)
Fth1 = JuBat.assemble_thermal_chemical_load(case.czm_mesh, E_eff0, ν_eff0, α_eff0, β_n0, β_p0, dT1, dsn1, dsp1)

@printf("T1 range: %.6f .. %.6f\n", minimum(T1), maximum(T1))
@printf("dT1 range: %.6f .. %.6f, norm = %.6e\n", minimum(dT1), maximum(dT1), norm(dT1))
@printf("dsn1 range: %.6f .. %.6f, norm = %.6e\n", minimum(dsn1), maximum(dsn1), norm(dsn1))
@printf("dsp1 range: %.6f .. %.6f, norm = %.6e\n", minimum(dsp1), maximum(dsp1), norm(dsp1))
@printf("Fth1 norm = %.6e, max|Fth1| = %.6e\n", norm(Fth1), maximum(abs.(Fth1)))

res1_zero, mesh1_zero = solve_czm_state(case, case.czm_mesh, dT1, dsn1, dsp1, zeros(Float64, 2 * case.czm_mesh.nnode))
res1_seed, mesh1_seed = solve_czm_state(case, case.czm_mesh, dT1, dsn1, dsp1, res0.displacement)
res1_loose, _ = solve_czm_state(case, case.czm_mesh, dT1, dsn1, dsp1, zeros(Float64, 2 * case.czm_mesh.nnode); tol=3e-4)
res1_very_loose, _ = solve_czm_state(case, case.czm_mesh, dT1, dsn1, dsp1, zeros(Float64, 2 * case.czm_mesh.nnode); tol=1e-3)

@printf("first-step basic (zero init) converged = %s, iterations = %d, residual = %.6e\n", string(res1_zero.converged), res1_zero.iterations, res1_zero.residual_norm)
@printf("first-step basic (seeded with t=0 solution) converged = %s, iterations = %d, residual = %.6e\n", string(res1_seed.converged), res1_seed.iterations, res1_seed.residual_norm)
@printf("first-step basic (tol=3e-4) converged = %s, iterations = %d, residual = %.6e\n", string(res1_loose.converged), res1_loose.iterations, res1_loose.residual_norm)
@printf("first-step basic (tol=1e-3) converged = %s, iterations = %d, residual = %.6e\n", string(res1_very_loose.converged), res1_very_loose.iterations, res1_very_loose.residual_norm)
@printf("first-step zero-init max|u| = %.6e, seeded max|u| = %.6e\n", maximum(abs.(res1_zero.displacement)), maximum(abs.(res1_seed.displacement)))

cache = JuBat.ensure_czm_cache(case, case.czm_mesh, E_eff0, ν_eff0)
u0 = zeros(Float64, 2 * case.czm_mesh.nnode)
K_total, f_int_total, separations, tractions = JuBat.assemble_coupled_system(
    case.czm_mesh, u0, E_eff0, ν_eff0, case.param_dim.cohesive;
    damage_states=case.czm_mesh.damage_states,
    K_bulk_cached=cache.K_bulk,
    geom_cache=cache.cohesive_geom,
    ws=cache.ws
)
R0 = Fth1 - f_int_total
JuBat.apply_czm_dirichlet!(R0, cache.bc_dofs, zeros(length(cache.bc_dofs)))
K_bc, R_bc = JuBat.apply_bc_czm(K_total, R0; bc_dofs=cache.bc_dofs, bc_vals=zeros(length(cache.bc_dofs)))
Δu = K_bc \ R_bc

u_plus = copy(u0)
u_plus .= u0 .+ Δu
JuBat.apply_czm_dirichlet!(u_plus, cache.bc_dofs, cache.bc_vals)
_, f_int_plus, _, _ = JuBat.assemble_coupled_system(
    case.czm_mesh, u_plus, E_eff0, ν_eff0, case.param_dim.cohesive;
    damage_states=case.czm_mesh.damage_states,
    K_bulk_cached=cache.K_bulk,
    geom_cache=cache.cohesive_geom,
    ws=cache.ws
)
R_plus = Fth1 - f_int_plus
JuBat.apply_czm_dirichlet!(R_plus, cache.bc_dofs, zeros(length(cache.bc_dofs)))

u_minus = copy(u0)
u_minus .= u0 .- Δu
JuBat.apply_czm_dirichlet!(u_minus, cache.bc_dofs, cache.bc_vals)
_, f_int_minus, _, _ = JuBat.assemble_coupled_system(
    case.czm_mesh, u_minus, E_eff0, ν_eff0, case.param_dim.cohesive;
    damage_states=case.czm_mesh.damage_states,
    K_bulk_cached=cache.K_bulk,
    geom_cache=cache.cohesive_geom,
    ws=cache.ws
)
R_minus = Fth1 - f_int_minus
JuBat.apply_czm_dirichlet!(R_minus, cache.bc_dofs, zeros(length(cache.bc_dofs)))

@printf("manual first-step residual at u=0 = %.6e\n", norm(R0))
@printf("manual first-step |Δu| = %.6e\n", norm(Δu))
@printf("manual first-step residual after u+Δu = %.6e\n", norm(R_plus))
@printf("manual first-step residual after u-Δu = %.6e\n", norm(R_minus))

free_dofs = setdiff(collect(1:length(R0)), cache.bc_dofs)
K_free = K_total[free_dofs, free_dofs]
R_free = R0[free_dofs]
Δu_free = K_free \ R_free
u_exact = zeros(Float64, length(R0))
u_exact[free_dofs] .= Δu_free
_, f_int_exact, _, _ = JuBat.assemble_coupled_system(
    case.czm_mesh, u_exact, E_eff0, ν_eff0, case.param_dim.cohesive;
    damage_states=case.czm_mesh.damage_states,
    K_bulk_cached=cache.K_bulk,
    geom_cache=cache.cohesive_geom,
    ws=cache.ws
)
R_exact = Fth1 - f_int_exact
JuBat.apply_czm_dirichlet!(R_exact, cache.bc_dofs, zeros(length(cache.bc_dofs)))
@printf("manual first-step residual after exact BC elimination = %.6e\n", norm(R_exact))

_, _, separations_exact, _ = JuBat.assemble_coupled_system(
    case.czm_mesh, u_exact, E_eff0, ν_eff0, case.param_dim.cohesive;
    damage_states=case.czm_mesh.damage_states,
    K_bulk_cached=cache.K_bulk,
    geom_cache=cache.cohesive_geom,
    ws=cache.ws
)
δn_exact = [s[1] for s in separations_exact]
@printf("exact-elim δn range = %.6e .. %.6e, positive fraction = %.3f\n", minimum(δn_exact), maximum(δn_exact), count(>(0.0), δn_exact) / length(δn_exact))

println("=== Manual Newton trace on first-step state ===")
let u_trace = copy(u0)
    for iter in 1:5
        K_trace, f_int_trace, _, _ = JuBat.assemble_coupled_system(
            case.czm_mesh, u_trace, E_eff0, ν_eff0, case.param_dim.cohesive;
            damage_states=case.czm_mesh.damage_states,
            K_bulk_cached=cache.K_bulk,
            geom_cache=cache.cohesive_geom,
            ws=cache.ws
        )
        R_trace = Fth1 - f_int_trace
        JuBat.apply_czm_dirichlet!(R_trace, cache.bc_dofs, zeros(length(cache.bc_dofs)))
        R_trace_norm = norm(R_trace)
        _, _, separations_trace, _ = JuBat.assemble_coupled_system(
            case.czm_mesh, u_trace, E_eff0, ν_eff0, case.param_dim.cohesive;
            damage_states=case.czm_mesh.damage_states,
            K_bulk_cached=cache.K_bulk,
            geom_cache=cache.cohesive_geom,
            ws=cache.ws
        )
        δn_trace = [s[1] for s in separations_trace]
        @printf("iter %d: residual = %.6e, max|u| = %.6e, δn range = %.6e .. %.6e, positive frac = %.3f\n", iter, R_trace_norm, maximum(abs.(u_trace)), minimum(δn_trace), maximum(δn_trace), count(>(0.0), δn_trace) / length(δn_trace))
        K_trace_bc, R_trace_bc = JuBat.apply_bc_czm(K_trace, R_trace; bc_dofs=cache.bc_dofs, bc_vals=zeros(length(cache.bc_dofs)))
        Δu_trace = K_trace_bc \ R_trace_bc
        @printf("iter %d: |Δu| = %.6e\n", iter, norm(Δu_trace))

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
            R_trial_norm = norm(R_trial)
            @printf("    ls %d: α = %.3f, residual = %.6e\n", ls, α, R_trial_norm)
            if R_trial_norm < R_trace_norm
                u_trace = u_trial
                accepted = true
                break
            end
            α *= 0.5
        end
        if !accepted
            println("    extra tiny-step scan:")
            for α_extra in (0.004, 0.002, 0.001, 0.0005, 0.00025)
                u_trial = u_trace .+ α_extra .* Δu_trace
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
                @printf("      α = %.5f, residual = %.6e\n", α_extra, norm(R_trial))
            end
            println("    line search failed; stopping trace")
            break
        end
    end
end
