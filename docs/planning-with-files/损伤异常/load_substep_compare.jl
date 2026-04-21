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
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=16, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
    return case
end

function build_first_step_state(case)
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    M0, K0, F0, vars0, yphi0 = JuBat.CallModel(case, y0, 0.0; jacobi="update")
    vc = 1:size(M0, 1)
    dt_init = 1e-8
    y_c0 = (M0 - K0 * dt_init) \ (M0 * y0[vc] + F0 * dt_init)
    y_old = vcat(y_c0, yphi0)

    dt = case.opt.dt[1] / case.param.scale.t0
    theta = 0.5
    M1, K1, F1, vars1, _ = JuBat.CallModel(case, y_old, dt; jacobi="update")
    Mt = M1 - theta * K1 * dt
    Kt = (1.0 - theta) * K0 * dt + M1
    Ft = theta * F1 * dt + (1.0 - theta) * F0 * dt
    y_c1 = Mt \ (Kt * y_old[vc] + Ft)
    T1 = y_c1[(end - case.layout.nT + 1):end]
    dT1, dsn1, dsp1 = JuBat.compute_czm_strain_inputs(case, vars1, case.czm_mesh, T1)
    return dT1, dsn1, dsp1
end

function run_compare()
    case = build_case()
    dT1, dsn1, dsp1 = build_first_step_state(case)
    E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
    u0 = zeros(Float64, 2 * case.czm_mesh.nnode)

    println("=== load_substep comparison on the first coupled CZM state ===")
    @printf("dT norm = %.6e, dsn norm = %.6e, dsp norm = %.6e\n", norm(dT1), norm(dsn1), norm(dsp1))
    @printf("default czm_tol = %.6e\n", case.opt.czm_tol)

    for (label, method, steps, tol) in [
        ("basic", "basic", 0, 1e-4),
        ("load_substep-2", "load_substep", 2, 1e-4),
        ("load_substep-5", "load_substep", 5, 1e-4),
        ("load_substep-10", "load_substep", 10, 1e-4),
        ("load_substep-20", "load_substep", 20, 1e-4),
    ]
        result, _ = JuBat.solve_czm_step(
            deepcopy(case.czm_mesh),
            zeros(Float64, 2 * case.czm_mesh.nnode),
            E_eff, ν_eff, case.param.cohesive, case.param, u0;
            α_eff=α_eff,
            β_n=β_n,
            β_p=β_p,
            dT_elem=dT1,
            Δsoc_n_elem=dsn1,
            Δsoc_p_elem=dsp1,
            max_iter=100,
            tol=tol,
            n_load_steps=steps == 0 ? 2 : steps,
            arc_length_alpha=1.0,
            iter_method=method,
        )
        @printf("%s: converged = %s, iterations = %d, residual = %.6e, max|u| = %.6e\n",
            label, string(result.converged), result.iterations, result.residual_norm, maximum(abs.(result.displacement)))
    end
end

run_compare()
