# tools/probe_c4lite_free_core.jl
#
# D10 复核探针：自由芯部（Γ_in,free）下扫描用户给定物理包络
# |Δsoc|≤1、|ΔT|≤20°C，并覆盖 2 个卷绕预应力倍率。仅调用力学/CZM
# 求解器，不接入电化学或热反馈；超范围载荷不作为能力验收点。
# 输出 → output/probe_c4lite_free_core/probe_results.csv（AGENTS §9.9）。

using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

const OUT = joinpath(@__DIR__, "..", "output", "probe_c4lite_free_core")

function build_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true
    opt.czm_geo_nonlinear = true
    opt.czm_j2_plasticity = true
    opt.czm_winding_prestress = true
    opt.czm_fix_inner = false

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(
        case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(
        mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    param_cache = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(
        case, case.czm_mesh, param_cache; fix_inner=case.opt.czm_fix_inner)
    layout = JuBat.CzmLayout(case.czm_mesh)
    layout.node_ref = copy(case.czm_mesh.node)
    layout.plastic_states = [JuBat.PlasticState()
        for _ in 1:size(case.czm_mesh.bulk_element, 1), _ in 1:4]
    return case, param_cache, cache, layout
end

function run_probe(
    prestress_scale::Float64,
    delta_soc::Float64;
    delta_T_C::Float64=0.0,
    solver::String="load_substep",
    plasticity::Bool=true,
)
    case, param_cache, cache, layout = build_fixture()
    mesh = case.czm_mesh
    ndof = 2 * mesh.nnode
    ne = size(mesh.bulk_element, 1)
    alpha = param_cache.by_interface[:PE_PCC].α
    beta_n = case.param.NE.Omega / 3
    beta_p = case.param.PE.Omega / 3
    eigenstrain = (
        α_eff=alpha,
        β_n=beta_n,
        β_p=beta_p,
        # α 已按 T_ref 归一，因此温差入口必须使用 ΔT*=ΔT[K]/T_ref。
        dT=fill(delta_T_C / case.param_dim.scale.T_ref, ne),
        Δsn=fill(delta_soc, ne),
        Δsp=fill(delta_soc, ne),
    )
    prestress_full = JuBat.winding_prestress_field(mesh, case.param)
    prestress = [(prestress_scale * a, prestress_scale * b, prestress_scale * c)
                 for (a, b, c) in prestress_full]

    try
        result, updated = JuBat.solve_czm_step(
            mesh, zeros(ndof), param_cache, case.param, layout.u_prev;
            α_eff=alpha,
            β_n=beta_n,
            β_p=beta_p,
            dT_elem=eigenstrain.dT,
            Δsoc_n_elem=eigenstrain.Δsn,
            Δsoc_p_elem=eigenstrain.Δsp,
            max_iter=200,
            tol=1e-8,
            n_load_steps=50,
            iter_method=solver,
            cache=cache,
            geo_nl=true,
            eigenstrain=eigenstrain,
            plasticity=plasticity,
            mech_state=plasticity ? layout.plastic_states : nothing,
            prestress=prestress,
        )
        result.converged || error(
            "mechanical solve returned converged=false (residual=$(result.residual_norm))")

        if plasticity
            JuBat.assemble_coupled_system(
                mesh, result.displacement, param_cache;
                geo_nl=true,
                eigenstrain=eigenstrain,
                plasticity=true,
                mech_state=layout.plastic_states,
                commit_plastic=true,
                prestress=prestress,
            )
        end
        _, delta_core = JuBat.core_ovalization(
            mesh, result.displacement, layout.node_ref)
        damage_max = maximum(state.D for state in updated.damage_states)
        kappa_max = plasticity ?
            maximum(state.kappa for state in layout.plastic_states) : missing
        normal_min = minimum(result.separation_n)
        normal_max = maximum(result.separation_n)
        initiation_ratio = maximum(
            max(result.separation_n[i], 0.0) /
            param_cache.by_interface[mesh.cohesive_elements[i].interface_type].δ_0_n
            for i in eachindex(result.separation_n))
        return (
            prestress_scale=prestress_scale,
            delta_soc=delta_soc,
            delta_T_C=delta_T_C,
            solver=solver,
            plasticity=plasticity,
            status="OK",
            delta_core=delta_core,
            damage_max=damage_max,
            kappa_max=kappa_max,
            normal_min=normal_min,
            normal_max=normal_max,
            initiation_ratio=initiation_ratio,
            iterations=result.iterations,
            residual=result.residual_norm,
            error="",
        )
    catch err
        return (
            prestress_scale=prestress_scale,
            delta_soc=delta_soc,
            delta_T_C=delta_T_C,
            solver=solver,
            plasticity=plasticity,
            status="ERROR",
            delta_core=missing,
            damage_max=missing,
            kappa_max=missing,
            normal_min=missing,
            normal_max=missing,
            initiation_ratio=missing,
            iterations=missing,
            residual=missing,
            error=sprint(showerror, err),
        )
    end
end

csv_field(value) = value === missing ? "NA" : string(value)
csv_field(value::AbstractString) = "\"" * replace(value, '"' => "\"\"") * "\""

function main()
    mkpath(OUT)
    rows = Any[]
    envelope = (
        (-1.0, -20.0), (-1.0, 0.0), (-1.0, 20.0),
        ( 0.0, -20.0),               ( 0.0, 20.0),
        ( 1.0, -20.0), ( 1.0, 0.0), ( 1.0, 20.0),
    )
    for prestress_scale in (0.2, 1.0), (delta_soc, delta_T_C) in envelope
        row = run_probe(prestress_scale, delta_soc; delta_T_C=delta_T_C)
        push!(rows, row)
        if row.status == "OK"
            @printf(
                "solver=%-12s plastic=%5s prestress=%3.1f Δsoc=%4.1f ΔT=%5.1f°C  Δ_core=%10.4e  D_max=%10.4e  κ_max=%10s  δn/δ0=%8.3f  res=%8.2e\n",
                row.solver, string(row.plasticity), row.prestress_scale,
                row.delta_soc, row.delta_T_C, row.delta_core,
                row.damage_max, string(row.kappa_max), row.initiation_ratio, row.residual)
        else
            @printf(
                "solver=%-12s plastic=%5s prestress=%3.1f Δsoc=%4.1f ΔT=%5.1f°C  ERROR: %s\n",
                row.solver, string(row.plasticity), row.prestress_scale,
                row.delta_soc, row.delta_T_C, row.error)
        end
    end
    if "--arc-followup" in ARGS
        # 仅对物理包络内 load-substep 失败点做独立路径跟踪；不是静默回退。
        failed_rows = filter(row -> row.status == "ERROR", rows)
        for failed in failed_rows
            row = run_probe(
                failed.prestress_scale,
                failed.delta_soc;
                delta_T_C=failed.delta_T_C,
                solver="arc_length",
            )
            push!(rows, row)
            if row.status == "OK"
                @printf(
                    "solver=%-12s plastic=%5s prestress=%3.1f Δsoc=%4.1f ΔT=%5.1f°C  Δ_core=%10.4e  D_max=%10.4e  κ_max=%10s  δn/δ0=%8.3f  res=%8.2e\n",
                    row.solver, string(row.plasticity), row.prestress_scale,
                    row.delta_soc, row.delta_T_C, row.delta_core, row.damage_max,
                    string(row.kappa_max), row.initiation_ratio, row.residual)
            else
                @printf(
                    "solver=%-12s plastic=%5s prestress=%3.1f Δsoc=%4.1f ΔT=%5.1f°C  ERROR: %s\n",
                    row.solver, string(row.plasticity), row.prestress_scale,
                    row.delta_soc, row.delta_T_C, row.error)
            end
        end
    end

    output_path = joinpath(OUT, "probe_results.csv")
    open(output_path, "w") do io
        println(io,
            "solver,plasticity,prestress_scale,delta_soc,delta_T_C,status,delta_core,damage_max,kappa_max,normal_min,normal_max,initiation_ratio,iterations,residual,error")
        for row in rows
            fields = (
                row.solver,
                row.plasticity,
                row.prestress_scale,
                row.delta_soc,
                row.delta_T_C,
                row.status,
                row.delta_core,
                row.damage_max,
                row.kappa_max,
                row.normal_min,
                row.normal_max,
                row.initiation_ratio,
                row.iterations,
                row.residual,
                row.error,
            )
            println(io, join(csv_field.(fields), ','))
        end
    end
    println("CSV → ", output_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
