# ========================================================================
# CZM Post-processing, Statistics & Damage Management
# ========================================================================
# 从 CzmSolve.jl 拆分出来的后处理、统计和损伤管理函数。
# 这些函数无状态、纯计算，独立于求解器主循环。
# 损伤状态自 2026-08-30 重构起存于 MechState.damage_states——
# 本文件函数直接接收 damage_states 向量（不再从 czm_mesh 读取）。
# ========================================================================

"""
    get_damage_statistics(damage_states)
"""
function get_damage_statistics(damage_states::AbstractVector{DamageState})
    n = length(damage_states)

    if n == 0
        return (max_D=0.0, mean_D=0.0, min_D=0.0, n_fractured=0, fraction_damaged=0.0, total_accumulated=0.0)
    end

    D_vals = [s.D for s in damage_states]
    n_fractured = count(s -> s.fractured, damage_states)
    n_damaged = count(d -> d > 0.01, D_vals)
    accumulated = [s.accumulated_damage for s in damage_states]

    return (
        max_D = maximum(D_vals),
        mean_D = mean(D_vals),
        min_D = minimum(D_vals),
        n_fractured = n_fractured,
        fraction_damaged = n_damaged / n,
        total_accumulated = sum(accumulated)
    )
end

"""
    check_fracture_criterion(damage_states; threshold=0.99)
"""
function check_fracture_criterion(damage_states::AbstractVector{DamageState}; threshold::Float64=0.99)
    stats = get_damage_statistics(damage_states)

    is_fractured_avg = stats.mean_D >= threshold
    is_fractured_count = (stats.n_fractured / max(1, length(damage_states))) > 0.5

    is_fractured = is_fractured_avg || is_fractured_count

    fracture_info = (
        is_fractured = is_fractured,
        criterion = is_fractured_avg ? :average_damage : (is_fractured_count ? :fractured_count : :none),
        stats = stats
    )

    return is_fractured, fracture_info
end

"""
    reset_damage_states!(ms)

原位重置 MechState 上的损伤状态。
"""
function reset_damage_states!(ms::MechState)
    if ms.gp_damage_states === nothing
        ms.mechanics_initialized && throw(ArgumentError("initialized mechanics lacks required GP history"))
        new_gp = nothing
    else
        size(ms.gp_damage_states,1) == length(ms.damage_states) ||
            throw(DimensionMismatch("GP and element damage history lengths differ"))
        new_gp = [DamageState() for _ in axes(ms.gp_damage_states,1), _ in axes(ms.gp_damage_states,2)]
    end
    new_elements = [DamageState() for _ in eachindex(ms.damage_states)]
    ms.damage_states = new_elements
    ms.gp_damage_states = new_gp
    return ms
end

function check_cycle_gp_history(gp_states::Matrix{DamageState}, n_elements::Int)
    size(gp_states,1) == n_elements || throw(DimensionMismatch(
        "GP and element damage history lengths differ"))
    # The cohesive assembler uses NCweight(order), whose order is the number
    # of points (2–5). The actual GP shape supplies this existing rule.
    weights = NCweight(size(gp_states,2))[1]
    for s in gp_states
        all(isfinite,(s.D,s.D_visc,s.δ_max_n,s.δ_max_t,s.δ_max_eff,s.accumulated_damage)) &&
            0.0 <= s.D_visc <= s.D <= 1.0 &&
            min(s.δ_max_n,s.δ_max_t,s.δ_max_eff,s.accumulated_damage) >= 0.0 ||
                throw(ArgumentError("invalid cycle GP damage history"))
        if s.fractured
            s.D == s.D_visc == 1.0 || throw(ArgumentError(
                "fractured GP history requires D = D_visc = 1"))
        elseif s.D == 1.0 && s.D_visc >= 1.0 - 1e-10
            throw(ArgumentError("noncanonical GP saturation history"))
        end
    end
    return weights
end

"""
    accumulate_cycle_damage!(ms, cycle_damage_increment)

对真实GP的 accumulated_damage 加给定非负增量；达到1则D=1。
此接口不含物理时间，故保留D_visc；仅D和D_visc均规范饱和时标记断裂。
零增量完整no-op；非零增量要求已有GP历史，不从元素均值重建。
"""
function accumulate_cycle_damage!(ms::MechState, cycle_damage_increment::Float64)
    isfinite(cycle_damage_increment) && cycle_damage_increment >= 0.0 ||
        throw(ArgumentError("cycle damage increment must be finite and nonnegative"))
    cycle_damage_increment == 0.0 && return ms
    weights = check_cycle_gp_history(ms.gp_damage_states,length(ms.damage_states))
    new_gp = clone_damage_states(ms.gp_damage_states)
    for s in new_gp
        if !s.fractured
            s.accumulated_damage += cycle_damage_increment
            if s.accumulated_damage >= 1.0
                s.D = 1.0
            end
            s.fractured = s.D == 1.0 && s.D_visc >= 1.0 - 1e-10
            s.fractured && (s.D_visc = 1.0)
        end
    end
    check_cycle_gp_history(new_gp,length(ms.damage_states))
    new_damage_states = aggregate_gp_damage_states(new_gp,weights)
    all(s->all(isfinite,(s.D,s.D_visc,s.δ_max_n,s.δ_max_t,s.δ_max_eff,s.accumulated_damage)),
        new_damage_states) || throw(ArgumentError("non-finite aggregated cycle damage history"))
    ms.damage_states = new_damage_states
    ms.gp_damage_states = new_gp
    return ms
end

"""
    czm_output_to_variables(result, variables; committed_states)

Export a successful step with its real committed element history. Equilibrium
damage alone does not determine whether viscous traction has fully vanished.
"""
function czm_output_to_variables(result::CZMResult, variables::Dict{String, Union{Array{Float64}, Float64}};
                                 committed_states::AbstractVector{DamageState})
    check_czm_committed_output(result, committed_states)
    new_variables = copy(variables)
    new_variables["czm displacement x"] = result.displacement[1:2:end]
    new_variables["czm displacement y"] = result.displacement[2:2:end]
    new_variables["czm damage"] = copy(result.damage)
    new_variables["czm traction normal"] = copy(result.traction_n)
    new_variables["czm traction tangent"] = copy(result.traction_t)
    new_variables["czm separation normal"] = copy(result.separation_n)
    new_variables["czm separation tangent"] = copy(result.separation_t)
    separation_eff = hypot.(max.(result.separation_n, 0.0), result.separation_t)
    new_variables["czm separation effective"] = separation_eff

    new_variables["czm D_max"] = maximum(result.damage)
    new_variables["czm D_mean"] = mean(result.damage)
    new_variables["czm δ_max_eff"] = [maximum(s.δ_max_eff for s in committed_states)]
    new_variables["czm n_fractured"] = Float64(count(result.fractured))

    return new_variables
end

"""Validate a successful physical output; failed CZMResult remains available for diagnostics."""
function check_czm_physical_output(result::Union{CZMResult,CZMSnapshot})
    result.converged || throw(ArgumentError("cannot publish an unconverged CZM physical result"))
    n = length(result.damage)
    n > 0 && !isempty(result.displacement) && iseven(length(result.displacement)) ||
        throw(DimensionMismatch("CZM physical output requires elements and interleaved displacement"))
    all(a -> length(a) == n, (result.damage_visc, result.fractured,
        result.separation_n, result.separation_t, result.traction_n, result.traction_t)) ||
        throw(DimensionMismatch("CZM physical output element fields have different lengths"))
    all(a -> all(isfinite, a), (result.displacement, result.damage, result.damage_visc,
        result.separation_n, result.separation_t, result.traction_n, result.traction_t)) ||
        throw(ArgumentError("CZM physical output contains non-finite fields"))
    isfinite(result.residual_norm) && result.residual_norm >= 0.0 && result.iterations >= 0 ||
        throw(ArgumentError("CZM physical output has invalid solver diagnostics"))
    for i in eachindex(result.damage)
        0.0 <= result.damage_visc[i] <= result.damage[i] <= 1.0 ||
            throw(ArgumentError("CZM physical output requires 0 <= D_visc <= D <= 1"))
        result.fractured[i] && !(result.damage[i] == result.damage_visc[i] == 1.0) &&
            throw(ArgumentError("fractured CZM output requires D = D_visc = 1"))
    end
    return result
end

function check_czm_committed_output(result::CZMResult, committed_states::AbstractVector{DamageState})
    check_czm_physical_output(result)
    length(committed_states) == length(result.damage) || throw(DimensionMismatch(
        "committed CZM history must match result damage length"))
    for (i, state) in enumerate(committed_states)
        all(isfinite, (state.D, state.D_visc, state.δ_max_n, state.δ_max_t,
                       state.δ_max_eff, state.accumulated_damage)) || throw(ArgumentError(
            "committed CZM history contains non-finite values at element $i"))
        state.D == result.damage[i] && state.D_visc == result.damage_visc[i] &&
            state.fractured == result.fractured[i] || throw(ArgumentError(
            "committed CZM D/D_visc/fractured differs from result at element $i"))
    end
    return result
end

"""
    CZMSnapshot(result, mech, mesh, param; time_s, mechanical_dt_s, cycle, phase, method, czm_model, plasticity, czm_workspace)

Copy the current accepted mechanical output. GP separations/tractions are captured
by the production COH loop at fixed accepted history and beta=0; no history advances.
Bulk J2 kappa keeps its own element/GP topology and is absent when plasticity is off.
"""
function CZMSnapshot(result::CZMResult, mech::MechState, mesh::CohesiveMesh, param::Params;
        time_s::Float64, mechanical_dt_s::Float64, cycle::Int, phase::String,
        method::String, czm_model::String, plasticity::Bool, czm_workspace::CZMWorkspace)
    check_czm_committed_output(result, mech.damage_states)
    mech.mechanics_initialized || throw(ArgumentError("snapshot requires an accepted mechanical step"))
    length(result.displacement) == 2size(mesh.node,1) ||
        throw(DimensionMismatch("snapshot displacement differs from the actual mechanical mesh"))
    result.displacement == mech.u_prev || throw(ArgumentError("snapshot displacement differs from accepted mechanics"))
    gp = clone_damage_states(mech.gp_damage_states)
    capture = (separation_n=zeros(size(gp)), separation_t=zeros(size(gp)),
        traction_n=zeros(size(gp)), traction_t=zeros(size(gp)))
    check_czm_workspace(mesh, czm_workspace)
    assemble_czm_system(mesh, result.displacement, param;
        czm_workspace=czm_workspace, gp_damage_states=gp, gp_trial_states=nothing,
        gp_capture=capture, visc_beta=0.0, czm_model=czm_model)
    bulk_kappa = if plasticity
        size(mech.plastic_states) == (size(mesh.bulk_element,1),4) ||
            throw(DimensionMismatch("J2 snapshot history must use actual bulk element/GP topology"))
        map(s -> s.kappa, mech.plastic_states)
    else
        nothing
    end
    snapshot = CZMSnapshot(time_s,cycle,phase,copy(result.displacement),copy(result.damage),
        copy(result.damage_visc),copy(result.fractured),copy(result.separation_n),copy(result.separation_t),
        copy(result.traction_n),copy(result.traction_t),result.converged,result.iterations,
        result.residual_norm,method,mechanical_dt_s,czm_model,gp,
        copy(NCweight(mesh.czm_submesh.mesh_bonded.gs.order)[2]),
        capture.separation_n,capture.separation_t,capture.traction_n,capture.traction_t,plasticity,bulk_kappa)
    return check_czm_snapshot(snapshot)
end

"""Public snapshot/CSV boundary, including real GP aggregation and explicit J2 mode."""
function check_czm_snapshot(snapshot::CZMSnapshot)
    check_czm_physical_output(snapshot)
    isfinite(snapshot.time_s) && snapshot.time_s >= 0.0 &&
        isfinite(snapshot.mechanical_dt_s) && snapshot.mechanical_dt_s > 0.0 ||
        throw(ArgumentError("snapshot requires finite local time and positive actual mechanical dt"))
    snapshot.czm_model in ("model1","mix") || throw(ArgumentError("invalid snapshot CZM model"))
    weights = check_cycle_gp_history(snapshot.gp_damage_states,length(snapshot.damage))
    snapshot.gp_points == NCweight(length(weights))[2] ||
        throw(ArgumentError("snapshot GP point order differs from its actual quadrature"))
    all(a -> size(a) == size(snapshot.gp_damage_states),
        (snapshot.gp_separation_n,snapshot.gp_separation_t,snapshot.gp_traction_n,snapshot.gp_traction_t)) ||
        throw(DimensionMismatch("snapshot GP current fields and history sizes differ"))
    all(a -> all(isfinite,a), (snapshot.gp_separation_n,snapshot.gp_separation_t,
        snapshot.gp_traction_n,snapshot.gp_traction_t)) ||
        throw(ArgumentError("snapshot GP output contains non-finite fields"))
    # Match the production COH loop's sequential weighted sums exactly.
    for (gp,element) in ((snapshot.gp_separation_n,snapshot.separation_n),
            (snapshot.gp_separation_t,snapshot.separation_t),
            (snapshot.gp_traction_n,snapshot.traction_n),(snapshot.gp_traction_t,snapshot.traction_t))
        for e in axes(gp,1)
            total = 0.0
            for g in eachindex(weights)
                total += weights[g] * gp[e,g]
            end
            total / sum(weights) == element[e] ||
                throw(ArgumentError("snapshot GP current response differs from its element output"))
        end
    end
    aggregate = aggregate_gp_damage_states(snapshot.gp_damage_states,weights)
    all(i -> aggregate[i].D == snapshot.damage[i] && aggregate[i].D_visc == snapshot.damage_visc[i] &&
        aggregate[i].fractured == snapshot.fractured[i], eachindex(aggregate)) ||
        throw(ArgumentError("snapshot GP aggregation differs from element output"))
    if snapshot.plasticity
        size(snapshot.bulk_kappa,1) > 0 && size(snapshot.bulk_kappa,2) == 4 ||
            throw(DimensionMismatch("J2 snapshot requires nonempty bulk history with four actual GPs"))
        all(x -> isfinite(x) && x >= 0.0, snapshot.bulk_kappa) ||
            throw(ArgumentError("snapshot has invalid accepted bulk kappa"))
    else
        snapshot.bulk_kappa === nothing || throw(ArgumentError("J2-off snapshot cannot publish bulk kappa"))
    end
    return snapshot
end

"""
    czm_max_separation_key(czm_model)

Return the physical result key corresponding to the separation measure that
drives the selected cohesive model.
"""
function czm_max_separation_key(czm_model::String)
    czm_model == "mix" && return "czm δ_max_eff [m]"
    czm_model == "model1" && return "czm δ_max_n [m]"
    throw(ArgumentError("unsupported CZM model: $czm_model"))
end
