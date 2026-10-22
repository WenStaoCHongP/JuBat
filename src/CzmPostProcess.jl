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

    D_vals = map(damage_states) do s
        isfinite(s.D) && 0<=s.D<=1 && isfinite(s.accumulated_damage) ||
            throw(ArgumentError("invalid damage statistics source"))
        s.D
    end
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
    ms.gp_damage_states === nothing && throw(ArgumentError(
        "nonzero cycle damage increment requires real GP history"))
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
    length(committed_states) == length(result.damage) || throw(DimensionMismatch(
        "committed CZM history must match result damage length"))
    for (i, state) in enumerate(committed_states)
        all(isfinite, (state.D, state.D_visc, state.δ_max_n, state.δ_max_t,
                       state.δ_max_eff, state.accumulated_damage)) || throw(ArgumentError(
            "committed CZM history contains non-finite values at element $i"))
        state.D == result.damage[i] || throw(ArgumentError(
            "committed CZM damage differs from result at element $i"))
    end
    new_variables = copy(variables)
    u_x = result.displacement[1:2:end]
    u_y = result.displacement[2:2:end]
    new_variables["czm displacement x"] = u_x
    new_variables["czm displacement y"] = u_y
    new_variables["czm damage"] = result.damage
    new_variables["czm traction normal"] = result.traction_n
    new_variables["czm traction tangent"] = result.traction_t
    new_variables["czm separation normal"] = result.separation_n
    new_variables["czm separation tangent"] = result.separation_t
    separation_eff = hypot.(max.(result.separation_n, 0.0), result.separation_t)
    new_variables["czm separation effective"] = separation_eff

    new_variables["czm D_max"] = maximum(result.damage)
    new_variables["czm D_mean"] = mean(result.damage)
    new_variables["czm δ_max_eff"] = [maximum(separation_eff)]
    new_variables["czm n_fractured"] = Float64(count(s -> s.fractured, committed_states))

    return new_variables
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
