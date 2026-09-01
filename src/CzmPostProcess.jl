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
function get_damage_statistics(damage_states::AbstractVector{<:AbstractDamageState})
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
function check_fracture_criterion(damage_states::AbstractVector{<:AbstractDamageState}; threshold::Float64=0.99)
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
    ms.damage_states = [DamageState() for _ in 1:length(ms.damage_states)]
    return ms
end

"""
    accumulate_cycle_damage!(ms, cycle_damage_increment)

原位累积循环损伤（未断裂单元 accumulated_damage += 增量；累积到 1 判断裂）。
"""
function accumulate_cycle_damage!(ms::MechState, cycle_damage_increment::Float64)
    new_damage_states = Vector{DamageState}(undef, length(ms.damage_states))
    for (i, state) in enumerate(ms.damage_states)
        new_state = DamageState()
        new_state.D = state.D
        new_state.D_visc = state.D_visc
        new_state.δ_max_n = state.δ_max_n
        new_state.δ_max_t = state.δ_max_t
        new_state.δ_max_eff = state.δ_max_eff
        new_state.fractured = state.fractured
        new_state.accumulated_damage = state.accumulated_damage

        if !new_state.fractured
            new_state.accumulated_damage += cycle_damage_increment
            if new_state.accumulated_damage >= 1.0
                new_state.D = 1.0
                new_state.fractured = true
            end
        end

        new_damage_states[i] = new_state
    end
    ms.damage_states = new_damage_states
    return ms
end

"""
    czm_output_to_variables(result, variables)
"""
function czm_output_to_variables(result::CZMResult, variables::Dict{String, Union{Array{Float64}, Float64}})
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
    new_variables["czm n_fractured"] = Float64(count(d -> d >= 0.99, result.damage))

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
