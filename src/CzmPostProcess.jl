# ========================================================================
# CZM Post-processing, Statistics & Damage Management
# ========================================================================
# 从 CzmSolve.jl 拆分出来的后处理、统计和损伤管理函数。
# 这些函数无状态、纯计算，独立于求解器主循环。
# 注: Statistics 已在 JuBat.jl module 级别 using，此处无需重复。
# ========================================================================

"""
    get_damage_statistics(czm_mesh)
"""
function get_damage_statistics(czm_mesh::CohesiveMesh)
    n = czm_mesh.n_cohesive

    if n == 0
        return (max_D=0.0, mean_D=0.0, min_D=0.0, n_fractured=0, fraction_damaged=0.0, total_accumulated=0.0)
    end

    D_vals = [s.D for s in czm_mesh.damage_states]
    n_fractured = count(s -> s.fractured, czm_mesh.damage_states)
    n_damaged = count(d -> d > 0.01, D_vals)
    accumulated = [s.accumulated_damage for s in czm_mesh.damage_states]

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
    check_fracture_criterion(czm_mesh; threshold=0.99)
"""
function check_fracture_criterion(czm_mesh::CohesiveMesh; threshold::Float64=0.99)
    stats = get_damage_statistics(czm_mesh)

    is_fractured_avg = stats.mean_D >= threshold
    is_fractured_count = (stats.n_fractured / max(1, czm_mesh.n_cohesive)) > 0.5

    is_fractured = is_fractured_avg || is_fractured_count

    fracture_info = (
        is_fractured = is_fractured,
        criterion = is_fractured_avg ? :average_damage : (is_fractured_count ? :fractured_count : :none),
        stats = stats
    )

    return is_fractured, fracture_info
end

"""
    reset_damage_states(czm_mesh)

Return new CZM mesh with reset damage states.
"""
function reset_damage_states(czm_mesh::CohesiveMesh)
    new_damage_states = [DamageState() for _ in 1:czm_mesh.n_cohesive]
    new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, new_damage_states)
    return new_czm_mesh
end

"""
    accumulate_cycle_damage(czm_mesh, cycle_damage_increment)

Return new CZM mesh with accumulated damage.
"""
function accumulate_cycle_damage(czm_mesh::CohesiveMesh, cycle_damage_increment::Float64)
    new_damage_states = Vector{DamageState}(undef, czm_mesh.n_cohesive)
    for (i, state) in enumerate(czm_mesh.damage_states)
        new_state = DamageState()
        new_state.D = state.D
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

    new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, new_damage_states)
    return new_czm_mesh
end

"""
    czm_output_to_variables(czm_mesh, result, variables)
"""
function czm_output_to_variables(czm_mesh::CohesiveMesh, result::CZMResult, variables::Dict{String, Union{Array{Float64}, Float64}})
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

    stats = get_damage_statistics(czm_mesh)
    new_variables["czm max damage"] = stats.max_D
    new_variables["czm mean damage"] = stats.mean_D
    new_variables["czm fractured elements"] = Float64(stats.n_fractured)

    return new_variables
end
