mutable struct CZMResult
    displacement::Vector{Float64}
    damage::Vector{Float64}
    traction_n::Vector{Float64}
    traction_t::Vector{Float64}
    separation_n::Vector{Float64}
    separation_t::Vector{Float64}
    converged::Bool
    iterations::Int64
    residual_norm::Float64
    
    CZMResult(ndof::Int, n_coh::Int) = new(
        zeros(ndof), zeros(n_coh), zeros(n_coh), zeros(n_coh),
        zeros(n_coh), zeros(n_coh), false, 0, Inf)
end

# ========================================================================
# 7. Newton-Raphson solver
# ========================================================================

"""
    newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param;α_eff=0.0, β_n=0.0, β_p=0.0,dT_elem=nothing, Δsoc_n_elem=nothing, Δsoc_p_elem=nothing,max_iter=50, tol=1e-8, u0=nothing, n_load_steps=10)

Newton-Raphson nonlinear solver with load substeps.

# Returns
- `result`: CZMResult
- `new_czm_mesh`: updated CZM mesh with damage states
"""
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, u0::Union{Vector{Float64},Nothing}=nothing, n_load_steps::Int=10)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)

    u = copy(u0)
    damage_states = czm_mesh.damage_states

    bc_nodes, inner_count, outer_count = identify_bc_nodes_czm(czm_mesh, param)
    bc_dofs = Int64[]
    bc_vals = Float64[]

    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1)
            push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node)
            push!(bc_vals, 0.0)
        end
    end

    F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)

    total_iter = 0

    for load_step in 1:n_load_steps
        load_factor = load_step / n_load_steps
        F_thermo_chem = load_factor * F_thermo_chem_total
        converged_substep = false

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)

            R = F_ext + F_thermo_chem - f_int_total

            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            substep_tol = tol * 10.0

            if R_norm < substep_tol
                converged_substep = true
                damage_states = update_damage(damage_states, separations, cohesive_params)
                break
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
            Δu = K_bc \ R_bc

            Δu_norm = norm(Δu)
            max_Δu = 1e-6
            if Δu_norm > max_Δu
                Δu = Δu * (max_Δu / Δu_norm)
            end

            u = u + Δu

            for (dof, val) in zip(bc_dofs, bc_vals)
                u[dof] = val
            end
        end

        if !converged_substep
            @debug "Load substep $load_step did not fully converge, continuing..."
        end
    end

    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)

    R = F_ext + F_thermo_chem_total - f_int_total
    for (dof, val) in zip(bc_dofs, bc_vals)
        R[dof] = val - u[dof]
    end
    R_norm = norm(R)

    final_tol = tol * 100.0

    result.converged = R_norm < final_tol
    result.iterations = total_iter
    result.residual_norm = R_norm
    result.displacement = u

    damage_states = update_damage(damage_states, separations, cohesive_params)

    for i in 1:n_coh
        result.damage[i] = damage_states[i].D
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end

    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = damage_states

    return result, new_czm_mesh
end

"""
    solve_czm_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev; ...)

Solve a single CZM step with selectable iteration method.
"""
function solve_czm_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, iter_method::String="load_substep")
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)
    u = copy(u_prev)
    damage_states = czm_mesh.damage_states

    bc_nodes, inner_count, outer_count = identify_bc_nodes_czm(czm_mesh, param)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1)
            push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node)
            push!(bc_vals, 0.0)
        end
    end

    F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)

    method = lowercase(iter_method)

    if method == "load_substep"
        result, new_czm_mesh = newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, u0=u, n_load_steps=n_load_steps)
        return result, new_czm_mesh
    elseif method == "basic"
        R_norm_0 = 1.0
        R_norm = 0.0

        for iter in 1:max_iter
            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)

            R = F_ext + F_thermo_chem - f_int_total
            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            if iter == 1
                R_norm_0 = max(R_norm, 1e-10)
            end
            rel_norm = R_norm / R_norm_0

            if R_norm < tol || rel_norm < tol
                result.converged = true
                result.iterations = iter
                result.residual_norm = R_norm
                result.displacement = u

                damage_states = update_damage(damage_states, separations, cohesive_params)

                for i in 1:n_coh
                    result.damage[i] = damage_states[i].D
                    result.separation_n[i] = separations[i][1]
                    result.separation_t[i] = separations[i][2]
                    result.traction_n[i] = tractions[i][1]
                    result.traction_t[i] = tractions[i][2]
                end

                new_czm_mesh = CohesiveMesh()
                new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
                new_czm_mesh.node = czm_mesh.node
                new_czm_mesh.nnode = czm_mesh.nnode
                new_czm_mesh.bulk_element = czm_mesh.bulk_element
                new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
                new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
                new_czm_mesh.n_layers = czm_mesh.n_layers
                new_czm_mesh.node_map = czm_mesh.node_map
                new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
                new_czm_mesh.damage_states = damage_states
                return result, new_czm_mesh
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
            Δu = K_bc \ R_bc
            u += Δu
            for (dof, val) in zip(bc_dofs, bc_vals)
                u[dof] = val
            end
        end

        @warn "Newton-Raphson (basic) did not converge" max_iter=max_iter residual=R_norm
        result.converged = false
        result.iterations = max_iter
        result.displacement = u
    elseif method == "arc_length" || method == "arclength" || method == "arc-length"
        F_total = F_ext + F_thermo_chem
        if norm(F_total) < 1e-20
            result.converged = true
            result.iterations = 0
            result.residual_norm = 0.0
            result.displacement = u
        else
            λ = 0.0
            total_iter = 0

            for step in 1:n_load_steps
                K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)

                K_bc, _ = apply_bc_czm(K_total, zeros(Float64, ndof); bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
                u_dot = K_bc \ F_total
                u_dot_norm2 = dot(u_dot, u_dot)

                Δλ_pred = 1.0 / max(1, n_load_steps)
                Δs = sqrt(u_dot_norm2 * Δλ_pred^2 + (arc_length_alpha * Δλ_pred)^2)

                u_step0 = copy(u)
                λ_step0 = λ

                u = u + Δλ_pred * u_dot
                λ = λ + Δλ_pred

                converged_step = false
                for iter in 1:max_iter
                    total_iter += 1

                    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)

                    R = λ * F_total - f_int_total
                    for (dof, val) in zip(bc_dofs, bc_vals)
                        R[dof] = val - u[dof]
                    end

                    R_norm = norm(R)
                    if R_norm < tol
                        converged_step = true
                        break
                    end

                    K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
                    Δu_R = K_bc \ R_bc
                    Δu_F = K_bc \ F_total

                    Δu_total = u - u_step0
                    Δλ_total = λ - λ_step0

                    A = dot(Δu_F, Δu_F) + arc_length_alpha^2
                    B = 2.0 * (dot(Δu_total, Δu_F) + arc_length_alpha^2 * Δλ_total + dot(Δu_R, Δu_F))
                    C = dot(Δu_total, Δu_total) + (arc_length_alpha * Δλ_total)^2 +
                        2.0 * dot(Δu_total, Δu_R) + dot(Δu_R, Δu_R) - Δs^2

                    disc = B^2 - 4.0 * A * C
                    if disc < 0
                        disc = 0.0
                    end
                    sqrt_disc = sqrt(disc)

                    δλ_1 = (-B + sqrt_disc) / (2.0 * A)
                    δλ_2 = (-B - sqrt_disc) / (2.0 * A)
                    δλ = abs(δλ_1 - Δλ_pred) <= abs(δλ_2 - Δλ_pred) ? δλ_1 : δλ_2

                    δu = Δu_R + δλ * Δu_F

                    u += δu
                    λ += δλ

                    for (dof, val) in zip(bc_dofs, bc_vals)
                        u[dof] = val
                    end
                end

                if !converged_step
                    @debug "Arc-length substep $step did not fully converge, continuing..."
                end

                damage_states = update_damage(damage_states, separations, cohesive_params)
            end

            K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)
            R = λ * F_total - f_int_total
            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end
            R_norm = norm(R)

            result.converged = R_norm < tol * 100.0
            result.iterations = total_iter
            result.residual_norm = R_norm
            result.displacement = u

            damage_states = update_damage(damage_states, separations, cohesive_params)

            for i in 1:n_coh
                result.damage[i] = damage_states[i].D
                result.separation_n[i] = separations[i][1]
                result.separation_t[i] = separations[i][2]
                result.traction_n[i] = tractions[i][1]
                result.traction_t[i] = tractions[i][2]
            end
        end
    else
        error("Unknown CZM iteration method: $(iter_method). Use 'basic', 'load_substep', or 'arc_length'.")
    end

    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = damage_states

    return result, new_czm_mesh
end

# ========================================================================
# 8. Damage statistics and fracture criteria
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

    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = new_damage_states

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

    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = new_damage_states

    return new_czm_mesh
end

# ========================================================================
# 9. Post-processing
# ========================================================================

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
