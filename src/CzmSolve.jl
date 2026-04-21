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

function clone_damage_states(damage_states::AbstractVector{<:AbstractDamageState})
    return [begin
        new_state = DamageState()
        new_state.D = s.D
        new_state.δ_max_n = s.δ_max_n
        new_state.δ_max_t = s.δ_max_t
        new_state.δ_max_eff = s.δ_max_eff
        new_state.fractured = s.fractured
        new_state.accumulated_damage = s.accumulated_damage
        new_state
    end for s in damage_states]
end

function clone_czm_mesh_with_damage(czm_mesh::CohesiveMesh, damage_states::AbstractVector{<:AbstractDamageState})
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
    return new_czm_mesh
end

function apply_czm_dirichlet!(u::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64}, bc_vals::AbstractVector{Float64})
    for (dof, val) in zip(bc_dofs, bc_vals)
        u[dof] = val
    end
    return u
end

function zero_czm_bc_entries!(v::AbstractVector{Float64}, bc_dofs::AbstractVector{Int64})
    for dof in bc_dofs
        v[dof] = 0.0
    end
    return v
end

function fill_czm_result!(result::CZMResult, u::Vector{Float64}, damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, tractions::Vector{Tuple{Float64, Float64}})
    result.displacement = u
    for i in eachindex(damage_states)
        result.damage[i] = damage_states[i].D
        result.separation_n[i] = separations[i][1]
        result.separation_t[i] = separations[i][2]
        result.traction_n[i] = tractions[i][1]
        result.traction_t[i] = tractions[i][2]
    end
    return result
end

function build_arc_length_augmented_matrix(K_bc::SparseMatrixCSC{Float64, Int64}, load_vector::Vector{Float64}, delta_u::Vector{Float64}, delta_lambda::Float64, arc_length_alpha::Float64)
    ndof = length(load_vector)
    A = spzeros(Float64, ndof + 1, ndof + 1)
    A[1:ndof, 1:ndof] = K_bc
    for i in 1:ndof
        A[i, ndof + 1] = -load_vector[i]
        A[ndof + 1, i] = 2.0 * delta_u[i]
    end
    A[ndof + 1, ndof + 1] = 2.0 * arc_length_alpha^2 * delta_lambda
    return A
end

function solve_czm_basic_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, cache::Union{Nothing, CZMAssemblyCache}=nothing)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(u_prev)
        u_start = copy(u)
        damage_states = czm_mesh.damage_states
        damage_start = clone_damage_states(damage_states)

        if cache !== nothing
            bc_dofs = cache.bc_dofs
            bc_vals = cache.bc_vals
        else
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
        end

        F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
        geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
        ws_basic = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(ndof, n_coh)

        total_iter = 0
        R_norm_0 = 1.0
        R_norm = Inf
        converged = false
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, E_eff, ν_eff, cohesive_params;
                damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws_basic)

            R = F_ext + F_thermo_chem - f_int_total
            apply_czm_dirichlet!(R, bc_dofs, zeros(length(bc_dofs)))

            R_norm = norm(R)
            if iter == 1
                R_norm_0 = max(R_norm, 1e-10)
            end
            rel_norm = R_norm / R_norm_0

            if R_norm < tol || rel_norm < tol
                converged = true
                damage_states = update_damage(damage_states, separations, cohesive_params)
                break
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))

            Δu = try
                K_bc \ R_bc
            catch
                break
            end

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            α = 1.0
            ls_accepted = false
            for _ in 1:8
                u_trial = u + α * Δu
                apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

                _, f_int_trial, _, _ = assemble_coupled_system(
                    czm_mesh, u_trial, E_eff, ν_eff, cohesive_params;
                    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                    geom_cache=geom_cache, ws=ws_basic)

                R_trial = F_ext + F_thermo_chem - f_int_trial
                apply_czm_dirichlet!(R_trial, bc_dofs, zeros(length(bc_dofs)))

                R_trial_norm = norm(R_trial)
                if !isnan(R_trial_norm) && R_trial_norm < R_norm
                    ls_accepted = true
                    break
                end

                α *= 0.5
            end

            if !ls_accepted
                break
            end

            u = u + α * Δu
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
        end

        if !converged
            u = u_start
            damage_states = damage_start
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(
            czm_mesh, u, E_eff, ν_eff, cohesive_params;
            damage_states=damage_states, K_bulk_cached=K_bulk_cached,
            geom_cache=geom_cache, ws=ws_basic)

        R = F_ext + F_thermo_chem - f_int_total
        apply_czm_dirichlet!(R, bc_dofs, zeros(length(bc_dofs)))
        R_norm = norm(R)

        final_tol = tol * 100.0

        result.converged = converged
        result.iterations = total_iter
        result.residual_norm = R_norm
        result.displacement = u
        fill_czm_result!(result, u, damage_states, separations, tractions)

        new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
        return result, new_czm_mesh
    end

    function solve_czm_arc_length_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, cache::Union{Nothing, CZMAssemblyCache}=nothing)
        nnode = czm_mesh.nnode
        ndof = 2 * nnode
        n_coh = czm_mesh.n_cohesive

        result = CZMResult(ndof, n_coh)
        u = copy(u_prev)
        damage_states = czm_mesh.damage_states

        if cache !== nothing
            bc_dofs = cache.bc_dofs
            bc_vals = cache.bc_vals
        else
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
        end

        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
        K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
        geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
        ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(ndof, n_coh)

        total_iter = 0
        load_progress = 0.0
        load_step = 0
        step_size = 1.0 / max(1, n_load_steps)
        step_size_min = step_size / 128.0
        step_size_max = step_size
        last_residual = Inf
        converged_substep = false
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        while load_progress < 1.0 - 1e-12
            load_step += 1
            load_start = load_progress
            target_progress = min(1.0, load_start + step_size)

            u_start = copy(u)
            damage_start = clone_damage_states(damage_states)
            converged_substep = false
            last_residual = Inf

            K_total, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, E_eff, ν_eff, cohesive_params;
                damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws)

            F_load_bc = copy(F_thermo_chem_total)
            zero_czm_bc_entries!(F_load_bc, bc_dofs)

            F_thermo_chem = load_start * F_thermo_chem_total
            R = F_ext + F_thermo_chem - f_int_total
            apply_czm_dirichlet!(R, bc_dofs, zeros(length(bc_dofs)))
            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))

            tangent = try
                K_bc \ F_load_bc
            catch
                nothing
            end

            if tangent === nothing || any(isnan, tangent) || any(isinf, tangent)
                break
            end

            delta_lambda_pred = target_progress - load_start
            if delta_lambda_pred <= 0.0
                delta_lambda_pred = step_size
            end

            delta_u_pred = tangent * delta_lambda_pred
            arc_target = sqrt(sum(abs2, delta_u_pred) + (arc_length_alpha * delta_lambda_pred)^2)
            if !isfinite(arc_target) || arc_target <= 0.0
                break
            end

            u = u_start + delta_u_pred
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            load_progress = target_progress

            for iter in 1:max_iter
                total_iter += 1

                F_thermo_chem = load_progress * F_thermo_chem_total
                K_total, f_int_total, separations, tractions = assemble_coupled_system(
                    czm_mesh, u, E_eff, ν_eff, cohesive_params;
                    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                    geom_cache=geom_cache, ws=ws)

                R = F_ext + F_thermo_chem - f_int_total
                apply_czm_dirichlet!(R, bc_dofs, zeros(length(bc_dofs)))

                delta_u = u - u_start
                delta_lambda = load_progress - load_start
                arc_constraint = dot(delta_u, delta_u) + (arc_length_alpha * delta_lambda)^2 - arc_target^2
                residual_norm = sqrt(norm(R)^2 + arc_constraint^2)
                last_residual = residual_norm

                substep_tol = tol * 10.0
                if norm(R) < substep_tol && abs(arc_constraint) < substep_tol
                    converged_substep = true
                    damage_states = update_damage(damage_states, separations, cohesive_params)
                    load_progress = min(load_progress, 1.0)
                    step_size = min(step_size * 1.25, step_size_max)
                    break
                end

                K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
                A = build_arc_length_augmented_matrix(K_bc, F_load_bc, delta_u, delta_lambda, arc_length_alpha)
                rhs = vcat(-R_bc, -arc_constraint)

                sol = try
                    A \ rhs
                catch
                    nothing
                end

                if sol === nothing
                    break
                end

                delta_u_corr = sol[1:ndof]
                delta_lambda_corr = sol[end]
                if any(isnan, delta_u_corr) || any(isinf, delta_u_corr) || isnan(delta_lambda_corr) || isinf(delta_lambda_corr)
                    break
                end

                α = 1.0
                accepted = false
                for _ in 1:8
                    u_trial = u + α * delta_u_corr
                    lambda_trial = load_progress + α * delta_lambda_corr

                    if lambda_trial < load_start - 1e-12 || lambda_trial > 1.0 + 1e-12
                        α *= 0.5
                        continue
                    end

                    apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

                    F_trial = lambda_trial * F_thermo_chem_total
                    _, f_int_trial, _, _ = assemble_coupled_system(
                        czm_mesh, u_trial, E_eff, ν_eff, cohesive_params;
                        damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                        geom_cache=geom_cache, ws=ws)

                    R_trial = F_ext + F_trial - f_int_trial
                    apply_czm_dirichlet!(R_trial, bc_dofs, zeros(length(bc_dofs)))
                    delta_u_trial = u_trial - u_start
                    delta_lambda_trial = lambda_trial - load_start
                    arc_trial = dot(delta_u_trial, delta_u_trial) + (arc_length_alpha * delta_lambda_trial)^2 - arc_target^2
                    merit_trial = sqrt(norm(R_trial)^2 + arc_trial^2)

                    if !isnan(merit_trial) && merit_trial < residual_norm
                        u = u_trial
                        load_progress = lambda_trial
                        accepted = true
                        break
                    end

                    α *= 0.5
                end

                if !accepted
                    break
                end
            end

            if !converged_substep
                u = u_start
                damage_states = damage_start
                load_progress = load_start
                step_size *= 0.5

                if step_size < step_size_min
                    @warn "CZM arc-length stepping stalled" load_progress=load_progress target_progress=target_progress residual=last_residual step_size=step_size
                    break
                end

                @debug "Arc-length substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_residual step_size=step_size
                continue
            end
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(
            czm_mesh, u, E_eff, ν_eff, cohesive_params;
            damage_states=damage_states, K_bulk_cached=K_bulk_cached,
            geom_cache=geom_cache, ws=ws)

        F_thermo_chem = load_progress * F_thermo_chem_total
        R = F_ext + F_thermo_chem - f_int_total
        apply_czm_dirichlet!(R, bc_dofs, zeros(length(bc_dofs)))
        R_norm = norm(R)

        final_tol = tol * 100.0
        result.converged = load_progress >= 1.0 - 1e-12 && converged_substep
        result.iterations = total_iter
        result.residual_norm = R_norm
        result.displacement = u
        fill_czm_result!(result, u, damage_states, separations, tractions)

        new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
        return result, new_czm_mesh
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
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, u0::Union{Vector{Float64},Nothing}=nothing, n_load_steps::Int=10, cache::Union{Nothing, CZMAssemblyCache}=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive

    result = CZMResult(ndof, n_coh)

    u = u0 === nothing ? zeros(Float64, ndof) : copy(u0)
    damage_states = czm_mesh.damage_states

    # 使用缓存的 BC 或重新计算
    if cache !== nothing
        bc_dofs = cache.bc_dofs
        bc_vals = cache.bc_vals
    else
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
    end

    F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)

    # 提取缓存
    K_bulk_cached = cache !== nothing ? cache.K_bulk : nothing
    geom_cache = cache !== nothing ? cache.cohesive_geom : nothing
    ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(2 * nnode, czm_mesh.n_cohesive)

    total_iter = 0
    load_progress = 0.0
    load_step = 0
    step_size = 1.0 / max(1, n_load_steps)
    step_size_min = step_size / 128.0
    step_size_max = step_size

    last_R_norm = Inf
    converged_substep = false

    while load_progress < 1.0 - 1e-12
        load_step += 1
        target_progress = min(1.0, load_progress + step_size)
        F_thermo_chem = target_progress * F_thermo_chem_total

        u_start = copy(u)
        damage_start = clone_damage_states(damage_states)
        converged_substep = false
        last_R_norm = Inf

        for iter in 1:max_iter
            total_iter += 1

            K_total, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, E_eff, ν_eff, cohesive_params;
                damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                geom_cache=geom_cache, ws=ws)

            R = F_ext + F_thermo_chem - f_int_total

            for (dof, val) in zip(bc_dofs, bc_vals)
                R[dof] = val - u[dof]
            end

            R_norm = norm(R)
            last_R_norm = R_norm
            substep_tol = tol * 10.0

            if R_norm < substep_tol
                converged_substep = true
                damage_states = update_damage(damage_states, separations, cohesive_params)
                load_progress = target_progress
                step_size = min(step_size * 1.25, step_size_max)
                break
            end

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=zeros(length(bc_dofs)))
            Δu = K_bc \ R_bc

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            α = 1.0
            ls_accepted = false
            for _ in 1:8
                u_trial = u + α * Δu
                for (dof, val) in zip(bc_dofs, bc_vals)
                    u_trial[dof] = val
                end

                _, f_int_trial, _, _ = assemble_coupled_system(
                    czm_mesh, u_trial, E_eff, ν_eff, cohesive_params;
                    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
                    geom_cache=geom_cache, ws=ws)

                R_trial = F_ext + F_thermo_chem - f_int_trial
                for (dof, val) in zip(bc_dofs, bc_vals)
                    R_trial[dof] = val - u_trial[dof]
                end

                R_trial_norm = norm(R_trial)
                if !isnan(R_trial_norm) && R_trial_norm < R_norm
                    ls_accepted = true
                    break
                end

                α *= 0.5
            end

            if !ls_accepted
                break
            end

            u = u + α * Δu

            for (dof, val) in zip(bc_dofs, bc_vals)
                u[dof] = val
            end
        end

        if !converged_substep
            u = u_start
            damage_states = damage_start
            step_size *= 0.5

            if step_size < step_size_min
                @warn "CZM adaptive load stepping stalled" load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
                break
            end

            @debug "Load substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
            continue
        end
    end

    K_total, f_int_total, separations, tractions = assemble_coupled_system(
        czm_mesh, u, E_eff, ν_eff, cohesive_params;
        damage_states=damage_states, K_bulk_cached=K_bulk_cached,
        geom_cache=geom_cache, ws=ws)

    F_thermo_chem = load_progress * F_thermo_chem_total
    R = F_ext + F_thermo_chem - f_int_total
    for (dof, val) in zip(bc_dofs, bc_vals)
        R[dof] = val - u[dof]
    end
    R_norm = norm(R)

    final_tol = tol * 100.0

    result.converged = load_progress >= 1.0 - 1e-12 && R_norm < final_tol
    result.iterations = total_iter
    result.residual_norm = R_norm
    result.displacement = u

    if result.converged
        damage_states = update_damage(damage_states, separations, cohesive_params)
    end

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
function solve_czm_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64}, E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive, param, u_prev::Vector{Float64}; α_eff::Float64=0.0, β_n::Float64=0.0, β_p::Float64=0.0, dT_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_n_elem::Union{Vector{Float64}, Nothing}=nothing, Δsoc_p_elem::Union{Vector{Float64}, Nothing}=nothing, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, iter_method::String="load_substep", cache::Union{Nothing, CZMAssemblyCache}=nothing)
    method = lowercase(iter_method)

    if method == "load_substep"
        return newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, u0=u_prev, n_load_steps=n_load_steps, cache=cache)
    elseif method == "basic"
        return solve_czm_basic_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, cache=cache)
    elseif method == "arc_length" || method == "arclength" || method == "arc-length"
        return solve_czm_arc_length_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
            max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, cache=cache)
    else
        error("Unknown CZM iteration method: $(iter_method). Use 'basic', 'load_substep', or 'arc_length'.")
    end
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

# ========================================================================
# CZM 损伤更新（由 Solve.jl 和 CycleSolver.jl 共用）
# ========================================================================

"""
    compute_czm_effective_params(case)

计算 CZM 求解所需的有效材料参数，全部基于 `SetParams.NormaliseParam`
后的归一化参数。

# 返回
- `E_eff`: 有效弹性模量 [-]
- `ν_eff`: 有效泊松比 [-]
- `α_eff`: 有效热膨胀系数 [-]
- `β_n`: 负极扩散应变系数 [-]
- `β_p`: 正极扩散应变系数 [-]
"""
function compute_czm_effective_params(case::Case)
    param = case.param

    # 有效弹性模量（厚度加权平均，均为归一化量）
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 有效泊松比（厚度加权平均）
    ν_eff = (param.NE.nu * param.NE.thickness + param.PE.nu * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 有效热膨胀系数（厚度加权平均，已按 T_ref 归一化）
    α_eff = (param.NE.alphaT * param.NE.thickness + param.PE.alphaT * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)

    # 扩散应变系数 β = (Ω * c_s,max) / 3 已在 SetParams 中完成归一化
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0

    return E_eff, ν_eff, α_eff, β_n, β_p
end

"""
    compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

计算 CZM 损伤计算所需的单元级应变输入。

# 返回
- `dT_elem`: 每个单元的温度变化 [K]
- `Δsoc_n_elem`: 每个单元的负极 SOC 变化 [-]
- `Δsoc_p_elem`: 每个单元的正极 SOC 变化 [-]
"""
function compute_czm_strain_inputs(case::Case, variables::Dict, czm_mesh, T_nodes_carry)
    ne = size(czm_mesh.bulk_element, 1)
    param = case.param

    # 参考 SOC（归一化值）
    soc_ref_n = param.NE.cs0
    soc_ref_p = param.PE.cs0

    # 初始化输出数组
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)

    # 提取温度场（无量纲温度 T* = T / T_ref）
    if length(T_nodes_carry) >= czm_mesh.nnode
        for e in 1:ne
            nodes = czm_mesh.bulk_element[e, :]
            T_elem_nd = 0.0
            valid_nodes = 0
            for n in nodes
                if n <= length(T_nodes_carry)
                    T_elem_nd += T_nodes_carry[n]
                    valid_nodes += 1
                end
            end
            if valid_nodes > 0
                T_elem_nd /= valid_nodes
                dT_elem[e] = T_elem_nd - param.cell.T0
            end
        end
    end

    # 提取 SOC 分布（如果 variables 中有）
    soc_n_elem = variables["thermal2D element soc_n"]
    soc_p_elem = variables["thermal2D element soc_p"]

    # 处理数组维度（可能是 ne×1 或 ne×num）
    if isa(soc_n_elem, AbstractMatrix)
        soc_n_elem = soc_n_elem[:, end]
        soc_p_elem = soc_p_elem[:, end]
    end

    for e in 1:min(ne, length(soc_n_elem))
        Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
        Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
    end

    return dT_elem, Δsoc_n_elem, Δsoc_p_elem
end

"""
    update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)

更新 CZM 网格的损伤状态。

使用牛顿-拉弗森迭代求解力学平衡方程，通过载荷子步法处理软化收敛问题。

# 参数
- `czm_mesh`: CZM 网格对象
- `czm_params`: CZM 参数（cohesive）
- `case`: Case 对象
- `variables`: 当前时间步的变量字典
- `T_nodes_carry`: 当前温度场
- `u_czm_prev`: 上一步的 CZM 位移场

# 返回
- `u_czm`: 更新后的 CZM 位移场
- `converged`: 是否收敛
"""
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
    param_dim = case.param_dim
    param = case.param

    # 同步CZM模型选项（model1 or mix）
    czm_params.czm_model = case.opt.czm_model

    # 计算有效材料参数
    E_eff, ν_eff, α_eff, β_n, β_p = compute_czm_effective_params(case)

    # 构建或复用 CZM 缓存
    cache = ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)

    # 计算应变输入
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

    # 外力向量（一般为零）
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)

    # 诊断：检测输入异常（帮助定位 NaN 来源）
    has_nan_T = any(isnan, T_nodes_carry)
    has_nan_soc_n = any(isnan, variables["thermal2D element soc_n"])
    has_nan_soc_p = any(isnan, variables["thermal2D element soc_p"])
    if has_nan_T || has_nan_soc_n || has_nan_soc_p || any(isnan, dT_elem) || any(isnan, Δsoc_n_elem) || any(isnan, Δsoc_p_elem)
        @warn "CZM inputs contain NaN" has_nan_T=has_nan_T has_nan_soc_n=has_nan_soc_n has_nan_soc_p=has_nan_soc_p n_nan_dT=count(isnan, dT_elem) n_nan_soc_n=count(isnan, Δsoc_n_elem) n_nan_soc_p=count(isnan, Δsoc_p_elem)
    end

    # 初始化位移（如果没有上一步的值）
    if u_czm_prev === nothing || length(u_czm_prev) != ndof
        u_czm_prev = zeros(Float64, ndof)
    elseif any(isnan, u_czm_prev)
        @warn "CZM u_czm_prev contains NaN, resetting to zeros"
        u_czm_prev = zeros(Float64, ndof)
    end

    # 调用 CZM 求解器（可选迭代方式）
    iter_method = case.opt.czm_iter_method
    max_iter = case.opt.czm_max_iter
    tol = case.opt.czm_tol
    n_load_steps = case.opt.czm_load_steps
    arc_length_alpha = case.opt.czm_arc_length_alpha

    result, updated_czm_mesh = solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, czm_params, param, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, iter_method=iter_method,
        cache=cache
    )

    # 诊断：检查求解结果是否异常
    has_nan_disp = any(isnan, result.displacement)
    has_nan_damage = any(ds -> isnan(ds.D), updated_czm_mesh.damage_states)
    if has_nan_disp || has_nan_damage || !result.converged
        @warn "CZM solve issue" converged=result.converged iterations=result.iterations residual=round(result.residual_norm; digits=4) has_nan_disp=has_nan_disp has_nan_damage=has_nan_damage
    end

    # Only commit damage states when the nonlinear solve converged.
    # This avoids propagating a partially converged or diverged state into the next time step.
    if result.converged
        czm_mesh.damage_states = updated_czm_mesh.damage_states
    end

    return result.displacement, result.converged
end
