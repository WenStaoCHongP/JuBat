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

"""
    extract_bc_dofs(czm_mesh, param; cache=nothing)

从 czm_mesh 提取 Dirichlet BC 的自由度列表和对应值。
优先使用缓存中的 bc_dofs/bc_vals，否则从 identify_bc_nodes_czm 重新计算。
"""
function extract_bc_dofs(czm_mesh::CohesiveMesh, param; cache::Union{Nothing, CZMAssemblyCache}=nothing)
    if cache !== nothing
        return cache.bc_dofs, cache.bc_vals
    end
    bc_nodes, _, _ = identify_bc_nodes_czm(czm_mesh, param)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        end
    end
    return bc_dofs, bc_vals
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

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

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

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

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

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)

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

    new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)

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
