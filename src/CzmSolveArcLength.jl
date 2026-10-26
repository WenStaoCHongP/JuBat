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

function spherical_arc_length_correction(delta_u_bar::AbstractVector{<:Real},
        delta_lambda_bar::Real, delta_u_R::AbstractVector{<:Real},
        delta_u_F::AbstractVector{<:Real}, arc_length_alpha::Real, arc_radius::Real)
    isfinite(arc_radius) && arc_radius > 0.0 || throw(ArgumentError(
        "spherical arc correction requires finite positive radius, got $arc_radius"))
    g = dot(delta_u_bar, delta_u_bar) +
        arc_length_alpha^2 * delta_lambda_bar^2 - arc_radius^2
    denominator = 2.0 * dot(delta_u_bar, delta_u_F) +
                  2.0 * arc_length_alpha^2 * delta_lambda_bar
    denominator_scale = 2.0 * (norm(delta_u_bar) * norm(delta_u_F) +
                                arc_length_alpha^2 * abs(delta_lambda_bar))
    abs(denominator) > eps(Float64) * max(1.0, denominator_scale) || error(
        "spherical arc correction is singular (denominator=$denominator, g=$g)")
    delta_lambda = (-g - 2.0 * dot(delta_u_bar, delta_u_R)) / denominator
    delta_u = delta_u_R .+ delta_lambda .* delta_u_F
    all(isfinite, delta_u) && isfinite(delta_lambda) || error(
        "spherical arc correction produced a non-finite update")
    return delta_u, delta_lambda, g
end

    # A cylinder controls displacement length. Its direction therefore uses
    # displacement cosine only, independently of the spherical lambda weight.
    # The normalized dot-product rounding bound is separate from Newton tol.
    cylindrical_arc_direction_tolerance(n::Int) = 8.0 * eps(Float64) * max(n, 1)

    function cylindrical_arc_direction_cosine(path::Vector{Float64}, direction::Vector{Float64})
        path_norm, direction_norm = norm(path), norm(direction)
        isfinite(path_norm) && isfinite(direction_norm) &&
            path_norm > 0.0 && direction_norm > 0.0 || return nothing
        return dot(path ./ path_norm, direction ./ direction_norm)
    end

    function cylindrical_arc_forward(path::Vector{Float64}, direction::Vector{Float64})
        cosine = cylindrical_arc_direction_cosine(path, direction)
        return cosine !== nothing && cosine > cylindrical_arc_direction_tolerance(length(path))
    end

    function cylindrical_arc_predictor(tangent::Vector{Float64}, load_start::Float64,
            step_size::Float64, previous_u_tangent::Union{Nothing,Vector{Float64}})
        delta_lambda = step_size
        if previous_u_tangent !== nothing
            cosine = cylindrical_arc_direction_cosine(tangent, previous_u_tangent)
            cosine !== nothing && abs(cosine) >
                cylindrical_arc_direction_tolerance(length(tangent)) || return nothing
            cosine < 0.0 && (delta_lambda = -delta_lambda)
        end
        delta_u = tangent * delta_lambda
        return (delta_lambda=delta_lambda, delta_u=delta_u,
                lambda=load_start + delta_lambda, radius_sq=dot(delta_u, delta_u))
    end

    function cylindrical_arc_converged(R_norm::Float64, delta_u_sq::Float64,
            radius_sq::Float64, tol::Float64)
        isfinite(radius_sq) && radius_sq > 0.0 || return false
        return isfinite(R_norm) && R_norm < 10.0 * tol &&
               abs(delta_u_sq - radius_sq) / radius_sq < tol
    end

    function cylindrical_arc_endpoint_action(R_norm::Float64,
            delta_u_sq::Float64, radius_sq::Float64, tol::Float64,
            load_progress::Float64)
        cylindrical_arc_converged(R_norm, delta_u_sq, radius_sq, tol) &&
            R_norm < tol || return :continue
        if load_progress > 1.0 + 1e-12
            return :target
        end
        return :accept
    end

    arc_wall_expired(start_ns::UInt64, limit_seconds::Float64) =
        limit_seconds > 0.0 && (time_ns() - start_ns) * 1e-9 >= limit_seconds

    function arc_observe(event::Function)
        observer = mechanical_step_observer[]
        observer === nothing || observer(event())
    end

    arc_merit(residual::Float64, constraint::Float64, radius_sq::Float64,
        tol::Float64) = 0.5 * ((residual / tol)^2 + (constraint / (tol * radius_sq))^2)

    function arc_gp_event_values(gaps, mesh, param, history, model)
        values = Array{Float64}(undef, size(history)..., 5)
        for i in axes(history,1), g in axes(history,2)
            ip = collector_params(param, mesh.cohesive_elements[i].interface_type)
            dn, dt = gaps.separation_n[i,g], gaps.separation_t[i,g]
            env = bilinear_opening_envelope(dn, dt, ip, model)
            D = history[i,g].D
            threshold = env.a * env.c / (env.c - D * (env.c - env.a))
            values[i,g,1] = env.e - threshold
            values[i,g,2] = env.e - env.a
            values[i,g,3] = env.e - env.c
            values[i,g,4] = dn
            values[i,g,5] = model == "mix" ? dt : 1.0
        end
        return values
    end

    function arc_gp_history_boundary(gaps, mesh, param, history, model)
        boundary = falses(size(history))
        q_n, q_t = zeros(size(history)), zeros(size(history))
        for i in axes(history,1), g in axes(history,2)
            ip = collector_params(param, mesh.cohesive_elements[i].interface_type)
            dn, dt = gaps.separation_n[i,g], gaps.separation_t[i,g]
            envelope = bilinear_opening_envelope(dn, dt, ip, model)
            if 0.0 < history[i,g].D < 1.0 && envelope.a < envelope.e < envelope.c &&
                    bilinear_envelope_damage(envelope) == history[i,g].D
                boundary[i,g] = true
                _, _, trial = cohesive_response(BilinearCohesiveLaw(ip, model), (dn, dt),
                    history[i,g], CohesiveStepData(1.0, true))
                q, _ = bilinear_damage_gradient(dn, dt, history[i,g], trial, ip, model, 1.0;
                    history_loading=true)
                q_n[i,g], q_t[i,g] = q
            end
        end
        return (;boundary, q_n, q_t)
    end

    arc_same_point(a, b) = a.lambda == b.lambda && a.u == b.u

    function arc_component_resolved(x::Float64, y::Float64, z::Float64)
        lo, hi = min(x,y,z), max(x,y,z)
        return lo == hi || nextfloat(lo) == hi
    end

    function arc_state_resolution(lo, mid, hi)
        return arc_component_resolved(lo.lambda,mid.lambda,hi.lambda) &&
            all(arc_component_resolved(lo.u[k],mid.u[k],hi.u[k]) for k in eachindex(lo.u))
    end

    # Events propose steps only; every proposal still needs a real assembly.
    function arc_event_fractions(path, start_gaps, end_gaps, mesh, param,
            history, model, max_fraction, linear_path,
            wall_started_ns, wall_limit; czm_workspace::CZMWorkspace)
        roots = Tuple{Float64,CartesianIndex{3}}[]
        samples = 0
        localized_bounds = nothing
        if linear_path && all(geom -> begin
                ip = collector_params(param, mesh.cohesive_elements[geom].interface_type)
                model == "model1" || (ip.δ_0 == ip.δ_0_t && ip.δ_c == ip.δ_c_t)
            end, eachindex(mesh.cohesive_elements))
            for i in axes(history,1), g in axes(history,2)
                arc_wall_expired(wall_started_ns, wall_limit) && return (;fractions=Float64[], key=nothing, status=:wall_time_limit, samples)
                ip = collector_params(param, mesh.cohesive_elements[i].interface_type)
                n0, t0 = start_gaps.separation_n[i,g], start_gaps.separation_t[i,g]
                nv = (end_gaps.separation_n[i,g] - n0) / max_fraction
                tv = (end_gaps.separation_t[i,g] - t0) / max_fraction
                D = history[i,g].D
                levels = (ip.δ_0 * ip.δ_c / (ip.δ_c - D * (ip.δ_c - ip.δ_0)), ip.δ_0, ip.δ_c)
                for (k, level) in enumerate(levels), opening in (false, true)
                    k == 1 && D == 1.0 && continue
                    aa = opening ? nv^2 : 0.0
                    bb = opening ? 2n0*nv : 0.0
                    cc = opening ? n0^2 : 0.0
                    if model == "mix"
                        aa += tv^2; bb += 2t0*tv; cc += t0^2
                    elseif !opening
                        continue
                    end
                    rr = arc_real_roots(aa, bb, cc-level^2)
                    rr === nothing && return (;fractions=Float64[], key=nothing, status=:nonfinite, samples)
                    for s in rr
                        0.0 < s < max_fraction || continue
                        (n0 + s*nv > 0.0) == opening || continue
                        push!(roots,(s,CartesianIndex(i,g,k)))
                    end
                end
                for (k, x, v) in ((4,n0,nv),(5,t0,tv))
                    k == 5 && model != "mix" && continue
                    v == 0.0 && continue
                    s = -x/v
                    0.0 < s < max_fraction && push!(roots,(s,CartesianIndex(i,g,k)))
                end
            end
        else
            previous_s = 0.0
            start_point = path(0.0)
            start_point === nothing && return (;fractions=Float64[], key=nothing, status=:invalid_path, samples)
            previous_point = start_point
            previous = arc_gp_event_values(start_gaps, mesh, param, history, model)
            for j in 1:16
                arc_wall_expired(wall_started_ns, wall_limit) && return (;fractions=Float64[], key=nothing, status=:wall_time_limit, samples)
                s = max_fraction * j / 16
                point = path(s)
                point === nothing && return (;fractions=Float64[], key=nothing, status=:invalid_path, samples)
                gaps = j == 16 ? end_gaps : cohesive_gp_separations(mesh, point.u, param; czm_workspace=czm_workspace)
                samples += j == 16 ? 0 : 1
                values = arc_gp_event_values(gaps, mesh, param, history, model)
                keys = [key for key in CartesianIndices(values)
                    if !(key[3] == 1 && history[key[1],key[2]].D == 1.0) &&
                       previous[key] != 0.0 &&
                       (values[key] == 0.0 || signbit(previous[key]) != signbit(values[key]))]
                if !isempty(keys)
                    lo, hi, values_lo = previous_s, s, previous
                    point_lo, point_hi = previous_point, point
                    while true
                        arc_wall_expired(wall_started_ns, wall_limit) && return (;fractions=Float64[], key=first(keys), status=:wall_time_limit, samples)
                        mid = lo + (hi-lo)/2
                        hi-lo <= max(64eps(max(abs(lo),abs(hi))), 1e-8*max(abs(lo),abs(hi))) && break
                        point = path(mid)
                        point === nothing && return (;fractions=Float64[], key=first(keys), status=:invalid_path, samples)
                        if arc_state_resolution(point_lo,point,point_hi)
                            if arc_same_point(point_lo,point_hi)
                                status = arc_same_point(point_lo,start_point) ?
                                    :no_state_progress : :inconsistent_event_state
                                return (;fractions=Float64[], key=first(keys), status, samples)
                            end
                            fractions = Float64[]
                            distinct_points = [start_point]
                            for (fraction, endpoint) in ((hi,point_hi),(lo,point_lo))
                                if !any(previous -> arc_same_point(endpoint,previous), distinct_points)
                                    push!(fractions,fraction)
                                    push!(distinct_points,endpoint)
                                end
                            end
                            status = isempty(fractions) ? :no_state_progress : :state_resolution
                            return (;fractions, key=first(keys), status, samples)
                        end
                        samples < 64 || return (;fractions=Float64[], key=first(keys), status=:localization_budget, samples)
                        values_mid = arc_gp_event_values(cohesive_gp_separations(mesh, point.u, param;
                            czm_workspace=czm_workspace), mesh, param, history, model)
                        samples += 1
                        left_keys = [key for key in keys if values_mid[key] == 0.0 ||
                            signbit(values_lo[key]) != signbit(values_mid[key])]
                        if isempty(left_keys)
                            lo, values_lo = mid, values_mid
                            point_lo = point
                        else
                            hi, keys = mid, left_keys
                            point_hi = point
                        end
                    end
                    hi-lo <= max(64eps(max(abs(lo),abs(hi))), 1e-8*max(abs(lo),abs(hi))) ||
                        return (;fractions=Float64[], key=first(keys), status=:localization_precision, samples)
                    push!(roots,((lo+hi)/2,first(keys)))
                    localized_bounds = (hi, lo)
                end
                isempty(roots) || break
                previous_s, previous = s, values
                previous_point = point
            end
        end
        isempty(roots) && return (;fractions=Float64[], key=nothing, status=:unbracketed, samples)
        sort!(roots; by=first)
        event, key = first(roots)
        # An analytic root is also a candidate, so Newton can reassemble
        # at the corner. Localized nonlinear paths retain their real bracket.
        fractions = localized_bounds === nothing ?
            unique([min(max_fraction, event*1.01), event, event*0.99]) :
            unique(collect(localized_bounds))
        return (;fractions, key, status=:located, samples)
    end

    function arc_real_roots(a::Float64, b::Float64, c::Float64)
        scale = max(abs(a), abs(b), abs(c))
        isfinite(scale) || return nothing
        scale == 0.0 && return ()
        a, b, c = a / scale, b / scale, c / scale
        a == 0.0 && return b == 0.0 ? () : (-c / b,)
        disc = b^2 - 4a*c
        roundoff = 8eps(Float64) * (b^2 + abs(4a*c))
        disc >= -roundoff || return ()
        root_disc = sqrt(max(disc, 0.0))
        root_disc == 0.0 && return (-b / (2a),)
        q = -0.5 * (b + copysign(root_disc, b))
        return (q / a, c / q)
    end

    # The split path exists only on the component containing s=0, eta=0.
    function cylindrical_arc_split_path(p::Vector{Float64}, q_R::Vector{Float64},
            q_F::Vector{Float64}, radius_sq::Float64)
        a = dot(q_F, q_F)
        b0 = 2dot(p, q_F)
        c0 = dot(p, p) - radius_sq
        isfinite(a) && a > 0.0 && isfinite(b0) && isfinite(c0) || return nothing
        abs(c0) <= 64eps(Float64) * radius_sq || return nothing
        abs(b0) > 64eps(Float64) * sqrt(a * radius_sq) || return nothing
        b0^2 - 4a*c0 >= 0.0 || return nothing
        eta0 = -2c0 / (b0 + copysign(sqrt(b0^2 - 4a*c0), b0))
        isfinite(eta0) && abs(eta0) * sqrt(a) <= 64eps(Float64) * sqrt(radius_sq) || return nothing
        b1 = 2dot(q_R, q_F)
        c1, c2 = 2dot(p, q_R), dot(q_R, q_R)
        d0, d1, d2 = b0^2 - 4a*c0, 2b0*b1 - 4a*c1, b1^2 - 4a*c2
        roots = arc_real_roots(d2, d1, d0)
        roots === nothing && return nothing
        max_fraction = 1.0
        for root in roots
            isfinite(root) || return nothing
            if 0.0 < root <= max_fraction
                max_fraction = prevfloat(root)
            end
        end
        max_fraction > 0.0 || return nothing
        return (; a, b0, b1, c0, c1, c2, sigma=sign(b0), max_fraction)
    end

    function cylindrical_arc_correction(delta_u::Vector{Float64},
            delta_u_R::Vector{Float64}, delta_u_F::Vector{Float64},
            radius_sq::Float64, direction_u::Vector{Float64};
            residual_fraction::Float64=1.0,
            path=cylindrical_arc_split_path(delta_u, delta_u_R, delta_u_F, radius_sq))
        path === nothing && return nothing
        0.0 < residual_fraction <= path.max_fraction || return nothing
        s = residual_fraction
        a, b, c = path.a, path.b0 + s * path.b1,
            path.c0 + s * path.c1 + s^2 * path.c2
        coefficient_scale = max(abs(a), abs(b), abs(c))
        isfinite(coefficient_scale) && coefficient_scale > 0.0 || return nothing
        a, b, c = a / coefficient_scale, b / coefficient_scale, c / coefficient_scale
        a > 0.0 || return nothing
        discriminant = b^2 - 4.0 * a * c
        # Only a negative discriminant within the arithmetic rounding bound
        # represents a repeated real root; a genuinely negative one fails.
        discriminant_roundoff = 8.0 * eps(Float64) * (b^2 + abs(4.0 * a * c))
        discriminant >= -discriminant_roundoff || return nothing
        discriminant = max(discriminant, 0.0)
        sqrt_disc = sqrt(discriminant)
        delta_lambda_corr = if sign(b) == path.sigma
            -2c / (b + path.sigma * sqrt_disc)
        else
            (-b + path.sigma * sqrt_disc) / (2a)
        end
        delta_u_corr = s * delta_u_R + delta_lambda_corr * delta_u_F
        all(isfinite, delta_u_corr) && isfinite(delta_lambda_corr) || return nothing
        cylindrical_arc_forward(delta_u + delta_u_corr, direction_u) || return nothing
        return (delta_u=delta_u_corr, delta_lambda=delta_lambda_corr,
            residual_fraction=s, root_sigma=path.sigma)
    end

    function cylindrical_arc_bordered_solve(K_bc::SparseMatrixCSC{Float64,Int64},
            F_load_bc::Vector{Float64}, delta_u::Vector{Float64},
            delta_lambda::Float64, arc_length_alpha::Float64,
            rhs_u::Vector{Float64}, rhs_constraint::Float64; bc_dofs=nothing)
        A = build_arc_length_augmented_matrix(K_bc, F_load_bc, delta_u,
            delta_lambda, arc_length_alpha)
        # 方案1（2026-09-26）约束行平衡：约束行元素（2Δu ~1e-6）与 K 行（罚刚度 ~1e10+）
        # 量级差可达 10 个数量级，∞-范数后向误差被 K 行主导而误拒可用解（Chan–Saad
        # bordered 求解的标准处理）。行缩放不改变线性系统的精确解，只改善条件与误差度量。
        constraint_row_scale = 2.0 * max(norm(delta_u),
            arc_length_alpha^2 * abs(delta_lambda))
        if isfinite(constraint_row_scale) && constraint_row_scale > 0.0
            A[size(A, 1), :] .*= inv(constraint_row_scale)
        end
        rhs = vcat(rhs_u, rhs_constraint *
            (constraint_row_scale > 0.0 ? inv(constraint_row_scale) : 1.0))
        if bc_dofs isa MechanicalBCDOFs && bc_dofs.gauge !== nothing
            z = bc_dofs.gauge
            abs(dot(z,rhs_u)) <= 1e-8*max(1.0,norm(rhs_u)) || return nothing
            column = sparse(reshape(vcat(z,0.0),:,1))
            A = [A column; transpose(column) spzeros(1,1)]
            rhs = vcat(rhs,0.0)
        end
        solution = try
            A \ rhs
        catch err
            err isa SingularException || rethrow()
            return nothing
        end
        all(isfinite, solution) || return nothing
        backward_scale = opnorm(A, Inf) * norm(solution, Inf) + norm(rhs, Inf)
        backward_scale > 0.0 || return nothing
        backward_error = norm(A * solution - rhs, Inf) / backward_scale
        isfinite(backward_error) && backward_error < 1e-8 || return nothing
        return (delta_u=solution[1:length(rhs_u)], delta_lambda=solution[length(rhs_u)+1],
            backward_error=backward_error)
    end

    function cylindrical_arc_bordered_predictor(K_bc::SparseMatrixCSC{Float64,Int64},
            F_load_bc::Vector{Float64}, previous_u::Vector{Float64},
            radius::Float64; bc_dofs=nothing)
        radius > 0.0 || return nothing
        tangent = cylindrical_arc_bordered_solve(K_bc, F_load_bc,
            previous_u, 0.0, 0.0,
            zeros(length(F_load_bc)), 2.0; bc_dofs=bc_dofs)
        tangent === nothing && return nothing
        tangent_norm = norm(tangent.delta_u)
        tangent_norm > 0.0 || return nothing
        scale = radius / tangent_norm
        delta_u = scale * tangent.delta_u
        cylindrical_arc_forward(delta_u, previous_u) || return nothing
        return (delta_u=delta_u, delta_lambda=scale * tangent.delta_lambda,
            radius_sq=dot(delta_u, delta_u),
            backward_error=tangent.backward_error)
    end

    function cylindrical_arc_bordered_correction(K_bc::SparseMatrixCSC{Float64,Int64},
            R_bc::Vector{Float64}, F_load_bc::Vector{Float64},
            delta_u::Vector{Float64}, radius_sq::Float64; bc_dofs=nothing)
        g = dot(delta_u, delta_u) - radius_sq
        return cylindrical_arc_bordered_solve(K_bc, F_load_bc,
            delta_u, 0.0, 0.0, R_bc, -g; bc_dofs=bc_dofs)
    end

    # A singular initial tangent has no previously accepted branch direction.
    # Later steps may use that real direction to make the augmented solve unique.
    function spherical_arc_tangent(K_bc::SparseMatrixCSC{Float64,Int64},
            f_hat::Vector{Float64}, previous_u::Union{Nothing,Vector{Float64}},
            previous_lambda::Float64, alpha::Float64; bc_dofs=nothing)
        delta_u = try
            K_bc \ f_hat
        catch err
            err isa SingularException || rethrow()
            nothing
        end
        delta_lambda = 1.0
        if delta_u === nothing
            previous_u === nothing && return nothing
            bordered = cylindrical_arc_bordered_solve(K_bc, f_hat, previous_u,
                previous_lambda, alpha, zeros(length(f_hat)), 2.0; bc_dofs=bc_dofs)
            bordered === nothing && return nothing
            delta_u, delta_lambda = bordered.delta_u, bordered.delta_lambda
        end
        all(isfinite,delta_u) && isfinite(delta_lambda) || error(
            "spherical arc tangent produced a non-finite direction")
        augmented_norm = hypot(norm(delta_u), alpha * delta_lambda)
        isfinite(augmented_norm) && augmented_norm > 0.0 || error(
            "spherical arc tangent has invalid norm $augmented_norm")
        if previous_u !== nothing
            previous_norm = hypot(norm(previous_u), alpha * previous_lambda)
            previous_norm > 0.0 || return nothing
            cosine = dot(delta_u ./ augmented_norm, previous_u ./ previous_norm) +
                (alpha * delta_lambda / augmented_norm) * (alpha * previous_lambda / previous_norm)
            abs(cosine) > cylindrical_arc_direction_tolerance(length(delta_u)+1) || return nothing
            if cosine < 0.0
                delta_u, delta_lambda = -delta_u, -delta_lambda
            end
        end
        return (;delta_u,delta_lambda,augmented_norm)
    end

    function spherical_arc_corrector(K_bc::SparseMatrixCSC{Float64,Int64},
            R_bc::Vector{Float64}, f_hat::Vector{Float64},
            delta_u::Vector{Float64}, delta_lambda::Float64,
            alpha::Float64, radius::Float64; bc_dofs=nothing)
        split = try
            directions = K_bc \ hcat(R_bc, f_hat)
            (directions[:, 1], directions[:, 2])
        catch err
            err isa SingularException || rethrow()
            nothing
        end
        if split !== nothing
            delta_u_R,delta_u_F = split
            all(isfinite,delta_u_R) && all(isfinite,delta_u_F) || error(
                "spherical arc split solve produced a non-finite correction")
            denominator = 2.0 * dot(delta_u,delta_u_F) + 2.0 * alpha^2 * delta_lambda
            denominator_scale = 2.0 * (norm(delta_u) * norm(delta_u_F) + alpha^2 * abs(delta_lambda))
            if abs(denominator) > eps(Float64) * max(1.0,denominator_scale)
                return spherical_arc_length_correction(delta_u,delta_lambda,
                    delta_u_R,delta_u_F,alpha,radius)
            end
        end
        g = dot(delta_u,delta_u) + alpha^2 * delta_lambda^2 - radius^2
        bordered = cylindrical_arc_bordered_solve(K_bc,f_hat,delta_u,delta_lambda,
            alpha,R_bc,-g;bc_dofs=bc_dofs)
        bordered === nothing && return nothing
        return bordered.delta_u,bordered.delta_lambda,g
    end

    # Each Newton candidate derives its damage from the last accepted substep.
    # A rejected candidate must not become the history for the next candidate.
    function arc_equilibrium_trial(czm_mesh::CohesiveMesh, u::Vector{Float64},
            param::Params, gp_damage_start::Matrix{DamageState},
            F_applied::Vector{Float64}, bc_dofs::AbstractVector{Int64},
            bc_vals::Vector{Float64}; czm_workspace::CZMWorkspace, visc_beta::Float64=1.0,
            czm_model::String="model1", gp_history_loading=nothing)
        weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
        gp_trial_states = similar(gp_damage_start)
        capture = (separation_n=zeros(size(gp_damage_start)),
            separation_t=zeros(size(gp_damage_start)),
            traction_n=zeros(size(gp_damage_start)), traction_t=zeros(size(gp_damage_start)))
        K_total, f_int, separations, tractions = assemble_coupled_system(
            czm_mesh, u, param; czm_workspace=czm_workspace,
            gp_damage_states=gp_damage_start,
            gp_trial_states=gp_trial_states,
            gp_capture=capture,
            gp_history_loading=gp_history_loading,
            visc_beta=visc_beta, czm_model=czm_model)
        trial_damage = aggregate_gp_damage_states(gp_trial_states, weights)
        R = F_applied - f_int
        set_czm_residual!(R, u, bc_dofs, bc_vals)
        return (damage_states=trial_damage,
            gp_damage_states=gp_trial_states, K_total=K_total, R=R,
            residual_norm=norm(R), separations=separations,
            tractions=tractions, gp_separations=capture)
    end

    # Fixed-load endpoint of the nongeometric arc path. Every line-search
    # candidate uses the same accepted history, including after a rejected
    # candidate; only a balanced candidate advances MechState.
    function solve_czm_arc_target_step(czm_mesh::CohesiveMesh,
            load::MechanicalLoad, param::Params, ms::MechState; czm_workspace::CZMWorkspace,
            max_iter::Int=50, tol::Float64=1e-8, visc_beta::Float64=1.0,
            czm_model::String="model1", fix_inner::Bool=true,
            initial_u::Union{Nothing,Vector{Float64}}=nothing,
            wall_started_ns::UInt64=time_ns(), wall_time_limit_seconds::Float64=0.0, boundary=nothing,
            costs=Dict(:assemblies=>0,:candidates=>0,:event_samples=>0,:linear_solves=>0,:corrections=>0))
        F_ext = load.F_ext
        u_start = copy(ms.u_prev)
        u = initial_u === nothing ? copy(u_start) : copy(initial_u)
        gp_damage_start = gp_history_for_step(ms, czm_mesh)
        weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
        damage_start = aggregate_gp_damage_states(gp_damage_start, weights)
        result = CZMResult(length(u), length(czm_mesh.cohesive_elements))
        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner, boundary=boundary)
        boundary === nothing || apply_czm_dirichlet!(u,bc_dofs,bc_vals)
        F_target = F_ext + assemble_thermal_chemical_load(czm_mesh, param,
            load.dT, load.Δsoc_n, load.Δsoc_p)
        check_czm_external_load(bc_dofs,F_target)

        last_residual = Inf
        corrections = 0

        for iter in 0:max_iter
            arc_wall_expired(wall_started_ns, wall_time_limit_seconds) && break
            trial = arc_equilibrium_trial(czm_mesh, u, param, gp_damage_start,
                F_target, bc_dofs, bc_vals; czm_workspace=czm_workspace,
                visc_beta=visc_beta,
                czm_model=czm_model)
            costs[:assemblies] += 1
            last_residual = trial.residual_norm
            if last_residual < tol
                committed = arc_equilibrium_trial(czm_mesh, u, param, trial.gp_damage_states,
                    F_target, bc_dofs, bc_vals; czm_workspace=czm_workspace,
                    visc_beta=0.0, czm_model=czm_model)
                costs[:assemblies] += 1
                committed.residual_norm < tol || break
                result.converged = true
                result.iterations = corrections
                result.residual_norm = committed.residual_norm
                fill_czm_result!(result, u, trial.damage_states,
                    committed.separations, committed.tractions)
                arc_observe(() -> (;phase=:arc_stop,reason=:balanced_target,costs=copy(costs),corrections=costs[:corrections]))
                commit_czm_state!(ms,czm_mesh,u,trial.damage_states,trial.gp_damage_states,
                    mechanical_eigenstrain(load))
                return result
            end
            iter == max_iter && break

            K_bc, R_bc = apply_bc_czm(trial.K_total, trial.R;
                bc_dofs=bc_dofs, bc_vals=bc_vals)
            Δu = try
                costs[:linear_solves] += 1
                czm_linear_solve(K_bc, R_bc, bc_dofs; czm_workspace=czm_workspace)
            catch err
                err isa Union{SingularException,ZeroPivotException} || rethrow()
                break
            end
            all(isfinite, Δu) || break
            corrections += 1
            costs[:corrections] += 1

            accepted = false
            fractions = [2.0^-j for j in 0:10]
            path = alpha -> begin
                alpha == 0.0 && return (u=copy(u),lambda=1.0)
                candidate = u + alpha * Δu
                apply_czm_dirichlet!(candidate, bc_dofs, bc_vals)
                (u=candidate, lambda=1.0)
            end
            for evaluation in 1:11
                arc_wall_expired(wall_started_ns, wall_time_limit_seconds) && break
                isempty(fractions) && break
                α = popfirst!(fractions)
                candidate = path(α).u
                candidate == u && break
                candidate_trial = arc_equilibrium_trial(czm_mesh, candidate,
                    param, gp_damage_start, F_target, bc_dofs, bc_vals; czm_workspace=czm_workspace,
                    visc_beta=visc_beta, czm_model=czm_model)
                costs[:assemblies] += 1
                costs[:candidates] += 1
                phi = 0.5 * (last_residual/tol)^2
                candidate_phi = 0.5 * (candidate_trial.residual_norm/tol)^2
                candidate_ok = isfinite(candidate_phi) && candidate_phi < phi &&
                    candidate_phi <= (1.0 - 2e-4*α)*phi
                arc_observe(() -> (;phase=:arc_candidate,attempt=0,accepted=candidate_ok,
                    lambda=1.0,residual_norm=candidate_trial.residual_norm,constraint_relative=0.0,
                    gp_history_start=clone_damage_states(gp_damage_start),mode=:target,fraction=α))
                if candidate_ok
                    u = candidate
                    accepted = true
                    break
                end
                if evaluation == 1
                    events = arc_event_fractions(path, trial.gp_separations,
                        candidate_trial.gp_separations, czm_mesh, param, gp_damage_start,
                        czm_model, 1.0, true, wall_started_ns, wall_time_limit_seconds; czm_workspace=czm_workspace)
                    costs[:event_samples] += events.samples
                    @debug "arc target event" status=events.status key=events.key samples=events.samples
                    events.status in (:wall_time_limit,:localization_budget,:localization_precision,:nonfinite,:invalid_path,
                        :no_state_progress,:inconsistent_event_state) && break
                    prepend!(fractions, events.fractions)
                end
            end
            accepted || break
        end

        _, f_int_start, separations, tractions = assemble_coupled_system(czm_mesh,
            u_start, param; czm_workspace=czm_workspace,
            gp_damage_states=gp_damage_start,
            visc_beta=0.0, czm_model=czm_model)
        costs[:assemblies] += 1
        R_start = F_target - f_int_start
        set_czm_residual!(R_start, u_start, bc_dofs, bc_vals)
        result.converged = false
        result.iterations = corrections
        result.residual_norm = norm(R_start)
        fill_czm_result!(result, u_start, damage_start, separations, tractions)
        arc_observe(() -> (;phase=:arc_stop,reason=:target_not_converged,costs=copy(costs),corrections=costs[:corrections]))
        return result
    end

    # 方案1（2026-09-26）实验口径：JUBAT_ARC_CORRECTOR=bordered 时增广系统作主校正器。
    # 环境变量门控以保持缺省行为与公开 API 不变；实验通过后应提升为正式选项。
    arc_bordered_primary_corrector() =
        lowercase(get(ENV, "JUBAT_ARC_CORRECTOR", "split")) == "bordered"

    function solve_czm_arc_length_step(czm_mesh::CohesiveMesh, load::MechanicalLoad, param, ms::MechState; czm_workspace::CZMWorkspace, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, arc_length_alpha::Float64=1.0, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, wall_time_limit_seconds::Float64=0.0, boundary=nothing)
        F_ext = load.F_ext
        isfinite(visc_beta) && visc_beta == 1.0 || throw(ArgumentError(
            "nongeometric arc_length requires rate-independent visc_beta=1"))
        isfinite(wall_time_limit_seconds) && wall_time_limit_seconds >= 0.0 ||
            throw(ArgumentError("arc wall-time limit must be finite and nonnegative"))
        wall_started_ns = time_ns()
        nnode = size(czm_mesh.node, 1)
        ndof = 2 * nnode
        n_coh = length(czm_mesh.cohesive_elements)

        result = CZMResult(ndof, n_coh)
        u = copy(ms.u_prev)
        gp_states = gp_history_for_step(ms, czm_mesh)
        any(s -> s.D_visc != s.D, gp_states) && throw(ArgumentError(
            "nongeometric arc_length cannot advance lagged GP viscous damage without physical-time relaxation"))
        weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
        damage_states = aggregate_gp_damage_states(gp_states, weights)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner, boundary=boundary)
        boundary === nothing || apply_czm_dirichlet!(u,bc_dofs,bc_vals)

        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param, load.dT, load.Δsoc_n, load.Δsoc_p)

        costs = Dict(:assemblies=>1,:candidates=>0,:event_samples=>0,:linear_solves=>0,:corrections=>0)
        # 增量载荷参考
        _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param; czm_workspace=czm_workspace,
            gp_damage_states=gp_states,
            visc_beta=0.0, czm_model=czm_model)
        F_target = F_ext + F_thermo_chem_total
        check_czm_external_load(bc_dofs,F_target)
        F_delta = F_target - f_int_ref
        F_load_bc = copy(F_delta)
        zero_czm_bc_entries!(F_load_bc, bc_dofs)
        if all(iszero, F_load_bc)
            return solve_czm_arc_target_step(czm_mesh, load, param, ms; czm_workspace=czm_workspace,
                max_iter=max_iter, tol=tol,
                visc_beta=visc_beta, czm_model=czm_model,
                fix_inner=fix_inner, boundary=boundary, wall_started_ns=wall_started_ns,
                wall_time_limit_seconds=wall_time_limit_seconds, costs=costs)
        end

        total_iter = 0
        load_progress = 0.0
        load_step = 0
        max_substep_attempts = max(100, 20 * n_load_steps)
        # A failed physical step must not consume an unbounded sequence of
        # expensive sparse solves while repeatedly shrinking arc substeps.
        max_total_iter = max(2 * max_iter, 8 * n_load_steps)
        initial_increment = 1.0 / max(1, n_load_steps)
        step_size = initial_increment
        step_size_min = initial_increment / 128.0
        target_neighborhood = initial_increment
        last_residual = Inf
        converged_substep = false
        previous_u_tangent = nothing
        previous_arc_control = nothing
        last_target_attempt_start = NaN
        stop_reason = :none
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)

        function attempt_target_step(u0, damage0, gp0, target_max_iter; initial_u=nothing)
            target_ms = MechState(czm_mesh)
            target_ms.u_prev = copy(u0)
            target_ms.damage_states = clone_damage_states(damage0)
            target_ms.gp_damage_states = clone_damage_states(gp0)
            target_result = solve_czm_arc_target_step(czm_mesh, load, param, target_ms;
                czm_workspace=czm_workspace,
                max_iter=target_max_iter, tol=tol, visc_beta=visc_beta,
                czm_model=czm_model, fix_inner=fix_inner, boundary=boundary,
                initial_u=initial_u, wall_started_ns=wall_started_ns,
                wall_time_limit_seconds=wall_time_limit_seconds, costs=costs)
            return target_result, target_ms
        end

        while load_progress < 1.0 - 1e-12
            if arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
                stop_reason = :wall_time_limit
                break
            end
            if total_iter >= max_total_iter
                stop_reason = :correction_budget
                break
            end
            load_step += 1
            if load_step > max_substep_attempts
                stop_reason = :substep_budget
                break
            end
            load_start = load_progress
            target_progress = min(1.0, load_start + step_size)
            # Keep the continuation scale independent of the remaining load;
            # an equilibrated overshoot is handled by the strict target solve.
            delta_nom_current = step_size

            u_start = copy(u)
            damage_start = clone_damage_states(damage_states)   # 子步试探回滚（ms 不受影响）
            gp_start = clone_damage_states(gp_states)
            converged_substep = false
            overshoot_u = nothing
            last_residual = Inf
            corrector_mode = arc_bordered_primary_corrector() ? :bordered : :split
            predictor_iterations = 0

            F_applied = f_int_ref + load_start * F_delta
            start_trial = arc_equilibrium_trial(czm_mesh, u, param, gp_start,
                F_applied, bc_dofs, bc_vals; czm_workspace=czm_workspace,
                visc_beta=0.0, czm_model=czm_model)
            K_total, R = start_trial.K_total, start_trial.R
            separations, tractions = start_trial.separations, start_trial.tractions
            costs[:assemblies] += 1

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            branch = arc_gp_history_boundary(start_trial.gp_separations,
                czm_mesh, param, gp_start, czm_model)
            # Without incoming damage growth, resume the prescribed load direction.
            incoming_loading = false
            if previous_u_tangent !== nothing && any(branch.boundary)
                incoming_gaps = cohesive_gp_separations(czm_mesh, previous_u_tangent, param;
                    czm_workspace=czm_workspace)
                incoming_gamma = branch.q_n .* incoming_gaps.separation_n +
                    branch.q_t .* incoming_gaps.separation_t
                incoming_loading = any(branch.boundary .& (incoming_gamma .> 0.0))
            end
            tangent = try
                costs[:linear_solves] += 1
                czm_linear_solve(K_bc, F_load_bc, bc_dofs; czm_workspace=czm_workspace)
            catch err
                err isa SingularException || rethrow()
                nothing
            end

            predictor = nothing
            if tangent !== nothing && all(isfinite, tangent)
                predictor = cylindrical_arc_predictor(tangent, load_start,
                    delta_nom_current, incoming_loading ? previous_u_tangent : nothing)
            elseif incoming_loading
                # At a fold K can be singular although the bordered arc
                # system remains invertible. Use the accepted direction.
                radius = previous_arc_control.radius *
                    (delta_nom_current / previous_arc_control.nominal_increment)
                bordered = cylindrical_arc_bordered_predictor(K_bc,
                    F_load_bc, previous_u_tangent, radius; bc_dofs=bc_dofs)
                costs[:linear_solves] += 1
                if bordered !== nothing
                    predictor = (delta_u=bordered.delta_u,
                        delta_lambda=bordered.delta_lambda,
                        lambda=load_start + bordered.delta_lambda,
                        radius_sq=bordered.radius_sq)
                end
            else
                stop_reason = :invalid_target_predictor_solve
                break
            end
            if predictor !== nothing
                if previous_arc_control !== nothing
                    radius = previous_arc_control.radius *
                        (delta_nom_current / previous_arc_control.nominal_increment)
                    factor = radius / sqrt(predictor.radius_sq)
                    du = factor * predictor.delta_u
                    dl = factor * predictor.delta_lambda
                    predictor = (delta_u=du,delta_lambda=dl,
                        lambda=load_start+dl,radius_sq=dot(du,du))
                end
                if any(branch.boundary)
                    seed_direction = incoming_loading ?
                        previous_u_tangent : predictor.delta_u
                    secant_gaps = cohesive_gp_separations(czm_mesh, seed_direction, param;
                        czm_workspace=czm_workspace)
                    gamma = branch.q_n .* secant_gaps.separation_n +
                        branch.q_t .* secant_gaps.separation_t
                    loading = branch.boundary .& (gamma .> 0.0)
                    seen = BitMatrix[]
                    pending = BitMatrix[copy(loading)]
                    radius = sqrt(predictor.radius_sq)
                    predictor = nothing
                    for branch_iter in 1:max_iter
                        if arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
                            stop_reason = :wall_time_limit
                            break
                        end
                        if total_iter >= max_total_iter
                            stop_reason = :correction_budget
                            break
                        end
                        if isempty(pending)
                            stop_reason = :history_branch_exhausted
                            break
                        end
                        loading = popfirst!(pending)
                        push!(seen,copy(loading))
                        total_iter += 1
                        predictor_iterations += 1
                        costs[:corrections] += 1
                        branch_trial = arc_equilibrium_trial(czm_mesh, u_start,
                            param, gp_start, F_applied, bc_dofs, bc_vals;
                            czm_workspace=czm_workspace,
                            visc_beta=1.0,
                            czm_model=czm_model, gp_history_loading=loading)
                        costs[:assemblies] += 1
                        branch_K, _ = apply_bc_czm(branch_trial.K_total, branch_trial.R;
                            bc_dofs=bc_dofs, bc_vals=bc_vals)
                        costs[:linear_solves] += 1
                        branch_tangent = cylindrical_arc_bordered_solve(branch_K,
                            F_load_bc, seed_direction, 0.0, 0.0,
                            zeros(ndof), 2.0; bc_dofs=bc_dofs)
                        if branch_tangent === nothing
                            stop_reason = :invalid_history_tangent
                            break
                        end
                        tangent_gaps = cohesive_gp_separations(czm_mesh, branch_tangent.delta_u, param;
                            czm_workspace=czm_workspace)
                        gamma = branch.q_n .* tangent_gaps.separation_n +
                            branch.q_t .* tangent_gaps.separation_t
                        orientation = 0.0
                        for orientation_sign in (1.0,-1.0)
                            consistent = all(!branch.boundary[k] ||
                                (loading[k] ? orientation_sign*gamma[k] >= 0.0 : orientation_sign*gamma[k] <= 0.0)
                                for k in eachindex(loading))
                            corner_forward = any(loading[k] && orientation_sign*gamma[k] > 0.0
                                for k in eachindex(loading))
                            if consistent && (orientation_sign > 0.0 || corner_forward)
                                orientation = orientation_sign
                                break
                            end
                        end
                        if orientation != 0.0
                            factor = orientation * radius / norm(branch_tangent.delta_u)
                            du = factor * branch_tangent.delta_u
                            dl = factor * branch_tangent.delta_lambda
                            predictor = (delta_u=du,delta_lambda=dl,
                                lambda=load_start+dl,radius_sq=dot(du,du))
                            @debug "arc history predictor" load_start loading_count=count(loading) orientation
                            break
                        end
                        for next_loading in (branch.boundary .& (gamma .> 0.0),
                                branch.boundary .& (gamma .< 0.0))
                            if !any(previous -> previous == next_loading, seen) &&
                                    !any(queued -> queued == next_loading, pending)
                                push!(pending,next_loading)
                            end
                        end
                    end
                    if predictor === nothing
                        stop_reason == :none && (stop_reason = :history_branch_budget)
                        break
                    end
                end
            end
            if predictor === nothing
                step_size *= 0.5
                if step_size < step_size_min
                    stop_reason = :invalid_predictor_solve
                    break
                end
                continue
            end
            delta_u_pred = predictor.delta_u
            arc_target_sq = predictor.radius_sq
            if !isfinite(arc_target_sq) || arc_target_sq <= 0.0
                stop_reason = :invalid_arc_radius
                break
            end

            u = u_start + delta_u_pred
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
            load_progress = predictor.lambda

            for iter in 1:(max_iter - predictor_iterations)
                if arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
                    stop_reason = :wall_time_limit
                    break
                end
                if total_iter >= max_total_iter
                    stop_reason = :correction_budget
                    break
                end
                total_iter += 1
                costs[:corrections] += 1

                F_applied = f_int_ref + load_progress * F_delta
                trial = arc_equilibrium_trial(czm_mesh, u, param,
                    gp_start, F_applied, bc_dofs, bc_vals; czm_workspace=czm_workspace,
                    visc_beta=visc_beta, czm_model=czm_model)
                costs[:assemblies] += 1
                K_total, R = trial.K_total, trial.R
                separations, tractions = trial.separations, trial.tractions

                delta_u = u - u_start
                delta_u_sq = dot(delta_u, delta_u)
                arc_constraint = delta_u_sq - arc_target_sq
                residual_norm = sqrt(trial.residual_norm^2 + arc_constraint^2)
                last_residual = residual_norm

                endpoint_action = cylindrical_arc_endpoint_action(
                    trial.residual_norm, delta_u_sq, arc_target_sq,
                    tol, load_progress)
                direction_u = delta_u_pred
                direction_ok = cylindrical_arc_forward(delta_u, direction_u)
                direction_ok || (endpoint_action = :continue)
                if endpoint_action != :continue
                    last_residual = trial.residual_norm
                    if endpoint_action == :accept
                        damage_states = trial.damage_states
                        gp_states = trial.gp_damage_states
                        converged_substep = true
                        arc_observe(() -> (;phase=:arc_substep,attempt=load_step,lambda=load_progress,
                            u=copy(u),gp_history_start=clone_damage_states(gp_start),
                            gp_damage_states=clone_damage_states(gp_states),residual_norm=trial.residual_norm,
                            constraint_relative=abs(arc_constraint)/arc_target_sq,direction_ok))
                        @debug "arc accepted" load_start load_progress u_norm=norm(u) damage_max=maximum(s.D for s in gp_states)
                        previous_u_tangent = u - u_start
                        previous_arc_control = (nominal_increment=delta_nom_current,
                            radius=sqrt(delta_u_sq))
                        step_size *= 1.25
                        break
                    elseif endpoint_action == :target
                        # The arc point is already equilibrated beyond the
                        # prescribed load. Repeating Newton cannot bring it
                        # back to λ=1; use only its displacement as the target
                        # guess, with history still anchored at u_start.
                        overshoot_u = copy(u)
                        break
                    end
                end

                K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

                # A bordered correction may leave the circle. This attempt
                # therefore never returns to split after entering that mode.
                direction_u = delta_u_pred
                phi = arc_merit(trial.residual_norm, arc_constraint, arc_target_sq, tol)
                candidate_accepted = false
                candidate_search_failed = false
                evaluations = 0
                while evaluations < 11
                    arc_wall_expired(wall_started_ns, wall_time_limit_seconds) && break
                    correction = nothing
                    split_path = nothing
                    if corrector_mode == :split
                        split = try
                            costs[:linear_solves] += 1
                            q_R = czm_linear_solve(K_bc, R_bc, bc_dofs; czm_workspace=czm_workspace)
                            costs[:linear_solves] += 1
                            (q_R, czm_linear_solve(K_bc, F_load_bc, bc_dofs; czm_workspace=czm_workspace))
                        catch err
                            err isa Union{SingularException,ZeroPivotException} || rethrow()
                            nothing
                        end
                        if split !== nothing && all(isfinite, split[1]) && all(isfinite, split[2])
                            split_path = cylindrical_arc_split_path(delta_u, split[1], split[2], arc_target_sq)
                        end
                        if split_path === nothing
                            corrector_mode = :bordered
                            continue
                        end
                        max_fraction = split_path.max_fraction
                        path = fraction -> begin
                            fraction == 0.0 && return (u=copy(u),lambda=load_progress)
                            corr = cylindrical_arc_correction(delta_u, split[1], split[2],
                                arc_target_sq, direction_u; residual_fraction=fraction, path=split_path)
                            corr === nothing && return nothing
                            candidate = u + corr.delta_u
                            apply_czm_dirichlet!(candidate, bc_dofs, bc_vals)
                            (u=candidate, lambda=load_progress+corr.delta_lambda)
                        end
                    else
                        costs[:linear_solves] += 1
                        correction = cylindrical_arc_bordered_correction(K_bc, R_bc,
                            F_load_bc, delta_u, arc_target_sq; bc_dofs=bc_dofs)
                        correction === nothing && break
                        max_fraction = 1.0
                        path = fraction -> begin
                            fraction == 0.0 && return (u=copy(u),lambda=load_progress)
                            candidate = u + fraction * correction.delta_u
                            apply_czm_dirichlet!(candidate, bc_dofs, bc_vals)
                            (u=candidate, lambda=load_progress+fraction*correction.delta_lambda)
                        end
                    end
                    fractions = [max_fraction * 2.0^-j for j in 0:10]
                    first_evaluation = true
                    while evaluations < 11 && !isempty(fractions)
                        if arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
                            candidate_search_failed = true
                            stop_reason = :wall_time_limit
                            break
                        end
                        fraction = popfirst!(fractions)
                        candidate = path(fraction)
                        if candidate === nothing
                            break
                        end
                        candidate.u == u && candidate.lambda == load_progress && break
                        candidate_trial = arc_equilibrium_trial(czm_mesh, candidate.u, param,
                            gp_start, f_int_ref + candidate.lambda * F_delta, bc_dofs, bc_vals;
                            czm_workspace=czm_workspace,
                            visc_beta=visc_beta, czm_model=czm_model)
                        costs[:assemblies] += 1
                        costs[:candidates] += 1
                        evaluations += 1
                        candidate_p = candidate.u - u_start
                        candidate_g = dot(candidate_p,candidate_p) - arc_target_sq
                        candidate_phi = arc_merit(candidate_trial.residual_norm, candidate_g, arc_target_sq, tol)
                        direction_ok = cylindrical_arc_forward(candidate_p, direction_u)
                        candidate_ok = direction_ok && isfinite(candidate_phi) && candidate_phi < phi &&
                            candidate_phi <= (1.0 - 2e-4*fraction)*phi
                        arc_observe(() -> (;phase=:arc_candidate,attempt=load_step,accepted=candidate_ok,
                            lambda=candidate.lambda,residual_norm=candidate_trial.residual_norm,
                            constraint_relative=abs(candidate_g)/arc_target_sq,
                            gp_history_start=clone_damage_states(gp_start),mode=corrector_mode,fraction))
                        if candidate_ok
                            u, load_progress = candidate.u, candidate.lambda
                            candidate_accepted = true
                            break
                        end
                        if first_evaluation
                            events = arc_event_fractions(path, trial.gp_separations,
                                candidate_trial.gp_separations, czm_mesh, param, gp_start,
                                czm_model, max_fraction, corrector_mode == :bordered,
                                wall_started_ns, wall_time_limit_seconds; czm_workspace=czm_workspace)
                            costs[:event_samples] += events.samples
                            @debug "arc correction event" attempt=load_step mode=corrector_mode status=events.status key=events.key samples=events.samples
                            if events.status in (:wall_time_limit,:localization_budget,:localization_precision,:nonfinite,:invalid_path,
                                    :no_state_progress,:inconsistent_event_state)
                                candidate_search_failed = true
                                stop_reason = events.status
                                break
                            end
                            prepend!(fractions, events.fractions)
                            first_evaluation = false
                        end
                    end
                    candidate_accepted && break
                    candidate_search_failed && break
                    corrector_mode == :bordered && break
                    corrector_mode = :bordered
                end
                candidate_accepted || break
            end

            if !converged_substep
                u = u_start
                damage_states = damage_start
                gp_states = gp_start
                load_progress = load_start
                if arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
                    stop_reason = :wall_time_limit
                    break
                end
                if (target_progress == 1.0 || overshoot_u !== nothing) &&
                   (load_start >= 1.0 - target_neighborhood || overshoot_u !== nothing) &&
                   load_start != last_target_attempt_start && total_iter < max_total_iter
                    # The arc constraint could not reach the prescribed load
                    # from this accepted state. Try its strict fixed-load
                    # equilibrium once before retrying smaller arc radii.
                    last_target_attempt_start = load_start
                    corrected, trial_ms = attempt_target_step(u_start, damage_start, gp_start,
                        min(max_iter, max_total_iter - total_iter); initial_u=overshoot_u)
                    total_iter += corrected.iterations
                    if corrected.converged
                        @debug "arc target correction accepted" load_start target_progress u_norm=norm(corrected.displacement) damage_max=maximum(s.D for s in trial_ms.gp_damage_states)
                        commit_czm_state!(ms,czm_mesh,trial_ms.u_prev,trial_ms.damage_states,
                            trial_ms.gp_damage_states,mechanical_eigenstrain(load))
                        corrected.iterations = total_iter
                        return corrected
                    end
                end
                if total_iter >= max_total_iter
                    stop_reason = :correction_budget
                    break
                end
                step_size *= 0.5

                if step_size < step_size_min
                    stop_reason = :minimum_step_size
                    break
                end

                @debug "Arc-length substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_residual step_size=step_size
                continue
            end
        end

        K_total, f_int_total, separations, tractions = assemble_coupled_system(
            czm_mesh, u, param; czm_workspace=czm_workspace,
            gp_damage_states=gp_states,
            visc_beta=0.0, czm_model=czm_model)
        costs[:assemblies] += 1

        # On failure report the residual against the requested target load;
        # the accepted continuation point may itself be perfectly balanced.
        R = F_target - f_int_total
        set_czm_residual!(R, u, bc_dofs, bc_vals)
        R_norm = norm(R)

        result.converged = false
        result.iterations = total_iter
        result.residual_norm = R_norm
        result.displacement = u

        fill_czm_result!(result, u, damage_states, separations, tractions)
        if load_progress >= 1.0 - target_neighborhood &&
           load_progress <= 1.0 + 1e-12 &&
           load_progress != last_target_attempt_start &&
           total_iter < max_total_iter &&
           !arc_wall_expired(wall_started_ns, wall_time_limit_seconds)
            # The arc constraint no longer applies at the prescribed final
            # load. A stalled final arc substep may still have reached the
            # target's neighborhood; only the strict fixed-load correction
            # may accept and commit that state.
            corrected, trial_ms = attempt_target_step(u, damage_states, gp_states,
                min(max_iter, max_total_iter - total_iter))
            corrected.iterations += total_iter
            if corrected.converged
                commit_czm_state!(ms,czm_mesh,trial_ms.u_prev,trial_ms.damage_states,
                    trial_ms.gp_damage_states,mechanical_eigenstrain(load))
            end
            return corrected
        end
        @warn "CZM arc-length path did not reach a balanced target load" reason=stop_reason load_progress=load_progress attempts=load_step corrections=total_iter residual=R_norm
        arc_observe(() -> (;phase=:arc_stop,reason=stop_reason,costs=copy(costs),corrections=total_iter))
        return result
    end

"""
    solve_czm_arc_geo_step(czm_mesh, load::MechanicalLoad, param, ms; ...) -> CZMResult

Crisfield 球面弧长 geo 路径（Theory §6.10；λ 缩放本征应变增量）。
约束为 `‖Δu‖² + α²Δλ² = Δl²`。平衡残差采用代码约定
`R = F_ext - f_int`，故修正分解为 `δu = K⁻¹R + δλ K⁻¹f̂`，
其中 `f̂ = ∂R/∂λ = -∂f_int/∂λ` 在每个当前迭代状态重新差分。
接受的率无关子步推进私有损伤/塑性历史；失败试算不晋级。
外部 MechState 只在最终 `λ=1` 平衡验收后发布一次。
"""
function solve_czm_arc_geo_step(czm_mesh::CohesiveMesh, load::MechanicalLoad,
        param, ms::MechState; czm_workspace::CZMWorkspace,
        max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10,
        arc_length_alpha::Float64=1.0,
        visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true,
        eigenstrain_ref=nothing,
        plasticity::Bool=false, prestress=nothing, boundary=nothing)
    F_ext = load.F_ext
    eigenstrain = mechanical_eigenstrain(load)
    isfinite(visc_beta) && visc_beta == 1.0 || throw(ArgumentError(
        "geometric arc_length requires rate-independent visc_beta=1"))
    isfinite(arc_length_alpha) && arc_length_alpha > 0.0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: arc_length_alpha must be finite and positive, got $arc_length_alpha"))
    max_iter > 0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: max_iter must be positive, got $max_iter"))
    n_load_steps > 0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: n_load_steps must be positive, got $n_load_steps"))
    isfinite(tol) && tol > 0.0 || throw(ArgumentError(
        "solve_czm_arc_geo_step: tol must be finite and positive, got $tol"))
    ndof = 2 * size(czm_mesh.node, 1)
    ne = size(czm_mesh.bulk_element, 1)
    n_coh = length(czm_mesh.cohesive_elements)
    result = CZMResult(ndof, n_coh)
    u = copy(ms.u_prev)
    gp_start = gp_history_for_step(ms, czm_mesh)
    any(s -> s.D_visc != s.D, gp_start) && throw(ArgumentError(
        "geometric arc_length cannot advance lagged GP viscous damage without physical-time relaxation"))
    gp_trial_states = similar(gp_start)
    gp_weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
    damage_states = aggregate_gp_damage_states(gp_start, gp_weights)
    committed_plastic_states, trial_plastic_states =
        prepare_plastic_state_buffers(ms, czm_mesh, plasticity)
    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner, boundary=boundary)
    validate_mechanical_boundary_mode(boundary,true)
    zero_bc_vals = zeros(Float64, length(bc_vals))

    successful_load = require_mechanical_load_history(ms)
    if ms.mechanics_initialized
        eigenstrain_ref === nothing || eigenstrain_ref == successful_load || throw(ArgumentError(
            "geometric arc reference differs from the last successful mechanical load"))
        eigenstrain_ref = successful_load
    end
    # A fresh explicit reference is an initialization input. Only a truly
    # fresh state without that input has the defined zero-load reference.
    if eigenstrain_ref === nothing
        eigenstrain_ref = (dT=zeros(ne), Δsn=zeros(ne), Δsp=zeros(ne))
    end
    mix(lam) = (dT=eigenstrain_ref.dT .+ lam .* (eigenstrain.dT .- eigenstrain_ref.dT),
                Δsn=eigenstrain_ref.Δsn .+ lam .* (eigenstrain.Δsn .- eigenstrain_ref.Δsn),
                Δsp=eigenstrain_ref.Δsp .+ lam .* (eigenstrain.Δsp .- eigenstrain_ref.Δsp))
    function assemble_at(ul, lam; gp=gp_start, plastic=committed_plastic_states,
            gp_trial=gp_trial_states, plastic_trial=trial_plastic_states, beta=visc_beta)
        return assemble_coupled_system(czm_mesh, ul, param; czm_workspace=czm_workspace,
            gp_damage_states=gp, gp_trial_states=gp_trial,
            visc_beta=beta, czm_model=czm_model,
            geo_nl=true, plasticity=plasticity, committed_plastic_states=plastic,
            trial_plastic_states=plastic_trial,prestress=prestress,eigenstrain=mix(lam))
    end
    function residual_at(f_int, ul)
        R = F_ext .- f_int
        for dof in bc_dofs
            R[dof] = 0.0
        end
        return R
    end
    function load_direction_at(ul, lam)
        hλ = cbrt(eps(Float64)) * max(1.0, abs(lam))
        f_plus = copy(assemble_at(ul, lam + hλ)[2])
        f_minus = copy(assemble_at(ul, lam - hλ)[2])
        f_hat = -(f_plus .- f_minus) ./ (2.0 * hλ)
        zero_czm_bc_entries!(f_hat, bc_dofs)
        all(isfinite, f_hat) || error(
            "solve_czm_arc_geo_step: non-finite load direction at λ=$lam")
        return f_hat
    end

    # λ=0 参考态可能含卷绕预应力；自由芯部下通常需要多次非线性 Newton，
    # 不能只做一次修正后就判失败。损伤/塑性在此仍为 trial，不提交历史状态。
    K0 = spzeros(Float64, ndof, ndof)
    f0 = zeros(Float64, ndof)
    R0 = fill(Inf, ndof)
    reference_iter = 0
    reference_converged = false
    zero_external = zeros(Float64, ndof)
    if ms.mechanics_initialized
        norm(u[bc_dofs] - bc_vals) <= tol || throw(ArgumentError(
            "geometric arc reference displacement is incompatible with the current Dirichlet boundary"))
        K0,f0,_,_ = assemble_at(u,0.0;gp_trial=nothing,beta=0.0)
        R0 = residual_at(f0,u)
        norm(R0) <= tol || throw(ArgumentError(
            "geometric arc accepted reference is not in equilibrium under the current load/boundary " *
            "(residual=$(norm(R0)), tol=$tol)"))
        if plasticity
            all(zip(trial_plastic_states,committed_plastic_states)) do (trial,committed)
                history_roundoff = 16eps(Float64) * max(1.0,abs(committed.kappa))
                abs(trial.kappa-committed.kappa) <= history_roundoff &&
                    norm(collect(trial.eps_p).-collect(committed.eps_p)) <= history_roundoff
            end || throw(ArgumentError("geometric arc accepted reference would advance plastic history"))
        end
        reference_converged = true
    else
        apply_czm_dirichlet!(u,bc_dofs,bc_vals)
    end
    for _ in 1:(reference_converged ? 0 : max_iter)
        K0, f0, _, _ = assemble_at(u, 0.0)
        R0 = residual_at(f0, u)
        R0_norm = norm(R0)
        isfinite(R0_norm) || error(
            "solve_czm_arc_geo_step: non-finite reference residual")
        if R0_norm <= tol
            reference_converged = true
            break
        end

        K0_bc, R0_bc = apply_bc_czm(
            K0, R0; bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        delta_u = K0_bc \ R0_bc
        all(isfinite, delta_u) || error(
            "solve_czm_arc_geo_step: non-finite reference-state Newton correction")
        u_trial, _, accepted, _, _ = backtrack_line_search!(
            u, delta_u, czm_mesh, param,
            F_ext, zero_external, R0_norm, bc_dofs, bc_vals;
            czm_workspace=czm_workspace,
            max_halvings=12, visc_beta=visc_beta, czm_model=czm_model, geo_nl=true,
            eigenstrain=mix(0.0), plasticity=plasticity,
            committed_plastic_states=committed_plastic_states,
            trial_plastic_states=trial_plastic_states, prestress=prestress,
            gp_damage_states=gp_start, gp_trial_states=gp_trial_states)
        accepted || error(
            "solve_czm_arc_geo_step: reference-state line search failed " *
            "(iteration=$(reference_iter + 1), residual=$R0_norm)")
        u .= u_trial
        reference_iter += 1
    end
    if !reference_converged
        K0, f0, _, _ = assemble_at(u, 0.0)
        R0 = residual_at(f0, u)
        reference_converged = norm(R0) <= tol
    end
    reference_converged || error(
        "solve_czm_arc_geo_step: reference state is not in equilibrium " *
        "after $max_iter iterations (residual=$(norm(R0)))")
    if !ms.mechanics_initialized
        # Explicit fresh initialization can itself create irreversible state.
        # Accept that equilibrated reference privately before following its path.
        assemble_at(u,0.0)
        reference_gp = clone_damage_states(gp_trial_states)
        reference_plastic = plasticity ? clone_plastic_states(trial_plastic_states) : nothing
        K0,f0,_,_ = assemble_at(u,0.0;gp=reference_gp,plastic=reference_plastic,
            gp_trial=nothing,plastic_trial=plasticity ? clone_plastic_states(reference_plastic) : nothing,beta=0.0)
        R0 = residual_at(f0,u)
        norm(R0) <= tol || error(
            "solve_czm_arc_geo_step: committed fresh reference is not in equilibrium (residual=$(norm(R0)))")
        gp_start = reference_gp
        committed_plastic_states = reference_plastic
        damage_states = aggregate_gp_damage_states(gp_start,gp_weights)
    end
    f_hat0 = load_direction_at(u, 0.0)
    K0_bc, _ = apply_bc_czm(K0, zeros(Float64, ndof);
                            bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
    tangent0 = spherical_arc_tangent(K0_bc,f_hat0,nothing,0.0,arc_length_alpha;bc_dofs=bc_dofs)
    tangent0 === nothing && error(
        "solve_czm_arc_geo_step: singular initial tangent has no accepted branch direction")
    tangent_norm0 = tangent0.augmented_norm

    λ = 0.0
    Δl0 = tangent_norm0 / n_load_steps
    Δl = Δl0
    Δl_min = Δl0 / 128.0
    previous_u_tangent = nothing
    previous_lambda_tangent = 0.0
    total_iter = reference_iter
    step_count = 0
    residual_history = Float64[]
    lambda_history = Float64[λ]
    step_history = Float64[]
    arc_step_tol = 10.0 * tol  # 中间路径点容差；最终 λ=1 平衡仍严格使用 tol

    while λ < 1.0 - 1e-10
        step_count += 1
        step_count <= 10000 || error(
            "solve_czm_arc_geo_step: exceeded 10000 arc steps " *
            "(λ_history=$lambda_history, residual_history=$residual_history)")
        u_start = copy(u)
        λ_start = λ

        K_start, _, _, _ = assemble_at(u_start, λ_start)
        f_hat_start = load_direction_at(u_start, λ_start)
        K_start_bc, _ = apply_bc_czm(K_start, zeros(Float64, ndof);
                                     bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        tangent = spherical_arc_tangent(K_start_bc,f_hat_start,previous_u_tangent,
            previous_lambda_tangent,arc_length_alpha;bc_dofs=bc_dofs)
        tangent === nothing && error(
            "solve_czm_arc_geo_step: no usable augmented predictor at λ=$λ_start")
        predictor_scale = Δl / tangent.augmented_norm
        dλ_pred = predictor_scale * tangent.delta_lambda
        if dλ_pred > 1.0 - λ_start
            predictor_scale *= (1.0 - λ_start) / dλ_pred
            dλ_pred = 1.0 - λ_start
        end
        Δl_step = predictor_scale * tangent.augmented_norm
        Δl_step > 0.0 || error(
            "solve_czm_arc_geo_step: zero predictor step at λ=$λ_start")
        u .= u_start .+ predictor_scale .* tangent.delta_u
        λ = λ_start + dλ_pred
        apply_czm_dirichlet!(u, bc_dofs, bc_vals)

        step_ok = false
        last_residual = Inf
        last_step_error = nothing
        for _ in 1:max_iter
            total_iter += 1
            K, f_int, _, _ = assemble_at(u, λ)
            R = residual_at(f_int, u)
            Δu_bar = u .- u_start
            Δλ_bar = λ - λ_start
            g = dot(Δu_bar, Δu_bar) +
                arc_length_alpha^2 * Δλ_bar^2 - Δl_step^2
            last_residual = norm(R)
            isfinite(last_residual) || error(
                "solve_czm_arc_geo_step: non-finite substep residual")
            push!(residual_history, last_residual)
            if last_residual <= arc_step_tol &&
               abs(g) <= arc_step_tol * max(Δl_step^2, eps(Float64))
                step_ok = true
                break
            end

            K_bc, R_bc = apply_bc_czm(
                K, R; bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
            f_hat = load_direction_at(u, λ)
            correction = spherical_arc_corrector(K_bc,R_bc,f_hat,Δu_bar,Δλ_bar,
                arc_length_alpha,Δl_step;bc_dofs=bc_dofs)
            if correction === nothing
                last_step_error = "singular augmented correction"
                break
            end
            Δu,dλ,_ = correction
            u .+= Δu
            λ += dλ
            apply_czm_dirichlet!(u, bc_dofs, bc_vals)
        end

        # Reaching/overshooting the endpoint supplies only a Newton guess.
        # Its history is accepted only after the strict final tolerance gate.
        step_ok && λ >= 1.0 && break

        if step_ok
            candidate_gp = clone_damage_states(gp_trial_states)
            candidate_plastic = plasticity ? clone_plastic_states(trial_plastic_states) : nothing
            _,f_committed,_,_ = assemble_at(u,λ;gp=candidate_gp,plastic=candidate_plastic,
                gp_trial=nothing,plastic_trial=plasticity ? clone_plastic_states(candidate_plastic) : nothing,
                beta=0.0)
            last_residual = norm(residual_at(f_committed,u))
            if last_residual <= arc_step_tol
                gp_start = candidate_gp
                committed_plastic_states = candidate_plastic
                damage_states = aggregate_gp_damage_states(gp_start,gp_weights)
            else
                step_ok = false
                last_step_error = "committed substep is not in equilibrium"
            end
        end
        if step_ok
            Δu_step = u .- u_start
            Δλ_step = λ - λ_start
            previous_norm = hypot(norm(Δu_step),arc_length_alpha * Δλ_step)
            previous_norm > 0.0 || error(
                "solve_czm_arc_geo_step: accepted a zero arc step at λ=$λ")
            previous_u_tangent = Δu_step ./ previous_norm
            previous_lambda_tangent = Δλ_step / previous_norm
            push!(lambda_history, λ)
            push!(step_history, Δl_step)
            observer = mechanical_step_observer[]
            if observer !== nothing
                observer((phase=:arc_geo_substep,lambda=λ,u=copy(u),
                    gp_damage_states=clone_damage_states(gp_start),
                    plastic_states=plasticity ? clone_plastic_states(committed_plastic_states) : nothing,
                    residual=last_residual,arc_radius=Δl_step,ms=ms))
            end
            continue
        end

        u .= u_start
        λ = λ_start
        Δl *= 0.5
        if Δl < Δl_min
            error("solve_czm_arc_geo_step: arc stepping stalled at λ=$λ " *
                  "(residual=$last_residual, Δl=$Δl, λ_history=$lambda_history, " *
                  "step_error=$(repr(last_step_error)), step_history=$step_history, " *
                  "residual_history=$residual_history)")
        end
    end

    # 精确落到 λ=1，并以固定 λ Newton 复算最终平衡；不能用弧长预测残差代替。
    λ = 1.0
    final_residual = Inf
    separations = Vector{Tuple{Float64,Float64}}(undef, n_coh)
    tractions = Vector{Tuple{Float64,Float64}}(undef, n_coh)
    for iter in 1:(max_iter + 1)
        iter <= max_iter && (total_iter += 1)
        K, f_int, separations, tractions = assemble_at(u, λ)
        R = residual_at(f_int, u)
        final_residual = norm(R)
        push!(residual_history, final_residual)
        final_residual <= tol && break
        iter > max_iter && break
        K_bc, R_bc = apply_bc_czm(K, R;
            bc_dofs=bc_dofs, bc_vals=zero_bc_vals)
        Δu = K_bc \ R_bc
        all(isfinite, Δu) || error(
            "solve_czm_arc_geo_step: non-finite final λ=1 correction")
        u .+= Δu
        apply_czm_dirichlet!(u, bc_dofs, bc_vals)
    end
    final_residual <= tol || error(
        "solve_czm_arc_geo_step: final λ=1 equilibrium did not converge " *
        "(residual=$final_residual, λ_history=$lambda_history, " *
        "step_history=$step_history, residual_history=$residual_history)")

    _, f_final, separations, tractions = assemble_at(u, 1.0)
    final_residual = norm(residual_at(f_final, u))
    final_residual <= tol || error(
        "solve_czm_arc_geo_step: final residual recomputation failed (residual=$final_residual)")
    # load_direction_at overwrites the shared trial buffer at λ±h; only the
    # actual final λ=1 assembly immediately above defines the accepted GP state.
    committed_gp = clone_damage_states(gp_trial_states)
    damage_states = aggregate_gp_damage_states(committed_gp, gp_weights)
    committed_plastic = plasticity ? clone_plastic_states(trial_plastic_states) : ms.plastic_states
    _, f_committed, separations, tractions = assemble_at(u,1.0;
        gp=committed_gp,plastic=plasticity ? committed_plastic : nothing,gp_trial=nothing,
        plastic_trial=plasticity ? clone_plastic_states(committed_plastic) : nothing,beta=0.0)
    final_residual = norm(residual_at(f_committed,u))
    final_residual <= tol || error(
        "solve_czm_arc_geo_step: committed damage is not in equilibrium (residual=$final_residual)")
    result.converged = true
    result.iterations = total_iter
    result.residual_norm = final_residual
    fill_czm_result!(result, u, damage_states, separations, tractions)
    commit_czm_state!(ms,czm_mesh,u,damage_states,committed_gp,mix(1.0);
        plastic_states=committed_plastic)
    return result
end

