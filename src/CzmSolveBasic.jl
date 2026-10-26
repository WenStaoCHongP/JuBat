function solve_czm_basic_step(czm_mesh::CohesiveMesh, load::MechanicalLoad, param, ms::MechState; czm_workspace::CZMWorkspace, max_iter::Int=50, tol::Float64=1e-8, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, plasticity::Bool=false, prestress=nothing, boundary=nothing)
        F_ext = load.F_ext
        eigenstrain = mechanical_eigenstrain(load)
        nnode = size(czm_mesh.node, 1)
        ndof = 2 * nnode
        n_coh = length(czm_mesh.cohesive_elements)

        result = CZMResult(ndof, n_coh)
        u = copy(ms.u_prev)
        u_start = copy(u)
        gp_states = gp_history_for_step(ms, czm_mesh)
        gp_trial_states = similar(gp_states)
        gp_weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
        damage_states = aggregate_gp_damage_states(gp_states, gp_weights)
        committed_plastic_states, trial_plastic_states =
            prepare_plastic_state_buffers(ms, czm_mesh, plasticity)

        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner, boundary=boundary)
        boundary === nothing || apply_czm_dirichlet!(u,bc_dofs,bc_vals)
        validate_mechanical_boundary_mode(boundary,geo_nl)

        # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
        if geo_nl
            F_thermo_chem = zeros(Float64, ndof)
        else
            F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, param, load.dT, load.Δsoc_n, load.Δsoc_p)
        end
        check_czm_external_load(bc_dofs,F_ext+F_thermo_chem)
        eig_kwargs = geo_nl ? (
            geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity,
            committed_plastic_states=committed_plastic_states,
            trial_plastic_states=trial_plastic_states,
            prestress=prestress) : ()

        total_iter = 0
        R_norm = Inf
        failure_R_norm = Inf
        converged = false
        converged_R_norm = Inf
        separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)
        cached_K = nothing
        cached_f_int = nothing
        cached_separations = nothing
        cached_tractions = nothing

        # One final assembly after the last Newton correction is needed to
        # judge that correction; it does not grant an extra Newton step.
        for iter in 1:(max_iter + 1)
            iter <= max_iter && (total_iter += 1)

            # 线搜索接受的试探装配在相同输入下与入口装配逐位一致，可直接复用；
            # 缓存生命周期内唯一的 ws 写者是线性求解的因子缓存字段（与装配
            # 输出别名 disjoint），任何 break 路径都不再消费缓存
            if cached_K === nothing
                K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param; czm_workspace=czm_workspace,gp_damage_states=gp_states, gp_trial_states=gp_trial_states, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)
            else
                K_total = cached_K
                f_int_total = cached_f_int
                separations = cached_separations
                tractions = cached_tractions
                cached_K = nothing
                cached_f_int = nothing
                cached_separations = nothing
                cached_tractions = nothing
            end

            R = F_ext + F_thermo_chem - f_int_total

            set_czm_residual!(R, u, bc_dofs, bc_vals)

            R_norm = norm(R)
            failure_R_norm = R_norm
            if R_norm < tol
                trial_gp_states = clone_damage_states(gp_trial_states)
                trial_damage_states = aggregate_gp_damage_states(trial_gp_states, gp_weights)
                _, f_int_committed, sep_committed, trac_committed =
                    assemble_coupled_system(czm_mesh, u, param; czm_workspace=czm_workspace,
                        gp_damage_states=trial_gp_states,
                        visc_beta=0.0,
                        czm_model=czm_model,
                        eig_kwargs...)
                R_committed = F_ext + F_thermo_chem - f_int_committed
                set_czm_residual!(R_committed, u, bc_dofs, bc_vals)
                committed_R_norm = norm(R_committed)
                failure_R_norm = committed_R_norm
                if committed_R_norm < tol
                    damage_states = trial_damage_states
                    gp_states = trial_gp_states
                    converged = true
                    converged_R_norm = committed_R_norm
                    separations, tractions = sep_committed, trac_committed
                    break
                end
                break
            end

            iter > max_iter && break
            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)

            Δu = try
                if geo_nl
                    solve_equilibrated(K_bc, R_bc)
                else
                    czm_linear_solve(K_bc, R_bc, bc_dofs; czm_workspace=czm_workspace)
                end
            catch err
                err isa Union{SingularException,ZeroPivotException} || rethrow()
                break
            end

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            u, R_norm, ls_accepted, α_used, ls_assembly = backtrack_line_search!(u, Δu, czm_mesh, param, F_ext, F_thermo_chem, R_norm,bc_dofs, bc_vals; czm_workspace=czm_workspace,visc_beta=visc_beta, czm_model=czm_model, geo_nl=geo_nl, eigenstrain=eigenstrain, plasticity=plasticity, committed_plastic_states=committed_plastic_states, trial_plastic_states=trial_plastic_states, prestress=prestress, gp_damage_states=gp_states, gp_trial_states=gp_trial_states)
            failure_R_norm = R_norm

            if !ls_accepted
                break
            end
            cached_K, cached_f_int, cached_separations, cached_tractions = ls_assembly
        end

        if !converged
            u = u_start
            damage_states = clone_damage_states(ms.damage_states)   # 未收敛不触碰 ms（试探态丢弃）
        end

        if !converged
            _, f_int_total, separations, tractions = assemble_coupled_system(
                czm_mesh, u, param; czm_workspace=czm_workspace,
                gp_damage_states=gp_states,
                visc_beta=0.0,
                czm_model=czm_model, eig_kwargs...)
            R = F_ext + F_thermo_chem - f_int_total
            set_czm_residual!(R, u, bc_dofs, bc_vals)
            R_norm = norm(R)
        end

        # 失败时位移/历史回滚，但报告最后试算残差；回滚态可能恰好
        # 在旧损伤下平衡，其小残差不能解释本次求解为何失败。
        final_R_norm = converged ? converged_R_norm : failure_R_norm

        result.converged = converged
        result.iterations = total_iter
        result.residual_norm = final_R_norm
        result.displacement = u
        fill_czm_result!(result, u, damage_states, separations, tractions)
        if converged
            commit_czm_state!(ms,czm_mesh,u,damage_states,gp_states,eigenstrain;
                plastic_states=plasticity ? trial_plastic_states : ms.plastic_states)
        end
        return result
    end

