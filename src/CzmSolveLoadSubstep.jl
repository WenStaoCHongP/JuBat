# ========================================================================
# 7. Newton-Raphson solver
# ========================================================================

"""
    newton_raphson_czm(czm_mesh, load::MechanicalLoad, param, ms; max_iter=50, tol=1e-8, n_load_steps=10)

Newton-Raphson nonlinear solver with load substeps.

# Returns
- `result`: CZMResult
- `new_czm_mesh`: updated CZM mesh with damage states
"""
function newton_raphson_czm(czm_mesh::CohesiveMesh, load::MechanicalLoad, param, ms::MechState; czm_workspace::CZMWorkspace, max_iter::Int=50, tol::Float64=1e-8, n_load_steps::Int=10, visc_beta::Float64=1.0, czm_model::String="model1", fix_inner::Bool=true, geo_nl::Bool=false, plasticity::Bool=false, prestress=nothing, boundary=nothing)
    F_ext = load.F_ext
    eigenstrain = mechanical_eigenstrain(load)
    nnode = size(czm_mesh.node, 1)
    ndof = 2 * nnode
    n_coh = length(czm_mesh.cohesive_elements)

    result = CZMResult(ndof, n_coh)

    u = copy(ms.u_prev)
    gp_start = gp_history_for_step(ms, czm_mesh)
    gp_trial_states = similar(gp_start)
    gp_weights = NCweight(czm_mesh.czm_submesh.mesh_bonded.gs.order)[1]
    damage_states = aggregate_gp_damage_states(gp_start, gp_weights)
    committed_plastic_states, trial_plastic_states =
        prepare_plastic_state_buffers(ms, czm_mesh, plasticity)

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; fix_inner=fix_inner, boundary=boundary)
        boundary === nothing || apply_czm_dirichlet!(u,bc_dofs,bc_vals)
    validate_mechanical_boundary_mode(boundary,geo_nl)

    # geo_nl（Batch 2，D-B2-1）：ε* 内嵌 f_int^GL，F_tc 不再外载；切线依赖 u，禁用缓存
    if geo_nl
        F_thermo_chem_total = zeros(Float64, ndof)
    else
        F_thermo_chem_total = assemble_thermal_chemical_load(czm_mesh, param, load.dT, load.Δsoc_n, load.Δsoc_p)
    end
    eig_kwargs = geo_nl ? (
        geo_nl=true, eigenstrain=eigenstrain, plasticity=plasticity,
        committed_plastic_states=committed_plastic_states,
        trial_plastic_states=trial_plastic_states,
        prestress=prestress) : ()

    # 增量载荷参考：u_prev 近似在上一时间步的平衡态
    # f_int(u_prev) ≈ 上一步外力，F_delta = 目标载荷 - 平衡内力（增量，通常很小）
    _, f_int_ref, _, _ = assemble_coupled_system(czm_mesh, u, param; czm_workspace=czm_workspace,gp_damage_states=gp_start, visc_beta=0.0, czm_model=czm_model, eig_kwargs...)
    F_target = F_ext + F_thermo_chem_total
    check_czm_external_load(bc_dofs,F_target)
    F_delta = F_target - f_int_ref

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
        # 从平衡态 f_int_ref 逐步推进到目标态 F_target
        F_applied = f_int_ref + target_progress * F_delta

        u_start = copy(u)
        converged_substep = false
        last_R_norm = Inf
        cached_assembly = nothing

        for iter in 1:(max_iter + 1)
            iter <= max_iter && (total_iter += 1)

            if cached_assembly === nothing
                K_total, f_int_total, separations, tractions = assemble_coupled_system(
                    czm_mesh, u, param; czm_workspace=czm_workspace,
                    gp_damage_states=gp_start, gp_trial_states=gp_trial_states,
                    visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)
            else
                K_total, f_int_total, separations, tractions = cached_assembly
                cached_assembly = nothing
            end

            R = F_applied - f_int_total

            set_czm_residual!(R, u, bc_dofs, bc_vals)

            R_norm = norm(R)
            last_R_norm = R_norm
            if R_norm < tol
                converged_substep = true
                load_progress = target_progress
                step_size = min(step_size * 1.25, step_size_max)
                break
            end

            # The last pass only judges the last permitted correction.
            iter > max_iter && break

            K_bc, R_bc = apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
            # R1-d：geo 路径与 basic 统一走对角均衡化（装配切线量级失衡，未均衡化
            # 分解填充 ~69×；非 geo 保持原路径）。失败语义不变：奇异即抛出
            Δu = czm_linear_solve(K_bc,R_bc,bc_dofs;equilibrated=geo_nl)

            if any(isnan, Δu) || any(isinf, Δu)
                break
            end

            α = 1.0
            ls_accepted = false
            for _ in 1:8
                u_trial = u + α * Δu
                apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

                K_trial, f_int_trial, sep_trial, trac_trial = assemble_coupled_system(
                    czm_mesh, u_trial, param; czm_workspace=czm_workspace,
                    gp_damage_states=gp_start, gp_trial_states=gp_trial_states,
                    visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

                R_trial = F_applied - f_int_trial
                set_czm_residual!(R_trial, u_trial, bc_dofs, bc_vals)

                R_trial_norm = norm(R_trial)
                if !isnan(R_trial_norm) && R_trial_norm < R_norm
                    # 同一子步内输入不变；下一轮在任何新装配前消费响应。
                    # GP/塑性 trial 缓冲也保留在此接受点，不提交为历史。
                    u = u_trial
                    cached_assembly = (K_trial, f_int_trial, sep_trial, trac_trial)
                    ls_accepted = true
                    break
                end

                α *= 0.5
            end

            if !ls_accepted
                break
            end
        end

        if !converged_substep
            u = u_start
            step_size *= 0.5

            if step_size < step_size_min
                @warn "CZM adaptive load stepping stalled" load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
                break
            end

            @debug "Load substep $load_step failed, reducing step size and retrying..." load_progress=load_progress target_progress=target_progress residual=last_R_norm step_size=step_size
            continue
        end
    end

    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, param; czm_workspace=czm_workspace,gp_damage_states=gp_start, gp_trial_states=gp_trial_states, visc_beta=visc_beta, czm_model=czm_model, eig_kwargs...)

    F_applied_final = f_int_ref + load_progress * F_delta
    R = F_applied_final - f_int_total
    set_czm_residual!(R, u, bc_dofs, bc_vals)
    R_norm = norm(R)

    result.converged = load_progress >= 1.0 - 1e-12 && R_norm < tol
    result.iterations = total_iter
    result.residual_norm = R_norm
    result.displacement = u

    # 子步仅推进位移；所有试算均来自同一外部接受史，仅发布最终端点的GP候选。
    if result.converged
        committed_gp = clone_damage_states(gp_trial_states)
        committed_damage = aggregate_gp_damage_states(committed_gp, gp_weights)
        _, committed_force, committed_sep, committed_trac = assemble_coupled_system(
            czm_mesh,u,param; czm_workspace=czm_workspace,gp_damage_states=committed_gp,
            visc_beta=0.0,czm_model=czm_model,eig_kwargs...)
        committed_R = F_target-committed_force
        set_czm_residual!(committed_R,u,bc_dofs,bc_vals)
        result.residual_norm = norm(committed_R)
        result.converged = result.residual_norm < tol
        if result.converged
            damage_states = committed_damage
            gp_start = committed_gp
            separations,tractions = committed_sep,committed_trac
        end
    end
    if !result.converged
        u = copy(ms.u_prev)
        damage_states = clone_damage_states(ms.damage_states)
        _,_,separations,tractions = assemble_coupled_system(czm_mesh,u,param; czm_workspace=czm_workspace,
            gp_damage_states=gp_start,
            visc_beta=0.0,czm_model=czm_model,eig_kwargs...)
        result.displacement = u
    end

    fill_czm_result!(result, u, damage_states, separations, tractions)

    if result.converged
        commit_czm_state!(ms,czm_mesh,u,damage_states,gp_start,eigenstrain;
            plastic_states=plasticity ? trial_plastic_states : ms.plastic_states)
    end
    return result
end

