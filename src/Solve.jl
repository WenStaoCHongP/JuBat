function Solve(case::Case;initial_state::Union{NamedTuple,Nothing}=nothing,return_final_state::Bool=false,thermal_variables::Union{Dict{String,Any},Nothing}=nothing,thermal_update_fn::Union{Function,Nothing}=nothing,thermal_record::Bool=false,polar_mesh_data::Any=nothing,czm_snapshots::Union{Vector{CZMSnapshot},Nothing}=nothing,czm_cycle::Int=1,czm_phase::String="unknown",phase_type::Union{Nothing,PhaseType}=nothing,chemical_step_observer::Union{Nothing,Function}=nothing,variables::Dict{String,Any}=Dict{String,Any}())
    check_continuous_feedback(case.opt.czm)
    if case.opt.model == "thermal"
            t0 = case.opt.time[1]/case.param.scale.t0
            t_end = case.opt.time[end]/case.param.scale.t0
            dt = case.opt.dt[1]/case.param.scale.t0
            vars = thermal_variables
            update_fn = thermal_update_fn
            record = thermal_record

            if case.opt.solveType == "Crank-Nicolson"
                theta = 0.5
            elseif case.opt.solveType == "forward"
                theta = 0.0
            elseif case.opt.solveType == "backward"
                theta = 1.0
            else
                error("Error: $(case.opt.solveType) difference scheme has not been implemented!\n")
            end

            mesh = case.mesh["thermal2D"]
            nnode = mesh.nlen
            T_nodes = vars["thermal2D temperature at nodes"]

            times = collect(range(t0, step=dt, length=Int(cld(t_end - t0, dt)) + 1))
            T_hist = if record
                zeros(Float64, nnode, length(times))
            else
                zeros(Float64, 0, 0)
            end
            record && (T_hist[:, 1] .= T_nodes)

            update_fn !== nothing && update_fn(t0, vars)
            if case.opt.thermalmodel == "ring2D_polar"
                MT_old, KT_old, FT_old = ThermalPolar2D_Ring(case, vars, polar_mesh_data)
            elseif case.opt.thermalmodel == "ring2D"
                MT_old, KT_old, FT_old = ThermalDistributed2D_Ring(case, vars)
                KT_old, FT_old = ThermalRing2D_BC(KT_old, FT_old, case, vars["thermal2D outer_nodes"], t0)
            elseif case.opt.thermalmodel == "distributed2D"
                MT_old, KT_old, FT_old = ThermalDistributed2D(case, vars)
                KT_old, FT_old = ThermalDistributed2D_BC(KT_old, FT_old, case, t0)
            else
                throw(ArgumentError("Unsupported pure thermal model: $(case.opt.thermalmodel)"))
            end

            for step in 1:(length(times) - 1)
                t = times[step + 1]
                update_fn !== nothing && update_fn(t, vars)

                if case.opt.thermalmodel == "ring2D_polar"
                    MT_new, KT_new, FT_new = ThermalPolar2D_Ring(case, vars, polar_mesh_data)
                elseif case.opt.thermalmodel == "ring2D"
                    MT_new, KT_new, FT_new = ThermalDistributed2D_Ring(case, vars)
                    KT_new, FT_new = ThermalRing2D_BC(KT_new, FT_new, case, vars["thermal2D outer_nodes"], t)
                elseif case.opt.thermalmodel == "distributed2D"
                    MT_new, KT_new, FT_new = ThermalDistributed2D(case, vars)
                    KT_new, FT_new = ThermalDistributed2D_BC(KT_new, FT_new, case, t)
                else
                    throw(ArgumentError("Unsupported pure thermal model: $(case.opt.thermalmodel)"))
                end

                A = MT_new - theta * dt * KT_new
                rhs = (MT_new + (1.0 - theta) * dt * KT_old) * T_nodes +
                      dt * (theta * FT_new + (1.0 - theta) * FT_old)
                T_nodes = A \ rhs
                vars["thermal2D temperature at nodes"] = T_nodes

                record && (T_hist[:, step + 1] .= T_nodes)

                KT_old = KT_new
                FT_old = FT_new
            end

            result = (time = times, T_nodes = T_nodes, T_hist = T_hist)
            return result
        end
    # Every phase reads the same physical-voltage window; never overwrite param.
    v_l, v_h = case.param.cell.v_l, case.param.cell.v_h
    isfinite(v_l) && isfinite(v_h) && v_l < v_h ||
        throw(ArgumentError("invalid param.cell.v_l/v_h voltage window"))
        # Electrochemical time still scaled by t0
        dt_min = case.opt.dt[1] / case.param.scale.t0
        dt_max = case.opt.dt[2] / case.param.scale.t0
        RunTime = case.opt.time / case.param.scale.t0
        t0 = RunTime[1]
        t_end = RunTime[end]

    # 判断是否启用多SPMe模式：与 CallModel 保持一致，由 per_element_spme 控制
    multi_spme_enabled = case.opt.per_element_spme

    if case.opt.czm.enabled
        case.mesh["czm"] # Enabled mechanics requires its real topology.
    end
    global_time = if initial_state === nothing
        case.opt.time[1]
    else
        initial_state.t_global
    end
    phase_time_origin=global_time
    y0 = if initial_state === nothing
        copy(if case.opt.per_element_spme
            ModelInitialisation_MultiSPMe(case)
        else
            ModelInitialisation(case)
        end)
    else
        copy(initial_state.y)
    end
    mech0 = if initial_state === nothing
        if haskey(case.mesh,"czm")
            MechState(case.mesh["czm"])
        else
            nothing
        end
    else
        deepcopy(initial_state.mech)
    end
    y0 isa Vector{Float64} || throw(ArgumentError("restored chemistry must be a Float64 vector"))
    all(isfinite,y0) || throw(ArgumentError("restored chemistry is non-finite"))
    isfinite(global_time) || throw(ArgumentError("non-finite restored global time"))
    if haskey(case.mesh,"czm")
        mech0===nothing && throw(ArgumentError("CZM topology requires mechanics"))
    else
        mech0===nothing || throw(ArgumentError("mechanics requires CZM topology"))
    end
    a=(y=y0,mech=mech0,t_global=global_time)
    mech0===nothing || require_mechanical_load_history(mech0)
    # Scheduling is phase-local; every accepted physical interval is consumed.
    pending_dt=0.0
    mechanical_step_count=0
    if case.opt.solveType == "Crank-Nicolson"
        theta = 0.5
    elseif case.opt.solveType == "forward"
        theta = 0
    elseif case.opt.solveType == "backward"
        theta = 1
    else
        error( "Error: $(case.opt.solveType) difference scheme has not been implemented!\n ")
    end

    dt = dt_min
    # 计算预期时间步数，添加安全限制避免内存溢出
    num_estimated = round(Int64, (t_end - t0)/dt * 1.5)
    # 限制最大预分配步数（避免内存溢出）
    max_steps = if multi_spme_enabled
        50000
    else
        100000
    end
    num = min(num_estimated+2, max_steps) # initial column and real BE interval

    if num_estimated > max_steps
        @warn "预期时间步数 $(num_estimated) 超过最大限制 $(max_steps)，将使用动态扩展策略"
        @warn "这可能导致性能下降。建议增大时间步长 dt_min = $(dt*case.param.scale.t0) 秒"
    end

    history = StandardVariables(case, num)
    czm_workspace = if haskey(case.mesh,"czm")
        CZMWorkspace(case.mesh["czm"],case.param)
    else
        nothing
    end

    t = t0
    vt = 2
    v = 1
    timing_totals = Dict{String,Float64}("spme" => 0.0,"branch" => 0.0,"thermal" => 0.0,"czm" => 0.0,)
    timing_call_count = 0
    # 5.7.4 返工 P2-2：独立候选/接受/eval 计数器（不复用输出列 v）
    chem_trial_idx = 0
    chem_accepted_idx = 0

    function accumulate_callmodel_timing!(totals::Dict{String,Float64}, vars::Dict{String, Union{Array{Float64},Float64}})
        totals["spme"] += get(vars, "timing spme solve [s]", 0.0)
        totals["branch"] += get(vars, "timing branch solver [s]", 0.0)
        totals["thermal"] += get(vars, "timing thermal distributed [s]", 0.0)
        totals["czm"] += get(vars, "timing czm model [s]", 0.0)
        return nothing
    end

    M_old, K_old, F_old, step_variables, y_phi= CallModel(case, y0, t; jacobi="update",mech=a.mech)
    timing_call_count += 1
    accumulate_callmodel_timing!(timing_totals, step_variables)
    # 打印初始信息
    V_init = step_variables["cell voltage"] * case.param.scale.phi
    println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")
    # 持久化热场（跨 CallModel 迭代携带）
    T_nodes_carry = if multi_spme_enabled
        T_nodes_initial = step_variables["thermal2D temperature at nodes"]
        T_nodes_initial = vec(T_nodes_initial)
        length(T_nodes_initial) == case.mesh["thermal2D"].nlen || throw(DimensionMismatch("initial thermal node count $(length(T_nodes_initial)) does not match layout count $(case.mesh["thermal2D"].nlen)"))
        copy(T_nodes_initial)
    else
        Float64[]
    end

    # CZM 与宏观应力状态。latest_macro_stress 只在实际力学求解后刷新，
    # 其间的输出列保持最近一次有效解，禁止把预分配零当作应力结果。
    czm_active = case.opt.czm.enabled
    macro_stress_active = czm_active && (case.opt.czm.j2_plasticity || !(case.opt.czm.geo_nonlinear || case.opt.czm.winding_prestress))
    latest_macro_stress = if !macro_stress_active
        nothing
    elseif !a.mech.mechanics_initialized
        nothing
    else
        compute_macro_stress(case, a.mech, step_variables, T_nodes_carry)
    end

    vc = 1:size(M_old,1)
    isfinite(step_variables["cell voltage"]*case.param.scale.phi) || error("non-finite initial terminal voltage")
    # The BE startup changes physical chemistry/temperature, so it is a real
    # first interval. Keep the original state/history at the actual start.
    dt_init = min(1e-8,RunTime[vt]-t0,t_end-t0)
    dt_init = min(t0+dt_init,RunTime[vt],t_end)-t0
    initial_trial = if dt_init > 0.0
        y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
        all(isfinite,y_c) || error("non-finite initial numerical state")
        if chemical_step_observer !== nothing
            chemical_step_observer((; kind=:init_BE,h=dt_init,t_eval=t,
                t_start=t,t_end=t+dt_init,t_load=t,
                M=M_old,K=K_old,F=F_old,y_before=y0[vc],y_after=y_c,step_variables))
        end
        vcat(y_c,y_phi)
    else
        nothing
    end
    latest_macro_stress===nothing || write_macro_stress!(step_variables,1,latest_macro_stress)
    candidate=(y=copy(y0),mech=a.mech,t_global=phase_time_origin)
    Variable_update!(history,step_variables,v;debug_coupling=case.opt.debug_coupling)
    variables["history"]=history
    czm_workspace===nothing || (variables["czm_workspace"]=czm_workspace)
    variables["accepted"]=candidate
    y_old=candidate.y
    first_step=true
    if case.opt.jacobi == "constant"
        RecordMatrix!(case, M_old, K_old)    # record system matrix information
    end

    print( "start to solve the problem \n")

    # 单元截止追踪变量
    first_cutoff_detected = false
    first_cutoff_time = 0.0
    first_cutoff_element = 0
    first_cutoff_ocv = 0.0
    total_cutoff_count = 0
    termination_reason = "time_limit"  # 默认终止原因

    # run the model
    while t < t_end
        # 电化学步
        # 5.7.4 返工 P2-2：候选开始（独立计数器）
        y_old=variables["accepted"].y
        t_start=t
        t_candidate=min(t+(first_step ? dt_init : dt),RunTime[vt],t_end)
        h_used=t_candidate-t_start
        t_eval=first_step ? t_start : t_candidate
        theta_used=first_step ? 1.0 : theta
        chem_trial_idx += 1
        if chemical_step_observer !== nothing
            chemical_step_observer((; kind = :attempt_begin, h = h_used, t_eval,
                t_start,t_end=t_candidate,t_load=t_start,
                trial_idx = chem_trial_idx))
        end
        if first_step
            M_new,K_new,F_new=M_old,K_old,F_old
        else
            M_new, K_new, F_new, step_variables, y_phi = CallModel(case, y_old, t_eval; jacobi="update",mech=variables["accepted"].mech)
            timing_call_count += 1
            accumulate_callmodel_timing!(timing_totals, step_variables)
        end
        # SOC/temperature are sampled from the accepted left state; forcing
        # uses t_eval. Preserve that explicit load before recording new T.
        T_nodes_czm_current = if multi_spme_enabled
            copy(vec(step_variables["thermal2D temperature at nodes"]))
        else
            T_nodes_carry
        end
        Mt = M_new - theta_used * K_new * h_used
        Kt = (1 - theta_used) * K_old * h_used + M_new
        Ft = theta_used * F_new * h_used + (1 - theta_used) * F_old * h_used
        y_c = first_step ? initial_trial[vc] : convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old[vc] + Ft)
        y_new = first_step ? initial_trial : vcat(y_c, y_phi)
        all(isfinite,y_new) || error("non-finite numerical state")
        isfinite(step_variables["cell voltage"]*case.param.scale.phi) || error("non-finite terminal voltage")

        # multi-SPMe：求解后提取温度自由度，用于记录/后处理
        if multi_spme_enabled
            nT = case.mesh["thermal2D"].nlen
            expected_multi_len = size(case.mesh["thermal2D"].element, 1) * case.index["SPMe chemical block length"] + nT
            length(y_c) == expected_multi_len || throw(DimensionMismatch("time-step state length $(length(y_c)) does not match multi-SPMe layout length $expected_multi_len"))
            T_nodes = get_thermal_dofs(y_c, case)
            step_variables["thermal2D temperature at nodes"] = T_nodes
        end
        error_y = ErrorEstimation(case, y_old, y_new, dt_min/h_used)
        if !first_step && error_y > 2 * case.opt.dtThreshold && case.opt.dtType == "auto" && t_candidate >= t_start + dt_min * 4
            # reduce dt to dt/2 and recaculate y_new
            # 5.7.4 返工 P2-2：拒绝（独立计数，不计入库存/接受账本）
            if chemical_step_observer !== nothing
                chemical_step_observer((; kind = :rejected, h_used = h_used,
                    t_eval,t_start,t_end=t_candidate,t_load=t_start,
                    voltage=step_variables["cell voltage"]*case.param.scale.phi,
                    error_y = error_y, trial_idx = chem_trial_idx))
            end
            dt = h_used/2
        else
            # Cutoff decisions belong to accepted candidates. Determine the
            # terminal interval before mechanics so its pending time is part
            # of this same candidate transaction.
            V_cell=step_variables["cell voltage"]*case.param.scale.phi
            area_diag=multi_spme_enabled && step_variables["thermal2D Vsolve status"]==4.0
            candidate_reason="time_limit"
            if multi_spme_enabled && !area_diag
                n_cutoff=Int(step_variables["thermal2D n_cutoff_elements"])
                ne_total=size(case.mesh["thermal2D"].element,1)
                if n_cutoff>0 && !first_cutoff_detected
                    first_cutoff_detected=true
                    first_cutoff_time=t_candidate*case.param.scale.t0
                    first_cutoff_element=Int(step_variables["thermal2D cutoff_elements"][1])
                    first_cutoff_ocv=step_variables["thermal2D cutoff_ocv"][1]
                end
                total_cutoff_count=n_cutoff
                if ne_total>0 && n_cutoff>=ne_total
                    candidate_reason="all_elements_cutoff"
                elseif ne_total>0 && 0<n_cutoff<ne_total &&
                       Int(step_variables["thermal2D n_active_elements"])==0
                    I_now=step_variables["cell current"]
                    if I_now!=0.0
                        candidate_reason=I_now>0 ? "voltage_cutoff_low" : "voltage_cutoff_high"
                    end
                end
            end
            if candidate_reason=="time_limit" && !area_diag
                if V_cell<v_l
                    candidate_reason="voltage_cutoff_low"
                elseif V_cell>v_h
                    candidate_reason="voltage_cutoff_high"
                elseif phase_type===PHASE_DISCHARGE && V_cell<=v_l
                    candidate_reason="phase_target_low"
                elseif phase_type===PHASE_CHARGE && V_cell>=v_h
                    candidate_reason="phase_target_high"
                end
            end
            terminal=t_candidate==t_end || candidate_reason!="time_limit"
            # Keep the accepted left-state load; the mechanical dt accounts
            # for elapsed accepted intervals, and its result is published at
            # this candidate endpoint. This is explicit staggered coupling.
            a=variables["accepted"]
            candidate_mech=a.mech
            candidate_stress=latest_macro_stress
            prepared_snapshot=nothing
            next_count=mechanical_step_count
            candidate_pending_dt=pending_dt
            if czm_active
                next_count += 1
                candidate_pending_dt += h_used
                if next_count % case.opt.czm.update_interval == 0 || terminal
                    candidate_mech=deepcopy(a.mech)
                    t_czm_ns = time_ns()
                    try
                        czm_result = update_czm_damage!(case, candidate_mech, step_variables, T_nodes_czm_current;
                            dt_seconds=candidate_pending_dt * case.param.scale.t0,
                            sample_time_seconds=a.t_global, czm_workspace=variables["czm_workspace"])
                        mechanical_dt_s=candidate_pending_dt * case.param.scale.t0
                        candidate_pending_dt=0.0
                        step_variables = czm_output_to_variables(czm_result, step_variables;
                            committed_states=candidate_mech.damage_states)
                        δ_max_n_vals = [s.δ_max_n for s in candidate_mech.damage_states]
                        step_variables["czm δ_max_n"] = [maximum(δ_max_n_vals)]
                        step_variables["czm δ_mean_n"] = [mean(δ_max_n_vals)]

                        if czm_snapshots !== nothing
                            prepared_snapshot = CZMSnapshot(czm_result,candidate_mech,case.mesh["czm"],case.param;
                                time_s=(t_candidate-t0) * case.param.scale.t0,
                                mechanical_dt_s=mechanical_dt_s,cycle=czm_cycle,phase=czm_phase,
                                method=case.opt.czm.iter_method,czm_model=case.opt.czm.model,
                                plasticity=case.opt.czm.j2_plasticity,czm_workspace=variables["czm_workspace"])
                        end
                    finally
                        timing_totals["czm"] += (time_ns() - t_czm_ns) * 1e-9
                    end
                    if macro_stress_active
                        candidate_stress = compute_macro_stress(
                            case, candidate_mech, step_variables, T_nodes_czm_current)
                    end
                end
            end

            step_variables["time"]=t_candidate
            record_history = case.opt.outputType == "auto" || t_candidate==RunTime[vt] || terminal
            if record_history
                candidate_stress===nothing || write_macro_stress!(step_variables,1,candidate_stress)
                if t_candidate==RunTime[vt]
                    vt=min(vt+1,length(RunTime))
                end
            end
            prepared_snapshot===nothing || sizehint!(czm_snapshots,length(czm_snapshots)+1)

            # adjust time incremental step dt
            if !first_step && t_candidate==t_start+dt && case.opt.dtType == "auto"
                if error_y < 0.5 * case.opt.dtThreshold
                    dt = min(dt * 2, dt_max)
                elseif error_y >= 1.5 * case.opt.dtThreshold
                    dt = dt_min
                elseif error_y > case.opt.dtThreshold
                    dt = max(dt / 2, dt_min)
                end
            end

            # Publish only after chemistry, mechanics, stress, history and snapshot preparation succeed.
            candidate=(y=copy(y_new),mech=candidate_mech,t_global=phase_time_origin+(t_candidate-t0)*case.param.scale.t0)
            prepared_K_old=copy(K_new);prepared_F_old=copy(F_new)
            prepared_T_nodes = if multi_spme_enabled
                get_thermal_dofs(candidate.y,case)
            else
                T_nodes_carry
            end
            if chemical_step_observer !== nothing
                y_before_ref=y_old[vc]
                y_after_ref=candidate.y[vc]
                K_old_ref=K_old; F_old_ref=F_old
            end
            if record_history
                Variable_update!(variables["history"],step_variables,v+1;debug_coupling=case.opt.debug_coupling)
                v += 1
            end
            variables["accepted"]=candidate
            pending_dt=candidate_pending_dt
            mechanical_step_count=next_count
            prepared_snapshot===nothing || push!(czm_snapshots,prepared_snapshot)
            latest_macro_stress=candidate_stress
            T_nodes_carry=prepared_T_nodes
            y_old=candidate.y
            K_old=prepared_K_old; F_old=prepared_F_old
            if chemical_step_observer !== nothing
                chem_accepted_idx += 1
                chemical_step_observer((;kind=:accepted_commit,h_used=h_used,dt_next=dt,theta=theta_used,
                    t_eval,t_start,t_end=t_candidate,t_load=t_start,
                    voltage=step_variables["cell voltage"]*case.param.scale.phi,
                    trial_idx=chem_trial_idx,accepted_idx=chem_accepted_idx,
                    M_new=M_new,K_old=K_old_ref,K_new=K_new,F_old=F_old_ref,F_new=F_new,
                    y_before=y_before_ref,y_after=y_after_ref))
            end
            t=t_candidate
            first_step=false
            termination_reason=candidate_reason
            terminal && break
        end

    end

    result = PostProcessing(case, variables["history"], v)
    # 汇总耗时统计（用于识别主要瓶颈）
    timing_total = timing_totals["spme"] + timing_totals["branch"] + timing_totals["thermal"] + timing_totals["czm"]

    result["timing SPMe solve total [s]"] = timing_totals["spme"]
    result["timing branch solver total [s]"] = timing_totals["branch"]
    result["timing thermal distributed total [s]"] = timing_totals["thermal"]
    result["timing CZM model total [s]"] = timing_totals["czm"]

    result["timing SPMe solve avg [ms]"] = 1000.0 * timing_totals["spme"] / timing_call_count
    result["timing branch solver avg [ms]"] = 1000.0 * timing_totals["branch"] / timing_call_count
    result["timing thermal distributed avg [ms]"] = 1000.0 * timing_totals["thermal"] / timing_call_count
    result["timing CZM model avg [ms]"] = 1000.0 * timing_totals["czm"] / timing_call_count
    result["timing CallModel calls"] = Float64(timing_call_count)

    if timing_total > 0
        result["timing SPMe solve ratio [%]"] = 100.0 * timing_totals["spme"] / timing_total
        result["timing branch solver ratio [%]"] = 100.0 * timing_totals["branch"] / timing_total
        result["timing thermal distributed ratio [%]"] = 100.0 * timing_totals["thermal"] / timing_total
        result["timing CZM model ratio [%]"] = 100.0 * timing_totals["czm"] / timing_total
    else
        result["timing SPMe solve ratio [%]"] = 0.0
        result["timing branch solver ratio [%]"] = 0.0
        result["timing thermal distributed ratio [%]"] = 0.0
        result["timing CZM model ratio [%]"] = 0.0
    end

    if case.opt.debug_coupling
        println("\n[Solve-Timing] CallModel 阶段累计耗时（用于优化定位）")
        println("  calls = $(timing_call_count)")
        @printf("  SPMe solve           : %.3f s (%.2f%%), avg %.3f ms/call\n",
            timing_totals["spme"], result["timing SPMe solve ratio [%]"], result["timing SPMe solve avg [ms]"])
        @printf("  Branch solver        : %.3f s (%.2f%%), avg %.3f ms/call\n",
            timing_totals["branch"], result["timing branch solver ratio [%]"], result["timing branch solver avg [ms]"])
        @printf("  Thermal distributed  : %.3f s (%.2f%%), avg %.3f ms/call\n",
            timing_totals["thermal"], result["timing thermal distributed ratio [%]"], result["timing thermal distributed avg [ms]"])
        @printf("  CZM model            : %.3f s (%.2f%%), avg %.3f ms/call\n",
            timing_totals["czm"], result["timing CZM model ratio [%]"], result["timing CZM model avg [ms]"])
    end

    # 添加截止信息到结果
    result["termination_reason"] = termination_reason
    result["first_cutoff_detected"] = first_cutoff_detected
    if first_cutoff_detected
        result["first_cutoff_time [s]"] = first_cutoff_time
        result["first_cutoff_element"] = first_cutoff_element
        result["first_cutoff_ocv [V]"] = first_cutoff_ocv
    end
    result["total_cutoff_count"] = total_cutoff_count

    # 附加热相关历史数据
    if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D"
        Tref = case.param_dim.scale.T_ref
        result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
        result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
    end
    print("finish the simulation\n")
    if return_final_state
        result["final_state"]=deepcopy(variables["accepted"])
    end
    return result
end

function RecordMatrix!(case::Case, M::SparseArrays.SparseMatrixCSC{Float64, Int64}, K::SparseArrays.SparseMatrixCSC{Float64, Int64})
    l_np= case.mesh["negative particle"].nlen
    l_pp= case.mesh["positive particle"].nlen
    case.param.NE.M_d = M[1:l_np, 1:l_np]
    case.param.NE.K_d = K[1:l_np, 1:l_np]
    case.param.PE.M_d = M[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
    case.param.PE.K_d = K[l_np+1:l_np+l_pp, l_np+1:l_np+l_pp]
    return case
end

function ErrorEstimation(case::Case, y_old::Array{Float64}, y_new::Array{Float64}, coeff::Float64)
    error_y = 0.0
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        error_y = norm(y_new - y_old) / norm(y_old) * coeff
    else
        v_c_np = case.index["negative particle lithium concentration"]
        v_c_pp = case.index["positive particle lithium concentration"]
        v_c_el = case.index["electrolyte lithium concentration"]
        v_phi_np = case.index["negative electrode potential"]
        v_phi_pp = case.index["positive electrode potential"]
        v_phi_el = case.index["electrolyte potential"]
        for i in [v_c_np, v_c_pp, v_c_el, v_phi_pp, v_phi_el]
            if norm(y_old[i])>0
                error_y = max(error_y, norm(y_new[i] - y_old[i]) / norm(y_old[i]) * coeff)
            end
        end

    end
    return error_y
end
