function Solve(case::Case;initial_state::Union{Dict{String,Any},Nothing}=nothing,return_final_state::Bool=false,thermal_variables::Union{Dict{String,Any},Nothing}=nothing,thermal_update_fn::Union{Function,Nothing}=nothing,thermal_record::Bool=false,polar_mesh_data::Any=nothing)
    if case.opt.model == "thermal"
            if thermal_variables === nothing
                error("Solve: model==\"thermal\" requires thermal_variables keyword argument")
            end
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
            T_hist = record ? zeros(Float64, nnode, length(times)) : zeros(Float64, 0, 0)
            record && (T_hist[:, 1] .= T_nodes)

            update_fn !== nothing && update_fn(t0, vars)
            if case.opt.thermalmodel == "ring2D_polar"
                MT_old, KT_old, FT_old = ThermalPolar2D_Ring(case, vars, polar_mesh_data)
            elseif case.opt.thermalmodel == "ring2D"
                MT_old, KT_old, FT_old = ThermalDistributed2D_Ring(case, vars)
                KT_old, FT_old = ThermalRing2D_BC(KT_old, FT_old, case, vars["thermal2D outer_nodes"], t0)
            else
                MT_old, KT_old, FT_old = ThermalDistributed2D(case, vars)
                KT_old, FT_old = ThermalDistributed2D_BC(KT_old, FT_old, case, t0)
            end

            for step in 1:(length(times) - 1)
                t = times[step + 1]
                update_fn !== nothing && update_fn(t, vars)

                if case.opt.thermalmodel == "ring2D_polar"
                    MT_new, KT_new, FT_new = ThermalPolar2D_Ring(case, vars, polar_mesh_data)
                elseif case.opt.thermalmodel == "ring2D"
                    MT_new, KT_new, FT_new = ThermalDistributed2D_Ring(case, vars)
                    KT_new, FT_new = ThermalRing2D_BC(KT_new, FT_new, case, vars["thermal2D outer_nodes"], t)
                else
                    MT_new, KT_new, FT_new = ThermalDistributed2D(case, vars)
                    KT_new, FT_new = ThermalDistributed2D_BC(KT_new, FT_new, case, t)
                end

                A = MT_new - theta * dt * KT_new
                rhs = (MT_new + (1.0 - theta) * dt * KT_old) * T_nodes +
                      dt * (theta * FT_new + (1.0 - theta) * FT_old)
                T_nodes = A \ rhs
                vars["thermal2D temperature at nodes"] = T_nodes

                record && (T_hist[:, step + 1] .= T_nodes)

                MT_old = MT_new
                KT_old = KT_new
                FT_old = FT_new
            end

            result = (time = times, T_nodes = T_nodes, T_hist = T_hist)
            return result
        end
        # Electrochemical time still scaled by t0
        dt_min = case.opt.dt[1] / case.param.scale.t0 
        dt_max = case.opt.dt[2] / case.param.scale.t0 
        RunTime = case.opt.time / case.param.scale.t0 
        t0 = RunTime[1] 
        t_end = RunTime[end]
    
    # 判断是否启用多SPMe模式：与 CallModel 保持一致，由 per_element_spme 控制
    multi_spme_enabled = case.opt.per_element_spme

    # 允许从外部状态继续求解（用于循环阶段衔接）
    y0_input = initial_state

    # 检查模式并初始化（优先级：外部状态 > 多SPMe > 标准）
    if y0_input === nothing
        if multi_spme_enabled
            y0 = ModelInitialisation_MultiSPMe(case)
        else
            y0 = ModelInitialisation(case)
        end
    else
        # initial_state 可能是向量或 Dict（来自循环求解器）
        if isa(y0_input, Dict)
            y_from_state = get(y0_input, "y", nothing)
            if y_from_state !== nothing
                y0 = vec(y_from_state)
            else
                # Dict 中没有 y，回退到模型初始化
                if multi_spme_enabled
                    y0 = ModelInitialisation_MultiSPMe(case)
                else
                    y0 = ModelInitialisation(case)
                end
            end
        else
            y0 = vec(y0_input)
        end
        if multi_spme_enabled
            ne = size(case.mesh["thermal2D"].element, 1)
            nT = case.mesh["thermal2D"].nlen
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Nel = case.mesh["electrolyte"].nlen
            n_chem = Nrn + Nrp + Nel
            expected_multi_len = ne * n_chem + nT

            if case.layout === nothing
                case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
            end

            if length(y0) != expected_multi_len
                @warn "外部状态长度与多SPMe布局不匹配，回退到模型初始化" got=length(y0) expected=expected_multi_len
                y0 = ModelInitialisation_MultiSPMe(case)
            end
        end
    end
    if case.opt.solveType == "Crank-Nicolson"
        theta = 0.5 
    elseif case.opt.solveType == "forward"
        theta = 0 
    elseif case.opt.solveType == "backward"
        theta = 1 
    else
        error( "Error: $(opt.solve_type) difference scheme has not been implemented!\n ") 
    end

    dt = dt_min
    dt_temp = 0
    dt_temp_flag = false
    # 计算预期时间步数，添加安全限制避免内存溢出
    num_estimated = round(Int64, (t_end - t0)/dt * 1.5)
    # 限制最大预分配步数（避免内存溢出）
    max_steps = multi_spme_enabled ? 50000 : 100000
    num = min(num_estimated, max_steps)
    
    if num_estimated > max_steps
        @warn "预期时间步数 $(num_estimated) 超过最大限制 $(max_steps)，将使用动态扩展策略"
        @warn "这可能导致性能下降。建议增大时间步长 dt_min = $(dt*case.param.scale.t0) 秒"
    end
    
    variables_hist = StandardVariables(case, num)
    errors = zeros(num, 1)

    t = t0
    vt = 2  
    v = 1
    timing_totals = Dict{String,Float64}(
        "spme" => 0.0,
        "branch" => 0.0,
        "thermal" => 0.0,
        "czm" => 0.0,
    )
    timing_call_count = 0

    function accumulate_callmodel_timing!(totals::Dict{String,Float64}, vars::Dict{String, Union{Array{Float64},Float64}})
        totals["spme"] += get(vars, "timing spme solve [s]", 0.0)
        totals["branch"] += get(vars, "timing branch solver [s]", 0.0)
        totals["thermal"] += get(vars, "timing thermal distributed [s]", 0.0)
        totals["czm"] += get(vars, "timing czm model [s]", 0.0)
        return nothing
    end

    M_old, K_old, F_old, variables, y_phi= CallModel(case, y0, t, jacobi="update")
    timing_call_count += 1
    accumulate_callmodel_timing!(timing_totals, variables)
    # 打印初始信息
    V_init = variables["cell voltage"] * case.param.scale.phi
    println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")
    # 持久化热场（跨 CallModel 迭代携带）
    T_nodes_carry = get(variables, "thermal2D temperature at nodes", Float64[])

    dt_init = 1e-8
    vc = 1:size(M_old,1)
    y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
    y_old = vcat(y_c, y_phi)
    Variable_update!(variables_hist, variables, v)
    t += dt 
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

    # CZM 损伤演化状态
    czm_active = case.opt.czm_enabled && case.czm_mesh !== nothing
    u_czm_prev = nothing       # 上一步 CZM 位移场
    czm_step_count = 0         # CZM 更新步计数器
    
    # run the model
    while t <= t_end
        # 电化学步
        M_new, K_new, F_new, variables, y_phi = CallModel(case, y_old, t, jacobi="update") 
        timing_call_count += 1
        accumulate_callmodel_timing!(timing_totals, variables)
        Mt = M_new - theta * K_new * dt 
        Kt = (1 - theta) * K_old * dt + M_new 
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt 
        y_c = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old[vc] + Ft) 
        y_new = vcat(y_c, y_phi)

        # multi-SPMe：求解后提取温度自由度，用于记录/后处理
        if multi_spme_enabled
            nT = case.layout.nT
            if length(y_c) >= nT
                T_nodes = y_c[(end - nT + 1):end]
                variables["thermal2D temperature at nodes"] = T_nodes
                T_nodes_carry = T_nodes
            end
        end
        error_y = ErrorEstimation(case, y_old, y_new, dt_min/dt) 
        errors[v] = error_y
        if error_y > 2 * case.opt.dtThreshold && case.opt.dtType == "auto" && dt >= dt_min * 4
            # reduce dt to dt/2 and recaculate y_new
            dt = dt  /2
            t -= dt
        else 
            # record the results
            if case.opt.outputType == "auto" || abs(t - RunTime[vt]) < 1e-7
                v = v + 1 
                Variable_update!(variables_hist, variables, v) 
                if abs(t - RunTime[vt]) < 1e-7
                    vt = min(vt + 1, length(RunTime)) 
                end
            end
            
            # adjust time incremental step dt
            if  case.opt.dtType == "auto" && dt_temp_flag == false
                if error_y < 0.5 * case.opt.dtThreshold
                    dt = min(dt * 2, dt_max) 
                elseif error_y >= 1.5 * case.opt.dtThreshold
                    dt = dt_min
                elseif error_y > case.opt.dtThreshold
                    dt = max(dt / 2, dt_min) 
                end
            elseif dt_temp_flag
                dt = dt_temp
                dt_temp_flag = false
            end
            if t + dt > RunTime[vt] && t < RunTime[vt]
                dt_temp = dt
                dt = abs(RunTime[vt] - t) 
                dt_temp_flag = true
            end

            # update system information
            y_old = copy(y_new)
            K_old = copy(K_new)
            F_old = copy(F_new)
            t += dt

            # CZM 损伤演化（按间隔更新）
            if czm_active
                czm_step_count += 1
                if czm_step_count % case.opt.czm_update_interval == 0
                    t_czm_ns = time_ns()
                    try
                        u_czm_new, czm_converged = update_czm_damage!(
                            case.czm_mesh, case.param.cohesive,
                            case, variables, T_nodes_carry, u_czm_prev)
                        if czm_converged
                            u_czm_prev = u_czm_new
                        end
                    catch e
                        @debug "CZM damage update failed at step $czm_step_count: $e"
                    end
                    timing_totals["czm"] += (time_ns() - t_czm_ns) * 1e-9
                end
            end
        end
        
        # ====================================================================
        # 精细化截止电压检测（单元级别）
        # ====================================================================
        V_cell = variables["cell voltage"] * case.param.scale.phi
        v_l = case.param.cell.v_l
        v_h = case.param.cell.v_h
        t_phys = t * case.param.scale.t0  # 物理时间 (s)
        
        # 检查是否有单元截止
        if multi_spme_enabled
            n_cutoff = Int(variables["thermal2D n_cutoff_elements"])
            
            # 记录首个截止单元信息
            if n_cutoff > 0 && !first_cutoff_detected
                first_cutoff_detected = true
                first_cutoff_time = t_phys
                
                # 获取截止单元详细信息
                first_cutoff_element = Int(variables["thermal2D cutoff_elements"][1])
                first_cutoff_ocv = variables["thermal2D cutoff_ocv"][1]
            end
            
            total_cutoff_count = n_cutoff
            
            # 获取总单元数
            ne_total = size(case.mesh["thermal2D"].element, 1)
            
            # 检查是否所有单元都截止
            if ne_total > 0 && n_cutoff >= ne_total
                termination_reason = "all_elements_cutoff"
                break
            end
        end
        
        # 整体电压截止检测（备用）
        if V_cell < v_l
            termination_reason = "voltage_cutoff_low"
            break
        elseif V_cell > v_h
            termination_reason = "voltage_cutoff_high"
            break
        end
    end
    
    # 记录终止原因和截止信息
    if t >= t_end
        termination_reason = "time_limit"
    end
    result = PostProcessing(case, variables_hist, v) 
    # 汇总耗时统计（用于识别主要瓶颈）
    timing_total = timing_totals["spme"] + timing_totals["branch"] + timing_totals["thermal"] + timing_totals["czm"]
    call_count_safe = max(timing_call_count, 1)

    result["timing SPMe solve total [s]"] = timing_totals["spme"]
    result["timing branch solver total [s]"] = timing_totals["branch"]
    result["timing thermal distributed total [s]"] = timing_totals["thermal"]
    result["timing CZM model total [s]"] = timing_totals["czm"]

    result["timing SPMe solve avg [ms]"] = 1000.0 * timing_totals["spme"] / call_count_safe
    result["timing branch solver avg [ms]"] = 1000.0 * timing_totals["branch"] / call_count_safe
    result["timing thermal distributed avg [ms]"] = 1000.0 * timing_totals["thermal"] / call_count_safe
    result["timing CZM model avg [ms]"] = 1000.0 * timing_totals["czm"] / call_count_safe
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
    try
        if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D" && !isempty(T_nodes_carry)
            Tref = case.param_dim.scale.T_ref
            result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
            result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
        end
    catch
        # non-fatal
    end
    print("finish the simulation\n") 
    errors = errors[1:v] 
    if return_final_state
        V_final = variables["cell voltage"] * case.param.scale.phi
        t_final = max(0.0, t * case.param.scale.t0)
        result["final_state"] = Dict{String, Any}(
            "y" => copy(y_old),
            "T_nodes" => copy(T_nodes_carry),
            "V" => V_final,
            "t_global" => t_final
        )
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
