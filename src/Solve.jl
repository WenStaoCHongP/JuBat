function Solve(case::Case; initial_state::Union{Dict{String,Any},Nothing}=nothing, return_final_state::Bool=false)
    if case.opt.model == "thermal"
            t0 = case.opt.time[1]/case.param.scale.t0
            t_end = case.opt.time[end]/case.param.scale.t0
            dt = case.opt.dt[1]/case.param.scale.t0
            vars = case.thermal_extras["thermal_variables"]
            update_fn = case.thermal_extras["thermal_update_fn"]
            record = case.thermal_extras["thermal_record"]

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
            T_nodes = vars["T_nodes"]

            times = collect(range(t0, step=dt, length=Int(cld(t_end - t0, dt)) + 1))
            T_hist = record ? zeros(Float64, nnode, length(times)) : zeros(Float64, 0, 0)
            record && (T_hist[:, 1] .= T_nodes)

            update_fn !== nothing && update_fn(t0, vars)
            if case.opt.thermalmodel == "ring2D_polar"
                mesh_data = case.thermal_extras["polar_mesh_data"]
                MT_old, KT_old, FT_old = ThermalPolar2D_Ring(case, vars, mesh_data)
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
                    mesh_data = case.thermal_extras["polar_mesh_data"]
                    MT_new, KT_new, FT_new = ThermalPolar2D_Ring(case, vars, mesh_data)
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
                vars["T_nodes"] = T_nodes

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
    
    # 判断是否启用多SPMe模式（当使用分布式2D热模型时自动启用）
    multi_spme_enabled = (case.opt.model == "SPMe" &&case.opt.thermalmodel == "distributed2D")
    
    # 允许从外部状态继续求解（用于循环阶段衔接）
    y0_input = initial_state === nothing ? nothing : get(initial_state, "y", nothing)
    T_nodes_input = initial_state === nothing ? nothing : get(initial_state, "T_nodes", nothing)

    # 检查模式并初始化（优先级：外部状态 > 多SPMe > 标准）
    if y0_input === nothing
        if multi_spme_enabled
            y0 = ModelInitialisation_MultiSPMe(case)
        else
            y0 = ModelInitialisation(case)
        end
    else
        y0 = vec(y0_input)
        if multi_spme_enabled
            ne = size(case.mesh["thermal2D"].element, 1)
            nT = case.mesh["thermal2D"].nlen
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Nel = case.mesh["electrolyte"].nlen
            n_chem = Nrn + Nrp + Nel
            expected_multi_len = ne * n_chem + nT

            if case.layout === nothing
                case.layout = MultiSPMeLayout(ne, n_chem, nT)
            end

            if length(y0) != expected_multi_len
                @warn "外部状态长度与多SPMe布局不匹配，回退到模型初始化" got=length(y0) expected=expected_multi_len
                y0 = ModelInitialisation_MultiSPMe(case)
            end
        elseif case.opt.thermalmodel == "distributed2D"
            nT = case.mesh["thermal2D"].nlen
            n_chem_base = case.mesh["negative particle"].nlen + case.mesh["positive particle"].nlen + case.mesh["electrolyte"].nlen
            if length(y0) == n_chem_base
                T_seed = (T_nodes_input !== nothing && length(T_nodes_input) == nT) ? vec(T_nodes_input) : fill(case.param.cell.T0, nT)
                y0 = [y0; T_seed]
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

    dt = deepcopy(dt_min)
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
    M_old, K_old, F_old, variables, y_phi= CallModel(case, y0, t, jacobi="update")
    # 打印初始信息
    V_init = variables["cell voltage"] * case.param.scale.phi
    println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")
    # Thermal-distributed init if enabled
    if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
        # 初始化热场，优先使用传入状态，其次使用状态向量中的热自由度。
        nnode_th = case.mesh["thermal2D"].nlen
        if T_nodes_input !== nothing && length(T_nodes_input) == nnode_th
            T_nodes = copy(vec(T_nodes_input))
        elseif length(y0) >= nnode_th
            T_nodes = copy(y0[(end - nnode_th + 1):end])
        else
            T_nodes = fill(case.param.cell.T0, nnode_th)
        end
        variables["T_nodes"] = T_nodes
    end
    # 持久化热场（跨 CallModel 迭代携带）
    T_nodes_carry = case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D" ? variables["T_nodes"] : Float64[]

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
    
    # run the model
    while t <= t_end
        # 1) 先用当前热场的均温影响动力学
        if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
            # T_nodes_carry is dimensionless (T/T_ref) under Scheme B thermal scaling
            Tm = thermal2D_volume_average_temperature(case.mesh["thermal2D"], T_nodes_carry)  # dimensionless
            case.param.cell.T0 = Tm  # SPMe expects dimensionless T
        end

        # 2) 电化学步
        M_new, K_new, F_new, variables, y_phi = CallModel(case, y_old, t, jacobi="update") 
        Mt = M_new - theta * K_new * dt 
        Kt = (1 - theta) * K_old * dt + M_new 
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt 
        y_c = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old[vc] + Ft) 
        y_new = vcat(y_c, y_phi)

        # 3) 分布式热：统一到 CallModel 的 M/K/F 内，不再在主循环单独做一步。
        # 仅在求解后提取温度自由度并回写到 variables 以用于记录/后处理。
        if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
            nT = case.mesh["thermal2D"].nlen
            n_tot = size(M_new, 1)
            # 热自由度位于化学自由度之后（按 blockdiag 顺序追加），属于 y_c 中最后 nT 个条目
            if length(y_c) == n_tot
                T_nodes = y_c[(end - nT + 1):end]
                variables["T_nodes"] = T_nodes
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
                    dt = deepcopy(dt_min)
                elseif error_y > case.opt.dtThreshold
                    dt = max(dt / 2, dt_min) 
                end
            elseif dt_temp_flag
                dt = deepcopy(dt_temp)
                dt_temp_flag = false    
            end
            if t + dt > RunTime[vt] && t < RunTime[vt]
                dt_temp = deepcopy(dt)
                dt = abs(RunTime[vt] - t) 
                dt_temp_flag = true
            end

            # update system information
            y_old = deepcopy(y_new)
            K_old = deepcopy(K_new)
            F_old = deepcopy(F_new)
            t += dt 
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
        if case.opt.thermalmodel == "distributed2D"
            for key in ["thermal2D element current", "thermal2D eta_n_e", "thermal2D eta_p_e", "thermal2D element soc_n", "thermal2D element soc_p",
                    "thermal2D q_rxn_ne", "thermal2D q_rev_ne", "thermal2D q_ohm_s_ne", "thermal2D q_ohm_e_ne",
                    "thermal2D q_sp", "thermal2D q_rxn_pe", "thermal2D q_rev_pe", "thermal2D q_ohm_s_pe", "thermal2D q_ohm_e_pe",
                    "thermal2D q_pcc", "thermal2D q_ncc",
                    "thermal2D element OCV", "thermal2D n_cutoff_elements", "thermal2D active_mask",
                    "thermal2D nearest_cutoff_element", "thermal2D nearest_cutoff_ocv", "thermal2D margin_to_cutoff"]
                result[key] = variables_hist[key][:, 1:v]
            end
            result["heat_source_fields"] = variables_hist["heat_source_fields"][:, 1:v]
            result["total heat source [W]"] = vec(variables_hist["total heat source"][1, 1:v])

            # 热源物理单位转换（无量纲 → 物理单位)
            q_scale = case.param_dim.scale.q
            result["thermal2D Q_rxn_NE [W/m3]"] = variables_hist["thermal2D q_rxn_ne"][:, 1:v] .* q_scale
            result["thermal2D Q_rev_NE [W/m3]"] = variables_hist["thermal2D q_rev_ne"][:, 1:v] .* q_scale
            result["thermal2D Q_ohm_s_NE [W/m3]"] = variables_hist["thermal2D q_ohm_s_ne"][:, 1:v] .* q_scale
            result["thermal2D Q_ohm_e_NE [W/m3]"] = variables_hist["thermal2D q_ohm_e_ne"][:, 1:v] .* q_scale
            result["thermal2D Q_SP [W/m3]"] = variables_hist["thermal2D q_sp"][:, 1:v] .* q_scale
            result["thermal2D Q_rxn_PE [W/m3]"] = variables_hist["thermal2D q_rxn_pe"][:, 1:v] .* q_scale
            result["thermal2D Q_rev_PE [W/m3]"] = variables_hist["thermal2D q_rev_pe"][:, 1:v] .* q_scale
            result["thermal2D Q_ohm_s_PE [W/m3]"] = variables_hist["thermal2D q_ohm_s_pe"][:, 1:v] .* q_scale
            result["thermal2D Q_ohm_e_PE [W/m3]"] = variables_hist["thermal2D q_ohm_e_pe"][:, 1:v] .* q_scale
            result["thermal2D Q_PCC [W/m3]"] = variables_hist["thermal2D q_pcc"][:, 1:v] .* q_scale
            result["thermal2D Q_NCC [W/m3]"] = variables_hist["thermal2D q_ncc"][:, 1:v] .* q_scale

            # 节点温度时间序列
            Tref = case.param_dim.scale.T_ref
            T_nodes_hist = variables_hist["T_nodes"][:, 1:v] .* Tref
            result["thermal2D temperature at nodes [K]"] = T_nodes_hist

            # 单元温度时间序列（节点平均）
            mesh_th = case.mesh["thermal2D"]
            ne = size(mesh_th.element, 1)
            n_t = size(T_nodes_hist, 2)
            T_elem_temp_hist = zeros(Float64, ne, n_t)
            for ti in 1:n_t
                T_nodes_t = variables_hist["T_nodes"][:, ti]
                T_elem_temp_hist[:, ti] = element_nodal_mean(mesh_th, T_nodes_t)
            end
            result["thermal2D temperature [K]"] = T_elem_temp_hist .* Tref
        end
        if case.opt.thermal_enabled
            if isa(T_nodes_carry, Array{Float64}) && length(T_nodes_carry) == case.mesh["thermal2D"].nlen
                Tref = case.param_dim.scale.T_ref
                result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
                result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
            end
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
"""
    CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)

多SPMe并行架构的 CallModel 实现。

此函数为每个热单元调用独立的 SPMe_element 求解器，然后将所有单元的矩阵和热学矩阵
全局装配成统一的系统。

# 状态向量结构
输入 yt:
    [yt_e[1]; yt_e[2]; ...; yt_e[ne]; T_nodes]
输出 M, K, F 对应的状态向量结构相同。

# 工作流程
1. 解析状态向量（提取每个单元的电化学状态和热场）
2. 计算元素均温和面积
3. 调用分流求解器获取逐单元电流 I_e
4. 并行调用 SPMe_element 求解每个单元的电化学响应
5. 计算逐单元热源（使用各单元的局部 η 和 dUdT）
6. 装配热学矩阵
7. 全局装配：blockdiag(M_e[1], ..., M_e[ne], MT)

# 返回
- M, K, F: 全局系统矩阵
- variables: 合并的变量字典（包含全局信息和逐单元信息）
- y_phi: 空向量（SPMe无电势自由度）
"""
function CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    # 验证前提条件
    if case.layout === nothing
        error("CallModel_MultiSPMe requires case.layout to be set. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.layout
    ne = layout.ne
    n_chem = layout.n_chem
    nT = layout.nT
    mesh_th = case.mesh["thermal2D"]
    param = case.param
    
    # 初始化 variables（提前，以便后续缓存/读取中使用）
    variables = StandardVariables(case, 1)
    
    # 1) 解析状态向量
    # 提取热场
    T_nodes = get_thermal_dofs(yt, case.layout)
    # 提取每个单元的电化学状态
    yt_chem = Vector{Vector{Float64}}(undef, ne)
    for e in 1:ne
        yt_chem[e] = vec(extract_element_state(yt, e, case.layout))
    end
    
    # 2) 计算元素面积和均温
    A = zeros(Float64, ne)
    ngs = length(mesh_th.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh_th.gs.ele[g]
        A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    variables["thermal2D element area"] = A
    areas = A
    
    Te_prev = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        Te_prev[e] = sum(T_nodes[nds]) / length(nds)
    end
    # 3) 分流求解（获取 I_e）
    # 使用与 SPMe 变量一致的电流无量纲尺度（param.scale.I_typ）
    I_total = case.opt.Current(t * case.param.scale.t0) / case.param.scale.I_typ
    
    # 更新 variables（用于传递给分流求解器）
    variables["cell current"] = I_total
    variables["T_nodes"] = T_nodes
    variables["thermal2D element area"] = areas
    
    # 使用缓存的 I_e 作为初值
    I_e_prev = hasproperty(case, :I_e_cache) ? case.I_e_cache : nothing
    
    # 分流求解需要代表性的全局状态（用于计算 prefactor）。
    # 这里使用所有单元的平均状态，并先生成与之匹配的电化学变量，确保标度一致。
    yt_representative = mean(yt_chem)
    T_rep = mean(Te_prev)
    vars_rep = SPMe_variables(case, yt_representative, t; I_app=I_total, T_e=T_rep)
    for (k, v) in vars_rep
        variables[k] = v
    end
    # 保留热相关场
    variables["T_nodes"] = T_nodes
    variables["thermal2D element area"] = areas

    # 获取CZM失效单元列表（仅在启用CZM时）
    if case.opt.czm_enabled
        deactivated_elements = Int64[]
        try
            deactivated_elements = convert(Vector{Int64}, variables["deactivated_elements"])
        catch
            deactivated_elements = Int64[]
        end
    else
        deactivated_elements = Int64[]
    end
    
    variables, I_e, Vc = solve_branch_currents_newton(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev; deactivated_elements=deactivated_elements)
    
    # 4) 并行求解每个单元的SPMe
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)
    
    # 可选：并行化（如果Julia配置了多线程）
    Threads.@threads for e in 1:ne
    #for e in 1:ne
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e;I_e = I_e[e],T_e = Te_prev[e],jacobi = jacobi)
        M_elems[e] = sparse(M_e)
        K_elems[e] = sparse(K_e)
        F_elems[e] = vec(F_e)
        variables_elems[e] = vars_e
    end
    
    # 5) 装配电化学全局矩阵
    M_chem = blockdiag(M_elems...)
    K_chem = blockdiag(K_elems...)
    F_chem = vcat(F_elems...)
    
    # 6) 计算逐单元热源（调用 ThermalDistributed.jl 中的统一函数)
    if case.opt.czm_enabled == true
        variables = compute_heat_sources_with_czm(case, variables, variables_elems, I_e, Te_prev, areas, czm_mesh, mesh_th)
    else
        variables = compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas; per_element_spme=true)
    end

    # 保存辅助变量（用于调试）
    for e in 1:ne
        vars_e = variables_elems[e]
        variables["thermal2D eta_n_e"][e] = vars_e["negative electrode overpotential"][1]
        variables["thermal2D eta_p_e"][e] = vars_e["positive electrode overpotential"][end]
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        variables["thermal2D dUdT_n_e"][e] = param.NE.dUdT(cn_surf_e)[1]
        variables["thermal2D dUdT_p_e"][e] = param.PE.dUdT(cp_surf_e)[1]
        csn_data = vars_e["negative particle lithium concentration"]
        csp_data = vars_e["positive particle lithium concentration"]
        variables["thermal2D element soc_n"][e] = mean(vec(csn_data))
        variables["thermal2D element soc_p"][e] = mean(vec(csp_data))
    end
    variables["thermal2D element current"] = I_e
    
    # 7) 装配热学矩阵
    MT, KT, FT = ThermalDistributed2D(case, variables)
    # 统一时间尺度：电化学和热模型使用相同时间尺度 t0，时间比为 1
    t_ratio = 1.0
    MT = MT .* t_ratio
    KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)
    
    # 8) 全局拼装
    M = blockdiag(M_chem, sparse(MT))
    K = blockdiag(K_chem, sparse(KT))
    F = [F_chem; FT]
    
    # 9) 合并 variables（保留关键全局信息）
    # 电压取公共电压 Vc
    variables["cell voltage"] = Vc
    variables["time"] = t
    variables["temperature"] = thermal2D_volume_average_temperature(case.mesh["thermal2D"], T_nodes)  # 体积平均温度
    variables["T_nodes"] = T_nodes

    # 可选：添加单元电压分布（用于诊断）
    variables["thermal2D element voltages"] = [variables_elems[e]["cell voltage"] for e in 1:ne]
    
    y_phi = Float64[]
    
    return M, K, F, variables, y_phi
end
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    # 判断是否应该启用多SPMe模式
    # 当使用分布式2D热模型时自动启用多SPMe
    should_use_multi_spme = (case.opt.model == "SPMe" && case.opt.thermalmodel == "distributed2D")
    
    # 如果应该使用多SPMe但布局为空，尝试从状态向量长度推断
    if should_use_multi_spme && case.layout === nothing
        # 计算期望的布局
        ne = size(case.mesh["thermal2D"].element, 1)
        nT = case.mesh["thermal2D"].nlen
        Nrn = case.mesh["negative particle"].nlen
        Nrp = case.mesh["positive particle"].nlen
        Nel = case.mesh["electrolyte"].nlen
        n_chem = Nrn + Nrp + Nel
        expected_multi_len = ne * n_chem + nT

        # 检查传入的状态向量长度
        if length(yt) == expected_multi_len
            # 状态向量是多SPMe格式，需要初始化布局
            @warn "CallModel: case.layout 为 nothing 但状态向量长度匹配多SPMe格式，自动初始化布局"
            case.layout = MultiSPMeLayout(ne, n_chem, nT)
        end
    end
    
    # 最终判断
    multi_spme_enabled = should_use_multi_spme && case.layout !== nothing
    
    if multi_spme_enabled
        return CallModel_MultiSPMe(case, yt, t, jacobi=jacobi)
    end
    
    # 原有逻辑（单SPMe模式）
    if case.opt.model == "SPM"
        M, K, F, variables = SPM(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "SPMe"
        M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "P2D"
        M, K, F, variables, y_phi = P2D(case, yt, t, jacobi=jacobi)
    elseif case.opt.model == "sP2D"
        M, K, F, variables, y_phi = sP2D(case, yt, t, jacobi=jacobi)
    else
        error( "Error: $(case.opt.model) model has not been implemented!\n ")
    end
    if case.opt.thermalmodel == "lumped"
        MT, FT = ThermalLumped(case, variables)        
        M = blockdiag(M, sparse(MT))
        K = blockdiag(K, sparse(zeros(1,1)))
        F = [F; FT]
    elseif case.opt.thermalmodel == "distributed2D"
        # 与 lumped 一致：在 CallModel 内根据最新电化学变量更新热源并装配热学 M/K/F，随后拼接。
        # 使用状态向量的热自由度作为当前温度场
        nT = case.mesh["thermal2D"].nlen
        variables["T_nodes"] = yt[(end - nT + 1):end]
        # 面积
        mesh_th = case.mesh["thermal2D"]
        ne_loc = size(mesh_th.element, 1)
        A_loc = zeros(Float64, ne_loc)
        ngs_loc = length(mesh_th.gs.detJ)
        @inbounds for g in 1:ngs_loc
            e = mesh_th.gs.ele[g]
            A_loc[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
        end
        variables["thermal2D element area"] = A_loc
        # 计算元素均温，准备非线性分流求解
        areas = A_loc
        T_nodes_loc = variables["T_nodes"]
        Te_prev = zeros(Float64, ne_loc)
        @inbounds for e in 1:ne_loc
            nds = mesh_th.element[e, :]
            Te_prev[e] = sum(T_nodes_loc[nds]) / length(nds)
        end
        # 总电流（无量纲，相对 I1C_total）
        Ival = variables["cell current"]
        I_total = isa(Ival, Float64) ? Ival : (ndims(Ival) == 1 ? Ival[1] : Ival[1,1])
        # 获取CZM失效单元列表（仅在启用CZM时）
        if case.opt.czm_enabled
            deactivated_elements_loc = Int64[]
            try
                deactivated_elements_loc = convert(Vector{Int64}, variables["deactivated_elements"])
            catch
                deactivated_elements_loc = Int64[]
            end
        else
            deactivated_elements_loc = Int64[]
        end
        # 使用非线性分流求解器求每单元电流（不进行面积分流回退）
        variables, _Ie, _Vc = solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, nothing; deactivated_elements=deactivated_elements_loc)
        # 更新热源（调用 ThermalDistributed.jl 中的统一函数）
        variables = compute_heat_sources(case, variables, nothing, _Ie, Te_prev, areas; per_element_spme=false)
        # 装配热学矩阵并施加边界条件
        MT, KT, FT = ThermalDistributed2D(case, variables)
        # 统一时间尺度：电化学和热模型使用相同时间尺度 t0，时间比为 1
        t_ratio = 1.0
        MT = MT .* t_ratio
        KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)
        # 拼接到主系统
        M = blockdiag(M, sparse(MT))
        K = blockdiag(K, sparse(KT))
        F = [F; FT]
    end
    return M, K, F, variables, y_phi
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
