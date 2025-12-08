function Solve(case::Case)
    # 将调试输出统一写入 output 目录下的日志文件（可通过 opt.debug_to_file 关闭）
    enable_file_log = true
    if hasproperty(case.opt, :debug_to_file)
        try
            enable_file_log = Bool(getfield(case.opt, :debug_to_file))
        catch
            enable_file_log = true
        end
    end
    log_io = nothing
    log_file = ""
    old_out = stdout
    old_err = stderr
    if enable_file_log
        log_dir = normpath(joinpath(@__DIR__, "..", "output"))
        try
            isdir(log_dir) || mkpath(log_dir)
            # 使用 epoch 秒作为时间戳，避免依赖 Dates 标准库的 import
            timestamp = string(round(Int, time()))
            log_file = joinpath(log_dir, "debug_$(timestamp).log")
            log_io = open(log_file, "w")
            println(log_io, "===== JuBat Debug Log $(timestamp) =====")
            flush(log_io)
            redirect_stdout(log_io)
            redirect_stderr(log_io)
        catch err
            # 如果日志初始化失败，则回退为标准控制台输出
            try
                log_io !== nothing && close(log_io)
            catch
            end
            log_io = nothing
        end
    end
    result = nothing
    try
        # Electrochemical time still scaled by t0
        dt_min = case.opt.dt[1] / case.param.scale.t0 
        dt_max = case.opt.dt[2] / case.param.scale.t0 
        RunTime = case.opt.time / case.param.scale.t0 
        t0 = RunTime[1] 
        t_end = RunTime[end]
    
    # initialisation（根据模式选择初始化函数）
    multi_spme_enabled = (
        case.opt.model == "SPMe" &&
        hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme &&
        case.opt.thermalmodel == "distributed2D" &&
        haskey(case.mesh, "thermal2D")
    )
    
    # 检查模式并初始化（优先级：简化耦合 > 多SPMe > 标准）
    if hasproperty(case.opt, :simple_thermal_coupling) && 
       case.opt.simple_thermal_coupling &&
       case.opt.thermalmodel == "distributed2D"
        y0 = ModelInitialisation_SimpleCoupling(case)
        if hasproperty(case.opt, :debug_simple_coupling) && case.opt.debug_simple_coupling
            println("[Solve] 简化耦合模式：状态向量维度 = $(length(y0))")
        end
    elseif multi_spme_enabled
        y0 = ModelInitialisation_MultiSPMe(case)
        if hasproperty(case.opt, :debug_multi_spme) && case.opt.debug_multi_spme
            println("[Solve] 多SPMe模式：状态向量维度 = $(length(y0))")
            layout = case.multi_spme_layout
            println("  ne = $(layout["ne"]), n_chem = $(layout["n_chem"]), nT = $(layout["nT"])")
        end
    else
        y0 = ModelInitialisation(case)
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
    # 对于多SPMe模式，由于每步需要存储ne个单元的数据，限制更严格
    max_steps = 100000  # 默认最大步数
    if hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme
        max_steps = 50000  # 多SPMe模式限制更严格
    end
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
    # DEBUG: 检查初始状态向量（只在有问题时打印）
    nan_in_y0 = sum(.!isfinite.(y0))
    if nan_in_y0 > 0
        @warn "初始状态向量包含 $(nan_in_y0) 个 NaN/Inf，长度 $(length(y0))"
    end
    M_old, K_old, F_old, variables, y_phi= CallModel(case, y0, t, jacobi="update")
    # 打印初始信息
    V_init = variables["cell voltage"] * case.param.scale.phi
    println("\n[Solve] 初始化完成: V=$V_init V, t_end=$(t_end * case.param.scale.t0)")
    # Thermal-distributed init if enabled
    if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
        # 初始化热场 T 为环境初温
        nnode_th = case.mesh["thermal2D"].nlen
        T_nodes = fill(case.param.cell.T0, nnode_th)
        variables["T_nodes"] = T_nodes
        # 若启用 collector-seeded 逻辑，则优先尝试从 mesh 读取内置的 layer_weights；如无，则回退采样计算
        if hasproperty(case.opt, :collector_seeded) && case.opt.collector_seeded
            try
                fks_mesh = try
                    jellyroll_get_layer_weights(case.mesh["thermal2D"])
                catch
                    nothing
                end
                if fks_mesh !== nothing
                    variables["thermal2D layer_weights"] = fks_mesh
                else
                    fks = jellyroll_element_layer_weights(case.mesh["thermal2D"], case.param_dim; nsamples_per_dim=4, logic=:spiral)
                    variables["thermal2D layer_weights"] = fks
                end
            catch err
                @warn "Failed to set layer_weights: $err"
            end
        end
    # 计算初始热应力（仅二维分布热）
        variables = thermal_stress(case, variables)
    end
    # 持久化热场（跨 CallModel 迭代携带）
    T_nodes_carry = if haskey(variables, "T_nodes") && isa(variables["T_nodes"], Array{Float64}) && !isempty(variables["T_nodes"])
        variables["T_nodes"]
    elseif haskey(case.mesh, "thermal2D")
        fill(case.param.cell.T0, case.mesh["thermal2D"].nlen)
    else
        Float64[]
    end

    # 跟踪元素温度（用于调试/作图）
    track_elem_index = 0
    T_elem_hist, time_hist = Float64[], Float64[]
    if case.opt.thermal_enabled && haskey(case.mesh, "thermal2D")
        ne_track = size(case.mesh["thermal2D"].element, 1)
        if ne_track > 0
            idx_env = try parse(Int, get(ENV, "JUBAT_TRACK_ELEM", "")) catch; nothing end
            track_elem_index = Int(clamp(idx_env !== nothing ? idx_env : round(ne_track/2), 1, ne_track))
            nodes_e0 = case.mesh["thermal2D"].element[track_elem_index, :]
            Te0 = length(T_nodes_carry) == case.mesh["thermal2D"].nlen ? sum(T_nodes_carry[nodes_e0]) / length(nodes_e0) : case.param.cell.T0
            push!(T_elem_hist, Te0)
            push!(time_hist, t0 * case.param.scale.t0)
        end
    end

    dt_init = 1e-8
    vc = 1:size(M_old,1)
    y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
    y_old = vcat(y_c, y_phi)
    # 检查初始求解步骤
    nan_in_yold = sum(.!isfinite.(y_old))
    if nan_in_yold > 0 || (multi_spme_enabled && maximum(abs.(y_old)) > 100.0)
        @warn "初始求解步骤异常: NaN=$(nan_in_yold), 范围=[$(minimum(y_old)), $(maximum(y_old))]"
    end
    Variable_update!(variables_hist, variables, v)
    t += dt 
    if case.opt.jacobi == "constant"
        RecordMatrix!(case, M_old, K_old)    # record system matrix information 
    end

    print( "start to solve the problem \n")

    # run the model
    while t <= t_end
        # 1) 先用当前热场的均温影响动力学
        if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D" && !isempty(T_nodes_carry)
            # T_nodes_carry is dimensionless (T/T_ref) under Scheme B thermal scaling
            Tm = mean(T_nodes_carry)  # dimensionless
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
        if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
            nT = case.mesh["thermal2D"].nlen
            n_tot = size(M_new, 1)
            # 热自由度位于化学自由度之后（按 blockdiag 顺序追加），属于 y_c 中最后 nT 个条目
            if length(y_c) == n_tot
                T_nodes = y_c[(end - nT + 1):end]
                variables["T_nodes"] = T_nodes
                T_nodes_carry = T_nodes
                variables = thermal_stress(case, variables)
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
                # 同步跟踪元素温度（在提交此时间步后、时间推进前记录）
                if track_elem_index > 0 && haskey(variables, "T_nodes") && isa(variables["T_nodes"], Array{Float64}) && haskey(case.mesh, "thermal2D")
                    nodes_e = case.mesh["thermal2D"].element[track_elem_index, :]
                    Tn_now = variables["T_nodes"]
                    if length(Tn_now) == case.mesh["thermal2D"].nlen
                        Te_now = sum(Tn_now[nodes_e]) / length(nodes_e)
                        push!(T_elem_hist, Te_now)
                        push!(time_hist, t * case.param.scale.t0)
                    end
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
        if variables["cell voltage"] * case.param.scale.phi < case.param.cell.v_l || variables["cell voltage"] * case.param.scale.phi > case.param.cell.v_h
            break
        end
    end
    result = PostProcessing(case, variables_hist, v) 
     # 附加热相关历史数据
    try
        for key in ["thermal2D element current", "thermal2D eta_n_e", "thermal2D eta_p_e"]
            haskey(variables_hist, key) && (result[key] = variables_hist[key][:, 1:v])
        end
        if haskey(variables_hist, "heat_source_fields") && size(variables_hist["heat_source_fields"], 1) > 0
            result["heat_source_fields"] = variables_hist["heat_source_fields"][:, 1:v]
        end
        if haskey(variables_hist, "average temperature") && size(variables_hist["average temperature"], 2) >= v
            result["average temperature [K]"] = vec(variables_hist["average temperature"][1, 1:v])
        end
        if haskey(variables_hist, "total heat source") && size(variables_hist["total heat source"], 2) >= v
            result["total heat source [W]"] = vec(variables_hist["total heat source"][1, 1:v])
        end
        if haskey(variables_hist, "thermal2D temperature") && size(variables_hist["thermal2D temperature"], 2) >= v
            result["thermal2D temperature [K]"] = variables_hist["thermal2D temperature"][:, 1:v]
        end
        if case.opt.thermal_enabled && haskey(case.mesh, "thermal2D")
            if isa(T_nodes_carry, Array{Float64}) && length(T_nodes_carry) == case.mesh["thermal2D"].nlen
                Tref = case.param_dim.scale.T_ref
                result["thermal2D T_nodes [K]"] = T_nodes_carry .* Tref
                result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
            end
            if !isempty(T_elem_hist) && length(T_elem_hist) == length(time_hist)
                result["thermal2D tracked element index"] = track_elem_index
                result["thermal2D tracked element time [s]"] = time_hist
                result["thermal2D tracked element T [K]"] = T_elem_hist .* case.param_dim.scale.T_ref
            end
        end
    catch
        # non-fatal
    end
    print("finish the simulation\n") 
    errors = errors[1:v] 
    # 延后返回，以便执行 finally 中的日志恢复
    finally
        if log_io !== nothing
            # 恢复标准输出/错误到控制台
            redirect_stdout(old_out)
            redirect_stderr(old_err)
            try
                close(log_io)
            catch
            end
            println("Debug log saved to $(log_file)")
        end
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
    if isempty(case.multi_spme_layout)
        error("CallModel_MultiSPMe requires populated multi_spme_layout. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    nT = layout["nT"]
    # 检查输入状态向量
    nan_count = sum(.!isfinite.(yt))
    if nan_count > 0
        chem_range = 1:(ne * n_chem)
        thermal_range = (ne * n_chem + 1):(ne * n_chem + nT)
        nan_chem = sum(.!isfinite.(yt[chem_range]))
        nan_thermal = sum(.!isfinite.(yt[thermal_range]))
        @warn "CallModel_MultiSPMe 收到 NaN 状态向量" t=t total_nan=nan_count chem_nan=nan_chem thermal_nan=nan_thermal
    end
    mesh_th = case.mesh["thermal2D"]
    param = case.param
    
    # 初始化 variables（提前，以便后续缓存/读取中使用）
    variables = StandardVariables(case, 1)
    
    # 1) 解析状态向量
    # 提取热场
    T_nodes = MultiSPMe_get_thermal_dofs(yt, case)
    # 温度场检查
    if length(T_nodes) > 0
        nan_count_here = sum(.!isfinite.(T_nodes))
        abnormal_count = sum(abs.(T_nodes) .> 10.0)
        large_deviation_count = sum(abs.(T_nodes .- 1.0) .> 5.0)
        
        if nan_count_here > 0 || abnormal_count > 0 || large_deviation_count > 0
            T_min, T_max = extrema(T_nodes)
            @warn "温度场异常" t=t range=(T_min,T_max) nan=nan_count_here abnormal=abnormal_count large_dev=large_deviation_count
        end
    end
    # 提取每个单元的电化学状态
    yt_chem = Vector{Vector{Float64}}(undef, ne)
    for e in 1:ne
        yt_chem[e] = vec(MultiSPMe_extract_element_state(yt, e, case))
    end
    
    # 2) 计算元素面积和均温
    areas = if haskey(variables, "thermal2D element area")
        variables["thermal2D element area"]
    else
        A = zeros(Float64, ne)
        ngs = length(mesh_th.gs.detJ)
        @inbounds for g in 1:ngs
            e = mesh_th.gs.ele[g]
            A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
        end
        variables["thermal2D element area"] = A
        A
    end
    
    Te_prev = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        Te_prev[e] = sum(T_nodes[nds]) / length(nds)
    end
      # 检查温度场
    nan_count_nodes = sum(.!isfinite.(T_nodes))
    nan_count_elem = sum(.!isfinite.(Te_prev))
    if nan_count_nodes > 0 || nan_count_elem > 0
        @warn "温度场包含 NaN/Inf" T_nodes_nan=nan_count_nodes Te_prev_nan=nan_count_elem
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

    variables, I_e, Vc = solve_branch_currents_newton(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev)
    
    # 4) 并行求解每个单元的SPMe
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)
    
    # 可选：并行化（如果Julia配置了多线程）
    # Threads.@threads for e in 1:ne
    for e in 1:ne
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
    
    # 6) 计算逐单元热源
    q_elem = zeros(Float64, ne)
    eta_n_e = zeros(Float64, ne)
    eta_p_e = zeros(Float64, ne)
    dUdT_n_e = zeros(Float64, ne)
    dUdT_p_e = zeros(Float64, ne)
    
    L_th = case.param_dim.scale.L_th
    
    # 获取 layer_weights（如果存在）
    fks = haskey(variables, "thermal2D layer_weights") ? variables["thermal2D layer_weights"] : nothing
    if fks === nothing
        mesh_loc = mesh_th
        fks = try
            jellyroll_get_layer_weights(mesh_loc)
        catch
            nothing
        end
        if fks === nothing
            fks = jellyroll_element_layer_weights(mesh_loc, case.param_dim; nsamples_per_dim=4, logic=:spiral)
        end
        variables["thermal2D layer_weights"] = fks
    end
    
    for e in 1:ne
        vars_e = variables_elems[e]
        # 提取逐单元变量
        eta_n_e[e] = vars_e["negative electrode overpotential"][1]
        eta_p_e[e] = vars_e["positive electrode overpotential"][end]
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        dUdT_n_e[e] = param.NE.dUdT(cn_surf_e)[1]
        dUdT_p_e[e] = param.PE.dUdT(cp_surf_e)[1]
        # 当前单元的温度和电流
        T_e = Te_prev[e]
        I_e_local = I_e[e]
        # 界面电流密度（从单元 variables_e 获取）
        j_n_e = vars_e["negative electrode interfacial current density"]
        j_p_e = vars_e["positive electrode interfacial current density"]
        # 比表面积
        as_n = param.NE.as
        as_p = param.PE.as
        # 有效电导率
        sig_n_eff = param.NE.sig * param.NE.eps_s
        sig_p_eff = param.PE.sig * param.PE.eps_s
        kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps ^ param.NE.brugg
        kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps ^ param.PE.brugg
        kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps ^ param.SP.brugg
        # 分层计算热源
        # 负极层
        Q_rxn_NE = as_n * abs(j_n_e) * abs(eta_n_e[e])  # 反应热：a_s * j * η
        Q_rev_NE = as_n * j_n_e * T_e * dUdT_n_e[e]     # 可逆热：a_s * j * T * dU/dT
        Q_ohm_s_NE = I_e_local^2 / (3.0 * sig_n_eff)    # 固相欧姆热：I²/(3σ)
        Q_ohm_e_NE = I_e_local^2 / (3.0 * kappa_ne)     # 液相欧姆热：I²/(3κ)
        Q_NE = Q_rxn_NE + Q_rev_NE + Q_ohm_s_NE + Q_ohm_e_NE
        # 隔膜层（SP）- 仅液相欧姆热
        Q_SP = I_e_local^2 / kappa_sp  # 隔膜欧姆热：I²/κ（无1/3因子）
        # 正极层（PE）
        Q_rxn_PE = as_p * abs(j_p_e) * abs(eta_p_e[e])  # 反应热
        Q_rev_PE = as_p * j_p_e * T_e * dUdT_p_e[e]     # 可逆热
        Q_ohm_s_PE = I_e_local^2 / (3.0 * sig_p_eff)    # 固相欧姆热
        Q_ohm_e_PE = I_e_local^2 / (3.0 * kappa_pe)     # 液相欧姆热
        Q_PE = Q_rxn_PE + Q_rev_PE + Q_ohm_s_PE + Q_ohm_e_PE
        # 集流体层
        Q_PCC = I_e_local^2 / (3.0 * param.PCC.sig) 
        Q_NCC = I_e_local^2 / (3.0 * param.NCC.sig) 
        # 按层权重聚合
        q_elem[e] = fks[e,1]*Q_NE + fks[e,2]*Q_SP + fks[e,3]*Q_PE + fks[e,4]*Q_PCC + fks[e,5]*Q_NCC
    end
    # 写入逐单元变量到 variables（用于调试和后处理）
    variables["thermal2D element current"] = I_e
    variables["thermal2D eta_n_e"] = eta_n_e
    variables["thermal2D eta_p_e"] = eta_p_e
    variables["thermal2D dUdT_n_e"] = dUdT_n_e
    variables["thermal2D dUdT_p_e"] = dUdT_p_e
    
    # 无量纲化热源
    # 注意: q_elem 使用电化学尺度计算，需要转换为热传导尺度
    if hasproperty(case.opt, :units_thermal) && case.opt.units_thermal == "SI"
        variables["heat_source_fields"] = q_elem
        variables["heat_source_units_code"] = 1.0
    else
        # 尺度转换: 电化学尺度 → 傅里叶尺度
        q_ec_scale = case.param_dim.scale.I_typ * case.param_dim.scale.phi / case.param_dim.cell.volume
        q_fourier_scale = case.param_dim.scale.k_th * case.param_dim.scale.T_ref / case.param_dim.scale.L_th^2
        scale_conversion = q_ec_scale / q_fourier_scale
        
        q_ref = case.param_dim.scale.q_th  # 傅里叶尺度
        variables["heat_source_fields"] = (q_elem .* scale_conversion) ./ q_ref
        variables["heat_source_units_code"] = 0.0
    end
    
    # 7) 装配热学矩阵
    MT, KT, FT = ThermalDistributed2D(case, variables)
    t_ratio =  case.param_dim.scale.t0 / case.param_dim.scale.t_th
    # 检查参数（只在初始调用时）
    if t < 1e-6
        (t_ratio < 0.001 || t_ratio > 1000.0) && @warn "时间尺度比异常" t_ratio=t_ratio
        if haskey(variables, "heat_source_fields")
            q_elem = variables["heat_source_fields"]
            q_max = maximum(abs.(q_elem))
            q_max > 100.0 && @warn "无量纲热源过大" q_max=q_max q_range=extrema(q_elem)
        end
    end
    MT = MT .* t_ratio
    ThermalDistributed2D_BC(KT, FT, case, t)
    
    # 8) 全局拼装
    M = blockdiag(M_chem, sparse(MT))
    K = blockdiag(K_chem, sparse(KT))
    F = [F_chem; FT]
    
    # 9) 合并 variables（保留关键全局信息）
    # 电压取公共电压 Vc
    variables["cell voltage"] = Vc
    variables["time"] = t
    variables["temperature"] = mean(T_nodes)  # 平均温度
    variables["T_nodes"] = T_nodes
    
    # 可选：添加单元电压分布（用于诊断）
    variables["thermal2D element voltages"] = [variables_elems[e]["cell voltage"] for e in 1:ne]
    
    y_phi = Float64[]
    
    return M, K, F, variables, y_phi
end
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    # 检查简化耦合模式
    if hasproperty(case.opt, :simple_thermal_coupling) && 
       case.opt.simple_thermal_coupling &&
       case.opt.thermalmodel == "distributed2D"
        return CallModel_SimpleCoupling(case, yt, t; jacobi=jacobi)
    end
    
    # 判断是否启用多SPMe模式
    multi_spme_enabled = (
        case.opt.model == "SPMe" &&
        hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme &&
        case.opt.thermalmodel == "distributed2D" &&
        haskey(case.mesh, "thermal2D") &&
        !isempty(case.multi_spme_layout)
    )
    
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
        if haskey(case.mesh, "thermal2D")
            # 确保有 T_nodes（用于热源/材料系数等），若缺失则填充为初温
            if !haskey(variables, "T_nodes") || (isa(variables["T_nodes"], Array{Float64}) && length(variables["T_nodes"]) == 0)
                nT = case.mesh["thermal2D"].nlen
                variables["T_nodes"] = fill(case.param.cell.T0, nT)
            end
            # 面积缓存
            mesh_th = case.mesh["thermal2D"]
            if !haskey(variables, "thermal2D element area")
                ne_loc = size(mesh_th.element, 1)
                A_loc = zeros(Float64, ne_loc)
                ngs_loc = length(mesh_th.gs.detJ)
                @inbounds for g in 1:ngs_loc
                    e = mesh_th.gs.ele[g]
                    A_loc[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
                end
                variables["thermal2D element area"] = A_loc
            end
            # 计算元素均温，准备非线性分流求解
            areas = variables["thermal2D element area"]
            ne_loc = length(areas)
            T_nodes_loc = variables["T_nodes"]
            Te_prev = zeros(Float64, ne_loc)
            @inbounds for e in 1:ne_loc
                nds = mesh_th.element[e, :]
                Te_prev[e] = sum(T_nodes_loc[nds]) / length(nds)
            end
            # 总电流（无量纲，相对 I1C_total）
            I_total = 0.0
            if haskey(variables, "cell current")
                Ival = variables["cell current"]
                if isa(Ival, Float64)
                    I_total = Ival
                elseif isa(Ival, Array{Float64})
                    I_total = (ndims(Ival) == 1 ? (length(Ival) > 0 ? Ival[1] : 0.0) : (size(Ival,1) > 0 ? Ival[1,1] : 0.0))
                end
            end
            # 使用非线性分流求解器求每单元电流（不进行面积分流回退）
            try
                variables, _Ie, _Vc = solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, nothing)
            catch err
                # 不回退到面积分流：直接抛出异常以便暴露问题
                error("solve_branch_currents_newton failed in CallModel: $(err)")
            end
            # 更新热源（统一在 CallModel 内完成）
            try
                variables = heatQ_Source(case, variables, t, yt)
            catch err
                @warn "heatQ_Source failed in CallModel, continue with zero heat" err
                variables["heat_source_fields"] = zeros(Float64, ne_loc)
                variables["heat_source_units_code"] = 0.0
            end
            # 装配热学矩阵并施加边界条件
            MT, KT, FT = ThermalDistributed2D(case, variables)
            # 时间尺度匹配：主求解器以 t0 为时间标尺，热模块以 t_th 为标尺，
            t_ratio = case.param_dim.scale.t0 / case.param_dim.scale.t_th
            MT = MT .* t_ratio
            ThermalDistributed2D_BC(KT, FT, case, t)
            # 拼接到主系统
            M = blockdiag(M, sparse(MT))
            K = blockdiag(K, sparse(KT))
            F = [F; FT]
        end
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