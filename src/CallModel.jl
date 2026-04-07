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
    variables["thermal2D temperature at nodes"] = T_nodes
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
    variables["thermal2D temperature at nodes"] = T_nodes
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
    
    t_branch_ns = time_ns()
    variables, I_e, Vc = solve_branch_currents(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev; deactivated_elements=deactivated_elements)
    t_branch_s = (time_ns() - t_branch_ns) * 1e-9
    
    # 4) 并行求解每个单元的SPMe
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)
    
    # 可选：并行化（如果Julia配置了多线程）
    t_spme_ns = time_ns()
    Threads.@threads for e in 1:ne
    #for e in 1:ne
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e;I_e = I_e[e],T_e = Te_prev[e],jacobi = jacobi)
        M_elems[e] = sparse(M_e)
        K_elems[e] = sparse(K_e)
        F_elems[e] = vec(F_e)
        variables_elems[e] = vars_e
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
    
    # 5) 装配电化学全局矩阵
    M_chem = blockdiag(M_elems...)
    K_chem = blockdiag(K_elems...)
    F_chem = vcat(F_elems...)
    
    # 6) 计算逐单元热源（调用 ThermalDistributed.jl 中的统一函数)
    t_thermal_ns = time_ns()
    t_czm_model_s = 0.0
    if case.opt.czm_enabled == true
        t_czm_ns = time_ns()
        variables = compute_heat_sources_with_czm(case, variables, variables_elems, I_e, Te_prev, areas, case.czm_mesh, mesh_th)
        t_czm_model_s = (time_ns() - t_czm_ns) * 1e-9
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
    t_thermal_s = (time_ns() - t_thermal_ns) * 1e-9
    
    # 8) 全局拼装
    M = blockdiag(M_chem, sparse(MT))
    K = blockdiag(K_chem, sparse(KT))
    F = [F_chem; FT]
    
    # 9) 合并 variables（保留关键全局信息）
    # 电压取公共电压 Vc
    variables["cell voltage"] = Vc
    variables["time"] = t
    variables["temperature"] = thermal2D_volume_average_temperature(case.mesh["thermal2D"], T_nodes)  # 体积平均温度
    variables["thermal2D temperature at nodes"] = T_nodes

    # 可选：添加单元电压分布（用于诊断）
    variables["thermal2D element voltages"] = [variables_elems[e]["cell voltage"] for e in 1:ne]
    variables["timing branch solver [s]"] = t_branch_s
    variables["timing spme solve [s]"] = t_spme_s
    variables["timing thermal distributed [s]"] = t_thermal_s
    variables["timing czm model [s]"] = t_czm_model_s
    
    y_phi = Float64[]
    
    return M, K, F, variables, y_phi
end
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    # 多SPMe模式：由 per_element_spme 控制，直接委托给 CallModel_MultiSPMe
    if case.opt.per_element_spme
        return CallModel_MultiSPMe(case, yt, t, jacobi=jacobi)
    end

    # 单模型模式：与 main 分支逻辑一致
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
    end
    return M, K, F, variables, y_phi
end