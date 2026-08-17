"""
	map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, ne_thermal::Int) -> Vector{Float64}

按 spec v2 §5.2 把 cohesive 单元损伤 max 归约到粗热单元粒度。
对每个粗热单元 e_thermal，扫描所有 cohesive_to_thermal[e_coh] == e_thermal 的 cohesive 单元，
取 damage_states[e_coh].D 的最大值；无 cohesive 覆盖则取 0。
"""
function map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, ne_thermal::Int)
	D_per_thermal = zeros(ne_thermal)
	for e_coh in 1:czm_mesh.n_cohesive
		e_thermal = czm_mesh.cohesive_to_thermal[e_coh]
		if 1 <= e_thermal <= ne_thermal
			D = czm_mesh.damage_states[e_coh].D
			if D > D_per_thermal[e_thermal]
				D_per_thermal[e_thermal] = D
			end
		end
	end
	return D_per_thermal
end

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
    # 提取每个单元的电化学状态（使用 @views 消除拷贝）
    yt_chem = Vector{SubArray{Float64,1}}(undef, ne)
    for e in 1:ne
        offset = (e - 1) * layout.n_chem
        yt_chem[e] = @view yt[(offset + 1):(offset + layout.n_chem)]
    end
    
    # 2) 使用缓存的元素面积（网格不变量）
    areas = layout.areas
    variables["thermal2D element area"] = areas
    
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
    # 当CZM单元损伤D >= 0.95时，对应热单元电流归零，总电流重分配
    if case.opt.czm_enabled && case.czm_mesh !== nothing
        fractured_czm = get_fractured_elements(case.czm_mesh)
        deactivated_elements = Int64[]
        geom = case.geometry
        for e in 1:ne
            for czm_idx in get(geom.czm_element_map, e, Int64[])
                if czm_idx in fractured_czm
                    push!(deactivated_elements, e)
                    break
                end
            end
        end
    else
        deactivated_elements = Int64[]
    end

    # 计算渐进式面积损失的 D 映射（仅在启用时）
    D_elem_area_loss = nothing
    if case.opt.czm_area_loss_enabled && case.czm_mesh !== nothing && case.czm_mesh.cohesive_to_thermal !== nothing
        D_elem_area_loss = map_czm_damage_to_thermal(case.czm_mesh, ne)
        # 调试：输出 D 映射统计
        D_above = filter(d -> d > case.opt.czm_area_loss_threshold, D_elem_area_loss)
        if !isempty(D_above) && case.opt.debug_coupling
            t_phys = round(t * case.param.scale.t0, digits=1)
            println("  [AreaLoss] t=$(t_phys)s | D_max=$(round(maximum(D_elem_area_loss), digits=4)) | 超阈值单元=$(length(D_above))/$(ne) | threshold=$(case.opt.czm_area_loss_threshold)")
        end
    end

    t_branch_ns = time_ns()
    variables, I_e, Vc = solve_branch_currents(case, variables, yt_representative, t, I_total, areas, Te_prev, nothing; deactivated_elements=deactivated_elements, D_elem=D_elem_area_loss)
    t_branch_s = (time_ns() - t_branch_ns) * 1e-9
    
    # 4) 并行求解每个单元的SPMe（使用线程本地精简工作区）
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)

    # 预分配线程本地工作区
    nthreads_avail = Threads.nthreads()
    ws_pool = [create_element_workspace(case) for _ in 1:nthreads_avail]

    # 可选：并行化（如果Julia配置了多线程）
    t_spme_ns = time_ns()
    Threads.@threads for e in 1:ne
        tid = Threads.threadid()
        ws_e = ws_pool[tid]
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e;I_e=I_e[e], T_e=Te_prev[e],jacobi=jacobi, workspace=ws_e)
        M_elems[e] = M_e
        K_elems[e] = K_e
        F_elems[e] = vec(F_e)
        variables_elems[e] = copy_element_results(vars_e)
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
    
    # 5) 装配电化学全局矩阵
    M_chem = blockdiag(M_elems...)
    K_chem = blockdiag(K_elems...)
    F_chem = vcat(F_elems...)
    
    # 6) 计算逐单元热源（调用 ThermalDistributed.jl 中的统一函数)
    t_thermal_ns = time_ns()
    t_czm_model_s = 0.0
    if case.opt.czm_enabled == true && case.czm_mesh !== nothing
        t_czm_ns = time_ns()
        variables = compute_heat_sources_with_czm(case, variables, variables_elems, I_e, Te_prev, areas, case.czm_mesh, mesh_th)
        t_czm_model_s = (time_ns() - t_czm_ns) * 1e-9
        # 记录 CZM 摘要统计
        stats = get_damage_statistics(case.czm_mesh)
        variables["czm D_max"] = [stats.max_D]
        variables["czm D_mean"] = [stats.mean_D]
        δ_max_n_vals = [s.δ_max_n for s in case.czm_mesh.damage_states]
        variables["czm δ_max_n"] = [maximum(δ_max_n_vals)]
        variables["czm δ_mean_n"] = [mean(δ_max_n_vals)]
        variables["czm n_fractured"] = [Float64(stats.n_fractured)]
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
    M = blockdiag(M_chem, MT)
    K = blockdiag(K_chem, KT)
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

"""
    copy_element_results(vars_e)

从 workspace Dict 中提取下游代码需要的键，返回轻量级独立 Dict。
状态提取键（通过 case.index 原位写入 workspace）必须 copy()。
计算结果键（通过 = 赋值，每次创建新数组/标量）可安全引用拷贝。
"""
function copy_element_results(vars_e)
    Dict{String, Union{Array{Float64},Float64}}(
        # ── 计算结果键：每次 = 赋值创建新对象，引用安全 ──
        "negative electrode overpotential"             => vars_e["negative electrode overpotential"],
        "positive electrode overpotential"             => vars_e["positive electrode overpotential"],
        "cell voltage"                                 => vars_e["cell voltage"],
        "negative electrode exchange current density"  => vars_e["negative electrode exchange current density"],
        "positive electrode exchange current density"  => vars_e["positive electrode exchange current density"],
        "negative electrode interfacial current density" => vars_e["negative electrode interfacial current density"],
        "positive electrode interfacial current density" => vars_e["positive electrode interfacial current density"],
        "negative electrode open circuit potential"    => vars_e["negative electrode open circuit potential"],
        "positive electrode open circuit potential"    => vars_e["positive electrode open circuit potential"],
        "temperature"                                  => vars_e["temperature"],
        "electrolyte lithium concentration at negative electrode Gauss point" => vars_e["electrolyte lithium concentration at negative electrode Gauss point"],
        "electrolyte lithium concentration at positive electrode Gauss point" => vars_e["electrolyte lithium concentration at positive electrode Gauss point"],
        "electrolyte lithium concentration at separator Gauss point" => vars_e["electrolyte lithium concentration at separator Gauss point"],
        # ── 状态提取键：替换赋值，workspace 不再持有引用，直接传递 ──
        "negative particle surface lithium concentration" => vars_e["negative particle surface lithium concentration"],
        "positive particle surface lithium concentration" => vars_e["positive particle surface lithium concentration"],
        "negative particle lithium concentration"      => vars_e["negative particle lithium concentration"],
        "positive particle lithium concentration"      => vars_e["positive particle lithium concentration"],
    )
end
