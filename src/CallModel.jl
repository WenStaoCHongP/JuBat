"""
	map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, damage_states, ne_thermal::Int) -> Vector{Float64}

按 spec v2 §5.2 把 cohesive 单元损伤 max 归约到粗热单元粒度。
对每个粗热单元 e_thermal，扫描所有 cohesive_to_thermal[e_coh] == e_thermal 的 cohesive 单元，
取 damage_states[e_coh].D 的最大值；无 cohesive 覆盖则取 0。
（2026-08-30 重构：damage_states 来自 MechState，显式传入。）
"""
function map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, damage_states::AbstractVector{DamageState}, ne_thermal::Int)
	thermal_elem_map = czm_mesh.czm_submesh.thermal_elem_map
	n_coh = length(czm_mesh.cohesive_elements)
	D_per_thermal = zeros(ne_thermal)
	for e_coh in 1:n_coh
		e_thermal = thermal_elem_map[czm_mesh.cohesive_elements[e_coh].host_outer_elem]
		D = damage_states[e_coh].D
		if D > D_per_thermal[e_thermal]
			D_per_thermal[e_thermal] = D
		end
	end
	return D_per_thermal
end

# 任务63 B2（R4/5.5）：可选分流消费观察点。默认 nothing 惰性；诊断入口在
# finally 中恢复，生产路径零开销。事件字段见 CallModel_MultiSPMe 内回调处。
const branch_consumption_observer = Ref{Any}(nothing)

# 任务63 R0-c/M0：球体积平均（FE/高斯口径，与 ElectrodeDiffusion 的 r² 质量装配同一
# 积分：q = gs.x²·weight·detJ，c_bar = Σ q_gp·c_h(r_gp)/Σ q_gp = 1ᵀMc/1ᵀM1）。
# 网格必须是 NE/PE 对应的 negative/positive particle 径向网格（Nrn/Nrp 决定），
# 不得使用电极厚度坐标网格（R1 返工：后者长度 Nn/Np+1，合法不同配置下 DimensionMismatch）。
function particle_volume_average(c::AbstractVector{Float64}, mesh::Mesh)
    c_gs = vec(sum(mesh.gs.Ni .* c[mesh.element[mesh.gs.ele, :]], dims=2))
    q = vec(mesh.gs.x).^2 .* mesh.gs.weight .* mesh.gs.detJ
    den = sum(q)
    c_bar = dot(c_gs, q) / den
    return c_bar
end

function CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String, mech::Union{Nothing,MechState})
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)
    n_chem = case.index["SPMe chemical block length"]
    nT = mesh_th.nlen
    param = case.param
    
    # 唯一完整输出容器；电化学量由各单元真实结果填入。
    variables = StandardVariables(case, 1)
    
    # 1) 解析状态向量
    # 提取热场
    T_nodes = get_thermal_dofs(yt, case)
    # 提取每个单元的电化学状态（使用 @views 消除拷贝）
    yt_chem = Vector{SubArray{Float64,1}}(undef, ne)
    for e in 1:ne
        offset = (e - 1) * n_chem
        yt_chem[e] = @view yt[(offset + 1):(offset + n_chem)]
    end
    
    # 2) 使用缓存的元素面积（网格不变量）
    areas = case.mesh["thermal2D element area"]
    
    Te_prev = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        Te_prev[e] = sum(T_nodes[nds]) / length(nds)
    end
    # 3) 分流求解（获取 I_e）
    # 使用与 SPMe 变量一致的电流无量纲尺度（param.scale.I_typ）
    I_applied = case.opt.Current(t * case.param.scale.t0)
    I_total = I_applied / case.param.scale.I_typ
    
    # 结果电流沿用 I1C 尺度；分流方程使用上面的 I_typ 尺度。
    variables["cell current"] = I_applied / case.param_dim.cell.I1C
    variables["thermal2D temperature at nodes"] = T_nodes
    variables["thermal2D element area"] = areas

    # 渐进式面积损失的 D 映射（唯一损伤消费入口；fractured 不再生成独立停流列表，任务59 P1）
    # R1 返工：反馈源已存在（挂了 czm_mesh）时结构必须完整合法——缺 mech/缺映射均显式失败，
    # 不再按 czm.enabled 分叉；完全没有反馈源（未挂 czm_mesh 且 CZM 关）才是合法几何权重路径。
    D_elem_area_loss = nothing
    if case.opt.czm.area_loss_enabled
        if haskey(case.mesh, "czm")
            D_elem_area_loss = map_czm_damage_to_thermal(case.mesh["czm"], mech.damage_states, ne)
        elseif case.opt.czm.enabled
            error("CallModel_MultiSPMe: area_loss_enabled 且 CZM 启用，但缺少 czm_mesh（非法状态，不得静默回退几何权重）")
        end
        # 其余：无反馈源（未挂 czm_mesh、CZM 关）→ 合法普通路径，按几何权重
    end

    t_branch_ns = time_ns()
    variables, I_e, Vc = solve_branch_currents(case, variables, t, I_total, areas, Te_prev, nothing; D_elem=D_elem_area_loss, yt_chem=yt_chem)
    t_branch_s = (time_ns() - t_branch_ns) * 1e-9

    # Voltage-exhausted positive area cannot carry the requested current.
    # Keep status4's legitimate zero-current diffusion/thermal mode intact.
    if variables["thermal2D Vsolve status"] == 5.0
        n_active = variables["thermal2D n_active_elements"]
        error("MultiSPMe voltage cutoff exhausted the active set: status=5, " *
            "n_active=$n_active, " *
            "I_total=$I_total (normalized), t=$t (normalized), " *
            "time_s=$(t * case.param.scale.t0); refusing SPMe/thermal assembly")
    end

    # 任务63 B2（R4/5.5）：可选分流消费观察点——每次实际消费后原位回调，惰性无开销。
    # 事件携带真实消费输入（D_elem、areas）与产出（I_e、w、A_eff、活动集），t 为归一化时间。
    obs_branch = branch_consumption_observer[]
    obs_branch === nothing || obs_branch((; t=t, I_total=I_total, areas=areas,
        D_elem=D_elem_area_loss, I_e=copy(I_e), Vc=Vc, Te_prev=copy(Te_prev),
        w=copy(variables["thermal2D branch weight"]),
        A_eff=copy(variables["thermal2D effective area"]),
        a=copy(variables["thermal2D branch area factor"]),
        gamma=variables["thermal2D total effective area factor"],
        I_branch=copy(variables["thermal2D element branch current [A]"]),
        active_mask=copy(variables["thermal2D active_mask"]),
        deactivated_mask=copy(variables["thermal2D deactivated_mask"]),
        inactive_reason=copy(variables["thermal2D inactive_reason"])))

    # active_elements：实际分流活动支路索引（面积活跃 ∩ 电压截止；任务59 Q6 裁决）
    variables["active_elements"] = Float64.(findall(==(1.0), variables["thermal2D active_mask"]))
    
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
        M_e, K_e, F_e, vars_e = SPMe(case, yt_chem[e], t; e=e, I_e=I_e[e], T_e=Te_prev[e],jacobi=jacobi, workspace=ws_e)
        M_elems[e] = M_e
        K_elems[e] = K_e
        F_elems[e] = vec(F_e)
        variables_elems[e] = copy_element_results(case, vars_e)
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
    
    # 5) 装配电化学全局矩阵
    M_chem = blockdiag(M_elems...)
    K_chem = blockdiag(K_elems...)
    F_chem = vcat(F_elems...)
    
    # 6) 计算逐单元热源（调用 ThermalDistributed.jl 中的统一函数)
    t_thermal_ns = time_ns()
    t_czm_model_s = 0.0
    if case.opt.czm.enabled == true && haskey(case.mesh, "czm")
        t_czm_ns = time_ns()
        variables = compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas)
        t_czm_model_s = (time_ns() - t_czm_ns) * 1e-9
        # 记录 CZM 摘要统计
        stats = get_damage_statistics(mech.damage_states)
        variables["czm D_max"] = [stats.max_D]
        if case.opt.czm.winding_prestress
            variables["winding prestress"] = [1.0]   # 浮点历史标记；最终结果再解释为功能开启
        end
        variables["czm D_mean"] = [stats.mean_D]
        δ_max_n_vals = [s.δ_max_n for s in mech.damage_states]
        δ_max_eff_vals = [s.δ_max_eff for s in mech.damage_states]
        variables["czm δ_max_n"] = [maximum(δ_max_n_vals)]
        variables["czm δ_max_eff"] = [maximum(δ_max_eff_vals)]
        variables["czm δ_mean_n"] = [mean(δ_max_n_vals)]
        variables["czm n_fractured"] = [Float64(stats.n_fractured)]
    else
        variables = compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas)
    end

    # 保存辅助变量（用于调试）；M0：体积平均改用 NE/PE 对应 particle 径向网格
    mesh_np = case.mesh["negative particle"]
    mesh_pp = case.mesh["positive particle"]
    element_output_keys = spme_element_output_keys(case)
    for e in 1:ne
        vars_e = variables_elems[e]
        for key in element_output_keys
            value = vars_e[key]
            if value isa Float64
                variables[key][e, 1] = value
            else
                n = size(variables[key], 1) ÷ ne
                variables[key][((e - 1) * n + 1):(e * n), 1] = vec(value)
            end
        end
        variables["thermal2D eta_n_e"][e] = vars_e["negative electrode overpotential"][1]
        variables["thermal2D eta_p_e"][e] = vars_e["positive electrode overpotential"][end]
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        variables["thermal2D dUdT_n_e"][e] = param.NE.dUdT(cn_surf_e)[1]
        variables["thermal2D dUdT_p_e"][e] = param.PE.dUdT(cp_surf_e)[1]
        csn_data = vars_e["negative particle lithium concentration"]
        csp_data = vars_e["positive particle lithium concentration"]
        # M1（2026-10-12 用户裁决）：生产化学膨胀消费切换为 FE 球体积平均（1ᵀMc/1ᵀM1，
        # 与扩散质量装配同一积分；静置零通量时库存守恒、平均严格不变）。
        # 旧节点均值保留为有界诊断键（*_nodeavg）。
        vol_n = particle_volume_average(vec(csn_data), mesh_np)
        vol_p = particle_volume_average(vec(csp_data), mesh_pp)
        variables["thermal2D element soc_n"][e] = vol_n
        variables["thermal2D element soc_p"][e] = vol_p
        variables["thermal2D element soc_n_vol"][e] = vol_n   # 同一正确值，不再各算一套
        variables["thermal2D element soc_p_vol"][e] = vol_p
        variables["thermal2D element soc_n_nodeavg"][e] = mean(vec(csn_data))
        variables["thermal2D element soc_p_nodeavg"][e] = mean(vec(csp_data))
    end
    variables["thermal2D element current"] = I_e
    
    # 7) 装配热学矩阵
    MT, KT, FT = ThermalDistributed2D(case, variables)
    # 统一时间尺度：电化学和热模型使用相同时间尺度 t0，时间比为 1
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
    active_idx = findall(==(1.0), variables["thermal2D active_mask"])
    if !isempty(active_idx)
        residual = maximum(abs(variables["thermal2D element voltages"][e] - Vc) for e in active_idx)
        residual <= 1e-8 || error("MultiSPMe local voltage mismatch: residual=$residual")
    end
    variables["timing branch solver [s]"] = t_branch_s
    variables["timing spme solve [s]"] = t_spme_s
    variables["timing thermal distributed [s]"] = t_thermal_s
    variables["timing czm model [s]"] = t_czm_model_s
    
    y_phi = Float64[]
    
    return M, K, F, variables, y_phi
end
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String, mech::Union{Nothing,MechState}=nothing)
    # 多SPMe模式：由 per_element_spme 控制，直接委托给 CallModel_MultiSPMe
    if case.opt.per_element_spme
        return CallModel_MultiSPMe(case, yt, t; jacobi=jacobi,mech=mech)
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
    copy_element_results(case, vars_e)

从 workspace Dict 中保留逐单元输出及热源所需的结果。
这些键均由 SPMe_variables / Mechanicaloutput 替换赋值，后续 workspace
复用不会原位修改已保存数组；包括 full 模式真实颗粒应力结果。
"""
function copy_element_results(case::Case, vars_e)
    result = Dict{String, Union{Array{Float64},Float64}}(
        key => vars_e[key] for key in spme_element_output_keys(case)
    )
    result["cell voltage"] = vars_e["cell voltage"]
    result["temperature"] = vars_e["temperature"]
    return result
end
