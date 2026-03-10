# ========================================================================
# solve_branch_currents_newton - 辅助函数
# ========================================================================

# 调试输出函数
"""检查并报告 NaN/Inf 值（调试用）"""
function _debug_check_prefactors(prefactor_n, prefactor_p, csn_av, csp_av, u_n_ref_val, u_p_ref_val, du_n_dT_val, du_p_dT_val, c_sigma, cn_surf, cp_surf, ce_n_gs, ce_p_gs)
	has_nan = !isfinite(prefactor_n) || !isfinite(prefactor_p) || !isfinite(csn_av) || !isfinite(csp_av) || !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val) || !isfinite(c_sigma)
    
	if !has_nan
		return false
	end
    
	println("\n" * "="^80)
	println("❌ [DEBUG] 预计算值包含 NaN/Inf")
	println("="^80)
	println("📊 预因子: prefactor_n=$prefactor_n, prefactor_p=$prefactor_p")
	println("📊 平均浓度: csn_av=$csn_av, csp_av=$csp_av")
	println("📊 开路电位: u_n=$u_n_ref_val, u_p=$u_p_ref_val")
	println("📊 温度系数: du_n_dT=$du_n_dT_val, du_p_dT=$du_p_dT_val")
	println("📊 电导: c_sigma=$c_sigma")
	println("\n💡 输入浓度前5个值:")
	println("  cn_surf: $(cn_surf[1:min(5,end)])")
	println("  cp_surf: $(cp_surf[1:min(5,end)])")
	println("  ce_n_gs: $(ce_n_gs[1:min(5,end)])")
	println("  ce_p_gs: $(ce_p_gs[1:min(5,end)])")
	println("="^80 * "\n")
    
	return true
end

"""检查单元系数是否有效"""
function _debug_check_coefficients(e, has_nan_prefactor, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n_val_T, u_p_val_T)
	if has_nan_prefactor || e > 1
		return  # 只打印第一个异常单元，且预计算值正常时
	end
    
	has_nan_coeff = !isfinite(C1) || !isfinite(C2) || !isfinite(alpha_p) || !isfinite(alpha_n) || !isfinite(C5)
    
	if !has_nan_coeff
		return
	end
    
	println("\n" * "="^80)
	println("❌ [DEBUG] 单元 $e 系数异常")
	println("="^80)
	println("📊 T_e=$T_e, j0_n=$j0_n, j0_p=$j0_p")
	println("📊 u_n(T)=$u_n_val_T, u_p(T)=$u_p_val_T")
	println("📊 系数: C1=$C1, C2=$C2")
	println("   alpha_p=$alpha_p, alpha_n=$alpha_n, C5=$C5")
	println("="^80 * "\n")
end

"""检查初始电压是否有效"""
function _debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, coeffs, I_total, ne)
	if has_nan_prefactor || isfinite(V)
		return
	end
    
	println("\n" * "="^80)
	println("❌ [DEBUG] 初始电压异常: V=$V")
	println("="^80)
	println("  I_total=$I_total, ne=$ne")
	println("  前3个异常单元:")
    
	count = 0
	for e in 1:ne
		if !isfinite(V_branches[e]) && count < 3
			count += 1
			println("  单元 $e: V=$(V_branches[e]), I=$(I_e[e])")
			println("    C1=$(coeffs[e].C1), C2=$(coeffs[e].C2)")
			println("    α_p=$(coeffs[e].alpha_p), α_n=$(coeffs[e].alpha_n)")
		end
	end
	println("="^80 * "\n")
end

# 电化学计算函数
"""标量化：将数组转为标量（取第一个元素）"""
_scalarize(x) = isa(x, Number) ? Float64(x) : Float64(x[1])

"""计算电化学预因子"""
function _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
	cn_surf = variables["negative particle surface lithium concentration"]
	cp_surf = variables["positive particle surface lithium concentration"]
	ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
	ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
    
	# 预因子计算
	prefactor_n = IntV(abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5, mesh_ne) / param.NE.thickness
	prefactor_p = IntV(abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5, mesh_pe) / param.PE.thickness
	csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
	csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    
	# 开路电位
	u_n_ref = param.NE.U(cn_surf)
	u_p_ref = param.PE.U(cp_surf)
	du_n_dT = param.NE.dUdT(cn_surf)
	du_p_dT = param.PE.dUdT(cp_surf)
    
	# 标量化
	u_n_ref_val = _scalarize(u_n_ref)
	u_p_ref_val = _scalarize(u_p_ref)
	du_n_dT_val = _scalarize(du_n_dT)
	du_p_dT_val = _scalarize(du_p_dT)
    
	# 固相电导
	c_sigma = (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig) / 3.0
    
	return (prefactor_n=prefactor_n, prefactor_p=prefactor_p,
			csn_av=csn_av, csp_av=csp_av,
			u_n_ref_val=u_n_ref_val, u_p_ref_val=u_p_ref_val,
			du_n_dT_val=du_n_dT_val, du_p_dT_val=du_p_dT_val,
			c_sigma=c_sigma,
			cn_surf=cn_surf, cp_surf=cp_surf, ce_n_gs=ce_n_gs, ce_p_gs=ce_p_gs)
end

"""计算单个单元的电化学系数"""
function _compute_element_coefficients(e, T_e, param, prefactors, T_ref, debug_mode=false)
	# 交换电流密度
	arr_n = Arrhenius(param.NE.Eac_k, T_e)
	arr_p = Arrhenius(param.PE.Eac_k, T_e)
	j0_n = param.NE.k * arr_n * prefactors.prefactor_n
	j0_p = param.PE.k * arr_p * prefactors.prefactor_p
    
	# 电解液电导率
	kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps^param.NE.brugg
	kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps^param.PE.brugg
	kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps^param.SP.brugg
	R_EL = param.NE.thickness / (3.0 * kappa_ne) + param.SP.thickness / kappa_sp + param.PE.thickness / (3.0 * kappa_pe)
    
	# 开路电位（温度相关）
	u_n = prefactors.u_n_ref_val + (T_e - T_ref) * prefactors.du_n_dT_val
	u_p = prefactors.u_p_ref_val + (T_e - T_ref) * prefactors.du_p_dT_val
    
	# 计算系数
	C1 = (u_p - u_n) + 2.0 * T_e * (1.0 - param.EL.tplus) * (prefactors.csp_av - prefactors.csn_av) / param.EL.ce0
	C2 = 2.0 * T_e
	alpha_p = -1.0 / (2.0 * j0_p * param.PE.as * param.PE.thickness)
	alpha_n = 1.0 / (2.0 * j0_n * param.NE.as * param.NE.thickness)
	C5 = R_EL + prefactors.c_sigma
    
	# 调试检查（仅第一个异常单元）
	if debug_mode
		_debug_check_coefficients(e, false, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n, u_p)
	end
    
	return (C1=C1, C2=C2, alpha_p=alpha_p, alpha_n=alpha_n, C5=C5)
end

"""批量计算所有单元的系数"""
function _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode=false)
	coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
	for e in 1:ne
		coeffs[e] = _compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref, debug_mode)
	end
	return coeffs
end

# 分支电压模型
"""计算分支电压 V = C1 + C2*(asinh(α_p*I) - asinh(α_n*I)) - C5*I"""
function _branch_voltage(coeff, I::Float64)
	apI = coeff.alpha_p * I
	anI = coeff.alpha_n * I
	return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
end

"""计算分支电压对电流的导数 dV/dI"""
function _branch_dVdI(coeff, I::Float64)
	apI = coeff.alpha_p * I
	anI = coeff.alpha_n * I
	denom_p = sqrt(1.0 + apI * apI)
	denom_n = sqrt(1.0 + anI * anI)
	return coeff.C2 * (coeff.alpha_p / denom_p - coeff.alpha_n / denom_n) - coeff.C5
end

# 初始化和边界检查
"""初始化单元电流猜测"""
function _initialize_currents(ne, w, I_total, x_prev)
	I_e = x_prev !== nothing && length(x_prev) == ne ? copy(x_prev) : (w .* I_total)
    
	# 归一化到满足总电流约束
	sx = sum(w .* I_e)
	if sx != 0.0
		I_e .*= (I_total / sx)
	else
		I_e .= w .* I_total
	end
    
	return I_e
end

"""检查电压边界"""
function _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="")
	V_phys = V * phi_scale
	V_MIN_phys = V_MIN * phi_scale
	V_MAX_phys = V_MAX * phi_scale
    
	if V < V_MIN
		# 放电截止
		return false, 1, V_phys, V_MIN_phys
	elseif V > V_MAX
		# 充电截止
		return false, 2, V_phys, V_MAX_phys
	end

	return true, 0, V_phys, 0.0
end

# 牛顿迭代求解器
"""
	_detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)

检测达到截止电压的单元。

# 返回
- `active_mask`: 布尔数组，true 表示单元活跃（可接受电流）
- `n_cutoff`: 达到截止的单元数
- `cutoff_info`: 截止单元详细信息字典

# 物理含义
- 充电时 (I_total < 0)：如果单元的 OCV >= V_MAX，该单元已满充
- 放电时 (I_total > 0)：如果单元的 OCV <= V_MIN，该单元已完全放电
- 静置时 (I_total ≈ 0)：所有单元都是活跃的
"""
function _detect_cutoff_elements(coeffs, ne::Int, V_MIN::Float64, V_MAX::Float64, I_total::Float64, phi_scale::Float64)
	active_mask = trues(ne)
    
	# 计算各单元的开路电压 (OCV = C1，当 I=0 时的电压)
	OCV = [coeffs[e].C1 for e in 1:ne]
	 # 记录截止单元的详细信息
	cutoff_elements = Int[]
	cutoff_ocv = Float64[]
	cutoff_type = Int[]  # 1 = 放电截止, 2 = 充电截止
	 # 严格检测：无容差
	if I_total < 0  # 充电
		# 充电时，如果单元 OCV >= V_MAX，该单元已满充
		for e in 1:ne
			if OCV[e] >= V_MAX
				active_mask[e] = false
				push!(cutoff_elements, e)
				push!(cutoff_ocv, OCV[e] * phi_scale)
				push!(cutoff_type, 2)
			end
		end
	elseif I_total > 0  # 放电
		# 放电时，如果单元 OCV <= V_MIN，该单元已完全放电
		for e in 1:ne
			if OCV[e] <= V_MIN
				active_mask[e] = false
				push!(cutoff_elements, e)
				push!(cutoff_ocv, OCV[e] * phi_scale)
				push!(cutoff_type, 1)
			end
		end
	end
	# 静置时 (I_total ≈ 0)，所有单元保持活跃
    
	n_cutoff = sum(.!active_mask)
	# 构建截止信息字典
	cutoff_info = Dict{String, Any}(
		"cutoff_elements" => cutoff_elements,
		"cutoff_ocv" => cutoff_ocv,
		"cutoff_type" => cutoff_type,
		"all_ocv" => OCV .* phi_scale,
		"V_MIN" => V_MIN * phi_scale,
		"V_MAX" => V_MAX * phi_scale
	)
    
	# 找出最接近截止的单元（用于预警）
	if I_total > 0  # 放电
		min_ocv_idx = argmin(OCV)
		cutoff_info["nearest_cutoff_element"] = min_ocv_idx
		cutoff_info["nearest_cutoff_ocv"] = OCV[min_ocv_idx] * phi_scale
		cutoff_info["margin_to_cutoff"] = (OCV[min_ocv_idx] - V_MIN) * phi_scale
	elseif I_total < 0  # 充电
		max_ocv_idx = argmax(OCV)
		cutoff_info["nearest_cutoff_element"] = max_ocv_idx
		cutoff_info["nearest_cutoff_ocv"] = OCV[max_ocv_idx] * phi_scale
		cutoff_info["margin_to_cutoff"] = (V_MAX - OCV[max_ocv_idx]) * phi_scale
	end

	return active_mask, n_cutoff, cutoff_info
end

"""牛顿迭代主循环（支持部分单元截止）

当所有单元都活跃时（n_cutoff=0），行为与原实现完全一致。
当有截止单元时，只对活跃单元求解，截止单元电流保持为 0。
"""
function _newton_iteration!(I_e, V, ne, w, I_total, coeffs; tol_V=1e-8, tol_I=1e-10, max_iters=25,active_mask::Union{Nothing, BitVector}=nothing)
	converged = false
	last_iter = 0
	F = zeros(Float64, ne)
	dFdI = similar(F)
	I_trial = similar(I_e)
	ΔI = similar(I_e)
    
	# 如果没有提供 active_mask，所有单元都是活跃的
	if active_mask === nothing
		active_mask = trues(ne)
	end
    
	# 获取活跃单元索引
	active_idx = findall(active_mask)
	n_active = length(active_idx)
	all_active = (n_active == ne)

	# 如果没有活跃单元，直接返回
	if n_active == 0
		return V, true, 0
	end

	for iter in 1:max_iters
		last_iter = iter
        
		# 计算残差和雅可比
		for e in 1:ne
			V_e = _branch_voltage(coeffs[e], I_e[e])
			F[e] = V_e - V
			dFdI[e] = _branch_dVdI(coeffs[e], I_e[e])
            
			# 防止奇异雅可比
			if abs(dFdI[e]) < 1e-12
				dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5
			end
		end
        
		# 检查收敛
		if all_active
			res_V = maximum(abs.(F))
		else
			res_V = maximum(abs.(F[active_idx]))
		end
		res_I = sum(w .* I_e) - I_total
        
		if res_V <= tol_V && abs(res_I) <= tol_I
			converged = true
			break
		end
        
		# 牛顿步
		if all_active
			# 所有单元活跃：使用原始向量操作
			denom = sum(w ./ dFdI)
			abs(denom) < 1e-12 && break
            
			num = -res_I + sum(w .* F ./ dFdI)
			ΔV = num / denom
			ΔI .= ((-F) .+ ΔV) ./ dFdI
		else
			# 有截止单元：只对活跃单元计算
			denom = sum(w[active_idx] ./ dFdI[active_idx])
			abs(denom) < 1e-12 && break
            
			num = -res_I + sum(w[active_idx] .* F[active_idx] ./ dFdI[active_idx])
			ΔV = num / denom
            
			# 计算电流增量（只对活跃单元）
			fill!(ΔI, 0.0)
			for e in active_idx
				ΔI[e] = ((-F[e]) + ΔV) / dFdI[e]
			end
		end
        
		# 线搜索
		λ, V_trial = _line_search(I_e, V, ΔI, ΔV, I_trial, ne)
		λ == 0.0 && break

		# 更新
		if all_active
			# 所有单元活跃：使用原始向量操作
			I_e .= I_trial
		else
			# 有截止单元：只更新活跃单元
			for e in active_idx
				I_e[e] = I_trial[e]
			end
			# 非活跃单元保持 I_e[e] = 0
		end
		V = V_trial
	end
    
	return V, converged, last_iter
end

"""线搜索：确保更新后的值有效"""
function _line_search(I_e, V, ΔI, ΔV, I_trial, ne; max_attempts=12)
	λ = 1.0
    
	for attempt in 1:max_attempts
		# 试探电压
		V_trial = V + λ * ΔV
		!isfinite(V_trial) && (λ *= 0.5; continue)
        
		# 试探电流
		ok = true
		@inbounds for e in 1:ne
			val = I_e[e] + λ * ΔI[e]
			if !isfinite(val) || abs(val) > 1e12
				ok = false
				break
			end
			I_trial[e] = val
		end
        
		ok && return λ, V_trial
		λ *= 0.5
	end
    
	return 0.0, V  # 失败
end

# ========================================================================
# 主函数：solve_branch_currents_newton（精简版）
# ========================================================================

"""
	solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev; deactivated_elements=nothing)

非线性分流求解器（精简版）

使用牛顿法求解电流分配问题：
- 未知量：各单元电流 I_e 和公共端电压 V
- 约束：V_e(I_e, T_e) = V (每个单元) 和 Σ(I_e*A_e) = I_total

# 参数
- `case`: 案例对象
- `variables`: 变量字典
- `yt`: 状态向量
- `t`: 时间
- `I_total`: 总电流（无量纲）
- `areas`: 各单元面积
- `Te_prev`: 各单元温度
- `x_prev`: 上一步的电流分配（可选）
- `deactivated_elements`: 失效单元索引列表（可选，CZM断裂导致的）

# 返回
- `variables`: 更新后的变量字典
- `I_e`: 各单元电流
- `V`: 公共端电压

# CZM断裂处理
当 `deactivated_elements` 不为空时，这些单元被视为永久退出电化学反应：
- 失效单元的电流固定为零
- 总电流由剩余活跃单元承担
- 失效单元不参与牛顿迭代
"""
function solve_branch_currents_newton(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing; deactivated_elements::Union{Nothing,Vector{Int64}}=nothing)
	# 1. 初始化
	ne = length(areas)
	A_global = sum(areas)
	w = areas ./ A_global  # 面积权重
	phi_scale = case.param.scale.phi
	V_MIN = case.param_dim.cell.v_l / phi_scale
	V_MAX = case.param_dim.cell.v_h / phi_scale
	# CZM失效单元处理：创建失效掩码
	deactivated_mask = zeros(Bool, ne)
	if deactivated_elements !== nothing && !isempty(deactivated_elements)
		for e in deactivated_elements
			if 1 <= e <= ne
				deactivated_mask[e] = true
			end
		end
	end
	n_deactivated = sum(deactivated_mask)

	# 2. 计算电化学预因子
	param = case.param
	mesh_ne = case.mesh["negative electrode"]
	mesh_pe = case.mesh["positive electrode"]
	prefactors = _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
    
	# 调试：检查预因子
	debug_mode = case.opt.debug_coupling
	has_nan_prefactor = _debug_check_prefactors(
		prefactors.prefactor_n, prefactors.prefactor_p, 
		prefactors.csn_av, prefactors.csp_av,
		prefactors.u_n_ref_val, prefactors.u_p_ref_val,
		prefactors.du_n_dT_val, prefactors.du_p_dT_val,
		prefactors.c_sigma,
		prefactors.cn_surf, prefactors.cp_surf, 
		prefactors.ce_n_gs, prefactors.ce_p_gs
	)
    
	# 3. 计算各单元系数
	T_ref = case.param.cell.T0
	coeffs = _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode)
    
	# 4. 检测达到截止电压的单元
	active_mask, n_cutoff, cutoff_info = _detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)
	# 4b. 合并CZM失效单元到活跃掩码
	# 失效单元不参与电化学反应
	for e in 1:ne
		if deactivated_mask[e]
			active_mask[e] = false
		end
	end
	n_inactive_total = sum(.!active_mask)

	# 5. 初始化电流猜测
	I_e = _initialize_currents(ne, w, I_total, x_prev)
    
	# 活跃单元索引
	active_idx = findall(active_mask)
	all_active = (n_inactive_total == 0)
    
	# 将非活跃单元的电流设为0
	if !all_active
		for e in 1:ne
			if !active_mask[e]
				I_e[e] = 0.0
			end
		end
	end
    
	# 计算初始电压
	if all_active
		# 所有单元活跃：使用原始逻辑
		V_branches = [_branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
		V = sum(V_branches) / ne
	elseif !isempty(active_idx)
		# 有截止单元：使用活跃单元的平均值
		V_branches = [_branch_voltage(coeffs[e], I_e[e]) for e in active_idx]
		V = sum(V_branches) / length(active_idx)
	else
		# 所有单元都达到截止，使用 OCV 的平均值作为公共电压
		V = sum(coeffs[e].C1 for e in 1:ne) / ne
	end
    
	# 调试：检查初始电压
	V_branches_all = [_branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
	_debug_check_initial_voltage(has_nan_prefactor, V, V_branches_all, I_e, coeffs, I_total, ne)
    
	# 6. 牛顿迭代求解
	if all_active
		# 所有单元活跃
		V, converged, last_iter = _newton_iteration!(I_e, V, ne, w, I_total, coeffs)
	elseif !isempty(active_idx)
		# 有非活跃单元（截止或失效）：传入 active_mask
		V, converged, last_iter = _newton_iteration!(I_e, V, ne, w, I_total, coeffs; active_mask=active_mask)
	else
		# 所有单元都非活跃
		converged = true
		last_iter = 0
	end
    
	# 7. 归一化确保总电流约束
	sx = sum(w .* I_e)
	if all_active
		# 所有单元活跃：使用原始向量操作
		if sx != 0.0
			I_e .*= (I_total / sx)
		end
	elseif !isempty(active_idx)
		# 有截止单元：只调整活跃单元
		if abs(sx) > 1e-12
			scale_factor = I_total / sx
			for e in active_idx
				I_e[e] *= scale_factor
			end
		elseif abs(I_total) > 1e-12
			# sx ≈ 0 但 I_total ≠ 0：按面积权重分配给活跃单元
			w_active_sum = sum(w[active_idx])
			if w_active_sum > 0
				for e in active_idx
					I_e[e] = I_total * w[e] / w_active_sum
				end
			end
		end
	end
    
	# 8. 边界检查
	voltage_in_bounds, cutoff_type_global, V_phys, V_limit = _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e)
    
	# 9. 写入结果
	variables["thermal2D element current"] = I_e
	variables["thermal2D element current A"] = case.param.scale.I_typ .* I_e
	variables["thermal2D common voltage"] = V
	variables["thermal2D common voltage V"] = V * phi_scale
	variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5
	variables["thermal2D Vsolve iters"] = Float64(last_iter)
	variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0
	# 截止状态信息（包括CZM失效单元）
	variables["thermal2D n_active_elements"] = Float64(sum(active_mask))
	variables["thermal2D n_cutoff_elements"] = Float64(n_cutoff)
	variables["thermal2D n_deactivated_elements"] = Float64(n_deactivated)  # CZM失效单元数
	variables["thermal2D active_mask"] = Float64.(active_mask)  # 1.0 = 活跃, 0.0 = 截止/失效
	variables["thermal2D deactivated_mask"] = Float64.(deactivated_mask)  # 1.0 = CZM失效
	# 失效原因编码: 0=活跃, 1=电压截止, 2=CZM断裂失效
	inactive_reason = zeros(Float64, ne)
	for e in 1:ne
		if deactivated_mask[e]
			inactive_reason[e] = 2.0  # CZM断裂失效
		elseif !active_mask[e]
			inactive_reason[e] = 1.0  # 电压截止
		end
		# 0.0 表示活跃
	end
	variables["thermal2D inactive_reason"] = inactive_reason
	variables["thermal2D voltage_in_bounds"] = voltage_in_bounds ? 1.0 : 0.0
	variables["thermal2D cutoff_type_global"] = Float64(cutoff_type_global)  # 0=正常, 1=放电截止, 2=充电截止
    
	# 各单元的开路电压（用于诊断）
	variables["thermal2D element OCV"] = cutoff_info["all_ocv"]
	# 截止单元详细信息
	variables["thermal2D cutoff_elements"] = Float64.(cutoff_info["cutoff_elements"])
	variables["thermal2D cutoff_ocv"] = cutoff_info["cutoff_ocv"]
	variables["thermal2D cutoff_type"] = Float64.(cutoff_info["cutoff_type"])
    
	# 最接近截止的单元信息（预警）
	if haskey(cutoff_info, "nearest_cutoff_element")
		variables["thermal2D nearest_cutoff_element"] = Float64(cutoff_info["nearest_cutoff_element"])
		variables["thermal2D nearest_cutoff_ocv"] = cutoff_info["nearest_cutoff_ocv"]
		variables["thermal2D margin_to_cutoff"] = cutoff_info["margin_to_cutoff"]
	end
	return variables, I_e, V
end
