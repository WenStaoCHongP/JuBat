# ========================================================================
# solve_branch_currents - 辅助函数
# ========================================================================

# 截止检测结果结构体（替代 Dict{String,Any}）
struct CutoffInfo
    active_mask::BitVector
    n_cutoff::Int
    cutoff_elements::Vector{Int}
    cutoff_ocv::Vector{Float64}
    cutoff_type::Vector{Int}
    all_ocv::Vector{Float64}
    nearest_element::Int
    nearest_ocv::Float64
    margin::Float64
end

# 电化学计算函数
"""标量化：将数组转为标量（取第一个元素）"""
scalarize(x) = if isa(x, Number)
	Float64(x)
else
	Float64(x[1])
end

"""计算电化学预因子"""
function compute_prefactors(variables, param, mesh_ne, mesh_pe, mesh_sp)
	cn_surf = variables["negative particle surface lithium concentration"]
	cp_surf = variables["positive particle surface lithium concentration"]
	ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
	ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
	ce_sp_gs = variables["electrolyte lithium concentration at separator Gauss point"]
    
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
	u_n_ref_val = scalarize(u_n_ref)
	u_p_ref_val = scalarize(u_p_ref)
	du_n_dT_val = scalarize(du_n_dT)
	du_p_dT_val = scalarize(du_p_dT)
    
	# 固相电导
	c_sigma = (param.NE.thickness / (param.NE.sig * param.NE.eps_s) + param.PE.thickness / (param.PE.sig * param.PE.eps_s)) / 3.0
    
	return (prefactor_n=prefactor_n, prefactor_p=prefactor_p,
			csn_av=csn_av, csp_av=csp_av,
			u_n_ref_val=u_n_ref_val, u_p_ref_val=u_p_ref_val,
			du_n_dT_val=du_n_dT_val, du_p_dT_val=du_p_dT_val,
			c_sigma=c_sigma,
			cn_surf=cn_surf, cp_surf=cp_surf, ce_n_gs=ce_n_gs, ce_p_gs=ce_p_gs,
			ce_sp_gs=ce_sp_gs, mesh_ne=mesh_ne, mesh_pe=mesh_pe, mesh_sp=mesh_sp)
end

"""计算单个单元的电化学系数"""
function compute_element_coefficients(e, T_e, param, prefactors, T_ref)
	# 交换电流密度
	arr_n = Arrhenius(param.NE.Eac_k, T_e)
	arr_p = Arrhenius(param.PE.Eac_k, T_e)
	j0_n = param.NE.k * arr_n * prefactors.prefactor_n
	j0_p = param.PE.k * arr_p * prefactors.prefactor_p
    
	# 与局部 SPMe 一致：在当前 Gauss 浓度和支路温度上评价后积分平均。
	kappa_ne = IntV(param.EL.kappa(prefactors.ce_n_gs, T_e) * param.NE.eps^param.NE.brugg, prefactors.mesh_ne) / param.NE.thickness
	kappa_pe = IntV(param.EL.kappa(prefactors.ce_p_gs, T_e) * param.PE.eps^param.PE.brugg, prefactors.mesh_pe) / param.PE.thickness
	kappa_sp = IntV(param.EL.kappa(prefactors.ce_sp_gs, T_e) * param.SP.eps^param.SP.brugg, prefactors.mesh_sp) / param.SP.thickness
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

	return (C1=C1, C2=C2, alpha_p=alpha_p, alpha_n=alpha_n, C5=C5)
end

"""批量计算所有单元的系数"""
function compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref)
	coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
	for e in 1:ne
		coeffs[e] = compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref)
	end
	return coeffs
end

# 分支电压模型
"""计算分支电压 V = C1 + C2*(asinh(α_p*I) - asinh(α_n*I)) - C5*I"""
function branch_voltage(coeff, I::Float64)
	apI = coeff.alpha_p * I
	anI = coeff.alpha_n * I
	return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
end

"""计算分支电压对电流的导数 dV/dI"""
function branch_dVdI(coeff, I::Float64)
	apI = coeff.alpha_p * I
	anI = coeff.alpha_n * I
	denom_p = sqrt(1.0 + apI * apI)
	denom_n = sqrt(1.0 + anI * anI)
	return coeff.C2 * (coeff.alpha_p / denom_p - coeff.alpha_n / denom_n) - coeff.C5
end

# 初始化和边界检查
"""初始化单元电流猜测"""
function initialize_currents(ne, w, I_total, x_prev; active_mask=trues(ne))
	I_e = if x_prev === nothing
		zeros(Float64, ne)
	else
		copy(x_prev)
	end
	length(I_e) == ne || throw(DimensionMismatch("branch current initial guess length"))
	all(isfinite, I_e) || error("nonfinite branch current initial guess")
	I_e[.!active_mask] .= 0.0
	active_idx = findall(active_mask)
	isempty(active_idx) && return I_e
	# Uniform density shift satisfies the linear constraint without erasing rest balancing currents.
	correction = (I_total - sum(w .* I_e)) / sum(w[active_idx])
	I_e[active_idx] .+= correction
	all(isfinite, I_e) || error("nonfinite branch current initial guess")
	return I_e
end

"""检查电压边界"""
function check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="")
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
	detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)

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
function detect_cutoff_elements(coeffs, ne::Int, V_MIN::Float64, V_MAX::Float64, I_total::Float64, phi_scale::Float64)
	active_mask = trues(ne)
	OCV = [coeffs[e].C1 for e in 1:ne]
	cutoff_elements = Int[]
	cutoff_ocv = Float64[]
	cutoff_type = Int[]

	# 检测截止单元：充电时 OCV ≥ V_MAX 满，放电时 OCV ≤ V_MIN 空
	for e in 1:ne
		cutoff = (I_total < 0 && OCV[e] >= V_MAX) || (I_total > 0 && OCV[e] <= V_MIN)
		if cutoff
			active_mask[e] = false
			push!(cutoff_elements, e)
			push!(cutoff_ocv, OCV[e] * phi_scale)
			push!(cutoff_type, if I_total < 0
				2
			else
				1
			end)
		end
	end

	n_cutoff = length(cutoff_elements)
	all_ocv = OCV .* phi_scale

	# 最近截止预警
	nearest_element, nearest_ocv, margin = 0, NaN, NaN
	if I_total > 0
		idx = argmin(OCV)
		nearest_element, nearest_ocv, margin = idx, OCV[idx] * phi_scale, (OCV[idx] - V_MIN) * phi_scale
	elseif I_total < 0
		idx = argmax(OCV)
		nearest_element, nearest_ocv, margin = idx, OCV[idx] * phi_scale, (V_MAX - OCV[idx]) * phi_scale
	end

	return CutoffInfo(active_mask, n_cutoff, cutoff_elements, cutoff_ocv, cutoff_type, all_ocv, nearest_element, nearest_ocv, margin)
end

"""牛顿迭代主循环（支持部分单元截止）

当所有单元都活跃时（n_cutoff=0），行为与原实现完全一致。
当有截止单元时，只对活跃单元求解，截止单元电流保持为 0。
"""
function newton_iteration(I_e, V, ne, w, I_total, coeffs; tol_V=1e-8, tol_I=1e-10, max_iters=25,active_mask::Union{Nothing, BitVector}=nothing)
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
			V_e = branch_voltage(coeffs[e], I_e[e])
			F[e] = V_e - V
			dFdI[e] = branch_dVdI(coeffs[e], I_e[e])
            
			# 防止奇异雅可比
			if abs(dFdI[e]) < 1e-12
				dFdI[e] = if sign(dFdI[e]) != 0.0
					sign(dFdI[e]) * 1e-12
				else
					-coeffs[e].C5
				end
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
		λ, V_trial = line_search(I_e, V, ΔI, ΔV, I_trial, ne)
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
    
	# 最后一次更新后的实际返回量必须同时满足两类残差；不再事后缩放电流。
	res_V = maximum(abs(branch_voltage(coeffs[e], I_e[e]) - V) for e in active_idx)
	res_I = sum(w .* I_e) - I_total
	converged = res_V <= tol_V && abs(res_I) <= tol_I
	return V, converged, last_iter
end

"""线搜索：确保更新后的值有效"""
function line_search(I_e, V, ΔI, ΔV, I_trial, ne; max_attempts=12)
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
# 主函数：solve_branch_currents（精简版）
# 任务59 P1：无阈值渐进面积损失（factor=1−D）统一损伤反馈；独立 fractured 停流已删除。
# ========================================================================

"""
	solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev; yt_chem, D_elem=nothing)

非线性分流求解器（精简版）。活动集 = 面积活跃 ∩ 电压截止活跃。

使用牛顿法求解电流分配问题：
- 未知量：各单元电流 I_e 和公共端电压 V
- 约束：V_e(I_e, yt_chem[e], T_e) = V (每个活跃单元) 和 Σ(a_e*I_e) = I_total（a=A_eff/Σareas）
- 活动面积非零时，静置沿用同一 Newton，仅将总电流设为零，允许内部均衡电流。
- yt_chem 为必需输入；逐支路直接提取浓度，不回退到平均状态。

# 面积损失（任务59 P1，无阈值）
`area_loss_enabled=true` 且传入 `D_elem` 时：A_eff = areas·(1−D_elem)，
a=A_eff/Σareas，γ=ΣA_eff/Σareas，w=A_eff/ΣA_eff；Newton目标为I_total/γ。
A_eff==0 的支路从 Newton 活动集排除、I_e 严格为 0。
NaN/越界 D 显式失败，不当作无损或零面积。

# 全失活（ΣA_eff == 0）
- I_total != 0：抛 ErrorException（面积全失活，无法承载外加电流）。
- I_total == 0：静置诊断分支——I_e 全 0、允许零电流扩散/导热，跳过 w 归一化与
  Newton；V=各支路 OCV(C1) 算术均值（诊断值，非求解所得公共端电压），
  写 Vsolve status=4.0 / converged=0.0 / iters=0.0，上层据此绕过电压截止。

# 参数
- `case`: 案例对象
- `variables`: 变量字典
- `yt`: 状态向量
- `t`: 时间
- `I_total`: 总电流（无量纲）
- `areas`: 各单元面积
- `Te_prev`: 各单元温度
- `x_prev`: 上一步的电流分配（可选）
- `D_elem`: 热单元 max(D) 映射（可选，仅 area_loss_enabled 时消费）
- `yt_chem`: 各支路电化学状态（必需，长度与单元数一致）

# 返回
- `variables`: 更新后的变量字典
- `I_e`: 各单元电流
- `V`: 公共端电压（全失活静置时为 OCV 均值诊断值）
"""
function solve_branch_currents(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing; yt_chem, D_elem::Union{Nothing,Vector{Float64}}=nothing)
	ne = length(areas)
	# 渐进式有效面积损失：无阈值 factor = 1 - D（仅开关开启时消费 D）
	if case.opt.czm.area_loss_enabled && D_elem !== nothing
		all(d -> isfinite(d) && 0.0 <= d <= 1.0, D_elem) || error(
			"solve_branch_currents: D_elem 含 NaN/越界值（合法域 [0,1]），面积损失路径拒绝消费")
		A_eff = areas .* (1.0 .- D_elem)
	else
		A_eff = areas
	end
	area_sum = sum(areas)
	all(x -> isfinite(x) && x > 0.0, areas) && isfinite(area_sum) && area_sum > 0.0 ||
		throw(ArgumentError("branch reference areas must be finite and positive"))
	a = A_eff ./ area_sum
	gamma = sum(A_eff) / area_sum
	variables["thermal2D branch area factor"] = a
	variables["thermal2D total effective area factor"] = gamma
	area_active = A_eff .> 0.0

	phi_scale = case.param.scale.phi
	# 截止集与 Solve 始终读取同一 param 电压窗口，充电/放电/静置不改写边界。
	# 伏特/phi_scale = 无量纲，单位一致；不引入第二套阈值。
	V_MIN, V_MAX = case.param.cell.v_l / phi_scale, case.param.cell.v_h / phi_scale

	# 沿用现有系数函数；生产路径逐支路消费与 SPMe_element 相同的状态、温度。
	mesh_ne, mesh_pe, mesh_sp = case.mesh["negative electrode"], case.mesh["positive electrode"], case.mesh["separator"]
	length(yt_chem) == ne || throw(DimensionMismatch("local chemistry branch count"))
	coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
	local_concentrations = Dict{String,Union{Array{Float64},Float64}}()
	for e in 1:ne
		yt_e = yt_chem[e]
		local_concentrations["negative particle surface lithium concentration"] = yt_e[case.index["negative particle surface lithium concentration"]]
		local_concentrations["positive particle surface lithium concentration"] = yt_e[case.index["positive particle surface lithium concentration"]]
		ce_n = yt_e[case.index["electrolyte lithium concentration in negative electrode"]]
		ce_p = yt_e[case.index["electrolyte lithium concentration in positive electrode"]]
		ce_sp = yt_e[case.index["electrolyte lithium concentration in separator"]]
		local_concentrations["electrolyte lithium concentration at negative electrode Gauss point"] = sum(mesh_ne.gs.Ni .* ce_n[mesh_ne.element[mesh_ne.gs.ele, :]], dims = 2)
		local_concentrations["electrolyte lithium concentration at positive electrode Gauss point"] = sum(mesh_pe.gs.Ni .* ce_p[mesh_pe.element[mesh_pe.gs.ele, :]], dims = 2)
		local_concentrations["electrolyte lithium concentration at separator Gauss point"] = sum(mesh_sp.gs.Ni .* ce_sp[mesh_sp.element[mesh_sp.gs.ele, :]], dims = 2)
		prefactors = compute_prefactors(local_concentrations, case.param, mesh_ne, mesh_pe, mesh_sp)
		coeffs[e] = compute_element_coefficients(e, Te_prev[e], case.param, prefactors, case.param.cell.T0)
	end

	deactivated_mask = .!area_active
	inactive_reason = zeros(Float64, ne)

	# 面积全失活：带载显式失败；静置走诊断分支（status=4）
	if !any(area_active)
		if I_total != 0.0
			error("面积全失活，无法承载外加电流：I_total=$(I_total)（无量纲，×I_typ=$(case.param.scale.I_typ) A = $(I_total * case.param.scale.I_typ) A），t=$(t)（无量纲，×t0=$(case.param.scale.t0) s）")
		end
		# R3 返工：诊断值发布前的合法性门——非法温度/非有限 OCV 直接抛错，不发布 NaN 诊断电压
		(all(isfinite, Te_prev) && minimum(Te_prev) > 0.0) || error(
			"solve_branch_currents: 全失活静置诊断拒绝非法温度（NaN/非正绝对温度），Te_prev 最小值=$(if isempty(Te_prev)
				NaN
			else
				minimum(Te_prev)
			end)")
		I_e = zeros(Float64, ne)
		V = mean(c.C1 for c in coeffs)   # OCV 算术均值诊断值，非求解所得公共端电压
		isfinite(V) || error(
			"solve_branch_currents: 全失活静置诊断 OCV 均值非有限（NaN/Inf），拒绝发布诊断电压")
		ci = detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)
		voltage_in_bounds, cutoff_type_global, _, _ = check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, zeros(ne), I_e)

		variables["thermal2D branch weight"] = zeros(Float64, ne)
		variables["thermal2D effective area"] = A_eff
		variables["thermal2D element current"] = I_e
		variables["thermal2D element current A"] = zeros(Float64, ne)
		variables["thermal2D element branch current [A]"] = zeros(Float64, ne)
		variables["thermal2D common voltage"] = V
		variables["thermal2D common voltage V"] = V * phi_scale
		variables["thermal2D Vsolve status"] = 4.0
		variables["thermal2D Vsolve iters"] = 0.0
		variables["thermal2D Vsolve converged"] = 0.0
		variables["thermal2D n_active_elements"] = 0.0
		variables["thermal2D n_cutoff_elements"] = Float64(ci.n_cutoff)
		variables["thermal2D n_deactivated_elements"] = Float64(ne)
		variables["thermal2D active_mask"] = zeros(Float64, ne)
		variables["thermal2D deactivated_mask"] = Float64.(deactivated_mask)
		fill!(inactive_reason, 2.0)
		variables["thermal2D inactive_reason"] = inactive_reason
		variables["thermal2D voltage_in_bounds"] = if voltage_in_bounds
			1.0
		else
			0.0
		end
		variables["thermal2D cutoff_type_global"] = Float64(cutoff_type_global)
		variables["thermal2D element OCV"] = ci.all_ocv
		variables["thermal2D cutoff_elements"] = Float64.(ci.cutoff_elements)
		variables["thermal2D cutoff_ocv"] = ci.cutoff_ocv
		variables["thermal2D cutoff_type"] = Float64.(ci.cutoff_type)
		variables["thermal2D nearest_cutoff_element"] = Float64(ci.nearest_element)
		variables["thermal2D nearest_cutoff_ocv"] = ci.nearest_ocv
		variables["thermal2D margin_to_cutoff"] = ci.margin
		return variables, I_e, V
	end

	w = A_eff ./ sum(A_eff)
	# 任务63 B2：暴露真实消费权重/有效面积（ON=A_eff/ΣA_eff，OFF=几何权重），供诊断观察点与结果导出
	variables["thermal2D branch weight"] = w
	variables["thermal2D effective area"] = A_eff
	I_density_target = I_total / gamma
	isfinite(I_density_target) || error("unrepresentable branch density target for effective area factor $gamma")

	# 截止检测 + 面积活跃合并
	ci = detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)
	active_mask = ci.active_mask .& area_active
	active_idx = findall(active_mask)

	# 初始化电流（非活跃单元置零）
	I_e = initialize_currents(ne, w, I_density_target, x_prev; active_mask=active_mask)

	# 初始电压（活跃单元平均，无活跃单元则用 OCV 均值——既有全截止路径，保留原语义）
	V = if isempty(active_idx)
		mean([coeffs[e].C1 for e in 1:ne])
	else
		mean([branch_voltage(coeffs[e], I_e[e]) for e in active_idx])
	end

	# 牛顿迭代（统一路径，active_mask 为空时直接跳过）
	converged, last_iter = true, 0
	if !isempty(active_idx)
		V, converged, last_iter = newton_iteration(I_e, V, ne, w, I_density_target, coeffs; active_mask=active_mask)
		converged = converged && abs(sum(a .* I_e) - I_total) <= 1e-10
	end
	converged || error("branch-current Newton solver failed to converge after $last_iter iterations")

	# 边界检查
	voltage_in_bounds, cutoff_type_global, _, _ = check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e)

	# 写入 variables
	variables["thermal2D element current"] = I_e
	variables["thermal2D element current A"] = case.param.scale.I_typ .* I_e # legacy equivalent-current alias
	variables["thermal2D element branch current [A]"] = case.param.scale.I_typ .* a .* I_e
	variables["thermal2D common voltage"] = V
	variables["thermal2D common voltage V"] = V * phi_scale
	variables["thermal2D Vsolve status"] = 3.0
	variables["thermal2D Vsolve iters"] = Float64(last_iter)
	variables["thermal2D Vsolve converged"] = 1.0
	variables["thermal2D n_active_elements"] = Float64(sum(active_mask))
	variables["thermal2D n_cutoff_elements"] = Float64(ci.n_cutoff)
	variables["thermal2D n_deactivated_elements"] = Float64(sum(deactivated_mask))
	variables["thermal2D active_mask"] = Float64.(active_mask)
	variables["thermal2D deactivated_mask"] = Float64.(deactivated_mask)
	for e in 1:ne
		inactive_reason[e] = if !area_active[e]
			2.0
		else
			if !active_mask[e]
				1.0
			else
				0.0
			end
		end
	end
	variables["thermal2D inactive_reason"] = inactive_reason
	variables["thermal2D voltage_in_bounds"] = if voltage_in_bounds
		1.0
	else
		0.0
	end
	variables["thermal2D cutoff_type_global"] = Float64(cutoff_type_global)
	variables["thermal2D element OCV"] = ci.all_ocv
	variables["thermal2D cutoff_elements"] = Float64.(ci.cutoff_elements)
	variables["thermal2D cutoff_ocv"] = ci.cutoff_ocv
	variables["thermal2D cutoff_type"] = Float64.(ci.cutoff_type)
	variables["thermal2D nearest_cutoff_element"] = Float64(ci.nearest_element)
	variables["thermal2D nearest_cutoff_ocv"] = ci.nearest_ocv
	variables["thermal2D margin_to_cutoff"] = ci.margin
	return variables, I_e, V
end
