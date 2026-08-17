"""
Materialmatrix.jl - constitutive and gap-conductance models.
"""

# ========================================================================
# Thermal material matrices
# ========================================================================

"""
	thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)

Compute per-Gauss-point capacity weights for jellyroll 2D thermal assembly.

网格已无量纲化，直接使用 wJ。
"""
function thermal_capacity_weights_2d(param::Params, fks::Matrix{Float64}, ele_of_gp::Vector{Int64}, wJ::Vector{Float64})
	ne = size(fks, 1)
	rho_c_e = zeros(Float64, ne)
	@inbounds for e in 1:ne
		rho_c_e[e] = fks[e, 1] * (param.NE.rho * param.NE.heat_Q) + fks[e, 2] * (param.SP.rho * param.SP.heat_Q) + fks[e, 3] * (param.PE.rho * param.PE.heat_Q) + fks[e, 4] * (param.PCC.rho * param.PCC.heat_Q) + fks[e, 5] * (param.NCC.rho * param.NCC.heat_Q)
	end
	return rho_c_e[ele_of_gp] .* wJ
end

"""
	thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy)

Compute per-Gauss-point anisotropic conductivity components (k_xx, k_xy, k_yy).
"""
function thermal_anisotropic_conductivity_2d(param::Params,fks::Matrix{Float64},ele_of_gp::Vector{Int64},gx::Vector{Float64},gy::Vector{Float64})
	ne = size(fks, 1)
	lam_r_e = zeros(Float64, ne)
	lam_t_e = zeros(Float64, ne)
	@inbounds for e in 1:ne
		f = @view fks[e, :]
		# 径向热导率：串联热阻模型，要求各层热导率 > 0
		denom = f[1] / param.NE.lambda + f[2] / param.SP.lambda + f[3] / param.PE.lambda + f[4] / param.PCC.lambda + f[5] / param.NCC.lambda
		denom > 0 || error("thermal_anisotropic_conductivity_2d: element $e has zero radial thermal conductivity (check lambda values)")
		lam_r_e[e] = 1.0 / denom
		lam_t_e[e] = f[1] * param.NE.lambda + f[2] * param.SP.lambda + f[3] * param.PE.lambda + f[4] * param.PCC.lambda + f[5] * param.NCC.lambda
	end

	ngs = length(ele_of_gp)
	k_xx = zeros(Float64, ngs)
	k_xy = zeros(Float64, ngs)
	k_yy = zeros(Float64, ngs)
	@inbounds for g in 1:ngs
		theta = atan(gy[g], gx[g])
		c, s = cos(theta), sin(theta)
		lr, lt = lam_r_e[ele_of_gp[g]], lam_t_e[ele_of_gp[g]]
		k_xx[g] = lr * c * c + lt * s * s
		k_xy[g] = (lt - lr) * s * c
		k_yy[g] = lr * s * s + lt * c * c
	end

	return k_xx, k_xy, k_yy
end

# ========================================================================
# CZM bilinear constitutive model
# ========================================================================

"""
	bilinear_traction_state(δ_n, δ_t, damage_state, params)

Compute bilinear traction and return updated damage state.
"""
function bilinear_traction_state(δ_n::Float64, δ_t::Float64, damage_state::DamageState, params::CzmInterfaceParams; visc_beta::Float64=1.0)
	K_n = params.K_n
	K_t = params.K_t
	δ_0_n = params.δ_0_n
	δ_c_n = params.δ_c_n
	δ_0_t = params.δ_0_t
	δ_c_t = params.δ_c_t
	η = params.η

	new_state = DamageState()
	new_state.D = damage_state.D
	new_state.D_visc = damage_state.D_visc
	new_state.δ_max_n = damage_state.δ_max_n
	new_state.δ_max_t = damage_state.δ_max_t
	new_state.δ_max_eff = damage_state.δ_max_eff
	new_state.fractured = damage_state.fractured
	new_state.accumulated_damage = damage_state.accumulated_damage

	czm_model = params.czm_model

	if damage_state.fractured
		new_state.D = 1.0
		new_state.D_visc = 1.0
		new_state.fractured = true
		if czm_model == "model1"
			T_n = 0.0
			T_t = K_t * δ_t
		else
			T_n = 0.0
			T_t = 0.0
		end
		return T_n, T_t, 1.0, new_state
	end
	δ_n_pos = max(0.0, δ_n)
	if czm_model == "model1"
		δ_eff = δ_n_pos
		δ_0_eff = δ_0_n
		δ_c_eff = δ_c_n
	else
		δ_eff = sqrt(δ_n_pos^2 + δ_t^2)
		if δ_eff > 1e-15
			β = abs(δ_t) / δ_eff
			δ_0_eff = sqrt(δ_0_n^2 + (δ_0_t^2 - δ_0_n^2) * β^η)
			δ_c_eff = sqrt(δ_c_n^2 + (δ_c_t^2 - δ_c_n^2) * β^η)
		else
			δ_0_eff = δ_0_n
			δ_c_eff = δ_c_n
		end
	end

	δ_max_hist = damage_state.δ_max_eff
	D_eq = damage_state.D

	if δ_eff > δ_max_hist
		if δ_eff <= δ_0_eff
			D_eq = 0.0
		elseif δ_eff >= δ_c_eff
			D_eq = 1.0
		else
			D_eq = δ_c_eff * (δ_eff - δ_0_eff) / (δ_eff * (δ_c_eff - δ_0_eff))
		end

		new_state.δ_max_eff = δ_eff
		new_state.δ_max_n = max(new_state.δ_max_n, δ_n_pos)
		if czm_model != "model1"
			new_state.δ_max_t = max(new_state.δ_max_t, abs(δ_t))
		end
		new_state.D = D_eq
		new_state.accumulated_damage = max(new_state.accumulated_damage, D_eq)
		if D_eq >= 1.0 - 1e-10
			new_state.fractured = true
		end
	end

	# Viscous damage: D_visc = D_visc_committed + visc_beta * (D_eq - D_visc_committed)
	D_visc = damage_state.D_visc + visc_beta * (D_eq - damage_state.D_visc)
	D_visc = max(damage_state.D_visc, D_visc)  # monotonicity
	new_state.D_visc = D_visc

	# Traction uses D_visc (not D_eq)
	if δ_n >= 0
		T_n = (1.0 - D_visc) * K_n * δ_n
	else
		T_n = K_n * δ_n
	end

	if czm_model == "model1"
		T_t = K_t * δ_t
	else
		T_t = (1.0 - D_visc) * K_t * δ_t
	end

    return T_n, T_t, D_eq, new_state
end

function bilinear_traction(δ_n::Float64, δ_t::Float64, damage_state::DamageState, params::CzmInterfaceParams; update::Bool=true, visc_beta::Float64=1.0)
	T_n, T_t, D, new_state = bilinear_traction_state(δ_n, δ_t, damage_state, params; visc_beta=visc_beta)
	if update
		damage_state.D = new_state.D
		damage_state.D_visc = new_state.D_visc
		damage_state.δ_max_n = new_state.δ_max_n
		damage_state.δ_max_t = new_state.δ_max_t
		damage_state.δ_max_eff = new_state.δ_max_eff
		damage_state.fractured = new_state.fractured
		damage_state.accumulated_damage = new_state.accumulated_damage
	end
	return T_n, T_t, D
end

"""
	bilinear_tangent(δ_n, δ_t, damage_state, params)

Compute bilinear tangent stiffness matrix.
"""
function bilinear_tangent(δ_n::Float64, δ_t::Float64, damage_state::DamageState, params::CzmInterfaceParams; visc_beta::Float64=1.0)
	K_n = params.K_n
	K_t = params.K_t
	δ_0_n = params.δ_0_n
	δ_c_n = params.δ_c_n
	δ_0_t = params.δ_0_t
	δ_c_t = params.δ_c_t
	η = params.η
	czm_model = params.czm_model
	dT_dδ = zeros(Float64, 2, 2)

	if damage_state.fractured
		dT_dδ[1, 1] = 1e-10 * K_n
		if czm_model == "model1"
			dT_dδ[2, 2] = K_t
		else
			dT_dδ[2, 2] = 1e-10 * K_t
		end
		return dT_dδ
	end

	δ_n_pos = max(0.0, δ_n)
	if czm_model == "model1"
		δ_eff = δ_n_pos
		δ_0_eff = δ_0_n
		δ_c_eff = δ_c_n
	else
		δ_eff = sqrt(δ_n_pos^2 + δ_t^2)
		if δ_eff > 1e-15
			β = abs(δ_t) / δ_eff
			δ_0_eff = sqrt(δ_0_n^2 + (δ_0_t^2 - δ_0_n^2) * β^η)
			δ_c_eff = sqrt(δ_c_n^2 + (δ_c_t^2 - δ_c_n^2) * β^η)
		else
			δ_0_eff = δ_0_n
			δ_c_eff = δ_c_n
		end
	end

	# Compute D_eq and D_visc (same logic as bilinear_traction_state for consistency)
	δ_max_hist = damage_state.δ_max_eff
	is_loading = (δ_eff > δ_max_hist - 1e-15)

	D_eq = damage_state.D
	if is_loading && δ_eff > δ_0_eff && δ_eff < δ_c_eff
		D_eq = δ_c_eff * (δ_eff - δ_0_eff) / (δ_eff * (δ_c_eff - δ_0_eff))
	elseif is_loading && δ_eff >= δ_c_eff
		D_eq = 1.0
	elseif is_loading && δ_eff <= δ_0_eff
		D_eq = 0.0
	end
	D_visc = damage_state.D_visc + visc_beta * (D_eq - damage_state.D_visc)
	D_visc = max(damage_state.D_visc, D_visc)  # monotonicity

	if czm_model == "model1"
		dT_dδ[2, 2] = K_t

		if δ_eff <= δ_0_eff || !is_loading
			if δ_n >= 0
				dT_dδ[1, 1] = (1.0 - D_visc) * K_n
			else
				dT_dδ[1, 1] = K_n
			end
		elseif δ_eff >= δ_c_eff
			dT_dδ[1, 1] = 1e-10 * K_n
		else
			if δ_n >= 0 && δ_eff > 1e-15
				dD_dδn = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
				# Key: dD/dδ multiplied by visc_beta for consistent linearization
				dT_dδ[1, 1] = (1.0 - D_visc) * K_n - K_n * δ_n * visc_beta * dD_dδn
			else
				dT_dδ[1, 1] = K_n
			end
		end
		dT_dδ[1, 2] = 0.0
		dT_dδ[2, 1] = 0.0
	else
		if δ_eff <= δ_0_eff || !is_loading
			if δ_n >= 0
				dT_dδ[1, 1] = (1.0 - D_visc) * K_n
			else
				dT_dδ[1, 1] = K_n
			end
			dT_dδ[2, 2] = (1.0 - D_visc) * K_t
		elseif δ_eff >= δ_c_eff
			dT_dδ[1, 1] = 1e-10 * K_n
			dT_dδ[2, 2] = 1e-10 * K_t
		else
			dD_dδeff = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
			if δ_n >= 0 && δ_eff > 1e-15
				dδeff_dδn = δ_n_pos / δ_eff
				dδeff_dδt = δ_t / δ_eff
				# Key: dD/dδ multiplied by visc_beta for consistent linearization
				dT_dδ[1, 1] = (1.0 - D_visc) * K_n - K_n * δ_n * visc_beta * dD_dδeff * dδeff_dδn
				dT_dδ[1, 2] = -K_n * δ_n * visc_beta * dD_dδeff * dδeff_dδt
				dT_dδ[2, 1] = -K_t * δ_t * visc_beta * dD_dδeff * dδeff_dδn
				dT_dδ[2, 2] = (1.0 - D_visc) * K_t - K_t * δ_t * visc_beta * dD_dδeff * dδeff_dδt
			else
				dT_dδ[1, 1] = K_n
				dT_dδ[2, 2] = (1.0 - D_visc) * K_t
			end
		end
	end

	return dT_dδ
end

"""
	update_damage(damage_states, separations, params)

Batch update of damage states.
"""
function update_damage(damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, params::CzmInterfaceParams; visc_beta::Float64=1.0)
	n = length(damage_states)
	@assert length(separations) == n "Mismatch in array lengths"

	new_states = Vector{DamageState}(undef, n)
	for i in 1:n
		state = damage_states[i]
		state isa DamageState || error("update_damage: expected DamageState, got $(typeof(state))")
		δ_n, δ_t = separations[i]
		new_state = last(bilinear_traction_state(δ_n, δ_t, state, params; visc_beta=visc_beta))
		new_states[i] = new_state
	end

	return new_states
end

# ========================================================================
# Gap conductance model
# ========================================================================

"""
	compute_gap_conductance(D, δ_n, params) -> h_eff

Compute effective interface conductance using a parallel thermal circuit model.

The heat transfer across the interface has two parallel paths:
  - Solid contact: h_contact = h_c0 * (1 - D)
  - Gap medium:    h_gap    = k_air / (δ + 2βλ_m)

Effective conductance: h_eff = h_contact + h_gap

单位契约（重设计 v2）：δ_n / δ_0_n / δ_c_n 以 scale.δ_czm 归一（分离空间），
而 h_c0 / k_air / lambda_m / threshold 以 scale.L 归一（热模型长度空间）。
入口处将分离量 ÷Λ（= ×δ_czm/L）转换到 L 空间后再运算。
旧方案 δ_czm = L（Λ = 1）时行为不变。
"""
function compute_gap_conductance(D::Float64, δ_n::Float64, params::CzmInterfaceParams)
	h_c0 = params.h_c0
	k_air = params.k_air
	lambda_m = params.lambda_m
	beta = params.beta
	threshold = params.threshold

	# 分离空间（δ_czm 归一）→ 热模型长度空间（L 归一）
	inv_Λ = 1.0 / params.Λ
	delta0 = params.δ_0_n * inv_Λ
	delta_c = params.δ_c_n * inv_Λ
	delta = max(δ_n, 0.0) * inv_Λ
	D_clamped = clamp(D, 0.0, 0.9999)
	two_beta_lambda = 2.0 * beta * lambda_m

	h_eff = if delta < delta0
		h_c0 + k_air / (delta + two_beta_lambda)
	elseif delta < threshold
		h_c0 * (1.0 - D_clamped) + k_air / (delta + two_beta_lambda)
	else
		h_c0 * (1.0 - D_clamped) + k_air / (delta + delta0)
	end

	h_eff > 0 || error("compute_gap_conductance: zero or negative conductance (h_c0=$h_c0, k_air=$k_air, delta=$delta)")
	return h_eff
end

"""
	compute_element_gap_conductance(czm_mesh, elem_idx, params) -> h_eff
"""
function compute_element_gap_conductance(czm_mesh::CohesiveMesh, elem_idx::Int64, params::CzmInterfaceParams)
	state = czm_mesh.damage_states[elem_idx]
	D = state.D
	δ_n = state.δ_max_n
	return compute_gap_conductance(D, δ_n, params)
end

"""
	get_fractured_elements(czm_mesh) -> Vector{Int64}
"""
function get_fractured_elements(czm_mesh::CohesiveMesh)
	fractured = Int64[]
	for (i, state) in enumerate(czm_mesh.damage_states)
		if state.fractured || state.D >= 0.99
			push!(fractured, i)
		end
	end
	return fractured
end

"""
	get_active_elements(czm_mesh, mesh_data) -> Vector{Int64}
"""
function get_active_elements(czm_mesh::CohesiveMesh, mesh_data::MeshGeometry)
	ne = length(mesh_data.element_layer)
	active = ones(Bool, ne)
	fractured_czm = get_fractured_elements(czm_mesh)

	for e in 1:ne
		if !mesh_data.is_inner_layer[e]
			continue
		end
		for czm_idx in get(mesh_data.czm_element_map, e, Int64[])
			if czm_idx in fractured_czm
				active[e] = false
				break
			end
		end
	end

	return findall(active)
end

"""
	compute_all_gap_conductances(czm_mesh, params) -> Vector{Float64}
"""
function compute_all_gap_conductances(czm_mesh::CohesiveMesh, params::CzmInterfaceParams)
	n_czm = czm_mesh.n_cohesive
	h_eff_all = zeros(Float64, n_czm)
	for i in 1:n_czm
		h_eff_all[i] = compute_element_gap_conductance(czm_mesh, i, params)
	end
	return h_eff_all
end

"""
	effective_area_factor(D::Float64, D_threshold::Float64) -> Float64

计算热单元的有效面积比例因子。

当 D ≤ D_threshold 时返回 1.0（无缩减）；
当 D > D_threshold 时线性缩减至 D=1.0 时为 0.0。

公式: factor = (1 - D) / (1 - D_threshold)
"""
function effective_area_factor(D::Float64, D_threshold::Float64)
	D ≤ D_threshold && return 1.0
	return (1.0 - D) / (1.0 - D_threshold)
end
