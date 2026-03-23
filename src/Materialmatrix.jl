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
	bilinear_traction_state(δ_n, δ_t, damage_state, cohesive_params)

Compute bilinear traction and return updated damage state.
"""
function bilinear_traction_state(δ_n::Float64, δ_t::Float64, damage_state::DamageState, cohesive_params::Cohesive)
	K_n = cohesive_params.K_n
	K_t = cohesive_params.K_t
	δ_0_n = cohesive_params.δ_0_n
	δ_c_n = cohesive_params.δ_c_n
	δ_0_t = cohesive_params.δ_0_t
	δ_c_t = cohesive_params.δ_c_t
	η = cohesive_params.eta

	new_state = DamageState()
	new_state.D = damage_state.D
	new_state.δ_max_n = damage_state.δ_max_n
	new_state.δ_max_t = damage_state.δ_max_t
	new_state.δ_max_eff = damage_state.δ_max_eff
	new_state.fractured = damage_state.fractured
	new_state.accumulated_damage = damage_state.accumulated_damage

	if damage_state.fractured
		new_state.D = 1.0
		new_state.fractured = true
		return 0.0, 0.0, 1.0, new_state
	end

	czm_model = cohesive_params.czm_model
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
	D = damage_state.D

	if δ_eff > δ_max_hist
		if δ_eff <= δ_0_eff
			D = 0.0
		elseif δ_eff >= δ_c_eff
			D = 1.0
		else
			D = δ_c_eff * (δ_eff - δ_0_eff) / (δ_eff * (δ_c_eff - δ_0_eff))
		end

		new_state.δ_max_eff = δ_eff
		new_state.δ_max_n = max(new_state.δ_max_n, δ_n_pos)
		if czm_model != "model1"
			new_state.δ_max_t = max(new_state.δ_max_t, abs(δ_t))
		end
		new_state.D = D
		new_state.accumulated_damage = max(new_state.accumulated_damage, D)
		if D >= 1.0 - 1e-10
			new_state.fractured = true
		end
	else
		D = new_state.D
	end

	if δ_n >= 0
		T_n = (1.0 - D) * K_n * δ_n
	else
		T_n = K_n * δ_n
	end

	if czm_model == "model1"
		T_t = K_t * δ_t
	else
		T_t = (1.0 - D) * K_t * δ_t
	end

    return T_n, T_t, D, new_state
end

function bilinear_traction(δ_n::Float64, δ_t::Float64, damage_state::DamageState, cohesive_params::Cohesive; update::Bool=true)
	T_n, T_t, D, new_state = bilinear_traction_state(δ_n, δ_t, damage_state, cohesive_params)
	if update
		damage_state.D = new_state.D
		damage_state.δ_max_n = new_state.δ_max_n
		damage_state.δ_max_t = new_state.δ_max_t
		damage_state.δ_max_eff = new_state.δ_max_eff
		damage_state.fractured = new_state.fractured
		damage_state.accumulated_damage = new_state.accumulated_damage
	end
	return T_n, T_t, D
end

"""
	bilinear_tangent(δ_n, δ_t, damage_state, cohesive_params)

Compute bilinear tangent stiffness matrix.
"""
function bilinear_tangent(δ_n::Float64, δ_t::Float64, damage_state::DamageState, cohesive_params::Cohesive)
	K_n = cohesive_params.K_n
	K_t = cohesive_params.K_t
	δ_0_n = cohesive_params.δ_0_n
	δ_c_n = cohesive_params.δ_c_n
	δ_0_t = cohesive_params.δ_0_t
	δ_c_t = cohesive_params.δ_c_t
	η = cohesive_params.eta
	czm_model = cohesive_params.czm_model
	dT_dδ = zeros(Float64, 2, 2)

	if damage_state.fractured
		dT_dδ[1, 1] = 1e-10 * K_n
		dT_dδ[2, 2] = 1e-10 * K_t
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

	D = damage_state.D
	δ_max_hist = damage_state.δ_max_eff
	is_loading = (δ_eff > δ_max_hist - 1e-15)

	if czm_model == "model1"
		dT_dδ[2, 2] = K_t

		if δ_eff <= δ_0_eff || !is_loading
			if δ_n >= 0
				dT_dδ[1, 1] = (1.0 - D) * K_n
			else
				dT_dδ[1, 1] = K_n
			end
		elseif δ_eff >= δ_c_eff
			dT_dδ[1, 1] = 1e-10 * K_n
		else
			if δ_n >= 0 && δ_eff > 1e-15
				dD_dδn = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
				dT_dδ[1, 1] = (1.0 - D) * K_n - K_n * δ_n * dD_dδn
			else
				dT_dδ[1, 1] = K_n
			end
		end
		dT_dδ[1, 2] = 0.0
		dT_dδ[2, 1] = 0.0
	else
		if δ_eff <= δ_0_eff || !is_loading
			if δ_n >= 0
				dT_dδ[1, 1] = (1.0 - D) * K_n
			else
				dT_dδ[1, 1] = K_n
			end
			dT_dδ[2, 2] = (1.0 - D) * K_t
		elseif δ_eff >= δ_c_eff
			dT_dδ[1, 1] = 1e-10 * K_n
			dT_dδ[2, 2] = 1e-10 * K_t
		else
			dD_dδeff = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
			if δ_n >= 0 && δ_eff > 1e-15
				dδeff_dδn = δ_n_pos / δ_eff
				dδeff_dδt = δ_t / δ_eff
				dT_dδ[1, 1] = (1.0 - D) * K_n - K_n * δ_n * dD_dδeff * dδeff_dδn
				dT_dδ[1, 2] = -K_n * δ_n * dD_dδeff * dδeff_dδt
				dT_dδ[2, 1] = -K_t * δ_t * dD_dδeff * dδeff_dδn
				dT_dδ[2, 2] = (1.0 - D) * K_t - K_t * δ_t * dD_dδeff * dδeff_dδt
			else
				dT_dδ[1, 1] = K_n
				dT_dδ[2, 2] = (1.0 - D) * K_t
			end
		end
	end

	return dT_dδ
end

"""
	update_damage(damage_states, separations, cohesive_params)

Batch update of damage states.
"""
function update_damage(damage_states::AbstractVector{<:AbstractDamageState}, separations::Vector{Tuple{Float64, Float64}}, cohesive_params::Cohesive)
	n = length(damage_states)
	@assert length(separations) == n "Mismatch in array lengths"

	new_states = Vector{DamageState}(undef, n)
	for i in 1:n
		state = damage_states[i]
		state isa DamageState || error("update_damage: expected DamageState, got $(typeof(state))")
		δ_n, δ_t = separations[i]
		new_state = last(bilinear_traction_state(δ_n, δ_t, state, cohesive_params))
		new_states[i] = new_state
	end

	return new_states
end

# ========================================================================
# Gap conductance model
# ========================================================================

"""
	compute_gap_conductance(D, δ_n, cohesive) -> h_eff

Compute effective interface conductance.
"""
function compute_gap_conductance(D::Float64, δ_n::Float64, cohesive)
	h_c0 = cohesive.h_c0
	k_air = cohesive.k_air
	lambda_m = cohesive.lambda_m
	beta = cohesive.beta
	threshold = cohesive.threshold

	delta0 = cohesive.δ_0_n
	delta_c = cohesive.δ_c_n
	delta = max(δ_n, 0.0)
	D_sep = if delta <= delta0
		0.0
	elseif delta < delta_c
		(delta - delta0) / (delta_c - delta0)
	else
		1.0
	end
	D_clamped = clamp(max(D, D_sep), 0.0, 0.9999)
	two_beta_lambda = 2.0 * beta * lambda_m

	denom = if delta < delta0
		h_c0 + k_air / (delta + two_beta_lambda)
	elseif delta < threshold
		h_c0 * (1.0 - D_clamped) + k_air / (delta + two_beta_lambda)
	else
		h_c0 * (1.0 - D_clamped) + k_air / (delta + delta0)
	end

	denom > 0 || error("compute_gap_conductance: zero or negative denominator (h_c0=$h_c0, k_air=$k_air, delta=$delta)")
	return 1.0 / denom
end

"""
	compute_element_gap_conductance(czm_mesh, elem_idx, cohesive) -> h_eff
"""
function compute_element_gap_conductance(czm_mesh::CohesiveMesh, elem_idx::Int64, cohesive)
	state = czm_mesh.damage_states[elem_idx]
	D = state.D
	δ_n = state.δ_max_n
	return compute_gap_conductance(D, δ_n, cohesive)
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
function get_active_elements(czm_mesh::CohesiveMesh, mesh_data)
	ne = mesh_data.ne
	active = ones(Bool, ne)
	fractured_czm = get_fractured_elements(czm_mesh)

	for e in 1:ne
		if !mesh_data.is_inner_layer[e]
			continue
		end
		for czm_idx in mesh_data.czm_element_map[e]
			if czm_idx in fractured_czm
				active[e] = false
				break
			end
		end
	end

	return findall(active)
end

"""
	compute_all_gap_conductances(czm_mesh, cohesive) -> Vector{Float64}
"""
function compute_all_gap_conductances(czm_mesh::CohesiveMesh, cohesive)
	n_czm = czm_mesh.n_cohesive
	h_eff_all = zeros(Float64, n_czm)
	for i in 1:n_czm
		h_eff_all[i] = compute_element_gap_conductance(czm_mesh, i, cohesive)
	end
	return h_eff_all
end
