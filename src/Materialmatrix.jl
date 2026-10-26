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
		k_xy[g] = (lr - lt) * s * c
		k_yy[g] = lr * s * s + lt * c * c
	end

	return k_xx, k_xy, k_yy
end

# ========================================================================
# CZM bilinear constitutive model
# ========================================================================

# Opening interpolation (not an energy BK law). This same mapping supplies
# the traction history and its full directional derivative.
function bilinear_opening_envelope(δ_n::Float64, δ_t::Float64,
        ip::CurrentCollector, czm_model::String)
    δn_pos = max(δ_n, 0.0)
    if czm_model == "model1"
        return (e=δn_pos, a=ip.δ_0, c=ip.δ_c,
            e_n=Float64(δ_n > 0.0), e_t=0.0,
            a_n=0.0, a_t=0.0, c_n=0.0, c_t=0.0)
    end
    czm_model == "mix" || throw(ArgumentError("unknown CZM model: $czm_model"))
    e = hypot(δn_pos, δ_t)
    if e == 0.0
        return (e=e, a=ip.δ_0, c=ip.δ_c,
            e_n=0.0, e_t=0.0, a_n=0.0, a_t=0.0, c_n=0.0, c_t=0.0)
    end
    β = abs(δ_t) / e
    a = sqrt(ip.δ_0^2 + (ip.δ_0_t^2 - ip.δ_0^2) * β^ip.eta)
    c = sqrt(ip.δ_c^2 + (ip.δ_c_t^2 - ip.δ_c^2) * β^ip.eta)
    e_n = δn_pos / e
    e_t = δ_t / e
    a_n = a_t = c_n = c_t = 0.0
    if δ_t != 0.0
        β_n = -β * e_n / e
        β_t = sign(δ_t) / e - β * e_t / e
        a_β = (ip.δ_0_t^2 - ip.δ_0^2) * ip.eta * β^(ip.eta - 1.0) / (2a)
        c_β = (ip.δ_c_t^2 - ip.δ_c^2) * ip.eta * β^(ip.eta - 1.0) / (2c)
        a_n, a_t = a_β * β_n, a_β * β_t
        c_n, c_t = c_β * β_n, c_β * β_t
    end
    return (; e, a, c, e_n, e_t, a_n, a_t, c_n, c_t)
end

function bilinear_envelope_damage(envelope)
    e, a, c = envelope.e, envelope.a, envelope.c
    return if e <= a
        0.0
    elseif e >= c
        1.0
    else
        c * (e - a) / (e * (c - a))
    end
end

"""Bilinear cohesive parameters from one frozen interface and its selected model."""
struct BilinearCohesiveLaw
    ip::CurrentCollector
    model::String
end

"""Step coefficient from real mechanical dt/τ at the solver boundary:
visc_beta=1 with viscosity off and 0 for frozen evaluation; no time is advanced here.
history_loading selects only the tangent side at a history corner."""
struct CohesiveStepData
    visc_beta::Float64
    history_loading::Bool
end

"""
    cohesive_response(law, gap, committed_state, step_data)
        -> (traction, tangent, trial_state)

`gap=(δ_n,δ_t)` uses the existing dimensionless separation coordinates.
The committed GP state is read-only. Traction and tangent use the same trial
state; only the solver may publish an accepted state.
"""
function cohesive_response(law::BilinearCohesiveLaw, gap::NTuple{2,Float64},
        damage_state::DamageState, step_data::CohesiveStepData)
    δ_n, δ_t = gap
    ip, czm_model = law.ip, law.model
    visc_beta = step_data.visc_beta
	envelope = bilinear_opening_envelope(δ_n, δ_t, ip, czm_model)
	new_state = DamageState()
	new_state.D = damage_state.D
	new_state.D_visc = damage_state.D_visc
	new_state.δ_max_n = damage_state.δ_max_n
	new_state.δ_max_t = damage_state.δ_max_t
	new_state.δ_max_eff = damage_state.δ_max_eff
	new_state.fractured = damage_state.fractured
	new_state.accumulated_damage = damage_state.accumulated_damage
	# Component envelopes record geometry independently of the damage driver.
	δ_n_pos = max(0.0, δ_n)
	new_state.δ_max_n = max(new_state.δ_max_n, δ_n_pos)
	new_state.δ_max_t = max(new_state.δ_max_t, abs(δ_t))

	e, a, c = envelope.e, envelope.a, envelope.c
	d = bilinear_envelope_damage(envelope)
	# Geometry peaks are diagnostics; irreversibility acts on damage itself.
	D_eq = max(damage_state.D, d)
	new_state.δ_max_eff = max(damage_state.δ_max_eff, e)
	new_state.D = D_eq
	new_state.accumulated_damage = max(new_state.accumulated_damage, D_eq)

	# Viscous damage: D_visc = D_visc_committed + visc_beta * (D_eq - D_visc_committed)
	D_visc = visc_beta == 1.0 ? D_eq :
		damage_state.D_visc + visc_beta * (D_eq - damage_state.D_visc)
	# Canonical saturation is part of this same traction map. A near-one D
	# still follows the true softening derivative and is never truncated.
	new_state.fractured = D_eq == 1.0 && D_visc >= 1.0 - 1e-10
	if new_state.fractured
		D_visc = 1.0
	end
	new_state.D_visc = D_visc

	# Traction uses D_visc (not D_eq)
	if δ_n >= 0
		T_n = (1.0 - D_visc) * ip.K_n * δ_n
	else
		T_n = ip.K_n * δ_n
	end

	if czm_model == "model1"
		T_t = ip.K_t * δ_t
	else
		T_t = (1.0 - D_visc) * ip.K_t * δ_t
	end

    tangent = _bilinear_tangent(δ_n, δ_t, damage_state, new_state,
        ip, czm_model, visc_beta; history_loading=step_data.history_loading)
    return (T_n, T_t), tangent, new_state
end

function bilinear_damage_gradient(δ_n::Float64, δ_t::Float64,
        previous::DamageState, trial::DamageState, ip::CurrentCollector,
        czm_model::String, visc_beta::Float64; history_loading::Bool=false)
    zero_grad = (0.0, 0.0)
    (previous.fractured || trial.fractured) && return zero_grad, zero_grad
    envelope = bilinear_opening_envelope(δ_n, δ_t, ip, czm_model)
    (; e, a, c, e_n, e_t, a_n, a_t, c_n, c_t) = envelope
    (a < e < c && (trial.D > previous.D ||
        (history_loading && bilinear_envelope_damage(envelope) == previous.D))) ||
        return zero_grad, zero_grad

    den = c - a
    D_e = c * a / (e^2 * den)
    D_a = c * (e - c) / (e * den^2)
    D_c = -a * (e - a) / (e * den^2)
    q_eq = (D_e * e_n + D_a * a_n + D_c * c_n,
            D_e * e_t + D_a * a_t + D_c * c_t)
    q_visc = (visc_beta * q_eq[1], visc_beta * q_eq[2])
    return q_eq, q_visc
end

function _bilinear_tangent(δ_n::Float64, δ_t::Float64,
        previous::DamageState, trial::DamageState, ip::CurrentCollector,
        czm_model::String, visc_beta::Float64; history_loading::Bool=false)
    dT_dδ = zeros(Float64, 2, 2)
    retained = 1.0 - trial.D_visc
    dT_dδ[1, 1] = if δ_n >= 0.0
        retained * ip.K_n
    else
        ip.K_n
    end
    dT_dδ[2, 2] = if czm_model == "model1"
        ip.K_t
    else
        retained * ip.K_t
    end

    _, q_visc = bilinear_damage_gradient(δ_n, δ_t, previous,
        trial, ip, czm_model, visc_beta; history_loading=history_loading)
    if δ_n > 0.0
        dT_dδ[1, 1] -= ip.K_n * δ_n * q_visc[1]
        dT_dδ[1, 2] -= ip.K_n * δ_n * q_visc[2]
    end
    if czm_model != "model1"
        dT_dδ[2, 1] -= ip.K_t * δ_t * q_visc[1]
        dT_dδ[2, 2] -= ip.K_t * δ_t * q_visc[2]
    end
    return dT_dδ
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

单位契约（重设计 v2）：δ_n / δ_0 / δ_c 以 scale.δ_czm 归一（分离空间），
而 h_c0 / k_air / lambda_m / threshold 以 scale.L 归一（热模型长度空间）。
入口处将分离量 ÷Λ（= ×δ_czm/L）转换到 L 空间后再运算；Λ 使用点内联
`param.scale.L / param.scale.δ_czm`（2026-08-30 重构，不再存字段）。
旧方案 δ_czm = L（Λ = 1）时行为不变。
"""
function compute_gap_conductance(D::Float64, δ_n::Float64, ip::CurrentCollector, param::Params)
	# 分离空间（δ_czm 归一）→ 热模型长度空间（L 归一）
	inv_Λ = param.scale.δ_czm / param.scale.L
	delta0 = ip.δ_0 * inv_Λ
	delta_c = ip.δ_c * inv_Λ
	delta = max(δ_n, 0.0) * inv_Λ
	D_clamped = clamp(D, 0.0, 0.9999)
	two_beta_lambda = 2.0 * ip.beta * ip.lambda_m

	h_eff = if delta < delta0
		ip.h_c0 + ip.k_air / (delta + two_beta_lambda)
	elseif delta < ip.threshold
		ip.h_c0 * (1.0 - D_clamped) + ip.k_air / (delta + two_beta_lambda)
	else
		ip.h_c0 * (1.0 - D_clamped) + ip.k_air / (delta + delta0)
	end

	h_eff > 0 || error("compute_gap_conductance: zero or negative conductance (h_c0=$(ip.h_c0), k_air=$(ip.k_air), delta=$delta)")
	return h_eff
end

"""
	compute_element_gap_conductance(damage_states, elem_idx, ip, param) -> h_eff
"""
function compute_element_gap_conductance(damage_states::AbstractVector{DamageState}, elem_idx::Int64, ip::CurrentCollector, param::Params)
	state = damage_states[elem_idx]
	D = state.D
	δ_n = state.δ_max_n
	return compute_gap_conductance(D, δ_n, ip, param)
end

"""
	get_fractured_elements(damage_states) -> Vector{Int64}

机械断裂统计接口（任务59 P1 后仅保留统计用途，不再接入电流/热源消费）。
"""
function get_fractured_elements(damage_states::AbstractVector{DamageState})
	fractured = Int64[]
	for (i, state) in enumerate(damage_states)
		if state.fractured
			push!(fractured, i)
		end
	end
	return fractured
end

"""
	compute_all_gap_conductances(czm_mesh, params) -> Vector{Float64}
"""
function compute_all_gap_conductances(damage_states::AbstractVector{DamageState}, ip::CurrentCollector, param::Params)
	n_czm = length(damage_states)
	h_eff_all = zeros(Float64, n_czm)
	for i in 1:n_czm
		h_eff_all[i] = compute_element_gap_conductance(damage_states, i, ip, param)
	end
	return h_eff_all
end

"""
	effective_area_factor(D::Float64) -> Float64

无阈值热单元有效面积比例因子（任务59 P1）：factor = 1 - D，D∈[0,1]。
D=0 → 1.0（无损）；D=1 → 0.0（该支路面积失活、电流严格为 0）。
无起始阈值，不保留双参数阈值重载。
"""
function effective_area_factor(D::Float64)
	return 1.0 - D
end
