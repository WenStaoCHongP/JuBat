# PCC/NCC 平面应力 J2 塑性（Batch 3，spec §3.3；Simo–Hughes 一致返回映射）
#
# 单位契约：调用方传入的 C、σ_y、H 与 e_mech 同单位（生产路径均为 σ_czm/L 归一系；
# 本构函数本身与单位无关，unit 测试直接用物理量纲）。

"""
    PlasticState

PCC/NCC 高斯点塑性状态（spec §4.2）。
`eps_p` 为塑性应变（工程剪切约定，与 E_vec = [E11, E22, γ12] 同构）；`kappa` 为等效塑性应变。
"""
mutable struct PlasticState
    eps_p::NTuple{3, Float64}
    kappa::Float64
end
PlasticState() = PlasticState((0.0, 0.0, 0.0), 0.0)

"""
    foil_params_of(param, mt) -> (E, ν, σ_y, H)

物理箔参数（D-B3-0）：`czm_j2_plasticity=true` 时 PCC/NCC 采用金属箔本构（σ_czm 归一系，
不经 `moduli_of` 的 E_coat 双缩放链）。其余材料类型返回 `nothing`。
"""
function foil_params_of(param, mt::Symbol)
    mt === :PCC && return (param.PCC.E_foil, param.PCC.nu_foil, param.PCC.sigma_y, param.PCC.H)
    mt === :NCC && return (param.NCC.E_foil, param.NCC.nu_foil, param.NCC.sigma_y, param.NCC.H)
    return nothing
end

"""
    plane_stress_C(E, ν) -> 3×3 平面应力弹性矩阵（[σ11,σ22,σ12] vs [E11,E22,γ12]）
"""
plane_stress_C(E::Float64, ν::Float64) =
    E / (1.0 - ν^2) * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1.0 - ν) / 2.0]

"""
    qbar(σ) -> σ̄（σ33≡0 的 J2 等效应力：σ̄² = σ11²+σ22²−σ11σ22+3σ12²）
"""
qbar(σ::AbstractVector{<:Real}) =
    sqrt(σ[1]^2 + σ[2]^2 - σ[1] * σ[2] + 3 * σ[3]^2)

const _A_PS = [1.0 -0.5 0.0; -0.5 1.0 0.0; 0.0 0.0 3.0]  # σ̄ = sqrt(σᵀAσ)

"""
    return_mapping_plane_stress(e_mech, C, σ_y, H, eps_p, κ) -> (σ, C_ep, Δeps_p, Δκ)

平面应力一致 J2 返回映射（关联流动、各向同性线性硬化 σ_y+Hκ）。

- 弹性步（f≤0）：`σ = C(e_mech − eps_p)`，`C_ep = C`。
- 塑性步：解 4×4 系统 `x=[σ11,σ22,σ12,Δγ]`：
  `R₁․₃ = σ − C(e_mech − eps_p − Δγ n(σ))`，`R₄ = σ̄(σ) − σ_y − H(κ+Δγ)`，
  `n = ∂σ̄/∂σ = Aσ/σ̄`；Newton（3D 径向返回初值，tol 1e-12，≤50 步）。
- 一致切线：收敛系统线性化 `J[dσ;dΔγ] = −∂R/∂e·de`，`∂R₁․₃/∂e = −C` ⟹
  `C_ep = dσ/de = (J⁻¹)[1:3,1:3]·C`（算法一致性切线，有限差分可精确复现）。
"""
function return_mapping_plane_stress(e_mech, C::Matrix{Float64}, σ_y::Float64,
                                     H::Float64, eps_p::NTuple{3, Float64}, κ::Float64)
    e_trial = [e_mech[1] - eps_p[1], e_mech[2] - eps_p[2], e_mech[3] - eps_p[3]]
    σ = C * e_trial
    f = qbar(σ) - (σ_y + H * κ)
    if f ≤ 0.0
        return σ, copy(C), (0.0, 0.0, 0.0), 0.0
    end

    G = C[3, 3]
    Δγ0 = f / (3.0 * G + H)
    x = [σ[1], σ[2], σ[3], Δγ0]
    J = zeros(4, 4)
    for _ in 1:50
        σ1, σ2, σ3, Δγ = x
        q = qbar([σ1, σ2, σ3])
        q > 1e-300 || error("return_mapping_plane_stress: σ̄→0 in plastic step (invalid state)")
        n = (_A_PS * [σ1, σ2, σ3]) / q
        σhat = C * (e_trial .- Δγ .* n)
        R = [σ1 - σhat[1], σ2 - σhat[2], σ3 - σhat[3], q - σ_y - H * (κ + Δγ)]
        Dn = (_A_PS - n * n') / q
        J[1:3, 1:3] .= (Δγ .* (C * Dn)) + I
        J[1:3, 4] .= C * n
        J[4, 1:3] .= n
        J[4, 4] = -H
        δ = -(J \ R)
        x .+= δ
        norm(R, Inf) < 1e-12 && norm(δ, Inf) < 1e-12 && break
    end
    σn = x[1:3]
    Δγ = x[4]
    Δγ ≥ 0.0 || error("return_mapping_plane_stress: negative plastic multiplier Δγ=$Δγ")
    q = qbar(σn)
    n = (_A_PS * σn) / q
    Δeps_p = (Δγ * n[1], Δγ * n[2], Δγ * n[3])
    # 一致切线：C_ep = (J⁻¹)[1:3,1:3]·C（J 为收敛点雅可比）
    X = J \ Matrix{Float64}(I, 4, 4)
    C_ep = X[1:3, 1:3] * C
    return σn, C_ep, Δeps_p, Δγ
end

"""
    clone_plastic_states(states::Matrix{PlasticState}) -> 深拷贝（试算/回滚用）
"""
clone_plastic_states(states::Matrix{PlasticState}) =
    PlasticState[PlasticState(s.eps_p, s.kappa) for s in states]
