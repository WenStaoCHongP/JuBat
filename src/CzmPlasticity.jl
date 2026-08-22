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
    foil_params_of(param, mt) -> (σ_y, H)

箔塑性屈服/硬化参数（Batch 3，`/σ_czm` 归一）。弹性模量/泊松比不在此取——
统一走 `moduli_of(param, mt)`（`param.PCC.E` 为 ÷E_coat 归一，须乘
`E_coat/σ_czm` 链才到 σ_czm 系；用户参数修正后 PCC.E/NCC.E 即物理箔模量）。
其余材料类型返回 `nothing`。
"""
function foil_params_of(param, mt::Symbol)
    mt === :PCC && return (param.PCC.sigma_y, param.PCC.H)
    mt === :NCC && return (param.NCC.sigma_y, param.NCC.H)
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

"""
    winding_prestress_field(czm_mesh, param) -> Vector{NTuple{3,Float64}}

卷绕预应力初始应力场 σ₀（Batch 2'，spec §3.7，D-B2'-1；σ_czm/L 归一系）。
逐单元（质心）全局系三分量 (σ_xx, σ_yy, σ_xy)：

1. **卷入张力（等应变/Voigt 分担）**：`ε_w,side = T_side/Ē_side`，
   `Ē_side = Σ(E_k·t_k)/Σt_k`（NE 侧 {NE,NCC}、PE 侧 {PE,PCC}，SP 取两侧均值），
   层 i 卷入环向张力 `σ_t,i = E_i·ε_w,side(i)`；
2. **对数累积压力**（外层缠绕压内层）：`p(r) = f·ln(R_end/r)`，
   `f = F_rev/t_repeat`，`F_rev = T_ne·t_ne_rev + T_pe·t_pe_rev`（每匝每轴长环向力）；
   `σ_θ0,i = σ_t,i − p(r̄_i)`，`σ_r0 = −p(r̄)`。

未设张力（T 全 0）或负张力 → `error`（AGENTS 9.4/9.7）。自平衡重分布由首个平衡求解完成（D-B2'-4）。
"""
function winding_prestress_field(czm_mesh::CohesiveMesh, param)
    T_ne = param.cell.winding_T_ne
    T_pe = param.cell.winding_T_pe
    (T_ne > 0.0 || T_pe > 0.0) || error(
        "winding_prestress_field: 卷绕张力未设置（cell.winding_T_ne/pe 均为 0）。" *
        "开启 czm_winding_prestress 前须给定参数，不默认、不置零（AGENTS 9.4/9.7）。")
    (T_ne >= 0.0 && T_pe >= 0.0) || error(
        "winding_prestress_field: 卷绕张力为负（T_ne=$T_ne, T_pe=$T_pe）物理非法。")

    Ew(mt) = moduli_of(param, mt)[1]
    t_ne = 2 * param.NE.thickness + param.NCC.thickness
    t_pe = 2 * param.PE.thickness + param.PCC.thickness
    E_ne = (2 * Ew(:NE) * param.NE.thickness + Ew(:NCC) * param.NCC.thickness) / t_ne
    E_pe = (2 * Ew(:PE) * param.PE.thickness + Ew(:PCC) * param.PCC.thickness) / t_pe
    ε_ne = T_ne / E_ne
    ε_pe = T_pe / E_pe
    ε_sp = 0.5 * (ε_ne + ε_pe)
    F_rev = T_ne * t_ne + T_pe * t_pe
    f = F_rev / param.cell.layer

    node = czm_mesh.node
    element = czm_mesh.bulk_element
    mt_all = czm_mesh.czm_submesh.material_type
    ne = size(element, 1)
    R_end = maximum(hypot.(node[:, 1], node[:, 2]))
    σ0 = Vector{NTuple{3, Float64}}(undef, ne)
    for e in 1:ne
        mt = mt_all[e]
        ε_w = mt in (:NE, :NCC) ? ε_ne : mt in (:PE, :PCC) ? ε_pe : ε_sp
        E_i = Ew(mt)
        xs = view(node, element[e, :], 1)
        ys = view(node, element[e, :], 2)
        xc = (xs[1] + xs[2] + xs[3] + xs[4]) / 4
        yc = (ys[1] + ys[2] + ys[3] + ys[4]) / 4
        r = hypot(xc, yc)
        pref = f * log(R_end / r)
        σ_θ = E_i * ε_w - pref
        σ_r = -pref
        c, s = xc / r, yc / r
        tx, ty, nx, ny = -s, c, c, s
        σ0[e] = (σ_θ * tx * tx + σ_r * nx * nx,
                 σ_θ * ty * ty + σ_r * ny * ny,
                 σ_θ * tx * ty + σ_r * nx * ny)
    end
    return σ0
end
