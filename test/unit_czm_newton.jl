using LinearAlgebra

function unit_czm_newton_step!(czm_mesh, u::Vector{Float64}, param_cache;
        bc_dofs::Vector{Int}, bc_vals::Vector{Float64},
        F_ext::Union{Nothing,Vector{Float64}}=nothing,
        F_thermo_chem::Union{Nothing,Vector{Float64}}=nothing,
        max_iter::Int=80, tol::Float64=1e-8)

    ndof = length(u)
    F_e = F_ext === nothing ? zeros(ndof) : F_ext
    F_tc = F_thermo_chem === nothing ? zeros(ndof) : F_thermo_chem
    damage_states = czm_mesh.damage_states
    # 将位移 BC 硬写入初值，避免首步残差主导
    for (dof, val) in zip(bc_dofs, bc_vals)
        u[dof] = val
    end

    separations = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    tractions = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    converged = false
    R_norm = Inf

    for iter in 1:max_iter
        K, f_int, separations, tractions = JuBat.assemble_coupled_system(
            czm_mesh, u, param_cache; damage_states=damage_states)
        R = F_e + F_tc - f_int
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        R_norm = norm(R)
        if R_norm < tol
            converged = true
            damage_states = JuBat.update_damage_per_interface(
                czm_mesh, damage_states, separations, param_cache)
            czm_mesh.damage_states = damage_states
            break
        end
        K_bc, R_bc = JuBat.apply_bc_czm(K, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
        Δu = K_bc \ R_bc
        any(!isfinite, Δu) && break
        u .+= Δu
        for (dof, val) in zip(bc_dofs, bc_vals)
            u[dof] = val
        end
    end
    return u, separations, tractions, converged, R_norm
end

"""解析双线性（单调加载，无历史损伤时）。δ,δ_0,δ_c 同空间。"""
function analytic_bilinear_T(δ::Float64, K::Float64, δ_0::Float64, δ_c::Float64)
    δ <= 0 && return K * δ          # 压缩罚（与本构一致的简化）
    δ <= δ_0 && return K * δ
    δ >= δ_c && return 0.0
    # 软化：T = σ_max * (δ_c - δ)/(δ_c - δ_0)，σ_max = K*δ_0
    return (K * δ_0) * (δ_c - δ) / (δ_c - δ_0)
end
