using LinearAlgebra

"""
单元条带增量 Newton。对位移 BC 使用正确的罚形式：
  K[dof,dof] += P,  R[dof] = P * (val - u[dof])
（生产 apply_bc_czm 对非零 val 写成 P*val，不适配残差 Newton。）

损伤仅在收敛后 commit（与生产 solve_czm_basic_step 一致）；装配内用 trial D_visc。
"""
function unit_czm_newton_step!(czm_mesh, u::Vector{Float64}, param_cache;
        bc_dofs::Vector{Int}, bc_vals::Vector{Float64},
        F_ext::Union{Nothing,Vector{Float64}}=nothing,
        F_thermo_chem::Union{Nothing,Vector{Float64}}=nothing,
        max_iter::Int=80, tol::Float64=1e-8,
        visc_beta::Float64=1.0)

    ndof = length(u)
    F_e = F_ext === nothing ? zeros(ndof) : F_ext
    F_tc = F_thermo_chem === nothing ? zeros(ndof) : F_thermo_chem
    damage_states = czm_mesh.damage_states
    for (dof, val) in zip(bc_dofs, bc_vals)
        u[dof] = val
    end

    separations = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    tractions = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    converged = false
    R_norm = Inf

    for iter in 1:max_iter
        K, f_int, separations, tractions = JuBat.assemble_coupled_system(
            czm_mesh, u, param_cache; damage_states=damage_states,
            visc_beta=visc_beta)
        R = F_e + F_tc - f_int

        dmax = size(K, 1) > 0 ? maximum(abs, diag(K)) : 0.0
        penalty = dmax > 0 ? 1e6 * dmax : 1e12
        K_bc = copy(K)
        R_bc = copy(R)
        for (dof, val) in zip(bc_dofs, bc_vals)
            K_bc[dof, dof] += penalty
            R_bc[dof] = penalty * (val - u[dof])
        end

        R_check = copy(R)
        for (dof, val) in zip(bc_dofs, bc_vals)
            R_check[dof] = val - u[dof]
        end
        R_norm = norm(R_check)
        if R_norm < tol
            converged = true
            damage_states = JuBat.update_damage_per_interface(
                czm_mesh, damage_states, separations, param_cache;
                visc_beta=visc_beta)
            czm_mesh.damage_states = damage_states
            break
        end

        Δu = try
            K_bc \ R_bc
        catch
            break
        end
        any(!isfinite, Δu) && break

        α = 1.0
        accepted = false
        for _ in 1:10
            u_trial = u .+ α .* Δu
            for (dof, val) in zip(bc_dofs, bc_vals)
                u_trial[dof] = val
            end
            _, f_try, sep_try, tr_try = JuBat.assemble_coupled_system(
                czm_mesh, u_trial, param_cache; damage_states=damage_states,
                visc_beta=visc_beta)
            R_try = F_e + F_tc - f_try
            for (dof, val) in zip(bc_dofs, bc_vals)
                R_try[dof] = val - u_trial[dof]
            end
            if norm(R_try) <= (1.0 + 1e-8) * R_norm || α < 1e-4
                u .= u_trial
                separations = sep_try
                tractions = tr_try
                accepted = true
                break
            end
            α *= 0.5
        end
        !accepted && break
    end
    return u, separations, tractions, converged, R_norm
end

"""解析双线性（单调加载，无历史损伤时）。δ,δ_0,δ_c 同空间。"""
function analytic_bilinear_T(δ::Float64, K::Float64, δ_0::Float64, δ_c::Float64)
    δ <= 0 && return K * δ
    δ <= δ_0 && return K * δ
    δ >= δ_c && return 0.0
    return (K * δ_0) * (δ_c - δ) / (δ_c - δ_0)
end
