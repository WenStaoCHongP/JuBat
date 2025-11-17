function SPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    variables = SPMe_variables(case, yt, t)
    if case.opt.mechanicalmodel == "full"
        variables = Mechanicaloutput(case,variables)
        theta_Mn = variables["negative particle stress coupling diffusion coefficient"][1]
        theta_Mp = variables["positive particle stress coupling diffusion coefficient"][1]
        else
        theta_Mn = 0.0
        theta_Mp = 0.0
    end
    csn_gs = variables["negative particle concentration at gauss point"]
    csp_gs = variables["positive particle concentration at gauss point"]
    param = case.param
    if jacobi == "constant" && param.NE.M_d != [] # no need to update M and K
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
    M_np = M_np .* param.scale.ts_n / case.param_dim.scale.t0
    M_pp = M_pp .* param.scale.ts_p / case.param_dim.scale.t0
    end
    mesh_el = case.mesh["electrolyte"]        
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)   
    M_el = M_el .* param.scale.te / case.param_dim.scale.t0 
    F = SPMe_BC(case, variables)
    M = blockdiag(M_np, M_pp, M_el)
    K = blockdiag(K_np, K_pp, K_el)

    return M, K, F, variables
end


function SPMe_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen, 1)
    flux_np[end] = - j_n * param.NE.rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen, 1)
    flux_pp[end] = - j_p * param.PE.rs^2

    # electrolyte source term
    mesh_el = case.mesh["electrolyte"]
    coeff = mesh_el.gs.weight .* mesh_el.gs.detJ
    v_ne = collect(1:case.opt.Nn * mesh_el.gs.order)
    v_sp = case.opt.Nn * mesh_el.gs.order .+ collect(1:case.opt.Ns * mesh_el.gs.order)
    v_pe = (case.opt.Nn + case.opt.Ns) * mesh_el.gs.order .+ collect(1:case.opt.Np * mesh_el.gs.order)
    coeff[v_ne] .*= (1 - param.EL.tplus) * param.NE.as .* j_n
    coeff[v_sp] .*= 0
    coeff[v_pe] .*= (1 - param.EL.tplus) * param.PE.as .* j_p

    Vi = mesh_el.element[mesh_el.gs.ele,:]
    flux_el = Assemble1D(Vi, mesh_el.gs.Ni, coeff, mesh_el.nlen)
    flux = [flux_np; flux_pp; flux_el]
    return flux
end

function SPMe_variables(case::Case, yt::Array{Float64}, t::Float64; I_app::Union{Nothing,Float64}=nothing, T_e::Union{Nothing,Float64}=nothing)
    param = case.param
    variables = StandardVariables(case, 1)
    # 移除与热模块相关的占位键，避免覆盖调用方设置的热场/热源（例如 T_nodes 由热网格提供）。
    # 这些键若由外部提前写入（如检查脚本或热耦合流程），应该以外部为准。
    if haskey(variables, "T_nodes"); pop!(variables, "T_nodes"); end
    if haskey(variables, "T_prev"); pop!(variables, "T_prev"); end
    if haskey(variables, "heat_source_fields"); pop!(variables, "heat_source_fields"); end
    # 允许外部覆盖无量纲电流与温度
    if isnothing(I_app)
        I_app = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    else
        I_app = Float64(I_app)
    end

    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = - I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]
    var_list = collect(keys(case.index))
    for i in var_list
        variables[i] = yt[case.index[i]]
    end
    if T_e === nothing
        if "temperature" in var_list
            T = yt[case.index["temperature"]]
        else
            T = case.param.cell.T0
        end
    else
        T = T_e
    end
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    ce_n = variables["electrolyte lithium concentration in negative electrode"]
    ce_p = variables["electrolyte lithium concentration in positive electrode"]
    ce_sp = variables["electrolyte lithium concentration in separator"]

    ce_n_gs = sum(mesh_ne.gs.Ni .* ce_n[mesh_ne.element[mesh_ne.gs.ele, :]], dims = 2)
    ce_p_gs = sum(mesh_pe.gs.Ni .* ce_p[mesh_pe.element[mesh_pe.gs.ele, :]], dims = 2)
    ce_sp_gs = sum(mesh_sp.gs.Ni .* ce_sp[mesh_sp.element[mesh_sp.gs.ele, :]], dims = 2)

    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5
    j0_n_av = IntV(j0_n_gs, mesh_ne) / param.NE.thickness
    j0_p_av = IntV(j0_p_gs, mesh_pe) / param.PE.thickness
    eta_n = 2.0 * T * asinh.(j_n / 2.0 / j0_n_av)
    eta_p = 2.0 * T * asinh.(j_p / 2.0 / j0_p_av)

    ## another implementation in pybamm
    dphi_S =  I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)    
    kappa_ne = param.EL.kappa(param.EL.ce0, T) * param.NE.eps ^ param.NE.brugg
    kappa_pe = param.EL.kappa(param.EL.ce0, T) * param.PE.eps ^ param.PE.brugg
    kappa_sp = param.EL.kappa(param.EL.ce0, T) * param.SP.eps ^ param.SP.brugg
    R_EL = param.NE.thickness / kappa_ne / 3.0  + param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe / 3.0
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    dphi_e = 2.0 * T * (1 - param.EL.tplus) * (csp_av - csn_av)/param.EL.ce0 .- I_app * R_EL .- dphi_S

    u_n = param.NE.U(cn_surf) .+ (T .- case.param.cell.T0) .* param.NE.dUdT(cn_surf)
    u_p = param.PE.U(cp_surf) .+ (T .- case.param.cell.T0) .* param.PE.dUdT(cp_surf)   
    V_cell = u_p - u_n + eta_p - eta_n + dphi_e
    variables["negative particle surface lithium concentration"] = cn_surf
    variables["positive particle surface lithium concentration"] = cp_surf
    variables["cell voltage"] = V_cell[1]
    variables["negative electrode exchange current density"] = j0_n_av
    variables["positive electrode exchange current density"] = j0_p_av
    variables["negative electrode interfacial current density"] = j_n
    variables["positive electrode interfacial current density"] = j_p
    variables["negative electrode overpotential"] = eta_n
    variables["positive electrode overpotential"] = eta_p
    variables["negative electrode open circuit potential"] = u_n
    variables["positive electrode open circuit potential"] = u_p
    variables["electrolyte lithium concentration at negative electrode Gauss point"] = ce_n_gs
    variables["electrolyte lithium concentration at positive electrode Gauss point"] = ce_p_gs
    variables["electrolyte lithium concentration at separator Gauss point"] = ce_sp_gs
    variables["time"] = t
    variables["temperature"] = T
    variables["cell current"] = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C
    thermal_distributed = hasproperty(case.opt, :thermal_enabled) && case.opt.thermal_enabled &&
                         hasproperty(case.opt, :thermalmodel) && case.opt.thermalmodel == "distributed2D"
    if !thermal_distributed
        variables["thermal2D element current"] = I_app * (param.scale.I_typ / case.param_dim.cell.I1C)
    end
    return variables
end

"""
    solve_branch_currents_newton!(case, variables, yt, t, I_total, areas, Te_prev, x_prev=nothing)

非线性分流求解器（放在 ThermalDistributed 模块中以与热相关逻辑靠近）。
外层未知量为公共端电压 V；对于给定 V，逐单元用 SPMe 的局部覆写求解 I_e 使 V_cell_e(x_e)=V。
满足约束 Σ I_e*A_e = I_total。函数返回更新后的 variables、I_e 与 V。
"""
function solve_branch_currents_newton(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing)
    # Goal: solve the coupled nonlinear system derived in 求分电流办法.md via Newton with analytical derivatives.
    # Unknowns: element currents I_e and common terminal voltage V satisfying
    #   V_e(I_e, T_e) - V = 0  for every element e
    #   Σ (I_e * A_e) - I_total = 0
    mesh_ok = haskey(case.mesh, "thermal2D")
    if !mesh_ok
        A_tot = sum(areas)
        w = areas ./ A_tot
        I_e = w .* I_total
        # store both nondimensional local currents (i_e) and physical branch currents (A)
        variables["thermal2D element current"] = I_e
        variables["thermal2D element current A"] = case.param_dim.cell.I1C .* w .* I_e
        variables["thermal2D common voltage"] = 0.0
        variables["thermal2D Vsolve status"] = 1.0
        variables["thermal2D Vsolve iters"] = 0.0
        variables["thermal2D Vsolve converged"] = 0.0
        return variables, I_e, 0.0
    end

    ne = length(areas)
    A_global = sum(areas)
    # area fractions (dimensionless) used to weight nondimensional local currents i_e
    w = areas ./ A_global
    # Voltage clamp limits: ensure they are in the same (nondimensional) units as branch_voltage
    # Case may provide normalized limits in case.param.cell.v_l / v_h (already nondim) or
    # physical limits in case.param_dim.cell.v_l / v_h (in V). Convert physical->nondim via phi scale.
    phi_scale = case.param.scale.phi
    # Use dimensional limits from param_dim and convert to nondimensional via phi_scale.
    # This avoids ambiguity where case.param may or may not have v_l/v_h normalized.
    V_MIN = hasproperty(case.param_dim.cell, :v_l) ? case.param_dim.cell.v_l / phi_scale : -Inf
    V_MAX = hasproperty(case.param_dim.cell, :v_h) ? case.param_dim.cell.v_h / phi_scale : Inf

    # (V_MIN, V_MAX) are nondimensional limits derived from param_dim via phi_scale
    param = case.param
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    # Helper to coerce potential arrays to scalars (use first entry by convention).
    scalarize(x) = isa(x, Number) ? Float64(x) : Float64(x[1])
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
    # Pre-compute temperature-independent electrochemical prefactors.
    prefactor_n = IntV(abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5, mesh_ne) / param.NE.thickness
    prefactor_p = IntV(abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5, mesh_pe) / param.PE.thickness
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    u_n_ref = param.NE.U(cn_surf)
    u_p_ref = param.PE.U(cp_surf)
    du_n_dT = param.NE.dUdT(cn_surf)
    du_p_dT = param.PE.dUdT(cp_surf)
    u_n_ref_val = scalarize(u_n_ref)
    u_p_ref_val = scalarize(u_p_ref)
    du_n_dT_val = scalarize(du_n_dT)
    du_p_dT_val = scalarize(du_p_dT)
    T_ref = case.param.cell.T0
    u_n_val(T) = u_n_ref_val + (T - T_ref) * du_n_dT_val
    u_p_val(T) = u_p_ref_val + (T - T_ref) * du_p_dT_val

    c_sigma = (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig) / 3.0

    coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
    for e in 1:ne
        T_e = Te_prev[e]
        # Exchange current densities (temperature dependent via Arrhenius law).
        arr_n = Arrhenius(param.NE.Eac_k, T_e)
        arr_p = Arrhenius(param.PE.Eac_k, T_e)
        j0_n = param.NE.k * arr_n * prefactor_n
        j0_p = param.PE.k * arr_p * prefactor_p
        # Prevent zero exchange current density from causing singular derivatives.
        j0_n = abs(j0_n) < 1e-16 ? 1e-16 : j0_n
        j0_p = abs(j0_p) < 1e-16 ? 1e-16 : j0_p

        kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps ^ param.NE.brugg
        kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps ^ param.PE.brugg
        kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps ^ param.SP.brugg
        kappa_ne = abs(kappa_ne) < 1e-16 ? 1e-16 : kappa_ne
        kappa_pe = abs(kappa_pe) < 1e-16 ? 1e-16 : kappa_pe
        kappa_sp = abs(kappa_sp) < 1e-16 ? 1e-16 : kappa_sp
        R_EL = param.NE.thickness / kappa_ne / 3.0 + param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe / 3.0

        C1 = (u_p_val(T_e) - u_n_val(T_e)) + 2.0 * T_e * (1.0 - param.EL.tplus) * (csp_av - csn_av) / param.EL.ce0
        C2 = 2.0 * T_e
        alpha_p = -1.0 / (2.0 * j0_p * param.PE.as * param.PE.thickness)
        alpha_n =  1.0 / (2.0 * j0_n * param.NE.as * param.NE.thickness)
        C5 = R_EL + c_sigma
        coeffs[e] = (C1=C1, C2=C2, alpha_p=alpha_p, alpha_n=alpha_n, C5=C5)
    end

    function branch_voltage(coeff, I::Float64)
        apI = coeff.alpha_p * I
        anI = coeff.alpha_n * I
        return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
    end

    function branch_dVdI(coeff, I::Float64)
        apI = coeff.alpha_p * I
        anI = coeff.alpha_n * I
        denom_p = sqrt(1.0 + apI * apI)
        denom_n = sqrt(1.0 + anI * anI)
        return coeff.C2 * (coeff.alpha_p / denom_p - coeff.alpha_n / denom_n) - coeff.C5
    end

    # Initial element current guess: previous step if available, else fraction-weighted split (i_e are nondimensional)
    I_e = x_prev !== nothing && length(x_prev) == ne ? copy(x_prev) : (w .* I_total)
    if sum(w .* I_e) != 0.0
        scale = I_total / sum(w .* I_e)
        I_e .*= scale
    else
        I_e .= (w .* I_total)
    end

    V = sum(branch_voltage(coeffs[e], I_e[e]) for e in 1:ne) / ne

    if abs(I_total) <= 1e-14
        I_e .= 0.0
        V = sum(coeffs[e].C1 for e in 1:ne) / ne
        # 检查无量纲电压是否在允许范围内
        if !(V_MIN <= V <= V_MAX)
            V_phys = V * phi_scale
            throw(ErrorException("thermal2D common voltage out of bounds at zero-current: V(nd)=$(V), V(V)=$(V_phys), allowed [$(V_MIN), $(V_MAX)] nd -> [$(V_MIN*phi_scale), $(V_MAX*phi_scale)] V"))
        end
        variables["thermal2D element current"] = I_e
        variables["thermal2D common voltage"] = V
        variables["thermal2D Vsolve status"] = 0.5
        variables["thermal2D Vsolve iters"] = 0.0
        variables["thermal2D Vsolve converged"] = 1.0
        return variables, I_e, V
    end

    tol_V = 1e-8
    tol_I = 1e-10
    max_iters = 25
    converged = false
    last_iter = 0
    F = zeros(Float64, ne)
    dFdI = similar(F)
    I_trial = similar(I_e)

    for iter in 1:max_iters
        last_iter = iter
        for e in 1:ne
            V_e = branch_voltage(coeffs[e], I_e[e])
            F[e] = V_e - V
            dFdI[e] = branch_dVdI(coeffs[e], I_e[e])
            if abs(dFdI[e]) < 1e-12
                dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5
            end
        end
        res_V = maximum(abs.(F))
    # constraint residual: weighted sum of nondimensional local currents minus total nondimensional current
    res_I = sum(w .* I_e) - I_total
        if res_V <= tol_V && abs(res_I) <= tol_I
            converged = true
            variables["thermal2D Vsolve iters"] = float(iter)
            break
        end

        denom = sum(w ./ dFdI)
        if abs(denom) < 1e-12
            break
        end
        num = -res_I + sum(w .* F ./ dFdI)
        ΔV = num / denom
        ΔI = ((-F) .+ ΔV) ./ dFdI

        λ = 1.0
        success = false
        V_trial = V
        for attempt in 1:12
            ok = true
            V_trial = V + λ * ΔV
            if !isfinite(V_trial)
                ok = false
            end
            if ok
                @inbounds for e in 1:ne
                    val = I_e[e] + λ * ΔI[e]
                    if !isfinite(val) || abs(val) > 1e12
                        ok = false
                        break
                    end
                    I_trial[e] = val
                end
            end
            if ok
                success = true
                break
            end
            λ *= 0.5
        end
        if !success
            break
        end

        I_e .= I_trial
        V = V_trial
    end

    if !converged
        I_e .= (areas ./ (A_global > 0 ? A_global : 1.0)) .* I_total
        V = sum(coeffs[e].C1 for e in 1:ne) / ne
        variables["thermal2D Vsolve iters"] = float(last_iter)
    end

    sx = sum(w .* I_e)
    if sx != 0.0
        I_e .*= (I_total / sx)
    end

    # 最终边界检查：无量纲电压必须处于 [V_MIN, V_MAX]
    if !(V_MIN <= V <= V_MAX)
        V_phys = V * phi_scale
        throw(ErrorException("thermal2D common voltage out of bounds: V(nd)=$(V), V(V)=$(V_phys), allowed [$(V_MIN), $(V_MAX)] nd -> [$(V_MIN*phi_scale), $(V_MAX*phi_scale)] V; I_total_nd=$(I_total), sum(w.*I_e)=$(sum(w .* I_e))"))
    end

    variables["thermal2D element current"] = I_e
    # also expose physical branch currents [A] for diagnostics
    variables["thermal2D element current A"] = case.param_dim.cell.I1C .* w .* I_e
    variables["thermal2D common voltage"] = V
    variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5
    variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0

    return variables, I_e, V
end