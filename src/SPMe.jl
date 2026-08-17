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

function SPMe_element(case::Case, yt_e, t::Float64, e::Int; I_e::Float64, T_e::Float64, jacobi::String="update",workspace::Union{Nothing, Dict{String, Union{Array{Float64},Float64}}}=nothing)
    # 1) 使用 workspace 或创建新 variables
    if workspace !== nothing
        variables_e = SPMe_variables!(workspace, case, yt_e, t; I_app=I_e, T_e=T_e)
    else
        variables_e = SPMe_variables(case, vec(yt_e), t; I_app=I_e, T_e=T_e)
    end
    
    # 2) 力学耦合（如果启用）
    if case.opt.mechanicalmodel == "full"
        variables_e = Mechanicaloutput(case, variables_e)
        theta_Mn = variables_e["negative particle stress coupling diffusion coefficient"][1]
        theta_Mp = variables_e["positive particle stress coupling diffusion coefficient"][1]
    else
        theta_Mn = 0.0
        theta_Mp = 0.0
    end
    
    # 3) 提取高斯点浓度（用于扩散系数计算）
    csn_gs = variables_e["negative particle concentration at gauss point"]
    csp_gs = variables_e["positive particle concentration at gauss point"]
    param = case.param
    
    # 4) 粒子扩散矩阵
    # 注意：在多线程环境下，应使用 jacobi="update" 避免竞争 param.NE.M_d
    if jacobi == "constant" && !isempty(param.NE.M_d) && !isempty(param.NE.K_d)
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
    end
    # 时间尺度归一化
    M_np = M_np .* param.scale.ts_n / case.param_dim.scale.t0
    M_pp = M_pp .* param.scale.ts_p / case.param_dim.scale.t0
    
    # 5) 电解液扩散矩阵
    mesh_el = case.mesh["electrolyte"]
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables_e)
    M_el = M_el .* param.scale.te / case.param_dim.scale.t0
    
    # 6) 边界条件（源项）
    F_e = SPMe_BC(case, variables_e)
    
    # 7) 装配该单元的局部系统矩阵
    M_e = blockdiag(M_np, M_pp, M_el)
    K_e = blockdiag(K_np, K_pp, K_el)
    
    # 8) 可选：添加单元编号到 variables_e（用于调试）
    variables_e["element index"] = Float64(e)
    
    return M_e, K_e, F_e, variables_e
end

"""
    SPMe_variables!(ws, case, yt, t; I_app, T_e)

SPMe_variables 的原位变体：直接写入预分配 workspace，不创建新 Dict。
计算逻辑与 SPMe_variables 完全一致，仅省去 StandardVariables 分配。
"""
function SPMe_variables!(ws::Dict{String, Union{Array{Float64},Float64}},case::Case, yt, t::Float64;I_app::Union{Nothing,Float64}=nothing,T_e::Union{Nothing,Float64}=nothing)
    param = case.param

    if isnothing(I_app)
        I_app = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    else
        I_app = Float64(I_app)
    end

    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = -I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]

    # 状态提取：仅覆写 ws 中已有的键（跳过 thermal2D 键）
    # 使用替换赋值（与原始 SPMe_variables 行为一致），保留原始类型
    var_list = collect(keys(case.index))
    if T_e !== nothing
        var_list = filter(k -> k != "temperature", var_list)
    end
    for i in var_list
        if haskey(ws, i)
            ws[i] = yt[case.index[i]]
        end
    end

    T = representative_temperature(case, yt; supplied=T_e)

    # ── 以下计算逻辑与 SPMe_variables 逐行一致 ──
    cn_surf = ws["negative particle surface lithium concentration"]
    cp_surf = ws["positive particle surface lithium concentration"]
    ce_n = ws["electrolyte lithium concentration in negative electrode"]
    ce_p = ws["electrolyte lithium concentration in positive electrode"]
    ce_sp = ws["electrolyte lithium concentration in separator"]

    ce_n_gs = sum(mesh_ne.gs.Ni .* ce_n[mesh_ne.element[mesh_ne.gs.ele, :]], dims = 2)
    ce_p_gs = sum(mesh_pe.gs.Ni .* ce_p[mesh_pe.element[mesh_pe.gs.ele, :]], dims = 2)
    ce_sp_gs = sum(mesh_sp.gs.Ni .* ce_sp[mesh_sp.element[mesh_sp.gs.ele, :]], dims = 2)

    j0_n_gs =  param.NE.k * Arrhenius(param.NE.Eac_k, T) .* abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5
    j0_p_gs =  param.PE.k * Arrhenius(param.PE.Eac_k, T) .* abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5
    j0_n_av = IntV(j0_n_gs, mesh_ne) / param.NE.thickness
    j0_p_av = IntV(j0_p_gs, mesh_pe) / param.PE.thickness
    eta_n = 2.0 * T * asinh.(j_n / 2.0 / j0_n_av)
    eta_p = 2.0 * T * asinh.(j_p / 2.0 / j0_p_av)

    dphi_S = I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)
    kappa_ne_gs = param.EL.kappa(ce_n_gs, T) * param.NE.eps ^ param.NE.brugg
    kappa_pe_gs = param.EL.kappa(ce_p_gs, T) * param.PE.eps ^ param.PE.brugg
    kappa_sp_gs = param.EL.kappa(ce_sp_gs, T) * param.SP.eps ^ param.SP.brugg
    kappa_ne_av = IntV(kappa_ne_gs, mesh_ne) / param.NE.thickness
    kappa_pe_av = IntV(kappa_pe_gs, mesh_pe) / param.PE.thickness
    kappa_sp_av = IntV(kappa_sp_gs, mesh_sp) / param.SP.thickness
    R_EL = param.NE.thickness / (3.0 * kappa_ne_av) + param.SP.thickness / kappa_sp_av + param.PE.thickness / (3.0 * kappa_pe_av)
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    dphi_e = 2.0 * T * (1 - param.EL.tplus) * (csp_av - csn_av) / param.EL.ce0 - I_app * R_EL - dphi_S

    u_n = param.NE.U(cn_surf) .+ (T .- case.param.cell.T0) .* param.NE.dUdT(cn_surf)
    u_p = param.PE.U(cp_surf) .+ (T .- case.param.cell.T0) .* param.PE.dUdT(cp_surf)
    V_cell = u_p - u_n + eta_p - eta_n + dphi_e

    # ── 写入 workspace ──
    ws["cell voltage"] = V_cell[1]
    ws["negative electrode exchange current density"] = j0_n_av
    ws["positive electrode exchange current density"] = j0_p_av
    ws["negative electrode interfacial current density"] = j_n
    ws["positive electrode interfacial current density"] = j_p
    ws["negative electrode overpotential"] = eta_n
    ws["positive electrode overpotential"] = eta_p
    ws["negative electrode open circuit potential"] = u_n
    ws["positive electrode open circuit potential"] = u_p
    ws["electrolyte lithium concentration at negative electrode Gauss point"] = ce_n_gs
    ws["electrolyte lithium concentration at positive electrode Gauss point"] = ce_p_gs
    ws["electrolyte lithium concentration at separator Gauss point"] = ce_sp_gs
    ws["time"] = t
    ws["temperature"] = T
    ws["cell current"] = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C

    return ws
end

function SPMe_BC(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    flux_np = zeros(Float64, case.mesh["negative particle"].nlen)
    flux_np[end] = - j_n * param.NE.rs^2
    flux_pp = zeros(Float64, case.mesh["positive particle"].nlen)
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
    flux = vcat(flux_np, flux_pp, flux_el)
    return flux
end

function SPMe_variables(case::Case, yt::AbstractVector{Float64}, t::Float64; I_app::Union{Nothing,Float64}=nothing, T_e::Union{Nothing,Float64}=nothing)
    param = case.param
    variables = StandardVariables(case, 1)
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
    if T_e !== nothing
        var_list = filter(k -> k != "temperature", var_list)
    end
    for i in var_list
        variables[i] = yt[case.index[i]]
    end
    T = representative_temperature(case, yt; supplied=T_e)
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
    dphi_S = I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)    
    kappa_ne_gs = param.EL.kappa(ce_n_gs, T) * param.NE.eps ^ param.NE.brugg
    kappa_pe_gs = param.EL.kappa(ce_p_gs, T) * param.PE.eps ^ param.PE.brugg
    kappa_sp_gs = param.EL.kappa(ce_sp_gs, T) * param.SP.eps ^ param.SP.brugg
    kappa_ne_av = IntV(kappa_ne_gs, mesh_ne) / param.NE.thickness
    kappa_pe_av = IntV(kappa_pe_gs, mesh_pe) / param.PE.thickness
    kappa_sp_av = IntV(kappa_sp_gs, mesh_sp) / param.SP.thickness
    R_EL = param.NE.thickness / (3.0 * kappa_ne_av) + param.SP.thickness / kappa_sp_av + param.PE.thickness / (3.0 * kappa_pe_av)
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    dphi_e = 2.0 * T * (1 - param.EL.tplus) * (csp_av - csn_av) / param.EL.ce0 - I_app * R_EL - dphi_S

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
    thermal_distributed = case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
    if !thermal_distributed
        variables["thermal2D element current"] = I_app * (param.scale.I_typ / case.param_dim.cell.I1C)
    end
    return variables
end
