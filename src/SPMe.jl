function SPMe(case::Case, yt, t::Float64; jacobi::String="update", e::Union{Nothing,Int}=nothing, I_e::Union{Nothing,Float64}=nothing, T_e::Union{Nothing,Float64}=nothing, workspace::Union{Nothing,Dict{String,Union{Array{Float64},Float64}}}=nothing)
    thermal_distributed = case.opt.thermalmodel == "distributed2D"
    if thermal_distributed
        variables = SPMe_variables(case, vec(yt), t; I_app=I_e, T_e=T_e, workspace=workspace)
    else
        variables = SPMe_variables(case, yt, t)
    end
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
    if thermal_distributed
        M_np = case.mesh["SPMe M_np"]
        M_pp = case.mesh["SPMe M_pp"]
        if case.opt.mechanicalmodel == "none"
            K_np = case.mesh["SPMe K_np none"]
            K_pp = case.mesh["SPMe K_pp none"]
        else
            mesh_np = case.mesh["negative particle"]
            mesh_pp = case.mesh["positive particle"]
            K_np = electrode_diffusion_stiffness(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
            K_pp = electrode_diffusion_stiffness(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
        end
    # 单 SPMe 缓存已归一化，直接复用；不与分布式单元缓存混用。
    elseif jacobi == "constant" && param.NE.M_d != []
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
    if thermal_distributed
        M_el = case.mesh["SPMe M_el"]
        K_el = electrolyte_diffusion_stiffness(param, mesh_el, mesh_el.nlen, variables)
    else
        M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)
        M_el = M_el .* param.scale.te / case.param_dim.scale.t0
    end
    F = SPMe_BC(case, variables)
    M = blockdiag(M_np, M_pp, M_el)
    K = blockdiag(K_np, K_pp, K_el)
    if thermal_distributed
        variables["element index"] = Float64(e)
    end
    return M, K, F, variables
end

function initialize_spme_matrices!(case::Case)
    # 冻结参数/网格下的只读装配缓存保存在 mesh；消费者不得原位修改。
    param = case.param
    mesh_np = case.mesh["negative particle"]
    mesh_pp = case.mesh["positive particle"]
    mesh_el = case.mesh["electrolyte"]
    case.mesh["SPMe M_np"] = electrode_diffusion_mass(mesh_np, mesh_np.nlen) .* param.scale.ts_n / case.param_dim.scale.t0
    case.mesh["SPMe M_pp"] = electrode_diffusion_mass(mesh_pp, mesh_pp.nlen) .* param.scale.ts_p / case.param_dim.scale.t0
    shape_ne = (length(case.mesh["negative electrode"].gs.ele), 1)
    shape_sp = (length(case.mesh["separator"].gs.ele), 1)
    shape_pe = (length(case.mesh["positive electrode"].gs.ele), 1)
    case.mesh["SPMe M_el"] = electrolyte_diffusion_mass(param, mesh_el, mesh_el.nlen, shape_ne, shape_sp, shape_pe) .* param.scale.te / case.param_dim.scale.t0
    # 预存无应力常数 K，支持 full case 的现有 deepcopy→none 对照；full 路径不消费它。
    case.mesh["SPMe K_np none"] = electrode_diffusion_stiffness(mesh_np, mesh_np.nlen, param.NE.Ds * Arrhenius(param.NE.Eac_D, 1.0))
    case.mesh["SPMe K_pp none"] = electrode_diffusion_stiffness(mesh_pp, mesh_pp.nlen, param.PE.Ds * Arrhenius(param.PE.Eac_D, 1.0))
    return case
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

"""
    SPMe_variables(case, yt, t; I_app=nothing, T_e=nothing, workspace=nothing)

计算 SPMe 变量；distributed2D 模式必须提供 workspace，直接覆写并返回该字典；其他热模式创建完整 StandardVariables 字典。
I_app 和 T_e 分别注入无量纲电流与单元温度；未注入温度时，关闭热模型取初温，否则读取唯一温度状态。
"""
function SPMe_variables(case::Case, yt::AbstractVector{Float64}, t::Float64; I_app::Union{Nothing,Float64}=nothing, T_e::Union{Nothing,Float64}=nothing, workspace::Union{Nothing,Dict{String,Union{Array{Float64},Float64}}}=nothing)
    param = case.param
    thermal_distributed = case.opt.thermalmodel == "distributed2D"
    if thermal_distributed
        variables = workspace
    else
        variables = StandardVariables(case, 1)
    end
    # 允许外部覆盖无量纲电流与温度
    if isnothing(I_app)
        I_app = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    end

    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = - I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]
    for i in ("negative particle lithium concentration",
        "positive particle lithium concentration",
        "negative particle surface lithium concentration",
        "positive particle surface lithium concentration",
        "electrolyte lithium concentration",
        "electrolyte lithium concentration in negative electrode",
        "electrolyte lithium concentration in positive electrode",
        "electrolyte lithium concentration in separator")
        variables[i] = yt[case.index[i]]
    end
    # 温度：单元温度注入优先；未注入时热模型关闭取环境初温、启用取唯一温度状态
    if T_e !== nothing
        T = T_e
    elseif case.opt.thermalmodel == "none"
        T = case.param.cell.T0
    else
        T = only(yt[case.index["temperature"]])
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
    dphi_S = I_app / 3 * (param.NE.thickness / (param.NE.sig * param.NE.eps_s) + param.PE.thickness / (param.PE.sig * param.PE.eps_s))
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
    if !thermal_distributed
        variables["thermal2D element current"] = I_app * (param.scale.I_typ / case.param_dim.cell.I1C)
    end
    return variables
end
