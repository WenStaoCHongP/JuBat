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

"""
    SPMe_element(case::Case, yt_e::Array{Float64}, t::Float64, e::Int; 
                 I_e::Float64, T_e::Float64, jacobi::String="update")

为单个热单元求解SPMe模型（多SPMe并行架构）。

此函数是 `SPMe` 的单元级版本，用于多SPMe并行模式。每个热单元对应一个独立的
电化学状态向量 yt_e，该状态根据该单元的分电流 I_e 和温度 T_e 独立演化。

# 参数
- `case::Case`: 全局案例对象（包含参数、网格等）
- `yt_e::Array{Float64}`: 该单元的局部电化学状态向量（接受向量或矩阵）
  - 结构: [cn_surf[1:Nrn]; cp_surf[1:Nrp]; ce[1:Nel]]
  - 长度: Nrn + Nrp + Nel（与全局SPMe的yt长度相同）
  - 注：自动转换为向量，兼容 ModelInitialisation 返回的矩阵形式
- `t::Float64`: 当前时间（无量纲，相对 t0）
- `e::Int`: 单元编号（用于调试和日志）
- `I_e::Float64`: 该单元的无量纲电流（由分流求解器提供，相对 I_typ）
- `T_e::Float64`: 该单元的无量纲温度（由热场提供，相对 T_ref）
- `jacobi::String`: 雅可比矩阵更新策略 ("constant" 或 "update")

# 返回
- `M_e::SparseMatrixCSC`: 该单元的质量矩阵（电化学部分）
- `K_e::SparseMatrixCSC`: 该单元的刚度矩阵（电化学部分）
- `F_e::Vector{Float64}`: 该单元的载荷向量（边界条件）
- `variables_e::Dict`: 该单元的电化学变量字典

# 注意事项
1. yt_e 是该单元的**局部**状态向量，不是全局向量
2. I_e 和 T_e 应为无量纲值（已按各自的特征尺度归一化）
3. 返回的 M_e, K_e, F_e 用于全局装配时的 blockdiag 操作
4. 如果 jacobi="constant"，会复用 case.param.NE.M_d 等缓存矩阵
   （注意：多线程并行时应使用 jacobi="update" 避免竞争）
5. yt_e 可以是向量或矩阵（列向量），函数内部自动转换

# 示例
```julia
# 初始化单元状态（可从全局初始化复制，可能是矩阵）
yt_e = ModelInitialisation(case)  # 可能返回 Matrix{Float64}

# 单元电流和温度（由外部提供）
I_e = 1.0  # 1C放电
T_e = 1.02  # 略高于参考温度

# 求解该单元的电化学响应（自动处理矩阵输入）
M_e, K_e, F_e, vars_e = SPMe_element(case, yt_e, 0.0, 1; I_e=I_e, T_e=T_e)

# 提取单元电压
V_e = vars_e["cell voltage"]
```

# 物理意义
在多SPMe架构中，每个单元维护独立的浓度场，根据该单元的电流和温度历史演化：
- 温度高的单元: 反应速率快，浓度梯度大，极化大
- 电流大的单元: 锂离子消耗快，浓度下降快
- 这种架构真实反映了空间异质性，适用于温度分布显著的场景
"""
function SPMe_element(case::Case, yt_e::Array{Float64}, t::Float64, e::Int; I_e::Float64, T_e::Float64, jacobi::String="update")
    # 0) 确保 yt_e 是向量（兼容 ModelInitialisation 返回的列向量矩阵）
    yt_e_vec = vec(yt_e)
    
    # 1) 调用 SPMe_variables，覆写 I_app 和 T_e
    # 这里 yt_e_vec 是该单元的局部状态向量，SPMe_variables 会从中提取浓度场
    variables_e = SPMe_variables(case, yt_e_vec, t; I_app=I_e, T_e=T_e)
    
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

# ========================================================================
# solve_branch_currents_newton - 辅助函数
# ========================================================================

# 调试输出函数
"""检查并报告 NaN/Inf 值（调试用）"""
function _debug_check_prefactors(prefactor_n, prefactor_p, csn_av, csp_av, 
                                  u_n_ref_val, u_p_ref_val, du_n_dT_val, du_p_dT_val, 
                                  c_sigma, cn_surf, cp_surf, ce_n_gs, ce_p_gs)
    has_nan = !isfinite(prefactor_n) || !isfinite(prefactor_p) || 
              !isfinite(csn_av) || !isfinite(csp_av) || 
              !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val) || 
              !isfinite(c_sigma)
    
    if !has_nan
        return false
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 预计算值包含 NaN/Inf")
    println("="^80)
    println("📊 预因子: prefactor_n=$prefactor_n, prefactor_p=$prefactor_p")
    println("📊 平均浓度: csn_av=$csn_av, csp_av=$csp_av")
    println("📊 开路电位: u_n=$u_n_ref_val, u_p=$u_p_ref_val")
    println("📊 温度系数: du_n_dT=$du_n_dT_val, du_p_dT=$du_p_dT_val")
    println("📊 电导: c_sigma=$c_sigma")
    println("\n💡 输入浓度前5个值:")
    println("  cn_surf: $(cn_surf[1:min(5,end)])")
    println("  cp_surf: $(cp_surf[1:min(5,end)])")
    println("  ce_n_gs: $(ce_n_gs[1:min(5,end)])")
    println("  ce_p_gs: $(ce_p_gs[1:min(5,end)])")
    println("="^80 * "\n")
    
    return true
end

"""检查单元系数是否有效"""
function _debug_check_coefficients(e, has_nan_prefactor, C1, C2, alpha_p, alpha_n, C5, 
                                    T_e, j0_n, j0_p, u_n_val_T, u_p_val_T)
    if has_nan_prefactor || e > 1
        return  # 只打印第一个异常单元，且预计算值正常时
    end
    
    has_nan_coeff = !isfinite(C1) || !isfinite(C2) || 
                    !isfinite(alpha_p) || !isfinite(alpha_n) || !isfinite(C5)
    
    if !has_nan_coeff
        return
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 单元 $e 系数异常")
    println("="^80)
    println("📊 T_e=$T_e, j0_n=$j0_n, j0_p=$j0_p")
    println("📊 u_n(T)=$u_n_val_T, u_p(T)=$u_p_val_T")
    println("📊 系数: C1=$C1, C2=$C2")
    println("   alpha_p=$alpha_p, alpha_n=$alpha_n, C5=$C5")
    println("="^80 * "\n")
end

"""检查初始电压是否有效"""
function _debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, 
                                      coeffs, I_total, ne)
    if has_nan_prefactor || isfinite(V)
        return
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 初始电压异常: V=$V")
    println("="^80)
    println("  I_total=$I_total, ne=$ne")
    println("  前3个异常单元:")
    
    count = 0
    for e in 1:ne
        if !isfinite(V_branches[e]) && count < 3
            count += 1
            println("  单元 $e: V=$(V_branches[e]), I=$(I_e[e])")
            println("    C1=$(coeffs[e].C1), C2=$(coeffs[e].C2)")
            println("    α_p=$(coeffs[e].alpha_p), α_n=$(coeffs[e].alpha_n)")
        end
    end
    println("="^80 * "\n")
end

# 电化学计算函数
"""标量化：将数组转为标量（取第一个元素）"""
_scalarize(x) = isa(x, Number) ? Float64(x) : Float64(x[1])

"""计算电化学预因子"""
function _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
    
    # 预因子计算
    prefactor_n = IntV(abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5, mesh_ne) / param.NE.thickness
    prefactor_p = IntV(abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5, mesh_pe) / param.PE.thickness
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    
    # 开路电位
    u_n_ref = param.NE.U(cn_surf)
    u_p_ref = param.PE.U(cp_surf)
    du_n_dT = param.NE.dUdT(cn_surf)
    du_p_dT = param.PE.dUdT(cp_surf)
    
    # 标量化
    u_n_ref_val = _scalarize(u_n_ref)
    u_p_ref_val = _scalarize(u_p_ref)
    du_n_dT_val = _scalarize(du_n_dT)
    du_p_dT_val = _scalarize(du_p_dT)
    
    # 固相电导
    c_sigma = (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig) / 3.0
    
    return (prefactor_n=prefactor_n, prefactor_p=prefactor_p,
            csn_av=csn_av, csp_av=csp_av,
            u_n_ref_val=u_n_ref_val, u_p_ref_val=u_p_ref_val,
            du_n_dT_val=du_n_dT_val, du_p_dT_val=du_p_dT_val,
            c_sigma=c_sigma,
            cn_surf=cn_surf, cp_surf=cp_surf, ce_n_gs=ce_n_gs, ce_p_gs=ce_p_gs)
end

"""计算单个单元的电化学系数"""
function _compute_element_coefficients(e, T_e, param, prefactors, T_ref, debug_mode=false)
    # 交换电流密度
    arr_n = Arrhenius(param.NE.Eac_k, T_e)
    arr_p = Arrhenius(param.PE.Eac_k, T_e)
    j0_n = param.NE.k * arr_n * prefactors.prefactor_n
    j0_p = param.PE.k * arr_p * prefactors.prefactor_p
    
    # 电解液电导率
    kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps^param.NE.brugg
    kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps^param.PE.brugg
    kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps^param.SP.brugg
    R_EL = param.NE.thickness / (3.0 * kappa_ne) + param.SP.thickness / kappa_sp + param.PE.thickness / (3.0 * kappa_pe)
    
    # 开路电位（温度相关）
    u_n = prefactors.u_n_ref_val + (T_e - T_ref) * prefactors.du_n_dT_val
    u_p = prefactors.u_p_ref_val + (T_e - T_ref) * prefactors.du_p_dT_val
    
    # 计算系数
    C1 = (u_p - u_n) + 2.0 * T_e * (1.0 - param.EL.tplus) * (prefactors.csp_av - prefactors.csn_av) / param.EL.ce0
    C2 = 2.0 * T_e
    alpha_p = -1.0 / (2.0 * j0_p * param.PE.as * param.PE.thickness)
    alpha_n = 1.0 / (2.0 * j0_n * param.NE.as * param.NE.thickness)
    C5 = R_EL + prefactors.c_sigma
    
    # 调试检查（仅第一个异常单元）
    if debug_mode
        _debug_check_coefficients(e, false, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n, u_p)
    end
    
    return (C1=C1, C2=C2, alpha_p=alpha_p, alpha_n=alpha_n, C5=C5)
end

"""批量计算所有单元的系数"""
function _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode=false)
    coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
    for e in 1:ne
        coeffs[e] = _compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref, debug_mode)
    end
    return coeffs
end

# 分支电压模型
"""计算分支电压 V = C1 + C2*(asinh(α_p*I) - asinh(α_n*I)) - C5*I"""
function _branch_voltage(coeff, I::Float64)
    apI = coeff.alpha_p * I
    anI = coeff.alpha_n * I
    return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
end

"""计算分支电压对电流的导数 dV/dI"""
function _branch_dVdI(coeff, I::Float64)
    apI = coeff.alpha_p * I
    anI = coeff.alpha_n * I
    denom_p = sqrt(1.0 + apI * apI)
    denom_n = sqrt(1.0 + anI * anI)
    return coeff.C2 * (coeff.alpha_p / denom_p - coeff.alpha_n / denom_n) - coeff.C5
end

# 初始化和边界检查
"""初始化单元电流猜测"""
function _initialize_currents(ne, w, I_total, x_prev)
    I_e = x_prev !== nothing && length(x_prev) == ne ? copy(x_prev) : (w .* I_total)
    
    # 归一化到满足总电流约束
    sx = sum(w .* I_e)
    if sx != 0.0
        I_e .*= (I_total / sx)
    else
        I_e .= w .* I_total
    end
    
    return I_e
end

"""检查电压边界

返回 true 如果在边界内，false 如果超出但在软边界内（允许继续），
只有严重超出时才抛出错误。

软边界设置为 ±5% 的容差，允许主循环正常处理截止条件。
"""
function _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="")
    # 正常范围内
    if V_MIN <= V <= V_MAX
        return true
    end
    
    V_phys = V * phi_scale
    V_MIN_phys = V_MIN * phi_scale
    V_MAX_phys = V_MAX * phi_scale
    
    # 计算软边界（允许 5% 的超出，让主循环处理截止）
    V_range = V_MAX - V_MIN
    soft_margin = 0.05 * V_range  # 5% 容差
    V_MIN_soft = V_MIN - soft_margin
    V_MAX_soft = V_MAX + soft_margin
    
    # 在软边界内：发出警告但继续
    if V_MIN_soft <= V <= V_MAX_soft
        # 电压略微超出正常范围，但在软边界内
        # 这通常发生在充电接近截止电压时，应该由主循环处理
        return false  # 返回 false 表示超出正常范围但可以继续
    end
    
    # 严重超出软边界：抛出错误
    error_msg = "thermal2D common voltage severely out of bounds$context: " *
                "V(nd)=$V, V(V)=$V_phys, " *
                "allowed [$V_MIN, $V_MAX] nd -> [$V_MIN_phys, $V_MAX_phys] V, " *
                "soft bounds [$V_MIN_soft, $V_MAX_soft] nd; " *
                "I_total_nd=$I_total, sum(w.*I_e)=$(sum(w .* I_e))"
    
    throw(ErrorException(error_msg))
end

# 牛顿迭代求解器
"""
    _detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total)

检测达到截止电压的单元。

# 返回
- `active_mask`: 布尔数组，true 表示单元活跃（可接受电流）
- `n_cutoff`: 达到截止的单元数

# 物理含义
- 充电时 (I_total < 0)：如果单元的 OCV >= V_MAX，该单元已满充
- 放电时 (I_total > 0)：如果单元的 OCV <= V_MIN，该单元已完全放电
- 静置时 (I_total ≈ 0)：所有单元都是活跃的
"""
function _detect_cutoff_elements(coeffs, ne::Int, V_MIN::Float64, V_MAX::Float64, I_total::Float64)
    active_mask = trues(ne)
    
    # 计算各单元的开路电压 (OCV = C1，当 I=0 时的电压)
    OCV = [coeffs[e].C1 for e in 1:ne]
    
    # 添加小容差避免数值问题
    tol = 0.001 * (V_MAX - V_MIN)  # 0.1% 容差
    
    if I_total < -1e-10  # 充电
        # 充电时，如果单元 OCV 已达到上限，无法继续充电
        for e in 1:ne
            if OCV[e] >= V_MAX - tol
                active_mask[e] = false
            end
        end
    elseif I_total > 1e-10  # 放电
        # 放电时，如果单元 OCV 已达到下限，无法继续放电
        for e in 1:ne
            if OCV[e] <= V_MIN + tol
                active_mask[e] = false
            end
        end
    end
    # 静置时 (I_total ≈ 0)，所有单元保持活跃
    
    n_cutoff = sum(.!active_mask)
    return active_mask, n_cutoff
end

"""牛顿迭代主循环（支持部分单元截止）"""
function _newton_iteration!(I_e, V, ne, w, I_total, coeffs; 
                           tol_V=1e-8, tol_I=1e-10, max_iters=25,
                           active_mask::Union{Nothing, BitVector}=nothing)
    converged = false
    last_iter = 0
    F = zeros(Float64, ne)
    dFdI = similar(F)
    I_trial = similar(I_e)
    
    # 如果没有提供 active_mask，所有单元都是活跃的
    if active_mask === nothing
        active_mask = trues(ne)
    end
    
    # 获取活跃单元索引
    active_idx = findall(active_mask)
    n_active = length(active_idx)
    
    # 如果没有活跃单元，直接返回
    if n_active == 0
        return V, true, 0
    end
    
    # 计算活跃单元的权重和（用于归一化）
    w_active_sum = sum(w[active_idx])
    
    for iter in 1:max_iters
        last_iter = iter
        
        # 计算残差和雅可比（只对活跃单元）
        for e in 1:ne
            if active_mask[e]
                V_e = _branch_voltage(coeffs[e], I_e[e])
                F[e] = V_e - V
                dFdI[e] = _branch_dVdI(coeffs[e], I_e[e])
                
                # 防止奇异雅可比
                if abs(dFdI[e]) < 1e-12
                    dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5
                end
            else
                # 非活跃单元：电流固定为 0，残差为 0
                F[e] = 0.0
                dFdI[e] = -1.0  # 任意非零值
            end
        end
        
        # 检查收敛（只检查活跃单元）
        res_V = n_active > 0 ? maximum(abs.(F[active_idx])) : 0.0
        res_I = sum(w .* I_e) - I_total
        
        if res_V <= tol_V && abs(res_I) <= tol_I
            converged = true
            break
        end
        
        # 牛顿步（只对活跃单元）
        denom = sum(w[active_idx] ./ dFdI[active_idx])
        abs(denom) < 1e-12 && break
        
        num = -res_I + sum(w[active_idx] .* F[active_idx] ./ dFdI[active_idx])
        ΔV = num / denom
        
        # 计算电流增量
        ΔI = zeros(Float64, ne)
        for e in active_idx
            ΔI[e] = ((-F[e]) + ΔV) / dFdI[e]
        end
        # 非活跃单元的 ΔI 保持为 0
        
        # 线搜索
        λ, V_trial = _line_search(I_e, V, ΔI, ΔV, I_trial, ne)
        λ == 0.0 && break
        
        # 更新（只更新活跃单元）
        for e in 1:ne
            if active_mask[e]
                I_e[e] = I_trial[e]
            end
            # 非活跃单元保持 I_e[e] = 0
        end
        V = V_trial
    end
    
    return V, converged, last_iter
end

"""线搜索：确保更新后的值有效"""
function _line_search(I_e, V, ΔI, ΔV, I_trial, ne; max_attempts=12)
    λ = 1.0
    
    for attempt in 1:max_attempts
        # 试探电压
        V_trial = V + λ * ΔV
        !isfinite(V_trial) && (λ *= 0.5; continue)
        
        # 试探电流
        ok = true
        @inbounds for e in 1:ne
            val = I_e[e] + λ * ΔI[e]
            if !isfinite(val) || abs(val) > 1e12
                ok = false
                break
            end
            I_trial[e] = val
        end
        
        ok && return λ, V_trial
        λ *= 0.5
    end
    
    return 0.0, V  # 失败
end

# ========================================================================
# 主函数：solve_branch_currents_newton（精简版）
# ========================================================================

"""
    solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev)

非线性分流求解器（精简版）

使用牛顿法求解电流分配问题：
- 未知量：各单元电流 I_e 和公共端电压 V
- 约束：V_e(I_e, T_e) = V (每个单元) 和 Σ(I_e*A_e) = I_total

# 参数
- `case`: 案例对象
- `variables`: 变量字典
- `yt`: 状态向量
- `t`: 时间
- `I_total`: 总电流（无量纲）
- `areas`: 各单元面积
- `Te_prev`: 各单元温度
- `x_prev`: 上一步的电流分配（可选）

# 返回
- `variables`: 更新后的变量字典
- `I_e`: 各单元电流
- `V`: 公共端电压
"""
function solve_branch_currents_newton(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing)
    # 1. 初始化
    ne = length(areas)
    A_global = sum(areas)
    w = areas ./ A_global  # 面积权重
    phi_scale = case.param.scale.phi
    V_MIN = hasproperty(case.param_dim.cell, :v_l) ? case.param_dim.cell.v_l / phi_scale : -Inf
    V_MAX = hasproperty(case.param_dim.cell, :v_h) ? case.param_dim.cell.v_h / phi_scale : Inf
    
    # 2. 计算电化学预因子
    param = case.param
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    prefactors = _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
    
    # 调试：检查预因子
    debug_mode = hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
    has_nan_prefactor = _debug_check_prefactors(
        prefactors.prefactor_n, prefactors.prefactor_p, 
        prefactors.csn_av, prefactors.csp_av,
        prefactors.u_n_ref_val, prefactors.u_p_ref_val,
        prefactors.du_n_dT_val, prefactors.du_p_dT_val,
        prefactors.c_sigma,
        prefactors.cn_surf, prefactors.cp_surf, 
        prefactors.ce_n_gs, prefactors.ce_p_gs
    )
    
    # 3. 计算各单元系数
    T_ref = case.param.cell.T0
    coeffs = _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode)
    
    # 4. 检测达到截止电压的单元
    active_mask, n_cutoff = _detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total)
    
    # 5. 初始化电流猜测
    I_e = _initialize_currents(ne, w, I_total, x_prev)
    
    # 将截止单元的电流设为 0
    for e in 1:ne
        if !active_mask[e]
            I_e[e] = 0.0
        end
    end
    
    # 计算初始电压（使用活跃单元的平均值）
    active_idx = findall(active_mask)
    if !isempty(active_idx)
        V_branches = [_branch_voltage(coeffs[e], I_e[e]) for e in active_idx]
        V = sum(V_branches) / length(active_idx)
    else
        # 所有单元都达到截止，使用 OCV 的平均值作为公共电压
        V = sum(coeffs[e].C1 for e in 1:ne) / ne
    end
    
    # 调试：检查初始电压
    V_branches_all = [_branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
    _debug_check_initial_voltage(has_nan_prefactor, V, V_branches_all, I_e, coeffs, I_total, ne)
    
    # 6. 牛顿迭代求解（只对活跃单元）
    if !isempty(active_idx)
        V, converged, last_iter = _newton_iteration!(I_e, V, ne, w, I_total, coeffs; 
                                                     active_mask=active_mask)
    else
        converged = true
        last_iter = 0
    end
    
    # 7. 归一化确保总电流约束（只对活跃单元）
    if !isempty(active_idx)
        sx = sum(w .* I_e)
        if abs(sx) > 1e-12
            # 按比例调整活跃单元的电流
            scale_factor = I_total / sx
            for e in active_idx
                I_e[e] *= scale_factor
            end
        elseif abs(I_total) > 1e-12
            # sx ≈ 0 但 I_total ≠ 0：按面积权重分配
            w_active_sum = sum(w[active_idx])
            if w_active_sum > 0
                for e in active_idx
                    I_e[e] = I_total * w[e] / w_active_sum
                end
            end
        end
    end
    
    # 8. 边界检查（软边界，允许主循环处理截止条件）
    voltage_in_bounds = _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e)
    
    # 9. 写入结果（包括边界状态和截止信息）
    variables["thermal2D element current"] = I_e
    variables["thermal2D element current A"] = case.param.scale.I_typ .* I_e
    variables["thermal2D common voltage"] = V
    variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5
    variables["thermal2D Vsolve iters"] = Float64(last_iter)
    variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0
    
    # 截止状态信息
    variables["thermal2D n_active_elements"] = Float64(sum(active_mask))
    variables["thermal2D n_cutoff_elements"] = Float64(n_cutoff)
    variables["thermal2D active_mask"] = Float64.(active_mask)  # 1.0 = 活跃, 0.0 = 截止
    
    # 各单元的开路电压（用于诊断）
    OCV_elements = [coeffs[e].C1 * phi_scale for e in 1:ne]  # 转换为物理单位 (V)
    variables["thermal2D element OCV"] = OCV_elements
    
    return variables, I_e, V
end