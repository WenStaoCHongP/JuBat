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
        # 使用与求解一致的无量纲尺度 I_typ 来恢复物理电流
        variables["thermal2D element current A"] = case.param.scale.I_typ .* w .* I_e
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
    
    # DEBUG: 检查关键预计算值
    has_nan_prefactor = !isfinite(prefactor_n) || !isfinite(prefactor_p) || !isfinite(csn_av) || !isfinite(csp_av) || 
                        !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val) || !isfinite(c_sigma)
    if has_nan_prefactor
        println("\n" * "="^80)
        println("❌ [DEBUG] 预计算值包含 NaN/Inf - 这是问题的根源！")
        println("="^80)
        println("📊 预因子和平均值:")
        println("  prefactor_n = $prefactor_n $(isfinite(prefactor_n) ? "✓" : "❌ NaN/Inf")")
        println("  prefactor_p = $prefactor_p $(isfinite(prefactor_p) ? "✓" : "❌ NaN/Inf")")
        println("  csn_av = $csn_av $(isfinite(csn_av) ? "✓" : "❌ NaN/Inf")")
        println("  csp_av = $csp_av $(isfinite(csp_av) ? "✓" : "❌ NaN/Inf")")
        
        println("\n📊 开路电位和温度系数:")
        println("  u_n_ref_val = $u_n_ref_val $(isfinite(u_n_ref_val) ? "✓" : "❌ NaN/Inf")")
        println("  u_p_ref_val = $u_p_ref_val $(isfinite(u_p_ref_val) ? "✓" : "❌ NaN/Inf")")
        println("  du_n_dT_val = $du_n_dT_val $(isfinite(du_n_dT_val) ? "✓" : "❌ NaN/Inf")")
        println("  du_p_dT_val = $du_p_dT_val $(isfinite(du_p_dT_val) ? "✓" : "❌ NaN/Inf")")
        println("  c_sigma = $c_sigma $(isfinite(c_sigma) ? "✓" : "❌ NaN/Inf")")
        
        println("\n📊 输入浓度数据 (前5个值):")
        println("  cn_surf: $(cn_surf[1:min(5,length(cn_surf))])")
        println("  cp_surf: $(cp_surf[1:min(5,length(cp_surf))])")
        println("  ce_n_gs: $(ce_n_gs[1:min(5,length(ce_n_gs))])")
        println("  ce_p_gs: $(ce_p_gs[1:min(5,length(ce_p_gs))])")
        
        println("\n💡 可能原因:")
        if !isfinite(prefactor_n) || !isfinite(prefactor_p)
            println("  • prefactor 计算出错 - 检查固相表面浓度和电解液浓度是否合理")
            println("    公式: prefactor = IntV(sqrt(|cs_surf*(1-cs_surf)*ce_gs|)) / thickness")
        end
        if !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val)
            println("  • 开路电位 U(cs_surf) 返回 NaN - 检查 cs_surf 是否在合理范围 [0,1]")
        end
        println("="^80 * "\n")
    end

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
        
        # DEBUG: 只在预计算值正常但单元系数异常时才打印（避免重复），且只打印第一个异常单元
        if !has_nan_prefactor && (!isfinite(C1) || !isfinite(C2) || !isfinite(alpha_p) || !isfinite(alpha_n) || !isfinite(C5)) && e == 1
            println("\n" * "="^80)
            println("❌ [DEBUG] 单元 $e 的系数包含 NaN/Inf (但预计算值正常)")
            println("="^80)
            println("📊 温度和 Arrhenius 因子:")
            println("  T_e = $T_e $(isfinite(T_e) ? "✓" : "❌ NaN/Inf")")
            println("  arr_n = $arr_n, arr_p = $arr_p")
            
            println("\n📊 交换电流密度:")
            println("  j0_n = $j0_n $(isfinite(j0_n) ? "✓" : "❌ NaN/Inf")")
            println("  j0_p = $j0_p $(isfinite(j0_p) ? "✓" : "❌ NaN/Inf")")
            
            println("\n📊 电导率:")
            println("  kappa_ne = $kappa_ne, kappa_pe = $kappa_pe, kappa_sp = $kappa_sp")
            println("  R_EL = $R_EL")
            
            println("\n📊 开路电位 (温度相关):")
            println("  u_n_val(T_e) = $(u_n_val(T_e))")
            println("  u_p_val(T_e) = $(u_p_val(T_e))")
            
            println("\n📊 计算的系数:")
            println("  C1 = $C1 $(isfinite(C1) ? "✓" : "❌ NaN/Inf")")
            println("  C2 = $C2 $(isfinite(C2) ? "✓" : "❌ NaN/Inf")")
            println("  alpha_p = $alpha_p $(isfinite(alpha_p) ? "✓" : "❌ NaN/Inf")")
            println("  alpha_n = $alpha_n $(isfinite(alpha_n) ? "✓" : "❌ NaN/Inf")")
            println("  C5 = $C5 $(isfinite(C5) ? "✓" : "❌ NaN/Inf")")
            println("="^80 * "\n")
        end
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

    # DEBUG: 计算初始V并检查 (只在之前没有检测到NaN时才详细打印)
    V_branches = [branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
    V = sum(V_branches) / ne
    
    if !has_nan_prefactor && !isfinite(V)
        println("\n" * "="^80)
        println("❌ [DEBUG] 初始电压 V 是 NaN/Inf (但预计算值和系数正常)")
        println("="^80)
        println("  I_total = $I_total, ne = $ne")
        println("  检查前3个异常单元:")
        printed_count = 0
        for e in 1:ne
            if !isfinite(V_branches[e]) && printed_count < 3
                printed_count += 1
                println("\n  单元 $e:")
                println("    V_branch = $(V_branches[e]), I_e = $(I_e[e])")
                println("    C1=$(coeffs[e].C1), C2=$(coeffs[e].C2)")
                println("    alpha_p=$(coeffs[e].alpha_p), alpha_n=$(coeffs[e].alpha_n), C5=$(coeffs[e].C5)")
            end
        end
        println("="^80 * "\n")
    end

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
        # 简化输出 - 关键信息已在上面打印
        if !has_nan_prefactor
            println("\n⚠️ [DEBUG] 电压超出边界，但前面的计算看起来正常？")
            println("  V (V) = $V_phys, 允许范围 [$(V_MIN*phi_scale), $(V_MAX*phi_scale)]")
            println("  收敛状态: converged=$converged, iters=$last_iter")
        end
        throw(ErrorException("thermal2D common voltage out of bounds: V(nd)=$(V), V(V)=$(V_phys), allowed [$(V_MIN), $(V_MAX)] nd -> [$(V_MIN*phi_scale), $(V_MAX*phi_scale)] V; I_total_nd=$(I_total), sum(w.*I_e)=$(sum(w .* I_e))"))
    end

    variables["thermal2D element current"] = I_e
    # also expose physical branch currents [A] for diagnostics (use I_typ for consistency)
    variables["thermal2D element current A"] = case.param.scale.I_typ .* w .* I_e
    variables["thermal2D common voltage"] = V
    variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5
    variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0

    return variables, I_e, V
end