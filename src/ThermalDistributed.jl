# ThermalDistributed.jl - 精简重构版
# 分布式热传导模型装配与边界条件
#
# 主要改进：
# 1. 提取重复的层参数获取逻辑
# 2. 统一单位转换处理
# 3. 向量化计算
# 4. 拆分长函数为子函数
# 5. 消除重复代码

# ========================================================================
# 辅助函数
# ========================================================================

"""计算各层无量纲热容 (ρc)* = (ρ·c)/ρc_ref"""
function _compute_layer_rho_c(param_dim, ρc_ref)
    layers = (:NE, :SP, :PE, :PCC, :NCC)
    return NamedTuple{layers}(
        (getfield(param_dim, layer).rho * getfield(param_dim, layer).heat_Q) / ρc_ref
        for layer in layers
    )
end

"""计算各层无量纲热导率 λ* = λ/k_ref"""
function _compute_layer_lambda(param_dim, k_ref)
    layers = (:NE, :SP, :PE, :PCC, :NCC)
    return NamedTuple{layers}(
        max(getfield(param_dim, layer).lambda, 0.0) / k_ref
        for layer in layers
    )
end

"""判断热源单位是否为 SI"""
function _is_SI_units(variables)
    if haskey(variables, "heat_source_units_code")
        code = variables["heat_source_units_code"]
        return (isa(code, Float64) && code > 0.5) || 
               (isa(code, AbstractVector) && length(code) > 0 && code[1] > 0.5)
    end
    return false
end

"""获取层权重矩阵 fks"""
function _get_layer_weights(variables, ne)
    if haskey(variables, "thermal2D layer_weights")
        fks = variables["thermal2D layer_weights"]
        return (isa(fks, AbstractArray) && size(fks, 1) == ne && size(fks, 2) >= 5) ? fks : nothing
    end
    return nothing
end

"""计算元素面积并缓存"""
function _compute_element_areas!(variables, mesh)
    if haskey(variables, "thermal2D element area")
        return variables["thermal2D element area"]
    end
    
    ne = size(mesh.element, 1)
    A = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    variables["thermal2D element area"] = A
    return A
end

"""向量化计算元素平均温度"""
function _compute_element_temperatures(T_nodes, elements)
    ne = size(elements, 1)
    T_e = zeros(Float64, ne)
    @inbounds for e in 1:ne
        T_e[e] = sum(@view T_nodes[elements[e, :]]) / 4  # Q4 单元固定4节点
    end
    return T_e
end

"""识别边界节点（缓存结果）"""
function _identify_boundary_nodes(mesh, param_dim, opt)
    nnode = mesh.nlen
    pgeo = jellyroll_spiral_params(param_dim)
    s_in = 0.0
    s_out = pgeo.t_repeat
    bval = max(pgeo.b, 1e-12)
    θ0_mesh = max(0.0, (pgeo.Rin - pgeo.a - s_in) / bval)
    θ1_mesh = min((pgeo.Rout - pgeo.a - s_out) / bval, (pgeo.Rout - pgeo.a) / bval)
    # 获取配置
    θ_in_range = hasproperty(opt, :boundary_inner_theta) ? opt.boundary_inner_theta : (θ0_mesh, min(θ0_mesh + 2.0*π, θ1_mesh))
    θ_out_range = hasproperty(opt, :boundary_outer_theta) ? opt.boundary_outer_theta : (max(θ1_mesh - 2.0*π, θ0_mesh), θ1_mesh)
    tol = hasproperty(opt, :boundary_tol) ? opt.boundary_tol : 1e-4
    
    # 向量化识别
    is_inner = [edge_boundary(mesh, i, param_dim; which=:inner, theta_range=θ_in_range, tol=tol) 
                for i in 1:nnode]
    is_outer = [edge_boundary(mesh, i, param_dim; which=:outer, theta_range=θ_out_range, tol=tol) 
                for i in 1:nnode]
    
    return is_inner, is_outer
end

# ========================================================================
# 主函数：1D 热传导（占位）
# ========================================================================

"""
    ThermalDistributed1D(case, variables)

一维分布式热传导方程装配（占位实现）
"""
function ThermalDistributed1D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    @assert haskey(case.mesh, "thermal1D") "thermal1D mesh is missing in case.mesh"
    mesh = case.mesh["thermal1D"]
    nnode = mesh.nlen
    
    # 占位：返回空矩阵
    MT = spzeros(nnode, nnode)
    KT = spzeros(nnode, nnode)
    FT = zeros(Float64, nnode)
    return MT, KT, FT
end

# ========================================================================
# 主函数：2D 热传导装配
# ========================================================================

"""
    ThermalDistributed2D(case, variables)

果冻卷 2D 热传导装配：(ρc) ∂T/∂t = ∇·(k ∇T) + q

**重要约定**：返回的 KT 已包含负号，即 M dT/dt = KT T + F

# 返回
- `MT`: 质量矩阵 (无量纲)
- `KT`: 刚度矩阵 (无量纲，包含负号)
- `FT`: 载荷向量 (无量纲)
"""
function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is missing"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "ThermalDistributed2D currently assumes Q4 mesh"
    
    # 提取尺度和基本几何
    scale = case.param_dim.scale
    ρc_ref, k_ref, q_ref, L_th = scale.rho_c_th, scale.k_th, scale.q_th, scale.L_th
    nnode = mesh.nlen
    ne = size(mesh.element, 1)
    
    # 高斯积分数据
    Ni, dNdx, dNdy = mesh.gs.Ni, mesh.gs.dNidx[:, 1:4], mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    Vi, Vj = mesh.element[mesh.gs.ele, :], mesh.element[mesh.gs.ele, :]
    
    # 装配质量矩阵
    MT = _assemble_mass_matrix(case, variables, ne, Vi, Vj, Ni, wJ, mesh.gs.ele, ρc_ref, L_th, nnode)
    
    # 装配刚度矩阵
    KT = _assemble_stiffness_matrix(case, variables, mesh, ne, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    
    # 装配载荷向量
    FT = _assemble_force_vector(variables, mesh, ne, Vi, Ni, wJ, q_ref, L_th, nnode)
    
    return MT, KT, FT
end

"""装配质量矩阵 MT = ∬ (ρc/ρc_ref) NᵀN dΩ*"""
function _assemble_mass_matrix(case, variables, ne, Vi, Vj, Ni, wJ, ele_of_gp, ρc_ref, L_th, nnode)
    fks = _get_layer_weights(variables, ne)
    
    if fks !== nothing
        # 使用层权重聚合
        ρc_layers = _compute_layer_rho_c(case.param_dim, ρc_ref)
        ρc_e = zeros(Float64, ne)
        @inbounds for e in 1:ne
            ρc_e[e] = fks[e,1]*ρc_layers.NE + fks[e,2]*ρc_layers.SP + 
                      fks[e,3]*ρc_layers.PE + fks[e,4]*ρc_layers.PCC + fks[e,5]*ρc_layers.NCC
        end
        ρc_weights = ρc_e[ele_of_gp] .* (wJ ./ L_th^2)
    else
        # Fallback：使用全局平均
        ρc_cell_nd = (case.param_dim.cell.rho * case.param_dim.cell.heat_Q) / ρc_ref
        ρc_weights = ρc_cell_nd .* (wJ ./ L_th^2)
    end
    
    return Assemble(Vi, Vj, Ni, Ni, ρc_weights, nnode)
end

"""装配刚度矩阵 KT = ∬ Bᵀ K B dΩ（各向异性/各向同性自适应）"""
function _assemble_stiffness_matrix(case, variables, mesh, ne, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    # 判断是否使用各向异性
    use_aniso = _should_use_anisotropic(case, variables, mesh)
    
    if use_aniso
        return _assemble_anisotropic_stiffness(case, variables, mesh, ne, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    else
        return _assemble_isotropic_stiffness(case, mesh, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    end
end

"""判断是否使用各向异性热导"""
function _should_use_anisotropic(case, variables, mesh)
    # 检查几何
    gx, gy = mesh.gs.x[:,1], mesh.gs.x[:,2]
    Rin = hasproperty(case.param_dim.cell, :Rin) ? case.param_dim.cell.Rin : minimum(hypot.(gx, gy))
    Rout = hasproperty(case.param_dim.cell, :Rout) ? case.param_dim.cell.Rout : maximum(hypot.(gx, gy))
    
    use_aniso = (Rout > Rin) && (case.opt.thermalmodel == "distributed2D")
    
    # 用户强制各向同性开关
    if haskey(variables, "force_isotropic_k")
        flag = variables["force_isotropic_k"]
        if isa(flag, Float64) && flag > 0.5
            use_aniso = false
        end
    end
    
    return use_aniso
end

"""装配各向异性刚度矩阵"""
function _assemble_anisotropic_stiffness(case, variables, mesh, ne, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    # 计算各元素的等效热导率
    lam_r_e, lam_t_e = _compute_effective_conductivity(case, variables, ne, k_ref)
    
    # 在高斯点旋转得到 Kxx/Kxy/Kyy
    ngs = size(mesh.gs.Ni, 1)
    gx, gy = mesh.gs.x[:,1], mesh.gs.x[:,2]
    e_of_g = mesh.gs.ele
    
    Kxx, Kxy, Kyy = zeros(Float64, ngs), zeros(Float64, ngs), zeros(Float64, ngs)
    @inbounds for g in 1:ngs
        θ = atan(gy[g], gx[g])
        c, s = cos(θ), sin(θ)
        lr, lt = lam_r_e[e_of_g[g]], lam_t_e[e_of_g[g]]
        Kxx[g] = lr*c*c + lt*s*s
        Kxy[g] = (lt - lr)*s*c
        Kyy[g] = lr*s*s + lt*c*c
    end
    
    # 加负号与电化学约定统一
    cxx, cxy, cyy = -Kxx .* (wJ ./ L_th^2), -Kxy .* (wJ ./ L_th^2), -Kyy .* (wJ ./ L_th^2)
    
    # 组装
    KT_xx = Assemble(Vi, Vj, dNdx, dNdx, cxx, nnode)
    KT_xy = Assemble(Vi, Vj, dNdx, dNdy, cxy, nnode)
    KT_yx = Assemble(Vi, Vj, dNdy, dNdx, cxy, nnode)
    KT_yy = Assemble(Vi, Vj, dNdy, dNdy, cyy, nnode)
    
    return KT_xx + KT_xy + KT_yx + KT_yy
end

"""计算各元素的等效热导率（径向+切向）"""
function _compute_effective_conductivity(case, variables, ne, k_ref)
    fks = _get_layer_weights(variables, ne)
    lam_r_e, lam_t_e = zeros(Float64, ne), zeros(Float64, ne)
    
    if fks !== nothing
        # 使用层权重聚合
        λ_layers = _compute_layer_lambda(case.param_dim, k_ref)
        ϵ = 1e-12
        @inbounds for e in 1:ne
            f = @view fks[e, :]
            # 径向：调和平均（串联）
            denom = sum(f[i] / max(λ_layers[i], ϵ) for i in 1:5)
            lam_r_e[e] = denom > 0 ? (1.0 / denom) : 0.0
            # 切向：算术平均（并联）
            lam_t_e[e] = sum(f[i] * λ_layers[i] for i in 1:5)
        end
    else
        # Fallback：全局等效
        pgeo = jellyroll_spiral_params(case.param_dim)
        lam_r_e .= pgeo.λ_r_eff / k_ref
        lam_t_e .= pgeo.λ_t_eff / k_ref
    end
    
    return lam_r_e, lam_t_e
end

"""装配各向同性刚度矩阵"""
function _assemble_isotropic_stiffness(case, mesh, Vi, Vj, dNdx, dNdy, wJ, k_ref, L_th, nnode)
    # 计算各向同性热导率
    pgeo = jellyroll_spiral_params(case.param_dim)
    λ_iso_nd = if pgeo.λ_r_eff > 0 && pgeo.λ_t_eff > 0
        ((pgeo.λ_r_eff + pgeo.λ_t_eff) / 2) / k_ref
    else
        1.0
    end
    
    # 加负号与电化学约定统一
    weights = -λ_iso_nd .* (wJ ./ L_th^2)
    KT_x = Assemble(Vi, Vj, dNdx, dNdx, weights, nnode)
    KT_y = Assemble(Vi, Vj, dNdy, dNdy, weights, nnode)
    
    return KT_x + KT_y
end

"""装配载荷向量 FT = ∬ q N dΩ"""
function _assemble_force_vector(variables, mesh, ne, Vi, Ni, wJ, q_ref, L_th, nnode)
    FT = zeros(Float64, nnode)
    
    if !haskey(variables, "heat_source_fields")
        return FT
    end
    
    q_elem = variables["heat_source_fields"]
    
    # 单位转换
    if _is_SI_units(variables)
        q_elem = q_elem ./ q_ref
    end
    
    # 组装
    if isa(q_elem, AbstractArray) && size(q_elem, 1) == ne
        qe = ndims(q_elem) == 1 ? q_elem : q_elem[:, 1]
        q_gs = qe[mesh.gs.ele]
        coeff_f = q_gs .* (wJ ./ L_th^2)
        FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)
    end
    
    return FT
end

# ========================================================================
# 边界条件
# ========================================================================

"""
    ThermalDistributed2D_BC(KT, FT, case, t)

应用热边界条件：
1. 外边界对流：-k ∂T/∂n = h (T - T_amb)
2. 内边界绝热（默认）
3. 极耳强制温度（惩罚法）
"""
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64=0.0)
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is missing"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" && mesh.dimension == 2 "BC requires Q4/2D mesh"
    
    # 识别边界节点
    is_inner, is_outer = _identify_boundary_nodes(mesh, case.param_dim, case.opt)
    
    # 应用外边界对流
    _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    
    # 应用极耳边界条件
    _apply_tab_bc!(KT, FT, mesh, case, t)
    
    return nothing
end

"""应用外边界对流边界条件"""
function _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    scale = case.param_dim.scale
    Bi = scale.h_th
    
    Bi == 0 && return  # 无对流
    
    L_th = scale.L_th
    T_amb = case.param_dim.cell.T_amb / scale.T_ref
    
    # 高斯积分点和权重
    s_vals = (-0.577350269189626, 0.577350269189626)
    w_vals = (1.0, 1.0)
    
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    
    for e in 1:ne
        nodes = mesh.element[e, :]
        # 检查四条边
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]), 
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            # 跳过非外边界边
            (is_outer[a] && is_outer[b]) || continue
            
            # 去重
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            
            # 边长和雅可比
            L = hypot(x[b] - x[a], y[b] - y[a])
            J = L / 2
            
            # 高斯积分
            ke11, ke12, ke22 = 0.0, 0.0, 0.0
            fe1, fe2 = 0.0, 0.0
            
            for (s, w) in zip(s_vals, w_vals)
                N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
                wt = Bi * w * (J / L_th)
                
                # 加负号与体内扩散项统一
                ke11 += -wt * N1 * N1
                ke12 += -wt * N1 * N2
                ke22 += -wt * N2 * N2
                fe1 += wt * T_amb * N1
                fe2 += wt * T_amb * N2
            end
            
            # 组装到全局矩阵
            KT[a, a] += ke11; KT[a, b] += ke12
            KT[b, a] += ke12; KT[b, b] += ke22
            FT[a] += fe1; FT[b] += fe2
        end
    end
end

"""应用极耳边界条件（惩罚法强制温度）"""
function _apply_tab_bc!(KT, FT, mesh, case, t)
    try
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        
        isempty(tab_nodes) && return
        
        # 调试信息
        if hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
            @info "[thermal BC] tab nodes" pos=length(pos_idx) neg=length(neg_idx)
        end
        
        # 极耳温度：线性升温
        rate_Ks = hasproperty(case.opt, :tab_heating_rate) ? case.opt.tab_heating_rate : 0.1
        penalty = hasproperty(case.opt, :tab_penalty) ? case.opt.tab_penalty : 1e12
        
        scale = case.param_dim.scale
        T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref
        T_tab_nd = T_amb_nd + (rate_Ks * t) / scale.T_ref
        
        # 惩罚法
        for n in tab_nodes
            KT[n, n] += penalty
            FT[n] += penalty * T_tab_nd
        end
    catch err
        @warn "Tab BC failed" err
    end
end

# ========================================================================
# 热源计算
# ========================================================================

"""
    heatQ_Source(case, variables, t, y_state)

计算分层热源并映射到热网格。

# 热源类型
- 反应热：Q_rxn = a_s · j · η
- 可逆热：Q_rev = a_s · j · T · dU/dT
- 固相欧姆热：Q_ohm_s = I²/(3σ)
- 液相欧姆热：Q_ohm_e = I²/(3κ)
- 集流体欧姆热：Q_cc = I²/(3σ_cc)
"""
function heatQ_Source(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, 
                      t::Float64, y_state)
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh missing"
    
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 预处理：获取必要数据
    areas = _compute_element_areas!(variables, mesh)
    T_nodes = haskey(variables, "T_nodes") ? variables["T_nodes"] : 
              fill(case.param.cell.T0, mesh.nlen)
    T_e = _compute_element_temperatures(T_nodes, mesh.element)
    I_e = variables["thermal2D element current"]
    fks = _get_layer_weights(variables, ne)
    fks = fks !== nothing ? fks : ones(Float64, ne, 5)  # Fallback
    
    # 计算热源
    q_elem = _compute_heat_sources(case, variables, T_e, I_e, fks, ne)
    
    # 单位转换并写入
    _write_heat_sources!(variables, q_elem, case)
    
    # 调试输出
    _debug_heat_sources(case, variables)
    
    return variables
end

"""计算各元素的热源"""
function _compute_heat_sources(case, variables, T_e, I_e, fks, ne)
    param = case.param
    q_elem = zeros(Float64, ne)
    
    for e in 1:ne
        Q_layers = _compute_layer_heat_sources(case, variables, T_e[e], I_e[e])
        q_elem[e] = sum(fks[e, i] * Q_layers[i] for i in 1:5)
    end
    
    return q_elem
end

"""计算各层热源 [Q_NE, Q_SP, Q_PE, Q_PCC, Q_NCC]"""
function _compute_layer_heat_sources(case, variables, T_e, I_e)
    param = case.param
    
    if !(case.opt.model in ("SPM", "SPMe"))
        return zeros(Float64, 5)
    end
    
    # 获取电化学变量
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    
    # 比表面积和电导率
    as_n, as_p = param.NE.as, param.PE.as
    sig_n_eff = param.NE.sig * param.NE.eps_s
    sig_p_eff = param.PE.sig * param.PE.eps_s
    kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps^param.NE.brugg
    kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps^param.PE.brugg
    kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps^param.SP.brugg
    
    # 负极
    Q_NE = (as_n * abs(j_n) * abs(eta_n) +                    # 反应热
            as_n * j_n * T_e * param.NE.dUdT(csn_surf) +      # 可逆热
            I_e^2 / (3.0 * sig_n_eff) +                       # 固相欧姆热
            I_e^2 / (3.0 * kappa_ne))                         # 液相欧姆热
    
    # 隔膜
    Q_SP = I_e^2 / kappa_sp
    
    # 正极
    Q_PE = (as_p * abs(j_p) * abs(eta_p) +
            as_p * j_p * T_e * param.PE.dUdT(csp_surf) +
            I_e^2 / (3.0 * sig_p_eff) +
            I_e^2 / (3.0 * kappa_pe))
    
    # 集流体
    σ_PCC = hasproperty(param, :PCC) && hasproperty(param.PCC, :sig) ? 
            max(param.PCC.sig, 1e-12) : 1e12
    σ_NCC = hasproperty(param, :NCC) && hasproperty(param.NCC, :sig) ? 
            max(param.NCC.sig, 1e-12) : 1e12
    t_PCC = hasproperty(param, :PCC) && hasproperty(param.PCC, :thickness) ? 
            param.PCC.thickness : 0.0
    t_NCC = hasproperty(param, :NCC) && hasproperty(param.NCC, :thickness) ? 
            param.NCC.thickness : 0.0
    
    Q_PCC = t_PCC > 0 ? I_e^2 / (3.0 * σ_PCC) : 0.0
    Q_NCC = t_NCC > 0 ? I_e^2 / (3.0 * σ_NCC) : 0.0
    
    return [Q_NE, Q_SP, Q_PE, Q_PCC, Q_NCC]
end

"""写入热源到 variables"""
function _write_heat_sources!(variables, q_elem, case)
    if hasproperty(case.opt, :units_thermal) && case.opt.units_thermal == "SI"
        variables["heat_source_fields"] = q_elem
        variables["heat_source_units_code"] = 1.0
    else
        q_ref = case.param_dim.scale.q_th
        variables["heat_source_fields"] = q_elem ./ q_ref
        variables["heat_source_units_code"] = 0.0
    end
end

"""调试输出热源信息"""
function _debug_heat_sources(case, variables)
    if !hasproperty(case.opt, :debug_coupling) || !case.opt.debug_coupling
        return
    end
    
    qe = variables["heat_source_fields"]
    units = _is_SI_units(variables) ? "SI W/m^3" : "nd"
    @info "[thermal] heat sources" units=units q_min=minimum(qe) q_max=maximum(qe) q_mean=mean(qe)
end

# ========================================================================
# 能量守恒诊断
# ========================================================================

"""
    energy_balance_log!(case, MT, T_prev, T_new, dt_th, variables)

能量守恒诊断：Q_gen - Q_conv = dE/dt
"""
function energy_balance_log!(case::Case, MT, T_prev::AbstractVector{<:Real}, 
                             T_new::AbstractVector{<:Real}, dt_th::Float64, 
                             variables::Dict{String,Union{Array{Float64},Float64}})
    # 仅在调试模式
    if !hasproperty(case.opt, :debug_coupling) || !case.opt.debug_coupling
        return
    end
    
    haskey(case.mesh, "thermal2D") || return
    mesh = case.mesh["thermal2D"]
    scale = case.param_dim.scale
    
    # 计算生成功率
    Q_gen_nd = _compute_generation_power(variables, mesh, scale)
    
    # 计算对流散热
    Q_conv_nd = _compute_convection_power(mesh, case, T_new, scale)
    
    # 计算储能变化率
    dE_dt_nd = _compute_storage_rate(MT, T_prev, T_new, dt_th)
    
    # 计算残差
    residual = (Q_gen_nd - Q_conv_nd) - dE_dt_nd
    rel = residual / max(abs(Q_gen_nd), 1e-12)
    
    @info "[thermal] Energy balance" Q_gen=Q_gen_nd Q_conv=Q_conv_nd dE_dt=dE_dt_nd residual=residual rel=rel
    
    # 写入变量
    try
        variables["thermal2D energy residual"] = residual
    catch
    end
end

"""计算生成功率"""
function _compute_generation_power(variables, mesh, scale)
    haskey(variables, "heat_source_fields") || return 0.0
    
    q_elem = variables["heat_source_fields"]
    q_nd = _is_SI_units(variables) ? (q_elem ./ scale.q_th) : q_elem
    
    # 计算面积
    ne = size(mesh.element, 1)
    A = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    for g in 1:ngs
        e = mesh.gs.ele[g]
        A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    
    A_nd = A ./ scale.L_th^2
    return sum(q_nd .* A_nd)
end

"""计算对流散热功率"""
function _compute_convection_power(mesh, case, T_new, scale)
    Bi = scale.h_th
    Bi == 0 && return 0.0
    
    L_th = scale.L_th
    T_amb = case.param_dim.cell.T_amb / scale.T_ref
    
    # 识别外边界
    is_outer, _ = _identify_boundary_nodes(mesh, case.param_dim, case.opt)
    
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    Q_conv = 0.0
    
    for e in 1:ne
        nodes = mesh.element[e, :]
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]), 
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            (is_outer[a] && is_outer[b]) || continue
            
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            
            L_nd = hypot(x[b] - x[a], y[b] - y[a]) / L_th
            T_bar = 0.5 * (T_new[a] + T_new[b])
            Q_conv += Bi * (T_bar - T_amb) * L_nd
        end
    end
    
    return Q_conv
end

"""计算储能变化率"""
function _compute_storage_rate(MT, T_prev, T_new, dt_th)
    M_lumped = dropdims(sum(Matrix(MT), dims=2), dims=2)
    return dot(M_lumped, T_new .- T_prev) / max(dt_th, 1e-16)
end