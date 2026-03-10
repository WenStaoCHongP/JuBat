# ThermalDistributed.jl - 分布式热传导模型
# 
# 主要功能：
# 1. 2D 热传导方程装配（各向异性）
# 2. 边界条件处理
# 3. 热源计算
#
# 设计原则：
# - 所有计算使用无量纲参数（在 SetParams.jl 中已归一化）
# - 单元面积和层权重从 Jellyrollmodel.jl 获取
# - 始终使用各向异性热导（径向/切向等效热导率）

# ========================================================================
# 1D 热传导（占位）
# ========================================================================

"""
    ThermalDistributed1D(case, variables)

一维分布式热传导方程装配（占位实现）
"""
function ThermalDistributed1D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    mesh = case.mesh["thermal1D"]
    nnode = mesh.nlen
    MT = spzeros(nnode, nnode)
    KT = spzeros(nnode, nnode)
    FT = zeros(Float64, nnode)
    return MT, KT, FT
end

# ========================================================================
# 2D 热传导装配
# ========================================================================

"""
    ThermalDistributed2D(case, variables)

果冻卷 2D 热传导装配：(ρc) ∂T/∂t = ∇·(K ∇T) + q

使用各向异性热导（径向+切向等效热导率）。

**约定**：返回的 KT 已包含负号，即 M dT/dt = KT T + F

# 返回
- `MT`: 质量矩阵 (无量纲)
- `KT`: 刚度矩阵 (无量纲，包含负号)
- `FT`: 载荷向量 (无量纲)
"""
function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    mesh = case.mesh["thermal2D"]
    param = case.param
    scale = case.param_dim.scale
    L_th = scale.L_th
    
    nnode = mesh.nlen
    ne = size(mesh.element, 1)
    
    # 高斯积分数据
    Ni = mesh.gs.Ni
    dNdx, dNdy = mesh.gs.dNidx[:, 1:4], mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    Vi, Vj = mesh.element[mesh.gs.ele, :], mesh.element[mesh.gs.ele, :]
    ele_of_gp = mesh.gs.ele
    
    # 获取层权重
    fks = jellyroll_element_properties(mesh, case.param_dim)[2]
    
    # ========== 质量矩阵 ==========
    # 计算各单元的等效热容
    ρc_e = zeros(Float64, ne)
    @inbounds for e in 1:ne
        ρc_e[e] = fks[e,1]*param.NE.rho + fks[e,2]*param.SP.rho + 
                  fks[e,3]*param.PE.rho + fks[e,4]*param.PCC.rho + fks[e,5]*param.NCC.rho
    end
    ρc_weights = ρc_e[ele_of_gp] .* (wJ ./ L_th^2)
    MT = Assemble(Vi, Vj, Ni, Ni, ρc_weights, nnode)
    
    # ========== 刚度矩阵（各向异性）==========
    # 计算各单元的等效热导率
    lam_r_e, lam_t_e = zeros(Float64, ne), zeros(Float64, ne)
    ϵ = 1e-12
    @inbounds for e in 1:ne
        f = @view fks[e, :]
        # 径向：调和平均（串联）
        denom = f[1]/max(param.NE.lambda, ϵ) + f[2]/max(param.SP.lambda, ϵ) + 
                f[3]/max(param.PE.lambda, ϵ) + f[4]/max(param.PCC.lambda, ϵ) + f[5]/max(param.NCC.lambda, ϵ)
        lam_r_e[e] = denom > 0 ? (1.0 / denom) : 0.0
        # 切向：算术平均（并联）
        lam_t_e[e] = f[1]*param.NE.lambda + f[2]*param.SP.lambda + 
                     f[3]*param.PE.lambda + f[4]*param.PCC.lambda + f[5]*param.NCC.lambda
    end
    
    # 在高斯点旋转得到 Kxx/Kxy/Kyy
    ngs = size(Ni, 1)
    gx, gy = mesh.gs.x[:,1], mesh.gs.x[:,2]
    
    Kxx, Kxy, Kyy = zeros(Float64, ngs), zeros(Float64, ngs), zeros(Float64, ngs)
    @inbounds for g in 1:ngs
        θ = atan(gy[g], gx[g])
        c, s = cos(θ), sin(θ)
        lr, lt = lam_r_e[ele_of_gp[g]], lam_t_e[ele_of_gp[g]]
        Kxx[g] = lr*c*c + lt*s*s
        Kxy[g] = (lt - lr)*s*c
        Kyy[g] = lr*s*s + lt*c*c
    end
    
    # 加负号与电化学约定统一
    cxx = -Kxx .* (wJ ./ L_th^2)
    cxy = -Kxy .* (wJ ./ L_th^2)
    cyy = -Kyy .* (wJ ./ L_th^2)
    
    KT_xx = Assemble(Vi, Vj, dNdx, dNdx, cxx, nnode)
    KT_xy = Assemble(Vi, Vj, dNdx, dNdy, cxy, nnode)
    KT_yx = Assemble(Vi, Vj, dNdy, dNdx, cxy, nnode)
    KT_yy = Assemble(Vi, Vj, dNdy, dNdy, cyy, nnode)
    
    KT = KT_xx + KT_xy + KT_yx + KT_yy
    
    # ========== 载荷向量 ==========
    FT = zeros(Float64, nnode)
    if haskey(variables, "heat_source_fields")
        q_elem = variables["heat_source_fields"]
        q_gs = q_elem[ele_of_gp]
        coeff_f = q_gs .* (wJ ./ L_th^2)
        FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)
    end
    
    return MT, KT, FT
end

# ========================================================================
# 边界条件
# ========================================================================

"""
    ThermalDistributed2D_BC(KT, FT, case, t)

应用热边界条件：
1. 外边界对流
2. Z方向冷却（surface 或 tab）
"""
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64=0.0)
    mesh = case.mesh["thermal2D"]
    
    # 识别边界节点
    is_inner, is_outer = _identify_boundary_nodes(mesh, case.param_dim, case.opt)
    
    # 应用外边界对流
    _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    
    # 应用Z方向冷却
    _apply_cool_method!(KT, FT, mesh, case, t)
    
    return nothing
end

"""识别边界节点"""
function _identify_boundary_nodes(mesh, param_dim, opt)
    nnode = mesh.nlen
    pgeo = jellyroll_spiral_params(param_dim)
    bval = max(pgeo.b, 1e-12)
    
    θ0_mesh = max(0.0, (pgeo.Rin - pgeo.a) / bval)
    θ1_mesh = min((pgeo.Rout - pgeo.a - pgeo.t_repeat) / bval, (pgeo.Rout - pgeo.a) / bval)
    
    θ_in_range = (θ0_mesh, min(θ0_mesh + 2.0*π, θ1_mesh))
    θ_out_range = (max(θ1_mesh - 2.0*π, θ0_mesh), θ1_mesh)
    tol = 1e-4
    
    is_inner = [edge_boundary(mesh, i, param_dim; which=:inner, theta_range=θ_in_range, tol=tol) for i in 1:nnode]
    is_outer = [edge_boundary(mesh, i, param_dim; which=:outer, theta_range=θ_out_range, tol=tol) for i in 1:nnode]
    
    return is_inner, is_outer
end

"""应用外边界对流边界条件"""
function _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    scale = case.param_dim.scale
    Bi = scale.h_th
    Bi == 0 && return
    
    L_th = scale.L_th
    T_amb = case.param_dim.cell.T_amb / scale.T_ref
    
    s_vals = (-0.577350269189626, 0.577350269189626)
    w_vals = (1.0, 1.0)
    
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    
    for e in 1:ne
        nodes = mesh.element[e, :]
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]), 
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            (is_outer[a] && is_outer[b]) || continue
            
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            
            L = hypot(x[b] - x[a], y[b] - y[a])
            J = L / 2
            
            ke11, ke12, ke22 = 0.0, 0.0, 0.0
            fe1, fe2 = 0.0, 0.0
            
            for (s, w) in zip(s_vals, w_vals)
                N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
                wt = Bi * w * (J / L_th)
                
                ke11 += -wt * N1 * N1
                ke12 += -wt * N1 * N2
                ke22 += -wt * N2 * N2
                fe1 += wt * T_amb * N1
                fe2 += wt * T_amb * N2
            end
            
            KT[a, a] += ke11; KT[a, b] += ke12
            KT[b, a] += ke12; KT[b, b] += ke22
            FT[a] += fe1; FT[b] += fe2
        end
    end
end

"""应用Z方向冷却"""
function _apply_cool_method!(KT, FT, mesh, case, t)
    cool_method = case.opt.cool_method
    
    if cool_method == "surface"
        _apply_cool_surface!(KT, FT, mesh, case)
    elseif cool_method == "tab"
        _apply_cool_tab!(KT, FT, mesh, case)
    end
end

"""整体表面冷却"""
function _apply_cool_surface!(KT, FT, mesh, case)
    h_surface = case.param_dim.cell.h
    H = case.param_dim.cell.width
    scale = case.param_dim.scale
    k_th, L_th, T_ref = scale.k_th, scale.L_th, scale.T_ref
    T_amb_nd = case.param_dim.cell.T_amb / T_ref
    
    conv_factor = 2.0 * h_surface / (H * k_th * L_th)
    
    ngs = length(mesh.gs.detJ)
    Ni = mesh.gs.Ni
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele = mesh.gs.ele
    nn_per_elem = size(mesh.element, 2)
    
    for g in 1:ngs
        nodes = mesh.element[ele[g], :]
        wt = conv_factor * wJ[g]
        
        for i in 1:nn_per_elem
            ni = nodes[i]
            Ni_g = Ni[g, i]
            for j in 1:nn_per_elem
                nj = nodes[j]
                Nj_g = Ni[g, j]
                KT[ni, nj] -= wt * Ni_g * Nj_g
            end
            FT[ni] += wt * T_amb_nd * Ni_g
        end
    end
end

"""极耳强化冷却"""
function _apply_cool_tab!(KT, FT, mesh, case)
    pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))
    isempty(tab_nodes) && return
    
    h_tab = case.param_dim.tab.h
    tab_area = case.param_dim.tab.area
    H = case.param_dim.cell.width
    scale = case.param_dim.scale
    k_th, L_th, T_ref = scale.k_th, scale.L_th, scale.T_ref
    T_amb_nd = case.param_dim.cell.T_amb / T_ref
    
    arc_lengths = _compute_tab_node_arc_lengths(mesh, tab_nodes)
    total_arc_length = sum(arc_lengths)
    total_arc_length < 1e-12 && return
    
    for (i, n) in enumerate(tab_nodes)
        weight = arc_lengths[i] / total_arc_length
        coeff = h_tab * tab_area * weight / (H * k_th * L_th)
        KT[n, n] -= coeff
        FT[n] += coeff * T_amb_nd
    end
end

"""计算极耳节点代表的弧长"""
function _compute_tab_node_arc_lengths(mesh, tab_nodes)
    n_nodes = length(tab_nodes)
    arc_lengths = zeros(Float64, n_nodes)
    
    n_nodes == 0 && return arc_lengths
    n_nodes == 1 && (arc_lengths[1] = 1.0; return arc_lengths)
    
    coords = [mesh.node[n, :] for n in tab_nodes]
    
    for i in 1:n_nodes
        if i == 1
            arc_lengths[i] = norm(coords[2] - coords[1]) / 2.0
        elseif i == n_nodes
            arc_lengths[i] = norm(coords[i] - coords[i-1]) / 2.0
        else
            arc_lengths[i] = (norm(coords[i] - coords[i-1]) + norm(coords[i+1] - coords[i])) / 2.0
        end
    end
    
    return arc_lengths
end

# ========================================================================
# 热源计算
# ========================================================================

"""
    compute_element_heat_sources(case, variables, I_e, T_e, ne)

计算各单元的热源（无量纲）。

# 参数
- `case`: Case 对象
- `variables`: 变量字典（需包含电化学变量）
- `I_e`: 各单元电流 Vector{Float64}(ne)
- `T_e`: 各单元温度 Vector{Float64}(ne)
- `ne`: 单元数

# 返回
- `q_elem`: 各单元热源 Vector{Float64}(ne)（无量纲）
"""
function compute_element_heat_sources(case::Case, variables::Dict, I_e::Vector{Float64}, T_e::Vector{Float64}, ne::Int)
    mesh = case.mesh["thermal2D"]
    param = case.param
    q_ref = case.param_dim.scale.q_th
    
    # 获取层权重
    fks = jellyroll_element_properties(mesh, case.param_dim)[2]
    
    # 获取电化学变量
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    j_n = variables["negative electrode interfacial current density"]
    j_p = variables["positive electrode interfacial current density"]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    
    # 材料参数
    as_n, as_p = param.NE.as, param.PE.as
    sig_n_eff = param.NE.sig * param.NE.eps_s
    sig_p_eff = param.PE.sig * param.PE.eps_s
    σ_PCC = max(param.PCC.sig, 1e-12)
    σ_NCC = max(param.NCC.sig, 1e-12)
    
    q_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        T = T_e[e]
        I = I_e[e]
        
        # 电导率（温度相关）
        kappa_ne = param.EL.kappa(param.EL.ce0, T) * param.NE.eps^param.NE.brugg
        kappa_pe = param.EL.kappa(param.EL.ce0, T) * param.PE.eps^param.PE.brugg
        kappa_sp = param.EL.kappa(param.EL.ce0, T) * param.SP.eps^param.SP.brugg
        
        # 各层热源
        Q_NE = as_n * abs(j_n) * abs(eta_n) + as_n * j_n * T * param.NE.dUdT(csn_surf) + I^2 / (3.0 * sig_n_eff) + I^2 / (3.0 * kappa_ne)
        Q_SP = I^2 / kappa_sp
        Q_PE = as_p * abs(j_p) * abs(eta_p) + as_p * j_p * T * param.PE.dUdT(csp_surf) + I^2 / (3.0 * sig_p_eff) + I^2 / (3.0 * kappa_pe)
        Q_PCC = I^2 / (3.0 * σ_PCC)
        Q_NCC = I^2 / (3.0 * σ_NCC)
        
        # 加权求和
        q_elem[e] = (fks[e,1]*Q_NE + fks[e,2]*Q_SP + fks[e,3]*Q_PE + fks[e,4]*Q_PCC + fks[e,5]*Q_NCC) / q_ref
    end
    
    return q_elem
end

"""
    heatQ_Source(case, variables, t, y_state)

计算热源并写入 variables（兼容旧接口）。
"""
function heatQ_Source(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, 
                      t::Float64, y_state)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 获取单元电流和温度
    I_e = variables["thermal2D element current"]
    T_nodes = get(variables, "T_nodes", fill(case.param.cell.T0, mesh.nlen))
    T_e = [mean(@view T_nodes[mesh.element[e, :]]) for e in 1:ne]
    
    # 计算热源
    q_elem = compute_element_heat_sources(case, variables, I_e, T_e, ne)
    
    # 写入 variables
    variables["heat_source_fields"] = q_elem
    
    return variables
end