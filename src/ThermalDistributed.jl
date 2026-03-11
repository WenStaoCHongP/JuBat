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
    rho_c_weights = thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ, L_th)
    MT = Assemble(Vi, Vj, Ni, Ni, rho_c_weights, nnode)
    
    # ========== 刚度矩阵（各向异性）==========
    # 在高斯点旋转得到 Kxx/Kxy/Kyy
    gx, gy = mesh.gs.x[:, 1], mesh.gs.x[:, 2]
    k_xx, k_xy, k_yy = thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy)
    
    # 加负号与电化学约定统一
    cxx = -k_xx .* (wJ ./ L_th^2)
    cxy = -k_xy .* (wJ ./ L_th^2)
    cyy = -k_yy .* (wJ ./ L_th^2)
    
    KT_xx = Assemble(Vi, Vj, dNdx, dNdx, cxx, nnode)
    KT_xy = Assemble(Vi, Vj, dNdx, dNdy, cxy, nnode)
    KT_yx = Assemble(Vi, Vj, dNdy, dNdx, cxy, nnode)
    KT_yy = Assemble(Vi, Vj, dNdy, dNdy, cyy, nnode)
    
    KT = KT_xx + KT_xy + KT_yx + KT_yy
    
    # ========== 载荷向量 ==========
    FT = zeros(Float64, nnode)
    q_elem = variables["heat_source_fields"]
    q_gs = q_elem[ele_of_gp]
    coeff_f = q_gs .* (wJ ./ L_th^2)
    FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)
    
    return MT, KT, FT
end

function apply_convection_bc(KT, FT, mesh, is_outer, case)
    K = copy(KT)
    F = copy(FT)

    scale = case.param_dim.scale
    Bi = scale.h_th
    if Bi == 0
        return K, F
    end

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

            K[a, a] += ke11; K[a, b] += ke12
            K[b, a] += ke12; K[b, b] += ke22
            F[a] += fe1; F[b] += fe2
        end
    end

    return K, F
end


function apply_cool_method(KT, FT, mesh, case)
    cool_method = case.opt.cool_method
    if cool_method == "surface"
        h_surface = case.param_dim.cell.h
        H = case.param_dim.cell.width
        scale = case.param_dim.scale
        k_th, L_th, T_ref = scale.k_th, scale.L_th, scale.T_ref
        T_amb_nd = case.param_dim.cell.T_amb / T_ref

        conv_factor = 2.0 * h_surface / (H * k_th)

        ngs = length(mesh.gs.detJ)
        Ni = mesh.gs.Ni
        wJ = mesh.gs.weight .* mesh.gs.detJ
        ele = mesh.gs.ele
        nn_per_elem = size(mesh.element, 2)

        K = copy(KT)
        F = copy(FT)
        for g in 1:ngs
            nodes = mesh.element[ele[g], :]
            wt = conv_factor * wJ[g]

            for i in 1:nn_per_elem
                ni = nodes[i]
                Ni_g = Ni[g, i]
                for j in 1:nn_per_elem
                    nj = nodes[j]
                    Nj_g = Ni[g, j]
                    K[ni, nj] -= wt * Ni_g * Nj_g
                end
                F[ni] += wt * T_amb_nd * Ni_g
            end
        end
        return K, F
    elseif cool_method == "tab"
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        if isempty(tab_nodes)
            return copy(KT), copy(FT)
        end

        h_tab = case.param_dim.tab.h
        tab_area = case.param_dim.tab.area
        H = case.param_dim.cell.width
        scale = case.param_dim.scale
        k_th, L_th, T_ref = scale.k_th, scale.L_th, scale.T_ref
        T_amb_nd = case.param_dim.cell.T_amb / T_ref

        n_nodes = length(tab_nodes)
        arc_lengths = zeros(Float64, n_nodes)
        if n_nodes == 1
            arc_lengths[1] = 1.0
        else
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
        end

        total_arc_length = sum(arc_lengths)
        if total_arc_length < 1e-12
            return copy(KT), copy(FT)
        end

        K = copy(KT)
        F = copy(FT)
        for (i, n) in enumerate(tab_nodes)
            weight = arc_lengths[i] / total_arc_length
            coeff = h_tab * tab_area * weight / (H * k_th)
            K[n, n] -= coeff
            F[n] += coeff * T_amb_nd
        end
        return K, F
    end

    return copy(KT), copy(FT)
end

function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    K = copy(KT)
    F = copy(FT)

    if case.opt.czm_enabled
        czm_mesh = get(case.multi_spme_layout, "czm_mesh", nothing)
        if czm_mesh === nothing
            try
                czm_mesh = create_czm_mesh(mesh, case.param_dim)
                case.multi_spme_layout["czm_mesh"] = czm_mesh
            catch err
                @warn "Failed to build CZM mesh for interface thermal resistance" exception=(err, catch_backtrace())
            end
        end
        if czm_mesh !== nothing
            scale = case.param_dim.scale
            L_th = scale.L_th
            k_th = scale.k_th
            for (elem_idx, czm_elem) in enumerate(czm_mesh.cohesive_elements)
                state = czm_mesh.damage_states[elem_idx]
                D = state.D
                δ_n = state.δ_max_n
                h_eff = compute_gap_conductance(D, δ_n, case.param_dim.cohesive)
                coeff = h_eff * czm_elem.length / (k_th * L_th)
                n_bot = czm_elem.nodes_bottom
                n_top = czm_elem.nodes_top
                for (nb, nt) in zip(n_bot, n_top)
                    K[nb, nb] -= coeff
                    K[nb, nt] += coeff
                    K[nt, nb] += coeff
                    K[nt, nt] -= coeff
                end
            end
        end
    end

    is_inner, is_outer = identify_boundary_nodes(mesh, case.param_dim)
    K, F = apply_convection_bc(K, F, mesh, is_outer, case)
    K, F = apply_cool_method(K, F, mesh, case)
    return K, F
end

function ThermalDistributed2D_Ring(case::Case, variables::Dict{String,Any})
    mesh = case.mesh["thermal2D"]
    scale = case.param_dim.scale
    L_th = scale.L_th
    k_th = scale.k_th
    rho_c_th = scale.rho_c_th
    nnode = mesh.nlen
    ne = size(mesh.element, 1)

    Ni = mesh.gs.Ni
    dNdx, dNdy = mesh.gs.dNidx[:, 1:4], mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    Vi, Vj = mesh.element[mesh.gs.ele, :], mesh.element[mesh.gs.ele, :]
    ele_of_gp = mesh.gs.ele

    rho = case.param_dim.cell.rho
    cp = case.param_dim.cell.heat_Q
    k_r = case.param_dim.cell.lambda_r
    k_t = case.param_dim.cell.lambda_t
    rho_c_nd = (rho * cp) / rho_c_th
    k_r_nd = k_r / k_th
    k_t_nd = k_t / k_th

    # Mass matrix
    coeff_m = rho_c_nd .* (wJ ./ L_th^2)
    MT = Assemble(Vi, Vj, Ni, Ni, coeff_m, nnode)

    # Anisotropic conductivity in polar form
    gx, gy = mesh.gs.x[:, 1], mesh.gs.x[:, 2]
    ngs = size(Ni, 1)
    dNdr = zeros(Float64, ngs, 4)
    dNdtheta = zeros(Float64, ngs, 4)
    @inbounds for g in 1:ngs
        theta = atan(gy[g], gx[g])
        c, s = cos(theta), sin(theta)
        dNdx_g = dNdx[g, :]
        dNdy_g = dNdy[g, :]
        dNdr[g, :] = c .* dNdx_g .+ s .* dNdy_g
        dNdtheta[g, :] = -s .* dNdx_g .+ c .* dNdy_g
    end

    cr = -k_r_nd .* (wJ ./ L_th^2)
    ct = -k_t_nd .* (wJ ./ L_th^2)

    KT_r = Assemble(Vi, Vj, dNdr, dNdr, cr, nnode)
    KT_t = Assemble(Vi, Vj, dNdtheta, dNdtheta, ct, nnode)
    KT = KT_r + KT_t

    FT = zeros(Float64, nnode)
    q_elem = variables["heat_source_fields"]
    q_gs = q_elem[ele_of_gp]
    coeff_f = q_gs .* (wJ ./ L_th^2)
    FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)

    return MT, KT, FT
end

function ThermalRing2D_BC(KT, FT, case::Case, outer_nodes, t::Float64)
    mesh = case.mesh["thermal2D"]
    is_outer = falses(mesh.nlen)
    for n in outer_nodes
        is_outer[n] = true
    end
    return apply_convection_bc(KT, FT, mesh, is_outer, case)
end

function heatQ_Source(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, t::Float64, y_state)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 获取单元电流和温度
    I_e = variables["thermal2D element current"]
    T_nodes = variables["T_nodes"]
    T_e = element_nodal_mean(mesh, T_nodes)
    
    # 计算热源
    q_elem = compute_element_heat_sources(case, variables, I_e, T_e, ne)
    
    # 写入 variables
    variables["heat_source_fields"] = q_elem
    
    return variables
end

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
    sigma_pcc = max(param.PCC.sig, 1e-12)
    sigma_ncc = max(param.NCC.sig, 1e-12)
    
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
        Q_PCC = I^2 / (3.0 * sigma_pcc)
        Q_NCC = I^2 / (3.0 * sigma_ncc)
        
        # 加权求和
        q_elem[e] = (fks[e,1]*Q_NE + fks[e,2]*Q_SP + fks[e,3]*Q_PE + fks[e,4]*Q_PCC + fks[e,5]*Q_NCC) / q_ref
    end
    
    return q_elem
end

function heatQ_Source_with_czm(
    case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, 
    t::Float64, y_state, czm_mesh, mesh_data
)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 获取单元电流和温度
    I_e = variables["thermal2D element current"]
    T_nodes = get(variables, "T_nodes", fill(case.param.cell.T0, mesh.nlen))
    T_e = element_nodal_mean(mesh, T_nodes)
    
    # 获取活跃单元（未断裂退出的）
    active_elements = get_active_elements(czm_mesh, mesh_data)
    
    q_elem = compute_element_heat_sources(case, variables, I_e, T_e, ne)
    is_active = zeros(Bool, ne)
    for e in active_elements
        is_active[e] = true
    end
    for e in 1:ne
        if !is_active[e]
            q_elem[e] = 0.0
        end
    end
    
    # 写入variables
    variables["heat_source_fields"] = q_elem
    variables["active_elements"] = active_elements
    
    return variables
end

