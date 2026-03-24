function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    mesh = case.mesh["thermal2D"]
    param = case.param  # 只使用无量纲参数
    nnode = mesh.nlen
    ne = size(mesh.element, 1)

    # 高斯积分数据（网格已无量纲，detJ 和 wJ 已是无量纲值）
    Ni = mesh.gs.Ni
    dNdx, dNdy = mesh.gs.dNidx[:, 1:4], mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    Vi, Vj = mesh.element[mesh.gs.ele, :], mesh.element[mesh.gs.ele, :]
    ele_of_gp = mesh.gs.ele

    # 获取层权重
    fks = jellyroll_element_properties(mesh, case.param)[2]

    # ========== 质量矩阵 ==========
    rho_c_weights = thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)
    MT = Assemble(Vi, Vj, Ni, Ni, rho_c_weights, nnode)

    # ========== 刚度矩阵（各向异性）==========
    # 在高斯点旋转得到 Kxx/Kxy/Kyy
    gx, gy = mesh.gs.x[:, 1], mesh.gs.x[:, 2]
    k_xx, k_xy, k_yy = thermal_anisotropic_conductivity_2d(param, fks, ele_of_gp, gx, gy)

    # 加负号与电化学约定统一
    # 网格已无量纲，直接使用 wJ
    cxx = -k_xx .* wJ
    cxy = -k_xy .* wJ
    cyy = -k_yy .* wJ

    KT_xx = Assemble(Vi, Vj, dNdx, dNdx, cxx, nnode)
    KT_xy = Assemble(Vi, Vj, dNdx, dNdy, cxy, nnode)
    KT_yx = Assemble(Vi, Vj, dNdy, dNdx, cxy, nnode)
    KT_yy = Assemble(Vi, Vj, dNdy, dNdy, cyy, nnode)

    KT = KT_xx + KT_xy + KT_yx + KT_yy

    # ========== 载荷向量 ==========
    FT = zeros(Float64, nnode)
    q_elem = variables["heat_source_fields"]
    q_gs = q_elem[ele_of_gp]
    coeff_f = q_gs .* wJ
    FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)

    return MT, KT, FT
end

function apply_convection_bc(KT, FT, mesh, is_outer, case)
    K = copy(KT)
    F = copy(FT)

    Bi = case.param_dim.scale.h * case.param.cell.lambda_r  # Biot 数（统一能量尺度）
    if Bi == 0
        return K, F
    end

    param = case.param
    T_amb = param.cell.T_amb  # 已无量纲

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

            # 网格已无量纲，边长也是无量纲的
            L_edge = hypot(x[b] - x[a], y[b] - y[a])
            J = L_edge / 2

            ke11, ke12, ke22 = 0.0, 0.0, 0.0
            fe1, fe2 = 0.0, 0.0

            for (s, w) in zip(s_vals, w_vals)
                N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
                # 无量纲边界积分：ds* = J * dξ, J 已无量纲
                wt = Bi * w * J

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
    if cool_method == "none"
        # 无冷却：直接返回原始矩阵，不添加任何边界条件
        return copy(KT), copy(FT)
    elseif cool_method == "surface"
        param = case.param
        Bi = case.param_dim.scale.h  # Biot 数（统一能量尺度）
        # conv_factor = 2 * Bi / H* 
        conv_factor = 2.0 * Bi / param.cell.width

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
                F[ni] += wt * param.cell.T_amb * Ni_g
            end
        end
        return K, F
    elseif cool_method == "tab"
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        if isempty(tab_nodes)
            return copy(KT), copy(FT)
        end

        param = case.param
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
            # 无量纲形式：Bi * weight
            coeff = param.tab.h * param.tab.area * weight / param.cell.width
            K[n, n] -= coeff
            F[n] += coeff * param.cell.T_amb
        end
        return K, F
    end

    return copy(KT), copy(FT)
end

function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    K = copy(KT)
    F = copy(FT)

    if case.opt.czm_enabled && haskey(case.multi_spme_layout, "czm_mesh")
        czm_mesh = case.multi_spme_layout["czm_mesh"]
        param = case.param
        for (elem_idx, czm_elem) in enumerate(czm_mesh.cohesive_elements)
            state = czm_mesh.damage_states[elem_idx]
            D = state.D
            δ_n = state.δ_max_n
            # 使用无量纲 cohesive 参数，返回无量纲 h_eff*
            h_eff_nd = compute_gap_conductance(D, δ_n, param.cohesive)
            # czm_elem.length 已无量纲（网格已归一化）
            # 系数直接为 h_eff* * L*
            coeff = h_eff_nd * czm_elem.length
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

    is_inner, is_outer = identify_boundary_nodes(mesh, case.param)
    K, F = apply_convection_bc(K, F, mesh, is_outer, case)
    K, F = apply_cool_method(K, F, mesh, case)
    return K, F
end

function ThermalDistributed2D_Ring(case::Case, variables::Dict{String,Any})
    mesh = case.mesh["thermal2D"]
    param = case.param  # 使用无量纲参数
    nnode = mesh.nlen
    ne = size(mesh.element, 1)

    Ni = mesh.gs.Ni
    dNdx, dNdy = mesh.gs.dNidx[:, 1:4], mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    Vi, Vj = mesh.element[mesh.gs.ele, :], mesh.element[mesh.gs.ele, :]
    ele_of_gp = mesh.gs.ele

    # C* = param.cell.heat_Q, V* = param.cell.volume, (ρc)* = C* / V*
    rho_c_nd = param.cell.heat_Q / param.cell.volume

    # Mass matrix（网格已无量纲，直接使用 wJ）
    coeff_m = rho_c_nd .* wJ
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

    # Stiffness matrix（网格已无量纲）
    cr = -param.cell.lambda_r .* wJ
    ct = -param.cell.lambda_t .* wJ

    KT_r = Assemble(Vi, Vj, dNdr, dNdr, cr, nnode)
    KT_t = Assemble(Vi, Vj, dNdtheta, dNdtheta, ct, nnode)
    KT = KT_r + KT_t

    # Load vector（热源已无量纲）
    FT = zeros(Float64, nnode)
    q_elem = variables["heat_source_fields"]
    q_gs = q_elem[ele_of_gp]
    coeff_f = q_gs .* wJ
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

function compute_heat_sources(case::Case, variables::Dict,variables_elems::Union{Vector{<:Dict}, Nothing},I_e::Vector{Float64}, T_e::Vector{Float64},areas::Vector{Float64}; per_element_spme::Bool=false)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    param = case.param

    # 获取层权重
    fks = jellyroll_element_properties(mesh, param)[2]

    # 从 variables 获取预分配的数组
    q_rxn_ne = variables["thermal2D q_rxn_ne"]
    q_rev_ne = variables["thermal2D q_rev_ne"]
    q_ohm_s_ne = variables["thermal2D q_ohm_s_ne"]
    q_ohm_e_ne = variables["thermal2D q_ohm_e_ne"]
    q_sp = variables["thermal2D q_sp"]
    q_rxn_pe = variables["thermal2D q_rxn_pe"]
    q_rev_pe = variables["thermal2D q_rev_pe"]
    q_ohm_s_pe = variables["thermal2D q_ohm_s_pe"]
    q_ohm_e_pe = variables["thermal2D q_ohm_e_pe"]
    q_pcc = variables["thermal2D q_pcc"]
    q_ncc = variables["thermal2D q_ncc"]

    # 材料参数（已无量纲）
    as_n, as_p = param.NE.as, param.PE.as
    sig_n_eff = param.NE.sig * param.NE.eps_s
    sig_p_eff = param.PE.sig * param.PE.eps_s
    sigma_pcc = max(param.PCC.sig, 1e-12)
    sigma_ncc = max(param.NCC.sig, 1e-12)

    q_total = zeros(Float64, ne)

    @inbounds for e in 1:ne
        # 获取电化学变量（根据 per_element_spme 判断）
        if per_element_spme && variables_elems !== nothing
            vars_e = variables_elems[e]
            eta_n = vars_e["negative electrode overpotential"][1]
            eta_p = vars_e["positive electrode overpotential"][end]
            j_n = vars_e["negative electrode interfacial current density"]
            j_p = vars_e["positive electrode interfacial current density"]
            cn_surf = vars_e["negative particle surface lithium concentration"][1]
            cp_surf = vars_e["positive particle surface lithium concentration"][end]
        else
            eta_n = variables["negative electrode overpotential"][1]
            eta_p = variables["positive electrode overpotential"][end]
            j_n = variables["negative electrode interfacial current density"]
            j_p = variables["positive electrode interfacial current density"]
            cn_surf = variables["negative particle surface lithium concentration"][1]
            cp_surf = variables["positive particle surface lithium concentration"][end]
        end

        T = T_e[e]
        I_local = I_e[e]

        # 电导率（温度相关，无量纲）
        kappa_ne = param.EL.kappa(param.EL.ce0, T) * param.NE.eps^param.NE.brugg
        kappa_pe = param.EL.kappa(param.EL.ce0, T) * param.PE.eps^param.PE.brugg
        kappa_sp = param.EL.kappa(param.EL.ce0, T) * param.SP.eps^param.SP.brugg

        # 计算各层热源分量（无量纲）
        # 负极层
        Q_rxn_NE = as_n * abs(j_n) * abs(eta_n)
        Q_rev_NE = as_n * j_n * T * param.NE.dUdT(cn_surf)[1]
        Q_ohm_s_NE = I_local^2 / (3.0 * sig_n_eff)
        Q_ohm_e_NE = I_local^2 / (3.0 * kappa_ne)

        # 隔膜层
        Q_SP = I_local^2 / kappa_sp

        # 正极层
        Q_rxn_PE = as_p * abs(j_p) * abs(eta_p)
        Q_rev_PE = as_p * j_p * T * param.PE.dUdT(cp_surf)[1]
        Q_ohm_s_PE = I_local^2 / (3.0 * sig_p_eff)
        Q_ohm_e_PE = I_local^2 / (3.0 * kappa_pe)

        # 集流体层
        Q_PCC = I_local^2 / (3.0 * sigma_pcc)
        Q_NCC = I_local^2 / (3.0 * sigma_ncc)

        # 按层权重分配并存储（无量纲）
        q_rxn_ne[e] = fks[e,1] * Q_rxn_NE
        q_rev_ne[e] = fks[e,1] * Q_rev_NE
        q_ohm_s_ne[e] = fks[e,1] * Q_ohm_s_NE
        q_ohm_e_ne[e] = fks[e,1] * Q_ohm_e_NE
        q_sp[e] = fks[e,2] * Q_SP
        q_rxn_pe[e] = fks[e,3] * Q_rxn_PE
        q_rev_pe[e] = fks[e,3] * Q_rev_PE
        q_ohm_s_pe[e] = fks[e,3] * Q_ohm_s_PE
        q_ohm_e_pe[e] = fks[e,3] * Q_ohm_e_PE
        q_pcc[e] = fks[e,4] * Q_PCC
        q_ncc[e] = fks[e,5] * Q_NCC

        # 总热源
        q_total[e] = q_rxn_ne[e] + q_rev_ne[e] + q_ohm_s_ne[e] + q_ohm_e_ne[e] + q_sp[e] + q_rxn_pe[e] + q_rev_pe[e] + q_ohm_s_pe[e] + q_ohm_e_pe[e] + q_pcc[e] + q_ncc[e]
    end

    # 写入 variables
    variables["heat_source_fields"] = q_total * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_rxn_ne"] = q_rxn_ne * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_rev_ne"] = q_rev_ne * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_ohm_s_ne"] = q_ohm_s_ne * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_ohm_e_ne"] = q_ohm_e_ne * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_sp"] = q_sp * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_rxn_pe"] = q_rxn_pe * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_rev_pe"] = q_rev_pe * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_ohm_s_pe"] = q_ohm_s_pe * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_ohm_e_pe"] = q_ohm_e_pe * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_pcc"] = q_pcc * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["thermal2D q_ncc"] = q_ncc * case.param_dim.scale.L^3 / case.param_dim.cell.volume
    variables["total heat source"] = [sum(q_total .* areas) * case.param_dim.scale.L^3 / case.param_dim.cell.volume]

    return variables
end

function compute_heat_sources_with_czm(case::Case, variables::Dict,variables_elems::Union{Vector{<:Dict}, Nothing},I_e::Vector{Float64}, T_e::Vector{Float64},areas::Vector{Float64}, czm_mesh, mesh_data)
    # 先计算所有单元的热源
    variables = compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme=true)

    # 获取活跃单元
    active_elements = get_active_elements(czm_mesh, mesh_data)
    ne = length(variables["heat_source_fields"])

    # 创建活跃掩码
    is_active = falses(ne)
    for e in active_elements
        if 1 <= e <= ne
            is_active[e] = true
        end
    end

    # 将非活跃单元的热源设为零
    q_total = variables["heat_source_fields"]
    for e in 1:ne
        if !is_active[e]
            q_total[e] = 0.0
        end
    end

    variables["heat_source_fields"] = q_total
    variables["active_elements"] = active_elements

    # 更新总功率（仅活跃单元）
    variables["total heat source"] = [sum(q_total .* areas)]

    return variables
end

