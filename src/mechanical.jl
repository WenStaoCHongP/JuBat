function Mechanicaloutput(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        mesh_n = case.mesh["negative particle"]
        mesh_p = case.mesh["positive particle"]
        c_n = variables["negative particle lithium concentration"]
        c_p = variables["positive particle lithium concentration"]
        eta_n = variables["negative electrode overpotential"]
        eta_p = variables["positive electrode overpotential"]
        V_cell = variables["cell voltage"] 
        T = variables["temperature"]
        stress_rn_center,stress_theta_n_surf,disp_surf_n,theta_Mn,csn_gs = Calstressdisp(param.NE, mesh_n, c_n, T)
        stress_rp_center,stress_theta_p_surf,disp_surf_p,theta_Mp,csp_gs = Calstressdisp(param.PE, mesh_p, c_p, T)
        eta_p_new = eta_p - (2/3) * stress_theta_p_surf * param.PE.Omega 
        eta_n_new = eta_n - (2/3) * stress_theta_n_surf * param.NE.Omega
        V_cell_new = V_cell  - (2/3) * stress_theta_p_surf * param.PE.Omega + (2/3) * stress_theta_n_surf * param.NE.Omega
        variables["negative particle center radial stress"] = stress_rn_center
        variables["positive particle center radial stress"] = stress_rp_center
        variables["negative particle surface tangential stress"] = stress_theta_n_surf
        variables["positive particle surface tangential stress"] = stress_theta_p_surf
        variables["negative particle surface displacement"] = disp_surf_n
        variables["positive particle surface displacement"] = disp_surf_p
        variables["negative particle concentration at gauss point"] = csn_gs
        variables["positive particle concentration at gauss point"] = csp_gs
        variables["negative particle stress coupling diffusion coefficient"] = theta_Mn
        variables["positive particle stress coupling diffusion coefficient"] = theta_Mp
        variables["negative electrode overpotential"] = eta_n_new
        variables["positive electrode overpotential"] = eta_p_new
        variables["cell voltage"] = V_cell_new[1]
    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        mesh_n = case.mesh["negative particle"]
        mesh_p = case.mesh["positive particle"]
        mesh_ne = case.mesh["negative electrode"]
        mesh_pe = case.mesh["positive electrode"]
        gs_ne = case.mesh["negative electrode"].gs
        gs_pe = case.mesh["positive electrode"].gs
        element_ne = case.mesh["negative electrode"].element
        element_pe = case.mesh["positive electrode"].element
        T = variables["temperature"][1]
        cs_n = variables["negative particle lithium concentration"]
        cs_p = variables["positive particle lithium concentration"]
        eta_n = variables["negative electrode overpotential"]
        eta_p = variables["positive electrode overpotential"]
        eta_n_gs = variables["negative electrode overpotential at Gauss point"]
        eta_p_gs = variables["positive electrode overpotential at Gauss point"]
        j0_n = variables["negative electrode exchange current density"]
        j0_p = variables["positive electrode exchange current density"]
        j0_n_gs = variables["negative electrode exchange current density at Gauss point"]
        j0_p_gs = variables["positive electrode exchange current density at Gauss point"]
        stress_rn_center = variables["negative particle center radial stress"]
        stress_rp_center = variables["positive particle center radial stress"]
        stress_theta_n_surf = variables["negative particle surface tangential stress"]
        stress_theta_p_surf = variables["positive particle surface tangential stress"]
        stress_theta_n_surf_gs = variables["negative particle surface tangential stress at gauss point"] 
        stress_theta_p_surf_gs = variables["positive particle surface tangential stress at gauss point"]
        disp_surf_n = variables["negative particle surface displacement"]
        disp_surf_p = variables["positive particle surface displacement"]
        csn_gs = variables["negative particle concentration at gauss point"]
        csp_gs = variables["positive particle concentration at gauss point"]
        theta_Mn = variables["negative particle stress coupling diffusion coefficient"]
        theta_Mp = variables["positive particle stress coupling diffusion coefficient"]

        meshnum_perparticle_n = (mesh_n.nlen/ mesh_ne.nlen)-1
        meshnum_perparticle_p = (mesh_p.nlen/ mesh_pe.nlen)-1
        n_n_gs = (mesh_ne.nlen-1) * mesh_ne.gs.order
        n_p_gs = (mesh_pe.nlen-1) * mesh_pe.gs.order
        for i = 1:mesh_ne.nlen
            mesh = PickElement(mesh_n, Int64.(collect((i-1)*meshnum_perparticle_n .+ (1:meshnum_perparticle_n))) )
            cs = cs_n[(i-1)* mesh.nlen .+ (1:mesh.nlen)]
            stress_rn_center[i],stress_theta_n_surf[i],disp_surf_n[i],theta_Mn[i],csn_gs[(i-1)*n_n_gs+1: i*n_n_gs]= Calstressdisp(param.NE, mesh, cs, T)
        end
        for i = 1:mesh_pe.nlen 
            mesh = PickElement(mesh_p, Int64.(collect((i-1)*meshnum_perparticle_p .+ (1:meshnum_perparticle_p))))
            cs = cs_p[(i-1)* mesh.nlen .+ (1:mesh.nlen)]
            stress_rp_center[i],stress_theta_p_surf[i],disp_surf_p[i],theta_Mp[i],csp_gs[(i-1)*n_p_gs+1: i*n_p_gs] = Calstressdisp(param.PE, mesh, cs, T)
        end
        eta_p_new = eta_p .- (2/3) * stress_theta_p_surf * param.PE.Omega  
        eta_n_new = eta_n .- (2/3) * stress_theta_n_surf * param.NE.Omega
        j_n = j0_n .* sinh.(0.5 .* eta_n_new ./ T) * 2.0
        j_p = j0_p .* sinh.(0.5 .* eta_p_new ./ T) * 2.0
        stress_theta_n_surf_gs =  sum(gs_ne.Ni .* stress_theta_n_surf[element_ne[gs_ne.ele,:]], dims=2)
        stress_theta_p_surf_gs =  sum(gs_pe.Ni .* stress_theta_p_surf[element_pe[gs_ne.ele,:]], dims=2)
        eta_p_gs_new = eta_p_gs .- (2/3) * stress_theta_p_surf_gs * param.PE.Omega  
        eta_n_gs_new = eta_n_gs .- (2/3) * stress_theta_n_surf_gs * param.NE.Omega
        j_n_gs = j0_n_gs .* sinh.(0.5 * eta_n_gs_new ./ T) * 2.0
        j_p_gs = j0_p_gs .* sinh.(0.5 * eta_p_gs_new ./ T) * 2.0
        variables["negative electrode interfacial current density"] = j_n
        variables["positive electrode interfacial current density"] = j_p
        variables["negative electrode interfacial current at Gauss point"] = j_n_gs
        variables["positive electrode interfacial current at Gauss point"] = j_p_gs
        variables["negative electrode overpotential"] = eta_n_new
        variables["positive electrode overpotential"] = eta_p_new 
        variables["negative electrode overpotential at Gauss point"] = eta_n_gs_new
        variables["positive electrode overpotential at Gauss point"] = eta_p_gs_new
        variables["negative particle center radial stress"] = stress_rn_center
        variables["positive particle center radial stress"] = stress_rp_center
        variables["negative particle surface tangential stress"] = stress_theta_n_surf
        variables["positive particle surface tangential stress"] = stress_theta_p_surf
        variables["negative particle surface tangential stress at gauss point"] = stress_theta_n_surf_gs
        variables["positive particle surface tangential stress at gauss point"] = stress_theta_p_surf_gs
        variables["negative particle surface displacement"] = disp_surf_n
        variables["positive particle surface displacement"] = disp_surf_p
        variables["negative particle concentration at gauss point"] = csn_gs
        variables["positive particle concentration at gauss point"] = csp_gs
        variables["negative particle stress coupling diffusion coefficient"] = theta_Mn
        variables["positive particle stress coupling diffusion coefficient"] = theta_Mp

    end
    return variables
end

# 独立热应力计算：根据温度场计算 1D/2D 热应力，避免与扩散应力耦合混杂
function thermal_stress(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    Tref = param.scale.T_ref
    T0 = hasproperty(param.cell, :T0) ? param.cell.T0 : 298.0 / Tref

    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        # 平均温度 → 两电极 1D 热应力（标量）
        Tval = variables["temperature"]
        T̄ = isa(Tval, Number) ? Tval : (sum(Tval) / max(1, length(Tval)))
        dT_K = (T̄ - T0) * Tref
        α_n = hasproperty(param.NE, :alphaT) ? getfield(param.NE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        α_p = hasproperty(param.PE, :alphaT) ? getfield(param.PE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        E_n = hasproperty(param.NE, :E) ? getfield(param.NE, :E) : 0.0
        E_p = hasproperty(param.PE, :E) ? getfield(param.PE, :E) : 0.0
        ν_n = hasproperty(param.NE, :nu) ? getfield(param.NE, :nu) : 0.0
        ν_p = hasproperty(param.PE, :nu) ? getfield(param.PE, :nu) : 0.0
        σ_th_n = E_n * α_n * dT_K / max(1e-12, (1.0 - ν_n))
        σ_th_p = E_p * α_p * dT_K / max(1e-12, (1.0 - ν_p))
        variables["negative electrode thermal stress (1D)"] = σ_th_n
        variables["positive electrode thermal stress (1D)"] = σ_th_p

    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        # 轴向 1D 向量热应力 + 可选 2D 单元热应力
        Tscalar = variables["temperature"][1]
        dT_K = (Tscalar - T0) * Tref
        mesh_ne = case.mesh["negative electrode"]
        mesh_pe = case.mesh["positive electrode"]
        α_n = hasproperty(param.NE, :alphaT) ? getfield(param.NE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        α_p = hasproperty(param.PE, :alphaT) ? getfield(param.PE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        E_n = hasproperty(param.NE, :E) ? getfield(param.NE, :E) : 0.0
        E_p = hasproperty(param.PE, :E) ? getfield(param.PE, :E) : 0.0
        ν_n = hasproperty(param.NE, :nu) ? getfield(param.NE, :nu) : 0.0
        ν_p = hasproperty(param.PE, :nu) ? getfield(param.PE, :nu) : 0.0
        σ_th_n_vec = fill(E_n * α_n * dT_K / max(1e-12, (1.0 - ν_n)), mesh_ne.nlen)
        σ_th_p_vec = fill(E_p * α_p * dT_K / max(1e-12, (1.0 - ν_p)), mesh_pe.nlen)
        variables["negative electrode thermal stress (1D)"] = σ_th_n_vec
        variables["positive electrode thermal stress (1D)"] = σ_th_p_vec

    end

    # 通用 2D 单元热应力：只要存在二维热网格与节点温度，就计算（独立于电化学模型）
    if haskey(case.mesh, "thermal2D") && haskey(variables, "T_nodes") && size(variables["T_nodes"], 1) > 0
        mesh_th = case.mesh["thermal2D"]
        T_nodes_any = variables["T_nodes"]
        Tn = isa(T_nodes_any, AbstractVector) ? T_nodes_any : T_nodes_any[:, end]
        ne_th = size(mesh_th.element, 1)
        σ_th_elem = zeros(Float64, ne_th)
        # 厚度加权等效材料参数
        α_n = hasproperty(param.NE, :alphaT) ? getfield(param.NE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        α_p = hasproperty(param.PE, :alphaT) ? getfield(param.PE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        E_n = hasproperty(param.NE, :E) ? getfield(param.NE, :E) : 0.0
        E_p = hasproperty(param.PE, :E) ? getfield(param.PE, :E) : 0.0
        ν_n = hasproperty(param.NE, :nu) ? getfield(param.NE, :nu) : 0.0
        ν_p = hasproperty(param.PE, :nu) ? getfield(param.PE, :nu) : 0.0
        wt_den = max(1e-12, (param.NE.thickness + param.PE.thickness))
        E_eff = (E_n * param.NE.thickness + E_p * param.PE.thickness) / wt_den
        ν_eff = (ν_n * param.NE.thickness + ν_p * param.PE.thickness) / wt_den
        α_eff = (α_n * param.NE.thickness + α_p * param.PE.thickness) / wt_den
        for e = 1:ne_th
            nodes = mesh_th.element[e, :]
            Te = sum(Tn[nodes]) / length(nodes)
            dT_e_K = (Te - T0) * Tref
            σ_th_elem[e] = E_eff * α_eff * dT_e_K / max(1e-12, (1.0 - ν_eff))
        end
        variables["thermal2D element thermal stress"] = σ_th_elem
    end
    return variables
end

function Calstressdisp(electrode::Electrode, mesh::Mesh, cs::Array{Float64}, T::Union{Float64, Array{Float64}})
    """
        fuction of diffusion-induced stress and its effect on diffusivity (in a particle)
            Input: electrode::Electrode -- electrode type NE/PE
                   mesh::Mesh -- particle mesh
                   cs::Array{Float64} -- particle lithium concentration distribution
                   rs::Float64 -- particle radius
                   T::Union{Float64, Array{Float64} -- temperature
            Output: stress_r_center -- radial diffusion-induced stress at the particle center
                   stress_theta_surf -- diffusion-induced tangential stress on particle surface
                   disp_surf -- particle surface displacement
                   theta_M -- Coupled stress diffusion coefficient
                   cs_gs -- lithium concentration at gauss point
    """
        rs = electrode.rs
        nu = electrode.nu
        Omega = electrode.Omega
        E = electrode.E
        cs_surf = cs[end]
        cs_center = cs[1]
        cs_gs =  sum(mesh.gs.Ni .* cs[mesh.element[mesh.gs.ele,:]], dims=2)
        cs_av = (3 /(4 * π * (rs ^ 3))) * IntV((cs_gs.*( 4* π .*(mesh.gs.x).^2)) , mesh)
        stress_r_center = (2 * Omega * E  * (cs_av - cs_center)) ./ (9 * (1 - nu))
        stress_theta_surf = (Omega * E * (cs_av - cs_surf)) ./ (3 * (1 - nu)) 
        disp_surf = (Omega * rs * cs_av) / 3 
        theta_M =  2 * E * (Omega^2) ./(T *(9 * (1 - nu)))
        return stress_r_center, stress_theta_surf, disp_surf, theta_M, cs_gs
    end

# ========================================================================
# 宏观扩散应力计算 - 2D有限元方法
# ========================================================================

"""
    diffusion_stress_2D(case, variables)

计算宏观层面的扩散应力（2D平面应力/应变问题）

基于锂浓度（SOC）变化引起的体积膨胀，结合有限元方法求解应力场和位移场。

# 输入
- `case`: 案例对象，包含网格、参数等
- `variables`: 包含温度场、SOC分布等的字典

# 选项
- `case.opt.plane_type`: 选择平面类型
  - `:stress` (默认): 平面应力假设，适用于薄板结构
  - `:strain`: 平面应变假设，适用于厚体结构或约束情况

# 输出
- 更新 `variables` 字典，添加以下场：
  - `"diffusion stress xx"`: x方向正应力 [Pa]
  - `"diffusion stress yy"`: y方向正应力 [Pa]
  - `"diffusion stress xy"`: 剪应力 [Pa]
  - `"diffusion stress vonMises"`: Von Mises等效应力 [Pa]
  - `"displacement x"`: x方向位移 [m]
  - `"displacement y"`: y方向位移 [m]
  - `"diffusion stress zz"`: z方向正应力 [Pa] (仅平面应变)

# 理论
参见文档: 
- docs/Diffusion_Stress_Macroscale_Theory.md
- docs/Plane_Stress_vs_Plane_Strain.md
"""
function diffusion_stress_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is required for 2D diffusion stress"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "diffusion_stress_2D requires Q4 mesh"
    
    param = case.param
    Tref = param.scale.T_ref
    T0 = hasproperty(param.cell, :T0) ? param.cell.T0 : 298.0 / Tref
    
    # 确定平面类型（应力或应变）
    plane_type = hasproperty(case.opt, :plane_type) ? case.opt.plane_type : :stress
    if !(plane_type in [:stress, :strain])
        @warn "未知的平面类型 $(plane_type)，使用默认值 :stress"
        plane_type = :stress
    end
    
    # 提取温度场和SOC分布
    T_nodes = haskey(variables, "T_nodes") ? variables["T_nodes"] : fill(T0, mesh.nlen)
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    
    # 计算单元级别的温度和SOC
    ne = size(mesh.element, 1)
    T_elem, SOC_elem = _compute_element_T_SOC(case, variables, mesh, ne, T_nodes)
    
    # 获取材料参数（按单元层权重加权）
    E_elem, nu_elem, alpha_elem, beta_c_elem, cs_max_elem = _get_mechanical_properties(case, variables, ne)
    
    # 平面应变时检查泊松比
    if plane_type == :strain
        max_nu = maximum(nu_elem)
        if max_nu > 0.45
            @warn "平面应变模式下泊松比过大 (max=$(max_nu))，可能导致数值不稳定"
        end
    end
    
    # 装配力学刚度矩阵
    K_mech = _assemble_mechanical_stiffness_2D(mesh, E_elem, nu_elem, plane_type)
    
    # 装配热-扩散载荷向量
    F_mech = _assemble_thermal_diffusion_load_2D(mesh, T_elem, SOC_elem, 
                                                  alpha_elem, beta_c_elem, cs_max_elem,
                                                  E_elem, nu_elem, T0, plane_type)
    
    # 施加边界条件
    K_mech, F_mech = _apply_mechanical_BC_2D(K_mech, F_mech, mesh, case)
    
    # 求解位移场
    U = _solve_mechanical_displacement_2D(K_mech, F_mech, mesh.nlen)
    
    # 恢复应力场
    if plane_type == :stress
        σ_xx, σ_yy, σ_xy, σ_vm = _recover_stress_2D(U, mesh, T_elem, SOC_elem,
                                                      alpha_elem, beta_c_elem, cs_max_elem,
                                                      E_elem, nu_elem, T0, plane_type)
        σ_zz = nothing
    else  # plane_type == :strain
        σ_xx, σ_yy, σ_xy, σ_zz, σ_vm = _recover_stress_2D(U, mesh, T_elem, SOC_elem,
                                                            alpha_elem, beta_c_elem, cs_max_elem,
                                                            E_elem, nu_elem, T0, plane_type)
    end
    
    # 写入结果（转换为有量纲）
    L_ref = hasproperty(param.scale, :L_th) ? param.scale.L_th : 1.0
    _write_mechanical_results!(variables, U, σ_xx, σ_yy, σ_xy, σ_vm, L_ref, σ_zz, plane_type)
    
    return variables
end

"""计算单元温度和SOC"""
function _compute_element_T_SOC(case, variables, mesh, ne, T_nodes)
    T_elem = zeros(Float64, ne)
    SOC_elem = zeros(Float64, ne)
    
    # 计算单元平均温度
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        T_elem[e] = sum(T_nodes[nodes]) / 4  # Q4单元4个节点
    end
    
    # 计算单元SOC
    if haskey(variables, "SOC_elem")
        # 如果已经计算过单元SOC
        SOC_elem = variables["SOC_elem"]
    elseif haskey(variables, "negative electrode lithium concentration") && 
           haskey(variables, "positive electrode lithium concentration")
        # 从颗粒浓度计算SOC
        SOC_elem = _estimate_SOC_from_concentration(case, variables, ne)
    else
        # 使用全局SOC估计
        SOC0 = hasproperty(case.param, :SOC0) ? case.param.SOC0 : 0.5
        SOC_elem .= SOC0
    end
    
    return T_elem, SOC_elem
end

"""从锂浓度估计单元SOC"""
function _estimate_SOC_from_concentration(case, variables, ne)
    param = case.param
    
    # 获取负极和正极的平均浓度
    cs_n = variables["negative electrode lithium concentration"]
    cs_p = variables["positive electrode lithium concentration"]
    
    # 计算平均值（简单策略）
    if isa(cs_n, AbstractArray) && length(cs_n) > 0
        cs_n_avg = sum(cs_n) / length(cs_n)
    else
        cs_n_avg = param.NE.cs_max * 0.5
    end
    
    # 使用负极SOC作为全局SOC估计
    cs_n_max = param.NE.cs_max
    SOC_global = cs_n_avg / cs_n_max
    
    # 所有单元使用相同SOC（简化）
    return fill(SOC_global, ne)
end

"""获取单元级别的力学参数"""
function _get_mechanical_properties(case, variables, ne)
    param = case.param
    
    # 检查是否有层权重矩阵
    fks = haskey(variables, "thermal2D layer_weights") ? 
          variables["thermal2D layer_weights"] : nothing
    
    E_elem = zeros(Float64, ne)
    nu_elem = zeros(Float64, ne)
    alpha_elem = zeros(Float64, ne)
    beta_c_elem = zeros(Float64, ne)
    cs_max_elem = zeros(Float64, ne)
    
    # 获取各层参数
    E_n = hasproperty(param.NE, :E) ? getfield(param.NE, :E) : 15e9  # 石墨: 15 GPa
    E_p = hasproperty(param.PE, :E) ? getfield(param.PE, :E) : 150e9  # NMC: 150 GPa
    E_s = hasproperty(param.SP, :E) ? getfield(param.SP, :E) : 1e9    # 隔膜: 1 GPa
    E_pcc = hasproperty(param, :PCC) && hasproperty(param.PCC, :E) ? param.PCC.E : 70e9  # Al: 70 GPa
    E_ncc = hasproperty(param, :NCC) && hasproperty(param.NCC, :E) ? param.NCC.E : 130e9 # Cu: 130 GPa
    
    nu_n = hasproperty(param.NE, :nu) ? getfield(param.NE, :nu) : 0.3
    nu_p = hasproperty(param.PE, :nu) ? getfield(param.PE, :nu) : 0.3
    nu_s = hasproperty(param.SP, :nu) ? getfield(param.SP, :nu) : 0.3
    nu_pcc = hasproperty(param, :PCC) && hasproperty(param.PCC, :nu) ? param.PCC.nu : 0.33
    nu_ncc = hasproperty(param, :NCC) && hasproperty(param.NCC, :nu) ? param.NCC.nu : 0.34
    
    alpha_n = hasproperty(param.NE, :alphaT) ? getfield(param.NE, :alphaT) : 1e-5
    alpha_p = hasproperty(param.PE, :alphaT) ? getfield(param.PE, :alphaT) : 1e-5
    alpha_s = hasproperty(param.SP, :alphaT) ? getfield(param.SP, :alphaT) : 1e-5
    alpha_pcc = hasproperty(param, :PCC) && hasproperty(param.PCC, :alphaT) ? param.PCC.alphaT : 23e-6
    alpha_ncc = hasproperty(param, :NCC) && hasproperty(param.NCC, :alphaT) ? param.NCC.alphaT : 17e-6
    
    # 浓度膨胀系数（新增参数）
    beta_c_n = hasproperty(param.NE, :beta_c) ? getfield(param.NE, :beta_c) : 0.03  # 石墨: ~3%体积变化
    beta_c_p = hasproperty(param.PE, :beta_c) ? getfield(param.PE, :beta_c) : 0.01  # NMC: ~1%体积变化
    
    cs_max_n = hasproperty(param.NE, :cs_max) ? param.NE.cs_max : 31507.0
    cs_max_p = hasproperty(param.PE, :cs_max) ? param.PE.cs_max : 51765.0
    
    # 定义材料参数向量
    E_vec = [E_n, E_s, E_p, E_pcc, E_ncc]
    nu_vec = [nu_n, nu_s, nu_p, nu_pcc, nu_ncc]
    alpha_vec = [alpha_n, alpha_s, alpha_p, alpha_pcc, alpha_ncc]
    beta_vec = [beta_c_n, 0.0, beta_c_p, 0.0, 0.0]  # 仅电极有扩散应变
    cs_max_vec = [cs_max_n, 0.0, cs_max_p, 0.0, 0.0]
    
    if fks !== nothing && size(fks, 1) == ne && size(fks, 2) >= 5
        # 使用层权重加权平均
        @inbounds for e in 1:ne
            f = @view fks[e, :]
            # 弹性模量：调和平均（串联）
            denom = sum(f[i] / max(E_vec[i], 1e6) for i in 1:5)
            E_elem[e] = denom > 0 ? (1.0 / denom) : E_n
            
            # 其他参数：算术平均
            nu_elem[e] = sum(f[i] * nu_vec[i] for i in 1:5)
            alpha_elem[e] = sum(f[i] * alpha_vec[i] for i in 1:5)
            beta_c_elem[e] = sum(f[i] * beta_vec[i] for i in 1:5)
            cs_max_elem[e] = sum(f[i] * cs_max_vec[i] for i in 1:5)
        end
    else
        # Fallback：使用全局参数
        E_avg = (E_n + E_p) / 2
        nu_avg = (nu_n + nu_p) / 2
        alpha_avg = (alpha_n + alpha_p) / 2
        beta_c_avg = (beta_c_n + beta_c_p) / 2
        cs_max_avg = (cs_max_n + cs_max_p) / 2
        
        E_elem .= E_avg
        nu_elem .= nu_avg
        alpha_elem .= alpha_avg
        beta_c_elem .= beta_c_avg
        cs_max_elem .= cs_max_avg
    end
    
    return E_elem, nu_elem, alpha_elem, beta_c_elem, cs_max_elem
end

"""装配2D力学刚度矩阵"""
function _assemble_mechanical_stiffness_2D(mesh, E_elem, nu_elem, plane_type::Symbol=:stress)
    nnode = mesh.nlen
    ne = size(mesh.element, 1)
    ndof = 2 * nnode  # 每个节点2个自由度
    
    # 高斯积分点数据
    Ni = mesh.gs.Ni
    dNdx = mesh.gs.dNidx[:, 1:4]
    dNdy = mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele_of_gp = mesh.gs.ele
    
    ngs = length(wJ)
    
    # 计算每个高斯点的弹性矩阵
    D11 = zeros(Float64, ngs)
    D12 = zeros(Float64, ngs)
    D33 = zeros(Float64, ngs)
    
    if plane_type == :stress
        # 平面应力：D = E/(1-ν²) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
        @inbounds for g in 1:ngs
            e = ele_of_gp[g]
            E = E_elem[e]
            ν = nu_elem[e]
            factor = E / max(1.0 - ν^2, 1e-12)
            D11[g] = factor
            D12[g] = factor * ν
            D33[g] = factor * (1.0 - ν) / 2.0
        end
    else  # plane_type == :strain
        # 平面应变：D = E/((1+ν)(1-2ν)) * [1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2]
        @inbounds for g in 1:ngs
            e = ele_of_gp[g]
            E = E_elem[e]
            ν = nu_elem[e]
            denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
            factor = E / denom
            D11[g] = factor * (1.0 - ν)
            D12[g] = factor * ν
            D33[g] = factor * (1.0 - 2.0*ν) / 2.0
        end
    end
    
    # 构造应变-位移矩阵的索引
    # 对于Q4单元，每个节点i有2个DOF：u_i, v_i
    # 节点编号 -> DOF编号映射：node i -> DOF 2i-1 (u), 2i (v)
    Vi_u = zeros(Int64, ngs, 4)
    Vi_v = zeros(Int64, ngs, 4)
    Vj_u = zeros(Int64, ngs, 4)
    Vj_v = zeros(Int64, ngs, 4)
    
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        for i in 1:4
            node = mesh.element[e, i]
            Vi_u[g, i] = 2 * node - 1
            Vi_v[g, i] = 2 * node
            Vj_u[g, i] = 2 * node - 1
            Vj_v[g, i] = 2 * node
        end
    end
    
    # 装配刚度矩阵的各个部分
    # K = ∫ B^T D B dΩ
    # B = [dN/dx 0; 0 dN/dy; dN/dy dN/dx]
    
    # K_uu: ∫ (dNi/dx * D11 * dNj/dx + dNi/dy * D33 * dNj/dy) dΩ
    wJ_D11 = wJ .* D11
    wJ_D33 = wJ .* D33
    K_uu_11 = Assemble(Vi_u, Vj_u, dNdx, dNdx, wJ_D11, ndof)
    K_uu_33 = Assemble(Vi_u, Vj_u, dNdy, dNdy, wJ_D33, ndof)
    K_uu = K_uu_11 + K_uu_33
    
    # K_vv: ∫ (dNi/dy * D11 * dNj/dy + dNi/dx * D33 * dNj/dx) dΩ
    K_vv_11 = Assemble(Vi_v, Vj_v, dNdy, dNdy, wJ_D11, ndof)
    K_vv_33 = Assemble(Vi_v, Vj_v, dNdx, dNdx, wJ_D33, ndof)
    K_vv = K_vv_11 + K_vv_33
    
    # K_uv: ∫ (dNi/dx * D12 * dNj/dy + dNi/dy * D33 * dNj/dx) dΩ
    wJ_D12 = wJ .* D12
    K_uv_12 = Assemble(Vi_u, Vj_v, dNdx, dNdy, wJ_D12, ndof)
    K_uv_33 = Assemble(Vi_u, Vj_v, dNdy, dNdx, wJ_D33, ndof)
    K_uv = K_uv_12 + K_uv_33
    
    # K_vu: ∫ (dNi/dy * D12 * dNj/dx + dNi/dx * D33 * dNj/dy) dΩ
    K_vu_12 = Assemble(Vi_v, Vj_u, dNdy, dNdx, wJ_D12, ndof)
    K_vu_33 = Assemble(Vi_v, Vj_u, dNdx, dNdy, wJ_D33, ndof)
    K_vu = K_vu_12 + K_vu_33
    
    # 总刚度矩阵
    K = K_uu + K_vv + K_uv + K_vu
    
    return K
end

"""装配热-扩散载荷向量"""
function _assemble_thermal_diffusion_load_2D(mesh, T_elem, SOC_elem, 
                                              alpha_elem, beta_c_elem, cs_max_elem,
                                              E_elem, nu_elem, T0, plane_type::Symbol=:stress)
    nnode = mesh.nlen
    ndof = 2 * nnode
    
    # 计算每个单元的初始应变 ε_0 = α*ΔT + β_c*c_s_max*ΔSOC
    ne = length(T_elem)
    epsilon_0_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        ΔT = T_elem[e] - T0
        ΔSOC = SOC_elem[e] - 0.5  # 假设初始SOC为0.5
        epsilon_0_elem[e] = alpha_elem[e] * ΔT + beta_c_elem[e] * cs_max_elem[e] * ΔSOC
    end
    
    # 高斯积分点数据
    dNdx = mesh.gs.dNidx[:, 1:4]
    dNdy = mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele_of_gp = mesh.gs.ele
    
    ngs = length(wJ)
    
    # 计算载荷系数
    # F = ∫ B^T D ε_0 dΩ
    coeff_u = zeros(Float64, ngs)
    coeff_v = zeros(Float64, ngs)
    
    if plane_type == :stress
        # 平面应力：D * ε_0 = E/(1-ν²) * [ε_0*(1+ν), ε_0*(1+ν), 0]^T
        @inbounds for g in 1:ngs
            e = ele_of_gp[g]
            E = E_elem[e]
            ν = nu_elem[e]
            ε_0 = epsilon_0_elem[e]
            factor = E / max(1.0 - ν^2, 1e-12) * ε_0 * (1.0 + ν) * wJ[g]
            coeff_u[g] = factor
            coeff_v[g] = factor
        end
    else  # plane_type == :strain
        # 平面应变：D * ε_0 = E/((1+ν)(1-2ν)) * [ε_0*(1+2ν), ε_0*(1+2ν), 0]^T
        @inbounds for g in 1:ngs
            e = ele_of_gp[g]
            E = E_elem[e]
            ν = nu_elem[e]
            ε_0 = epsilon_0_elem[e]
            denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
            factor = E / denom * ε_0 * (1.0 + 2.0*ν) * wJ[g]
            coeff_u[g] = factor
            coeff_v[g] = factor
        end
    end
    
    # 构造DOF索引
    Vi_u = zeros(Int64, ngs, 4)
    Vi_v = zeros(Int64, ngs, 4)
    
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        for i in 1:4
            node = mesh.element[e, i]
            Vi_u[g, i] = 2 * node - 1
            Vi_v[g, i] = 2 * node
        end
    end
    
    # 装配载荷向量
    F = zeros(Float64, ndof)
    
    # F_u = ∫ dNi/dx * D11 * ε_0 dΩ
    F_u = Assemble1D(Vi_u, dNdx, coeff_u, ndof)
    
    # F_v = ∫ dNi/dy * D11 * ε_0 dΩ
    F_v = Assemble1D(Vi_v, dNdy, coeff_v, ndof)
    
    F = F_u + F_v
    
    return F
end

"""施加力学边界条件"""
function _apply_mechanical_BC_2D(K, F, mesh, case)
    nnode = mesh.nlen
    ndof = 2 * nnode
    
    # 识别边界节点
    bc_nodes = _identify_mechanical_bc_nodes(mesh, case)
    
    # 施加固定约束（惩罚法）
    penalty = 1e12
    
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_x
            # 固定x方向位移
            dof = 2 * node - 1
            K[dof, dof] += penalty
            F[dof] = 0.0
        elseif bc_type == :fixed_y
            # 固定y方向位移
            dof = 2 * node
            K[dof, dof] += penalty
            F[dof] = 0.0
        elseif bc_type == :fixed_xy
            # 固定x和y方向位移
            dof_x = 2 * node - 1
            dof_y = 2 * node
            K[dof_x, dof_x] += penalty
            K[dof_y, dof_y] += penalty
            F[dof_x] = 0.0
            F[dof_y] = 0.0
        end
    end
    
    return K, F
end

"""识别需要施加边界条件的节点"""
function _identify_mechanical_bc_nodes(mesh, case)
    nnode = mesh.nlen
    bc_nodes = Dict{Int, Symbol}()
    
    # 简单策略：固定一个节点防止刚体位移
    # 选择最接近原点的节点
    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    r = hypot.(x, y)
    fixed_node = argmin(r)
    
    bc_nodes[fixed_node] = :fixed_xy
    
    # 可选：添加对称边界条件等
    # 这里采用最小约束，只固定一个点
    
    return bc_nodes
end

"""求解位移场"""
function _solve_mechanical_displacement_2D(K, F, nnode)
    ndof = 2 * nnode
    
    # 求解线性方程组 K*U = F
    try
        U = K \ F
    catch e
        @warn "Mechanical solve failed, using zero displacement" e
        U = zeros(Float64, ndof)
    end
    
    return U
end

"""恢复应力场"""
function _recover_stress_2D(U, mesh, T_elem, SOC_elem,
                             alpha_elem, beta_c_elem, cs_max_elem,
                             E_elem, nu_elem, T0, plane_type::Symbol=:stress)
    ne = size(mesh.element, 1)
    
    # 初始化应力数组
    σ_xx = zeros(Float64, ne)
    σ_yy = zeros(Float64, ne)
    σ_xy = zeros(Float64, ne)
    σ_zz = plane_type == :strain ? zeros(Float64, ne) : nothing
    σ_vm = zeros(Float64, ne)
    
    # 计算每个单元的初始应变
    epsilon_0_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        ΔT = T_elem[e] - T0
        ΔSOC = SOC_elem[e] - 0.5
        epsilon_0_elem[e] = alpha_elem[e] * ΔT + beta_c_elem[e] * cs_max_elem[e] * ΔSOC
    end
    
    # 在单元中心恢复应力
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        
        # 单元中心的形函数导数（ξ=0, η=0）
        dNdxi = [-0.25, 0.25, 0.25, -0.25]
        dNdeta = [-0.25, -0.25, 0.25, 0.25]
        
        # 计算雅可比矩阵
        x_nodes = mesh.node[nodes, 1]
        y_nodes = mesh.node[nodes, 2]
        
        J11 = sum(dNdxi[i] * x_nodes[i] for i in 1:4)
        J12 = sum(dNdxi[i] * y_nodes[i] for i in 1:4)
        J21 = sum(dNdeta[i] * x_nodes[i] for i in 1:4)
        J22 = sum(dNdeta[i] * y_nodes[i] for i in 1:4)
        
        detJ = J11 * J22 - J12 * J21
        
        if abs(detJ) < 1e-12
            continue
        end
        
        # 形函数导数（物理坐标）
        invdetJ = 1.0 / detJ
        dNdx = [(J22 * dNdxi[i] - J12 * dNdeta[i]) * invdetJ for i in 1:4]
        dNdy = [(-J21 * dNdxi[i] + J11 * dNdeta[i]) * invdetJ for i in 1:4]
        
        # 计算应变
        ε_xx = sum(dNdx[i] * U[2*nodes[i]-1] for i in 1:4)
        ε_yy = sum(dNdy[i] * U[2*nodes[i]] for i in 1:4)
        γ_xy = sum(dNdy[i] * U[2*nodes[i]-1] + dNdx[i] * U[2*nodes[i]] for i in 1:4)
        
        # 弹性应变
        ε_0 = epsilon_0_elem[e]
        ε_elastic_xx = ε_xx - ε_0
        ε_elastic_yy = ε_yy - ε_0
        ε_elastic_xy = γ_xy  # 剪应变没有初始值
        
        # 计算应力
        E = E_elem[e]
        ν = nu_elem[e]
        
        if plane_type == :stress
            # 平面应力
            factor = E / max(1.0 - ν^2, 1e-12)
            σ_xx[e] = factor * (ε_elastic_xx + ν * ε_elastic_yy)
            σ_yy[e] = factor * (ε_elastic_yy + ν * ε_elastic_xx)
            σ_xy[e] = factor * (1.0 - ν) / 2.0 * ε_elastic_xy
            
            # Von Mises应力（平面应力）
            σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 - σ_xx[e]*σ_yy[e] + 3.0*σ_xy[e]^2)
        else  # plane_type == :strain
            # 平面应变
            denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
            factor = E / denom
            σ_xx[e] = factor * ((1.0 - ν) * ε_elastic_xx + ν * ε_elastic_yy)
            σ_yy[e] = factor * ((1.0 - ν) * ε_elastic_yy + ν * ε_elastic_xx)
            σ_xy[e] = factor * (1.0 - 2.0*ν) / 2.0 * ε_elastic_xy
            
            # z方向应力（平面应变特有）
            σ_zz[e] = factor * ν * (ε_elastic_xx + ε_elastic_yy)
            
            # Von Mises应力（平面应变，包含 σ_zz）
            σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 + σ_zz[e]^2 - 
                          σ_xx[e]*σ_yy[e] - σ_yy[e]*σ_zz[e] - σ_zz[e]*σ_xx[e] + 
                          3.0*σ_xy[e]^2)
        end
    end
    
    if plane_type == :stress
        return σ_xx, σ_yy, σ_xy, σ_vm
    else
        return σ_xx, σ_yy, σ_xy, σ_zz, σ_vm
    end
end

"""写入力学计算结果"""
function _write_mechanical_results!(variables, U, σ_xx, σ_yy, σ_xy, σ_vm, L_ref, 
                                     σ_zz=nothing, plane_type::Symbol=:stress)
    nnode = length(U) ÷ 2
    
    # 提取位移
    u_x = U[1:2:end]
    u_y = U[2:2:end]
    
    # 写入变量（无量纲 -> 有量纲）
    variables["displacement x"] = u_x .* L_ref
    variables["displacement y"] = u_y .* L_ref
    variables["diffusion stress xx"] = σ_xx
    variables["diffusion stress yy"] = σ_yy
    variables["diffusion stress xy"] = σ_xy
    variables["diffusion stress vonMises"] = σ_vm
    
    # 平面应变时，额外写入 z 方向应力
    if plane_type == :strain && σ_zz !== nothing
        variables["diffusion stress zz"] = σ_zz
    end
    
    # 记录平面类型
    variables["plane_type"] = String(plane_type)
    
    return nothing
end

"""
    thermal_diffusion_stress_coupled_2D(case, variables)

热-扩散应力耦合计算（统一框架）

同时考虑热应力和扩散应力，两者作为初始应变叠加。
"""
function thermal_diffusion_stress_coupled_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    # 直接调用 diffusion_stress_2D，因为它已经包含了热应力和扩散应力
    return diffusion_stress_2D(case, variables)
end
