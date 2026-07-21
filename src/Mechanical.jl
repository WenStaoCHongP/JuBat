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
    thermal_diffusion_stress_2D(case, variables)

计算宏观层面的热和扩散应力（2D平面应力问题）

基于锂浓度（SOC）变化引起的体积膨胀，结合有限元方法求解应力场和位移场。

# 输入
- `case`: 案例对象，包含网格、参数等
- `variables`: 包含温度场、SOC分布等的字典

# 输出
- 更新 `variables` 字典，添加以下场：
  - `"diffusion stress xx"`: x方向正应力 [Pa]
  - `"diffusion stress yy"`: y方向正应力 [Pa]
  - `"diffusion stress xy"`: 剪应力 [Pa]
  - `"diffusion stress vonMises"`: Von Mises等效应力 [Pa]
  - `"displacement x"`: x方向位移 [m]
  - `"displacement y"`: y方向位移 [m]

"""
function thermal_diffusion_stress_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    # === 入口断言：参数集必须定义 E_coat 才能启用宏观力学 ===
    @assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0 "宏观力学分析需要 PE/NE.E_coat > 0；当前参数集未定义极片模量（E_coat=0）。请在参数文件中补全 PE.E_coat/PE.nu_coat/NE.E_coat/NE.nu_coat，或禁用 mechanicalmodel=\"full\"。"

    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "diffusion_stress_2D requires Q4 mesh"

    param = case.param
    Tref = param.scale.T_ref
    T0 = param.cell.T0

    # 提取温度场和SOC分布
    T_nodes = variables["T_nodes"]
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    soc_n_elem = variables["thermal2D element soc_n"]
    soc_p_elem = variables["thermal2D element soc_p"]
    soc_ref_n = param.NE.cs0
    soc_ref_p = param.PE.cs0
    # 获取材料参数
    # TODO Chunk 4 Task 4.x: compute_effective_coating_modulus 已移除；
    # 此处 thermal_diffusion_stress_2D 是非 CZM 路径，暂用 PE_PCC 占位。
    czm_param_cache = compute_czm_params_per_interface(case)
    E_eff = czm_param_cache.by_interface[:PE_PCC].E_eff
    ν_eff = czm_param_cache.by_interface[:PE_PCC].ν
    α_eff = czm_param_cache.by_interface[:PE_PCC].α
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0

    # 计算单元级别的温度和SOC
    ne = size(mesh.element, 1)
    T_elem = zeros(Float64, ne)
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)
    T_elem .= element_nodal_mean(mesh, T_nodes)
    
    @inbounds for e in 1:ne
        dT_elem[e] = T_elem[e] - T0
        Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
        Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
    end

    # 装配力学刚度矩阵
    nnode = mesh.nlen
    ndof = 2 * nnode  # 每个节点2个自由度
    dNdx = mesh.gs.dNidx[:, 1:4]
    dNdy = mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele_of_gp = mesh.gs.ele
    ngs = length(wJ)
    
    # 平面应力：D = E/(1-ν^2) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    D11 = fill(E_eff / (1.0 - ν_eff^2), ngs)
    D12 = fill(E_eff / (1.0 - ν_eff^2) * ν_eff, ngs)
    D33 = fill(E_eff / (1.0 - ν_eff^2) * (1.0 - ν_eff) / 2.0, ngs)
    
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
    
    wJ_D11 = wJ .* D11
    wJ_D33 = wJ .* D33
    K_uu = Assemble(Vi_u, Vj_u, dNdx, dNdx, wJ_D11, ndof) + Assemble(Vi_u, Vj_u, dNdy, dNdy, wJ_D33, ndof)
    K_vv = Assemble(Vi_v, Vj_v, dNdy, dNdy, wJ_D11, ndof) + Assemble(Vi_v, Vj_v, dNdx, dNdx, wJ_D33, ndof)
    wJ_D12 = wJ .* D12
    K_uv = Assemble(Vi_u, Vj_v, dNdx, dNdy, wJ_D12, ndof) + Assemble(Vi_u, Vj_v, dNdy, dNdx, wJ_D33, ndof)
    K_vu = Assemble(Vi_v, Vj_u, dNdy, dNdx, wJ_D12, ndof) + Assemble(Vi_v, Vj_u, dNdx, dNdy, wJ_D33, ndof)
    K_mech = K_uu + K_vv + K_uv + K_vu
    
    # 装配热-扩散载荷向量
    epsilon_0_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + β_n * Δsoc_n_elem[e] + β_p * Δsoc_p_elem[e]
    end
    
    coeff_u = zeros(Float64, ngs)
    coeff_v = zeros(Float64, ngs)
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        ε_0 = epsilon_0_elem[e]
        factor = E_eff / (1.0 - ν_eff^2) * ε_0 * (1.0 + ν_eff) * wJ[g]
        coeff_u[g] = factor
        coeff_v[g] = factor
    end
    F_u = Assemble1D(Vi_u, dNdx, coeff_u, ndof)
    F_v = Assemble1D(Vi_v, dNdy, coeff_v, ndof)
    F_mech = F_u + F_v
    
    # 施加边界条件
    bc_nodes = Dict{Int, Symbol}()
    is_inner, is_outer = identify_boundary_nodes(mesh, case.param, case.opt)
    for i in 1:nnode
        if is_inner[i] || is_outer[i]
            bc_nodes[i] = :fixed_xy
        end
    end
    
    penalty = 1e12
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_x
            dof = 2 * node - 1
            K_mech[dof, dof] += penalty
            F_mech[dof] = 0.0
        elseif bc_type == :fixed_y
            dof = 2 * node
            K_mech[dof, dof] += penalty
            F_mech[dof] = 0.0
        elseif bc_type == :fixed_xy
            dof_x = 2 * node - 1
            dof_y = 2 * node
            K_mech[dof_x, dof_x] += penalty
            K_mech[dof_y, dof_y] += penalty
            F_mech[dof_x] = 0.0
            F_mech[dof_y] = 0.0
        end
    end
    
    # 求解位移场
    U_M = try
        K_mech \ F_mech
    catch e
        @warn "Mechanical solve failed, using zero displacement" e
        zeros(Float64, ndof)
    end
    
    # 恢复应力场
    σ_xx = zeros(Float64, ne)
    σ_yy = zeros(Float64, ne)
    σ_xy = zeros(Float64, ne)
    σ_vm = zeros(Float64, ne)
    σ_thermal_vm = zeros(Float64, ne)
    σ_diffusion_vm = zeros(Float64, ne)
    epsilon_thermal_elem = zeros(Float64, ne)
    epsilon_diffusion_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        ε_thermal = α_eff * dT_elem[e]
        ε_diff_n = β_n * Δsoc_n_elem[e]
        ε_diff_p = β_p * Δsoc_p_elem[e]
        epsilon_thermal_elem[e] = ε_thermal
        epsilon_diffusion_elem[e] = ε_diff_n + ε_diff_p
        epsilon_0_elem[e] = ε_thermal + ε_diff_n + ε_diff_p
    end
    
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        grads = q4_center_gradients(mesh.node, nodes)
        grads === nothing && continue
        dNdx_local, dNdy_local, _ = grads
        ε_xx = sum(dNdx_local[i] * U_M[2 * nodes[i] - 1] for i in 1:4)
        ε_yy = sum(dNdy_local[i] * U_M[2 * nodes[i]] for i in 1:4)
        γ_xy = sum(dNdy_local[i] * U_M[2 * nodes[i] - 1] + dNdx_local[i] * U_M[2 * nodes[i]] for i in 1:4)
        ε_0 = epsilon_0_elem[e]
        ε_elastic_xx = ε_xx - ε_0
        ε_elastic_yy = ε_yy - ε_0
        ε_elastic_xy = γ_xy
        factor = E_eff / (1.0 - ν_eff^2)
        σ_xx[e] = factor * (ε_elastic_xx + ν_eff * ε_elastic_yy)
        σ_yy[e] = factor * (ε_elastic_yy + ν_eff * ε_elastic_xx)
        σ_xy[e] = factor * (1.0 - ν_eff) / 2.0 * ε_elastic_xy
        σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 - σ_xx[e] * σ_yy[e] + 3.0 * σ_xy[e]^2)
        if abs(ε_0) > 1e-15
            ratio_thermal = epsilon_thermal_elem[e] / ε_0
            ratio_diffusion = epsilon_diffusion_elem[e] / ε_0
            σ_thermal_vm[e] = abs(ratio_thermal) * σ_vm[e]
            σ_diffusion_vm[e] = abs(ratio_diffusion) * σ_vm[e]
        else
            σ_thermal_vm[e] = 0.0
            σ_diffusion_vm[e] = 0.0
        end
    end
    
    # 写入结果（转换为有量纲）
    L_ref = param.scale.L  # 使用统一长度尺度
    new_variables = copy(variables)
    new_variables["displacement x"] = U_M[1:2:end] .* L_ref
    new_variables["displacement y"] = U_M[2:2:end] .* L_ref
    new_variables["diffusion stress xx"] = σ_xx
    new_variables["diffusion stress yy"] = σ_yy
    new_variables["diffusion stress xy"] = σ_xy
    new_variables["diffusion stress vonMises"] = σ_vm
    new_variables["thermal stress vonMises"] = σ_thermal_vm
    new_variables["diffusion stress vonMises only"] = σ_diffusion_vm
    
    return new_variables
end
