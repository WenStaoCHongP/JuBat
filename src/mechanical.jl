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
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is required for 2D diffusion stress"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "diffusion_stress_2D requires Q4 mesh"
    
    param = case.param
    Tref = param.scale.T_ref
    T0 = hasproperty(param.cell, :T0) ? param.cell.T0 : 298.0 / Tref
    
    # 提取温度场和SOC分布
    T_nodes = haskey(variables, "T_nodes") ? variables["T_nodes"] : fill(T0, mesh.nlen)
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    soc_n_elem = variables["thermal2D element soc_n"]
    soc_p_elem = variables["thermal2D element soc_p"]
    soc_ref_n = param.NE.cs0 
    soc_ref_p = param.PE.cs0
    # 获取材料参数
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    ν_eff = (param.NE.nu * param.NE.thickness + param.PE.nu * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    α_eff = (param.NE.alphaT * param.NE.thickness + param.PE.alphaT * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    β_n = param.NE.Omega / 3.0 
    β_p = param.PE.Omega / 3.0 

    # 计算单元级别的温度和SOC
    ne = size(mesh.element, 1)
    T_elem = zeros(Float64, ne)
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        T_elem[e] = sum(T_nodes[nodes]) / length(nodes)
        dT_elem[e] = T_elem[e] - T0
        Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
        Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
    end

    # 装配力学刚度矩阵
    K_mech = _assemble_mechanical_stiffness_2D(mesh, E_eff, ν_eff)
    
    # 装配热-扩散载荷向量
    F_mech = _assemble_thermal_diffusion_load_2D(mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    
    # 施加边界条件
    K_mech, F_mech = _apply_mechanical_BC_2D(K_mech, F_mech, mesh, case)
    
    # 求解位移场
    U_M = _solve_mechanical_displacement_2D(K_mech, F_mech, mesh.nlen)
    
    # 恢复应力场
    σ_xx, σ_yy, σ_xy, σ_vm = _recover_stress_2D(U_M, mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    
    # 写入结果（转换为有量纲）
    L_ref = hasproperty(param.scale, :L_th) ? param.scale.L_th : 1.0
    _write_mechanical_results!(variables, U_M, σ_xx, σ_yy, σ_xy, σ_vm, L_ref)
    
    return variables
end


"""装配2D力学刚度矩阵"""
function _assemble_mechanical_stiffness_2D(mesh, E_eff, ν_eff)
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
    # 平面应力：D = E/(1-ν²) * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    D11 = zeros(Float64, ngs)
    D12 = zeros(Float64, ngs)
    D33 = zeros(Float64, ngs)
    E = E_eff
    ν = ν_eff
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        D11[g] = E / (1.0 - ν^2)
        D12[g] = (E / (1.0 - ν^2)) * ν
        D33[g] = (E / (1.0 - ν^2)) * (1.0 - ν) / 2.0
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
    K_M = K_uu + K_vv + K_uv + K_vu
    
    return K_M
end

"""装配热-扩散载荷向量"""
function _assemble_thermal_diffusion_load_2D(mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    nnode = mesh.nlen
    ndof = 2 * nnode
    E = E_eff
    ν = ν_eff
    # 计算每个单元的初始应变 ε_0 = α*ΔT + β_c*c_s_max*ΔSOC
    ne = length(dT_elem)
    epsilon_0_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + β_n * Δsoc_n_elem[e] + β_p * Δsoc_p_elem[e]
    end
    
    # 高斯积分点数据
    dNdx = mesh.gs.dNidx[:, 1:4]
    dNdy = mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele_of_gp = mesh.gs.ele
    
    ngs = length(wJ)
    
    # 计算载荷系数
    # F = ∫ B^T D ε_0 dΩ
    # 其中 ε_0 = [ε_0, ε_0, 0]^T
    # D * ε_0 = E/(1-ν²) * [ε_0*(1+ν), ε_0*(1+ν), 0]^T
    
    coeff_u = zeros(Float64, ngs)
    coeff_v = zeros(Float64, ngs)
    
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        ε_0 = epsilon_0_elem[e]
        factor = E / (1.0 - ν^2) * ε_0 * (1.0 + ν) * wJ[g]
        coeff_u[g] = factor
        coeff_v[g] = factor
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
    F_M = zeros(Float64, ndof)
    
    # F_u = ∫ dNi/dx * D11 * ε_0 dΩ
    F_u = Assemble1D(Vi_u, dNdx, coeff_u, ndof)
    
    # F_v = ∫ dNi/dy * D11 * ε_0 dΩ
    F_v = Assemble1D(Vi_v, dNdy, coeff_v, ndof)
    
    F_M = F_u + F_v
    
    return F_M
end

"""施加力学边界条件"""
function _apply_mechanical_BC_2D(K_M, F_M, mesh, case)
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
            K_M[dof, dof] += penalty
            F_M[dof] = 0.0
        elseif bc_type == :fixed_y
            # 固定y方向位移
            dof = 2 * node
            K_M[dof, dof] += penalty
            F_M[dof] = 0.0
        elseif bc_type == :fixed_xy
            # 固定x和y方向位移
            dof_x = 2 * node - 1
            dof_y = 2 * node
            K_M[dof_x, dof_x] += penalty
            K_M[dof_y, dof_y] += penalty
            F_M[dof_x] = 0.0
            F_M[dof_y] = 0.0
        end
    end
    
    return K_M, F_M
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
function _solve_mechanical_displacement_2D(K_M::SparseMatrixCSC{Float64,Int}, F_M::Vector{Float64}, nnode::Int)
    ndof = 2 * nnode
    U_M = zeros(Float64, ndof)
    # 求解线性方程组 K*U = F
        U_M = try
            K_M \ F_M
    catch e
        @warn "Mechanical solve failed, using zero displacement" e
        U_M = zeros(Float64, ndof)
    end
    
    return U_M
end

"""恢复应力场"""
function _recover_stress_2D(U_M, mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    ne = size(mesh.element, 1)
    E = E_eff
    ν = ν_eff
    # 初始化应力数组
    σ_xx = zeros(Float64, ne)
    σ_yy = zeros(Float64, ne)
    σ_xy = zeros(Float64, ne)
    σ_vm = zeros(Float64, ne)
    
    # 计算每个单元的初始应变
    epsilon_0_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + β_n * Δsoc_n_elem[e] + β_p * Δsoc_p_elem[e]
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
        ε_xx = sum(dNdx[i] * U_M[2*nodes[i]-1] for i in 1:4)
        ε_yy = sum(dNdy[i] * U_M[2*nodes[i]] for i in 1:4)
        γ_xy = sum(dNdy[i] * U_M[2*nodes[i]-1] + dNdx[i] * U_M[2*nodes[i]] for i in 1:4)
        
        # 弹性应变
        ε_0 = epsilon_0_elem[e]
        ε_elastic_xx = ε_xx - ε_0
        ε_elastic_yy = ε_yy - ε_0
        ε_elastic_xy = γ_xy  # 剪应变没有初始值
        
        # 计算应力（平面应力）
        factor = E / (1.0 - ν^2)
        
        σ_xx[e] = factor * (ε_elastic_xx + ν * ε_elastic_yy)
        σ_yy[e] = factor * (ε_elastic_yy + ν * ε_elastic_xx)
        σ_xy[e] = factor * (1.0 - ν) / 2.0 * ε_elastic_xy
        
        # Von Mises应力
        σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 - σ_xx[e]*σ_yy[e] + 3.0*σ_xy[e]^2)
    end
    
    return σ_xx, σ_yy, σ_xy, σ_vm
end

"""写入力学计算结果"""
function _write_mechanical_results!(variables, U_M, σ_xx, σ_yy, σ_xy, σ_vm, L_ref)
    nnode = length(U_M) ÷ 2
    
    # 提取位移
    u_x = U_M[1:2:end]
    u_y = U_M[2:2:end]
    
    # 写入变量（无量纲 -> 有量纲）
    variables["displacement x"] = u_x .* L_ref
    variables["displacement y"] = u_y .* L_ref
    variables["diffusion stress xx"] = σ_xx
    variables["diffusion stress yy"] = σ_yy
    variables["diffusion stress xy"] = σ_xy
    variables["diffusion stress vonMises"] = σ_vm
    
    return nothing
end
