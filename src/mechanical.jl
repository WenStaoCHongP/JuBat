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
