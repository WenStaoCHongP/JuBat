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
    
    # 多尺度均匀化：将颗粒扩散应力映射到2D网格（如果启用力学模块）
    if case.opt.mechanical_enabled && haskey(case.mesh, "thermal2D")
        variables = homogenize_particle_stress_to_2D(case, variables)
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

"""
    homogenize_particle_stress_to_2D(case::Case, variables::Dict)

多尺度均匀化：将颗粒尺度扩散应力映射到宏观2D Jellyroll网格

# 理论背景
将微观颗粒应力（~10 μm尺度）均匀化到宏观单元应力（~1 mm尺度）

均匀化公式：
    σ_macro = Σ(f_layer · ε_active · σ_particle)

其中：
- f_layer = 层权重（负极/隔膜/正极的体积分数）
- ε_active = 活性材料体积分数（通常0.6-0.8）
- σ_particle = 颗粒表面切向应力

# 尺度桥接
```
微观（颗粒）         均匀化          宏观（单元）
~10 μm          =========>        ~1 mm
σ_particle                        σ_macro
```

# 方法选择
1. 如果有layer_weights：精细映射（推荐用于Jellyroll）
2. 否则：简单厚度加权平均

# 参数
- `case::Case` - JuBat案例对象
- `variables::Dict` - 包含颗粒应力的变量字典

# 返回
更新后的`variables`，新增：
- `"thermal2D element diffusion stress (homogenized)"` - 均匀化扩散应力 [Pa]
- `"thermal2D element total stress"` - 总应力（热应力+扩散应力） [Pa]

# 示例
```julia
# 在thermal_stress之后调用
variables = thermal_stress(case, variables)
variables = homogenize_particle_stress_to_2D(case, variables)

# 提取总应力
σ_total = variables["thermal2D element total stress"]
```
"""
function homogenize_particle_stress_to_2D(case::Case, variables::Dict)
    # 检查是否有2D网格
    if !haskey(case.mesh, "thermal2D")
        return variables
    end
    
    # 检查是否有颗粒应力
    if !haskey(variables, "negative particle surface tangential stress") ||
       !haskey(variables, "positive particle surface tangential stress")
        return variables
    end
    
    mesh_th = case.mesh["thermal2D"]
    ne_th = size(mesh_th.element, 1)
    param = case.param
    param_dim = case.param_dim
    
    # 提取颗粒应力（微观尺度）
    σ_particle_n = variables["negative particle surface tangential stress"]
    σ_particle_p = variables["positive particle surface tangential stress"]
    
    # 活性材料体积分数（考虑孔隙）
    ε_active_n = hasproperty(param_dim.NE, :epsilon_s) ? param_dim.NE.epsilon_s : 0.65
    ε_active_p = hasproperty(param_dim.PE, :epsilon_s) ? param_dim.PE.epsilon_s : 0.60
    
    # 有效应力 = 颗粒应力 × 体积分数
    σ_eff_n = isa(σ_particle_n, Number) ? σ_particle_n * ε_active_n : σ_particle_n .* ε_active_n
    σ_eff_p = isa(σ_particle_p, Number) ? σ_particle_p * ε_active_p : σ_particle_p .* ε_active_p
    
    # 获取layer_weights（如果可用）
    fks = nothing
    if haskey(variables, "thermal2D layer_weights")
        fks = variables["thermal2D layer_weights"]
    end
    
    # 均匀化到宏观单元
    σ_element = zeros(Float64, ne_th)
    
    if fks !== nothing && size(fks, 2) >= 3
        # 方法1：精细映射（使用layer_weights）
        @inbounds for e in 1:ne_th
            f_neg = fks[e, 1]  # 负极层体积分数
            f_pos = fks[e, 3]  # 正极层体积分数
            # 隔膜层（f_sep = fks[e,2]）不贡献扩散应力
            
            if isa(σ_eff_n, Number)
                # SPM/SPMe: 标量应力，所有颗粒相同
                σ_element[e] = f_neg * σ_eff_n + f_pos * σ_eff_p
            else
                # P2D: 向量应力，沿x轴分布
                # 简化：使用平均值（可改进为空间插值）
                σ_element[e] = f_neg * mean(σ_eff_n) + f_pos * mean(σ_eff_p)
            end
        end
        
    else
        # 方法2：简单厚度加权（没有layer_weights时的退化方案）
        t_n = param_dim.NE.thickness
        t_p = param_dim.PE.thickness
        t_total = t_n + t_p
        
        if isa(σ_eff_n, Number)
            σ_macro = (σ_eff_n * t_n + σ_eff_p * t_p) / t_total
            σ_element .= σ_macro
        else
            σ_macro = (mean(σ_eff_n) * t_n + mean(σ_eff_p) * t_p) / t_total
            σ_element .= σ_macro
        end
    end
    
    # 存储均匀化后的扩散应力
    variables["thermal2D element diffusion stress (homogenized)"] = σ_element
    
    # 叠加热应力（如果已计算）
    if haskey(variables, "thermal2D element thermal stress")
        σ_thermal = variables["thermal2D element thermal stress"]
        variables["thermal2D element total stress"] = σ_thermal .+ σ_element
    else
        # 如果只有扩散应力
        variables["thermal2D element total stress"] = σ_element
    end
    
    return variables
end
