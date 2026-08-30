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
# 层分辨宏观应力（热 + 扩散）：共享恢复核 + 耦合在线收割 + 固体按需求解
# ========================================================================

"""
    recover_bulk_stress(node, element, material_type, u, ε0, param)

逐单元恢复层分辨宏观应力（σ_czm 归一空间）。平面应力 σ = D·(ε − ε0·[1,1,0])，
D 由 `moduli_of(param, material_type[e])` 逐层给出；`u` 为 L 归一化位移
（按 [u1x,u1y,u2x,u2y,...] 排布），`ε0` 为逐单元特征应变。
"""
function recover_bulk_stress(node, element, material_type, u, ε0, param)
    ne = size(element, 1)
    sigma_xx = zeros(Float64, ne)
    sigma_yy = zeros(Float64, ne)
    sigma_xy = zeros(Float64, ne)
    sigma_vm = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nodes = element[e, :]
        dNdx, dNdy, _ = q4_center_gradients(node, nodes)
        eps_xx = sum(dNdx[i] * u[2 * nodes[i] - 1] for i in 1:4)
        eps_yy = sum(dNdy[i] * u[2 * nodes[i]] for i in 1:4)
        gam_xy = sum(dNdy[i] * u[2 * nodes[i] - 1] + dNdx[i] * u[2 * nodes[i]] for i in 1:4)
        E_e, nu_e = moduli_of(param, material_type[e])
        factor = E_e / (1.0 - nu_e^2)
        e_xx = eps_xx - ε0[e]
        e_yy = eps_yy - ε0[e]
        sigma_xx[e] = factor * (e_xx + nu_e * e_yy)
        sigma_yy[e] = factor * (e_yy + nu_e * e_xx)
        sigma_xy[e] = factor * (1.0 - nu_e) / 2.0 * gam_xy
        sigma_vm[e] = sqrt(sigma_xx[e]^2 + sigma_yy[e]^2 - sigma_xx[e] * sigma_yy[e] + 3.0 * sigma_xy[e]^2)
    end
    return sigma_xx, sigma_yy, sigma_xy, sigma_vm
end

"""
    macro_eigenstrain(case, variables, T_nodes) -> Vector{Float64}

逐力学体单元特征应变，与在线 CZM 热-化学载荷同源（`eigenstrain_of`：
alphaT(mt)·dT + Ω(mt)/3·Δsoc(mt)，α/β 分层分辨率）；dT 按父热单元，
Δsoc 按材料分发（集流体/隔膜只有热应变，电极膨胀只作用于本层涂层）。
"""
function macro_eigenstrain(case, variables, T_nodes)
    param = case.param
    strain_in = compute_czm_strain_inputs(case, variables, T_nodes)
    mt = case.czm_mesh.czm_submesh.material_type
    return [eigenstrain_of(param, mt[e], strain_in.dT_czm[e],
                           strain_in.Δsoc_n_czm[e], strain_in.Δsoc_p_czm[e])
            for e in eachindex(mt)]
end

"""
    compute_macro_stress(case, variables, T_nodes)

以当前已收敛 `case.czm_layout.u_prev` 和同一时间层载荷恢复四个层分辨应力分量。
返回归一化到 `scale.σ_czm` 的 NamedTuple；调用方负责决定历史记录时刻。
"""
function compute_macro_stress(case, variables, T_nodes)
    czm_mesh = case.czm_mesh
    ε0 = macro_eigenstrain(case, variables, T_nodes)
    sigma_xx, sigma_yy, sigma_xy, sigma_vm = recover_bulk_stress(
        czm_mesh.node, czm_mesh.bulk_element, czm_mesh.czm_submesh.material_type,
        case.czm_layout.u_prev, ε0, case.param)
    return (xx=sigma_xx, yy=sigma_yy, xy=sigma_xy, vonMises=sigma_vm)
end

function write_macro_stress!(variables_hist, v, stress)
    variables_hist["diffusion stress xx"][:, v] = stress.xx
    variables_hist["diffusion stress yy"][:, v] = stress.yy
    variables_hist["diffusion stress xy"][:, v] = stress.xy
    variables_hist["diffusion stress vonMises"][:, v] = stress.vonMises
    return variables_hist
end

"""
    export_macro_stress(case, variables, variables_hist, v, T_nodes)

恢复当前已收敛层分辨应力并写入第 `v` 个历史列。Solve 主循环另行保存最近一次
实际 CZM 更新的恢复结果，使非更新步保持该有效状态而不是保留预分配零。
非线性路径不分配历史键，此处按既有门控原样返回。
"""
function export_macro_stress(case, variables, variables_hist, v, T_nodes)
    (case.opt.czm_enabled && case.czm_mesh !== nothing) || return variables_hist
    haskey(variables_hist, "diffusion stress xx") || return variables_hist
    stress = compute_macro_stress(case, variables, T_nodes)
    write_macro_stress!(variables_hist, v, stress)
    return variables_hist
end

"""
    thermal_diffusion_stress_2D(case, variables)

仅固体力学流程（opt.czm_enabled=false）的层分辨宏观应力工具函数：在
`czm_submesh.mesh_bonded`（Φ 合并、无内聚力单元）上按逐层刚度（`moduli_of`）
与逐层特征应变求解线性弹性静力平衡，返回有量纲结果。

# 输入
- `case.czm_mesh`：需由 `jellyroll_collector_seed_mesh(czm_enabled=true)` +
  `create_czm_mesh` 构建（`opt.czm_enabled` 可保持 false，在线 CZM 路径不受影响）
- `variables["T_nodes"]`：归一化节点温度（T/T_ref；向量或历史矩阵取末列）
- `variables["thermal2D element soc_n/soc_p"]`：归一化化学计量比
- 边界：外圈固定；内圈按 `opt.czm_fix_inner`（默认 true = 内外均固定）

# 输出（返回 copy(variables) 增键，均为有量纲）
- `"diffusion stress xx/yy/xy/vonMises [Pa]"`（σ_czm 空间恢复后 ×scale.σ_czm）
- `"displacement x/y [m]"`（L 归一化位移 ×scale.L）
"""
function thermal_diffusion_stress_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    # === 入口断言：参数集必须定义 E_coat 才能启用宏观力学 ===
    @assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0 "宏观力学分析需要 PE/NE.E_coat > 0；当前参数集未定义极片模量（E_coat=0）。请在参数文件中补全 PE.E_coat/PE.nu_coat/NE.E_coat/NE.nu_coat，或禁用 mechanicalmodel=\"full\"。"

    param = case.param
    submesh = case.czm_mesh.czm_submesh
    mesh = submesh.mesh_bonded
    T_nodes = variables["T_nodes"]
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    ε0 = macro_eigenstrain(case, variables, T_nodes)

    ne = size(mesh.element, 1)
    nnode = mesh.nlen
    ndof = 2 * nnode
    dNdx = mesh.gs.dNidx[:, 1:4]
    dNdy = mesh.gs.dNidx[:, 5:8]
    wJ = mesh.gs.weight .* mesh.gs.detJ
    ele_of_gp = mesh.gs.ele
    ngs = length(wJ)

    # 逐高斯点平面应力矩阵系数（材料按 material_type 逐层）
    D11 = Vector{Float64}(undef, ngs)
    D12 = Vector{Float64}(undef, ngs)
    D33 = Vector{Float64}(undef, ngs)
    @inbounds for g in 1:ngs
        E_e, ν_e = moduli_of(param, submesh.material_type[ele_of_gp[g]])
        D11[g] = E_e / (1.0 - ν_e^2)
        D12[g] = D11[g] * ν_e
        D33[g] = D11[g] * (1.0 - ν_e) / 2.0
    end

    Vi_u = zeros(Int64, ngs, 4)
    Vi_v = zeros(Int64, ngs, 4)
    Vj_u = zeros(Int64, ngs, 4)
    Vj_v = zeros(Int64, ngs, 4)
    @inbounds for g in 1:ngs
        for i in 1:4
            node = mesh.element[ele_of_gp[g], i]
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

    coeff_u = Vector{Float64}(undef, ngs)
    coeff_v = Vector{Float64}(undef, ngs)
    @inbounds for g in 1:ngs
        e = ele_of_gp[g]
        E_e, ν_e = moduli_of(param, submesh.material_type[e])
        factor = E_e / (1.0 - ν_e^2) * ε0[e] * (1.0 + ν_e) * wJ[g]
        coeff_u[g] = factor
        coeff_v[g] = factor
    end
    F_mech = Assemble1D(Vi_u, dNdx, coeff_u, ndof) + Assemble1D(Vi_v, dNdy, coeff_v, ndof)

    # 边界：外圈固定；内圈按 opt.czm_fix_inner（默认 true = 内外均固定），相对罚
    is_inner, is_outer = identify_boundary_nodes(mesh, param, case.opt)
    penalty = 1e6 * maximum(abs, diag(K_mech))
    @inbounds for i in 1:nnode
        if is_outer[i] || (case.opt.czm_fix_inner && is_inner[i])
            K_mech[2 * i - 1, 2 * i - 1] += penalty
            K_mech[2 * i, 2 * i] += penalty
            F_mech[2 * i - 1] = 0.0
            F_mech[2 * i] = 0.0
        end
    end

    U = K_mech \ F_mech

    sigma_xx, sigma_yy, sigma_xy, sigma_vm = recover_bulk_stress(
        mesh.node, mesh.element, submesh.material_type, U, ε0, param)

    L_ref = param.scale.L
    σ_czm = param.scale.σ_czm
    new_variables = copy(variables)
    new_variables["diffusion stress xx [Pa]"] = sigma_xx .* σ_czm
    new_variables["diffusion stress yy [Pa]"] = sigma_yy .* σ_czm
    new_variables["diffusion stress xy [Pa]"] = sigma_xy .* σ_czm
    new_variables["diffusion stress vonMises [Pa]"] = sigma_vm .* σ_czm
    new_variables["displacement x [m]"] = U[1:2:end] .* L_ref
    new_variables["displacement y [m]"] = U[2:2:end] .* L_ref
    return new_variables
end
