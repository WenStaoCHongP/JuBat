function PostProcessing(case::Case, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    result = Dict()
    result["time [s]"]= variables["time"][1:v] * case.param.scale.t0
    result["cell voltage [V]"]= variables["cell voltage"][1:v] * case.param.scale.phi
    result["cell current [A]"]= variables["cell current"][1:v] * case.param_dim.cell.I1C
    result["temperature [K]"] = variables["temperature"][1:v] * case.param_dim.scale.T_ref
    result["negative particle center radial stress[Pa]"] = variables["negative particle center radial stress"][:,1:v] * case.param.scale.E_n
    result["positive particle center radial stress[Pa]"] = variables["positive particle center radial stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface tangential stress[Pa]"] = variables["negative particle surface tangential stress"][:,1:v] * case.param.scale.E_n
    result["positive particle surface tangential stress[Pa]"] = variables["positive particle surface tangential stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface displacement[m]"] = variables["negative particle surface displacement"][:,1:v] * case.param.scale.r0
    result["positive particle surface displacement[m]"] = variables["positive particle surface displacement"][:,1:v] * case.param.scale.r0
    if case.opt.model == "SPM"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
    elseif case.opt.model == "SPMe"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][:,1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][:,1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][:,1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][:,1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
        result["negative electrode potential [V]"] = variables["negative electrode potential"][:,1:v] * case.param.scale.phi
        result["positive electrode potential [V]"] = variables["positive electrode potential"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in negative electrode [V]"] = variables["electrolyte potential in negative electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in positive electrode [V]"] = variables["electrolyte potential in positive electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential [V]"] = variables["electrolyte potential"][:,1:v] * case.param.scale.phi
        result["negative electrode open circuit potential [V]"] = variables["negative electrode open circuit potential"] * case.param.scale.phi
        result["positive electrode open circuit potential [V]"] = variables["positive electrode open circuit potential"] * case.param.scale.phi
        result["negative electrode interfacial current density [A/m^2]"]  = variables["negative electrode interfacial current density"] * case.param.scale.j
        result["positive electrode interfacial current density [A/m^2]"]  = variables["positive electrode interfacial current density"] * case.param.scale.j
    end

    if case.opt.thermalmodel == "lumped"
        result["thermal lumped internal heat [W/m^3]"] = vec(variables["thermal lumped internal heat"][1, 1:v]) * case.param.scale.q
    elseif case.opt.thermalmodel == "distributed2D"
        # ── 热源物理单位还原 ──
        result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_s_NE [W/m3]"] = variables["thermal2D q_ohm_s_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_e_NE [W/m3]"] = variables["thermal2D q_ohm_e_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_SP [W/m3]"] = variables["thermal2D q_sp"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rxn_PE [W/m3]"] = variables["thermal2D q_rxn_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rev_PE [W/m3]"] = variables["thermal2D q_rev_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_s_PE [W/m3]"] = variables["thermal2D q_ohm_s_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_e_PE [W/m3]"] = variables["thermal2D q_ohm_e_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_PCC [W/m3]"] = variables["thermal2D q_pcc"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_NCC [W/m3]"] = variables["thermal2D q_ncc"][:, 1:v] * case.param.scale.q

        # ── 温度还原 ──
        result["thermal2D temperature at nodes [K]"] = variables["thermal2D temperature at nodes"][:, 1:v] * case.param_dim.scale.T_ref

        # ── 单元温度（从节点平均计算）──
        mesh_th = case.mesh["thermal2D"]
        ne = size(mesh_th.element, 1)
        n_t = v
        Tref = case.param_dim.scale.T_ref
        T_elem_hist = zeros(Float64, ne, n_t)
        for ti in 1:n_t
            T_nodes_t = variables["thermal2D temperature at nodes"][:, ti]
            T_elem_hist[:, ti] = element_nodal_mean(mesh_th, T_nodes_t)
        end
        result["thermal2D temperature [K]"] = T_elem_hist .* Tref

        # ── 热源场与总热源 ──
        result["heat_source_fields"] = variables["heat_source_fields"][:, 1:v]
        result["total heat source [W]"] = vec(variables["total heat source"][1, 1:v])

        # ── 单元级变量（无量纲直传）──
        for key in ["thermal2D element current", "thermal2D eta_n_e", "thermal2D eta_p_e",
                    "thermal2D element soc_n", "thermal2D element soc_p",
                    "thermal2D element voltages", "thermal2D element OCV",
                    "thermal2D dUdT_n_e", "thermal2D dUdT_p_e"]
            result[key] = variables[key][:, 1:v]
        end

        # ── 截止与激活信息 ──
        result["thermal2D active_mask"] = variables["thermal2D active_mask"][:, 1:v]
        result["thermal2D n_cutoff_elements"] = variables["thermal2D n_cutoff_elements"][1, 1:v]
        result["thermal2D nearest_cutoff_element"] = variables["thermal2D nearest_cutoff_element"][1, 1:v]
        result["thermal2D nearest_cutoff_ocv"] = variables["thermal2D nearest_cutoff_ocv"][1, 1:v]
        result["thermal2D margin_to_cutoff"] = variables["thermal2D margin_to_cutoff"][1, 1:v]
    end
    if case.opt.czm_enabled == true
        result["collapse_approx"] = "phi_perfect_bond"
        if case.opt.czm_winding_prestress
            result["winding prestress"] = vec(variables["winding prestress"][1, 1:v])
        end
        result["czm D_max"] = vec(variables["czm D_max"][1, 1:v])
        result["czm D_mean"] = vec(variables["czm D_mean"][1, 1:v])
        # δ_max_n 存储于分离空间（scale.δ_czm 归一，重设计 v2；修正原误用 scale.L）
        result["czm δ_max_n [m]"] = vec(variables["czm δ_max_n"][1, 1:v]) * case.param_dim.scale.δ_czm
        result["czm δ_mean_n [m]"] = vec(variables["czm δ_mean_n"][1, 1:v]) * case.param_dim.scale.δ_czm
        result["czm n_fractured"] = vec(variables["czm n_fractured"][1, 1:v])
        result["czm damage [0-1]"] = variables["czm damage"][:, 1:v]
        result["czm displacement x [m]"] = variables["czm displacement x"][:, 1:v] * case.param.scale.L
        result["czm displacement y [m]"] = variables["czm displacement y"][:, 1:v] * case.param.scale.L
        # 牵引力以 scale.σ_czm 归一（重设计 v2；修正原误用颗粒化学压尺度 E_n/E_p）
        result["czm traction normal [Pa]"] = variables["czm traction normal"][:, 1:v] * case.param.scale.σ_czm
        result["czm traction tangent [Pa]"] = variables["czm traction tangent"][:, 1:v] * case.param.scale.σ_czm
        # 分离位移以 scale.δ_czm 归一（重设计 v2；修正原误用 scale.r0 颗粒半径尺度）
        result["czm separation normal [m]"] = variables["czm separation normal"][:, 1:v] * case.param.scale.δ_czm
        result["czm separation tangent [m]"] = variables["czm separation tangent"][:, 1:v] * case.param.scale.δ_czm
    end
    return result
end
