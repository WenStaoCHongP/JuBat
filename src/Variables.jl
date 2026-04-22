function StandardVariables(case::Case, num::Int64)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        Nn = 1
        Np = 1
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
    elseif case.opt.model == "P2D"
        Nn = case.mesh["negative electrode"].nlen
        Np = case.mesh["positive electrode"].nlen
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
    end

    variables = Dict{String, Union{Array{Float64}, Float64}}(
        "negative particle lithium concentration" => zeros(Float64, Nrn, num),
        "positive particle lithium concentration" => zeros(Float64, Nrp, num),
        "negative particle averaged lithium concentration" => zeros(Float64, Nn, num),
        "positive particle averaged lithium concentration" => zeros(Float64, Np, num),
        "negative particle surface lithium concentration" => zeros(Float64, Nn, num),
        "positive particle surface lithium concentration" => zeros(Float64, Np, num),    
        "negative electrode porosity" => zeros(Float64, Nn, num),
        "positive electrode porosity" => zeros(Float64, Np, num),
        "negative electrode temperature" => zeros(Float64, Nn, num),
        "positive electrode temperature" => zeros(Float64, Np, num),  
        "negative electrode exchange current density" => zeros(Float64, Nn, num),
        "positive electrode exchange current density" => zeros(Float64, Np, num), 
        "negative electrode interfacial current density" => zeros(Float64, Nn, num),
        "positive electrode interfacial current density" => zeros(Float64, Np, num), 
        "negative electrode overpotential" => zeros(Float64, Nn, num),
        "positive electrode overpotential" => zeros(Float64, Np, num), 
        "negative electrode open circuit potential" => zeros(Float64, Nn, num),
        "positive electrode open circuit potential" => zeros(Float64, Np, num),
        "negative particle center radial stress" => zeros(Float64, Nn, num),
        "positive particle center radial stress" => zeros(Float64, Np, num),
        "negative particle surface tangential stress" => zeros(Float64, Nn, num),
        "positive particle surface tangential stress" => zeros(Float64, Np, num),
        "negative particle surface displacement" => zeros(Float64, Nn, num),
        "positive particle surface displacement" => zeros(Float64, Np, num),
        "negative particle concentration at gauss point" => zeros(Float64,  Nn* case.opt.Nrn* case.opt.gsorder, num),
        "positive particle concentration at gauss point" => zeros(Float64,  Np* case.opt.Nrp* case.opt.gsorder, num),
        "negative particle surface tangential stress at gauss point" => zeros(Float64, Ne_ngs, num),
        "positive particle surface tangential stress at gauss point" => zeros(Float64, Ne_pgs, num),
        "negative particle stress coupling diffusion coefficient" => zeros(Float64, Nn, num),
        "positive particle stress coupling diffusion coefficient" => zeros(Float64, Np, num),
        "cell voltage" => zeros(Float64, 1, num),
        "time" => zeros(Float64, 1, num),      
        "cell current" => zeros(Float64, 1, num),      
    )
    # additional variables for SPMe and P2D
    if case.opt.model == "SPMe" || case.opt.model == "P2D"
        Ne = case.mesh["electrolyte"].nlen
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        variables["electrolyte lithium concentration"] = zeros(Float64, Ne, num)
        variables["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n, num)
        variables["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p, num)
        variables["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp, num)
        variables["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs, num)
    end
    # extra variables for P2D
    if case.opt.model == "P2D"
        variables["negative electrode potential"] = zeros(Float64, Nn, num)
        variables["positive electrode potential"] = zeros(Float64, Np, num)
        variables["electrolyte potential"] = zeros(Float64, Ne, num)
        variables["electrolyte potential in negative electrode"] = zeros(Float64, Ne_n, num)
        variables["electrolyte potential in positive electrode"] = zeros(Float64, Ne_p, num)
        variables["electrolyte potential in separator"] = zeros(Float64, Ne_sp, num)
        variables["negative electrode interfacial current at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode interfacial current at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode overpotential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode overpotential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode exchange current density at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode exchange current density at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative particle surface lithium concentration at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive particle surface lithium concentration at Gauss point"] = zeros(Float64, Ne_pgs, num)
    end

    variables["temperature"] = zeros(Float64, length(case.index["temperature"]), num)
    
    if case.opt.thermalmodel == "lumped"
        variables["thermal lumped internal heat"] = zeros(Float64, 1, num)
    end

    if case.opt.thermalmodel == "distributed2D"
        ne = size(case.mesh["thermal2D"].element, 1)
        nT = case.mesh["thermal2D"].nlen
        variables["thermal2D temperature at nodes"] = zeros(Float64, nT, num)
        variables["thermal2D temperature history"] = zeros(Float64, ne, num)
        variables["heat_source_fields"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_sp"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_pcc"] = zeros(Float64, ne, num)
        variables["thermal2D q_ncc"] = zeros(Float64, ne, num)
        variables["thermal2D element current"] = zeros(Float64, ne, num)
        variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p"] = zeros(Float64, ne, num)
        variables["thermal2D element voltages"] = zeros(Float64, ne, num)
        variables["thermal2D element OCV"] = zeros(Float64, ne, num)
        variables["thermal2D active_mask"] = zeros(Float64, ne, num)
        variables["thermal2D n_cutoff_elements"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_element"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_ocv"] = zeros(Float64, 1, num)
        variables["thermal2D margin_to_cutoff"] = zeros(Float64, 1, num)
        variables["total heat source"] = zeros(Float64, 1, num)
    end
    if case.opt.czm_enabled == true
        variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
        variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
        variables["czm D_max"] = zeros(Float64, 1, num)
        variables["czm D_mean"] = zeros(Float64, 1, num)
        variables["czm δ_max_n"] = zeros(Float64, 1, num)
        variables["czm δ_mean_n"] = zeros(Float64, 1, num)
        variables["czm n_fractured"] = zeros(Float64, 1, num)
    end 
    return variables
end

"""
    create_element_workspace(case)

创建精简型单元工作区 Dict，仅包含 SPMe_element 调用链实际需要的键。
排除 distributed2D 专用的 ~30 个 thermal2D 键（由 CallModel_MultiSPMe
在单元循环外独立管理）。比 StandardVariables(case, 1) 减少约 60% 数组分配。
"""
function create_element_workspace(case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nn = 1; Np = 1  # SPMe 模式
    Ne_ngs = case.opt.Nn * case.opt.gsorder
    Ne_pgs = case.opt.Np * case.opt.gsorder

    ws = Dict{String, Union{Array{Float64}, Float64}}()

    # ── 状态提取键（case.index 对应，SPMe_variables! 从 yt 原位覆写）──
    ws["negative particle lithium concentration"] = zeros(Float64, Nrn, 1)
    ws["positive particle lithium concentration"] = zeros(Float64, Nrp, 1)
    ws["negative particle surface lithium concentration"] = zeros(Float64, Nn, 1)
    ws["positive particle surface lithium concentration"] = zeros(Float64, Np, 1)

    # ── SPMe 专用键 ──
    if case.opt.model == "SPMe"
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        ws["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n, 1)
        ws["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p, 1)
        ws["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp, 1)
        # Gauss 点计算结果（ElectrolyteDiffusion 读取）
        ws["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs, 1)
    end

    ws["temperature"] = 0.0

    # ── SPMe_variables! 计算结果键 ──
    ws["cell voltage"] = 0.0
    ws["time"] = 0.0
    ws["cell current"] = 0.0
    ws["negative electrode exchange current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode exchange current density"] = zeros(Float64, Np, 1)
    ws["negative electrode interfacial current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode interfacial current density"] = zeros(Float64, Np, 1)
    ws["negative electrode overpotential"] = zeros(Float64, Nn, 1)
    ws["positive electrode overpotential"] = zeros(Float64, Np, 1)
    ws["negative electrode open circuit potential"] = zeros(Float64, Nn, 1)
    ws["positive electrode open circuit potential"] = zeros(Float64, Np, 1)

    # ── 高斯点浓度键（SPMe_element L52-53 无条件读取）──
    # 注意：用 opt.Nrn/Nrp（单元数）而非 mesh.nlen（节点数）来计算高斯点数
    ws["negative particle concentration at gauss point"] = zeros(Float64, Nn * case.opt.Nrn * case.opt.gsorder, 1)
    ws["positive particle concentration at gauss point"] = zeros(Float64, Np * case.opt.Nrp * case.opt.gsorder, 1)

    # ── Mechanicaloutput 结果键（条件）──
    if case.opt.mechanicalmodel == "full"
        ws["negative particle center radial stress"] = zeros(Float64, Nn, 1)
        ws["positive particle center radial stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface tangential stress"] = zeros(Float64, Nn, 1)
        ws["positive particle surface tangential stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface displacement"] = zeros(Float64, Nn, 1)
        ws["positive particle surface displacement"] = zeros(Float64, Np, 1)
        ws["negative particle surface tangential stress at gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["positive particle surface tangential stress at gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["negative particle stress coupling diffusion coefficient"] = zeros(Float64, Nn, 1)
        ws["positive particle stress coupling diffusion coefficient"] = zeros(Float64, Np, 1)
    end

    # ── CZM 键（条件）──
    if case.opt.czm_enabled
        ws["negative electrode cohesive zone damage"] = zeros(Float64, Nn, 1)
        ws["positive electrode cohesive zone damage"] = zeros(Float64, Np, 1)
    end

    return ws
end

function Variable_update!(variables_hist::Dict{String, Union{Array{Float64},Float64}}, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    # 检查是否需要扩展数组（动态增长）
    for k in keys(variables_hist)
        if isa(variables_hist[k], Array{Float64}) && ndims(variables_hist[k]) == 2
            current_size = size(variables_hist[k], 2)
            if v > current_size
                # 需要扩展
                expansion_size = max(1000, current_size ÷ 2)
                @warn "时间步 $(v) 超过预分配 $(current_size)，扩展 $(expansion_size) 步（变量: $(k)）"
                n_rows = size(variables_hist[k], 1)
                new_cols = zeros(Float64, n_rows, expansion_size)
                variables_hist[k] = hcat(variables_hist[k], new_cols)
            end
        end
    end
    
    hist_keys = Set(keys(variables_hist))
    for (k, val) in pairs(variables)
        k in hist_keys || continue
        hist_val = variables_hist[k]
        if isa(hist_val, Array{Float64})
            nrows = size(hist_val, 1)
            if isa(val, Array{Float64})
                col = ndims(val) == 1 ? val : val[:, 1]
                if length(col) == nrows
                    hist_val[:, v] = col
                elseif nrows == 1 && !isempty(col)
                    hist_val[1, v] = col[1]
                end
            elseif isa(val, Float64) && nrows == 1
                hist_val[1, v] = val
            end
        elseif isa(hist_val, Float64)
            if isa(val, Float64)
                variables_hist[k] = val
            elseif isa(val, Array{Float64})
                col = ndims(val) == 1 ? val : val[:, 1]
                isempty(col) || (variables_hist[k] = col[1])
            end
        end
    end
    return variables_hist
end