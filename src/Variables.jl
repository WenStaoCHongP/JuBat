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

    if "temperature" in collect(keys(case.index))
        variables["temperature"] = zeros(Float64, length(case.index["temperature"]), num)
    else
        variables["temperature"] = zeros(Float64, 1, num)
    end

    variables["T_nodes"] = zeros(Float64, 0, num)
    variables["T_prev"] = zeros(Float64, 0, num)
    variables["heat_source_fields"] = zeros(Float64, 0, num)
    if case.opt.model == "SPMe" && case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
        ne = size(case.mesh["thermal2D"].element, 1)
        nT = case.mesh["thermal2D"].nlen
        # 覆盖热源历史为逐单元尺寸
        variables["heat_source_fields"] = zeros(Float64, ne, num)
        variables["thermal2D element current"] = zeros(Float64, ne, num)
        variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p"] = zeros(Float64, ne, num)
        variables["thermal2D element voltages"] = zeros(Float64, ne, num)
        variables["thermal2D element thermal stress"] = zeros(Float64, ne, num)
        variables["thermal2D element diffusion stress"] = zeros(Float64, ne, num)
        variables["thermal2D element total stress"] = zeros(Float64, ne, num)
        variables["thermal2D element diffusion strain"] = zeros(Float64, ne, num)
        variables["thermal2D element thermal strain"] = zeros(Float64, ne, num)
        variables["thermal2D temperature"] = zeros(Float64, nT, num)
    end
    
    if hasproperty(case.opt, :simple_thermal_coupling) && 
       case.opt.simple_thermal_coupling &&
       case.opt.thermalmodel == "distributed2D" && 
       haskey(case.mesh, "thermal2D")
        nT = case.mesh["thermal2D"].nlen
        # 温度场历史
        variables["thermal2D temperature"] = zeros(Float64, nT, num)
        # 全局变量
        variables["average temperature"] = zeros(Float64, 1, num)
        variables["total heat source"] = zeros(Float64, 1, num)
    end

    if case.opt.thermalmodel == "lumped"
        variables["thermal lumped internal heat"] = zeros(Float64, 1, num)
        variables["thermal lumped internal heat [W]"] = zeros(Float64, 1, num)
    end

    if case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
        nT = case.mesh["thermal2D"].nlen
        variables["thermal2D displacement x"] = zeros(Float64, nT, num)
        variables["thermal2D displacement y"] = zeros(Float64, nT, num)
    end

    return variables
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
    
    # 仅更新历史中已存在的键，避免临时/额外键尺寸不匹配
    for k in keys(variables_hist)
        # 支持标量历史（1行）和向量历史（n行）
        if isa(variables_hist[k], Array{Float64})
            # 目标历史为矩阵 (nrows x num)
            nrows = size(variables_hist[k], 1)
            if haskey(variables, k)
                val = variables[k]
                if isa(val, Array{Float64})
                    # 取列向量/首列并截断/填充
                    col = ndims(val) == 1 ? val : val[:,1]
                    if length(col) == nrows
                        variables_hist[k][:, v] = col
                    elseif nrows == 1 && length(col) >= 1
                        variables_hist[k][1, v] = col[1]
                    end
                elseif isa(val, Float64)
                    if nrows == 1
                        variables_hist[k][1, v] = val
                    end
                end
            end
        elseif isa(variables_hist[k], Float64)
            # 历史为标量
            if haskey(variables, k)
                val = variables[k]
                variables_hist[k] = isa(val, Float64) ? val : (isa(val, Array{Float64}) ? (ndims(val) == 1 ? (length(val) > 0 ? val[1] : variables_hist[k]) : (size(val,1) > 0 ? val[1,1] : variables_hist[k])) : variables_hist[k])
            end
        end
    end
    return variables_hist
end