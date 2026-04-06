function ModelInitialisation(case::Case)
    if isempty(case.opt.y0)
        if case.opt.model == "SPM"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            y0 = [csn0;  csp0]
        elseif case.opt.model == "SPMe"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) *  case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) *  case.param.PE.cs0
            Ne = case.mesh["electrolyte"].nlen
            ce0 = ones(Float64, Ne, 1) *  case.param.EL.ce0
            y0 = [csn0;  csp0; ce0]
        elseif case.opt.model == "P2D"
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Ne = case.mesh["electrolyte"].nlen
            Nn = case.mesh["negative electrode"].nlen
            Np = case.mesh["positive electrode"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            ce0 = ones(Float64, Ne, 1) * case.param.EL.ce0
            phie0 = - ones(Float64, Ne, 1) * case.param.NE.U(case.param.NE.cs0)
            phis_p =  ones(Float64, Nn, 1) * case.param.PE.U(case.param.PE.cs0) .+ phie0[1] # guessed values are not used
            phis_n = zeros(Float64, Np, 1)
            y0 = [csn0;  csp0; ce0]
        else
            error( "Error: $(case.opt.model{1}) model has not been implemented!\n ")
        end
        if case.opt.thermalmodel == "lumped"
            y0 =[y0; case.param.cell.T0]
        elseif case.opt.thermalmodel == "distributed2D"
            # 将二维分布式热的节点温度自由度追加到主状态向量
            nT = case.mesh["thermal2D"].nlen
            T0_nodes = fill(case.param.cell.T0, nT)
            y0 = [y0; T0_nodes]
        end
        if case.opt.model == "P2D"
            y0 =[y0; phis_n; phis_p; phie0]
        end
    else
        y0 = case.opt.y0 
    end
    return y0
end
"""
        ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution=nothing)

为多SPMe并行架构初始化扩展状态向量。
"""
function ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution::Union{Nothing, Vector{Float64}}=nothing)
    # 1) 获取单元数和网格信息
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    
    # 3) 创建单个单元的电化学初始状态（标准SPMe初始化，不含热场）
    # 临时修改 thermalmodel 以获取纯电化学部分
    original_thermalmodel = case.opt.thermalmodel
    case.opt.thermalmodel = "none"
    y0_single_chem = ModelInitialisation(case)
    case.opt.thermalmodel = original_thermalmodel
    
    n_chem = length(y0_single_chem)
    
    # 4) 为每个单元创建独立的电化学初始状态
    y0_chem_all = zeros(Float64, ne * n_chem)
    
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nel = case.mesh["electrolyte"].nlen
    
    for e in 1:ne
        offset = (e - 1) * n_chem
        
        if initial_soc_distribution === nothing
            # 均匀初始化：简单复制
            y0_chem_all[(offset + 1):(offset + n_chem)] .= y0_single_chem
        else
            # 非均匀初始化：根据SOC调整粒子浓度
            soc_e = initial_soc_distribution[e]
            
            # 负极浓度（SOC↑ → cn_surf↑）
            # 假设线性关系：cs = cs_min + SOC × (cs_max - cs_min)
            cn_surf_e = case.param.NE.cs0 * soc_e
            csn0_e = ones(Float64, Nrn) * cn_surf_e
            
            # 正极浓度（SOC↑ → cp_surf↓，因为锂从正极移出）
            # 简化模型：cp_surf = cs_max × (1 - SOC)
            cp_surf_e = case.param.PE.cs0 * (1.0 - soc_e)
            csp0_e = ones(Float64, Nrp) * cp_surf_e
            
            # 电解液浓度（假设均匀）
            ce0_e = ones(Float64, Nel) * case.param.EL.ce0
            
            # 组装单元状态
            yt_e = [csn0_e; csp0_e; ce0_e]
            y0_chem_all[(offset + 1):(offset + n_chem)] .= yt_e
        end
    end
    
    # 5) 添加热场初始温度
    T0_nodes = fill(case.param.cell.T0, nT)
    
    # 6) 组装全局状态向量
    y0 = [y0_chem_all; T0_nodes]
    # 7) 缓存布局信息
    case.layout = MultiSPMeLayout(ne, n_chem, nT)

    return y0
end


"""
    extract_element_state(y, e, layout)

从多SPMe全局状态向量中提取单个单元的电化学状态。
"""
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end


"""
    get_thermal_dofs(y, layout)

从多SPMe全局状态向量中提取热场节点温度。
"""
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end


"""
    update_state(y, layout; element_index, element_state, thermal_nodes)

更新多SPMe全局状态向量（返回新向量）。
"""
function update_state(y::AbstractVector, layout::MultiSPMeLayout;
                      element_index::Union{Nothing,Int}=nothing,
                      element_state::Union{Nothing,Vector{Float64}}=nothing,
                      thermal_nodes::Union{Nothing,Vector{Float64}}=nothing)
    y_new = copy(y)
    if element_index !== nothing
        @assert 1 <= element_index <= layout.ne "element_index $element_index out of range [1, $(layout.ne)]"
        @assert length(element_state) == layout.n_chem "element_state length $(length(element_state)) != n_chem $(layout.n_chem)"
        offset = (element_index - 1) * layout.n_chem
        y_new[(offset + 1):(offset + layout.n_chem)] .= element_state
    end
    if thermal_nodes !== nothing
        @assert length(thermal_nodes) == layout.nT "thermal_nodes length $(length(thermal_nodes)) != nT $(layout.nT)"
        y_new[layout.thermal_range] .= thermal_nodes
    end
    return y_new
end
