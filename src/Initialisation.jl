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
    # 7) 缓存状态向量结构信息（用于后续提取）
    empty!(case.multi_spme_layout)
    case.multi_spme_layout["ne"] = ne
    case.multi_spme_layout["n_chem"] = n_chem
    case.multi_spme_layout["nT"] = nT
    case.multi_spme_layout["n_total"] = length(y0)
    case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
    case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)
    
    return y0
end


"""
    MultiSPMe_extract_element_state(y::Vector{Float64}, e::Int, case::Case) -> Vector{Float64}

从多SPMe全局状态向量中提取单个单元的电化学状态。

# 参数
- `y::Vector{Float64}`: 全局状态向量（由 ModelInitialisation_MultiSPMe 生成）
- `e::Int`: 单元编号（1-based）
- `case::Case`: 案例对象（必须包含 multi_spme_layout 信息）

# 返回
- `yt_e::Vector{Float64}`: 单元 e 的局部电化学状态向量
  - 结构: [cn_surf; cp_surf; ce]
  - 长度: n_chem = Nrn + Nrp + Nel

# 示例
```julia
# 初始化多SPMe状态
y0 = ModelInitialisation_MultiSPMe(case)

# 提取单元5的状态
yt_5 = MultiSPMe_extract_element_state(y0, 5, case)

# 用于求解该单元
M_e, K_e, F_e, vars_e = SPMe_element(case, yt_5, t, 5; I_e=I_e[5], T_e=T_e[5])
```
"""
function MultiSPMe_extract_element_state(y::Array{Float64}, e::Int, case::Case)
    # 自动转换为向量（兼容矩阵输入）
    y_vec = vec(y)
    
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    
    offset = (e - 1) * n_chem
    yt_e = y_vec[(offset + 1):(offset + n_chem)]
    
    return yt_e
end


"""
    MultiSPMe_get_thermal_dofs(y::Vector{Float64}, case::Case) -> Vector{Float64}

从多SPMe全局状态向量中提取热场节点温度。

# 参数
- `y::Vector{Float64}`: 全局状态向量
- `case::Case`: 案例对象（必须包含 multi_spme_layout 信息）

# 返回
- `T_nodes::Vector{Float64}`: 热场节点温度（无量纲）
  - 长度: nT

# 示例
```julia
y0 = ModelInitialisation_MultiSPMe(case)
T_nodes = MultiSPMe_get_thermal_dofs(y0, case)
```
"""
function MultiSPMe_get_thermal_dofs(y::Array{Float64}, case::Case)
    # 自动转换为向量（兼容矩阵输入）
    y_vec = vec(y)
    
    layout = case.multi_spme_layout
    thermal_range = layout["thermal_range"]
    
    T_nodes = y_vec[thermal_range]
    
    return T_nodes
end


"""
    MultiSPMe_update_state(y::Vector{Float64}, case::Case; element_index=nothing, element_state=nothing, thermal_nodes=nothing)

更新多SPMe全局状态向量（返回新向量）。

可选更新：
- `element_index` 与 `element_state`：写回单元电化学状态
- `thermal_nodes`：写回热场节点温度
"""
function MultiSPMe_update_state(y::Vector{Float64},case::Case;element_index::Union{Nothing, Int}=nothing,element_state::Union{Nothing, Vector{Float64}}=nothing,thermal_nodes::Union{Nothing, Vector{Float64}}=nothing)
    if isempty(case.multi_spme_layout)
        error("case.multi_spme_layout is empty. Did you call ModelInitialisation_MultiSPMe?")
    end

    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    nT = layout["nT"]
    thermal_range = layout["thermal_range"]

    y_new = copy(y)

    if element_index !== nothing || element_state !== nothing
        if element_index === nothing || element_state === nothing
            error("element_index and element_state must be provided together")
        end
        if element_index < 1 || element_index > ne
            error("Element index $element_index out of range [1, $ne]")
        end
        if length(element_state) != n_chem
            error("element_state length ($(length(element_state))) must equal n_chem ($n_chem)")
        end
        offset = (element_index - 1) * n_chem
        y_new[(offset + 1):(offset + n_chem)] .= element_state
    end

    if thermal_nodes !== nothing
        if length(thermal_nodes) != nT
            error("thermal_nodes length ($(length(thermal_nodes))) must equal nT ($nT)")
        end
        y_new[thermal_range] .= thermal_nodes
    end

    return y_new
end