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
        elseif case.opt.thermalmodel == "distributed1D"
            # 预留：若未来加入 1D 分布式热网格，可在此追加 DOF
            # 当前无实现
        elseif case.opt.thermalmodel == "distributed2D"
            # 将二维分布式热的节点温度自由度追加到主状态向量
            if haskey(case.mesh, "thermal2D")
                nT = case.mesh["thermal2D"].nlen
                T0_nodes = fill(case.param.cell.T0, nT)
                y0 = [y0; T0_nodes]
            end
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

此函数构造包含所有热单元的独立电化学状态和热场节点温度的全局状态向量：
    y0 = [yt_e[1]; yt_e[2]; ...; yt_e[ne]; T_nodes]

其中每个 yt_e[e] 包含单元 e 的独立电化学状态（粒子和电解液浓度）。

# 参数
- `case::Case`: 案例对象（必须包含 thermal2D 网格）
- `initial_soc_distribution::Union{Nothing, Vector{Float64}}`: 可选的逐单元初始SOC分布
  - 如果为 `nothing`（默认）：所有单元使用相同的初始SOC（从 case.param 读取）
  - 如果为向量：长度必须等于热单元数 ne，每个单元使用对应的SOC初始化
  - SOC范围：[0, 1]（0=完全放电，1=完全充电）

# 返回
- `y0::Vector{Float64}`: 扩展状态向量
  - 结构: [电化学部分(ne个单元); 热场部分(nT个节点)]
  - 长度: ne × n_chem + nT
  - 其中 n_chem = Nrn + Nrp + Nel（单个SPMe单元的电化学自由度）

# 前提条件
1. `case.opt.model` 必须为 "SPMe"
2. `case.opt.thermalmodel` 必须为 "distributed2D"
3. `case.mesh` 必须包含 "thermal2D" 网格
4. `case.opt.per_element_spme` 应为 true（建议，非强制）

# 示例
```julia
# 均匀初始SOC（所有单元相同）
case = SetCase(param_dim, opt)
case.opt.per_element_spme = true
y0 = ModelInitialisation_MultiSPMe(case)

# 非均匀初始SOC（如热失控扩散研究）
ne = size(case.mesh["thermal2D"].element, 1)
soc_dist = range(0.8, 1.0, length=ne)  # 从80%到100%线性分布
y0 = ModelInitialisation_MultiSPMe(case; initial_soc_distribution=soc_dist)
```

# 状态向量布局
假设 ne=3 个热单元，Nrn=10, Nrp=10, Nel=40（总共 n_chem=60），nT=200 个热节点：

```
y0 = [
    # 电化学部分（3×60 = 180个自由度）
    cn_surf_1[1:10]; cp_surf_1[1:10]; ce_1[1:40];   # 单元1
    cn_surf_2[1:10]; cp_surf_2[1:10]; ce_2[1:40];   # 单元2
    cn_surf_3[1:10]; cp_surf_3[1:10]; ce_3[1:40];   # 单元3
    
    # 热场部分（200个自由度）
    T_nodes[1:200];                                   # 所有热节点
]
总长度: 180 + 200 = 380
```

# 注意事项
1. **内存占用**: 状态向量长度为 `ne × n_chem + nT`
2. **初始化策略**: 默认简单复制单个单元的初始化到所有单元
   - 未来可扩展支持空间变化的初始条件（如温度梯度、SOC分布）
3. **与单SPMe兼容**: 单SPMe模式仍使用 `ModelInitialisation`，不受影响
4. **状态提取**: 使用 `MultiSPMe_extract_element_state` 提取单个单元状态

# 相关函数
- `ModelInitialisation`: 单SPMe模式初始化
- `MultiSPMe_extract_element_state`: 从全局向量提取单元状态
- `MultiSPMe_get_thermal_dofs`: 提取热场自由度
"""
function ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution::Union{Nothing, Vector{Float64}}=nothing)
    # 1) 验证前提条件
    if case.opt.model != "SPMe"
        error("ModelInitialisation_MultiSPMe only supports SPMe model, got $(case.opt.model)")
    end
    
    if case.opt.thermalmodel != "distributed2D"
        error("ModelInitialisation_MultiSPMe requires thermalmodel='distributed2D', got $(case.opt.thermalmodel)")
    end
    
    if !haskey(case.mesh, "thermal2D")
        error("ModelInitialisation_MultiSPMe requires thermal2D mesh")
    end
    
    # 2) 获取单元数和网格信息
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    
    # 3) 创建单个单元的电化学初始状态（标准SPMe初始化，不含热场）
    # 临时修改 thermalmodel 以获取纯电化学部分
    original_thermalmodel = case.opt.thermalmodel
    case.opt.thermalmodel = "none"
    y0_single_chem = ModelInitialisation(case)
    case.opt.thermalmodel = original_thermalmodel
    
    n_chem = length(y0_single_chem)
    
    # 4) 处理逐单元初始SOC分布
    if initial_soc_distribution !== nothing
        # 验证SOC分布向量长度
        if length(initial_soc_distribution) != ne
            error("initial_soc_distribution length ($(length(initial_soc_distribution))) must equal number of elements ($ne)")
        end
        
        # 验证SOC范围
        if any(initial_soc_distribution .< 0) || any(initial_soc_distribution .> 1)
            error("initial_soc_distribution values must be in [0, 1], got range [$(minimum(initial_soc_distribution)), $(maximum(initial_soc_distribution))]")
        end
    end
    
    # 5) 为每个单元创建独立的电化学初始状态
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
    
    # 6) 添加热场初始温度
    # DEBUG: 检查 T0 是否有效
    if !hasproperty(case.param.cell, :T0)
        error("❌ case.param.cell.T0 未定义！请在参数初始化时设置初始温度。")
    end
    
    T0_val = case.param.cell.T0
    if !isfinite(T0_val)
        error("❌ case.param.cell.T0 = $T0_val (NaN/Inf)！初始温度必须是有限的正数。")
    end
    
    if T0_val <= 0 || T0_val > 500
        @warn "⚠️  case.param.cell.T0 = $T0_val 看起来不合理（期望 273-373 K 或类似范围）"
    end
    
    T0_nodes = fill(T0_val, nT)
    
    # 7) 组装全局状态向量
    y0 = [y0_chem_all; T0_nodes]
    # DEBUG: 验证初始化结果
    if any(!isfinite, y0)
        nan_in_chem = sum(.!isfinite.(y0_chem_all))
        nan_in_T = sum(.!isfinite.(T0_nodes))
        println("\n" * "="^80)
        println("❌ [DEBUG] 初始状态向量包含 NaN/Inf！")
        println("="^80)
        println("  化学部分 (长度 $(length(y0_chem_all))): $nan_in_chem 个 NaN/Inf")
        println("  温度部分 (长度 $(length(T0_nodes))): $nan_in_T 个 NaN/Inf")
        println("  T0_val = $T0_val")
        if nan_in_chem > 0
            println("\n💡 化学部分有 NaN - 检查 SOC 分布或单元初始化")
        end
        if nan_in_T > 0
            println("\n💡 温度部分有 NaN - case.param.cell.T0 有问题")
        end
        println("="^80 * "\n")
        error("ModelInitialisation_MultiSPMe: 初始状态向量包含 NaN/Inf")
    end
    # 8) 缓存状态向量结构信息（用于后续提取）
    if isempty(case.multi_spme_layout)
        # 初始化或清空后填充布局信息
        empty!(case.multi_spme_layout)
    end
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
    if isempty(case.multi_spme_layout)
        error("case.multi_spme_layout is empty. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    
    if e < 1 || e > ne
        error("Element index $e out of range [1, $ne]")
    end
    
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
    if isempty(case.multi_spme_layout)
        error("case.multi_spme_layout is empty. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.multi_spme_layout
    thermal_range = layout["thermal_range"]
    
    T_nodes = y_vec[thermal_range]
    
    return T_nodes
end


"""
    MultiSPMe_update_element_state!(y::Vector{Float64}, e::Int, yt_e::Vector{Float64}, case::Case)

将单个单元的电化学状态写回全局状态向量（原地修改）。

# 参数
- `y::Vector{Float64}`: 全局状态向量（将被修改）
- `e::Int`: 单元编号（1-based）
- `yt_e::Vector{Float64}`: 单元 e 的新电化学状态
- `case::Case`: 案例对象

# 示例
```julia
# 提取单元状态
yt_e = MultiSPMe_extract_element_state(y_old, e, case)

# 时间推进（求解该单元）
yt_e_new = time_step_solve(yt_e, ...)

# 写回全局向量
MultiSPMe_update_element_state!(y_new, e, yt_e_new, case)
```
"""
function MultiSPMe_update_element_state!(y::Vector{Float64}, e::Int, yt_e::Vector{Float64}, case::Case)
    if isempty(case.multi_spme_layout)
        error("case.multi_spme_layout is empty. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    
    if e < 1 || e > ne
        error("Element index $e out of range [1, $ne]")
    end
    
    if length(yt_e) != n_chem
        error("yt_e length ($(length(yt_e))) must equal n_chem ($n_chem)")
    end
    
    offset = (e - 1) * n_chem
    y[(offset + 1):(offset + n_chem)] .= yt_e
    
    return nothing
end


"""
    MultiSPMe_update_thermal_dofs!(y::Vector{Float64}, T_nodes::Vector{Float64}, case::Case)

将热场节点温度写回全局状态向量（原地修改）。

# 参数
- `y::Vector{Float64}`: 全局状态向量（将被修改）
- `T_nodes::Vector{Float64}`: 新的热场节点温度
- `case::Case`: 案例对象

# 示例
```julia
T_nodes = MultiSPMe_get_thermal_dofs(y_old, case)
T_nodes_new = thermal_solve(T_nodes, ...)
MultiSPMe_update_thermal_dofs!(y_new, T_nodes_new, case)
```
"""
function MultiSPMe_update_thermal_dofs!(y::Vector{Float64}, T_nodes::Vector{Float64}, case::Case)
    if isempty(case.multi_spme_layout)
        error("case.multi_spme_layout is empty. Did you call ModelInitialisation_MultiSPMe?")
    end
    
    layout = case.multi_spme_layout
    nT = layout["nT"]
    thermal_range = layout["thermal_range"]
    
    if length(T_nodes) != nT
        error("T_nodes length ($(length(T_nodes))) must equal nT ($nT)")
    end
    
    y[thermal_range] .= T_nodes
    
    return nothing
end