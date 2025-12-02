# ModelInitialisation_SimpleCoupling.jl
# 简化耦合模型初始化：单SPMe + 2D热传导

"""
    ModelInitialisation_SimpleCoupling(case::Case)

初始化简化耦合模型的状态向量。

# 状态向量布局
y = [
    # 电化学部分（单个SPMe）
    csn_all,    # 负极颗粒浓度 (Nrn)
    csp_all,    # 正极颗粒浓度 (Nrp)
    ce_all,     # 电解质浓度 (Ne) - 仅SPMe模型
    
    # 热传导部分（2D网格）
    T_all       # 温度场 (nT)
]

# 耦合方式
- SPMe使用平均温度
- 热源从SPMe结果均匀分配到所有单元
"""
function ModelInitialisation_SimpleCoupling(case::Case)
    @assert case.opt.model in ["SPM", "SPMe"] "简化耦合模式仅支持SPM/SPMe模型"
    @assert case.opt.thermalmodel == "distributed2D" "简化耦合模式需要2D热传导"
    @assert haskey(case.mesh, "thermal2D") "未找到thermal2D网格"
    
    # ========================================================================
    # 1. 电化学部分初始化（单个SPMe）
    # ========================================================================
    
    # 颗粒网格
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    
    # 初始固相浓度
    csn0 = ones(Float64, Nrn) * case.param.NE.cs0
    csp0 = ones(Float64, Nrp) * case.param.PE.cs0
    
    if case.opt.model == "SPM"
        # SPM: 只有固相浓度
        y_chem = vcat(csn0, csp0)
        n_chem = Nrn + Nrp
        
    elseif case.opt.model == "SPMe"
        # SPMe: 固相浓度 + 电解质浓度
        Ne = case.mesh["electrolyte"].nlen
        ce0 = ones(Float64, Ne) * case.param.EL.ce0
        y_chem = vcat(csn0, csp0, ce0)
        n_chem = Nrn + Nrp + Ne
    end
    
    # ========================================================================
    # 2. 热传导部分初始化（2D温度场）
    # ========================================================================
    
    nT = case.mesh["thermal2D"].nlen
    T_init = case.param.T0 / case.param.scale.T0  # 无量纲化
    y_T = fill(T_init, nT)
    
    # ========================================================================
    # 3. 组合状态向量
    # ========================================================================
    
    y0 = vcat(y_chem, y_T)
    
    # ========================================================================
    # 4. 保存布局信息（用于提取状态）
    # ========================================================================
    
    case.simple_coupling_layout = Dict{String, Any}(
        "n_chem" => n_chem,
        "nT" => nT,
        "total" => length(y0),
        "Nrn" => Nrn,
        "Nrp" => Nrp,
        "model" => case.opt.model
    )
    
    if case.opt.model == "SPMe"
        case.simple_coupling_layout["Ne"] = Ne
    end
    
    # 打印初始化信息
    if hasproperty(case.opt, :debug_simple_coupling) && case.opt.debug_simple_coupling
        println("[SimpleCoupling初始化]")
        println("  模型: $(case.opt.model)")
        println("  化学状态维度: $n_chem")
        println("    - 负极颗粒节点: $Nrn")
        println("    - 正极颗粒节点: $Nrp")
        if case.opt.model == "SPMe"
            println("    - 电解质节点: $Ne")
        end
        println("  温度场维度: $nT")
        println("  总状态向量维度: $(length(y0))")
        println("  初始温度: $(case.param.T0) K")
    end
    
    return y0
end


"""
    extract_states_simple_coupling(case, y)

从状态向量中提取化学状态和温度场。

# 返回
- `y_chem`: 化学状态向量
- `y_T`: 温度场向量
"""
function extract_states_simple_coupling(case, y)
    layout = case.simple_coupling_layout
    n_chem = layout["n_chem"]
    
    y_chem = y[1:n_chem]
    y_T = y[n_chem+1:end]
    
    return y_chem, y_T
end


"""
    compute_average_temperature(case, y_T)

计算温度场的平均值（用于反馈给SPMe）。

# 参数
- `case`: 案例对象
- `y_T`: 无量纲温度场

# 返回
- `T_avg`: 有量纲平均温度 (K)
"""
function compute_average_temperature(case, y_T)
    T_avg_scaled = mean(y_T)
    T_avg = T_avg_scaled * case.param.scale.T0
    return T_avg
end
