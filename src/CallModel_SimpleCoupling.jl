# CallModel_SimpleCoupling.jl
# 简化耦合模型核心求解器：单SPMe + 2D热传导

using SparseArrays: blockdiag
using Statistics: mean

"""
    CallModel_SimpleCoupling(case::Case, y::Array{Float64}, t::Float64; jacobi::String="update")

简化耦合模型的主求解函数。

# 模型特点
- 电化学：单个SPMe（总电流，全局平均）
- 热传导：2D分布式（空间温度场）
- 耦合：
  - 正向：SPMe热源 → 均匀分配 → 热传导
  - 反向：温度平均 → SPMe参数

# 状态向量布局
y = [y_chem; y_T]
  - y_chem: SPMe化学状态 (n_chem)
  - y_T: 温度场 (nT)

# 返回
- M: 质量矩阵
- K: 刚度矩阵
- F: 载荷向量
- variables: 变量字典
- y: 状态向量（用于phi）
"""
function CallModel_SimpleCoupling(case::Case, y::Array{Float64}, t::Float64; jacobi::String="update")
    # ========================================================================
    # 0. 验证与提取状态
    # ========================================================================
    
    if isempty(case.simple_coupling_layout)
        error("CallModel_SimpleCoupling: simple_coupling_layout未初始化，请先调用ModelInitialisation_SimpleCoupling")
    end
    
    layout = case.simple_coupling_layout
    n_chem = layout["n_chem"]
    nT = layout["nT"]
    
    # 提取化学状态和温度场
    y_vec = vec(y)
    y_chem = y_vec[1:n_chem]
    y_T = y_vec[n_chem+1:end]
    
    # ========================================================================
    # 1. 计算平均温度并反馈给SPMe
    # ========================================================================
    
    T_avg = compute_average_temperature(case, y_T)
    
    # 临时更新参数温度（用于SPMe计算）
    T0_original = case.param.T0
    case.param.T0 = T_avg
    
    # ========================================================================
    # 2. 求解单个SPMe（使用总电流）
    # ========================================================================
    
    if case.opt.model == "SPM"
        M_chem, K_chem, F_chem, variables_chem = SPM(case, y_chem, t; jacobi=jacobi)
    elseif case.opt.model == "SPMe"
        M_chem, K_chem, F_chem, variables_chem = SPMe(case, y_chem, t; jacobi=jacobi)
    else
        error("SimpleCoupling仅支持SPM/SPMe模型，当前模型：$(case.opt.model)")
    end
    
    # 恢复原始温度参数
    case.param.T0 = T0_original
    
    # ========================================================================
    # 3. 计算总内热源
    # ========================================================================
    
    Q_total = compute_total_heat_source_simple(case, variables_chem, y_T)
    
    # ========================================================================
    # 4. 将热源均匀分配到所有单元
    # ========================================================================
    
    q_volumetric = distribute_heat_source_uniform(case, Q_total)
    
    # ========================================================================
    # 5. 求解2D热传导
    # ========================================================================
    
    # 构建热变量字典（包含体积热源）
    variables_thermal = Dict{String, Any}(
        "volumetric heat source" => q_volumetric,
        "average temperature" => T_avg,
        "total heat source" => Q_total
    )
    
    # 调用ThermalDistributed求解
    M_T, K_T, F_T = solve_thermal2D_simple_coupling(case, y_T, variables_thermal, t, jacobi)
    
    # ========================================================================
    # 6. 组装全局系统矩阵
    # ========================================================================
    
    # 块对角组装
    M = blockdiag(M_chem, M_T)
    K = blockdiag(K_chem, K_T)
    F = vcat(F_chem, F_T)
    
    # ========================================================================
    # 7. 组装输出变量
    # ========================================================================
    
    # 合并化学和热变量
    variables = merge(variables_chem, variables_thermal)
    
    # 添加温度场（有量纲）
    variables["thermal2D temperature"] = y_T * case.param.scale.T0
    variables["thermal2D temperature field"] = y_T * case.param.scale.T0
    
    # 调试信息
    if hasproperty(case.opt, :debug_simple_coupling) && case.opt.debug_simple_coupling
        println("[SimpleCoupling @ t=$(t*case.param.scale.t0)s]")
        println("  T_avg = $(T_avg) K")
        println("  Q_total = $(Q_total) W")
        println("  q_vol (avg) = $(mean(q_volumetric)) W/m³")
        println("  V_cell = $(variables_chem["cell voltage"]) V")
    end
    
    return M, K, F, variables, y
end


"""
    compute_total_heat_source_simple(case, variables_chem, y_T)

从SPMe结果计算总内热源（有量纲）。

# 热源组成
- Q_rxn: 反应热（过电位×电流）
- Q_ohm: 欧姆热（I²R）
- Q_rev: 熵热（dU/dT × I × T）

# 返回
- Q_total: 总热源功率 (W)
"""
function compute_total_heat_source_simple(case, variables_chem, y_T)
    # 提取关键变量
    I_app = case.opt.I_app  # 总电流 (A)
    T_avg = mean(y_T) * case.param.scale.T0  # 平均温度 (K)
    
    # 初始化各热源项
    Q_rxn = 0.0
    Q_ohm = 0.0
    Q_rev = 0.0
    
    # ========================================================================
    # 1. 反应热（过电位）
    # ========================================================================
    # 负极反应热
    if haskey(variables_chem, "negative electrode overpotential")
        η_n = variables_chem["negative electrode overpotential"]
        # 如果是数组，取第一个值；如果是标量，直接使用
        η_n_val = isa(η_n, Array) ? η_n[1] : η_n
        
        # 负极反应热 = I × |η_n|
        Q_rxn_n = abs(I_app) * abs(η_n_val)
        Q_rxn += Q_rxn_n
    end
    
    # 正极反应热
    if haskey(variables_chem, "positive electrode overpotential")
        η_p = variables_chem["positive electrode overpotential"]
        η_p_val = isa(η_p, Array) ? η_p[1] : η_p
        
        # 正极反应热 = I × |η_p|
        Q_rxn_p = abs(I_app) * abs(η_p_val)
        Q_rxn += Q_rxn_p
    end
    
    # ========================================================================
    # 2. 欧姆热（电解质电阻）
    # ========================================================================
    # 对于SPMe，电解质欧姆极化可以从电解质电位降计算
    # Q_ohm = I × ΔΦ_e（电解质电位降）
    # 这部分通常已经包含在过电位中，避免重复计算
    # 暂时保留为0，避免过度估计
    Q_ohm = 0.0
    
    # ========================================================================
    # 3. 熵热（可逆热）
    # ========================================================================
    # 熵热 = I × T × (dU_p/dT - dU_n/dT)
    dUdT_n_val = 0.0
    dUdT_p_val = 0.0
    
    if haskey(variables_chem, "negative electrode entropic coefficient")
        dUdT_n = variables_chem["negative electrode entropic coefficient"]
        dUdT_n_val = isa(dUdT_n, Array) ? dUdT_n[1] : dUdT_n
    end
    
    if haskey(variables_chem, "positive electrode entropic coefficient")
        dUdT_p = variables_chem["positive electrode entropic coefficient"]
        dUdT_p_val = isa(dUdT_p, Array) ? dUdT_p[1] : dUdT_p
    end
    
    # 可逆热
    Q_rev = I_app * T_avg * (dUdT_p_val - dUdT_n_val)
    
    # ========================================================================
    # 总热源
    # ========================================================================
    Q_total = Q_rxn + Q_ohm + Q_rev
    
    # 调试输出（可选）
    if hasproperty(case.opt, :debug_simple_coupling) && case.opt.debug_simple_coupling
        if abs(Q_total) > 0.1  # 避免初始零值输出
            @printf("    Heat sources: Q_rxn=%.2f W, Q_ohm=%.2f W, Q_rev=%.2f W\n", 
                    Q_rxn, Q_ohm, Q_rev)
        end
    end
    
    return Q_total
end


"""
    distribute_heat_source_uniform(case, Q_total)

将总热源均匀分配到所有单元。

# 参数
- case: 案例对象
- Q_total: 总热源功率 (W)

# 返回
- q_volumetric: 每个单元的体积热源密度 (W/m³)，向量长度为ne
"""
function distribute_heat_source_uniform(case, Q_total)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 计算总体积
    V_total = compute_total_volume(case, mesh)
    
    # 均匀体积热源密度
    q_uniform = Q_total / V_total  # W/m³
    
    # 返回所有单元相同的热源
    return fill(q_uniform, ne)
end


"""
    compute_total_volume(case, mesh)

计算网格总体积。

# 参数
- case: 案例对象
- mesh: 2D网格

# 返回
- V_total: 总体积 (m³)
"""
function compute_total_volume(case, mesh)
    ne = size(mesh.element, 1)
    ngs = length(mesh.gs.detJ)
    
    # 计算总面积（2D）
    A_total = 0.0
    for g in 1:ngs
        A_total += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    
    # 体积 = 面积 × 高度
    H = case.param.cell.height
    V_total = A_total * H
    
    return V_total
end


"""
    solve_thermal2D_simple_coupling(case, y_T, variables, t, jacobi)

求解2D热传导（简化耦合模式）。

直接调用ThermalDistributed2D，无需重新实现。

# 参数
- case: 案例对象
- y_T: 温度场状态向量（无量纲）
- variables: 变量字典（包含体积热源）
- t: 时间
- jacobi: 雅可比更新策略

# 返回
- M_T: 质量矩阵
- K_T: 刚度矩阵
- F_T: 载荷向量（包含热源和边界条件）
"""
function solve_thermal2D_simple_coupling(case, y_T, variables, t, jacobi)
    # 准备variables字典（ThermalDistributed2D需要的格式）
    # 需要包含"volumetric heat source"和层权重信息
    
    # 添加层权重（如果可用）
    mesh = case.mesh["thermal2D"]
    layer_weights = jellyroll_get_layer_weights(mesh)
    if layer_weights !== nothing
        variables["layer_weights"] = layer_weights
    end
    
    # 调用ThermalDistributed2D
    M_T, K_T, F_T = ThermalDistributed2D(case, variables)
    
    return M_T, K_T, F_T
end
