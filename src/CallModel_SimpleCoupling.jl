# CallModel_SimpleCoupling.jl
# 简化耦合模型核心求解器：单SPMe + 2D热传导

using SparseArrays: blockdiag
using Statistics: mean

# 导入ThermalDistributed的内部函数（用于热源计算）
# _compute_layer_heat_sources 已在 ThermalDistributed.jl 中完整实现
# 包括：反应热、可逆热、固相/液相欧姆热、集流体欧姆热

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
    # 3. 计算体积热源密度（调用已有函数）
    # ========================================================================
    
    # 直接调用ThermalDistributed的层热源计算函数
    # 注意：这里使用单个SPMe的结果，假设所有单元条件相同
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 计算单个单元的热源密度（W/m³）- 直接调用已有函数
    # 获取总电流（有量纲，A）
    I_app = case.opt.Current(t * case.param.scale.t0)
    Q_layers = _compute_layer_heat_sources(case, variables_chem, T_avg, I_app)
    
    # 获取层权重
    layer_weights = jellyroll_get_layer_weights(mesh)
    
    if layer_weights !== nothing
        # 按层权重计算每个单元的热源密度
        q_volumetric = zeros(Float64, ne)
        for e in 1:ne
            q_volumetric[e] = sum(layer_weights[e, i] * Q_layers[i] for i in 1:5)
        end
    else
        # Fallback：使用简单平均（均匀假设）
        q_uniform = mean(Q_layers)
        q_volumetric = fill(q_uniform, ne)
    end
    
    # 计算总热源功率（用于统计，W）
    V_total = compute_total_volume(case, mesh)
    Q_total = mean(q_volumetric) * V_total
    
    # ========================================================================
    # 4. 求解2D热传导
    # ========================================================================
    
    # 构建热变量字典（ThermalDistributed2D需要的格式）
    variables_thermal = Dict{String, Any}(
        "heat_source_fields" => q_volumetric,  # W/m³
        "heat_source_units_code" => 1.0,        # SI单位
        "average temperature" => T_avg,
        "total heat source" => Q_total
    )
    
    # 添加层权重
    if layer_weights !== nothing
        variables_thermal["layer_weights"] = layer_weights
    end
    
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


# ========================================================================
# 热源计算：直接调用ThermalDistributed的已有函数
# ========================================================================

# 注意：这里导入内部函数（仅用于简化耦合模式）
# _compute_layer_heat_sources 在 ThermalDistributed.jl 中定义
# 它计算5层的体积热源密度 [Q_NE, Q_SP, Q_PE, Q_PCC, Q_NCC]（单位：W/m³）
# 包括：反应热 + 可逆热 + 固相欧姆热 + 液相欧姆热 + 集流体欧姆热


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
