# CallModel_SimpleCoupling.jl
# 简化耦合模型核心求解器：单SPMe + 2D热传导

using SparseArrays: blockdiag, sparse
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
    isempty(case.simple_coupling_layout) &&
        error("CallModel_SimpleCoupling: simple_coupling_layout未初始化，请先调用ModelInitialisation_SimpleCoupling")

    # 拆分状态向量
    y_chem, y_T = extract_states_simple_coupling(case, y)
    T_avg_nd, T_avg_K = compute_average_temperature(case, y_T)
    case.shared_spme_temperature = T_avg_nd

    # 临时使用平均温度更新 SPMe 参数
    T0_original = case.param.cell.T0
    case.param.cell.T0 = T_avg_nd
    if case.opt.model == "SPM"
        M_chem, K_chem, F_chem, variables_chem = SPM(case, y_chem, t; jacobi=jacobi)
    elseif case.opt.model == "SPMe"
        M_chem, K_chem, F_chem, variables_chem = SPMe(case, y_chem, t; jacobi=jacobi)
    else
        error("SimpleCoupling仅支持SPM/SPMe模型，当前模型：$(case.opt.model)")
    end
    case.param.cell.T0 = T0_original

    # 计算统一热源（假设各热单元条件相同）
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)
    I_total_nd = case.opt.Current(t * case.param.scale.t0) / case.param.scale.I_typ
    Q_layers_ec = try
        _compute_layer_heat_sources(case, variables_chem, T_avg_nd, I_total_nd)
    catch err
        @warn "_compute_layer_heat_sources 失败，热源将置零" err
        zeros(Float64, 5)
    end
    
    # 尺度转换: Q_layers_ec 是无量纲的（相对于电化学尺度 I×φ/V）
    # 需要转换为相对于傅里叶尺度 (k×T/L²) 的无量纲值
    q_ec_scale = case.param_dim.scale.I_typ * case.param_dim.scale.phi / case.param_dim.cell.volume
    q_ref = case.param_dim.scale.q_th  # 傅里叶尺度
    
    # 转换: 电化学无量纲 → 物理值 → 傅里叶无量纲
    Q_layers = Q_layers_ec .* (q_ec_scale / q_ref)

    layer_weights = try
        jellyroll_get_layer_weights(mesh_th)
    catch
        nothing
    end
    if layer_weights === nothing
        layer_weights = jellyroll_element_layer_weights(mesh_th, case.param_dim; nsamples_per_dim=4, logic=:spiral)
    end

    q_elem = zeros(Float64, ne)
    if layer_weights !== nothing && size(layer_weights, 1) == ne && size(layer_weights, 2) >= 5
        @inbounds for e in 1:ne
            q_elem[e] = sum(layer_weights[e, i] * Q_layers[i] for i in 1:5)
        end
    else
        q_mean = mean(Q_layers)
        fill!(q_elem, q_mean)
    end
    
    # q_ref 已在上面定义
    if q_ref <= 0
        @warn "simple coupling heat scaling fallback" q_ref=q_ref
        q_ref = 1.0
    end
    elem_volumes = compute_element_volumes(case, mesh_th)
    total_volume = sum(elem_volumes)
    total_heat_W = q_ref * sum(q_elem .* elem_volumes)
    q_avg_Wm3 = total_volume > 0 ? (total_heat_W / total_volume) : 0.0

    units_SI = hasproperty(case.opt, :units_thermal) && case.opt.units_thermal == "SI"
    if units_SI
        heat_fields = q_elem .* q_ref
        heat_units_code = 1.0
    else
        heat_fields = q_elem
        heat_units_code = 0.0
    end

    thermal_vars = Dict{String, Union{Array{Float64}, Float64}}(
        "T_nodes" => y_T,
        "heat_source_fields" => heat_fields,
        "heat_source_units_code" => heat_units_code,
    )
    thermal_vars["average temperature"] = T_avg_K
    thermal_vars["total heat source"] = total_heat_W
    if layer_weights !== nothing
        thermal_vars["thermal2D layer_weights"] = layer_weights
    end

    M_T, K_T, F_T, thermal_vars = solve_thermal2D_simple_coupling(case, thermal_vars, t)

    # 拼装全局系统
    M = blockdiag(M_chem, sparse(M_T))
    K = blockdiag(K_chem, sparse(K_T))
    F = [F_chem; F_T]

    # 汇总并补充变量
    variables = Dict{String, Union{Array{Float64}, Float64}}()
    merge!(variables, variables_chem)
    merge!(variables, thermal_vars)
    variables["T_nodes"] = y_T
    variables["temperature"] = T_avg_nd
    variables["average temperature"] = T_avg_K
    variables["total heat source"] = total_heat_W
    variables["thermal2D temperature"] = y_T .* case.param_dim.scale.T_ref
    variables["thermal2D temperature field"] = variables["thermal2D temperature"]
    variables["time"] = t

    if hasproperty(case.opt, :debug_simple_coupling) && case.opt.debug_simple_coupling
        println("[SimpleCoupling @ t=$(t * case.param.scale.t0)s]")
        println("  T_avg = $(T_avg_K) K")
        println("  Q_total = $(total_heat_W) W")
        println("  q_avg = $(q_avg_Wm3) W/m^3")
        Vc_val = begin
            Vc_nd = variables_chem["cell voltage"]  # 无量纲电压
            Vc_nd_scalar = isa(Vc_nd, AbstractArray) ? (length(Vc_nd) > 0 ? Vc_nd[1] : 0.0) : Vc_nd
            Vc_nd_scalar * case.param.scale.phi  # 转换为有量纲 [V]
        end
        println("  V_cell = $(Vc_val) V")
    end

    return M, K, F, variables, Float64[]
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
function compute_element_volumes(case, mesh)
    ne = size(mesh.element, 1)
    areas = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    area_total = sum(areas)
    height = hasproperty(case.param_dim.cell, :width) ? case.param_dim.cell.width : case.param_dim.cell.volume / max(area_total, eps())
    return areas .* height
end

function compute_total_volume(case, mesh)
    volumes = compute_element_volumes(case, mesh)
    return sum(volumes)
end


"""
    solve_thermal2D_simple_coupling(case, variables, t)

求解2D热传导（简化耦合模式），并施加时间尺度与边界条件。

# 参数
- case: 案例对象
- variables: 变量字典（包含热源、温度场等）
- t: 无量纲时间

# 返回
- M_T: 质量矩阵
- K_T: 刚度矩阵
- F_T: 载荷向量（包含热源和边界条件）
- vars_local: 更新后的变量字典（含可能的缓存项）
"""
function solve_thermal2D_simple_coupling(case, variables, t)
    vars_local = copy(variables)
    M_T, K_T, F_T = ThermalDistributed2D(case, vars_local)
    t_ratio = case.param_dim.scale.t0 / max(case.param_dim.scale.t_th, eps())
    M_T = M_T .* t_ratio
    ThermalDistributed2D_BC(K_T, F_T, case, t)
    return M_T, K_T, F_T, vars_local
end