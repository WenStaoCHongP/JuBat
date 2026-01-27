"""
    CycleSolver.jl - 充放电循环求解器

实现充放电循环的核心逻辑：
- 充电 → 静置 → 放电 → 静置 → (重复)
- 支持电压/时间截止条件
- 损伤状态跨循环累积
- 温度场循环内继承，循环间可重置

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Printf

# ========================================================================
# 1. 结果结构定义
# ========================================================================

"""
    PhaseResult - 单阶段求解结果

存储一个阶段（充电/静置/放电）的结果。
"""
mutable struct PhaseResult
    phase_type::PhaseType               # 阶段类型
    t_start::Float64                    # 起始时间 (s)
    t_end::Float64                      # 结束时间 (s)
    duration::Float64                   # 实际持续时间 (s)
    
    # 电化学
    V_start::Float64                    # 起始电压 (V)
    V_end::Float64                      # 结束电压 (V)
    capacity::Float64                   # 容量 (Ah)
    terminated_by::Symbol               # 截止原因 :time 或 :voltage
    
    # 热
    T_max::Float64                      # 最高温度 (K)
    T_mean_end::Float64                 # 结束时平均温度 (K)
    
    # 损伤
    D_max::Float64                      # 最大损伤
    D_mean::Float64                     # 平均损伤
    ΔD_max::Float64                     # 损伤增量（相对于阶段开始）
    
    # 最终状态（用于传递到下一阶段）
    final_state::Union{Dict, Nothing}
end

# 默认构造函数
function PhaseResult()
    PhaseResult(
        PHASE_REST, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, :none,
        0.0, 0.0,
        0.0, 0.0, 0.0,
        nothing
    )
end

"""
    CycleResult - 单个循环的结果
"""
mutable struct CycleResult
    cycle_idx::Int                      # 循环编号
    
    # 各阶段结果
    charge::PhaseResult
    rest1::PhaseResult
    discharge::PhaseResult
    rest2::PhaseResult
    
    # 循环汇总
    capacity_charge::Float64            # 充电容量 (Ah)
    capacity_discharge::Float64         # 放电容量 (Ah)
    coulombic_efficiency::Float64       # 库伦效率 (%)
    
    # 损伤汇总
    D_max_end::Float64                  # 循环结束时最大损伤
    D_mean_end::Float64                 # 循环结束时平均损伤
    n_fractured::Int                    # 断裂单元数
    
    # 热汇总
    T_max::Float64                      # 循环内最高温度
end

"""
    CyclingResult - 全部循环的结果汇总
"""
mutable struct CyclingResult
    n_cycles::Int                       # 完成的循环数
    
    # 每循环指标（向量）
    cycle_idx::Vector{Int}
    capacity_charge::Vector{Float64}
    capacity_discharge::Vector{Float64}
    coulombic_efficiency::Vector{Float64}
    D_max::Vector{Float64}
    D_mean::Vector{Float64}
    n_fractured::Vector{Int}
    T_max::Vector{Float64}
    
    # 详细结果（可选）
    cycle_results::Vector{CycleResult}
    
    # 最终状态
    final_czm_mesh::Any                 # 最终的CZM网格（包含损伤状态）
end

function CyclingResult(n_cycles::Int)
    CyclingResult(
        0,
        Int[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Int[], Float64[],
        CycleResult[],
        nothing
    )
end


# ========================================================================
# 2. 单阶段求解器
# ========================================================================

"""
    solve_phase(case, phase_type, t_max, I_current, V_limit, initial_state;
                czm_mesh=nothing, czm_params=nothing)

求解单个阶段（充电/静置/放电）。

# 参数
- `case`: JuBat Case对象
- `phase_type`: 阶段类型 (PHASE_CHARGE, PHASE_REST, PHASE_DISCHARGE)
- `t_max`: 最大持续时间 (s)
- `I_current`: 电流 (A)，充电为负，放电为正，静置为0
- `V_limit`: 电压限制 (V)，充电为上限，放电为下限，静置时无效
- `initial_state`: 初始状态字典

# 返回
- `PhaseResult`: 阶段结果
"""
function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64, 
                     I_current::Float64, V_limit::Float64,
                     initial_state::Dict;
                     czm_mesh=nothing, czm_params=nothing,
                     dt_range::Vector{Float64}=[1.0, 10.0])
    
    result = PhaseResult()
    result.phase_type = phase_type
    result.t_start = get(initial_state, "t_global", 0.0)
    
    # 设置电流
    case.opt.Current = x -> I_current
    
    # 时间设置
    case.opt.time = [0.0, t_max]
    case.opt.dt = dt_range
    
    # 获取初始状态
    y0 = get(initial_state, "y", nothing)
    T_nodes = get(initial_state, "T_nodes", nothing)
    V_init = get(initial_state, "V", nothing)  # 初始电压可能未知
    
    # V_start 将在实际求解后更新
    result.V_start = V_init !== nothing ? V_init : NaN
    
    # 根据模式选择求解
    # 简化：调用内部求解循环
    phase_result_data = _solve_phase_internal(
        case, phase_type, t_max, I_current, V_limit,
        y0, T_nodes, czm_mesh, czm_params, dt_range)
    
    # 填充结果
    result.t_end = result.t_start + phase_result_data["duration"]
    result.duration = phase_result_data["duration"]
    result.V_start = phase_result_data["V_start"]  # 使用实际计算的初始电压
    result.V_end = phase_result_data["V_end"]
    result.capacity = phase_result_data["capacity"]
    result.terminated_by = phase_result_data["terminated_by"]
    result.T_max = phase_result_data["T_max"]
    result.T_mean_end = phase_result_data["T_mean_end"]
    result.D_max = phase_result_data["D_max"]
    result.D_mean = phase_result_data["D_mean"]
    result.ΔD_max = phase_result_data["ΔD_max"]
    result.final_state = phase_result_data["final_state"]
    
    # 静置阶段：添加锂扩散信息到 final_state
    if phase_type == PHASE_REST && haskey(phase_result_data, "diffusion_active")
        result.final_state["diffusion_active"] = phase_result_data["diffusion_active"]
        result.final_state["cs_relaxation_n"] = get(phase_result_data, "cs_relaxation_n", 0.0)
        result.final_state["cs_relaxation_p"] = get(phase_result_data, "cs_relaxation_p", 0.0)
    end
    
    return result
end

"""
    _solve_phase_internal(case, phase_type, t_max, I_current, V_limit, 
                          y0, T_nodes, czm_mesh, czm_params, dt_range)

内部阶段求解实现。

对于静置阶段（PHASE_REST）：
- 电化学状态完全继承上一步的最终状态
- 电流设为 0，但继续进行锂扩散计算
- 颗粒内锂浓度会因扩散而趋向均匀化
"""
function _solve_phase_internal(case::Case, phase_type::PhaseType, 
                               t_max::Float64, I_current::Float64, V_limit::Float64,
                               y0, T_nodes, czm_mesh, czm_params, dt_range)
    
    # 静置阶段特殊处理：确保电流为0
    if phase_type == PHASE_REST
        I_current = 0.0  # 强制设为0
    end
    
    # 时间缩放
    t0_scale = case.param.scale.t0
    dt_min = dt_range[1] / t0_scale
    dt_max = dt_range[2] / t0_scale
    t_end_nd = t_max / t0_scale
    
    # 检测是否为多SPMe模式（当使用分布式2D热模型时自动启用）
    multi_spme = case.opt.model == "SPMe" && 
                 case.opt.thermalmodel == "distributed2D" &&
                 haskey(case.mesh, "thermal2D")
    
    # 初始化或验证状态向量
    # 静置阶段特别说明：y0 必须从上一阶段继承，不能重新初始化
    is_rest_phase = (phase_type == PHASE_REST)
    
    if y0 === nothing
        if is_rest_phase
            @warn "[CycleSolver] 静置阶段应继承上一步状态，但y0为空，将使用初始状态"
        end
        # 使用标准初始化
        if multi_spme
            y0 = ModelInitialisation_MultiSPMe(case)
        else
            y0 = ModelInitialisation(case)
        end
    else
        # 状态向量从上一阶段传递（静置阶段的正常行为）
        y0 = vec(y0)  # 确保是向量格式
        
        if multi_spme
            # 确保多SPMe布局已初始化
            if isempty(case.multi_spme_layout)
                _ensure_multi_spme_layout!(case)
            end
            
            # 验证状态向量长度
            expected_len = case.multi_spme_layout["n_total"]
            if length(y0) != expected_len
                println("⚠️  [CycleSolver] 状态向量长度不匹配: got=$(length(y0)), expected=$expected_len, phase=$phase_type")
                println("⚠️  [CycleSolver] 状态将被重新初始化为初始SOC！")
                y0 = ModelInitialisation_MultiSPMe(case)
            end
        else
            # 非多SPMe模式：检查是否需要追加热场自由度
            if case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
                nT = case.mesh["thermal2D"].nlen
                # 基础电化学自由度（不含热场）
                n_chem_base = case.mesh["negative particle"].nlen + 
                              case.mesh["positive particle"].nlen + 
                              case.mesh["electrolyte"].nlen
                
                if length(y0) == n_chem_base
                    # 状态向量缺少热场自由度，需要追加
                    T_init = T_nodes !== nothing ? T_nodes : fill(case.param.cell.T0, nT)
                    y0 = [y0; T_init]
                end
            end
        end
    end
    
    # 初始化温度场
    # T_nodes 有三种情况：
    # 1. 有效数组（非空）：使用传入的温度场
    # 2. nothing 或空数组，且 y0 包含热场：从 y0 中提取（继承）
    # 3. nothing 或空数组，且 y0 不包含热场：使用初始温度
    T_nodes_valid = T_nodes !== nothing && !isempty(T_nodes)
    
    if T_nodes_valid && haskey(case.mesh, "thermal2D")
        # 情况1：使用传入的有效温度场
        T_nodes_carry = copy(T_nodes)
        
        # 多SPMe模式：同步温度场到状态向量中
        if multi_spme && !isempty(case.multi_spme_layout)
            nT = case.multi_spme_layout["nT"]
            if length(T_nodes_carry) == nT
                thermal_range = case.multi_spme_layout["thermal_range"]
                y0[thermal_range] .= T_nodes_carry
            end
        end
    elseif multi_spme && !isempty(case.multi_spme_layout) && haskey(case.mesh, "thermal2D")
        # 情况2：多SPMe模式，T_nodes为nothing，从状态向量y0中提取温度场（继承）
        thermal_range = case.multi_spme_layout["thermal_range"]
        T_nodes_carry = y0[thermal_range]
    elseif haskey(case.mesh, "thermal2D")
        nT_mesh = case.mesh["thermal2D"].nlen
        
        # 检查是否应该重置温度场到初始温度
        # 条件：T_nodes 被显式设为 nothing（通常表示用户请求重置）
        # 且状态向量 y0 已经包含有效的温度场数据
        should_reset_thermal = (T_nodes === nothing)
        
        if multi_spme && !isempty(case.multi_spme_layout)
            # 多SPMe模式：检查状态向量是否包含温度场
            thermal_range = case.multi_spme_layout["thermal_range"]
            has_valid_thermal = (length(y0) >= thermal_range[end])
            
            if has_valid_thermal && should_reset_thermal
                # 重置温度场到初始温度
                T_nodes_carry = fill(case.param.cell.T0, nT_mesh)
                y0[thermal_range] .= T_nodes_carry
            elseif has_valid_thermal
                # 从状态向量提取（继承）
                T_nodes_carry = y0[thermal_range]
            else
                # 状态向量不包含温度场，使用初始温度
                T_nodes_carry = fill(case.param.cell.T0, nT_mesh)
            end
        else
            # 非多SPMe模式：尝试从 y0 末尾提取温度（如果状态向量包含热场）
            n_chem_base = case.mesh["negative particle"].nlen + 
                          case.mesh["positive particle"].nlen + 
                          case.mesh["electrolyte"].nlen
            if length(y0) > n_chem_base
                # y0 包含热场，从末尾提取
                T_nodes_carry = y0[(end - nT_mesh + 1):end]
            else
                # 使用初始温度
                T_nodes_carry = fill(case.param.cell.T0, nT_mesh)
            end
        end
    else
        T_nodes_carry = Float64[]
    end
    
    # Crank-Nicolson 参数
    theta = case.opt.solveType == "Crank-Nicolson" ? 0.5 : 
            (case.opt.solveType == "forward" ? 0.0 : 1.0)
    
    # 初始调用
    t = 0.0
    dt = dt_min
    
    # 调试：验证 CallModel 前的布局状态
    if multi_spme
        layout_empty = isempty(case.multi_spme_layout)
        if layout_empty
            println("⚠️  [CycleSolver] CallModel 前布局为空，尝试初始化, phase=$phase_type")
            _ensure_multi_spme_layout!(case)
        end
    end
    
    M_old, K_old, F_old, variables, y_phi = CallModel(case, y0, t, jacobi="update")
    
    # 验证 M 矩阵大小与 y0 一致（关键检查！）
    M_size = size(M_old, 1)
    y0_len = length(y0)
    if M_size != y0_len
        println("⚠️  [CycleSolver] M矩阵大小与y0不匹配: M_size=$M_size, y0_len=$y0_len, phase=$phase_type")
        println("⚠️  [CycleSolver] multi_spme_layout: empty=$(isempty(case.multi_spme_layout))")
        
        # 尝试修复：如果布局为空，重新初始化
        if multi_spme && isempty(case.multi_spme_layout)
            println("⚠️  [CycleSolver] 尝试重新初始化布局...")
            _ensure_multi_spme_layout!(case)
        end
        
        # 如果仍然不匹配，报错
        if M_size != y0_len
            error("无法解决 M 矩阵与 y0 的维度不匹配问题。请检查 multi_spme_layout 是否正确设置。" *
                  " M_size=$M_size, y0_len=$y0_len, layout_empty=$(isempty(case.multi_spme_layout))")
        end
    end
    
    # 初始化变量
    V_current = variables["cell voltage"] * case.param.scale.phi
    V_start_actual = V_current  # 保存实际的初始电压
    capacity = 0.0
    T_max_phase = haskey(variables, "T_nodes") ? 
                  maximum(variables["T_nodes"]) * case.param_dim.scale.T_ref : 
                  case.param_dim.cell.T0
    
    # 损伤初始状态
    D_max_init = czm_mesh !== nothing ? maximum(s.D for s in czm_mesh.damage_states) : 0.0
    D_mean_init = czm_mesh !== nothing ? mean(s.D for s in czm_mesh.damage_states) : 0.0
    
    # 静置阶段：记录初始锂浓度分布（用于分析扩散过程）
    cs_n_init_var = 0.0  # 负极锂浓度方差
    cs_p_init_var = 0.0  # 正极锂浓度方差
    if is_rest_phase
        # 提取初始锂浓度分布
        Nrn = case.mesh["negative particle"].nlen
        Nrp = case.mesh["positive particle"].nlen
        
        if multi_spme && !isempty(case.multi_spme_layout)
            # 多SPMe模式：提取所有单元的平均锂浓度
            ne = case.multi_spme_layout["ne"]
            n_chem = case.multi_spme_layout["n_chem"]
            
            cs_n_all = Float64[]
            cs_p_all = Float64[]
            for e in 1:ne
                offset = (e - 1) * n_chem
                cs_n_e = y0[(offset + 1):(offset + Nrn)]
                cs_p_e = y0[(offset + Nrn + 1):(offset + Nrn + Nrp)]
                push!(cs_n_all, mean(cs_n_e))
                push!(cs_p_all, mean(cs_p_e))
            end
            cs_n_init_var = var(cs_n_all)
            cs_p_init_var = var(cs_p_all)
        else
            # 单SPMe模式
            cs_n = y0[1:Nrn]
            cs_p = y0[(Nrn+1):(Nrn+Nrp)]
            cs_n_init_var = var(cs_n)
            cs_p_init_var = var(cs_p)
        end
    end
    
    # 初始求解步
    vc = 1:size(M_old, 1)
    dt_init = 1e-8
    y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
    y_old = vcat(y_c, y_phi)
    
    terminated_by = :time
    t_actual = 0.0
    
    # CZM 位移场初始化
    u_czm_prev = nothing
    if czm_mesh !== nothing
        ndof_czm = 2 * czm_mesh.nnode
        u_czm_prev = zeros(Float64, ndof_czm)
    end
    
    # CZM 更新计数器（不需要每一步都更新，每 N 步更新一次即可）
    czm_update_interval = 10  # 每 10 个时间步更新一次 CZM
    step_count = 0
    
    # 主循环
    while t < t_end_nd
        step_count += 1
        
        # 更新温度影响
        if case.opt.thermal_enabled && !isempty(T_nodes_carry)
            case.param.cell.T0 = mean(T_nodes_carry)
        end
        
        # 电化学步
        M_new, K_new, F_new, variables, y_phi = CallModel(case, y_old, t, jacobi="update")
        Mt = M_new - theta * K_new * dt
        Kt = (1 - theta) * K_old * dt + M_new
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt
        y_c = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old[vc] + Ft)
        y_new = vcat(y_c, y_phi)
        
        # 提取温度场
        if case.opt.thermal_enabled && haskey(case.mesh, "thermal2D")
            nT = case.mesh["thermal2D"].nlen
            if length(y_c) >= nT
                T_nodes_carry = y_c[(end - nT + 1):end]
                T_max_current = maximum(T_nodes_carry) * case.param_dim.scale.T_ref
                T_max_phase = max(T_max_phase, T_max_current)
            end
        end
        
        # ============================================================
        # CZM 损伤计算（周期性更新）
        # ============================================================
        if czm_mesh !== nothing && czm_params !== nothing && (step_count % czm_update_interval == 0)
            # 每100步输出一次调试信息
            debug_czm = (step_count % (czm_update_interval * 10) == 0)
            
            try
                u_czm_prev, czm_converged = _update_czm_damage!(
                    czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev;
                    debug=debug_czm
                )
                if !czm_converged
                    @debug "CZM solver did not converge at t=$(t * t0_scale)s"
                end
            catch e
                @debug "CZM update failed at t=$(t * t0_scale)s: $e"
            end
        end
        # ============================================================
        
        # 更新电压和容量
        V_current = variables["cell voltage"] * case.param.scale.phi
        dt_dim = dt * t0_scale
        capacity += abs(I_current) * dt_dim / 3600.0  # Ah
        
        # 检查截止条件
        if phase_type == PHASE_CHARGE && V_current >= V_limit
            terminated_by = :voltage
            break
        elseif phase_type == PHASE_DISCHARGE && V_current <= V_limit
            terminated_by = :voltage
            break
        end
        
        # 更新状态
        y_old = y_new
        K_old = K_new
        F_old = F_new
        t += dt
        t_actual = t * t0_scale
        
        # 自适应时间步长
        if case.opt.dtType == "auto"
            error_y = ErrorEstimation(case, y_old, y_new, dt_min/dt)
            if error_y < 0.5 * case.opt.dtThreshold
                dt = min(dt * 2, dt_max)
            elseif error_y > case.opt.dtThreshold
                dt = max(dt / 2, dt_min)
            end
        end
        
        # 确保不超过最大时间
        if t + dt > t_end_nd
            dt = t_end_nd - t
        end
    end
    
    # ============================================================
    # 最终 CZM 损伤计算（确保最后一步也计算损伤）
    # ============================================================
    D_max_end = D_max_init
    D_mean_end = D_mean_init
    if czm_mesh !== nothing && czm_params !== nothing && t_actual > 0
        # 最后一次 CZM 更新（如果循环中没有刚更新）
        if step_count % czm_update_interval != 0
            try
                u_czm_prev, _ = _update_czm_damage!(
                    czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev
                )
            catch e
                @debug "Final CZM update failed: $e"
            end
        end
        
        # 获取最终损伤统计
        stats = get_damage_statistics(czm_mesh)
        D_max_end = stats.max_D
        D_mean_end = stats.mean_D
    end
    
    # 静置阶段：计算最终锂浓度分布并分析扩散效果
    cs_relaxation_n = 0.0  # 负极锂浓度松弛率（方差减少百分比）
    cs_relaxation_p = 0.0  # 正极锂浓度松弛率
    if is_rest_phase && t_actual > 0
        # 提取最终锂浓度分布
        Nrn = case.mesh["negative particle"].nlen
        Nrp = case.mesh["positive particle"].nlen
        
        cs_n_final_var = 0.0
        cs_p_final_var = 0.0
        
        if multi_spme && !isempty(case.multi_spme_layout)
            # 多SPMe模式
            ne = case.multi_spme_layout["ne"]
            n_chem = case.multi_spme_layout["n_chem"]
            
            cs_n_all = Float64[]
            cs_p_all = Float64[]
            for e in 1:ne
                offset = (e - 1) * n_chem
                cs_n_e = y_old[(offset + 1):(offset + Nrn)]
                cs_p_e = y_old[(offset + Nrn + 1):(offset + Nrn + Nrp)]
                push!(cs_n_all, mean(cs_n_e))
                push!(cs_p_all, mean(cs_p_e))
            end
            cs_n_final_var = var(cs_n_all)
            cs_p_final_var = var(cs_p_all)
        else
            # 单SPMe模式
            cs_n = y_old[1:Nrn]
            cs_p = y_old[(Nrn+1):(Nrn+Nrp)]
            cs_n_final_var = var(cs_n)
            cs_p_final_var = var(cs_p)
        end
        
        # 计算松弛率（方差减少百分比）
        if cs_n_init_var > 1e-12
            cs_relaxation_n = 100.0 * (1.0 - cs_n_final_var / cs_n_init_var)
        end
        if cs_p_init_var > 1e-12
            cs_relaxation_p = 100.0 * (1.0 - cs_p_final_var / cs_p_init_var)
        end
    end
    
    # 构建最终状态
    # 确保保存的是最新的电压值（从 variables 获取，如果有的话）
    V_final = V_current
    if haskey(variables, "cell voltage")
        V_final = variables["cell voltage"] * case.param.scale.phi
    end
    
    final_state = Dict{String, Any}(
        "y" => copy(y_old),  # 使用 copy 避免引用问题
        "T_nodes" => copy(T_nodes_carry),
        "V" => V_final,
        "t_global" => t_actual
    )
    
    # 构建返回结果
    result_dict = Dict(
        "duration" => t_actual,
        "V_start" => V_start_actual,  # 实际计算的初始电压
        "V_end" => V_current,
        "capacity" => capacity,
        "terminated_by" => terminated_by,
        "T_max" => T_max_phase,
        "T_mean_end" => !isempty(T_nodes_carry) ? 
                        mean(T_nodes_carry) * case.param_dim.scale.T_ref : 
                        case.param_dim.cell.T0,
        "D_max" => D_max_end,
        "D_mean" => D_mean_end,
        "ΔD_max" => D_max_end - D_max_init,
        "final_state" => final_state
    )
    
    # 静置阶段额外信息
    if is_rest_phase
        result_dict["cs_relaxation_n"] = cs_relaxation_n
        result_dict["cs_relaxation_p"] = cs_relaxation_p
        result_dict["diffusion_active"] = true  # 标记扩散过程已执行
    end
    
    return result_dict
end


# ========================================================================
# 3. 循环求解器
# ========================================================================

"""
    solve_cycling(case, cycle_opt, czm_mesh; verbose=true)

执行完整的充放电循环仿真。

# 参数
- `case`: JuBat Case对象
- `cycle_opt`: CycleOption 循环参数
- `czm_mesh`: CohesiveMesh 内聚力网格（可选）

# 返回
- `CyclingResult`: 循环结果汇总
"""
function solve_cycling(case::Case, cycle_opt::CycleOption, czm_mesh=nothing;
                       verbose::Bool=true, save_detailed::Bool=false)
    
    n_cycles = cycle_opt.n_cycles
    result = CyclingResult(n_cycles)
    
    # 应用初始SOC设置
    soc_init = cycle_opt.SOC_init
    if soc_init >= 0.0 && soc_init <= 1.0
        # 获取param_dim（从case.param_dim，如果存在）
        if hasproperty(case, :param_dim) && case.param_dim !== nothing
            cs0_NE, cs0_PE = apply_initial_soc!(case, case.param_dim, soc_init)
            if verbose
                @printf("  初始SOC: %.1f%%\n", soc_init * 100)
                @printf("    → 负极cs0: %.1f mol/m³\n", cs0_NE)
                @printf("    → 正极cs0: %.1f mol/m³\n", cs0_PE)
            end
        else
            # 回退：直接计算归一化参数（不修改param_dim）
            theta_n = case.param.NE.theta_0 + soc_init * (case.param.NE.theta_100 - case.param.NE.theta_0)
            theta_p = case.param.PE.theta_0 - soc_init * (case.param.PE.theta_0 - case.param.PE.theta_100)
            case.param.NE.cs0 = theta_n
            case.param.PE.cs0 = theta_p
            if verbose
                @printf("  初始SOC: %.1f%% (归一化模式)\n", soc_init * 100)
                @printf("    → 负极θ_n: %.4f\n", theta_n)
                @printf("    → 正极θ_p: %.4f\n", theta_p)
            end
        end
    end
    
    # 初始化
    if verbose
        println("="^60)
        println("开始充放电循环仿真")
        println("="^60)
        @printf("  循环次数: %d\n", n_cycles)
        @printf("  初始SOC: %.1f%%\n", soc_init * 100)
        println("  循环顺序: 放电 → 静置 → 充电 → 静置")
        @printf("  放电: %.0fs (%.1fC), 截止 %.2fV\n", 
                cycle_opt.t_discharge, cycle_opt.I_discharge/5.0, cycle_opt.V_lower)
        @printf("  充电: %.0fs (%.1fC), 截止 %.2fV\n", 
                cycle_opt.t_charge, cycle_opt.I_charge/5.0, cycle_opt.V_upper)
        @printf("  静置: %.0fs + %.0fs\n", cycle_opt.t_rest1, cycle_opt.t_rest2)
    end
    
    # 初始状态
    current_state = Dict{String, Any}(
        "y" => nothing,
        "T_nodes" => nothing,
        "V" => 3.7,
        "t_global" => 0.0
    )
    
    # CZM参数
    czm_params = hasproperty(case.param_dim, :cohesive) ? case.param_dim.cohesive : nothing
    
    # 循环主体
    for cycle in 1:n_cycles
        if verbose
            println("\n" * "-"^40)
            @printf("循环 %d/%d\n", cycle, n_cycles)
            println("-"^40)
        end
        
        cycle_result = CycleResult(
            cycle,
            PhaseResult(), PhaseResult(), PhaseResult(), PhaseResult(),
            0.0, 0.0, 0.0,
            0.0, 0.0, 0,
            0.0
        )
        
        # 重置温度（如果需要）
        if cycle_opt.reset_T_each_cycle && cycle > 1
            current_state["T_nodes"] = nothing
            if verbose
                println("  温度场已重置")
            end
        end
        
        # 损伤不重置，跨循环累积
        D_max_cycle_start = czm_mesh !== nothing ? 
                            maximum(s.D for s in czm_mesh.damage_states) : 0.0
        
        # ============ 阶段1: 放电 ============
        # 循环顺序：放电 → 静置 → 充电 → 静置
        # （适配Jellyroll参数的初始高SOC状态）
        if verbose
            print("  [放电] ")
        end
        
        discharge_result = solve_phase(
            case, PHASE_DISCHARGE,
            cycle_opt.t_discharge,
            cycle_opt.I_discharge,  # 正电流表示放电
            cycle_opt.V_lower,
            current_state;
            czm_mesh=czm_mesh,
            czm_params=czm_params,
            dt_range=cycle_opt.dt_cycle
        )
        cycle_result.discharge = discharge_result
        current_state = discharge_result.final_state
        
        if verbose
            @printf("%.1fs, %.3fV→%.3fV, %.3fAh (%s)\n",
                    discharge_result.duration, discharge_result.V_start,
                    discharge_result.V_end, discharge_result.capacity,
                    discharge_result.terminated_by)
            # 调试：打印放电阶段后的状态信息
            y_out = get(current_state, "y", nothing)
            V_out = get(current_state, "V", NaN)
            @printf("    → V_out=%.3fV, y_len=%d\n", V_out, y_out === nothing ? 0 : length(y_out))
        end
        
        # ============ 阶段2: 静置1 ============
        # 静置阶段：电化学状态继承上一步，以0电流继续锂扩散
        if cycle_opt.t_rest1 > 0
            # 静置时间 > 0：执行静置阶段，继续锂扩散过程
            if verbose
                print("  [静置1] 状态继承+锂扩散 ")
            end
            
            rest1_result = solve_phase(
                case, PHASE_REST,
                cycle_opt.t_rest1,
                0.0,  # 零电流：仅进行锂扩散
                0.0,  # 无电压限制
                current_state;  # 继承上一步状态
                czm_mesh=czm_mesh,
                czm_params=czm_params,
                dt_range=cycle_opt.dt_cycle
            )
            cycle_result.rest1 = rest1_result
            current_state = rest1_result.final_state
            
            if verbose
                @printf("%.1fs, T_max=%.2fK\n", rest1_result.duration, rest1_result.T_max)
                y_out = get(current_state, "y", nothing)
                V_out = get(current_state, "V", NaN)
                @printf("    → V_out=%.3fV, y_len=%d", V_out, y_out === nothing ? 0 : length(y_out))
                # 显示锂浓度松弛信息（如果有）
                if haskey(rest1_result.final_state, "cs_relaxation_n") || 
                   hasproperty(rest1_result, :cs_relaxation_n)
                    # 从 solve_phase 结果中提取松弛信息
                    # 注意：这些信息可能在内部结果字典中
                end
                println()
            end
        else
            # 静置时间 = 0：跳过静置阶段，直接继承上一阶段状态
            cycle_result.rest1 = PhaseResult()
            cycle_result.rest1.duration = 0.0
            cycle_result.rest1.V_start = get(current_state, "V", NaN)
            cycle_result.rest1.V_end = cycle_result.rest1.V_start
            cycle_result.rest1.final_state = current_state
            # current_state 保持不变（无扩散过程）
            
            if verbose
                println("  [静置1] 跳过 (t=0，无扩散)")
            end
        end
        
        # ============ 阶段3: 充电 ============
        # 充电前温度场重置（可选）
        if hasproperty(cycle_opt, :reset_T_before_charge) && cycle_opt.reset_T_before_charge
            # 重置温度场到初始温度，但保留电化学状态
            current_state["T_nodes"] = nothing
            if verbose
                println("    (温度场已重置)")
            end
        end
        
        if verbose
            print("  [充电] ")
            # 调试：打印传入充电阶段的状态
            y_in = get(current_state, "y", nothing)
            V_in = get(current_state, "V", NaN)
            @printf("(V_in=%.3fV, y_len=%d) ", V_in, y_in === nothing ? 0 : length(y_in))
        end
        
        charge_result = solve_phase(
            case, PHASE_CHARGE, 
            cycle_opt.t_charge,
            -cycle_opt.I_charge,  # 负电流表示充电
            cycle_opt.V_upper,
            current_state;
            czm_mesh=czm_mesh, 
            czm_params=czm_params,
            dt_range=cycle_opt.dt_cycle
        )
        cycle_result.charge = charge_result
        current_state = charge_result.final_state
        
        if verbose
            @printf("%.1fs, %.3fV→%.3fV, %.3fAh (%s)\n",
                    charge_result.duration, charge_result.V_start, 
                    charge_result.V_end, charge_result.capacity,
                    charge_result.terminated_by)
        end
        
        # ============ 阶段4: 静置2 ============
        # 静置阶段：电化学状态继承上一步，以0电流继续锂扩散
        if cycle_opt.t_rest2 > 0
            # 静置时间 > 0：执行静置阶段，继续锂扩散过程
            if verbose
                print("  [静置2] 状态继承+锂扩散 ")
            end
            
            rest2_result = solve_phase(
                case, PHASE_REST,
                cycle_opt.t_rest2,
                0.0,  # 零电流：仅进行锂扩散
                0.0,
                current_state;  # 继承上一步状态
                czm_mesh=czm_mesh,
                czm_params=czm_params,
                dt_range=cycle_opt.dt_cycle
            )
            cycle_result.rest2 = rest2_result
            current_state = rest2_result.final_state
            
            if verbose
                @printf("%.1fs, T_max=%.2fK\n", rest2_result.duration, rest2_result.T_max)
            end
        else
            # 静置时间 = 0：跳过静置阶段，无扩散过程
            cycle_result.rest2 = PhaseResult()
            cycle_result.rest2.duration = 0.0
            cycle_result.rest2.V_start = get(current_state, "V", NaN)
            cycle_result.rest2.V_end = cycle_result.rest2.V_start
            cycle_result.rest2.final_state = current_state
            # current_state 保持不变（无扩散过程）
            
            if verbose
                println("  [静置2] 跳过 (t=0，无扩散)")
            end
        end
        
        # ============ 循环汇总 ============
        cycle_result.capacity_charge = charge_result.capacity
        cycle_result.capacity_discharge = discharge_result.capacity
        cycle_result.coulombic_efficiency = cycle_result.capacity_charge > 0 ?
            100.0 * cycle_result.capacity_discharge / cycle_result.capacity_charge : 0.0
        
        # 损伤汇总
        if czm_mesh !== nothing
            stats = get_damage_statistics(czm_mesh)
            cycle_result.D_max_end = stats.max_D
            cycle_result.D_mean_end = stats.mean_D
            cycle_result.n_fractured = stats.n_fractured
        end
        
        cycle_result.T_max = max(charge_result.T_max, rest1_result.T_max,
                                 discharge_result.T_max, rest2_result.T_max)
        
        # 记录到结果
        push!(result.cycle_idx, cycle)
        push!(result.capacity_charge, cycle_result.capacity_charge)
        push!(result.capacity_discharge, cycle_result.capacity_discharge)
        push!(result.coulombic_efficiency, cycle_result.coulombic_efficiency)
        push!(result.D_max, cycle_result.D_max_end)
        push!(result.D_mean, cycle_result.D_mean_end)
        push!(result.n_fractured, cycle_result.n_fractured)
        push!(result.T_max, cycle_result.T_max)
        
        if save_detailed
            push!(result.cycle_results, cycle_result)
        end
        
        result.n_cycles = cycle
        
        # 进度输出
        if verbose
            @printf("  → 循环%d完成: 充%.3fAh, 放%.3fAh, CE=%.1f%%, D_max=%.2f%%\n",
                    cycle, cycle_result.capacity_charge, cycle_result.capacity_discharge,
                    cycle_result.coulombic_efficiency, cycle_result.D_max_end * 100)
        end
        
        # 检查是否有单元完全断裂
        if czm_mesh !== nothing && cycle_result.n_fractured > 0.5 * czm_mesh.n_cohesive
            if verbose
                @warn "超过50%的内聚力单元断裂，提前终止循环"
            end
            break
        end
    end
    
    result.final_czm_mesh = czm_mesh
    
    if verbose
        println("\n" * "="^60)
        println("循环仿真完成")
        println("="^60)
        @printf("  完成循环数: %d\n", result.n_cycles)
        if result.n_cycles > 0
            @printf("  最终容量: 充%.3fAh, 放%.3fAh\n", 
                    result.capacity_charge[end], result.capacity_discharge[end])
            @printf("  容量保持率: %.1f%%\n", 
                    100 * result.capacity_discharge[end] / result.capacity_discharge[1])
            @printf("  最终损伤: D_max=%.2f%%, D_mean=%.2f%%\n",
                    result.D_max[end] * 100, result.D_mean[end] * 100)
        end
    end
    
    return result
end


# ========================================================================
# 4. 辅助函数
# ========================================================================

"""
    compute_cs0_from_soc(param_dim, soc::Float64) -> (cs0_NE, cs0_PE)

根据目标SOC计算正负极的初始锂浓度。

# 参数
- `param_dim`: 包含电极参数的维度参数结构体
- `soc::Float64`: 目标SOC (0~1)

# 返回
- `cs0_NE::Float64`: 负极初始锂浓度 (mol/m³)
- `cs0_PE::Float64`: 正极初始锂浓度 (mol/m³)

# 计算公式
负极：θ_n = θ_0_n + SOC × (θ_100_n - θ_0_n)
      cs_n = θ_n × cs_max_n

正极：θ_p = θ_0_p - SOC × (θ_0_p - θ_100_p)
      cs_p = θ_p × cs_max_p

# 示例
```julia
param_dim = JuBat.ChooseCell("Jellyroll")
cs0_NE, cs0_PE = compute_cs0_from_soc(param_dim, 0.8)  # 80% SOC
param_dim.NE.cs0 = cs0_NE
param_dim.PE.cs0 = cs0_PE
```
"""
function compute_cs0_from_soc(param_dim, soc::Float64)
    # 验证SOC范围
    if soc < 0.0 || soc > 1.0
        error("SOC must be in [0, 1], got $soc")
    end
    
    # 负极：SOC↑ → θ_n↑ → cs_n↑
    # θ_n = θ_0_n + SOC × (θ_100_n - θ_0_n)
    theta_n = param_dim.NE.theta_0 + soc * (param_dim.NE.theta_100 - param_dim.NE.theta_0)
    cs0_NE = theta_n * param_dim.NE.cs_max
    
    # 正极：SOC↑ → θ_p↓ → cs_p↓（锂从正极脱出）
    # θ_p = θ_0_p - SOC × (θ_0_p - θ_100_p)
    theta_p = param_dim.PE.theta_0 - soc * (param_dim.PE.theta_0 - param_dim.PE.theta_100)
    cs0_PE = theta_p * param_dim.PE.cs_max
    
    return cs0_NE, cs0_PE
end

"""
    apply_initial_soc!(case::Case, param_dim, soc::Float64)

应用初始SOC设置到case的参数中。

# 参数
- `case::Case`: JuBat Case对象
- `param_dim`: 维度参数结构体
- `soc::Float64`: 目标SOC (0~1)

# 说明
此函数会修改 `param_dim.NE.cs0` 和 `param_dim.PE.cs0`，
然后重新归一化参数到 `case.param` 中。

# 示例
```julia
apply_initial_soc!(case, param_dim, 0.8)  # 设置初始SOC为80%
```
"""
function apply_initial_soc!(case::Case, param_dim, soc::Float64)
    # 计算对应SOC的锂浓度
    cs0_NE, cs0_PE = compute_cs0_from_soc(param_dim, soc)
    
    # 更新维度参数
    param_dim.NE.cs0 = cs0_NE
    param_dim.PE.cs0 = cs0_PE
    
    # 更新归一化参数（重新归一化）
    case.param.NE.cs0 = param_dim.NE.cs0 / param_dim.NE.cs_max
    case.param.PE.cs0 = param_dim.PE.cs0 / param_dim.PE.cs_max
    
    return cs0_NE, cs0_PE
end

"""
    _compute_czm_effective_params(case, param_dim)

计算 CZM 求解所需的有效材料参数。

# 返回
- `E_eff`: 有效弹性模量 [Pa]
- `ν_eff`: 有效泊松比 [-]
- `α_eff`: 有效热膨胀系数 [1/K]
- `β_n`: 负极扩散应变系数 [-]
- `β_p`: 正极扩散应变系数 [-]
"""
function _compute_czm_effective_params(case::Case, param_dim)
    # 有效弹性模量（厚度加权平均）
    E_eff = (param_dim.NE.E * param_dim.NE.thickness + param_dim.PE.E * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    
    # 有效泊松比（厚度加权平均）
    ν_eff = (param_dim.NE.nu * param_dim.NE.thickness + param_dim.PE.nu * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    
    # 有效热膨胀系数（厚度加权平均）
    α_eff = (param_dim.NE.alphaT * param_dim.NE.thickness + param_dim.PE.alphaT * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    
    # 扩散应变系数 β = Ω/3（部分摩尔体积/3）
    β_n = param_dim.NE.Omega / 3.0
    β_p = param_dim.PE.Omega / 3.0
    
    return E_eff, ν_eff, α_eff, β_n, β_p
end

"""
    _compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

计算 CZM 损伤计算所需的单元级应变输入。

# 返回
- `dT_elem`: 每个单元的温度变化 [K]
- `Δsoc_n_elem`: 每个单元的负极 SOC 变化 [-]
- `Δsoc_p_elem`: 每个单元的正极 SOC 变化 [-]
"""
function _compute_czm_strain_inputs(case::Case, variables::Dict, czm_mesh, T_nodes_carry)
    ne = size(czm_mesh.bulk_element, 1)
    param_dim = case.param_dim
    
    # 参考温度（维度值）
    T_ref = param_dim.cell.T0
    T_ref_scale = param_dim.scale.T_ref
    
    # 参考 SOC（归一化值）
    soc_ref_n = case.param.NE.cs0
    soc_ref_p = case.param.PE.cs0
    
    # 初始化输出数组
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)
    
    # 提取温度场（转换为维度值）
    if !isempty(T_nodes_carry) && length(T_nodes_carry) >= czm_mesh.nnode
        # T_nodes_carry 是无量纲温度（T/T_ref）
        for e in 1:ne
            nodes = czm_mesh.bulk_element[e, :]
            # 计算单元平均温度（无量纲）
            T_elem_nd = 0.0
            valid_nodes = 0
            for n in nodes
                if n <= length(T_nodes_carry)
                    T_elem_nd += T_nodes_carry[n]
                    valid_nodes += 1
                end
            end
            if valid_nodes > 0
                T_elem_nd /= valid_nodes
                # 转换为有量纲温度变化 [K]
                T_elem_dim = T_elem_nd * T_ref_scale
                dT_elem[e] = T_elem_dim - T_ref
            end
        end
    end
    
    # 提取 SOC 分布（如果 variables 中有）
    if haskey(variables, "thermal2D element soc_n") && haskey(variables, "thermal2D element soc_p")
        soc_n_elem = variables["thermal2D element soc_n"]
        soc_p_elem = variables["thermal2D element soc_p"]
        
        # 处理数组维度（可能是 ne×1 或 ne×num）
        if isa(soc_n_elem, AbstractMatrix)
            soc_n_elem = soc_n_elem[:, end]
            soc_p_elem = soc_p_elem[:, end]
        end
        
        for e in 1:min(ne, length(soc_n_elem))
            Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
            Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
        end
    end
    
    return dT_elem, Δsoc_n_elem, Δsoc_p_elem
end

"""
    _update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev; debug=false)

更新 CZM 网格的损伤状态。

# 参数
- `czm_mesh`: CZM 网格对象
- `czm_params`: CZM 参数（cohesive）
- `case`: Case 对象
- `variables`: 当前时间步的变量字典
- `T_nodes_carry`: 当前温度场
- `u_czm_prev`: 上一步的 CZM 位移场
- `debug`: 是否输出调试信息

# 返回
- `u_czm`: 更新后的 CZM 位移场
- `converged`: 是否收敛
"""
function _update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev; debug::Bool=false)
    param_dim = case.param_dim
    
    # 计算有效材料参数
    E_eff, ν_eff, α_eff, β_n, β_p = _compute_czm_effective_params(case, param_dim)
    
    # 计算应变输入
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = _compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)
    
    # 调试输出
    if debug
        @printf("  [CZM Debug] dT: max=%.2f K, min=%.2f K\n", maximum(dT_elem), minimum(dT_elem))
        @printf("  [CZM Debug] Δsoc_n: max=%.4f, min=%.4f\n", maximum(Δsoc_n_elem), minimum(Δsoc_n_elem))
        @printf("  [CZM Debug] Δsoc_p: max=%.4f, min=%.4f\n", maximum(Δsoc_p_elem), minimum(Δsoc_p_elem))
        
        # 计算等效应变
        ε_max = maximum(abs.(α_eff .* dT_elem .+ β_n .* Δsoc_n_elem .+ β_p .* Δsoc_p_elem))
        @printf("  [CZM Debug] ε_max=%.2e, α_eff=%.2e, β_n=%.2e, β_p=%.2e\n", ε_max, α_eff, β_n, β_p)
    end
    
    # 外力向量（一般为零）
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)
    
    # 初始化位移（如果没有上一步的值）
    if u_czm_prev === nothing || length(u_czm_prev) != ndof
        u_czm_prev = zeros(Float64, ndof)
    end
    
    # 调用 CZM 求解器
    result = solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, czm_params, param_dim, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=30, tol=1e-6
    )
    
    # 调试输出：分离位移统计
    if debug && result.converged
        δ_n_max = maximum(result.separation_n)
        δ_t_max = maximum(abs.(result.separation_t))
        δ_0_n = czm_params.δ_0_n
        δ_0_t = czm_params.δ_0_t
        @printf("  [CZM Debug] δ_n_max=%.2e (%.1f%% of δ_0_n=%.2e)\n", δ_n_max, 100*δ_n_max/δ_0_n, δ_0_n)
        @printf("  [CZM Debug] δ_t_max=%.2e (%.1f%% of δ_0_t=%.2e)\n", δ_t_max, 100*δ_t_max/δ_0_t, δ_0_t)
        
        D_max = maximum(result.damage)
        @printf("  [CZM Debug] D_max=%.4f%%\n", D_max * 100)
    end
    
    return result.displacement, result.converged
end

"""
    _ensure_multi_spme_layout!(case::Case)

确保多SPMe布局信息已初始化（不改变状态向量）。

当状态向量从上一阶段传递时，case.multi_spme_layout 可能为空。
此函数计算并设置必要的布局信息，使 CallModel_MultiSPMe 能正常工作。
"""
function _ensure_multi_spme_layout!(case::Case)
    if !haskey(case.mesh, "thermal2D")
        error("_ensure_multi_spme_layout! requires thermal2D mesh")
    end
    
    # 获取维度信息
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    
    # 计算单个单元的电化学自由度数
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nel = case.mesh["electrolyte"].nlen
    n_chem = Nrn + Nrp + Nel
    
    # 设置布局信息
    if isempty(case.multi_spme_layout)
        empty!(case.multi_spme_layout)
    end
    case.multi_spme_layout["ne"] = ne
    case.multi_spme_layout["n_chem"] = n_chem
    case.multi_spme_layout["nT"] = nT
    case.multi_spme_layout["n_total"] = ne * n_chem + nT
    case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
    case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)
    
    return nothing
end

"""
    plot_cycling_results(result; save_path="output/")

绘制循环结果图。
"""
function plot_cycling_results(result::CyclingResult; save_path::String="output/")
    isdir(save_path) || mkdir(save_path)
    
    cycles = result.cycle_idx
    
    # 图1: 容量衰减
    p1 = plot(cycles, result.capacity_discharge,
              xlabel="Cycle Number", ylabel="Discharge Capacity (Ah)",
              label="Discharge", linewidth=2, marker=:circle,
              title="Capacity Fade")
    plot!(p1, cycles, result.capacity_charge,
          label="Charge", linewidth=2, marker=:square, linestyle=:dash)
    savefig(p1, joinpath(save_path, "cycling_capacity.png"))
    
    # 图2: 损伤演化
    p2 = plot(cycles, result.D_max .* 100,
              xlabel="Cycle Number", ylabel="Damage (%)",
              label="D_max", linewidth=2, color=:red,
              title="Damage Evolution")
    plot!(p2, cycles, result.D_mean .* 100,
          label="D_mean", linewidth=2, color=:blue, linestyle=:dash)
    savefig(p2, joinpath(save_path, "cycling_damage.png"))
    
    # 图3: 库伦效率
    p3 = plot(cycles, result.coulombic_efficiency,
              xlabel="Cycle Number", ylabel="Coulombic Efficiency (%)",
              label="CE", linewidth=2, marker=:circle,
              title="Coulombic Efficiency",
              ylims=(95, 105))
    savefig(p3, joinpath(save_path, "cycling_efficiency.png"))
    
    # 图4: 温度
    p4 = plot(cycles, result.T_max,
              xlabel="Cycle Number", ylabel="T_max (K)",
              label="Max Temperature", linewidth=2, marker=:circle,
              title="Temperature History")
    savefig(p4, joinpath(save_path, "cycling_temperature.png"))
    
    # 组合图
    p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900))
    savefig(p_all, joinpath(save_path, "cycling_summary.png"))
    
    println("✓ 循环结果图已保存至 $save_path")
    
    return p_all
end