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
    V_init = get(initial_state, "V", 3.7)
    
    result.V_start = V_init
    
    # 根据模式选择求解
    # 简化：调用内部求解循环
    phase_result_data = _solve_phase_internal(
        case, phase_type, t_max, I_current, V_limit,
        y0, T_nodes, czm_mesh, czm_params, dt_range)
    
    # 填充结果
    result.t_end = result.t_start + phase_result_data["duration"]
    result.duration = phase_result_data["duration"]
    result.V_end = phase_result_data["V_end"]
    result.capacity = phase_result_data["capacity"]
    result.terminated_by = phase_result_data["terminated_by"]
    result.T_max = phase_result_data["T_max"]
    result.T_mean_end = phase_result_data["T_mean_end"]
    result.D_max = phase_result_data["D_max"]
    result.D_mean = phase_result_data["D_mean"]
    result.ΔD_max = phase_result_data["ΔD_max"]
    result.final_state = phase_result_data["final_state"]
    
    return result
end

"""
    _solve_phase_internal(case, phase_type, t_max, I_current, V_limit, 
                          y0, T_nodes, czm_mesh, czm_params, dt_range)

内部阶段求解实现。
"""
function _solve_phase_internal(case::Case, phase_type::PhaseType, 
                               t_max::Float64, I_current::Float64, V_limit::Float64,
                               y0, T_nodes, czm_mesh, czm_params, dt_range)
    
    # 时间缩放
    t0_scale = case.param.scale.t0
    dt_min = dt_range[1] / t0_scale
    dt_max = dt_range[2] / t0_scale
    t_end_nd = t_max / t0_scale
    
    # 初始化
    if y0 === nothing
        # 使用标准初始化
        multi_spme = case.opt.model == "SPMe" && 
                     hasproperty(case.opt, :per_element_spme) && 
                     case.opt.per_element_spme
        if multi_spme
            y0 = ModelInitialisation_MultiSPMe(case)
        else
            y0 = ModelInitialisation(case)
        end
    end
    
    # 初始化温度场
    if T_nodes !== nothing && haskey(case.mesh, "thermal2D")
        # 使用传入的温度场
        T_nodes_carry = copy(T_nodes)
    elseif haskey(case.mesh, "thermal2D")
        T_nodes_carry = fill(case.param.cell.T0, case.mesh["thermal2D"].nlen)
    else
        T_nodes_carry = Float64[]
    end
    
    # Crank-Nicolson 参数
    theta = case.opt.solveType == "Crank-Nicolson" ? 0.5 : 
            (case.opt.solveType == "forward" ? 0.0 : 1.0)
    
    # 初始调用
    t = 0.0
    dt = dt_min
    M_old, K_old, F_old, variables, y_phi = CallModel(case, y0, t, jacobi="update")
    
    # 初始化变量
    V_current = variables["cell voltage"] * case.param.scale.phi
    capacity = 0.0
    T_max_phase = haskey(variables, "T_nodes") ? 
                  maximum(variables["T_nodes"]) * case.param_dim.scale.T_ref : 
                  case.param_dim.cell.T0
    
    # 损伤初始状态
    D_max_init = czm_mesh !== nothing ? maximum(s.D for s in czm_mesh.damage_states) : 0.0
    D_mean_init = czm_mesh !== nothing ? mean(s.D for s in czm_mesh.damage_states) : 0.0
    
    # 初始求解步
    vc = 1:size(M_old, 1)
    dt_init = 1e-8
    y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
    y_old = vcat(y_c, y_phi)
    
    terminated_by = :time
    t_actual = 0.0
    
    # 电压容差（处理数值边界问题）
    V_tolerance = 0.05  # 50mV 容差
    
    # 检查初始电压是否已接近截止条件
    if phase_type == PHASE_CHARGE && V_current >= (V_limit - V_tolerance)
        @info "充电阶段：初始电压已接近截止" V_current=V_current V_limit=V_limit
        terminated_by = :voltage
        return Dict(
            "duration" => 0.0, "V_end" => V_current, "capacity" => 0.0,
            "terminated_by" => :voltage,
            "T_max" => T_max_phase,
            "T_mean_end" => !isempty(T_nodes_carry) ? mean(T_nodes_carry) * case.param_dim.scale.T_ref : case.param_dim.cell.T0,
            "D_max" => D_max_init, "D_mean" => D_mean_init, "ΔD_max" => 0.0,
            "final_state" => Dict("y" => y_old, "T_nodes" => T_nodes_carry, "V" => V_current, "t_global" => 0.0)
        )
    elseif phase_type == PHASE_DISCHARGE && V_current <= (V_limit + V_tolerance)
        @info "放电阶段：初始电压已接近截止" V_current=V_current V_limit=V_limit
        terminated_by = :voltage
        return Dict(
            "duration" => 0.0, "V_end" => V_current, "capacity" => 0.0,
            "terminated_by" => :voltage,
            "T_max" => T_max_phase,
            "T_mean_end" => !isempty(T_nodes_carry) ? mean(T_nodes_carry) * case.param_dim.scale.T_ref : case.param_dim.cell.T0,
            "D_max" => D_max_init, "D_mean" => D_mean_init, "ΔD_max" => 0.0,
            "final_state" => Dict("y" => y_old, "T_nodes" => T_nodes_carry, "V" => V_current, "t_global" => 0.0)
        )
    end
    
    # 主循环
    while t < t_end_nd
        # 更新温度影响
        if case.opt.thermal_enabled && !isempty(T_nodes_carry)
            case.param.cell.T0 = mean(T_nodes_carry)
        end
        
        # 检查电压是否接近截止（在调用 CallModel 前）
        if phase_type == PHASE_CHARGE && V_current >= (V_limit - V_tolerance)
            terminated_by = :voltage
            break
        elseif phase_type == PHASE_DISCHARGE && V_current <= (V_limit + V_tolerance)
            terminated_by = :voltage
            break
        end
        
        # 电化学步 - 使用 try-catch 捕获电压越界错误
        local M_new, K_new, F_new, y_phi_new
        try
            M_new, K_new, F_new, variables, y_phi_new = CallModel(case, y_old, t, jacobi="update")
        catch e
            # 如果是电压越界错误，优雅终止
            if occursin("voltage out of bounds", string(e))
                @warn "电压越界，提前终止阶段" phase=phase_type t=t*t0_scale
                terminated_by = :voltage
                break
            else
                rethrow(e)
            end
        end
        y_phi = y_phi_new
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
    
    # CZM损伤计算（如果有）
    D_max_end = D_max_init
    D_mean_end = D_mean_init
    if czm_mesh !== nothing && czm_params !== nothing && t_actual > 0
        # 计算最终损伤状态
        stats = get_damage_statistics(czm_mesh)
        D_max_end = stats.max_D
        D_mean_end = stats.mean_D
    end
    
    # 构建最终状态
    final_state = Dict(
        "y" => y_old,
        "T_nodes" => T_nodes_carry,
        "V" => V_current,
        "t_global" => t_actual
    )
    
    return Dict(
        "duration" => t_actual,
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
    
    # 初始化
    if verbose
        println("="^60)
        println("开始充放电循环仿真")
        println("="^60)
        @printf("  循环次数: %d\n", n_cycles)
        println("  循环顺序: 放电 → 静置 → 充电 → 静置")
        @printf("  放电: %.0fs (%.1fC), 截止 %.2fV\n", 
                cycle_opt.t_discharge, cycle_opt.I_discharge/5.0, cycle_opt.V_lower)
        @printf("  充电: %.0fs (%.1fC), 截止 %.2fV\n", 
                cycle_opt.t_charge, cycle_opt.I_charge/5.0, cycle_opt.V_upper)
        @printf("  静置: %.0fs + %.0fs\n", cycle_opt.t_rest1, cycle_opt.t_rest2)
    end
    
    # 注意：初始SOC使用参数文件（如Jellyroll.jl）中的设置
    # 不在此处覆盖，以保持与参数文件的一致性
    
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
        end
        
        # ============ 阶段2: 静置1 ============
        if verbose
            print("  [静置1] ")
        end
        
        rest1_result = solve_phase(
            case, PHASE_REST,
            cycle_opt.t_rest1,
            0.0,  # 无电流
            0.0,  # 无电压限制
            current_state;
            czm_mesh=czm_mesh,
            czm_params=czm_params,
            dt_range=cycle_opt.dt_cycle
        )
        cycle_result.rest1 = rest1_result
        current_state = rest1_result.final_state
        
        if verbose
            @printf("%.1fs, T_max=%.2fK\n", rest1_result.duration, rest1_result.T_max)
        end
        
        # ============ 阶段3: 充电 ============
        if verbose
            print("  [充电] ")
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
        if verbose
            print("  [静置2] ")
        end
        
        rest2_result = solve_phase(
            case, PHASE_REST,
            cycle_opt.t_rest2,
            0.0,
            0.0,
            current_state;
            czm_mesh=czm_mesh,
            czm_params=czm_params,
            dt_range=cycle_opt.dt_cycle
        )
        cycle_result.rest2 = rest2_result
        current_state = rest2_result.final_state
        
        if verbose
            @printf("%.1fs, T_max=%.2fK\n", rest2_result.duration, rest2_result.T_max)
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
