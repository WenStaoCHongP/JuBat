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
    soh::Vector{Float64}                # SOH (State of Health)
    # 详细结果（可选）
    cycle_results::Vector{CycleResult}
    
    # 最终状态
    final_czm_mesh::Any                 # 最终的CZM网格（包含损伤状态）
    initial_capacity::Float64           # 初始容量（第一个循环的放电容量）
    soh_terminated::Bool                # 是否因SOH低于阈值终止
end

function CyclingResult(n_cycles::Int)
    CyclingResult(
        0,
        Int[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Int[], Float64[], Float64[],
        CycleResult[],
        nothing,          # final_czm_mesh
        0.0,              # initial_capacity
        false             # soh_terminated
    )
end


# ========================================================================
# 2. 单阶段求解器
# ========================================================================

function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64, I_current::Float64, V_limit::Float64, initial_state::Dict; czm_mesh=nothing, czm_params=nothing, dt_range::Vector{Float64}=[1.0, 10.0])
    result = PhaseResult()
    result.phase_type = phase_type
    result.t_start = get(initial_state, "t_global", 0.0)

    if phase_type == PHASE_REST
        I_current = 0.0
    end

    # 记录阶段开始损伤
    D_max_init = czm_mesh !== nothing ? maximum(s.D for s in czm_mesh.damage_states) : 0.0
    D_mean_init = czm_mesh !== nothing ? mean(s.D for s in czm_mesh.damage_states) : 0.0

    # 配置本阶段求解参数
    case.opt.Current = _ -> I_current
    case.opt.time = [0.0, t_max]
    case.opt.dt = dt_range

    old_v_l = case.param.cell.v_l
    old_v_h = case.param.cell.v_h
    try
        if phase_type == PHASE_DISCHARGE
            case.param.cell.v_l = V_limit
            case.param.cell.v_h = 1.0e6
        elseif phase_type == PHASE_CHARGE
            case.param.cell.v_l = -1.0e6
            case.param.cell.v_h = V_limit
        else
            case.param.cell.v_l = -1.0e6
            case.param.cell.v_h = 1.0e6
        end

        solve_result = Solve(case; initial_state=initial_state, return_final_state=true)

        duration = begin
            time_hist = get(solve_result, "time [s]", Float64[])
            isempty(time_hist) ? t_max : time_hist[end]
        end

        final_state = get(solve_result, "final_state", Dict{String, Any}())
        final_state["t_global"] = result.t_start + duration

        # 阶段末调用 CZM 求解器（CzmSolve.jl）更新损伤，避免在 CycleSolver 重复定义求解循环。
        if czm_mesh !== nothing && czm_params !== nothing && haskey(final_state, "y")
            y_end = final_state["y"]
            T_nodes_nd = get(final_state, "T_nodes", Float64[])
            t_end_nd = duration / case.param.scale.t0
            _, _, _, vars_end, _ = CallModel(case, y_end, t_end_nd, jacobi="update")
            try
                _update_czm_damage!(czm_mesh, czm_params, case, vars_end, T_nodes_nd, nothing)
            catch e
                @debug "Phase-end CZM update failed: $e"
            end
        end

        phase_data = _postprocess_phase_result(
            case, phase_type, solve_result, initial_state,
            I_current, result.t_start,
            D_max_init, D_mean_init, czm_mesh
        )

        result.t_end = result.t_start + phase_data["duration"]
        result.duration = phase_data["duration"]
        result.V_start = phase_data["V_start"]
        result.V_end = phase_data["V_end"]
        result.capacity = phase_data["capacity"]
        result.terminated_by = phase_data["terminated_by"]
        result.T_max = phase_data["T_max"]
        result.T_mean_end = phase_data["T_mean_end"]
        result.D_max = phase_data["D_max"]
        result.D_mean = phase_data["D_mean"]
        result.ΔD_max = phase_data["ΔD_max"]
        result.final_state = phase_data["final_state"]
    finally
        case.param.cell.v_l = old_v_l
        case.param.cell.v_h = old_v_h
    end

    return result
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
function solve_cycling(case::Case, cycle_opt::CycleOption, czm_mesh=nothing;verbose::Bool=true, save_detailed::Bool=false)
    
    n_cycles = cycle_opt.n_cycles
    result = CyclingResult(n_cycles)

    # 应用初始SOC设置
    soc_init = cycle_opt.SOC_init
    if soc_init >= 0.0 && soc_init <= 1.0
        cs0_NE, cs0_PE = apply_initial_soc!(case, case.param_dim, soc_init)
        if verbose
            @printf("  初始SOC: %.1f%%\n", soc_init * 100)
            @printf("    → 负极cs0: %.1f mol/m³\n", cs0_NE)
            @printf("    → 正极cs0: %.1f mol/m³\n", cs0_PE)
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
    czm_params = case.param_dim.cohesive
    
    # SOH监控参数
    soh_threshold = case.opt.czm_soh_threshold
    initial_capacity = 0.0  # 初始放电容量（第一个循环结束后设置）
    soh_terminated = false  # SOH终止标志

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
                # no-op
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
        if cycle_opt.reset_T_before_charge
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
        
        # ============ 循环汇总与后处理 ============
        _postprocess_cycle_result!(cycle_result, charge_result, discharge_result, rest1_result, rest2_result, czm_mesh)
        _append_cycle_result!(result, cycle, cycle_result; save_detailed=save_detailed)
        initial_capacity, current_soh = _update_soh_and_capacity!(result, cycle, cycle_result, initial_capacity)

        if verbose
            _print_cycle_summary(cycle, cycle_result, current_soh)
        end

        should_stop, soh_hit = _check_cycle_termination(cycle, cycle_result, czm_mesh, current_soh, soh_threshold; verbose=verbose)
        if soh_hit
            soh_terminated = true
            result.soh_terminated = true
        end
        if should_stop
            break
        end
    end
    
    result.final_czm_mesh = czm_mesh
    if verbose
        _print_cycling_summary(result, initial_capacity, soh_terminated)
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
    if length(T_nodes_carry) >= czm_mesh.nnode
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
    
    return dT_elem, Δsoc_n_elem, Δsoc_p_elem
end

"""
    _update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)

更新 CZM 网格的损伤状态。

使用牛顿-拉弗森迭代求解力学平衡方程，通过载荷子步法处理软化收敛问题。

# 参数
- `czm_mesh`: CZM 网格对象
- `czm_params`: CZM 参数（cohesive）
- `case`: Case 对象
- `variables`: 当前时间步的变量字典
- `T_nodes_carry`: 当前温度场
- `u_czm_prev`: 上一步的 CZM 位移场

# 返回
- `u_czm`: 更新后的 CZM 位移场
- `converged`: 是否收敛
"""
function _update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
    param_dim = case.param_dim

    # 同步CZM模型选项（model1 or mix）
    czm_params.czm_model = case.opt.czm_model
    
    # 计算有效材料参数
    E_eff, ν_eff, α_eff, β_n, β_p = _compute_czm_effective_params(case, param_dim)
    
    # 计算应变输入
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = _compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)
    
    # 外力向量（一般为零）
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)
    
    # 初始化位移（如果没有上一步的值）
    if u_czm_prev === nothing || length(u_czm_prev) != ndof
        u_czm_prev = zeros(Float64, ndof)
    end
    
    # 调用 CZM 求解器（可选迭代方式）
    iter_method = case.opt.czm_iter_method
    max_iter = case.opt.czm_max_iter
    tol = case.opt.czm_tol
    n_load_steps = case.opt.czm_load_steps
    arc_length_alpha = case.opt.czm_arc_length_alpha

    result, updated_czm_mesh = solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, czm_params, param_dim, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha, iter_method=iter_method
    )

    # Keep caller mesh damage states in sync for subsequent thermo-mechanical coupling.
    czm_mesh.damage_states = updated_czm_mesh.damage_states
    
    return result.displacement, result.converged
end

"""
    _ensure_multi_spme_layout!(case::Case)

确保多SPMe布局信息已初始化（不改变状态向量）。

当状态向量从上一阶段传递时，case.multi_spme_layout 可能为空。
此函数计算并设置必要的布局信息，使 CallModel_MultiSPMe 能正常工作。
"""
function _ensure_multi_spme_layout!(case::Case)
    # 获取维度信息
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    
    # 计算单个单元的电化学自由度数
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nel = case.mesh["electrolyte"].nlen
    n_chem = Nrn + Nrp + Nel
    
    # 设置布局信息
    empty!(case.multi_spme_layout)
    case.multi_spme_layout["ne"] = ne
    case.multi_spme_layout["n_chem"] = n_chem
    case.multi_spme_layout["nT"] = nT
    case.multi_spme_layout["n_total"] = ne * n_chem + nT
    case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
    case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)
    
    return nothing
end
