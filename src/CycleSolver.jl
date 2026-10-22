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
    terminated_by::Symbol               # :time / :voltage目标 / 实际异常原因
    
    # 热
    T_max::Float64                      # 最高温度 (K)
    T_mean_end::Float64                 # 结束时平均温度 (K)
    
    # 损伤
    D_max::Float64                      # 最大损伤
    D_mean::Float64                     # 平均损伤
    ΔD_max::Float64                     # 损伤增量（相对于阶段开始）
    
    # 最终状态（用于传递到下一阶段）
    final_state::Union{NamedTuple, Nothing}

    # Raw solve result (for CSV export, only when save_detailed=true)
    solve_result::Any   # Dict from Solve(), or nothing
end

# 默认构造函数
function PhaseResult()
    PhaseResult(
        PHASE_REST, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, :none,
        0.0, 0.0,
        0.0, 0.0, 0.0,
        nothing,
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
    rest1::Union{PhaseResult, Nothing}
    discharge::PhaseResult
    rest2::Union{PhaseResult, Nothing}
    
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

    # CZM快照（用于CSV导出）
    czm_snapshots::Vector{CZMSnapshot}

    # 最终状态
    final_mech::Any                     # 最终的 MechState（包含损伤/位移状态）
    initial_capacity::Float64           # 初始容量（第一个循环的放电容量）
    soh_terminated::Bool                # 是否因SOH低于阈值终止
end

function CyclingResult(n_cycles::Int)
    CyclingResult(
        0,
        Int[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Int[], Float64[], Float64[],
        CycleResult[],
        CZMSnapshot[],    # czm_snapshots
        nothing,          # final_mech
        0.0,              # initial_capacity
        false             # soh_terminated
    )
end


# ========================================================================
# 2. 单阶段求解器
# ========================================================================

function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64, initial_state::Union{NamedTuple,Nothing}; dt_range::Vector{Float64}=[1.0, 10.0], czm_snapshots::Union{Vector{CZMSnapshot},Nothing}=nothing, czm_cycle::Int=1, czm_phase::String=string(phase_type), chemical_step_observer::Union{Nothing,Function}=nothing,variables::Dict{String,Any}=Dict{String,Any}())
    configured_current=case.opt.Current
    configured_dt=case.opt.dt
    configured_time=case.opt.time
    if phase_type == PHASE_DISCHARGE
        case.opt.Current=configured_current
    elseif phase_type == PHASE_CHARGE
        case.opt.Current=t -> -configured_current(t)
    elseif phase_type == PHASE_REST
        case.opt.Current=t -> 0.0
    end
    try
    case.opt.dt=dt_range
    case.opt.time=[0.0 t_max]
    result = PhaseResult()
    result.phase_type = phase_type
    result.t_start = initial_state === nothing ? 0.0 : initial_state.t_global


    phase_state=initial_state
    D_max_init=phase_state===nothing || phase_state.mech===nothing ? 0.0 : maximum(s.D for s in phase_state.mech.damage_states)
    D_mean_init=phase_state===nothing || phase_state.mech===nothing ? 0.0 : mean(s.D for s in phase_state.mech.damage_states)
    solve_result = Solve(case; initial_state=phase_state,variables=variables, return_final_state=true,
                         czm_snapshots=czm_snapshots, czm_cycle=czm_cycle,
                         czm_phase=czm_phase, phase_type=phase_type,
                         chemical_step_observer=chemical_step_observer)

    reason = solve_result["termination_reason"]
    reason in ("time_limit", "phase_target_low", "phase_target_high") ||
        error("$(phase_type) stopped with $(reason); refusing to continue cycling")

    # CZM 损伤已由 Solve 主循环每步更新（Solve.jl:270-285），此处不再冗余调用。
    # 旧版在此处调用 update_czm_damage!(..., u_czm_prev=nothing) 会导致：
    #   1. 位移场从零重解，产生不一致的应力-位移解
    #   2. 每阶段额外浪费 ~12s 计算

    phase_data = postprocess_phase_result(
        case, phase_type, solve_result, initial_state,
        result.t_start,
        D_max_init, D_mean_init
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
    variables["accepted"]=deepcopy(phase_data["final_state"])
    result.final_state = deepcopy(variables["accepted"])
    result.solve_result = solve_result

    return result
    finally
        case.opt.Current=configured_current
        case.opt.dt=configured_dt
        case.opt.time=configured_time
    end
end


# ========================================================================
# 3. 循环求解器
# ========================================================================

"""
    solve_cycling(case, cycle_opt, mechanical_seed; verbose=true)

执行完整的充放电循环仿真。

# 参数
- `case`: JuBat Case对象
- `cycle_opt`: CycleOption 循环参数
- `mechanical_seed`: 仅用于新会话初始化的只读力学种子；入口复制，推进只读 accepted

# 返回
- `CyclingResult`: 循环结果汇总
"""
function solve_cycling(case::Case, cycle_opt::CycleOption, mechanical_seed=nothing;verbose::Bool=true, save_detailed::Bool=false,variables::Dict{String,Any}=Dict{String,Any}())

    if case.opt.czm.enabled
        case.mesh["czm"] # Enabled mechanics requires its real topology.
    end
    n_cycles = cycle_opt.n_cycles
    result = CyclingResult(n_cycles)

    # 2026-09-26 逐元素全历史导出：JUBAT_SNAPSHOT_CYCLES="1,5,10,30,50,100" 指定
    # 采集逐机械步 CZMSnapshot 的循环号；其余循环不采集（仅逐循环/逐相位汇总 CSV）。
    # 与 save_detailed 解耦：后者仍全程采集并额外保留 solve_result（内存不可行于百循环）。
    snapshot_cycles = czm_snapshot_cycles()
    want_snapshots = haskey(case.mesh,"czm") && (save_detailed || !isempty(snapshot_cycles))

    # CZM snapshots vector (shared across all phases)
    czm_snaps = want_snapshots ? CZMSnapshot[] : nothing

    # 应用初始SOC设置
    soc_init = cycle_opt.SOC_init
    cs0_NE, cs0_PE = apply_initial_soc!(case, case.param_dim, soc_init)
    if verbose
        @printf("  初始SOC: %.1f%%\n", soc_init * 100)
        @printf("    → 负极cs0: %.1f mol/m³\n", cs0_NE)
        @printf("    → 正极cs0: %.1f mol/m³\n", cs0_PE)
    end

    # 初始化
    if verbose
        println("="^60)
        println("开始充放电循环仿真")
        println("="^60)
        @printf("  循环次数: %d\n", n_cycles)
        @printf("  初始SOC: %.1f%%\n", soc_init * 100)
        println("  循环顺序: 放电 → 静置 → 充电 → 静置")
        @printf("  放电: %.0fs, 截止 %.2fV\n", cycle_opt.t_discharge, case.param.cell.v_l)
        @printf("  充电: %.0fs, 截止 %.2fV\n", cycle_opt.t_charge, case.param.cell.v_h)
        @printf("  静置: %.0fs + %.0fs\n", cycle_opt.t_rest1, cycle_opt.t_rest2)
    end
    
    # 初始状态
    current_state=(y=copy(case.opt.per_element_spme ? ModelInitialisation_MultiSPMe(case) : ModelInitialisation(case)),
        mech=mechanical_seed===nothing ? (haskey(case.mesh,"czm") ? MechState(case.mesh["czm"]) : nothing) : deepcopy(mechanical_seed),t_global=0.0)
    variables["accepted"]=deepcopy(current_state)
    
    
    # SOH监控参数
    soh_threshold = case.opt.czm.soh_threshold
    initial_capacity = 0.0  # 初始放电容量（第一个循环结束后设置）
    soh_terminated = false  # SOH终止标志

    # 循环主体（cycle_snaps 按循环号门控：仅选定循环采集全历史快照）
    for cycle in 1:n_cycles
        cycle_snaps = czm_snaps === nothing ? nothing :
            (save_detailed || cycle in snapshot_cycles ? czm_snaps : nothing)
        if verbose
            println("\n" * "-"^40)
            @printf("循环 %d/%d\n", cycle, n_cycles)
            println("-"^40)
        end
        
        cycle_result = CycleResult(
            cycle,
            PhaseResult(), nothing, PhaseResult(), nothing,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0,
            0.0
        )
        
        # ============ 阶段1: 放电 ============
        # 循环顺序：放电 → 静置 → 充电 → 静置
        # （适配Jellyroll参数的初始高SOC状态）
        if verbose
            print("  [放电] ")
        end
        
        discharge_result = solve_phase(
            case, PHASE_DISCHARGE,
            cycle_opt.t_discharge,
            current_state;
            dt_range=cycle_opt.dt_cycle,variables=variables,
            czm_snapshots=cycle_snaps, czm_cycle=cycle, czm_phase="discharge"
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
        # 静置阶段：电化学状态继承上一步，以0电流继续锂扩散
        if cycle_opt.t_rest1 > 0
            # 静置时间 > 0：执行静置阶段，继续锂扩散过程
            if verbose
                print("  [静置1] 状态继承+锂扩散 ")
            end
            
            rest1_result = solve_phase(
                case, PHASE_REST,
                cycle_opt.t_rest1,
                current_state;  # 继承上一步状态
                    dt_range=cycle_opt.dt_cycle,variables=variables,
                czm_snapshots=cycle_snaps, czm_cycle=cycle, czm_phase="rest1"
            )
            cycle_result.rest1 = rest1_result
            current_state = rest1_result.final_state
            
            if verbose
                @printf("%.1fs, T_max=%.2fK\n", rest1_result.duration, rest1_result.T_max)
            end
        else
            # 静置时间 = 0：跳过静置阶段，直接继承上一阶段状态
            # current_state 保持不变（无扩散过程）
            
            if verbose
                println("  [静置1] 跳过 (t=0，无扩散)")
            end
        end
        
        # ============ 阶段3: 充电 ============
        if verbose
            print("  [充电] ")
        end
        
        charge_result = solve_phase(
            case, PHASE_CHARGE,
            cycle_opt.t_charge,
            current_state;
            dt_range=cycle_opt.dt_cycle,variables=variables,
            czm_snapshots=cycle_snaps, czm_cycle=cycle, czm_phase="charge"
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
                current_state;  # 继承上一步状态
                    dt_range=cycle_opt.dt_cycle,variables=variables,
                czm_snapshots=cycle_snaps, czm_cycle=cycle, czm_phase="rest2"
            )
            cycle_result.rest2 = rest2_result
            current_state = rest2_result.final_state
            if verbose
                @printf("%.1fs, T_max=%.2fK\n", rest2_result.duration, rest2_result.T_max)
            end
        else
            # 静置时间 = 0：跳过静置阶段，无扩散过程
            # current_state 保持不变（无扩散过程）
            if verbose
                println("  [静置2] 跳过 (t=0，无扩散)")
            end
        end
        
        # ============ 循环汇总与后处理 ============
        postprocess_cycle_result!(cycle_result, charge_result, discharge_result,
                                  cycle_result.rest1, cycle_result.rest2, current_state.mech)
        append_cycle_result!(result, cycle, cycle_result; save_detailed=save_detailed)
        initial_capacity, current_soh = update_soh_and_capacity!(result, cycle, cycle_result, initial_capacity)

        if verbose
            print_cycle_summary(cycle, cycle_result, current_soh)
        end

        should_stop, soh_hit = check_cycle_termination(cycle, cycle_result, current_state.mech, current_soh, soh_threshold; verbose=verbose)
        if soh_hit
            soh_terminated = true
            result.soh_terminated = true
        end
        if should_stop
            break
        end
    end

    # Attach CZM snapshots to result
    if czm_snaps !== nothing
        result.czm_snapshots = czm_snaps
    end

    result.final_mech = deepcopy(current_state.mech)
    if verbose
        print_cycling_summary(result, initial_capacity, soh_terminated)
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

# 2026-09-26 逐元素全历史导出：解析 JUBAT_SNAPSHOT_CYCLES（逗号分隔循环号，空=不采集）
function czm_snapshot_cycles()
    raw = get(ENV, "JUBAT_SNAPSHOT_CYCLES", "")
    isempty(raw) && return Set{Int}()
    return Set(parse(Int, strip(x)) for x in split(raw, ","))
end
