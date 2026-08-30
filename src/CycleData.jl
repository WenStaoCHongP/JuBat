# ========================================================================
# 循环数据在线采集（为 CSV 导出提供温度/SOC 等结构化快照）
# ========================================================================
"""
    TimeStepData - 单个时间步的数据
"""
struct TimeStepData
    time::Float64               # 物理时间 (s)
    phase::PhaseType            # 阶段类型
    V::Float64                  # 电压 (V)
    I::Float64                  # 电流 (A)
    T_nodes::Vector{Float64}    # 节点温度 (K)
    T_max::Float64              # 最高温度 (K)
    T_mean::Float64             # 平均温度 (K)
    soc_n::Vector{Float64}      # 负极各单元SOC
    soc_p::Vector{Float64}      # 正极各单元SOC
    soc_mean::Float64           # 平均SOC
end

"""
    CycleExportData - 完整循环的导出数据
"""
struct CycleExportData
    cycle_idx::Int
    timesteps::Vector{TimeStepData}
    node_coords::Matrix{Float64}    # 节点坐标 (nT × 2)
    element_connectivity::Matrix{Int} # 单元连接 (ne × 4)
    ne::Int                         # 单元数
    nT::Int                         # 节点数
end

function solve_phase_with_export(case::Case, phase_type::PhaseType, t_max::Float64, I_current::Float64, V_limit::Float64, initial_state::Union{Nothing, Dict};czm_mesh=nothing, czm_params=nothing,dt_range::Vector{Float64}=[1.0, 10.0],export_interval::Int=1)
    if case.opt.czm.enabled || czm_mesh !== nothing || czm_params !== nothing
        throw(ArgumentError(
            "solve_phase_with_export does not implement CZM damage evolution; use solve_phase or solve_cycling for CZM simulations"
        ))
    end
    
    timestep_data = TimeStepData[]
    
    # 时间缩放
    t0_scale = case.param.scale.t0
    dt_min = dt_range[1] / t0_scale
    dt_max = dt_range[2] / t0_scale
    t_end_nd = t_max / t0_scale
    T_ref = case.param_dim.scale.T_ref
    phi_scale = case.param.scale.phi
    I_scale = case.param.scale.I_typ
    
    # 静置阶段特殊处理
    if phase_type == PHASE_REST
        I_current = 0.0
    end
    # 检测多SPMe模式
    multi_spme = case.opt.model == "SPMe" &&case.opt.per_element_spme &&case.opt.thermalmodel == "distributed2D"
    
    # 仅显式传入 nothing 时初始化；状态字典必须包含有效状态向量
    if initial_state === nothing
        T_nodes_init = nothing
        if multi_spme
            y0 = ModelInitialisation_MultiSPMe(case)
        else
            y0 = ModelInitialisation(case)
        end
    else
        haskey(initial_state, "y") || throw(ArgumentError("initial_state must contain key \"y\""))
        y_from_state = initial_state["y"]
        y_from_state isa AbstractArray || throw(ArgumentError("initial_state[\"y\"] must be an array"))
        y0 = vec(y_from_state)
        T_nodes_init = get(initial_state, "T_nodes", nothing)

        if multi_spme
            ne = size(case.mesh["thermal2D"].element, 1)
            nT = case.mesh["thermal2D"].nlen
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Nel = case.mesh["electrolyte"].nlen
            n_chem = Nrn + Nrp + Nel
            expected_multi_len = ne * n_chem + nT

            if case.layout === nothing
                case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
            end

            length(y0) == expected_multi_len || throw(DimensionMismatch(
                "external state length $(length(y0)) does not match multi-SPMe layout length $expected_multi_len"
            ))
        end
    end

    # 初始化温度场
    if T_nodes_init !== nothing
        T_nodes_carry = copy(T_nodes_init)
        if multi_spme
            thermal_range = case.layout.thermal_range
            y0[thermal_range] .= T_nodes_carry
        end
    elseif multi_spme
        thermal_range = case.layout.thermal_range
        T_nodes_carry = y0[thermal_range]
    else
        nT_mesh = case.mesh["thermal2D"].nlen
        T_nodes_carry = y0[(end - nT_mesh + 1):end]
    end
    
    # 时间离散格式必须显式受支持
    if case.opt.solveType == "Crank-Nicolson"
        theta = 0.5
    elseif case.opt.solveType == "forward"
        theta = 0.0
    elseif case.opt.solveType == "backward"
        theta = 1.0
    else
        throw(ArgumentError("unsupported solveType: $(repr(case.opt.solveType))"))
    end
    
    # 设置电流
    case.opt.Current = x -> I_current
    
    # 初始调用
    t = 0.0
    dt = dt_min
    
    M_old, K_old, F_old, variables, y_phi = CallModel(case, y0, t, jacobi="update")
    
    V_current = variables["cell voltage"] * phi_scale
    V_start_actual = V_current
    capacity = 0.0
    T_max_phase = maximum(T_nodes_carry) * T_ref
    
    # 直接使用一致的初始状态，避免导出首步出现非物理温度尖峰。
    vc = 1:size(M_old, 1)
    y_old = vcat(copy(y0[vc]), y_phi)
    
    terminated_by = :time
    t_actual = 0.0
    step_count = 0
    if initial_state === nothing
        global_time = 0.0
    else
        haskey(initial_state, "t_global") || throw(ArgumentError("initial_state must contain key \"t_global\""))
        t_global = initial_state["t_global"]
        t_global isa Real || throw(ArgumentError("initial_state[\"t_global\"] must be a real number"))
        isfinite(t_global) || throw(ArgumentError("initial_state[\"t_global\"] must be finite"))
        global_time = Float64(t_global)
    end
    T_nodes_prev_export = copy(T_nodes_carry)
    
    # 获取热网格信息
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    # Precompute per-element areas for physically meaningful temperature averaging.
    elem_area = zeros(Float64, ne)
    for g in eachindex(mesh_th.gs.weight)
        e = mesh_th.gs.ele[g]
        elem_area[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    area_sum = sum(elem_area)
    
    # 主循环
    while t < t_end_nd
        step_count += 1
        
        # 更新温度影响
        if case.opt.thermal_enabled
            case.param.cell.T0 = mean(T_nodes_carry)
        end
        
        # 电化学步
        M_new, K_new, F_new, variables, y_phi = CallModel(case, y_old, t, jacobi="update")
        Mt = M_new - theta * K_new * dt
        Kt = (1 - theta) * K_old * dt + M_old
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt
        y_c = convert(SparseMatrixCSC{Float64,Int}, Mt) \ (Kt * y_old[vc] + Ft)
        y_new = vcat(y_c, y_phi)

        # 以步后状态重新计算变量，避免导出时混用步前/步后状态导致的锯齿。
        t_next = t + dt
        _, _, _, variables_new, _ = CallModel(case, y_new, t_next, jacobi="update")
        
        # 提取温度场
        if case.opt.thermal_enabled
            if multi_spme
                T_nodes_carry = get_thermal_dofs(y_new, case.layout)
            else
                nT = mesh_th.nlen
                T_nodes_carry = y_c[(end - nT + 1):end]
            end
            T_max_current = maximum(T_nodes_carry) * T_ref
            T_max_phase = max(T_max_phase, T_max_current)
        end
        
        # 更新电压和容量
        V_current = variables_new["cell voltage"] * phi_scale
        dt_dim = dt * t0_scale
        capacity += abs(I_current) * dt_dim / 3600.0
        t_actual = t_next * t0_scale
        
        # 导出数据（按间隔）
        if step_count % export_interval == 0
            # 提取SOC数据
            soc_n_raw = variables_new["thermal2D element soc_n"]
            soc_p_raw = variables_new["thermal2D element soc_p"]
            soc_n = isa(soc_n_raw, AbstractMatrix) ? copy(vec(soc_n_raw[:, end])) : copy(vec(soc_n_raw))
            soc_p = isa(soc_p_raw, AbstractMatrix) ? copy(vec(soc_p_raw[:, end])) : copy(vec(soc_p_raw))
            
            # 使用步间中点温度导出，抑制 Crank-Nicolson 常见的奇偶步伪振荡。
            T_nodes_out = step_count == 1 ? T_nodes_carry : 0.5 .* (T_nodes_prev_export .+ T_nodes_carry)

            # 转换温度为有量纲值
            T_nodes_K = T_nodes_out .* T_ref
            T_max_K = maximum(T_nodes_K)
            T_elem_K = [mean(@view T_nodes_K[mesh_th.element[e, :]]) for e in 1:ne]
            T_mean_K = dot(T_elem_K, elem_area) / area_sum
            
            # 创建时间步数据
            ts_data = TimeStepData(
                global_time + t_actual,
                phase_type,
                V_current,
                I_current * I_scale,
                T_nodes_K,
                T_max_K,
                T_mean_K,
                soc_n,
                soc_p,
                mean(soc_n)
            )
            push!(timestep_data, ts_data)
        end

        T_nodes_prev_export = copy(T_nodes_carry)
        
        # 检查截止条件
        if phase_type == PHASE_CHARGE && V_current >= V_limit
            terminated_by = :voltage
            break
        elseif phase_type == PHASE_DISCHARGE && V_current <= V_limit
            terminated_by = :voltage
            break
        end
        
        # 自适应时间步长
        if case.opt.dtType == "auto"
            error_y = ErrorEstimation(case, y_old, y_new, dt_min/dt)
            if error_y < 0.5 * case.opt.dtThreshold
                dt = min(dt * 2, dt_max)
            elseif error_y > case.opt.dtThreshold
                dt = max(dt / 2, dt_min)
            end
        end

        # 更新状态
        y_old = y_new
        K_old = K_new
        F_old = F_new
        t = t_next
        
        if t + dt > t_end_nd
            dt = t_end_nd - t
        end
    end
    
    # 构建阶段结果
    result = PhaseResult()
    result.phase_type = phase_type
    result.t_start = global_time
    result.t_end = global_time + t_actual
    result.duration = t_actual
    result.V_start = V_start_actual
    result.V_end = V_current
    result.capacity = capacity
    result.terminated_by = terminated_by
    result.T_max = T_max_phase
    result.T_mean_end = mean(T_nodes_carry) * T_ref
    result.D_max = 0.0
    result.D_mean = 0.0
    result.ΔD_max = 0.0
    
    # 构建最终状态
    result.final_state = Dict{String, Any}(
        "y" => copy(y_old),
        "T_nodes" => copy(T_nodes_carry),
        "V" => V_current,
        "t_global" => global_time + t_actual
    )
    
    return result, timestep_data
end

"""
    solve_cycle_with_export(case, cycle_opt; verbose=true, export_interval=1)

执行单次充放电循环并导出每个时间步的数据。

# 返回
- `CycleExportData`: 包含所有时间步数据的导出对象
"""
function solve_cycle_with_export(case::Case, cycle_opt::CycleOption;verbose::Bool=true, export_interval::Int=1)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    if verbose
        println("="^60)
        println("电化学-热耦合单循环数据导出")
        println("="^60)
        @printf("  导出间隔: 每 %d 个时间步\n", export_interval)
        @printf("  网格: %d 单元, %d 节点\n", ne, nT)
    end
    
    # 应用初始SOC
    soc_init = cycle_opt.SOC_init
    apply_initial_soc!(case, case.param_dim, soc_init)
    if verbose
        @printf("  初始SOC: %.1f%%\n", soc_init * 100)
    end
    
    # 初始状态
    current_state = nothing
    
    all_timestep_data = TimeStepData[]
    
    # ============ 阶段1: 放电 ============
    if verbose
        println("\n[放电阶段]")
    end
    
    discharge_result, discharge_data = solve_phase_with_export(
        case, PHASE_DISCHARGE,
        cycle_opt.t_discharge,
        cycle_opt.I_discharge,
        cycle_opt.V_lower,
        current_state;
        dt_range=cycle_opt.dt_cycle,
        export_interval=export_interval
    )
    current_state = discharge_result.final_state
    append!(all_timestep_data, discharge_data)
    
    if verbose
        @printf("  完成: %.1fs, %d 个数据点, %.3fV→%.3fV\n",
                discharge_result.duration, length(discharge_data),
                discharge_result.V_start, discharge_result.V_end)
    end
    
    # ============ 阶段2: 静置1 ============
    if cycle_opt.t_rest1 > 0
        if verbose
            println("\n[静置1阶段]")
        end
        
        rest1_result, rest1_data = solve_phase_with_export(
            case, PHASE_REST,
            cycle_opt.t_rest1,
            0.0,
            0.0,
            current_state;
            dt_range=cycle_opt.dt_cycle,
            export_interval=export_interval
        )
        current_state = rest1_result.final_state
        append!(all_timestep_data, rest1_data)
        
        if verbose
            @printf("  完成: %.1fs, %d 个数据点\n",
                    rest1_result.duration, length(rest1_data))
        end
    end
    
    # ============ 阶段3: 充电 ============
    if verbose
        println("\n[充电阶段]")
    end
    
    charge_result, charge_data = solve_phase_with_export(
        case, PHASE_CHARGE,
        cycle_opt.t_charge,
        -cycle_opt.I_charge,
        cycle_opt.V_upper,
        current_state;
        dt_range=cycle_opt.dt_cycle,
        export_interval=export_interval
    )
    current_state = charge_result.final_state
    append!(all_timestep_data, charge_data)
    
    if verbose
        @printf("  完成: %.1fs, %d 个数据点, %.3fV→%.3fV\n",
                charge_result.duration, length(charge_data),
                charge_result.V_start, charge_result.V_end)
    end
    
    # ============ 阶段4: 静置2 ============
    if cycle_opt.t_rest2 > 0
        if verbose
            println("\n[静置2阶段]")
        end
        
        rest2_result, rest2_data = solve_phase_with_export(
            case, PHASE_REST,
            cycle_opt.t_rest2,
            0.0,
            0.0,
            current_state;
            dt_range=cycle_opt.dt_cycle,
            export_interval=export_interval
        )
        append!(all_timestep_data, rest2_data)
        
        if verbose
            @printf("  完成: %.1fs, %d 个数据点\n",
                    rest2_result.duration, length(rest2_data))
        end
    end
    
    # 构建导出数据
    export_data = CycleExportData(
        1,
        all_timestep_data,
        mesh_th.node,
        mesh_th.element,
        ne,
        nT
    )
    
    if verbose
        println("\n" * "="^60)
        @printf("导出完成: 共 %d 个时间步数据点\n", length(all_timestep_data))
        println("="^60)
    end
    
    return export_data
end
