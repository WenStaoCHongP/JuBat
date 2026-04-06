function PostProcessing(case::Case, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    result = Dict()
    result["time [s]"]= variables["time"][1:v] * case.param.scale.t0 
    result["cell voltage [V]"]= variables["cell voltage"][1:v] * case.param.scale.phi
    result["cell current [A]"]= variables["cell current"][1:v] * case.param_dim.cell.I1C
    result["temperature [K]"] = variables["temperature"][1:v] * case.param_dim.scale.T_ref
    result["negative particle center radial stress[Pa]"] = variables["negative particle center radial stress"][:,1:v] * case.param.scale.E_n
    result["positive particle center radial stress[Pa]"] = variables["positive particle center radial stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface tangential stress[Pa]"] = variables["negative particle surface tangential stress"][:,1:v] * case.param.scale.E_n
    result["positive particle surface tangential stress[Pa]"] = variables["positive particle surface tangential stress"][:,1:v] * case.param.scale.E_p
    result["negative particle surface displacement[m]"] = variables["negative particle surface displacement"][:,1:v] * case.param.scale.r0
    result["positive particle surface displacement[m]"] = variables["positive particle surface displacement"][:,1:v] * case.param.scale.r0
    if case.opt.model == "SPM"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
    elseif case.opt.model == "SPMe"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        result["negative particle lithium concentration [mol/m^3]"] = variables["negative particle lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle lithium concentration [mol/m^3]"] = variables["positive particle lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative particle surface lithium concentration [mol/m^3]"] = variables["negative particle surface lithium concentration"][:,1:v] * case.param.scale.cn_max
        result["positive particle surface lithium concentration [mol/m^3]"] = variables["positive particle surface lithium concentration"][:,1:v] * case.param.scale.cp_max
        result["negative electrode exchange current density [A/m^2]"] = variables["negative electrode exchange current density"][:,1:v] * case.param.scale.j
        result["positive electrode exchange current density [A/m^2]"] = variables["positive electrode exchange current density"][:,1:v] * case.param.scale.j
        result["negative electrode overpotential [V]"]= variables["negative electrode overpotential"][:,1:v] * case.param.scale.phi
        result["positive electrode overpotential [V]"]= variables["positive electrode overpotential"][:,1:v] * case.param.scale.phi
        result["electrolyte lithium concentration [mol/m^3]"] = variables["electrolyte lithium concentration"][:,1:v] * case.param.scale.ce
        result["negative electrode potential [V]"] = variables["negative electrode potential"][:,1:v] * case.param.scale.phi
        result["positive electrode potential [V]"] = variables["positive electrode potential"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in negative electrode [V]"] = variables["electrolyte potential in negative electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential in positive electrode [V]"] = variables["electrolyte potential in positive electrode"][:,1:v] * case.param.scale.phi
        result["electrolyte potential [V]"] = variables["electrolyte potential"][:,1:v] * case.param.scale.phi
        result["negative electrode open circuit potential [V]"] = variables["negative electrode open circuit potential"] * case.param.scale.phi
        result["positive electrode open circuit potential [V]"] = variables["positive electrode open circuit potential"] * case.param.scale.phi
        result["negative electrode interfacial current density [A/m^2]"]  = variables["negative electrode interfacial current density"] * case.param.scale.j
        result["positive electrode interfacial current density [A/m^2]"]  = variables["positive electrode interfacial current density"] * case.param.scale.j
    end

    if case.opt.thermalmodel == "lumped"
        result["thermal lumped internal heat [W/m^3]"] = vec(variables["thermal lumped internal heat"][1, 1:v]) * case.param.scale.q    
    elseif case.opt.thermalmodel == "distributed2D"
        result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_s_NE [W/m3]"] = variables["thermal2D q_ohm_s_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_e_NE [W/m3]"] = variables["thermal2D q_ohm_e_ne"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_SP [W/m3]"] = variables["thermal2D q_sp"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rxn_PE [W/m3]"] = variables["thermal2D q_rxn_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_rev_PE [W/m3]"] = variables["thermal2D q_rev_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_s_PE [W/m3]"] = variables["thermal2D q_ohm_s_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_ohm_e_PE [W/m3]"] = variables["thermal2D q_ohm_e_pe"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_PCC [W/m3]"] = variables["thermal2D q_pcc"][:, 1:v] * case.param.scale.q
        result["thermal2D Q_NCC [W/m3]"] = variables["thermal2D q_ncc"][:, 1:v] * case.param.scale.q
        result["thermal2D temperature at nodes [K]"] = variables["T_nodes"][:, 1:v] * case.param_dim.scale.T_ref
    end
    """
    if case.opt.czm_enabled == true
        result["czm damage [0-1]"] = variables["czm damage"][:, 1:v]
        result["czm displacement x [m]"] = variables["czm displacement x"][:, 1:v] * case.param.scale.r0
        result["czm displacement y [m]"] = variables["czm displacement y"][:, 1:v] * case.param.scale.r0
        result["czm traction normal [Pa]"] = variables["czm traction normal"][:, 1:v] * case.param.scale.E_n
        result["czm traction tangent [Pa]"] = variables["czm traction tangent"][:, 1:v] * case.param.scale.E_p
        result["czm separation normal [m]"] = variables["czm separation normal"][:, 1:v] * case.param.scale.r0
        result["czm separation tangent [m]"] = variables["czm separation tangent"][:, 1:v] * case.param.scale.r0
    end
    """
    return result
end

function _phase_termination_symbol(phase_type::PhaseType, reason::String)
    if phase_type == PHASE_REST
        return :time
    end
    if reason == "voltage_cutoff_low" || reason == "voltage_cutoff_high" || reason == "all_elements_cutoff"
        return :voltage
    end
    return :time
end

function _state_concentration_variance(case::Case, y_state)
    if y_state === nothing
        return 0.0, 0.0
    end

    y = vec(y_state)
    multi_spme = case.opt.model == "SPMe" && case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D"
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen

    if multi_spme
        ne = case.layout.ne
        n_chem = case.layout.n_chem
        cs_n_all = Float64[]
        cs_p_all = Float64[]
        for e in 1:ne
            offset = (e - 1) * n_chem
            cs_n_e = y[(offset + 1):(offset + Nrn)]
            cs_p_e = y[(offset + Nrn + 1):(offset + Nrn + Nrp)]
            push!(cs_n_all, mean(cs_n_e))
            push!(cs_p_all, mean(cs_p_e))
        end
        return var(cs_n_all), var(cs_p_all)
    end

    cs_n = y[1:Nrn]
    cs_p = y[(Nrn + 1):(Nrn + Nrp)]
    return var(cs_n), var(cs_p)
end

function _postprocess_phase_result(case::Case, phase_type::PhaseType, solve_result::Dict, initial_state::Dict,I_current::Float64, t_start::Float64,D_max_init::Float64, D_mean_init::Float64, czm_mesh)
    time_hist = get(solve_result, "time [s]", Float64[])
    voltage_hist = get(solve_result, "cell voltage [V]", Float64[])
    duration = isempty(time_hist) ? get(case.opt.time, 2, 0.0) : time_hist[end]

    V_start = get(initial_state, "V", isempty(voltage_hist) ? NaN : voltage_hist[1])
    V_end = isempty(voltage_hist) ? V_start : voltage_hist[end]

    T_nodes_K = get(solve_result, "thermal2D T_nodes [K]", Float64[])
    T_hist_K = get(solve_result, "thermal2D temperature [K]", zeros(Float64, 0, 0))
    T_max = if !isempty(T_hist_K)
        maximum(T_hist_K)
    elseif !isempty(T_nodes_K)
        maximum(T_nodes_K)
    else
        case.param_dim.cell.T0
    end
    T_mean_end = !isempty(T_nodes_K) ? mean(T_nodes_K) : case.param_dim.cell.T0

    final_state = get(solve_result, "final_state", Dict{String, Any}())
    final_state["t_global"] = t_start + duration

    D_max_end = D_max_init
    D_mean_end = D_mean_init
    if czm_mesh !== nothing
        stats = get_damage_statistics(czm_mesh)
        D_max_end = stats.max_D
        D_mean_end = stats.mean_D
    end

    terminated_by = _phase_termination_symbol(phase_type, get(solve_result, "termination_reason", "time_limit"))
    capacity = abs(I_current) * duration / 3600.0

    if phase_type == PHASE_REST
        y_start = get(initial_state, "y", nothing)
        y_end = get(final_state, "y", nothing)
        cs_n_init_var, cs_p_init_var = _state_concentration_variance(case, y_start)
        cs_n_final_var, cs_p_final_var = _state_concentration_variance(case, y_end)
        cs_relaxation_n = cs_n_init_var > 1e-12 ? 100.0 * (1.0 - cs_n_final_var / cs_n_init_var) : 0.0
        cs_relaxation_p = cs_p_init_var > 1e-12 ? 100.0 * (1.0 - cs_p_final_var / cs_p_init_var) : 0.0
        final_state["diffusion_active"] = true
        final_state["cs_relaxation_n"] = cs_relaxation_n
        final_state["cs_relaxation_p"] = cs_relaxation_p
    end

    return Dict(
        "duration" => duration,
        "V_start" => V_start,
        "V_end" => V_end,
        "T_max" => T_max,
        "T_mean_end" => T_mean_end,
        "capacity" => capacity,
        "terminated_by" => terminated_by,
        "D_max" => D_max_end,
        "D_mean" => D_mean_end,
        "ΔD_max" => D_max_end - D_max_init,
        "final_state" => final_state
    )
end

function _postprocess_cycle_result!(cycle_result, charge_result, discharge_result, rest1_result, rest2_result, czm_mesh)
    cycle_result.capacity_charge = charge_result.capacity
    cycle_result.capacity_discharge = discharge_result.capacity
    cycle_result.coulombic_efficiency = cycle_result.capacity_charge > 0 ?
        100.0 * cycle_result.capacity_discharge / cycle_result.capacity_charge : 0.0

    if czm_mesh !== nothing
        stats = get_damage_statistics(czm_mesh)
        cycle_result.D_max_end = stats.max_D
        cycle_result.D_mean_end = stats.mean_D
        cycle_result.n_fractured = stats.n_fractured
    end

    cycle_result.T_max = max(charge_result.T_max, rest1_result.T_max, discharge_result.T_max, rest2_result.T_max)
    return nothing
end

function _append_cycle_result!(result, cycle, cycle_result; save_detailed::Bool=false)
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
    return nothing
end

function _update_soh_and_capacity!(result, cycle, cycle_result, initial_capacity::Float64)
    if cycle == 1 && cycle_result.capacity_discharge > 0
        initial_capacity = cycle_result.capacity_discharge
        result.initial_capacity = initial_capacity
    end

    current_soh = initial_capacity > 0 ? cycle_result.capacity_discharge / initial_capacity : 1.0
    push!(result.soh, current_soh)
    return initial_capacity, current_soh
end

function _print_cycle_summary(cycle::Int, cycle_result, current_soh::Float64)
    @printf("  → 循环%d完成: 充%.3fAh, 放%.3fAh, CE=%.1f%%, D_max=%.2f%%, SOH=%.1f%%\n",
            cycle, cycle_result.capacity_charge, cycle_result.capacity_discharge,
            cycle_result.coulombic_efficiency, cycle_result.D_max_end * 100, current_soh * 100)
end

function _check_cycle_termination(cycle::Int, cycle_result, czm_mesh, current_soh::Float64, soh_threshold::Float64; verbose::Bool=true)
    soh_terminated = false

    if current_soh <= soh_threshold && cycle > 1
        if verbose
            @warn "SOH降至$(round(current_soh*100, digits=1))%，低于阈值$(round(soh_threshold*100, digits=1))%，终止循环"
        end
        soh_terminated = true
        return true, soh_terminated
    end

    if czm_mesh !== nothing && cycle_result.n_fractured > 0.5 * czm_mesh.n_cohesive
        if verbose
            @warn "超过50%的内聚力单元断裂，提前终止循环"
        end
        return true, soh_terminated
    end

    return false, soh_terminated
end

function _print_cycling_summary(result, initial_capacity::Float64, soh_terminated::Bool)
    final_soh = initial_capacity > 0 && result.n_cycles > 0 ? result.capacity_discharge[end] / initial_capacity : 1.0
    println("\n" * "="^60)
    println("循环仿真完成")
    println("="^60)
    @printf("  完成循环数: %d\n", result.n_cycles)
    if result.n_cycles > 0
        @printf("  初始容量: %.3fAh\n", initial_capacity)
        @printf("  最终容量: 充%.3fAh, 放%.3fAh\n", result.capacity_charge[end], result.capacity_discharge[end])
        @printf("  容量保持率(SOH): %.1f%%\n", final_soh * 100)
        @printf("  最终损伤: D_max=%.2f%%, D_mean=%.2f%%\n", result.D_max[end] * 100, result.D_mean[end] * 100)
        if soh_terminated
            println("  终止原因: SOH低于阈值")
        end
    end
    return nothing
end

function plot_cycling_results(result; save_path::String="output/")
    isdir(save_path) || mkdir(save_path)

    cycles = result.cycle_idx

    p1 = plot(cycles, result.capacity_discharge,
              xlabel="Cycle Number", ylabel="Discharge Capacity (Ah)",
              label="Discharge", linewidth=2, marker=:circle,
              title="Capacity Fade")
    plot!(p1, cycles, result.capacity_charge,
          label="Charge", linewidth=2, marker=:square, linestyle=:dash)
    savefig(p1, joinpath(save_path, "cycling_capacity.png"))

    p2 = plot(cycles, result.D_max .* 100,
              xlabel="Cycle Number", ylabel="Damage (%)",
              label="D_max", linewidth=2, color=:red,
              title="Damage Evolution")
    plot!(p2, cycles, result.D_mean .* 100,
          label="D_mean", linewidth=2, color=:blue, linestyle=:dash)
    savefig(p2, joinpath(save_path, "cycling_damage.png"))

    p3 = plot(cycles, result.coulombic_efficiency,
              xlabel="Cycle Number", ylabel="Coulombic Efficiency (%)",
              label="CE", linewidth=2, marker=:circle,
              title="Coulombic Efficiency",
              ylims=(95, 105))
    savefig(p3, joinpath(save_path, "cycling_efficiency.png"))

    p4 = plot(cycles, result.T_max,
              xlabel="Cycle Number", ylabel="T_max (K)",
              label="Max Temperature", linewidth=2, marker=:circle,
              title="Temperature History")
    savefig(p4, joinpath(save_path, "cycling_temperature.png"))

    p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900))
    savefig(p_all, joinpath(save_path, "cycling_summary.png"))

    println("✓ 循环结果图已保存至 $save_path")

    return p_all
end
