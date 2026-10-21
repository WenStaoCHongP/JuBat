#=
任务63（review1 返工版）：选定基准边界（fixed_xy_S0_E1_omit_a）多循环真实耦合工况入口。

返工内容（review_round1.md B1/B2 + 计划 5.5）：
- B1/R1 已在 src 修复 endpoint_variant 门控（CzmSolve/Mechanical），本入口显式四字段配置。
- B1/R2 消费生产循环契约：逐圈构造 CycleResult 并调用 update_soh_and_capacity! /
  check_cycle_termination（CyclePostProcess.jl 生产函数，SOH≤阈值第二圈起生效、
  >50% fractured 提前终止）；记录逐圈充/放容量、SOH、Q_abs、EFC 与停止原因。
  入口校验：--case 未知名称报错；默认仅运行 BRIDGE_1，主线工况须显式指定。
- B2/R3 电流单位：I_e 为无量纲契约字段（Σw_e·I_e=I_total），物理等效 I_eq_A=I_e·I_typ，
  权重贡献 I_contrib_A=w_e·I_e·I_typ（5.5 流式记录真实生产消费的 w/A_eff）。
- B2/R4 证据链：mechanical_step_observer 流式接受步（dt/残差/收敛/D），
  失败时保存当步输入与最后提交状态；GP 相位末快照补几何/界面/匝号/历史最大分离/
  线积分权重与末步单元分离/牵引；manifest 补 runner/配置身份、完成/退出状态、
  热重置前后温度统计；输出入 <case>/review1_rework/，旧目录不覆盖。
- 5.5 选定单元 thermal_elem=1：branch_consumption_observer 流式记录每次实际分流消费
  （D_used/A_eff/w/I_e/I_contrib/守恒误差/活动集/单元温度），供损伤—分电流配对图。

用法：julia -t8 --startup-file=no --project=. example/check_selected_boundary_multicycle.jl [--case=BRIDGE_1 ...]
=#
include(joinpath(@__DIR__, "check_mechanical_boundary_branches.jl"))

const SBM_OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "check_selected_boundary_multicycle"))
const SBM_SELECTED_THERMAL_ELEM = 1

const SBM_CASES = [
    (; id="BRIDGE_1",    cycles=1,  area_loss=true,  reset_T_each_cycle=false, reset_T_before_charge=true),
    (; id="BRIDGE_OFF_1",cycles=1,  area_loss=false, reset_T_each_cycle=false, reset_T_before_charge=true),
    (; id="MAIN_ON_5",   cycles=5,  area_loss=true,  reset_T_each_cycle=false, reset_T_before_charge=false),
    (; id="CTRL_OFF_5",  cycles=5,  area_loss=false, reset_T_each_cycle=false, reset_T_before_charge=false),
    (; id="MAIN_ON_12",  cycles=12, area_loss=true,  reset_T_each_cycle=false, reset_T_before_charge=false),
    (; id="CTRL_OFF_12", cycles=12, area_loss=false, reset_T_each_cycle=false, reset_T_before_charge=false),
]

mutable struct SohCarrier
    soh::Vector{Float64}
    initial_capacity::Float64
end

function region_of(r, rmin, rmax)
    f = (r - rmin) / (rmax - rmin)
    f < 1 / 3 ? "inner" : f < 2 / 3 ? "middle" : "outer"
end

function export_selected_element_identity(case, outdir)
    tm = case.mesh["thermal2D"]
    e = SBM_SELECTED_THERMAL_ELEM
    nodes = tm.element[e, :]
    L = case.param.scale.L
    xs = [tm.node[n, 1] for n in nodes] .* L
    ys = [tm.node[n, 2] for n in nodes] .* L
    cm = case.mesh["czm"]
    coh_ids = findall(==(e), [cm.czm_submesh.thermal_elem_map[c.host_outer_elem] for c in cm.cohesive_elements])
    open(joinpath(outdir, "selected_element_identity.csv"), "w") do f
        println(f, "thermal_elem,node1,node2,node3,node4,x_min_m,x_max_m,y_min_m,y_max_m,centroid_r_m")
        println(f, "$e,$(nodes[1]),$(nodes[2]),$(nodes[3]),$(nodes[4]),$(minimum(xs)),$(maximum(xs)),$(minimum(ys)),$(maximum(ys)),$(hypot(mean(xs), mean(ys)))")
        println(f, "cohesive_map,elem,iface,material_layer,winding_turn,x_m,y_m,r_m,length_m")
        for c in coh_ids
            ce = cm.cohesive_elements[c]
            ns = vcat(ce.nodes_bottom, ce.nodes_top)
            x = sum(cm.node[ns, 1]) / length(ns) * L
            y = sum(cm.node[ns, 2]) / length(ns) * L
            println(f, "cohesive_map,$c,$(ce.interface_type),$(cm.czm_submesh.material_type[ce.host_outer_elem]),$(cm.czm_submesh.winding_turn[ce.host_outer_elem]),$x,$y,$(hypot(x, y)),$(ce.length)")
        end
    end
    return length(coh_ids)
end

function export_gp_phaseend(case, outdir, cycle, label, last_result)
    ms = case.mech
    G = ms.gp_damage_states
    cm = case.mesh["czm"]
    L, δ_czm, σ_czm = case.param.scale.L, case.param.scale.δ_czm, case.param.scale.σ_czm
    rmin, rmax = extrema(hypot.(cm.node[:,1], cm.node[:,2])) .* L
    open(joinpath(outdir, "gp_snapshots_phaseend.csv"), "a") do f
        for k in eachindex(cm.cohesive_elements)
            ce = cm.cohesive_elements[k]
            ns = vcat(ce.nodes_bottom, ce.nodes_top)
            x = sum(cm.node[ns, 1]) / length(ns) * L
            y = sum(cm.node[ns, 2]) / length(ns) * L
            θ = atan(y, x) < 0 ? atan(y, x) + 2π : atan(y, x)
            turn = cm.czm_submesh.winding_turn[ce.host_outer_elem]
            sep_n = last_result === nothing ? "" : string(last_result.separation_n[k] * δ_czm)
            sep_t = last_result === nothing ? "" : string(last_result.separation_t[k] * δ_czm)
            tr_n = last_result === nothing ? "" : string(last_result.traction_n[k] * σ_czm)
            tr_t = last_result === nothing ? "" : string(last_result.traction_t[k] * σ_czm)
            for g in 1:size(G,2)
                s = G[k,g]
                println(f, "$cycle,$label,$k,$g,$(ce.interface_type),$turn,$x,$y,$(hypot(x,y)),$θ,$(region_of(hypot(x,y), rmin, rmax))," *
                        "$(s.D),$(s.D_visc),$(s.fractured),$(s.δ_max_n*δ_czm),$(s.δ_max_t*δ_czm),$(s.δ_max_eff*δ_czm)," *
                        "$sep_n,$sep_t,$tr_n,$tr_t,$(ce.length)")
            end
        end
    end
    return nothing
end

function run_case(spec, outdir)
    mkpath(outdir)
    t_all = time()
    case = build_collection_case(outdir)
    czm = case.opt.czm
    czm.fix_inner, czm.fix_start, czm.fix_end = false, false, true
    czm.endpoint_variant = :omit_a
    czm.area_loss_enabled = spec.area_loss
    boundary = JuBat.resolve_mechanical_bc(case.mesh["czm"], case.param, czm; opt=case.opt)
    length(boundary.nodes) == 92 || error("boundary node count $(length(boundary.nodes)) != 92 (nθ=80 baseline)")
    haskey(boundary.nodes, boundary.a) && error("baseline requires a free")
    write_csv(joinpath(outdir, "constraint_nodes.csv"), ("node","mode","sources"),
        ((n,m,join(get(boundary.provenance,n,Symbol[]),"|")) for (n,m) in sort(collect(boundary.nodes))))
    sources = source_identity()
    write_csv(joinpath(outdir, "source_identity.csv"), ("source_relative_path","sha256"),
        (r for r in vcat(sources, [(path="example/check_selected_boundary_multicycle.jl",
            sha256=bytes2hex(sha256(read(@__FILE__))))])))
    n_coh_map = export_selected_element_identity(case, outdir)
    I_typ = case.param.scale.I_typ
    t0 = case.param.scale.t0

    soh_carrier = SohCarrier(Float64[], 0.0)
    manifest_extra = Dict{String,Any}()
    # --- 双观察点：机械接受步流式 + 选定单元分流消费流式 ---
    mech_io = open(joinpath(outdir, "mech_step_stream.csv"), "w")
    println(mech_io, "commit_idx,cycle,phase,t_phase_s,t_cycle_s,dt_s,residual_norm,converged,D_max_elem,D_max_gp")
    branch_io = open(joinpath(outdir, "selected_element_history.csv"), "w")
    println(branch_io, "case,cycle,phase,step,t_phase_s,t_cycle_s,thermal_elem,last_mech_commit,area_loss_enabled," *
        "has_damage_input,D_state_at_eval,D_used,D_committed_after,A_geom,A_eff_used,w_used,I_e,I_eq_A,I_contrib_A," *
        "I_total_A,sum_I_contrib_A,conservation_abs_err,active_mask,inactive_reason,T_elem_K")
    commit_idx = 0
    commit_times = Float64[]
    commit_D = Float64[]
    pending_input = Ref{Any}(nothing)
    last_result = Ref{Any}(nothing)
    failure_dumped = Ref(false)
    cur_cycle, cur_phase, phase_offset, phase_elapsed = Ref(0), Ref(""), Ref(0.0), Ref(0.0)
    branch_count = 0
    prev_branch_row = Ref{Any}(nothing)

    mech_observer = function(event)
        if event.phase == :before
            pending_input[] = event
        elseif event.phase == :after
            last_result[] = event.result
            if event.result.converged
                commit_idx += 1
                dt = event.dt_seconds === nothing ? NaN : Float64(event.dt_seconds)
                phase_elapsed[] += dt
                t_phase = phase_elapsed[]
                t_cycle = phase_offset[] + t_phase
                push!(commit_times, t_cycle)
                Gp = event.ms.gp_damage_states
                dgp = maximum(Gp[e,g].D for e in 1:size(Gp,1), g in 1:size(Gp,2))
                de = maximum(s.D for s in event.ms.damage_states)
                println(mech_io, "$commit_idx,$(cur_cycle[]),$(cur_phase[]),$t_phase,$t_cycle,$dt," *
                        "$(event.result.residual_norm),true,$de,$dgp")
                push!(commit_D, de)
                flush(mech_io)
                pending_input[] = nothing
            else
                println(mech_io, "NA,$(cur_cycle[]),$(cur_phase[]),$(phase_elapsed[]),$(phase_offset[]+phase_elapsed[]),NA," *
                        "$(event.result.residual_norm),false,NA,NA")
            end
        elseif event.phase == :error
            if !failure_dumped[]
                failure_dumped[] = true
                fdir = joinpath(outdir, "failure")
                mkpath(fdir)
                open(joinpath(fdir, "failure.txt"), "w") do f
                    println(f, sprint(showerror, event.error))
                    println(f, "cycle=$(cur_cycle[]) phase=$(cur_phase[]) t_phase_s=$(phase_elapsed[]) commit_idx=$commit_idx")
                end
                pending_input[] === nothing || save_binary(joinpath(fdir, "failed_input.bin"),
                    (; F_ext=copy(pending_input[].F_ext), dt_seconds=pending_input[].dt_seconds,
                       dT_elem=copy(pending_input[].dT_elem), Δsoc_n_elem=copy(pending_input[].Δsoc_n_elem),
                       Δsoc_p_elem=copy(pending_input[].Δsoc_p_elem), ms_before=deepcopy(pending_input[].ms)))
                save_binary(joinpath(fdir, "last_committed_state.bin"), deepcopy(case.mech))
            end
        end
        return nothing
    end

    branch_observer = function(ev)
        branch_count += 1
        te = ev.t * t0
        te_global = phase_offset[] + te
        last_commit = isempty(commit_times) ? 0 : searchsortedlast(commit_times, te_global + 1e-9)
        d_state = ev.D_elem === nothing ?
            JuBat.map_czm_damage_to_thermal(case.mesh["czm"], case.mech.damage_states,
                length(ev.areas))[SBM_SELECTED_THERMAL_ELEM] : ev.D_elem[SBM_SELECTED_THERMAL_ELEM]
        d_used = ev.D_elem === nothing ? "" : string(ev.D_elem[SBM_SELECTED_THERMAL_ELEM])
        ie = ev.I_e[SBM_SELECTED_THERMAL_ELEM]
        w1 = ev.w[SBM_SELECTED_THERMAL_ELEM]
        contrib = w1 * ie * I_typ
        total = ev.I_total * I_typ
        sumc = sum(ev.w[e] * ev.I_e[e] * I_typ for e in eachindex(ev.I_e))
        row = (cycle=cur_cycle[], phase=cur_phase[], step=branch_count, t_phase=te,
            t_cycle=phase_offset[] + te, elem=SBM_SELECTED_THERMAL_ELEM, last_commit=last_commit,
            has_damage=ev.D_elem !== nothing, d_state=d_state, d_used=d_used,
            a_geom=ev.areas[SBM_SELECTED_THERMAL_ELEM], a_eff=ev.A_eff[SBM_SELECTED_THERMAL_ELEM],
            w=w1, ie=ie, ityp=I_typ, contrib=contrib, total=total, sumc=sumc,
            cons_err=abs(sumc - total), active=ev.active_mask[SBM_SELECTED_THERMAL_ELEM],
            inactive_reason=ev.inactive_reason[SBM_SELECTED_THERMAL_ELEM],
            T_elem=ev.Te_prev[SBM_SELECTED_THERMAL_ELEM])
        prev = prev_branch_row[]
        if prev !== nothing  # D_committed_after = 该消费之后下一次消费所见的 D_state
            write_branch_row(branch_io, spec.id, prev, string(row.d_state))
        end
        prev_branch_row[] = row
        return nothing
    end

    function flush_last_branch_row()
        prev = prev_branch_row[]
        prev === nothing && return nothing
        write_branch_row(branch_io, spec.id, prev, "")
        prev_branch_row[] = nothing
        return nothing
    end

    # --- 逐相位运行（含热重置记录、SOH/终止契约） ---
    steps_io = open(joinpath(outdir, "steps.csv"), "w")
    println(steps_io, "cycle,phase,step,time_s,voltage_V,D_max")
    phase_io = open(joinpath(outdir, "phase_summary.csv"), "w")
    println(phase_io, "cycle,phase,t_start_s,t_end_s,duration_s,terminated_by,V_start,V_end,T_max_K,D_max,D_mean,wall_s,reset_T_before_phase,T_mean_pre_reset_K,T_mean_post_reset_K")
    open(joinpath(outdir, "phase_end_currents.csv"), "w") do f
        println(f, "cycle,phase,elem,I_e,I_eq_A") end
    open(joinpath(outdir, "gp_snapshots_phaseend.csv"), "w") do f
        println(f, "cycle,phase,elem,gp,iface,winding_turn,x_m,y_m,r_m,theta_rad,region,D,D_visc,fractured," *
            "delta_max_n_m,delta_max_t_m,delta_max_eff_m,sep_n_m,sep_t_m,traction_n_Pa,traction_t_Pa,length_weight_m") end
    open(joinpath(outdir, "regional_damage_phaseend.csv"), "w") do f
        println(f, "cycle,phase,iface,region,elements,D_max,D_mean,damaged_elements") end
    open(joinpath(outdir, "cycle_summary.csv"), "w") do f
        println(f, "case,cycle,capacity_charge_Ah,capacity_discharge_Ah,coulombic_efficiency,D_max_end,n_fractured,T_max,SOH,Q_abs_Ah,EFC,stop_reason") end

    cycle_opt = diagnostic_cycle_options()
    phase_specs = (("discharge",JuBat.PHASE_DISCHARGE,cycle_opt.t_discharge,cycle_opt.I_discharge),
        ("rest1",JuBat.PHASE_REST,cycle_opt.t_rest1,0.0),
        ("charge",JuBat.PHASE_CHARGE,cycle_opt.t_charge,-cycle_opt.I_charge),
        ("rest2",JuBat.PHASE_REST,cycle_opt.t_rest2,0.0))
    state = nothing
    completed_cycles, stop_reason, error_text, last_phase_done = 0, "all_cycles_complete", "", ""
    Q_abs_Ah = 0.0
    old_mech_obs = JuBat.mechanical_step_observer[]
    old_branch_obs = JuBat.branch_consumption_observer[]
    JuBat.mechanical_step_observer[] = mech_observer
    JuBat.branch_consumption_observer[] = branch_observer
    try
        for cycle in 1:spec.cycles
            cycle_phases = Dict{String,Any}()
            for (label,phase,dur,cur) in phase_specs
                cur_cycle[], cur_phase[] = cycle, label
                # t_cycle 修正（2026-10-10）：全局偏移取上一相位末电化学状态（无状态=首相位 0），
                # 相位内累计复位；观察点写出的 t_phase_s/t_cycle_s 由此恢复各自语义。
                phase_offset[] = state === nothing ? 0.0 : Float64(state["t_global"])
                phase_elapsed[] = 0.0
                do_reset = label == "charge" && spec.reset_T_before_charge && state !== nothing
                T_pre = T_post = ""
                if do_reset
                    T_pre = string(mean(state["T_nodes"]) * case.param.scale.T_ref)
                    JuBat.reset_cycle_temperature!(case, state)
                    T_post = string(mean(state["T_nodes"]) * case.param.scale.T_ref)
                    println("[$(spec.id)] cycle=$cycle temperature reset before charge: mean T $T_pre -> $T_post K")
                end
                t0p = time()
                pr = JuBat.solve_phase(case, phase, dur, cur, state;
                    ms=case.mech, dt_range=cycle_opt.dt_cycle, czm_cycle=cycle)
                wall = time() - t0p
                state = pr.final_state
                cycle_phases[label] = pr
                last_phase_done = label
                Q_abs_Ah += abs(cur) * pr.duration / 3600.0
                sr = pr.solve_result
                tvec, V, D = sr["time [s]"], sr["cell voltage [V]"], sr["czm D_max"]
                for k in eachindex(tvec)
                    println(steps_io, "$cycle,$label,$k,$(pr.t_start + tvec[k] - tvec[1]),$(V[k]),$(D[k])")
                end
                Ie = sr["thermal2D element current"][:, end]
                open(joinpath(outdir, "phase_end_currents.csv"), "a") do f
                    for e in eachindex(Ie)
                        println(f, "$cycle,$label,$e,$(Ie[e]),$(Ie[e]*I_typ)")
                    end
                end
                Dm = maximum(s.D for s in case.mech.damage_states)
                Dmean = mean(s.D for s in case.mech.damage_states)
                println(phase_io, "$cycle,$label,$(pr.t_start),$(pr.t_end),$(pr.duration),$(pr.terminated_by)," *
                        "$(pr.V_start),$(pr.V_end),$(pr.T_max),$Dm,$Dmean,$(round(wall;digits=1)),$do_reset,$T_pre,$T_post")
                flush_last_branch_row()
                export_gp_phaseend(case, outdir, cycle, label, last_result[])
                # 区域损伤（半径三分）
                cm = case.mesh["czm"]
                rs = [hypot(cm.node[c.nodes_bottom[1],1], cm.node[c.nodes_bottom[1],2]) * case.param.scale.L
                      for c in cm.cohesive_elements]
                rmin, rmax = extrema(hypot.(cm.node[:,1], cm.node[:,2])) .* case.param.scale.L
                buckets = Dict{Tuple{Symbol,String},Vector{Float64}}()
                for (k,c) in enumerate(cm.cohesive_elements)
                    push!(get!(buckets, (c.interface_type, region_of(rs[k], rmin, rmax)), Float64[]),
                        case.mech.damage_states[k].D)
                end
                open(joinpath(outdir, "regional_damage_phaseend.csv"), "a") do f
                    for key in sort(collect(keys(buckets)); by = k -> (string(k[1]), k[2]))
                        v = buckets[key]
                        println(f, "$cycle,$label,$(key[1]),$(key[2]),$(length(v)),$(maximum(v)),$(sum(v)/length(v)),$(count(>(1e-8), v))")
                    end
                end
                println("[$(spec.id)] CYCLE=$cycle PHASE=$label t_end=$(pr.t_end) term=$(pr.terminated_by) D_max=$Dm wall=$(round(wall;digits=1))s")
            end
            # B1/R2：生产循环契约——SOH/容量更新与终止检查（复用 CyclePostProcess 生产函数）
            cycle_result = JuBat.CycleResult(cycle, cycle_phases["charge"], cycle_phases["rest1"],
                cycle_phases["discharge"], cycle_phases["rest2"],
                cycle_phases["charge"].capacity, cycle_phases["discharge"].capacity,
                cycle_phases["discharge"].capacity / max(cycle_phases["charge"].capacity, eps()) * 100.0,
                maximum(s.D for s in case.mech.damage_states),
                mean(s.D for s in case.mech.damage_states),
                count(s.fractured for s in case.mech.damage_states),
                maximum(cycle_phases[p].T_max for p in ("discharge","rest1","charge","rest2")))
            initial_capacity, current_soh = JuBat.update_soh_and_capacity!(soh_carrier, cycle, cycle_result, soh_carrier.initial_capacity)
            efc = Q_abs_Ah / (2.0 * 5.0)
            should_stop, soh_hit = JuBat.check_cycle_termination(cycle, cycle_result, case.mech,
                current_soh, case.opt.czm.soh_threshold)
            stop_label = should_stop ? (soh_hit ? "soh_threshold" : "fractured_majority") : "continue"
            open(joinpath(outdir, "cycle_summary.csv"), "a") do f
                println(f, "$(spec.id),$cycle,$(cycle_result.capacity_charge),$(cycle_result.capacity_discharge)," *
                        "$(cycle_result.coulombic_efficiency),$(cycle_result.D_max_end),$(cycle_result.n_fractured)," *
                        "$(cycle_result.T_max),$current_soh,$Q_abs_Ah,$efc,$stop_label")
            end
            completed_cycles = cycle
            if should_stop
                stop_reason = stop_label
                break
            end
        end
    catch err
        stop_reason = "failed"
        error_text = sprint(showerror, err)
        println("[$(spec.id)] FAILED after cycle $completed_cycles (last phase done: $last_phase_done): $error_text")
    finally
        JuBat.mechanical_step_observer[] = old_mech_obs
        JuBat.branch_consumption_observer[] = old_branch_obs
        flush_last_branch_row()
        close(mech_io); close(branch_io); close(steps_io); close(phase_io)
        manifest_extra = Dict("completed_cycles"=>completed_cycles, "stop_reason"=>stop_reason,
            "error"=>error_text, "last_phase_done"=>last_phase_done,
            "mech_commits"=>commit_idx, "branch_consumptions"=>branch_count,
            "total_wall_s"=>round(time() - t_all; digits=1),
            "Q_abs_Ah"=>Q_abs_Ah, "EFC"=>Q_abs_Ah / 10.0)
    end
    # 选定单元身份核对（5.5）：损伤映射位置须仍在此热单元
    D_elem_check = JuBat.map_czm_damage_to_thermal(case.mesh["czm"], case.mech.damage_states,
        size(case.mesh["thermal2D"].element,1))
    open(joinpath(outdir, "run_manifest.toml"), "w") do io
        TOML.print(io, merge(Dict("schema"=>"JuBat selected-boundary multicycle v1 (review1 rework)",
            "case"=>spec.id, "created_at"=>string(now()), "cycles_requested"=>spec.cycles,
            "area_loss_enabled"=>spec.area_loss, "reset_T_each_cycle"=>spec.reset_T_each_cycle,
            "reset_T_before_charge"=>spec.reset_T_before_charge,
            "boundary"=>Dict("outer_bc"=>"fixed_xy","fix_inner"=>false,"fix_start"=>false,
                "fix_end"=>true,"endpoint_variant"=>"omit_a"),
            "constrained_nodes"=>length(boundary.nodes), "outer_count"=>boundary.outer_count,
            "a_node"=>boundary.a, "a_fixed"=>haskey(boundary.nodes,boundary.a),
            "b_node"=>boundary.b, "b_fixed"=>haskey(boundary.nodes,boundary.b),
            "ntheta"=>80, "thermal_elements"=>size(case.mesh["thermal2D"].element,1),
            "cohesive_elements"=>length(case.mesh["czm"].cohesive_elements),
            "bulk_elements"=>size(case.mesh["czm"].bulk_element,1),
            "source_hash"=>bytes2hex(sha256(string(sources))),
            "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
            "config_hash"=>serialized_hash(options_identity(czm)),
            "review1_rework"=>true,
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
            "selected_thermal_elem"=>SBM_SELECTED_THERMAL_ELEM,
            "selected_elem_cohesive_count"=>n_coh_map,
            "selected_elem_D_final"=>D_elem_check[SBM_SELECTED_THERMAL_ELEM],
            "I_typ_A"=>I_typ, "current_units"=>"I_e dimensionless (sum(w*I_e)=I_total); I_eq_A=I_e*I_typ; I_contrib_A=w*I_eq_A"),
            manifest_extra))
    end
    open(joinpath(outdir, "report.md"), "w") do f
        println(f, "# $spec.id — review1 rework run\n")
        println(f, "- Boundary: fixed_xy / S0 / E1 / omit_a; constrained nodes $(length(boundary.nodes)) (a free, b via end).")
        println(f, "- Cycles: $completed_cycles/$(spec.cycles); stop=$stop_reason; wall $(round(time()-t_all;digits=1)) s; "
                  * "mech commits=$commit_idx, branch consumptions=$branch_count.")
        println(f, "- area_loss=$(spec.area_loss), reset_T_before_charge=$(spec.reset_T_before_charge); nθ=80; gp_basic+τ=5 s.")
        println(f, "- Selected thermal elem $SBM_SELECTED_THERMAL_ELEM: $n_coh_map cohesive elements mapped; final D_elem=$(D_elem_check[SBM_SELECTED_THERMAL_ELEM]).")
        if error_text != ""
            println(f, "- Error: ```$(error_text)```")
        end
        println(f, "\n证据：run_manifest / constraint_nodes / source_identity(+runner) / selected_element_identity / mech_step_stream / "
                  * "selected_element_history / steps / phase_summary / phase_end_currents(I_e,I_eq_A) / gp_snapshots_phaseend / "
                  * "regional_damage / cycle_summary(容量/SOH/Q_abs/EFC) / failure/（如失败）。")
    end
    println("[$(spec.id)] DONE cycles=$completed_cycles/$(spec.cycles) stop=$stop_reason wall=$(round(time()-t_all;digits=1))s")
    return stop_reason == "all_cycles_complete"
end

function write_branch_row(io, case_id, r, d_after)
    println(io, "$case_id,$(r.cycle),$(r.phase),$(r.step),$(r.t_phase),$(r.t_cycle),$(r.elem),$(r.last_commit)," *
            "$(r.has_damage ? "true" : "false"),$(r.has_damage ? "true" : "false"),$(r.d_state),$(r.d_used),$d_after," *
            "$(r.a_geom),$(r.a_eff),$(r.w),$(r.ie),$(r.ie*r.ityp),$(r.contrib),$(r.total),$(r.sumc)," *
            "$(r.cons_err),$(r.active),$(r.inactive_reason),$(r.T_elem)")
    return nothing
end

function main(args=ARGS)
    wanted = String[]
    out_tag = "review1_rework"
    for arg in args
        startswith(arg,"--case=") ? push!(wanted,split(arg,"=";limit=2)[2]) :
        startswith(arg,"--out=") ? (out_tag = split(arg,"=";limit=2)[2]) :
        throw(ArgumentError("unknown argument: $arg"))
    end
    isempty(wanted) && (wanted = ["BRIDGE_1"])
    known = [c.id for c in SBM_CASES]
    unknown = setdiff(wanted, known)
    isempty(unknown) || throw(ArgumentError("unknown case(s) $(join(unknown, ",")); valid: $(join(known, ", ")))"))
    all_ok = true
    for spec in SBM_CASES
        spec.id in wanted || continue
        all_ok &= run_case(spec, joinpath(SBM_OUTPUT, spec.id, out_tag))
    end
    return all_ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() || exit(1)
end
