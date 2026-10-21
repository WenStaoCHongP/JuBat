#=
任务63 R0-b（rest1_root_cause.md §6 / task_plan 5.7.2）：静置1根因诊断 v4。

协议与 BRIDGE_1 完全一致（基准边界 omit_a、gp_basic+τ5、面积 ON、充电前温度重置、nθ=80、gsorder=2）。

v4（5.7.2 返工）要点：
- ① CSV 列错配修复：averages 12 列满行（补 PE 末值），write_csv_checked 硬门行宽==列数；
  BASE 故障分支与正常分支同 schema（不再重复写窄表头）。
- ② 真实 GP 逐点重建：committed_gp_separations 按生产算法（czm.jl 装配：B_global→R 标架→Λ）
  分别重建 GP1/GP2 的即时 δ_n/δ_t（NC order=2 两端点各对应一对节点），不做端点平均。
  即时值/历史包络（δ_max_*）/单元 GP 加权均值（*_elemavg_*）三种口径各自命名。
- ④ 诊断补齐：radial_profiles.csv（NE/PE elem1 径向 c/r/FE 权重，相位末）；
  averages_stepwise.csv（全场逐步两种平均 min/mean/max）。
- ⑤ 状态闭环：try/catch 回写 failed 状态+相位+最后接受步；BC 违约 = 约束节点位移范数。
- ⑥ 报告引用真实文件，列出共同起点/实际末/冻结末/增量/差值与分离比较。

输出 output/rest1_root_diag/fe_mean_rework_r2/<attempt>/。旧 fe_mean_rework_r1 的
错列表已标 DIAGNOSTIC_INVALID.md，保留为历史证据。
Julia 只导数据；绘图/深度后处理按用户架构约定由 Python 在 Desktop/output/postprocess 完成。
=#
include(joinpath(@__DIR__, "check_mechanical_boundary_branches.jl"))

const DIAG_ROOT = normpath(joinpath(@__DIR__, "..", "output", "rest1_root_diag", "fe_mean_rework_r4"))
const DIAG_OUT = joinpath(DIAG_ROOT, "attempt_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
const DIAG_GP_ELEM = 3365   # findings 已核实的唯一受损记录 elem3365/gp1（NE_NCC 内侧 θ=0）

"""写 CSV 硬门：每行元组长度必须等于表头列数（5.7.2-①）。"""
function write_csv_checked(path, header, rows)
    n = length(header)
    for (k, r) in enumerate(rows)
        length(r) == n || error("write_csv_checked $path row $k: $(length(r)) fields vs $n columns")
    end
    write_csv(path, header, rows)
    # 写后读回核验列数
    open(path) do f
        head = split(readline(f), ',')
        length(head) == n || error("write_csv_checked $path: header readback $(length(head)) != $n")
        k = 0
        for line in eachline(f)
            k += 1
            # 简单计数：字段中不含引号包裹的逗号（本脚本全部数值/短标签）
            length(split(line, ',')) == n || error("write_csv_checked $path readback row $k width mismatch")
        end
    end
    return nothing
end

"""② 真实 GP 逐点分离重建（生产算法只读复刻，不推进损伤）：
节点序 n1,n2=bottom, n4,n3=top；B_global 行/列与 czm.jl 装配一致；
δ_local = R·B_global·u_e；δ_n = Λ·δ_local[1]、δ_t = Λ·δ_local[2]（δ_czm 归一）。
返回每个 GP 的 (δ_n, δ_t)（无量纲），物理米 = ×scale.δ_czm。"""
function committed_gp_separations(czm_mesh, u, elem_id, param)
    elem = czm_mesh.cohesive_elements[elem_id]
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top          # 生产节点序：dofs=[n1,n2,n3,n4]
    _, _, _, R = JuBat.cohesive_local_frame(czm_mesh, elem)
    order = czm_mesh.czm_submesh.mesh_bonded.gs.order
    _, pts = JuBat.NCweight(order)   # 5.7.3-①：NCweight 返回 (w, q)，旧版误接为 (pts, w)
    order == 2 && (pts == [-1.0, 1.0] || error("committed_gp_separations: NC order2 pts != [-1,1]"))
    Λ = param.scale.L / param.scale.δ_czm
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
    u_e = [u[d] for d in dofs]
    out = Vector{Tuple{Float64,Float64}}(undef, length(pts))
    for (g, ξ) in enumerate(pts)
        N1 = 0.5 * (1.0 - ξ)
        N2 = 0.5 * (1.0 + ξ)
        B = zeros(2, 8)
        B[1,1] = -N1; B[2,2] = -N1; B[1,3] = -N2; B[2,4] = -N2
        B[1,5] =  N2; B[2,6] =  N2; B[1,7] =  N1; B[2,8] =  N1
        δ_local = R * (B * u_e)
        out[g] = (Λ * δ_local[1], Λ * δ_local[2])
    end
    return out
end

"""⑤ BC 违约：约束节点（fixed_xy）的位移范数（同次约束值恒为 0，违约即非零）。"""
function bc_violation(u, boundary)
    v = 0.0
    for (n, mode) in boundary.nodes
        mode === :fixed_xy && (v = max(v, abs(u[2*n-1]), abs(u[2*n])))
    end
    return v
end

gp_key_quantities(ms) = (D = ms.gp_damage_states[DIAG_GP_ELEM,1].D,
    D_visc = ms.gp_damage_states[DIAG_GP_ELEM,1].D_visc,
    δ_max_n = ms.gp_damage_states[DIAG_GP_ELEM,1].δ_max_n,
    δ_max_t = ms.gp_damage_states[DIAG_GP_ELEM,1].δ_max_t,
    δ_max_eff = ms.gp_damage_states[DIAG_GP_ELEM,1].δ_max_eff,
    D_elem_mean = ms.damage_states[DIAG_GP_ELEM].D)

function write_manifest(status; extra = Dict{String,Any}())
    d = Dict("schema" => "rest1_root_diag v5 (fe_mean_rework_r4)",
        "status" => status, "created_at" => string(now()),
        "runner_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "attempt_dir" => DIAG_OUT)
    merge!(d, extra)
    open(joinpath(DIAG_OUT, "capture_manifest.toml"), "w") do io
        TOML.print(io, d)
    end
end

function main()
    mkpath(DIAG_OUT)
    t_all = time()
    write_manifest("running")
    write_csv(joinpath(DIAG_OUT, "source_identity.csv"), ("source_relative_path","sha256"), source_identity())
    last_stage = Ref("capture")       # 5.7.3-④：capture/BASE/HOLD/report 分阶段跟踪
    last_phase = Ref("")
    last_accepted_step = Ref(0)
    try
        run_diag(t_all, last_stage, last_phase, last_accepted_step)
    catch err
        write_manifest("failed"; extra = Dict("failed_stage" => last_stage[], "failed_phase" => last_phase[],
            "last_accepted_step" => last_accepted_step[], "error" => sprint(showerror, err)))
        rethrow()
    end
    return true
end

function run_diag(t_all, last_stage, last_phase, last_accepted_step)
    case = build_collection_case(DIAG_OUT)
    czm = case.opt.czm
    czm.fix_inner, czm.fix_start, czm.fix_end = false, false, true
    czm.endpoint_variant = :omit_a
    czm.area_loss_enabled = true
    boundary = JuBat.resolve_mechanical_bc(case.mesh["czm"], case.param, czm; opt=case.opt)
    cm, param = case.mesh["czm"], case.param
    δ_czm, σ_czm = param.scale.δ_czm, param.scale.σ_czm

    steps = Any[]
    phase_end_ms = Dict{String,Any}()
    phase_end_state = Dict{String,Any}()
    phase_ref = Ref("")
    # 5.7.4：化学步观察点事件采集
    chem_events = Any[]
    chem_observer_enabled = Ref(true)
    capture = function(event)
        if event.phase == :after && event.result.converged
            r = event.result
            gp = committed_gp_separations(cm, event.ms.u_prev, DIAG_GP_ELEM, param)
            push!(steps, (phase = phase_ref[], dt = event.dt_seconds === nothing ? 0.0 : Float64(event.dt_seconds),
                F_ext = copy(event.F_ext), dT = copy(event.dT_elem),
                dsn = copy(event.Δsoc_n_elem), dsp = copy(event.Δsoc_p_elem),
                after = gp_key_quantities(event.ms),
                resid = r.residual_norm, iters = r.iterations,
                bcviol = bc_violation(event.ms.u_prev, boundary),
                sep_n_elemavg = r.separation_n[DIAG_GP_ELEM],
                sep_t_elemavg = r.separation_t[DIAG_GP_ELEM],
                trac_n_elemavg = r.traction_n[DIAG_GP_ELEM],
                trac_t_elemavg = r.traction_t[DIAG_GP_ELEM],
                g1_dn = gp[1][1], g1_dt = gp[1][2],
                g2_dn = gp[2][1], g2_dt = gp[2][2]))
            phase_end_ms[phase_ref[]] = deepcopy(event.ms)
            last_accepted_step[] += 1
        end
        return nothing
    end
    old_obs = JuBat.mechanical_step_observer[]
    JuBat.mechanical_step_observer[] = capture
    cycle_opt = diagnostic_cycle_options()
    specs = (("discharge",JuBat.PHASE_DISCHARGE,cycle_opt.t_discharge,cycle_opt.I_discharge),
        ("rest1",JuBat.PHASE_REST,cycle_opt.t_rest1,0.0),
        ("charge",JuBat.PHASE_CHARGE,cycle_opt.t_charge,-cycle_opt.I_charge),
        ("rest2",JuBat.PHASE_REST,cycle_opt.t_rest2,0.0))
    state = nothing
    soc_avgs = Any[]      # ① 12 列满行
    radial_rows = Any[]   # ③ 径向节点表
    gp_rows = Any[]       # ③ 完整 GP 表
    stepwise_avgs = Any[] # ④ 全场逐步两种平均
    try
        for (label,phase,dur,cur) in specs
            phase_ref[] = label
            last_phase[] = label
            label == "charge" && JuBat.reset_cycle_temperature!(case, state)
            # 5.7.4：化学步观察点透传（库存预算采集，正式回调非空回调）
            chem_obs = if chem_observer_enabled[]
                ev -> push!(chem_events, (phase = label, kind = ev.kind,
                    h_used = get(ev, :h_used, nothing), h = get(ev, :h, nothing),
                    dt_next = get(ev, :dt_next, nothing),
                    trial_idx = get(ev, :trial_idx, nothing),
                    accepted_idx = get(ev, :accepted_idx, nothing)))
            else
                nothing
            end
            pr = JuBat.solve_phase(case, phase, dur, cur, state;
                ms=case.mech, dt_range=cycle_opt.dt_cycle, chemical_step_observer=chem_obs)
            state = pr.final_state
            phase_end_state[label] = deepcopy(state)
            sr = pr.solve_result
            if haskey(sr, "thermal2D element soc_n_vol")
                sn, snv = sr["thermal2D element soc_n"], sr["thermal2D element soc_n_vol"]
                sp, spv = sr["thermal2D element soc_p"], sr["thermal2D element soc_p_vol"]
                ne = size(sn, 1)
                push!(soc_avgs, (label, pr.t_end,
                    sn[1,1], snv[1,1], sp[1,1], spv[1,1],
                    sn[1,end], snv[1,end], sp[1,end], spv[1,end],
                    maximum(sn) - minimum(sn), maximum(snv) - minimum(snv)))
                for k in axes(sn, 2)
                    push!(stepwise_avgs, (label, k,
                        minimum(sn[:,k]), sum(sn[:,k])/ne, maximum(sn[:,k]),
                        minimum(snv[:,k]), sum(snv[:,k])/ne, maximum(snv[:,k]),
                        minimum(sp[:,k]), sum(sp[:,k])/ne, maximum(sp[:,k]),
                        minimum(spv[:,k]), sum(spv[:,k])/ne, maximum(spv[:,k])))
                end
            end
            # ④ 径向剖面：elem1 的 NE/PE 粒子径向 c/r/FE 权重（相位末，来自 final_state y）
            y = state["y"]
            nc_chem = case.index["SPMe chemical block length"]
            off1 = 0 * nc_chem   # elem1 偏移
            idx_n = case.index["negative particle lithium concentration"]
            idx_p = case.index["positive particle lithium concentration"]
            mesh_n, mesh_p = case.mesh["negative particle"], case.mesh["positive particle"]
            cn1 = y[off1 .+ idx_n]; cp1 = y[off1 .+ idx_p]
            for (side, mesh, c) in (("NE", mesh_n, cn1), ("PE", mesh_p, cp1))
                r = vec(mesh.node)
                for i in eachindex(r)   # 5.7.3-③：节点表只存节点
                    push!(radial_rows, (label, side, i, r[i], c[i]))
                end
                q = vec(mesh.gs.x).^2 .* mesh.gs.weight .* mesh.gs.detJ
                c_gs = vec(sum(mesh.gs.Ni .* c[mesh.element[mesh.gs.ele,:]], dims=2))
                for g in eachindex(q)  # 5.7.3-③：完整 GP 表（全部积分点+宿主元素）
                    push!(gp_rows, (label, side, g, mesh.gs.ele[g], mesh.gs.x[g], q[g], c_gs[g]))
                end
            end
            println("[DIAG] phase=$label done t_end=$(pr.t_end) GP3365 D=$(round(gp_key_quantities(case.mech).D; digits=10))")
        end
    finally
        JuBat.mechanical_step_observer[] = old_obs
    end

    # capture_steps（即时 GP1/GP2 分离 + 历史包络 + 单元均值 + BC 违约，全带口径命名）
    hdr_cap = ("phase","step","dt_s","residual","iterations","bc_violation_max",
        "D3365_gp1","D_visc3365_gp1","delta_max_eff3365_gp1_hist","delta_max_n3365_gp1_hist_m","delta_max_t3365_gp1_hist_m",
        "D_elem3365_mean",
        "sep_n3365_elemavg_dimless","sep_t3365_elemavg_dimless","trac_n3365_elemavg_dimless","trac_t3365_elemavg_dimless",
        "gp1_dn_instant_dimless","gp1_dt_instant_dimless","gp1_dn_instant_m","gp1_dt_instant_m",
        "gp2_dn_instant_dimless","gp2_dt_instant_dimless","gp2_dn_instant_m","gp2_dt_instant_m")
    rows_cap = [ (s.phase, k, s.dt, s.resid, s.iters, s.bcviol,
        s.after.D, s.after.D_visc, s.after.δ_max_eff, s.after.δ_max_n*δ_czm, s.after.δ_max_t*δ_czm,
        s.after.D_elem_mean,
        s.sep_n_elemavg, s.sep_t_elemavg, s.trac_n_elemavg, s.trac_t_elemavg,
        s.g1_dn, s.g1_dt, s.g1_dn*δ_czm, s.g1_dt*δ_czm,
        s.g2_dn, s.g2_dt, s.g2_dn*δ_czm, s.g2_dt*δ_czm)
        for (k,s) in enumerate(steps) ]
    write_csv_checked(joinpath(DIAG_OUT, "capture_steps.csv"), hdr_cap, rows_cap)

    # ① averages 12 列满行 + 硬门
    write_csv_checked(joinpath(DIAG_OUT, "averages_t_eval.csv"),
        ("phase","t_end_s",
         "soc_n_nodeavg_e1_first_teval","soc_n_volavg_e1_first_teval",
         "soc_p_nodeavg_e1_first_teval","soc_p_volavg_e1_first_teval",
         "soc_n_nodeavg_e1_last_teval","soc_n_volavg_e1_last_teval",
         "soc_p_nodeavg_e1_last_teval","soc_p_volavg_e1_last_teval",
         "soc_n_nodeavg_span_allfield_alltime","soc_n_volavg_span_allfield_alltime"), soc_avgs)
    write_csv_checked(joinpath(DIAG_OUT, "radial_node_profiles.csv"),
        ("phase","electrode","node_i","r_normalized","c_normalized"), radial_rows)
    write_csv_checked(joinpath(DIAG_OUT, "radial_gp_profiles.csv"),
        ("phase","electrode","gp_i","host_element","r_gp","q_r2wdetj","c_gp"), gp_rows)
    write_csv_checked(joinpath(DIAG_OUT, "averages_stepwise.csv"),
        ("phase","step_k","soc_n_node_min","soc_n_node_mean","soc_n_node_max",
         "soc_n_vol_min","soc_n_vol_mean","soc_n_vol_max",
         "soc_p_node_min","soc_p_node_mean","soc_p_node_max",
         "soc_p_vol_min","soc_p_vol_mean","soc_p_vol_max"), stepwise_avgs)

    # 5.7.4：化学步事件表（候选/接受/拒绝计数与步长来源）
    isempty(chem_events) || write_csv_checked(joinpath(DIAG_OUT, "solver_step_events.csv"),
        ("phase","kind","h_used_or_h","dt_next","trial_idx","accepted_idx"),
        [(e.phase, string(e.kind), e.h_used === nothing ? e.h : e.h_used,
          e.dt_next === nothing ? "" : e.dt_next,
          e.trial_idx === nothing ? "" : e.trial_idx,
          e.accepted_idx === nothing ? "" : e.accepted_idx) for e in chem_events])
    println("[DIAG] chem events: $(count(e->e.kind===:init_BE, chem_events)) init_BE, " *
        "$(count(e->e.kind===:attempt_begin, chem_events)) attempts, " *
        "$(count(e->e.kind===:rejected, chem_events)) rejected, " *
        "$(count(e->e.kind===:accepted_commit, chem_events)) accepted")

    dis_steps = [s for s in steps if s.phase == "discharge"]
    rest1_steps = [s for s in steps if s.phase == "rest1"]
    isempty(rest1_steps) && error("no rest1 steps captured")
    last_dis = dis_steps[end]
    println("[DIAG] captured steps=$(length(steps)) rest1=$(length(rest1_steps))")

    save_binary(joinpath(DIAG_OUT, "diag_packet.bin"),
        (; schema="rest1_root_diag v4", boundary,
           state_discharge_end = phase_end_state["discharge"],
           ms_discharge_end = phase_end_ms["discharge"],
           rest1_inputs = [(; dt = s.dt, F_ext = s.F_ext, dT = s.dT, dsn = s.dsn, dsp = s.dsp) for s in rest1_steps],
           last_discharge_input = (; dt = last_dis.dt, F_ext = last_dis.F_ext, dT = last_dis.dT,
               dsn = last_dis.dsn, dsp = last_dis.dsp),
           gp_elem = DIAG_GP_ELEM, created_at = string(now())))

    if get(ENV, "DIAG_FAULT_BASE", "0") == "1"
        perturbed = merge(rest1_steps[3].after, (D = rest1_steps[3].after.D + 1e-9,))
        rest1_steps[3] = merge(rest1_steps[3], (after = perturbed,))
        println("[DIAG] FAULT INJECTED: captured D at step 3 perturbed +1e-9")
    end

    # --- B: R_BASE（与 capture 同 schema 的差异表）---
    last_stage[] = "BASE"
    last_phase[] = "rest1"           # 5.7.3-④：重放阶段为 rest1
    ms_base = deepcopy(phase_end_ms["discharge"])
    workspace_base = JuBat.CZMWorkspace(cm,param)
    base_rows = Any[]
    for (k,s) in enumerate(rest1_steps)
        last_accepted_step[] = k     # 5.7.3-④：BASE 自己的索引，非 capture 计数
        r = JuBat.solve_czm_step(cm, ms_base, param, s.F_ext, czm; czm_workspace=workspace_base,
            dT_elem = s.dT, Δsoc_n_elem = s.dsn, Δsoc_p_elem = s.dsp,
            dt_seconds = s.dt, boundary = boundary)
        r.converged || error("R_BASE step $k did not converge")
        q = gp_key_quantities(ms_base)
        gp = committed_gp_separations(cm, ms_base.u_prev, DIAG_GP_ELEM, param)
        push!(base_rows, (k, s.dt, s.after.D, q.D, abs(q.D - s.after.D),
            s.after.D_visc, q.D_visc, abs(q.D_visc - s.after.D_visc),
            s.after.δ_max_eff, q.δ_max_eff, abs(q.δ_max_eff - s.after.δ_max_eff),
            r.residual_norm, r.iterations, bc_violation(ms_base.u_prev, boundary),
            gp[1][1], gp[1][2], gp[1][1]*δ_czm, gp[1][2]*δ_czm,
            gp[2][1], gp[2][2], gp[2][1]*δ_czm, gp[2][2]*δ_czm))
    end
    hdr_base = ("step","dt_s","D_cap","D_rep","dD","D_visc_cap","D_visc_rep","dD_visc",
        "dmax_eff_cap","dmax_eff_rep","ddmax_eff","residual_rep","iters_rep","bc_violation_rep",
        "gp1_dn_rep_dimless","gp1_dt_rep_dimless","gp1_dn_rep_m","gp1_dt_rep_m",
        "gp2_dn_rep_dimless","gp2_dt_rep_dimless","gp2_dn_rep_m","gp2_dt_rep_m")
    write_csv_checked(joinpath(DIAG_OUT, "r_base_compare.csv"), hdr_base, base_rows)
    max_dD = maximum(r[5] for r in base_rows)
    max_dDv = maximum(r[8] for r in base_rows)
    max_dδ = maximum(r[11] for r in base_rows)
    println("[DIAG] R_BASE: $(length(base_rows)) steps, max|dD|=$max_dD max|dD_visc|=$max_dDv max|dδ_max_eff|=$max_dδ")
    if !(max_dD == 0.0 && max_dDv == 0.0 && max_dδ == 0.0)
        error("R_BASE replay mismatch: max|dD|=$max_dD max|dD_visc|=$max_dDv max|dδ|=$max_dδ — 因果诊断停止（r_base_compare.csv 已保存）")
    end

    # --- C: R_HOLD_ALL ---
    last_stage[] = "HOLD"
    ms_hold = deepcopy(phase_end_ms["discharge"])
    workspace_hold = JuBat.CZMWorkspace(cm,param)
    hold_rows = Any[]
    t_acc = 0.0
    for (k,s) in enumerate(rest1_steps)
        last_accepted_step[] = k     # 5.7.3-④：HOLD 自己的索引
        r = JuBat.solve_czm_step(cm, ms_hold, param, last_dis.F_ext, czm; czm_workspace=workspace_hold,
            dT_elem = last_dis.dT, Δsoc_n_elem = last_dis.dsn, Δsoc_p_elem = last_dis.dsp,
            dt_seconds = s.dt, boundary = boundary)
        r.converged || error("R_HOLD_ALL step $k did not converge")
        t_acc += s.dt
        q = gp_key_quantities(ms_hold)
        gp = committed_gp_separations(cm, ms_hold.u_prev, DIAG_GP_ELEM, param)
        push!(hold_rows, (k, t_acc, q.D, q.D_visc, q.δ_max_n, q.δ_max_t, q.δ_max_eff,
            r.residual_norm, r.iterations, bc_violation(ms_hold.u_prev, boundary),
            gp[1][1], gp[1][2], gp[1][1]*δ_czm, gp[1][2]*δ_czm,
            gp[2][1], gp[2][2], gp[2][1]*δ_czm, gp[2][2]*δ_czm))
    end
    write_csv_checked(joinpath(DIAG_OUT, "r_hold_all.csv"),
        ("step","t_s","D3365_gp1","D_visc3365_gp1","dmax_n_hist_dimless","dmax_t_hist_dimless","dmax_eff_hist_dimless",
         "residual_rep","iters_rep","bc_violation_rep",
         "gp1_dn_dimless","gp1_dt_dimless","gp1_dn_m","gp1_dt_m",
         "gp2_dn_dimless","gp2_dt_dimless","gp2_dn_m","gp2_dt_m"), hold_rows)
    hold_end = hold_rows[end]
    cap_end = rest1_steps[end].after
    dis_end = last_dis.after
    hold_gp1 = (hold_rows[end][11], hold_rows[end][12])
    cap_gp1 = (rest1_steps[end].g1_dn, rest1_steps[end].g1_dt)
    println("[DIAG] R_HOLD_ALL: end D=$(hold_end[3]) (actual $(cap_end.D), start $(dis_end.D))")

    last_stage[] = "report"
    opt_dump = Dict{String,Any}(string(f) => (v = getfield(case.opt, f); v isa Function ? "callable" : string(v))
        for f in fieldnames(typeof(case.opt)))
    cyc_dump = Dict{String,Any}("I_charge"=>cycle_opt.I_charge, "I_discharge"=>cycle_opt.I_discharge,
        "t_charge"=>cycle_opt.t_charge, "t_discharge"=>cycle_opt.t_discharge, "t_rest1"=>cycle_opt.t_rest1,
        "t_rest2"=>cycle_opt.t_rest2, "V_upper"=>case.param.cell.v_h, "V_lower"=>case.param.cell.v_l,
        "SOC_init"=>cycle_opt.SOC_init, "n_cycles"=>cycle_opt.n_cycles,
        "reset_T_each_cycle"=>cycle_opt.reset_T_each_cycle, "reset_T_before_charge"=>cycle_opt.reset_T_before_charge,
        "dt_cycle"=>cycle_opt.dt_cycle)
    write_manifest("completed"; extra = Dict(
        "source_hash" => bytes2hex(sha256(string(source_identity()))),
        "config_hash" => serialized_hash(options_identity(czm)),
        "packet_sha256" => bytes2hex(sha256(read(joinpath(DIAG_OUT, "diag_packet.bin")))),
        "boundary_nodes" => length(boundary.nodes), "gp_elem" => DIAG_GP_ELEM,
        "delta_czm_m" => δ_czm, "sigma_czm_Pa" => σ_czm,
        "I_typ_A" => param.scale.I_typ, "t0_s" => param.scale.t0,
        "captured_steps" => length(steps), "rest1_steps" => length(rest1_steps),
        "ntheta" => 80, "gsorder" => 2,
        "phase_completion" => "4/4", "wall_seconds" => round(time() - t_all; digits=1),
        "option_dump" => opt_dump, "cycle_option_dump" => cyc_dump,
        "mesh_ne_nlen" => case.mesh["negative particle"].nlen,
        "mesh_pe_nlen" => case.mesh["positive particle"].nlen,
        "mesh_ne_gsorder" => case.mesh["negative particle"].gs.order,
        "n_bulk" => size(cm.bulk_element, 1), "n_cohesive" => length(cm.cohesive_elements),
        "initial_SOC" => 0.65))

    # ⑥ 报告：引用真实文件，列起点/末点/增量/差值与分离比较
    open(joinpath(DIAG_OUT, "report.md"), "w") do f
        println(f, "# rest1 root-cause diagnostics v5 (5.7.3)\n")
        println(f, "- Files: capture_steps.csv / averages_t_eval.csv / averages_stepwise.csv / radial_node_profiles.csv / radial_gp_profiles.csv /",
                " r_base_compare.csv / r_hold_all.csv / diag_packet.bin（write_csv_checked 硬门）。")
        println(f, "- Capture: $(length(steps)) steps (rest1 $(length(rest1_steps))); GP elem $DIAG_GP_ELEM; ",
                "即时 GP1/GP2 分离由 committed_gp_separations（生产 B_global→R→Λ 算法只读重建）。")
        println(f, "- R_BASE: max|ΔD|=$max_dD, max|ΔD_visc|=$max_dDv, max|Δδ_max_eff|=$max_dδ（逐位一致门通过）。")
        println(f, "\n## 实际 vs 冻结输入（rest1 600 s）\n")
        # GP1 共同起点：放电末最后一步的即时分离
        dis_gp1 = (last_dis.g1_dn, last_dis.g1_dt)
        println(f, "| 量 | 共同起点（放电末） | 实际序列末 | 冻结输入末 | 实际增量 | 冻结增量 | 末差 |")
        println(f, "|---|---:|---:|---:|---:|---:|---:|")
        println(f, "| D (GP1) | $(dis_end.D) | $(cap_end.D) | $(hold_end[3]) | $(cap_end.D - dis_end.D) | $(hold_end[3] - dis_end.D) | $(cap_end.D - hold_end[3]) |")
        println(f, "| δ_max_eff (hist, dimless) | $(dis_end.δ_max_eff) | $(cap_end.δ_max_eff) | $(hold_end[7]) | $(cap_end.δ_max_eff - dis_end.δ_max_eff) | $(hold_end[7] - dis_end.δ_max_eff) | $(cap_end.δ_max_eff - hold_end[7]) |")
        println(f, "| GP1 δ_n instant (m) | $(dis_gp1[1]*δ_czm) | $(cap_gp1[1]*δ_czm) | $(hold_gp1[1]*δ_czm) | $(cap_gp1[1]*δ_czm - dis_gp1[1]*δ_czm) | $(hold_gp1[1]*δ_czm - dis_gp1[1]*δ_czm) | $(cap_gp1[1]*δ_czm - hold_gp1[1]*δ_czm) |")
        println(f, "| GP1 δ_t instant (m) | $(dis_gp1[2]*δ_czm) | $(cap_gp1[2]*δ_czm) | $(hold_gp1[2]*δ_czm) | $(cap_gp1[2]*δ_czm - dis_gp1[2]*δ_czm) | $(hold_gp1[2]*δ_czm - dis_gp1[2]*δ_czm) | $(cap_gp1[2]*δ_czm - hold_gp1[2]*δ_czm) |")
        println(f, "\n- 解释边界：本试验未分别隔离黏性与储能；固定输入仍出现显著增长仅说明该模型、",
                "旧生产平均形成的共同放电末态与 τ=5 s 下，冻结输入足以产生显著后续增长。",
                "不把差异拆分为温度/化学贡献百分比，不称必要驱动。")
        println(f, "- 平均口径记录见 averages_t_eval.csv（12 列含 PE 末值）与 averages_stepwise.csv（全场逐步）。",
                "radial_node_profiles.csv / radial_gp_profiles.csv 提供 elem1 NE/PE 径向节点/完整 GP 表。")
        println(f, "- Wall $(round(time() - t_all; digits=1)) s.")
    end
    println("[DIAG] DONE wall=$(round(time() - t_all; digits=1))s")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() || exit(1)
end
