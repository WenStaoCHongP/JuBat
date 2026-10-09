#=
任务62：真实一次完整循环输入采集及独立力学边界重放。

单次入口（Julia 1.11.2，建议 -t1，GKSwstype=100）：
  julia --startup-file=no --project=. example/check_mechanical_boundary_branches.jl

默认采集完整一次循环：1800s放电→600s静置→1800s充电→600s静置。
仅调试时使用：--max-steps=1 或 --window-seconds=60；这两项明确属于部分窗口。
可选命令行参数：--collect-only
  --packet=<本脚本此前保存的 input_packet.bin> 只重放已采集包，核对源码身份。
  --max-steps 控制采集的实际已接受机械步数；不改真实自适应 dt。

所有数据写入 @__DIR__/../output/check_mechanical_boundary_branches/。
该实验隔离同一真实循环载荷下的力学边界效应；不属于容量或寿命验证。
=#
using Printf, Serialization, SHA, Dates, LinearAlgebra, Statistics, TOML, Plots
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat

const BC_DIAGNOSTIC_ROOT = normpath(joinpath(@__DIR__, ".."))
const BC_DIAGNOSTIC_OUTPUT = joinpath(BC_DIAGNOSTIC_ROOT, "output", "check_mechanical_boundary_branches")

struct MechanicalCaptureComplete <: Exception end

function diagnostic_arguments(args)
    max_steps, window_seconds, packet_path, collect_only = nothing, nothing, nothing, false
    for arg in args
        if startswith(arg, "--max-steps=")
            max_steps = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--window-seconds=")
            window_seconds = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--packet=")
            packet_path = abspath(split(arg, "="; limit=2)[2])
        elseif arg == "--collect-only"
            collect_only = true
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    max_steps === nothing || max_steps > 0 || throw(ArgumentError("--max-steps must be positive"))
    window_seconds === nothing || (isfinite(window_seconds) && window_seconds > 0) ||
        throw(ArgumentError("--window-seconds must be finite and positive"))
    return (; max_steps, window_seconds, packet_path, collect_only)
end

function serialized_hash(value)
    io = IOBuffer()
    serialize(io, value)
    return bytes2hex(sha256(take!(io)))
end

function save_binary(path, value)
    open(path, "w") do io
        serialize(io, value)
    end
    return bytes2hex(sha256(read(path)))
end

function csv_cell(value)
    value === nothing && return ""
    text = string(value)
    occursin(r"[,\"\r\n]", text) && return "\"" * replace(text, "\"" => "\"\"") * "\""
    return text
end

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(csv_cell.(header), ","))
        for row in rows
            println(io, join(csv_cell.(Tuple(row)), ","))
        end
    end
end

function source_identity()
    paths = sort([joinpath(dir,file) for (dir,_,files) in walkdir(joinpath(BC_DIAGNOSTIC_ROOT,"src"))
        for file in files if endswith(file,".jl")])
    append!(paths, [joinpath(BC_DIAGNOSTIC_ROOT, "Project.toml"),
        joinpath(@__DIR__, "czm_cycling.jl"), @__FILE__])
    return [(path=relpath(p, BC_DIAGNOSTIC_ROOT), sha256=bytes2hex(sha256(read(p)))) for p in paths]
end

damage_identity(s) = Tuple(getfield(s, f) for f in fieldnames(typeof(s)))
function mechanical_state_identity(ms)
    return (; u=ms.u_prev, element=map(damage_identity, ms.damage_states),
        gp=ms.gp_damage_states === nothing ? nothing : map(damage_identity, ms.gp_damage_states),
        plastic=ms.plastic_states, prestress=ms.winding_prestress, node_ref=ms.node_ref,
        contact=ms.contact)
end

function static_mesh_identity(mesh)
    sub = mesh.czm_submesh
    return (; nodes=mesh.node, bulk=mesh.bulk_element,
        cohesive=[(; id=e.id, nodes=e.nodes, bottom=e.nodes_bottom, top=e.nodes_top,
            iface=e.interface_type, inner=e.host_inner_elem, outer=e.host_outer_elem) for e in mesh.cohesive_elements],
        thermal_map=sub.thermal_elem_map, material=sub.material_type,
        winding_turn=sub.winding_turn, phi_pairs=sub.phi_pairs, phi_keep=sub.phi_keep,
        bonded_nodes=sub.mesh_bonded.node, bonded_elements=sub.mesh_bonded.element,
        gsorder=mesh.bulk_mesh.gs.order)
end

function clone_static_mesh(mesh)
    result = deepcopy(mesh)
    result.K_bulk = nothing
    result.cohesive_geom = nothing
    result.bulk_gp_geom = nothing
    result.bulk_asm = nothing
    result.ws = nothing
    return result
end

function options_identity(opt)
    return Dict(string(f) => (getfield(opt, f) === nothing ? "legacy" : string(getfield(opt, f)))
        for f in fieldnames(typeof(opt)))
end

function build_collection_case(outdir)
    # 原脚本在 SetCase 前进行这些物理参数缩放；归一化后的 param 冻结。
    pd = JuBat.ChooseCell("Jellyroll")
    pd.cell.v_l, pd.cell.v_h = 2.5, 4.2
    for ip in (pd.PCC, pd.NCC)
        ip.σ_max *= 0.09; ip.τ_max *= 0.09
        ip.K_n *= 0.1; ip.K_t *= 0.1
        ip.G_c *= 2.0; ip.G_c_t *= 2.0
        ip.δ_0 = ip.σ_max / ip.K_n
        ip.δ_0_t = ip.τ_max / ip.K_t
        ip.δ_c = 2.0 * ip.G_c / ip.σ_max
        ip.δ_c_t = 2.0 * ip.G_c_t / ip.τ_max
    end
    pd.scale.σ_czm = pd.PCC.σ_max
    pd.scale.δ_czm = 2.0 * pd.PCC.G_c / pd.PCC.σ_max
    pd.scale.G_czm = pd.scale.σ_czm * pd.scale.δ_czm
    pd.scale.K_czm = pd.scale.σ_czm / pd.scale.δ_czm
    opt = JuBat.Option()
    opt.Current = _ -> 5.0
    opt.model = "SPMe"
    opt.Nn, opt.Ns, opt.Np, opt.Nrn, opt.Nrp = 10, 5, 10, 10, 10
    opt.gsorder, opt.dimension = 2, 1
    opt.mechanicalmodel = "none"
    opt.dt = [0.5, 10.0] # 原示例原值；solve_phase 将按 dt_cycle=[1,10] 设置。
    opt.dtType, opt.jacobi, opt.solveType = "auto", "update", "Crank-Nicolson"
    opt.thermal_enabled, opt.thermalmodel, opt.thermal_dim = true, "distributed2D", "2D"
    opt.cool_method, opt.per_element_spme = "surface", true
    opt.debug_coupling = false
    opt.debug_log_path = joinpath(outdir, "collection_coupling.log")
    opt.czm.enabled, opt.czm.model, opt.czm.fix_inner = true, "mix", false
    opt.czm.iter_method, opt.czm.max_iter = "gp_basic", 100
    opt.czm.load_steps, opt.czm.tol, opt.czm.update_interval = 10, 1e-3, 1
    opt.czm.geo_nonlinear, opt.czm.j2_plasticity = false, false
    opt.czm.viscous_enabled, opt.czm.viscous_tau = true, 5.0
    opt.czm.area_loss_enabled = true
    case = JuBat.SetCase(pd, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    JuBat.apply_initial_soc!(case, case.param_dim, 0.65)
    return case
end

function diagnostic_cycle_options()
    return JuBat.CycleOption(n_cycles=1,I_charge=5.0,I_discharge=5.0,t_charge=1800.0,
        t_discharge=1800.0,t_rest1=600.0,t_rest2=600.0,V_upper=4.2,V_lower=2.5,
        SOC_init=0.65,reset_T_each_cycle=false,reset_T_before_charge=true)
end

function collect_input_packet(settings, outdir, sources)
    isdefined(JuBat, :mechanical_step_observer) || error("internal mechanical_step_observer seam is required")
    case = build_collection_case(outdir)
    inputs, baseline_after, baseline_results = Any[], Any[], Any[]
    pending = Ref{Any}(nothing)
    elapsed = Ref(0.0)
    phase_elapsed,phase_offset = Ref(0.0),Ref(0.0)
    phase_label = Ref("discharge")
    captured_param, captured_options = Ref{Any}(nothing), Ref{Any}(nothing)
    failed = Ref{Any}(nothing)
    old_observer = JuBat.mechanical_step_observer[]
    old_observer === nothing || error("another mechanical observer is already active")
    capture_error = Ref("")
    observer = function(event)
        if event.phase == :before
            settings.max_steps===nothing || length(inputs) < settings.max_steps || throw(MechanicalCaptureComplete())
            event.dt_seconds === nothing && error("capture received no physical dt_seconds")
            dt = Float64(event.dt_seconds)
            isfinite(dt) && dt > 0 || error("capture received invalid physical dt_seconds=$dt")
            if captured_param[] === nothing
                captured_param[] = deepcopy(event.param)
                captured_options[] = deepcopy(event.czm_opt)
            end
            before = deepcopy(event.ms)
            if !isempty(baseline_after) && isequal(mechanical_state_identity(before),mechanical_state_identity(baseline_after[end]))
                before = baseline_after[end] # 同一不可修改提交快照共享存储，减少完整循环包内存。
            end
            pending[] = (; index=length(inputs) + 1, cycle=1, phase=phase_label[],
                time_s=phase_offset[] + phase_elapsed[] + dt, phase_time_s=phase_elapsed[] + dt,
                dt_seconds=dt, F_ext=copy(event.F_ext),
                dT_elem=copy(event.dT_elem), Δsoc_n_elem=copy(event.Δsoc_n_elem),
                Δsoc_p_elem=copy(event.Δsoc_p_elem), ms_before=before,
                eigenstrain=deepcopy(event.eigenstrain), prestress=deepcopy(event.prestress))
        elseif event.phase == :after
            pending[] === nothing && error("capture :after without :before")
            if event.result.converged
                push!(inputs, pending[])
                push!(baseline_after, deepcopy(event.ms))
                push!(baseline_results, deepcopy(event.result))
                elapsed[] = pending[].time_s
                phase_elapsed[] = pending[].phase_time_s
            else
                capture_error[] = "baseline mechanical solve did not converge: $(event.result.residual_norm)"
                failed[] = pending[]
            end
            pending[] = nothing
        elseif event.phase == :error
            event.error isa MechanicalCaptureComplete || (capture_error[] = sprint(showerror, event.error))
        else
            error("unknown mechanical observer event $(event.phase)")
        end
    end
    stop_reason = "complete_cycle"
    completed_cycles = 0
    phases = Any[]
    chemical_thermal_state = nothing
    cycle_opt = diagnostic_cycle_options()
    try
        JuBat.mechanical_step_observer[] = observer
        specs = (("discharge",JuBat.PHASE_DISCHARGE,cycle_opt.t_discharge,cycle_opt.I_discharge,cycle_opt.V_lower),
            ("rest1",JuBat.PHASE_REST,cycle_opt.t_rest1,0.0,0.0),
            ("charge",JuBat.PHASE_CHARGE,cycle_opt.t_charge,-cycle_opt.I_charge,cycle_opt.V_upper),
            ("rest2",JuBat.PHASE_REST,cycle_opt.t_rest2,0.0,0.0))
        for (label,phase,duration,current,limit) in specs
            if settings.window_seconds !== nothing && label != "discharge"
                stop_reason = "explicit_partial_discharge_window"
                break
            end
            if label=="charge" && cycle_opt.reset_T_before_charge && case.opt.thermalmodel != "none"
                JuBat.reset_cycle_temperature!(case,chemical_thermal_state)
            end
            phase_label[],phase_elapsed[] = label,0.0
            phase_offset[] = chemical_thermal_state===nothing ? 0.0 : Float64(chemical_thermal_state["t_global"])
            t_max = settings.window_seconds===nothing ? duration : settings.window_seconds
            previous_state_hash = serialized_hash(chemical_thermal_state)
            phase_result = JuBat.solve_phase(case,phase,t_max,current,limit,chemical_thermal_state;
                ms=case.mech,dt_range=cycle_opt.dt_cycle)
            chemical_thermal_state = phase_result.final_state
            push!(phases,(; phase=label,start_s=phase_result.t_start,end_s=phase_result.t_end,
                duration_s=phase_result.duration,current_A=current,terminated_by=string(phase_result.terminated_by),
                initial_chemical_thermal_hash=previous_state_hash,final_chemical_thermal_hash=serialized_hash(chemical_thermal_state),
                case_initial_SOC_setting=cycle_opt.SOC_init,reset_T_before_phase=(label=="charge" && cycle_opt.reset_T_before_charge)))
            println("CAPTURE_PHASE=$label COMPLETED_PHASES=$(length(phases))/4 TIME_S=$(phase_result.t_end) ACCEPTED_MECHANICAL_STEPS=$(length(inputs))")
        end
        completed_cycles = length(phases)==4 ? 1 : 0
    catch err
        if err isa MechanicalCaptureComplete
            stop_reason = "accepted_step_limit"
        else
            capture_error[] = sprint(showerror, err)
            stop_reason = "collection_failed"
        end
    finally
        JuBat.mechanical_step_observer[] = old_observer
    end
    # 失败输入单独保留，不能冒充已接受时间窗口。
    failed_input = failed[] === nothing ? pending[] : failed[]
    mesh = clone_static_mesh(case.czm_mesh)
    param = captured_param[] === nothing ? deepcopy(case.param) : captured_param[]
    options = captured_options[] === nothing ? deepcopy(case.opt.czm) : captured_options[]
    config = (; origin="example/czm_cycling.jl one complete cycle", ntheta=80, current_A=5.0,
        SOC_init=0.65, dt_cycle=cycle_opt.dt_cycle, phase_duration_s=(1800.0,600.0,1800.0,600.0),
        requested_window_s=settings.window_seconds,
        max_accepted_steps=settings.max_steps, strength_scale=0.09, stiffness_scale=0.1,
        fracture_energy_scale=2.0, area_loss_enabled=true, reset_T_each_cycle=false,reset_T_before_charge=true,
        czm_options=options_identity(options))
    packet = (; schema="JuBat mechanical boundary packet v1", created_at=string(now()),
        sources, source_hash=serialized_hash(sources), config, config_hash=serialized_hash(config),
        mesh, mesh_hash=serialized_hash(static_mesh_identity(mesh)), param,
        param_hash=serialized_hash(param), options, inputs, baseline_after, baseline_results,
        failed_input, last_committed=deepcopy(case.mech), actual_end_time_s=elapsed[],
        completed_cycles,completed_phases=phases,final_chemical_thermal_state=chemical_thermal_state,
        stop_reason, collection_error=capture_error[])
    sha = save_binary(joinpath(outdir, "input_packet.bin"), packet)
    write_csv(joinpath(outdir, "input_steps.csv"),
        ("step", "cycle", "phase", "time_s", "phase_time_s", "dt_seconds", "F_ext_hash", "dT_hash", "delta_soc_n_hash",
            "delta_soc_p_hash", "before_state_hash", "committed_state_hash", "GP_count", "D_max_GP", "D_max_element"),
        ((x.index, x.cycle, x.phase, x.time_s, x.phase_time_s, x.dt_seconds, serialized_hash(x.F_ext),
            serialized_hash(x.dT_elem), serialized_hash(x.Δsoc_n_elem), serialized_hash(x.Δsoc_p_elem),
            serialized_hash(mechanical_state_identity(x.ms_before)),
            serialized_hash(mechanical_state_identity(packet.baseline_after[k])),
            length(packet.baseline_after[k].gp_damage_states),
            maximum(s.D for s in packet.baseline_after[k].gp_damage_states),
            maximum(s.D for s in packet.baseline_after[k].damage_states)) for (k,x) in enumerate(inputs)))
    println("INPUT_PACKET_SHA256=$sha")
    write_csv(joinpath(outdir,"capture_phases.csv"),("phase","start_s","end_s","duration_s","current_A","terminated_by",
        "initial_chemical_thermal_hash","final_chemical_thermal_hash","case_initial_SOC_setting","reset_T_before_phase"),(Tuple(p) for p in phases))
    println("CAPTURE_ACCEPTED_STEPS=$(length(inputs)) ACTUAL_END_TIME_S=$(elapsed[]) COMPLETED_CYCLES=$completed_cycles STOP=$stop_reason")
    return packet, sha
end

function validate_packet(packet, sources)
    packet.schema == "JuBat mechanical boundary packet v1" || error("unsupported packet schema")
    packet.source_hash == serialized_hash(sources) || error("packet source identity differs from current source")
    packet.mesh_hash == serialized_hash(static_mesh_identity(packet.mesh)) || error("packet mesh hash mismatch")
    packet.param_hash == serialized_hash(packet.param) || error("packet param hash mismatch")
    packet.config_hash == serialized_hash(packet.config) || error("packet config hash mismatch")
    length(packet.inputs) == length(packet.baseline_after) == length(packet.baseline_results) ||
        error("packet input/state/result lengths differ")
    isempty(packet.inputs) && error("no accepted mechanical inputs were captured")
    for (k,x) in enumerate(packet.inputs)
        x.index == k || error("noncontiguous packet indices")
        k == 1 || isequal(mechanical_state_identity(x.ms_before),
            mechanical_state_identity(packet.baseline_after[k-1])) || error("baseline committed-state chain differs at step $k")
        isfinite(x.dt_seconds) && x.dt_seconds > 0 || error("invalid packet dt at step $k")
        all(isfinite, x.F_ext) && all(isfinite, x.dT_elem) &&
            all(isfinite, x.Δsoc_n_elem) && all(isfinite, x.Δsoc_p_elem) || error("nonfinite packet input at step $k")
        packet.baseline_after[k].gp_damage_states === nothing && error("accepted gp_basic state lacks actual GP history")
    end
    return nothing
end

function candidate_matrix(packet)
    candidates = Any[]
    # 旧调用先重放，逐步检查其与采集状态完全一致。
    opt = deepcopy(packet.options)
    push!(candidates, (; id="legacy", outer="fixed_xy", start="legacy", finish="legacy", variant="keep", options=opt))
    for outer in (:fixed_xy, :radial_slide), start in (false,true), finish in (false,true), variant in (:keep,:omit_a,:omit_b,:omit_ab)
        opt = deepcopy(packet.options)
        opt.fix_inner = false
        opt.outer_bc, opt.fix_start, opt.fix_end = outer, start, finish
        id = "$(outer)_S$(Int(start))_E$(Int(finish))_$(variant)"
        push!(candidates, (; id, outer=string(outer), start=string(start), finish=string(finish), variant=string(variant), options=opt))
    end
    return candidates
end

function boundary_operator_identity(boundary)
    # 后端签名包含实际投影算子及数值规范，不以配置名称或来源标签去重。
    return JuBat.mechanical_boundary_signature(boundary)
end

enforcement_tag(candidate) = candidate.id=="legacy" ?
    "legacy_cartesian_penalty" : "directional_orthogonal_elimination"

function unfolded_angles(mesh, param)
    nseg = maximum(mesh.czm_submesh.thermal_elem_map)
    theta = Vector{Float64}(undef, nseg+1)
    nfirst = mesh.bulk_element[1,1]
    b = param.cell.layer / (2π)
    theta[1] = (hypot(mesh.node[nfirst,1], mesh.node[nfirst,2]) - param.cell.Rin) / b
    previous = atan(mesh.node[nfirst,2], mesh.node[nfirst,1])
    for s in 1:nseg
        node = mesh.bulk_element[s,4]
        current = atan(mesh.node[node,2], mesh.node[node,1])
        increment = mod(current - previous, 2π)
        0 < increment < π || error("nonmonotone inherited angular topology at segment $s")
        theta[s+1] = theta[s] + increment
        previous = current
    end
    return theta
end

function gp_fields(mesh, param, ms, options, theta; only_damaged=false)
    ms.gp_damage_states === nothing && error("GP fields require committed GP states")
    geom = JuBat.cohesive_geometry(mesh)
    nseg = length(theta)-1
    theta_span = theta[end]-theta[1]
    fields = Any[]
    for (e,elem) in enumerate(mesh.cohesive_elements)
        only_damaged && all(s -> s.D <= 1e-8 && s.D_visc <= 1e-8,@view(ms.gp_damage_states[e,:])) && continue
        g = geom[e]
        ue = ms.u_prev[g.dofs]
        segment = mesh.czm_submesh.thermal_elem_map[elem.host_inner_elem]
        face = (elem.host_inner_elem - 1) ÷ nseg + 1
        ip = JuBat.collector_params(param, elem.interface_type)
        for (gp,(xi,weight)) in enumerate(zip(g.gauss_pts,g.gauss_wts))
            state = ms.gp_damage_states[e,gp]
            only_damaged && state.D <= 1e-8 && state.D_visc <= 1e-8 && continue
            N1, N2 = (1-xi)/2, (1+xi)/2
            jump = [N1*(ue[7]-ue[1]) + N2*(ue[5]-ue[3]),
                N1*(ue[8]-ue[2]) + N2*(ue[6]-ue[4])]
            delta = (param.scale.L / param.scale.δ_czm) .* (g.R * jump)
            tn,tt,_,_ = JuBat.bilinear_traction_state(delta[1],delta[2],state,ip,options.model;visc_beta=0.0)
            n1,n2 = elem.nodes_bottom
            x = (N1*mesh.node[n1,1]+N2*mesh.node[n2,1]) * param.scale.L
            y = (N1*mesh.node[n1,2]+N2*mesh.node[n2,2]) * param.scale.L
            angle = N1*theta[segment] + N2*theta[segment+1]
            fraction = (angle-theta[1])/theta_span
            region = fraction < 1/3 ? "inner" : fraction < 2/3 ? "middle" : "outer"
            push!(fields, (; element=e, gp, iface=string(elem.interface_type), face, segment,
                winding_turn=mesh.czm_submesh.winding_turn[elem.host_inner_elem], theta=angle,
                region, x, y, radius=hypot(x,y), D=state.D, D_visc=state.D_visc,
                delta_n_m=delta[1]*param.scale.δ_czm, delta_t_m=delta[2]*param.scale.δ_czm,
                traction_n_Pa=tn*param.scale.σ_czm, traction_t_Pa=tt*param.scale.σ_czm,
                history_n_m=state.δ_max_n*param.scale.δ_czm,
                history_t_m=state.δ_max_t*param.scale.δ_czm,
                history_eff_m=state.δ_max_eff*param.scale.δ_czm,
                fractured=state.fractured, accumulated_damage=state.accumulated_damage,
                length_weight_m=g.length*weight*0.5*param.scale.L))
        end
    end
    return fields
end

function compact_GP_metrics(mesh,ms,theta)
    geom = JuBat.cohesive_geometry(mesh)
    region_max = zeros(3)
    weighted,total_length = 0.0,0.0
    d_max,dv_max,damaged,fractured = 0.0,0.0,0,0
    for (e,el) in enumerate(mesh.cohesive_elements)
        segment = mesh.czm_submesh.thermal_elem_map[el.host_inner_elem]
        for (g,s) in enumerate(@view(ms.gp_damage_states[e,:]))
            xi = geom[e].gauss_pts[g]
            angle = (1-xi)/2*theta[segment] + (1+xi)/2*theta[segment+1]
            fraction = (angle-theta[1])/(theta[end]-theta[1])
            region = fraction < 1/3 ? 1 : fraction < 2/3 ? 2 : 3
            region_max[region] = max(region_max[region],s.D)
            weight = geom[e].length*geom[e].gauss_wts[g]/2
            weighted += s.D*weight
            total_length += weight
            d_max = max(d_max,s.D)
            dv_max = max(dv_max,s.D_visc)
            damaged += s.D>1e-8
            fractured += s.fractured
        end
    end
    return (; d_max,dv_max,mean_length_weighted=weighted/total_length,damaged,fractured,
        inner_max=region_max[1],middle_max=region_max[2],outer_max=region_max[3],
        element_max=maximum(s.D for s in ms.damage_states))
end

function committed_equilibrium(mesh, param, ms, options, input, boundary)
    _,fint,_,_ = JuBat.assemble_coupled_system(mesh, ms.u_prev,param;
        damage_states=ms.damage_states,gp_damage_states=ms.gp_damage_states,
        K_bulk_cached=JuBat.bulk_stiffness(mesh,param),geom_cache=JuBat.cohesive_geometry(mesh),
        ws=JuBat.assembly_workspace(mesh),visc_beta=0.0,czm_model=options.model)
    thermal = JuBat.assemble_thermal_chemical_load(mesh,param,input.dT_elem,input.Δsoc_n_elem,input.Δsoc_p_elem)
    raw_residual = input.F_ext + thermal - fint
    return raw_residual, JuBat.mechanical_boundary_diagnostics(ms.u_prev,raw_residual,boundary)
end

function export_geometry(mesh, param, theta, outdir)
    write_csv(joinpath(outdir,"mesh_nodes.csv"),("node","x_m","y_m"),
        ((n,mesh.node[n,1]*param.scale.L,mesh.node[n,2]*param.scale.L) for n in 1:mesh.nnode))
    write_csv(joinpath(outdir,"mesh_bulk_elements.csv"),("element","n1","n2","n3","n4","material","parent_thermal","winding_turn"),
        ((e,mesh.bulk_element[e,:]...,mesh.czm_submesh.material_type[e],mesh.czm_submesh.thermal_elem_map[e],mesh.czm_submesh.winding_turn[e]) for e in axes(mesh.bulk_element,1)))
    write_csv(joinpath(outdir,"mesh_cohesive_elements.csv"),("element","interface","bottom1","bottom2","top1","top2","host_inner","host_outer","parent_thermal","length_m"),
        ((e,el.interface_type,el.nodes_bottom...,el.nodes_top...,el.host_inner_elem,el.host_outer_elem,
            mesh.czm_submesh.thermal_elem_map[el.host_inner_elem],el.length*param.scale.L) for (e,el) in enumerate(mesh.cohesive_elements)))
    write_csv(joinpath(outdir,"inherited_theta.csv"),("angular_node","theta_unfolded_rad"),enumerate(theta))
end

function export_branch_fields(branch_dir, fields, element_states)
    if isempty(fields)
        write(joinpath(branch_dir,"GP_unavailable.txt"),"No accepted GP solve in this branch; the saved initial state has gp_damage_states=nothing. No GP history was reconstructed.\n")
    else
        header = string.(propertynames(first(fields)))
        write_csv(joinpath(branch_dir,"committed_gp.csv"),header,(Tuple(f) for f in fields))
    end
    write_csv(joinpath(branch_dir,"committed_element_damage.csv"),("element","D_weighted_average","D_visc_weighted_average","fractured_all_GP"),
        ((i,s.D,s.D_visc,s.fractured) for (i,s) in enumerate(element_states)))
    isempty(fields) && return
    region_rows = Any[]
    for iface in ("PE_PCC","NE_NCC"), region in ("inner","middle","outer")
        fs = filter(f -> f.iface == iface && f.region == region,fields)
        isempty(fs) && error("empty geometric region $iface/$region")
        total_length = sum(f.length_weight_m for f in fs)
        damaged = filter(f -> f.D > 1e-8,fs)
        push!(region_rows,(iface,region,length(fs),maximum(f.D for f in fs),
            sum(f.D*f.length_weight_m for f in fs)/total_length,
            sum((f.length_weight_m for f in damaged);init=0.0),length(damaged),
            count(f -> f.fractured,fs),isempty(damaged) ? nothing : minimum(f.theta for f in damaged),
            isempty(damaged) ? nothing : maximum(f.theta for f in damaged)))
    end
    write_csv(joinpath(branch_dir,"regional_damage.csv"),("interface","region","GP_count","D_max_GP","D_length_weighted_mean",
        "damaged_GP_length_m","damaged_GP_count_D_gt_1e_8","fractured_GP_count","damaged_theta_min_rad","damaged_theta_max_rad"),region_rows)
end

function shared_field_plots(results,outdir)
    results = filter(r -> !isempty(r.fields), results)
    isempty(results) && return
    # GP峰值与单元均值分离。统一[0,1]损伤色标；分离/牵引使用所有独立分支的同一范围。
    all_fields = reduce(vcat,[r.fields for r in results])
    channels = ((:D,"Committed GP D",(0.0,1.0)),(:D_visc,"Committed GP D_visc",(0.0,1.0)),
        (:delta_n_m,"Current GP delta_n [m]",extrema(f.delta_n_m for f in all_fields)),
        (:delta_t_m,"Current GP delta_t [m]",extrema(f.delta_t_m for f in all_fields)),
        (:traction_n_Pa,"Current GP traction_n [Pa]",extrema(f.traction_n_Pa for f in all_fields)),
        (:traction_t_Pa,"Current GP traction_t [Pa]",extrema(f.traction_t_Pa for f in all_fields)))
    for r in results
        plots = Any[]
        for (field,label,limits) in channels
            lo,hi = limits
            lo == hi && ((lo,hi) = (lo-1e-15,hi+1e-15))
            p = scatter([f.x*1e3 for f in r.fields],[f.y*1e3 for f in r.fields];
                marker_z=[getfield(f,field) for f in r.fields],clims=(lo,hi),color=:viridis,
                markersize=1.7,markerstrokewidth=0,aspect_ratio=:equal,legend=false,colorbar=true,
                xlabel="x [mm]",ylabel="y [mm]",title=label,xlims=(-10.5,10.5),ylims=(-10.5,10.5))
            push!(plots,p)
        end
        figure = plot(plots...;layout=(2,3),size=(1800,1200),margin=7*Plots.mm,
            titlefontsize=11,guidefontsize=10,tickfontsize=8,
            plot_title="$(r.id): last committed t=$(r.end_time)s ($(r.status))")
        savefig(figure,joinpath(outdir,"branches",r.id,"committed_gp_fields.png"))
        p = scatter([f.theta/(2π) for f in r.fields],[f.D for f in r.fields];
            marker_z=[f.D for f in r.fields],clims=(0,1),color=:viridis,markersize=2,
            markerstrokewidth=0,legend=false,colorbar=true,xlabel="Unfolded angle / 2pi",ylabel="Committed GP D",
            ylims=(0,1),title="$(r.id) - actual GP locations",size=(1100,400))
        savefig(p,joinpath(outdir,"branches",r.id,"unfolded_GP_damage.png"))
    end
end

function diagnostic_report(packet,packet_sha,candidates,rows,results,outdir)
    unique_count = length(results)
    completed = count(r -> r.status == "completed",results)
    failed = count(r -> r.status != "completed",results)
    onset = maximum(s.D for s in packet.baseline_after[end].gp_damage_states)
    open(joinpath(outdir,"report.md"),"w") do io
        println(io,"# Task 62 mechanical boundary replay\n")
        println(io,"Captured **$(length(packet.inputs)) accepted mechanical steps**, ending at **$(packet.actual_end_time_s) s**. Completed real coupled cycles: **$(packet.completed_cycles)/1**; completed phases: **$(length(packet.completed_phases))/4**. Stop reason: `$(packet.stop_reason)`.")
        println(io,"\nDefault collection runs 1800 s discharge → 600 s rest → 1800 s charge → 600 s rest with carried chemical/thermal state and MechState. Temperature is reset immediately before charge exactly as in the original example. Any --max-steps or --window-seconds run is an explicitly partial diagnostic and cannot satisfy the one-cycle verification request.")
        println(io,"\nConfiguration: nθ=80, +5 A, SOC_init=0.65, gp_basic, τ=5 s, geo/J2=false, area_loss_enabled=true, strength/stiffness/Gc scale=0.09/0.1/2. Physical dt is captured from each real solve call. dT is Kelvin; electrode concentration changes are mol/m³.")
        println(io,"\nCandidates: $(length(candidates)) (32 explicit + 1 legacy). Unique replay cases including constraint operator, numerical rotation gauge and enforcement method: $unique_count. Completed unique replays: $completed; failed: $failed. Duplicate candidates are aliases in summary.csv and do not count as independent numerical evidence.")
        println(io,"\nThe legacy reference uses its original finite-penalty Cartesian enforcement; explicit configurations use orthogonal elimination. They remain distinct replay cases even when their resolved physical constraint projectors match. candidate_options.csv records both the physical operator hash and the full solve identity hash.")
        println(io,"\nBaseline committed GP D_max: $onset. ",onset > 1e-8 ? "Damage started within the captured cycle/window." : "The captured data contain no damage onset; damage redistribution cannot be inferred.")
        legacy = first(results)
        println(io,"\nThe legacy reference is checked against every captured committed displacement and complete GP/element/other MechState field. Legacy verification status: **$(legacy.status)**. Each unique replay starts from its own copy of the first captured state and consumes exactly the same input sequence. Replay damage cannot modify the captured electrochemical or thermal loads.")
        println(io,"\nEach branch saves its last committed MechState even on failure. Failed trial states are never exported as accepted solutions. Free residual and constraint violation refer to the committed state under the last successfully accepted load; the solver failure residual is reported separately.")
        println(io,"\n`committed_gp.csv` contains real GP D/D_visc, current signed separations, tractions, history maxima and fracture flags. `committed_element_damage.csv` contains GP-weighted element means. D>0.1 is softening, not complete fracture. Region membership uses inherited unfolded angle thirds, with winding turn and radius included; θ modulo 2π is not used to identify inner/outer turns.")
        println(io,"\nGP export sampling: every successfully finished phase end, first damage onset, global GP peak and last committed state contain all GP data. `compact_GP_history.csv` contains metrics for every accepted step. `damaged_GP_history.csv` contains full per-step mode/history data only for actual GP with D or D_visc > 1e-8; it is not a table of all undamaged GP at every step. Complete baseline committed GP histories are retained in input_packet.bin.")
        println(io,"\nAll damage plots use [0,1]. Current separation/traction plots share ranges across independent branches. Branches that failed earlier are labelled with their actual last committed time, so these fields must not be interpreted as a common terminal-time comparison.")
        println(io,"\nReactions use the assembled raw normalized residual and are split into local radial/tangential directions. Constraint provenance and numerical rotation gauge are exported separately. A rotation gauge is a numerical convention and cannot replace physical support or carry an incompatible applied torque.")
        println(io,"\nThis is a fixed-load mechanical causality experiment. The real baseline cycle completion count is reported above; branch replays do not re-run electrochemical or thermal feedback. The experiment does not establish physical support conditions, contact/friction, inner-core buckling, full-coupling current/capacity changes, 12-cycle completion or lifetime validity. Selecting a/b rules requires user review of the reaction and damage evidence; hotspot movement is not an acceptance condition.")
        println(io,"\n| Unique branch | Status | Accepted steps | Last committed time [s] | GP D_max | Element D_max |")
        println(io,"|---|---|---:|---:|---:|---:|")
        for r in results
            gp_max = isempty(r.fields) ? "unavailable" : maximum(f.D for f in r.fields)
            println(io,"| $(r.id) | $(r.status) | $(r.accepted) | $(r.end_time) | $gp_max | $(maximum(s.D for s in r.ms.damage_states)) |")
        end
        println(io,"\nPacket SHA-256: `$packet_sha`\n\nSource: `$(packet.source_hash)`\n\nConfiguration: `$(packet.config_hash)`\n\nMesh: `$(packet.mesh_hash)`\n\nParameters: `$(packet.param_hash)`")
        isempty(packet.collection_error) || println(io,"\nCollection failure: `$(replace(packet.collection_error,'`'=>'\''))`")
        println(io,"\nFiles: summary.csv, candidate_options.csv, constraint_contributions.csv, constraint_nodes.csv, capture_phases.csv, input_steps.csv, source_identity.csv, identity.toml, input_packet.bin, mesh CSVs, branches/<id>/ and this report.")
    end
    println("CANDIDATES=$(length(candidates)) UNIQUE_REPLAYS=$unique_count COMPLETED=$completed FAILED=$failed")
    return failed == 0 && isempty(packet.collection_error)
end

function export_candidate_boundaries(candidates,boundaries,owners,operators,mesh,param,outdir)
    option_rows, contribution_rows, node_rows, topology_rows = Any[], Any[], Any[], Any[]
    for (candidate,boundary,owner,operator) in zip(candidates,boundaries,owners,operators)
        opts = options_identity(candidate.options)
        resolved_start = candidate.options.fix_start === nothing ? candidate.options.fix_inner : candidate.options.fix_start
        resolved_end = candidate.options.fix_end === nothing ?
            (candidate.options.fix_inner ? "full_end" : "end_except_b") : string(candidate.options.fix_end)
        dofs,_ = JuBat.mechanical_boundary_dofs(boundary)
        gauge = dofs.gauge
        count_constraints = sum(mode == :fixed_xy ? 2 : 1 for mode in values(boundary.nodes))
        # 后端已经检查三刚体模态；只有已验证的纯转动秩2情形建立规范。
        rigid_body_rank = gauge===nothing ? 3 : 2
        push!(option_rows,(candidate.id,owner,candidate.outer,candidate.options.fix_inner,candidate.start,candidate.finish,
            resolved_start,resolved_end,candidate.variant,length(boundary.nodes),count_constraints,
            rigid_body_rank,gauge !== nothing,boundary.inner_count,boundary.outer_count,
            enforcement_tag(candidate),serialized_hash(boundary_operator_identity(boundary)),operator,serialized_hash(opts)))
        roles = Dict(:geometric_outer=>boundary.geometric_outer,:geometric_inner=>boundary.geometric_inner,
            :start=>boundary.start_nodes,:end=>boundary.end_nodes,:dual_a=>Set([boundary.a]),:dual_b=>Set([boundary.b]))
        important_nodes = sort(collect(union(values(roles)...)))
        for node in important_nodes
            x,y = mesh.node[node,1]*param.scale.L,mesh.node[node,2]*param.scale.L
            sources = get(boundary.provenance,node,Symbol[])
            mode = get(boundary.nodes,node,:free)
            r = hypot(x,y)
            rank = mode == :fixed_xy ? 2 : mode == :radial ? 1 : 0
            push!(node_rows,(candidate.id,node,x,y,r,node==boundary.a,node==boundary.b,
                node in boundary.geometric_outer,node in boundary.geometric_inner,node in boundary.start_nodes,
                node in boundary.end_nodes,mode,rank,r==0 ? nothing : x/r,r==0 ? nothing : y/r,join(string.(sources),"|")))
            for role in (:geometric_outer,:geometric_inner,:start,:end,:dual_a,:dual_b)
                node in roles[role] || continue
                retained = role in sources
                contribution_mode = role in (:geometric_outer,:dual_a) && candidate.outer=="radial_slide" ? "radial" : "fixed_xy"
                push!(contribution_rows,(candidate.id,node,role,retained,contribution_mode,x,y))
            end
            if node==boundary.a || node==boundary.b
                bulk = [e for e in axes(mesh.bulk_element,1) if node in @view(mesh.bulk_element[e,:])]
                coh = [e for (e,el) in enumerate(mesh.cohesive_elements) if node in el.nodes]
                push!(topology_rows,(candidate.id,node==boundary.a ? "a" : "b",node,x,y,
                    join(bulk,"|"),join(string.(mesh.czm_submesh.material_type[bulk]),"|"),join(coh,"|")))
            end
        end
        if gauge !== nothing
            gauge_dir = joinpath(outdir,"numerical_gauges")
            mkpath(gauge_dir)
            write_csv(joinpath(gauge_dir,candidate.id*".csv"),("node","gauge_x","gauge_y"),
                ((n,gauge[2n-1],gauge[2n]) for n in 1:mesh.nnode))
        end
    end
    write_csv(joinpath(outdir,"candidate_options.csv"),("candidate","unique_owner","outer_bc_raw","fix_inner_raw","fix_start_raw",
        "fix_end_raw","fix_start_resolved","fix_end_resolved","endpoint_variant","constrained_nodes","independent_physical_constraints",
        "rigid_body_rank","rotation_gauge","inner_count","outer_count","enforcement_method",
        "physical_operator_hash","solve_identity_hash","raw_options_hash"),option_rows)
    write_csv(joinpath(outdir,"constraint_contributions.csv"),("candidate","node","source","retained","contribution_mode","x_m","y_m"),contribution_rows)
    write_csv(joinpath(outdir,"constraint_nodes.csv"),("candidate","node","x_m","y_m","radius_m","is_a","is_b",
        "geometric_outer","geometric_inner","start","end","resolved_mode","independent_rank","radial_nx","radial_ny","retained_sources"),node_rows)
    write_csv(joinpath(outdir,"special_point_topology.csv"),("candidate","special_point","node","x_m","y_m","bulk_elements",
        "connected_materials","cohesive_elements"),topology_rows)
end

function export_reactions(branch_dir,mesh,param,ms,raw,diagnostics,boundary)
    # raw=F_ext+F_thermochemical-f_internal；物理支承反力方向取 -Q*raw。
    support = -diagnostics.reaction
    write_csv(joinpath(branch_dir,"committed_node_reactions.csv"),("node","x_m","y_m","ux_m","uy_m",
        "support_rx_normalized","support_ry_normalized","support_radial_normalized","support_tangent_normalized",
        "raw_residual_x","raw_residual_y","sources"),
        ((n,mesh.node[n,1]*param.scale.L,mesh.node[n,2]*param.scale.L,
            ms.u_prev[2n-1]*param.scale.L,ms.u_prev[2n]*param.scale.L,
            support[2n-1],support[2n],
            dot(support[2n-1:2n],mesh.node[n,:])/hypot(mesh.node[n,1],mesh.node[n,2]),
            dot(support[2n-1:2n],[-mesh.node[n,2],mesh.node[n,1]])/hypot(mesh.node[n,1],mesh.node[n,2]),
            raw[2n-1],raw[2n],join(string.(get(boundary.provenance,n,Symbol[])),"|")) for n in 1:mesh.nnode))
end

function replay_branch(candidate,boundary,packet,theta,outdir)
    branch_dir = joinpath(outdir,"branches",candidate.id)
    mkpath(branch_dir)
    mesh = clone_static_mesh(packet.mesh)
    param,ms,options = packet.param,deepcopy(packet.inputs[1].ms_before),deepcopy(candidate.options)
    accepted,end_time,status,error_text,last_residual = 0,0.0,"completed","",nothing
    step_rows,history_rows,compact_rows = Any[],Any[],Any[]
    onset_snapshot,peak_snapshot = nothing,nothing
    peak_D = -Inf
    last_diagnostics,last_raw = nothing,nothing
    t_start = time_ns()
    # 关闭采集回调，分支只调用力学求解，绝不回写共同输入包。
    JuBat.mechanical_step_observer[] === nothing || error("observer unexpectedly active during replay")
    for (k,input) in enumerate(packet.inputs)
        before = deepcopy(ms)
        result = nothing
        local solve_error = nothing
        try
            result = JuBat.solve_czm_step(mesh,ms,param,copy(input.F_ext),options;
                dT_elem=copy(input.dT_elem),Δsoc_n_elem=copy(input.Δsoc_n_elem),Δsoc_p_elem=copy(input.Δsoc_p_elem),
                eigenstrain=input.eigenstrain,prestress=input.prestress,dt_seconds=input.dt_seconds,
                boundary=candidate.id=="legacy" ? nothing : boundary)
        catch err
            solve_error = err
        end
        if solve_error !== nothing || !result.converged
            status = "failed_solve"
            error_text = solve_error === nothing ? "mechanical result.converged=false" : sprint(showerror,solve_error)
            unchanged = isequal(mechanical_state_identity(ms),mechanical_state_identity(before))
            if !unchanged
                save_binary(joinpath(branch_dir,"failure_illegally_mutated_state.bin"),ms)
                status = "failed_state_contract"
                error_text *= "; solver changed MechState on failure"
                ms = before
            end
            last_residual = result === nothing ? nothing : result.residual_norm
            push!(step_rows,(k,input.phase,input.time_s,input.phase_time_s,input.dt_seconds,false,result===nothing ? nothing : result.iterations,
                last_residual,nothing,nothing,nothing,nothing,unchanged,nothing,error_text))
            save_binary(joinpath(branch_dir,"failed_input.bin"),(;input,before_state=before,error=error_text,result))
            break
        end
        accepted,end_time,last_residual = k,input.time_s,result.residual_norm
        raw,diagnostics = committed_equilibrium(mesh,param,ms,options,input,boundary)
        last_diagnostics,last_raw = diagnostics,raw
        exact = candidate.id=="legacy" ?
            isequal(mechanical_state_identity(ms),mechanical_state_identity(packet.baseline_after[k])) : nothing
        push!(step_rows,(k,input.phase,input.time_s,input.phase_time_s,input.dt_seconds,true,result.iterations,result.residual_norm,
            diagnostics.free_residual,diagnostics.constraint_violation,diagnostics.gauge_coordinate,
            diagnostics.gauge_moment,nothing,exact,""))
        metrics = compact_GP_metrics(mesh,ms,theta)
        push!(compact_rows,(k,input.phase,input.time_s,input.phase_time_s,Tuple(metrics)...))
        fields = gp_fields(mesh,param,ms,options,theta;only_damaged=true)
        # 保存损伤GP的逐步模式/历史，保留热点驱动力而不生成虚构均值GP。
        for field in fields
            if field.D > 1e-8 || field.D_visc > 1e-8
                push!(history_rows,(k,input.phase,input.time_s,input.phase_time_s,Tuple(field)...))
            end
        end
        if metrics.d_max > peak_D
            peak_D = metrics.d_max
            peak_snapshot = (;index=k,phase=input.phase,time_s=input.time_s,ms=deepcopy(ms))
        end
        if onset_snapshot===nothing && metrics.d_max>1e-8
            onset_snapshot = (;index=k,phase=input.phase,time_s=input.time_s,ms=deepcopy(ms))
        end
        phase_finished = k < length(packet.inputs) ? packet.inputs[k+1].phase!=input.phase :
            any(p -> p.phase==input.phase,packet.completed_phases)
        if phase_finished
            snapshot_dir = joinpath(branch_dir,"snapshots","phase_final_"*input.phase)
            mkpath(snapshot_dir)
            save_binary(joinpath(snapshot_dir,"committed_state.bin"),(;index=k,phase=input.phase,time_s=input.time_s,ms=deepcopy(ms)))
            export_branch_fields(snapshot_dir,gp_fields(mesh,param,ms,options,theta),ms.damage_states)
        end
        if exact === false
            status,error_text = "failed_legacy_identity","legacy replay differs from captured committed MechState at step $k"
            break
        end
        if !(isfinite(diagnostics.free_residual) && diagnostics.free_residual < options.tol &&
             isfinite(diagnostics.constraint_violation) && diagnostics.constraint_violation < options.tol)
            status,error_text = "failed_equilibrium_contract","committed free residual or constraint violation exceeds tolerance at step $k"
            break
        end
    end
    wall_seconds = (time_ns()-t_start)*1e-9
    save_binary(joinpath(branch_dir,"last_committed_state.bin"),(;ms,last_accepted_input_index=accepted,
        last_accepted_time_s=end_time,status,error=error_text,options))
    fields = ms.gp_damage_states === nothing ? Any[] : gp_fields(mesh,param,ms,options,theta)
    export_branch_fields(branch_dir,fields,ms.damage_states)
    for (label,snapshot) in (("damage_onset",onset_snapshot),("global_GP_peak",peak_snapshot))
        snapshot===nothing && continue
        snapshot_dir = joinpath(branch_dir,"snapshots",label)
        mkpath(snapshot_dir)
        save_binary(joinpath(snapshot_dir,"committed_state.bin"),snapshot)
        export_branch_fields(snapshot_dir,gp_fields(mesh,param,snapshot.ms,options,theta),snapshot.ms.damage_states)
    end
    write_csv(joinpath(branch_dir,"step_metrics.csv"),("step","phase","input_time_s","phase_time_s","dt_seconds","accepted","iterations",
        "solver_residual","committed_free_residual","constraint_violation","gauge_coordinate","gauge_moment",
        "failure_state_unchanged","legacy_committed_state_exact","error"),step_rows)
    if !isempty(fields)
        write_csv(joinpath(branch_dir,"damaged_GP_history.csv"),("step","phase","time_s","phase_time_s",string.(propertynames(first(fields)))...),history_rows)
        write_csv(joinpath(branch_dir,"compact_GP_history.csv"),("step","phase","time_s","phase_time_s",
            "D_max_GP","D_visc_max_GP","D_length_weighted_mean","damaged_GP_count","fractured_GP_count",
            "inner_D_max_GP","middle_D_max_GP","outer_D_max_GP","D_max_element_mean"),compact_rows)
    end
    last_raw===nothing || export_reactions(branch_dir,mesh,param,ms,last_raw,last_diagnostics,boundary)
    isempty(error_text) || write(joinpath(branch_dir,"failure.txt"),error_text*"\n")
    println("BRANCH=$(candidate.id) STATUS=$status ACCEPTED=$accepted/$(length(packet.inputs)) LAST_COMMITTED_TIME_S=$end_time WALL_SECONDS=$wall_seconds")
    return (;id=candidate.id,status,accepted,end_time,ms,fields,last_residual,
        diagnostics=last_diagnostics,wall_seconds,error=error_text)
end

function replay_all_branches(packet,packet_sha,outdir)
    mesh,param = packet.mesh,packet.param
    theta = unfolded_angles(mesh,param)
    export_geometry(mesh,param,theta,outdir)
    candidates = candidate_matrix(packet)
    boundaries = [JuBat.resolve_mechanical_bc(mesh,param,c.options;endpoint_variant=Symbol(c.variant)) for c in candidates]
    operators = [serialized_hash((;operator=boundary_operator_identity(b),enforcement=enforcement_tag(c)))
        for (c,b) in zip(candidates,boundaries)]
    owners,seen = String[],Dict{String,String}()
    for (c,key) in zip(candidates,operators)
        owner = get!(seen,key,c.id)
        push!(owners,owner)
    end
    export_candidate_boundaries(candidates,boundaries,owners,operators,mesh,param,outdir)
    results = Any[]
    for (candidate,boundary,owner) in zip(candidates,boundaries,owners)
        candidate.id==owner || continue
        push!(results,replay_branch(candidate,boundary,packet,theta,outdir))
    end
    by_id = Dict(r.id=>r for r in results)
    summary_rows = Any[]
    for (candidate,owner,operator) in zip(candidates,owners,operators)
        r = by_id[owner]
        ds = r.diagnostics
        gp_max = isempty(r.fields) ? nothing : maximum(f.D for f in r.fields)
        gp_visc_max = isempty(r.fields) ? nothing : maximum(f.D_visc for f in r.fields)
        weighted = isempty(r.fields) ? nothing : sum(f.D*f.length_weight_m for f in r.fields)/sum(f.length_weight_m for f in r.fields)
        push!(summary_rows,(candidate.id,owner,candidate.id!=owner,operator,candidate.outer,
            candidate.start,candidate.finish,candidate.variant,r.status,r.accepted,length(packet.inputs),
            r.end_time,gp_max,gp_visc_max,maximum(s.D for s in r.ms.damage_states),weighted,
            isempty(r.fields) ? nothing : count(f->f.D>1e-8,r.fields),
            isempty(r.fields) ? nothing : count(f->f.D>=0.99,r.fields),
            isempty(r.fields) ? nothing : count(f->f.fractured,r.fields),
            r.last_residual,ds===nothing ? nothing : ds.free_residual,
            ds===nothing ? nothing : ds.constraint_violation,ds===nothing ? nothing : ds.gauge_moment,
            r.wall_seconds,r.error))
    end
    write_csv(joinpath(outdir,"summary.csv"),("candidate","unique_owner","duplicate_solve_case","solve_identity_hash",
        "outer_bc","fix_start_raw","fix_end_raw","endpoint_variant","status","accepted_steps","requested_input_steps",
        "last_committed_time_s","D_max_GP","D_visc_max_GP","D_max_element_mean","D_length_weighted_mean",
        "damaged_GP_count_D_gt_1e_8","GP_count_D_ge_0_99","fractured_GP_count","last_solver_residual",
        "last_committed_free_residual","constraint_violation","gauge_moment","unique_replay_wall_seconds","error"),summary_rows)
    shared_field_plots(results,outdir)
    return diagnostic_report(packet,packet_sha,candidates,summary_rows,results,outdir)
end

function main(args=ARGS)
    settings = diagnostic_arguments(args)
    stored_packet = settings.packet_path===nothing ? nothing : open(deserialize,settings.packet_path)
    is_partial = stored_packet===nothing ?
        (settings.max_steps!==nothing || settings.window_seconds!==nothing) : stored_packet.completed_cycles!=1
    outdir = is_partial ? joinpath(BC_DIAGNOSTIC_OUTPUT,"partial_validation") : BC_DIAGNOSTIC_OUTPUT
    mkpath(outdir)
    sources = source_identity()
    write_csv(joinpath(outdir,"source_identity.csv"),("source_relative_path","sha256"),
        ((s.path,s.sha256) for s in sources))
    packet,packet_sha = if settings.packet_path === nothing
        collect_input_packet(settings,outdir,sources)
    else
        (stored_packet,bytes2hex(sha256(read(settings.packet_path))))
    end
    validate_packet(packet,sources)
    open(joinpath(outdir,"identity.toml"),"w") do io
        TOML.print(io,Dict("schema"=>packet.schema,"source_sha256"=>packet.source_hash,
            "configuration_sha256"=>packet.config_hash,"mesh_sha256"=>packet.mesh_hash,
            "param_sha256"=>packet.param_hash,"packet_sha256"=>packet_sha,
            "accepted_input_steps"=>length(packet.inputs),"actual_end_time_s"=>packet.actual_end_time_s,
            "completed_cycles"=>packet.completed_cycles,"completed_phases"=>length(packet.completed_phases),
            "stop_reason"=>packet.stop_reason,"collection_error"=>packet.collection_error))
    end
    settings.collect_only && return isempty(packet.collection_error)
    return replay_all_branches(packet,packet_sha,outdir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() || exit(1)
end
