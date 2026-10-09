#=
任务63前置诊断（用户2026-10-10指令）：多种力学边界 × 真实耦合一圈，
采集"带载相位损伤最大时刻"的支路电流分布，供电流云图绘制。

协议与任务62采集口径一致（=任务63 BRIDGE_1）：SPMe 逐单元、distributed2D、
surface 冷却、SOC 0.65、±5 A、1800/600/1800/600 s、dt_cycle=[1,10]、
gp_basic+粘性τ=5 s、geo/J2 关、界面强度×0.09/刚度×0.1/Gc×2、
area_loss_enabled=true（损伤→面积→分流电流真实耦合）、reset_T_before_charge=true。
机械边界经公开 CzmOptions 字段（endpoint_variant 等）进入耦合链路；
每工况导出解析后的约束集合作入口核验证据。

机械重放不含电流（载荷为重放输入），因此电流分布必须来自本真实耦合运行；
不同边界→不同损伤→面积反馈→不同电流重分布。

用法：
  julia -t2 --startup-file=no --project=. example/run_boundary_current_clouds.jl [--case=<id>]
默认运行全部 6 个工况；输出 output/run_boundary_current_clouds/<case>/。
"损伤最大时刻"定义：带载相位（放电/充电）末步的 D_max 最大者（损伤单调不减，
带载末即带载窗口内损伤峰值；充电中途失败时取最后完成的带载相位末并如实标注）。
=#
include(joinpath(@__DIR__, "check_mechanical_boundary_branches.jl"))

const BCC_OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "run_boundary_current_clouds"))

const BOUNDARY_CURRENT_CASES = [
    (; id="fixed_xy_S0_E1_omit_a", outer=:fixed_xy, fix_inner=false, fix_start=false, fix_end=true, variant=:omit_a),
    (; id="fixed_xy_S0_E1_keep",   outer=:fixed_xy, fix_inner=false, fix_start=false, fix_end=true, variant=:keep),
    (; id="fixed_xy_S0_E1_omit_b", outer=:fixed_xy, fix_inner=false, fix_start=false, fix_end=true, variant=:omit_b),
    (; id="fixed_xy_S1_E1_keep",   outer=:fixed_xy, fix_inner=false, fix_start=true,  fix_end=true, variant=:keep),
    (; id="radial_slide_S0_E1_keep", outer=:radial_slide, fix_inner=false, fix_start=false, fix_end=true, variant=:keep),
    (; id="fixed_xy_S0_E0_keep",   outer=:fixed_xy, fix_inner=false, fix_start=false, fix_end=false, variant=:keep),
]

function export_case_geometry(case, outdir)
    tm = case.mesh["thermal2D"]
    nodes = tm.node
    radii = hypot.(nodes[:,1], nodes[:,2])
    # 热网格坐标若为归一化量纲（半径 ~O(1)），乘 scale.L 还原米
    coords_are_meters = maximum(radii) < 0.5
    factor = coords_are_meters ? 1.0 : case.param.scale.L
    write_csv(joinpath(outdir,"thermal_nodes.csv"),("node","x_m","y_m"),
        ((i, nodes[i,1]*factor, nodes[i,2]*factor) for i in 1:size(nodes,1)))
    write_csv(joinpath(outdir,"thermal_elements.csv"),("element","n1","n2","n3","n4"),
        ((e, tm.element[e,1], tm.element[e,2], tm.element[e,3], tm.element[e,4])
         for e in 1:size(tm.element,1)))
    cm = case.czm_mesh
    write_csv(joinpath(outdir,"cohesive_elements.csv"),
        ("element","iface","bottom1","bottom2","host_inner_bulk","host_outer_bulk"),
        ((c.id, c.interface_type, c.nodes_bottom[1], c.nodes_bottom[2],
          c.host_inner_elem, c.host_outer_elem) for c in cm.cohesive_elements))
    return size(tm.element,1), coords_are_meters, factor
end

function export_thermal_damage(case, outdir, tag)
    # 热单元损伤 = 其映射 cohesive 单元 D 均值的最大值（生产分流消费同口径的粗化）
    tm = case.mesh["thermal2D"]
    n_t = size(tm.element,1)
    td = zeros(n_t)
    for (k,c) in enumerate(case.czm_mesh.cohesive_elements)
        D = case.mech.damage_states[k].D
        D > 0 || continue
        for bulk in (c.host_inner_elem, c.host_outer_elem)
            te = case.czm_mesh.czm_submesh.thermal_elem_map[bulk]
            D > td[te] && (td[te] = D)
        end
    end
    write_csv(joinpath(outdir,"thermal_element_damage_$tag.csv"),("element","D_element_mean_max"),
        ((e, td[e]) for e in 1:n_t))
    return td
end

function export_gp_damage(case, outdir, tag)
    ms = case.mech
    ms.gp_damage_states === nothing && error("gp_basic state missing at $tag")
    G = ms.gp_damage_states
    size(G,1) == length(ms.damage_states) ||
        error("gp damage matrix rows $(size(G,1)) != cohesive elements $(length(ms.damage_states))")
    write_csv(joinpath(outdir,"gp_damage_$tag.csv"),("element","gp","D","D_visc","fractured"),
        ((e,g,G[e,g].D,G[e,g].D_visc,G[e,g].fractured)
         for e in 1:size(G,1) for g in 1:size(G,2)))
    return nothing
end

function run_boundary_case(spec, outdir)
    mkpath(outdir)
    case = build_collection_case(outdir)
    czm = case.opt.czm
    czm.fix_inner, czm.fix_start, czm.fix_end = spec.fix_inner, spec.fix_start, spec.fix_end
    czm.outer_bc, czm.endpoint_variant = spec.outer, spec.variant
    boundary = JuBat.resolve_mechanical_bc(case.czm_mesh, case.param, czm; opt=case.opt)
    write_csv(joinpath(outdir,"constraint_nodes.csv"),("node","mode","sources"),
        ((n,m,join(get(boundary.provenance,n,Symbol[]),"|")) for (n,m) in sort(collect(boundary.nodes))))
    ne_thermal, coords_meters, coord_factor = export_case_geometry(case, outdir)

    cycle_opt = diagnostic_cycle_options()
    specs = (("discharge",JuBat.PHASE_DISCHARGE,cycle_opt.t_discharge,cycle_opt.I_discharge,cycle_opt.V_lower),
        ("rest1",JuBat.PHASE_REST,cycle_opt.t_rest1,0.0,0.0),
        ("charge",JuBat.PHASE_CHARGE,cycle_opt.t_charge,-cycle_opt.I_charge,cycle_opt.V_upper),
        ("rest2",JuBat.PHASE_REST,cycle_opt.t_rest2,0.0,0.0))
    state = nothing
    rows, snapshots = Any[], Any[]
    status, error_text = "completed_4_phases", ""
    for (label,phase,duration,current,limit) in specs
        if label=="charge" && cycle_opt.reset_T_before_charge && case.opt.thermalmodel != "none"
            JuBat.reset_cycle_temperature!(case,state)
        end
        println("[$(spec.id)] phase=$label start t_max=$duration current=$current")
        try
            pr = JuBat.solve_phase(case,phase,duration,current,limit,state;
                ms=case.mech,dt_range=cycle_opt.dt_cycle)
            state = pr.final_state
            sr = pr.solve_result
            t = sr["time [s]"]; V = sr["cell voltage [V]"]; D = sr["czm D_max"]
            I = sr["thermal2D element current"]
            size(I,1) == ne_thermal || error("current rows $(size(I,1)) != thermal elements $ne_thermal")
            for k in eachindex(t)
                col = view(I,:,k)
                push!(rows,(label,k,t[k],V[k],D[k],minimum(col),maximum(col),
                    sum(abs,col)/length(col)))
            end
            # 相位末快照：带载相位末同时保留电流列与损伤场
            Iend = I[:,end]
            snapshot = (; phase=label, t_end=pr.t_end, D_end=D[end], terminated_by=string(pr.terminated_by),
                currents=copy(Iend), damage=copy(case.mech.damage_states))
            push!(snapshots,snapshot)
            write_csv(joinpath(outdir,"current_phase_end_$label.csv"),("element","I_A"),
                ((e,Iend[e]) for e in 1:length(Iend)))
            export_thermal_damage(case,outdir,label)
            println("[$(spec.id)] CAPTURE_PHASE=$label t_end=$(pr.t_end) terminated_by=$(pr.terminated_by) D_end=$(D[end]) steps=$(length(t))")
        catch err
            status = "failed_in_$label"
            error_text = sprint(showerror,err)
            println("[$(spec.id)] FAILED in $label: $error_text")
            break
        end
    end
    write_csv(joinpath(outdir,"steps.csv"),
        ("phase","step","time_s","voltage_V","D_max","I_min_A","I_max_A","I_mean_abs_A"),rows)
    export_gp_damage(case,outdir,"final_state")

    # 带载相位损伤最大时刻：放电/充电相位末 D_end 最大者
    loaded = [s for s in snapshots if s.phase in ("discharge","charge")]
    chosen = isempty(loaded) ? nothing : loaded[argmax([s.D_end for s in loaded])]
    manifest = Dict("case"=>spec.id,"outer_bc"=>string(spec.outer),"fix_inner"=>spec.fix_inner,
        "fix_start"=>spec.fix_start,"fix_end"=>spec.fix_end,"endpoint_variant"=>string(spec.variant),
        "constrained_nodes"=>length(boundary.nodes),"outer_count"=>boundary.outer_count,
        "a_node"=>boundary.a,"a_fixed"=>haskey(boundary.nodes,boundary.a),
        "b_node"=>boundary.b,"b_fixed"=>haskey(boundary.nodes,boundary.b),
        "thermal_elements"=>ne_thermal,"coords_are_meters"=>coords_meters,"coord_factor_m"=>coord_factor,
        "status"=>status,"error"=>error_text,
        "chosen_phase"=>chosen===nothing ? "" : chosen.phase,"chosen_t_end_s"=>chosen===nothing ? 0.0 : chosen.t_end,
        "chosen_D_max"=>chosen===nothing ? 0.0 : chosen.D_end,
        "chosen_terminated_by"=>chosen===nothing ? "" : chosen.terminated_by,
        "note"=>"chosen instant = loaded phase end with max D (damage monotone); current field in current_phase_end_<phase>.csv",
        "protocol"=>"task62 capture recipe: SPMe per-element, gp_basic, tau=5s, area_loss on, reset_T_before_charge=true, ntheta=80",
        "created_at"=>string(now()))
    open(joinpath(outdir,"manifest.toml"),"w") do io
        TOML.print(io,manifest)
    end
    ok = status == "completed_4_phases" && chosen !== nothing
    println("[$(spec.id)] DONE status=$status chosen=$(chosen===nothing ? "none" : chosen.phase*"@t="*string(chosen.t_end))")
    return ok
end

function main(args=ARGS)
    wanted = String[]
    for arg in args
        startswith(arg,"--case=") || throw(ArgumentError("unknown argument: $arg"))
        push!(wanted,split(arg,"=";limit=2)[2])
    end
    mkpath(BCC_OUTPUT)
    all_ok = true
    for spec in BOUNDARY_CURRENT_CASES
        (isempty(wanted) || spec.id in wanted) || continue
        all_ok &= run_boundary_case(spec, joinpath(BCC_OUTPUT,spec.id))
    end
    return all_ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() || exit(1)
end
