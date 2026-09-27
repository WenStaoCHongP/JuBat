"""
Jellyroll 电池长循环仿真：多 SPMe 并行电化学-热-CZM 循环（文字结果 + CSV 导出）

自由设定循环次数（JUBAT_CYCLES 环境变量，默认 5），一次运行产出后处理所需的
全部 CSV（逐循环汇总 / 逐相位明细 / 选定循环逐元素全历史 / 选定循环场数据长表）。

界面参数口径（2026-09-27 定版，KB jubat-czm-dev-playbook）：
基线文件不动（σ_max 82e6/92e6、K_n 2.4e17/1.2e17、G_c 25.3/6.2），脚本内
缩放 ×0.09（强度）/ ×0.1（刚度）/ ×2（断裂能）→ 等效 σ_max 7.38/8.28 MPa、
K_n 2.4e16/1.2e16 Pa/m、G_c 50.6/12.4 J/m²。geo/J2 成对关闭（D-B3-1）。

环境变量：
  JUBAT_CYCLES=<N>            总循环次数（默认 5）
  JUBAT_CZM_METHOD=<name>     求解方法（basic / gp_basic / arc_length）
  JUBAT_CZM_TAU_SECONDS=<s>   物理时间粘性 τ（设为 5 时 gp_basic + 粘性 = 稳健口径）
  JUBAT_SNAPSHOT_CYCLES=<csv> 选定循环号逗号列表（如 1,5,10,30,50,100）→ 采全历史
  JUBAT_EXPORT_FIELD_DATA=1   同时导出温度/位移/电流/损伤场数据长表（需 save_detailed）
  JUBAT_RUN_TAG=<name>        输出目录子目录名（默认 debug_a0.9_gc2）

网格信息（节点坐标/单元连接 CSV）不随本脚本导出——网格静态，需时用
`tools/check_collector_mesh.jl`（KB 恢复后运行）单独生成到
`output/网格信息/theta<N>/`（见 KB jubat-czm-dev-playbook 网格工具节）。

推荐运行（100 循环 + 全历史导出）：
  JUBAT_CYCLES=100 JUBAT_CZM_METHOD=gp_basic JUBAT_CZM_TAU_SECONDS=5 \
  JUBAT_SNAPSHOT_CYCLES=1,5,10,30,50,100 JUBAT_EXPORT_FIELD_DATA=1 \
  JUBAT_RUN_TAG=long_run julia -t8 example/czm_cycling.jl

日期：2026-09-14（2026-09-27 合并场数据导出、改为自由循环次数）
"""

using Printf
using Profile
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# JUBAT_SNAPSHOT_CYCLES 解析（空=不采集全历史）
function czm_snapshot_cycles()
    raw = get(ENV, "JUBAT_SNAPSHOT_CYCLES", "")
    isempty(raw) && return Set{Int}()
    return Set(parse(Int, strip(x)) for x in split(raw, ","))
end

function main()
    run_tag = get(ENV, "JUBAT_RUN_TAG", "debug_a0.9_gc2")
    outdir = joinpath(@__DIR__, "..", "output", "czm_cycling", run_tag)
    mkpath(outdir)

    n_cycles = parse(Int, get(ENV, "JUBAT_CYCLES", "5"))
    snapshot_cycles = czm_snapshot_cycles()
    export_field = get(ENV, "JUBAT_EXPORT_FIELD_DATA", "0") == "1"

    println("="^80)
    @printf("Jellyroll电池多SPMe并行电化学-热耦合循环仿真（SOC 0.65 / %d次循环 / 调试参数）\n", n_cycles)
    println("="^80)

    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/4] 参数设置...")

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    # 界面参数口径（定版）：基线文件 × 脚本缩放 ×0.09/×0.1/×2
    strength_scale, stiffness_scale, gc_scale = 0.09, 0.1, 2.0
    for ip in (param_dim.PCC, param_dim.NCC)
        ip.σ_max *= strength_scale
        ip.τ_max *= strength_scale
        ip.K_n *= stiffness_scale
        ip.K_t *= stiffness_scale
        ip.G_c *= gc_scale
        ip.G_c_t *= gc_scale
        ip.δ_0 = ip.σ_max / ip.K_n
        ip.δ_0_t = ip.τ_max / ip.K_t
        ip.δ_c = 2.0 * ip.G_c / ip.σ_max
        ip.δ_c_t = 2.0 * ip.G_c_t / ip.τ_max
    end
    param_dim.scale.σ_czm = param_dim.PCC.σ_max
    param_dim.scale.δ_czm = 2.0 * param_dim.PCC.G_c / param_dim.PCC.σ_max
    param_dim.scale.G_czm = param_dim.scale.σ_czm * param_dim.scale.δ_czm
    param_dim.scale.K_czm = param_dim.scale.σ_czm / param_dim.scale.δ_czm

    opt = JuBat.Option()
    i = 5.0
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"

    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.debug_coupling = false
    opt.debug_log_path = joinpath(outdir, "simple_coupling_debug.log")
    opt.czm.enabled = true
    opt.czm.model = "mix"
    opt.czm.fix_inner = false
    opt.czm.iter_method = get(ENV, "JUBAT_CZM_METHOD", "basic")
    opt.czm.max_iter = parse(Int, get(ENV, "JUBAT_CZM_MAX_ITER", string(opt.czm.max_iter)))
    if lowercase(opt.czm.iter_method) in ("arc_length", "arclength", "arc-length")
        opt.czm.wall_time_limit_seconds = parse(Float64,
            get(ENV, "JUBAT_CZM_WALL_LIMIT_SECONDS", "60"))
        opt.czm.snapshot_path = joinpath(outdir, "arc_failure.bin")
    end
    opt.czm.load_steps = 10
    opt.czm.tol = 1e-3
    opt.czm.area_loss_enabled = get(ENV, "JUBAT_AREA_LOSS", "0") == "1"  # 任务 58 双面连续面积反馈
    opt.czm.geo_nonlinear = false   # 调试口径：geo/J2 成对关闭（D-B3-1）
    opt.czm.j2_plasticity = false
    opt.czm.update_interval = 1
    tau_seconds = get(ENV, "JUBAT_CZM_TAU_SECONDS", "")
    if !isempty(tau_seconds)
        opt.czm.viscous_enabled = true
        opt.czm.viscous_tau = parse(Float64, tau_seconds)
    end

    cycle_opt = JuBat.CycleOption(
        n_cycles = n_cycles,
        I_charge = i,
        I_discharge = i,
        t_charge = 1800.0,
        t_discharge = 1800.0,
        t_rest1 = 600.0,
        t_rest2 = 600.0,
        V_upper = 4.2,
        V_lower = 2.5,
        SOC_init = 0.65,
        reset_T_each_cycle = false,
        reset_T_before_charge = true,
    )

    println("OK: 参数设置完成")
    @printf("  电流: %.2f A（充放同幅）\n", i)
    @printf("  循环: %d 次（放电 1800 s → 静置 600 s → 充电 1800 s → 静置 600 s）\n", n_cycles)
    @printf("  初始SOC: 0.65（solve_cycling 内部施加）\n")
    @printf("  等效界面参数: σ_max=%.2f/%.2f MPa, K=%.1e/%.1e, G_c=%.1f/%.1f J/m²\n",
        param_dim.PCC.σ_max * 1e-6, param_dim.NCC.σ_max * 1e-6,
        param_dim.PCC.K_n, param_dim.NCC.K_n,
        param_dim.PCC.G_c, param_dim.NCC.G_c)
    @printf("  CZM: %s / geo=%s / J2=%s / max_iter=%d / tol=%.1e / τ=%g s / area_loss=%s\n",
        opt.czm.iter_method, opt.czm.geo_nonlinear, opt.czm.j2_plasticity,
        opt.czm.max_iter, opt.czm.tol, opt.czm.viscous_tau, string(opt.czm.area_loss_enabled))
    @printf("  全历史采集循环: %s\n", isempty(snapshot_cycles) ? "无" : sort(collect(snapshot_cycles)) |> x -> join(x, ","))
    @printf("  场数据长表: %s\n", export_field ? "导出" : "不导出")

    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/4] 创建案例和Jellyroll网格...")

    case = JuBat.SetCase(param_dim, opt)

    n_theta = 80
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)

    ne = size(mesh_th.element, 1)
    println("OK: Jellyroll网格创建完成")
    @printf("  周向单元数 n_theta: %d / 总单元数 ne: %d / 总节点数 nT: %d\n", n_theta, ne, mesh_th.nlen)

    # === 设定全局环境（快照循环门控） ===
    if !isempty(snapshot_cycles)
        ENV["JUBAT_SNAPSHOT_CYCLES"] = join(sort(collect(snapshot_cycles)), ",")
    else
        delete!(ENV, "JUBAT_SNAPSHOT_CYCLES")
    end

    # ========================================================================
    # 3. 循环求解
    # ========================================================================
    println("\n[3/4] 运行循环求解器...")

    t_wall_start = time_ns()
    profile_on = get(ENV, "JUBAT_PROFILE", "0") == "1"
    local result
    try
        solve_fn = () -> JuBat.solve_cycling(case, cycle_opt, case.mech;
            save_detailed = export_field)
        if profile_on
            Profile.init(n = 10_000_000, delay = 0.01)
            result = Profile.@profile solve_fn()
        else
            result = solve_fn()
        end
    catch err
        println("\n[诊断] solve_cycling 失败，原位状态转储：")
        println("  错误: " * first(sprint(showerror, err), 200))
        mech = case.mech
        u = mech.u_prev
        if mech.gp_damage_states !== nothing
            gp = mech.gp_damage_states
            @printf("  GP历史: D_max=%.6f, D_visc_max=%.6f, D>=0.99=%d, fractured=%d / %d\n",
                maximum(s.D for s in gp), maximum(s.D_visc for s in gp),
                count(s -> s.D >= 0.99, gp), count(s -> s.fractured, gp), length(gp))
        end
        for (iface,) in ((:PE_PCC,), (:NE_NCC,))
            idx = findall(e -> e.interface_type == iface, case.czm_mesh.cohesive_elements)
            states = mech.damage_states[idx]
            @printf("  %s: fractured=%d, D>0.9=%d, D>0.99=%d, D_max=%.6f / %d\n",
                string(iface), count(s -> s.fractured, states),
                count(s -> s.D > 0.9, states), count(s -> s.D > 0.99, states),
                maximum(s.D for s in states), length(states))
        end
        rethrow(err)
    end
    t_wall_s = (time_ns() - t_wall_start) * 1e-9
    if profile_on
        prof_path = joinpath(outdir, "profile_flat.txt")
        open(prof_path, "w") do io
            Profile.print(io, format = :flat)
        end
        println("profile 已保存: $prof_path")
    end

    println("OK: 循环求解完成")
    @printf("  完成循环数: %d / 计划循环数: %d / 总墙钟: %.1f s\n",
        length(result.cycle_idx), cycle_opt.n_cycles, t_wall_s)

    # ========================================================================
    # 4. 损伤分析与 CSV 导出
    # ========================================================================
    println("\n[4/4] 损伤分析 + CSV 导出")
    println("-"^78)
    @printf("  %-6s %-12s %-12s %-12s %-12s %-10s\n",
        "循环", "放电容量[Ah]", "D_max", "D_mean", "T_max[K]", "SOH")
    for k in 1:length(result.cycle_idx)
        @printf("  %-6d %-12.4f %-12.4e %-12.4e %-12.2f %-10.4f\n",
            result.cycle_idx[k], result.capacity_discharge[k],
            result.D_max[k], result.D_mean[k], result.T_max[k], result.soh[k])
    end

    println()
    if !isempty(result.cycle_results)
        @printf("  阶段损伤增量分解（ΔD_max，相对各阶段起点）:\n")
        @printf("  %-6s %-10s %-14s %-14s %-14s\n", "循环", "阶段", "时长[s]", "D_max(末)", "ΔD_max")
        for cr in result.cycle_results
            for (label, ph) in (("放电", cr.discharge), ("静置1", cr.rest1), ("充电", cr.charge), ("静置2", cr.rest2))
                ph === nothing && continue
                @printf("  %-6d %-10s %-14.0f %-14.4e %-14.4e\n",
                    cr.cycle_idx, label, ph.duration, ph.D_max, ph.ΔD_max)
            end
        end
    end

    # --- CSV 导出：逐循环汇总 + 逐相位明细 ---
    open(joinpath(outdir, "cycle_summary.csv"), "w") do io
        println(io, "cycle_idx,capacity_charge_Ah,capacity_discharge_Ah,coulombic_efficiency,D_max,D_mean,n_fractured,T_max_K,soh")
        for k in eachindex(result.cycle_idx)
            println(io, join((
                result.cycle_idx[k], result.capacity_charge[k], result.capacity_discharge[k],
                result.coulombic_efficiency[k], result.D_max[k], result.D_mean[k],
                result.n_fractured[k], result.T_max[k], result.soh[k]), ","))
        end
    end
    open(joinpath(outdir, "phase_summary.csv"), "w") do io
        println(io, "cycle_idx,phase,t_start_s,t_end_s,duration_s,V_start_V,V_end_V,capacity_Ah,terminated_by,T_max_K,T_mean_end_K,D_max,D_mean,dD_max")
        for cr in result.cycle_results
            for (label, ph) in (("discharge", cr.discharge), ("rest1", cr.rest1),
                                ("charge", cr.charge), ("rest2", cr.rest2))
                ph === nothing && continue
                println(io, join((
                    cr.cycle_idx, label, ph.t_start, ph.t_end, ph.duration,
                    ph.V_start, ph.V_end, ph.capacity, string(ph.terminated_by),
                    ph.T_max, ph.T_mean_end, ph.D_max, ph.D_mean, ph.ΔD_max), ","))
            end
        end
    end
    println("\n  cycle_summary.csv / phase_summary.csv 已导出")

    # --- 逐元素全历史导出（选定循环）：element_map + damage/sep_n/sep_t + f/stress 宽表 ---
    if !isempty(result.czm_snapshots)
        n_coh = length(result.czm_snapshots[1].damage)
        open(joinpath(outdir, "element_map.csv"), "w") do io
            println(io, "elem_id,interface,x_m,y_m")
            for (j, elem) in enumerate(case.czm_mesh.cohesive_elements)
                ns = vcat(elem.nodes_bottom, elem.nodes_top)
                x = sum(case.czm_mesh.node[ns, 1]) / length(ns) * case.param.scale.L
                y = sum(case.czm_mesh.node[ns, 2]) / length(ns) * case.param.scale.L
                println(io, "$j,$(string(elem.interface_type)),$x,$y")
            end
        end
        δ_czm = case.param.scale.δ_czm
        σ_czm = case.param.scale.σ_czm
        snap_cycles_present = sort(unique([s.cycle for s in result.czm_snapshots]))
        header = "t_s," * join(("e$j" for j in 1:n_coh), ",")
        for cyc in snap_cycles_present
            snaps = sort(filter(s -> s.cycle == cyc, result.czm_snapshots), by = s -> s.time_s)
            for (fname, getcol) in (
                    ("damage_history_cyc$cyc.csv", s -> s.damage),
                    ("sep_n_history_cyc$cyc.csv", s -> s.separation_n .* δ_czm),
                    ("sep_t_history_cyc$cyc.csv", s -> s.separation_t .* δ_czm))
                open(joinpath(outdir, fname), "w") do io
                    println(io, header)
                    for s in snaps
                        println(io, string(s.time_s) * "," * join(getcol(s), ","))
                    end
                end
            end
            @printf("  全历史已导出: 循环 %d（%d 步 × %d 单元）\n", cyc, length(snaps), n_coh)
        end

        # f_n/f_p 逐机械步宽表（任务 58 面积反馈：从 phase solve_result 聚合选定圈）
        ne_th = size(mesh_th.element, 1)
        fn_header = "t_s," * join(("e$j" for j in 1:ne_th), ",")
        for cyc in snap_cycles_present
            open(joinpath(outdir, "area_fraction_history_cyc$cyc.csv"), "w") do io
                println(io, "t_s,phase,f_n_min,f_p_min")
                for cr in result.cycle_results
                    cr.cycle_idx == cyc || continue
                    for (pname, ph) in (("discharge", cr.discharge), ("rest1", cr.rest1),
                                        ("charge", cr.charge), ("rest2", cr.rest2))
                        ph === nothing && continue
                        raw = ph.solve_result === nothing ? nothing :
                            get(ph.solve_result, "thermal2D effective area fraction n", nothing)
                        raw === nothing && continue
                        raw_p = get(ph.solve_result, "thermal2D effective area fraction p", nothing)
                        tvec = get(ph.solve_result, "time [s]", nothing)
                        for k in axes(raw, 2)
                            t = tvec === nothing ? k : ph.t_start + tvec[k] - tvec[1]
                            println(io, "$t,$pname,$(minimum(raw[:, k])),$(minimum(raw_p[:, k]))")
                        end
                    end
                end
            end
        end

        # 应力逐机械步宽表（选定圈，从 diffusion stress 在线导出键聚合）
        for cyc in snap_cycles_present
            open(joinpath(outdir, "stress_history_cyc$cyc.csv"), "w") do io
                println(io, "t_s,phase,e_id,sigma_xx_MPa,sigma_yy_MPa,sigma_xy_MPa")
                wrote_any = false
                for cr in result.cycle_results
                    cr.cycle_idx == cyc || continue
                    for (pname, ph) in (("discharge", cr.discharge), ("rest1", cr.rest1),
                                        ("charge", cr.charge), ("rest2", cr.rest2))
                        ph === nothing && continue
                        raw_x = ph.solve_result === nothing ? nothing :
                            get(ph.solve_result, "diffusion stress xx [Pa]", nothing)
                        raw_x === nothing && continue
                        raw_y = get(ph.solve_result, "diffusion stress yy [Pa]", nothing)
                        raw_xy = get(ph.solve_result, "diffusion stress xy [Pa]", nothing)
                        tvec = get(ph.solve_result, "time [s]", nothing)
                        n_mech = size(raw_x, 1)
                        for k in axes(raw_x, 2)
                            t = tvec === nothing ? k : ph.t_start + tvec[k] - tvec[1]
                            for e in 1:n_mech
                                println(io, "$t,$pname,$e,$(raw_x[e, k] * 1e-6),$(raw_y[e, k] * 1e-6),$(raw_xy[e, k] * 1e-6)")
                            end
                            wrote_any = true
                        end
                    end
                end
                wrote_any || println("  （循环 $cyc 无层分辨应力键——需在线导出激活）")
            end
        end
        println("  f_n/f_p 与应力历史已导出（选定圈）")
    end

    # --- 场数据长表导出（JUBAT_EXPORT_FIELD_DATA=1，需 save_detailed 求解） ---
    if export_field && !isempty(result.cycle_results)
        L = case.param.scale.L
        δ_czm = case.param.scale.δ_czm

        # node_temperature.csv（热节点，全相位逐时间步）
        open(joinpath(outdir, "node_temperature.csv"), "w") do f
            println(f, "cycle,phase,time_s,node_id,T_K")
            for cr in result.cycle_results, (pname, ph) in
                    (("discharge", cr.discharge), ("rest1", cr.rest1),
                     ("charge", cr.charge), ("rest2", cr.rest2))
                ph === nothing && continue
                raw = ph.solve_result === nothing ? nothing :
                    get(ph.solve_result, "thermal2D temperature at nodes [K]", nothing)
                raw === nothing && continue
                tvec = ph.solve_result === nothing ? nothing :
                    get(ph.solve_result, "time [s]", nothing)
                for k in axes(raw, 2)
                    t = tvec === nothing ? k * (ph.duration / size(raw, 2)) :
                        ph.t_start + tvec[k] - tvec[1]
                    for n in 1:size(raw, 1)
                        println(f, "$(cr.cycle_idx),$pname,$t,$n,$(raw[n, k])")
                    end
                end
            end
        end

        # node_displacement.csv（力学节点，快照时刻）
        if !isempty(result.czm_snapshots)
            open(joinpath(outdir, "node_displacement.csv"), "w") do f
                println(f, "cycle,phase,time_s,node_id,ux,uy")
                for s in result.czm_snapshots
                    nn = case.czm_mesh.nnode
                    for n in 1:nn
                        println(f, "$(s.cycle),$(s.phase),$(s.time_s),$n,$(s.displacement[2n-1]*L),$(s.displacement[2n]*L)")
                    end
                end
            end
        end

        # element_currents.csv（热单元电流）
        open(joinpath(outdir, "element_currents.csv"), "w") do f
            println(f, "cycle,phase,time_s,elem_id,I_e")
            for cr in result.cycle_results, (pname, ph) in
                    (("discharge", cr.discharge), ("rest1", cr.rest1),
                     ("charge", cr.charge), ("rest2", cr.rest2))
                ph === nothing && continue
                raw = ph.solve_result === nothing ? nothing :
                    get(ph.solve_result, "thermal2D element current", nothing)
                raw === nothing && continue
                tvec = ph.solve_result === nothing ? nothing :
                    get(ph.solve_result, "time [s]", nothing)
                for k in axes(raw, 2)
                    t = tvec === nothing ? k * (ph.duration / size(raw, 2)) :
                        ph.t_start + tvec[k] - tvec[1]
                    for e in 1:size(raw, 1)
                        println(f, "$(cr.cycle_idx),$pname,$t,$e,$(raw[e, k])")
                    end
                end
            end
        end

        # cohesive_damage.csv（plot_czm 损伤云图契约：长表）
        if !isempty(result.czm_snapshots)
            open(joinpath(outdir, "cohesive_damage.csv"), "w") do f
                println(f, "cycle,phase,time_s,coh_id,D,sep_n_m,sep_t_m")
                for s in result.czm_snapshots
                    for j in 1:length(s.damage)
                        println(f, "$(s.cycle),$(s.phase),$(s.time_s),$j,$(s.damage[j]),$(s.separation_n[j]*δ_czm),$(s.separation_t[j]*δ_czm)")
                    end
                end
            end
        end

        println("  场数据长表已导出: node_temperature / node_displacement / element_currents / cohesive_damage")
    end

    println()
    D_max_all = isempty(result.D_max) ? 0.0 : maximum(result.D_max)
    n_frac_all = isempty(result.n_fractured) ? 0 : maximum(result.n_fractured)
    @printf("  全程 D_max = %.4e ｜ 断裂单元峰值 = %d\n", D_max_all, n_frac_all)

    fm = result.final_mech === nothing ? case.mech : result.final_mech
    if fm.gp_damage_states !== nothing
        gp = fm.gp_damage_states
        @printf("  最终GP: D_max=%.4e, D_visc_max=%.4e, D>=0.99点数=%d, fractured点数=%d / %d\n",
            maximum(s.D for s in gp), maximum(s.D_visc for s in gp),
            count(s -> s.D >= 0.99, gp), count(s -> s.fractured, gp), length(gp))
    end
    if fm.plastic_states !== nothing
        kappa_max = maximum(s.kappa for s in fm.plastic_states)
        @printf("  最终等效塑性应变 KAPPA_MAX = %.4e\n", kappa_max)
    else
        println("  （调试口径 J2 关闭：无塑性状态）")
    end

    for (iface, ip) in ((:PE_PCC, case.param.PCC), (:NE_NCC, case.param.NCC))
        idx = findall(e -> e.interface_type == iface, case.czm_mesh.cohesive_elements)
        states = fm.damage_states[idx]
        @printf("  %s 最终: D>0单元数=%d, D>=0.99单元数=%d, D_mean=%.4e\n",
            string(iface), count(s -> s.D > 1e-8, states),
            count(s -> s.D >= 0.99, states),
            sum(s.D for s in states) / length(states))
    end

    println()
    if D_max_all > 0.0
        @printf("  结论：损伤已启动（D_max = %.4e > 0）\n", D_max_all)
    else
        @printf("  结论：已完成的 %d 次循环 D_max 恒为 0\n", length(result.cycle_idx))
    end

    println("\n" * "="^80)
    println("全部完成")
    println("="^80)
end

main()
