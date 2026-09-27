"""
SOC 0.65 五次短循环：Jellyroll电池多SPMe并行电化学-热-CZM循环仿真（纯文字结果）

以 `example/testexample_soc065_1800s.jl` 为基，改用 `solve_cycling` 跑 5 次短循环
（放电 1800 s → 静置 600 s → 充电 1800 s → 静置 600 s），回答一个问题：
**循环载荷下 CZM 损伤是否启动**（单次 1800 s 放电全程 D=0，但 κ 已累积至 6.5e-3）。

与 1800 s 单次工况的差异：
- 循环求解（solve_cycling，循环序=放电→静置→充电→静置），n_cycles=5
- 充放电极幅同 5 A；截止 4.2/2.5 V；SOC_init=0.65（由 solve_cycling 内部施加）
- reset_T_each_cycle=false：保留跨循环热累积（更利于损伤发展的保守选择）
- 其余（mix、geo+J2、basic、tol=1e-3、nθ=80、表面冷却）与 1800 s 工况一致

2026-09-24 弧长修复复测：在当前生产参数基础上，界面强度乘 0.09、
刚度乘 0.1、断裂能乘 2，得到 PE/NE 强度 7.38/8.28 MPa、
刚度 2.4/1.2e16 Pa/m、断裂能 50.6/12.4 J/m²；
geo/J2 成对关闭（D-B3-1），回答：**损伤已可启动的参数下，D 是否随循环增长**。
输出 `output/testexample_soc065_5cycles/debug_a0.9_gc2/`（不覆盖任务 51 证据）。
日期：2026-09-14
"""

using Printf
using Profile
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    run_tag = get(ENV, "JUBAT_RUN_TAG", "debug_a0.9_gc2")
    outdir = joinpath(@__DIR__, "..", "output", "testexample_soc065_5cycles", run_tag)
    mkpath(outdir)

    println("="^80)
    println("Jellyroll电池多SPMe并行电化学-热耦合循环仿真（SOC 0.65 / 5次短循环 / 调试参数）")
    println("="^80)

    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/4] 参数设置...")

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    # 当前生产默认值上的同参复测：强度×0.09、刚度×0.1、断裂能×2。
    # ChooseCell 后缩放并同步派生量与 CZM 锚（与 testexample_soc065_1800s.jl 扫描口径一致）。
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
    opt.czm.geo_nonlinear = false   # 调试口径：geo/J2 成对关闭（D-B3-1）
    opt.czm.j2_plasticity = false
    opt.czm.update_interval = 1
    tau_seconds = get(ENV, "JUBAT_CZM_TAU_SECONDS", "")
    if !isempty(tau_seconds)
        opt.czm.viscous_enabled = true
        opt.czm.viscous_tau = parse(Float64, tau_seconds)
    end
    # 2026-09-22 诊断口径：关闭粘性（回到原始崩溃配置），失败时转储 mech 状态

    cycle_opt = JuBat.CycleOption(
        n_cycles = parse(Int, get(ENV, "JUBAT_CYCLES", "5")),
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
    @printf("  循环: %d 次（放电 1800 s → 静置 600 s → 充电 1800 s → 静置 600 s）\n", cycle_opt.n_cycles)
    @printf("  初始SOC: 0.65（solve_cycling 内部施加）\n")
    @printf("  跨循环温度: 累积（reset_T_each_cycle=false）\n")
    @printf("  充电前温度: 重置（reset_T_before_charge=true）\n")
    @printf("  CZM: %s / geo=%s / J2=%s / %s / max_iter=%d / tol=%.1e / τ=%g s\n",
        opt.czm.model, opt.czm.geo_nonlinear, opt.czm.j2_plasticity,
        opt.czm.iter_method, opt.czm.max_iter, opt.czm.tol, opt.czm.viscous_tau)

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

    # ========================================================================
    # 3. 循环求解
    # ========================================================================
    println("\n[3/4] 运行循环求解器...")

    t_wall_start = time_ns()
    profile_on = get(ENV, "JUBAT_PROFILE", "0") == "1"
    local result
    try
        if profile_on
            Profile.init(n = 10_000_000, delay = 0.01)
            result = Profile.@profile JuBat.solve_cycling(case, cycle_opt, case.mech)
        else
            result = JuBat.solve_cycling(case, cycle_opt, case.mech)
        end
    catch err
        # 阶段0诊断：失败点原位状态转储（fractured 置位？受损单元承压 or 受拉？剪切主导？）
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
            n_open = n_closed = n_shear_dom = 0
            max_dn = 0.0
            for (k, j) in enumerate(idx)
                states[k].D > 1e-8 || continue
                elem = case.czm_mesh.cohesive_elements[j]
                _, _, _, R = JuBat.cohesive_local_frame(case.czm_mesh, elem)
                n1, n2 = elem.nodes_bottom
                n4, n3 = elem.nodes_top
                dx = 0.5 * (u[2*n4-1] - u[2*n1-1]) + 0.5 * (u[2*n3-1] - u[2*n2-1])
                dy = 0.5 * (u[2*n4] - u[2*n1]) + 0.5 * (u[2*n3] - u[2*n2])
                δn = R[1, 1] * dx + R[1, 2] * dy
                δt = R[2, 1] * dx + R[2, 2] * dy
                δn >= 0 ? (n_open += 1) : (n_closed += 1)
                abs(δt) > abs(δn) && (n_shear_dom += 1)
                max_dn = max(max_dn, abs(δn))
            end
            @printf("  %s 受损单元(D>0): δn>=0张开=%d, δn<0承压=%d, 剪切主导=%d, max|δn|=%.3e(归一)\n",
                string(iface), n_open, n_closed, n_shear_dom, max_dn)
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
    # 4. 损伤是否启动：逐循环与阶段分解
    # ========================================================================
    println("\n[4/4] 损伤启动分析")
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

    # CSV 导出（后续后处理）：逐循环汇总 + 逐相位明细，写入 run 目录
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
    println("\n  CSV 已导出: $(joinpath(outdir, "cycle_summary.csv")) / phase_summary.csv")

    # 逐元素全历史导出（JUBAT_SNAPSHOT_CYCLES 选定循环）：element_map + 每圈
    # damage/sep_n/sep_t 宽表（行=机械步，列=cohesive 单元，id 对应 element_map 行序）
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
        println("  单元坐标系: $(joinpath(outdir, "element_map.csv"))")
    end

    println()
    D_max_all = isempty(result.D_max) ? 0.0 : maximum(result.D_max)
    n_frac_all = isempty(result.n_fractured) ? 0 : maximum(result.n_fractured)
    @printf("  全程 D_max = %.4e ｜ 断裂单元峰值 = %d\n", D_max_all, n_frac_all)

    # solve_cycling 内部经 update_czm_damage! 演化 case.mech（ms 参数仅控制快照）
    fm = result.final_mech === nothing ? case.mech : result.final_mech
    if fm.gp_damage_states !== nothing
        gp = fm.gp_damage_states
        @printf("  最终GP: D_max=%.4e, D_visc_max=%.4e, D>=0.99点数=%d, fractured点数=%d / %d\n",
            maximum(s.D for s in gp), maximum(s.D_visc for s in gp),
            count(s -> s.D >= 0.99, gp), count(s -> s.fractured, gp), length(gp))
    end
    if fm.plastic_states !== nothing
        kappa_max = maximum(s.kappa for s in fm.plastic_states)
        n_yielded = count(s.kappa > 0 for s in fm.plastic_states)
        @printf("  最终等效塑性应变 KAPPA_MAX = %.4e（对照：单次 1800 s 放电为 6.4801e-3）\n", kappa_max)
        @printf("  屈服过的高斯点数 YIELDED = %d / %d\n", n_yielded, length(fm.plastic_states))
    else
        println("  （调试口径 J2 关闭：无塑性状态）")
    end

    # 分界面最终损伤统计（调试参数下双界面顺序启动的对照）
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
