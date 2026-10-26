"""
Jellyroll 电池极耳多工况研究：nθ=20 × 1 次循环 × 9 工况极耳设计矩阵

仿照 example/czm_cycling.jl 的耦合调试口径（多 SPMe 并行 + 分布式热 + CZM，
gp_basic + 粘性 τ=5 s，geo/J2 关闭），在 nθ=20 粗网格上按极耳设计矩阵逐工况
跑 1 个完整循环（放电 1800 s → 静置 600 s → 充电 1800 s → 静置 600 s），
研究极耳位置/数量/尺寸对温度场与界面损伤的影响。

════════════════════════════════════════════════════════════════════════
极耳工况设计矩阵（9 工况，独立 7 个：Case 4 ≡ 2 ≡ 8 基准重复）
（W 基准 4 mm、L 基准 0.75×9.9 mm=7.425 mm、面积 W×L=29.7 mm²；rng=螺旋累计角范围）

  组别     Case  n_pos  θ_pos                        W[mm]  L[mm]    说明
  位置     1     1      10%·rng（螺旋起始，内层）      4.0    7.425   单耳位置扫描
           2     1      50%·rng（螺旋中间，基准）      4.0    7.425   ≡ Case 4/8
           3     1      90%·rng（螺旋末端，外层）      4.0    7.425
  数量     4     1      50%·rng                       4.0    7.425   ≡ Case 2
           5     3      50%·rng ± 25%·rng, 50%·rng    4.0    L/3     等间距覆盖 rng/2
           6     5      50%·rng + k·rng/8 (k=−2..2)   4.0    L/5     等间距覆盖 rng/2
  尺寸     7     1      50%·rng                       2.0    14.85   细长（L/W=7.43）
           8     1      50%·rng                       4.0    7.425   ≡ Case 2
           9     1      50%·rng                       8.0    3.7125  粗短（L/W=0.46）

变量控制：位置组固定 n/W/L 只变 θ；数量组固定总长 L_total=W×L 与 W、等间距覆盖
rng/2 只变 n（Δθ=rng/8）；尺寸组固定 n/θ/面积 29.7 mm² 只变长宽比。
所有工况 tab.area = 2·W·L·n_pos ≡ 59.4 mm² 恒定 → 总极耳冷却面积不变，
研究的是冷却/产热的"分布"而非总量。

════════════════════════════════════════════════════════════════════════
物理激活路径（src/ThermalDistributed.jl apply_cool_method "tab" 分支）：
  - opt.cool_method = "tab"（极耳体冷却替代表面冷却——极耳是本工况唯一的
    位置/数量/尺寸敏感通道，"surface" 分支完全不消费极耳参数）；
  - 极耳电阻产热（任务 50）：sig_pos=3.55e7（Al）/ sig_neg=5.96e7（镀镍铜）/
    thickness=165e-6 m 三者>0 启用；同极性 N 耳并联均分，每耳 (I/N)²·R，
    R = L/(sig·W·t) 逐耳注入自己的 θ 窗口节点组；
  - tab.h 默认匹配表面冷却总容量：h_tab = cell.h·cooling_surface/tab.area
    ≈ 800 W/m²K（保证与 surface 基线可比的热平衡温度；纯 tab 冷却 h=10 时
    冷却导纳仅 ~6e-4 W/K，1C 产热下热失控，不可用）。

已知简化（单 tab 结构限制）：tab.width/length 为正负极共享字段，数量/尺寸
工况中负极耳几何随正极耳联动（负耳 I²R 为二阶小量，~8 mW）；theta_neg 恒为
默认 [44π]（最外节点线），单因素设计只动正极耳。

界面参数直接以物理量赋给 param_dim.PCC（PE-PCC）和 param_dim.NCC（NE-NCC）：
σ_max/τ_max 为 6.97/7.82 MPa，K_n/K_t 为 2.4e16/1.2e16 Pa/m，
G_c/G_c_t 为 50.6/12.4 J/m²；沿用原脚本默认数值，分离量及归一化尺度由主参数计算。

环境变量：
  JUBAT_TAB_CASES=<csv>      运行的工况号（默认 "1,2,3,5,6,7,9" 独立集；"all"=1:9）
  JUBAT_TAB_H=<f>            极耳等效换热系数（默认按表面冷却容量匹配 ≈800）
  JUBAT_CZM_METHOD=<name>    求解方法（默认 gp_basic；basic 在此参数下不收敛）
  JUBAT_CZM_TAU_SECONDS=<s>  物理时间粘性 τ（默认 5；设为空串关闭粘性）
  JUBAT_RUN_TAG=<name>       输出子目录名（默认 tabstudy_a0.82_n20）

输出（output/czm_cycling_n20/<tag>/）：tab_case_matrix.csv（跨工况对比总表）+
case<k>/ 子目录（cycle_summary / phase_summary / element_map / 循环 1 全历史
damage / sep_n / sep_t / traction_n / traction_t / 云图契约 mesh_nodes /
mesh_elements / tab_nodes / cloud_damage / cloud_temperature / cloud_current /
case_manifest：损伤最大时刻瞬态 + 损伤起始时刻/位置，供
postprocess/plot_tab_study_clouds.py 绘图）。

日期：2026-10-12
"""

using Printf
# CSV 时间使用相位起点加真实局部时间；相位身份来自求解器，保留重复时间记录。
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
include(joinpath(@__DIR__, "cycling_csv_time.jl"))

# ── 工况矩阵定义（W/L 单位 m；theta_kind 在求值期以 rng 解释） ──
struct TabCase
    id::Int
    group::String
    n_pos::Int
    W::Float64
    L::Float64          # 每耳长度（数量组 L_total/n）
    theta_kind::Symbol  # :t10 :t50 :t90 :n3 :n5
    dup_of::Union{Int, Nothing}
    note::String
end

const L_BASE = 0.75 * 9.9e-3     # 7.425 mm 基准总长
const TAB_CASES = [
    TabCase(1, "位置", 1, 4e-3, L_BASE,       :t10, nothing, "单耳@螺旋起始（内层10%）"),
    TabCase(2, "位置", 1, 4e-3, L_BASE,       :t50, nothing, "单耳@螺旋中间（基准）"),
    TabCase(3, "位置", 1, 4e-3, L_BASE,       :t90, nothing, "单耳@螺旋末端（外层90%）"),
    TabCase(4, "数量", 1, 4e-3, L_BASE,       :t50, 2,       "基准（≡Case 2）"),
    TabCase(5, "数量", 3, 4e-3, L_BASE / 3,   :n3,  nothing, "三耳等间距覆盖rng/2"),
    TabCase(6, "数量", 5, 4e-3, L_BASE / 5,   :n5,  nothing, "五耳等间距覆盖rng/2（Δθ=rng/8）"),
    TabCase(7, "尺寸", 1, 2e-3, 29.7e-6 / 2e-3,  :t50, nothing, "细长 L/W=7.43（面积恒定）"),
    TabCase(8, "尺寸", 1, 4e-3, L_BASE,       :t50, 2,       "基准形（≡Case 2）"),
    TabCase(9, "尺寸", 1, 8e-3, 29.7e-6 / 8e-3,  :t50, nothing, "粗短 L/W=0.46（面积恒定）"),
]

theta_pos_for(kind::Symbol, rng::Float64) = kind == :t10 ? [0.10rng] :
    kind == :t50 ? [0.50rng] :
    kind == :t90 ? [0.90rng] :
    kind == :n3 ? [0.50rng - 0.25rng, 0.50rng, 0.50rng + 0.25rng] :
    kind == :n5 ? [0.50rng + k * rng / 8 for k in -2:2] :
    error("未知 theta_kind=$kind")

# 螺旋累计角范围（无量纲）：外表面最外节点的 theta_cum = (Rout−a)/b 上界，
# 与网格生成器同参数式；θ_pos 按其分数取值
function spiral_theta_range(param_dim)
    b = param_dim.cell.layer / (2π)
    return (param_dim.cell.Rout - param_dim.cell.Rin) / b
end

function apply_tab_design!(param_dim, tc::TabCase, thetas::Vector{Float64}, h_tab::Float64)
    t = param_dim.tab
    t.theta_pos = thetas
    t.width = tc.W
    t.length = tc.L                       # 每耳长度（R 逐耳用）
    t.area = 2.0 * tc.W * tc.L * tc.n_pos # 总极耳面积恒定 59.4 mm²（正负合计口径）
    t.h = h_tab
    # 任务 50 极耳电阻产热参考值（Al / 镀镍铜 / 165 µm）
    t.sig_pos = 3.55e7
    t.sig_neg = 5.96e7
    t.thickness = 165e-6
end

function main()
    iter_method = get(ENV, "JUBAT_CZM_METHOD", "gp_basic")
    tau_seconds = get(ENV, "JUBAT_CZM_TAU_SECONDS", "5")

    cases_raw = get(ENV, "JUBAT_TAB_CASES", "1,2,3,5,6,7,9")
    case_ids = lowercase(strip(cases_raw)) == "all" ? collect(1:9) :
        [parse(Int, strip(x)) for x in split(cases_raw, ",")]
    all(id -> 1 <= id <= 9, case_ids) || error("JUBAT_TAB_CASES 工况号须在 1..9: $cases_raw")

    run_tag = get(ENV, "JUBAT_RUN_TAG", "tabstudy_a0.85_n20")
    outdir = joinpath(@__DIR__, "..", "output", "czm_cycling_n20", run_tag)
    mkpath(outdir)

    n_theta = 20
    n_cycles = 1

    println("="^80)
    @printf("Jellyroll 极耳多工况研究：nθ=%d / %d 次循环 / %d 个工况 / SOC 0.65\n",
        n_theta, n_cycles, length(case_ids))
    println("="^80)

    # 共享几何：所有工况只差 tab 字段，热/CZM 网格完全相同，构建一次复用
    probe_dim = JuBat.ChooseCell("Jellyroll")
    rng = spiral_theta_range(probe_dim)
    tab_area_ref = 2.0 * 4e-3 * L_BASE            # 59.4 mm²（全部工况恒定）
    h_tab = get(ENV, "JUBAT_TAB_H", "") |> s -> isempty(s) ?
        probe_dim.cell.h * probe_dim.cell.cooling_surface / tab_area_ref :
        parse(Float64, s)

    @printf("  螺旋累计角范围 rng = %.2f rad = %.2fπ（θ_pos 按分数取值）\n", rng, rng / π)
    @printf("  极耳总面积恒定 %.2f mm² → 匹配 tab.h = %.1f W/m²K（表面冷却容量等价）\n",
        tab_area_ref * 1e6, h_tab)
    @printf("  CZM: %s / τ=%s s\n", iter_method,
        isempty(tau_seconds) ? "关" : tau_seconds)
    println("  工况: " * join(case_ids, ","))

    # 单循环全历史采集（对每个 solve_cycling 调用生效）
    ENV["JUBAT_SNAPSHOT_CYCLES"] = "1"

    mesh_data = nothing
    summaries = TabCaseSummary[]
    for id in case_ids
        tc = TAB_CASES[id]
        thetas = theta_pos_for(tc.theta_kind, rng)
        summary, mesh_data = run_tab_case(tc, thetas, h_tab, outdir, mesh_data;
            n_theta, n_cycles,
            iter_method, tau_seconds)
        push!(summaries, summary)
    end

    # ── 跨工况对比总表 ──
    println("\n" * "="^80)
    println("极耳工况对比总表")
    println("="^80)
    @printf("  %-4s %-4s %-5s %-22s %-6s %-7s %-6s %-10s %-10s %-8s %-8s\n",
        "Case", "组别", "n", "θ_pos/π", "W[mm]", "L[mm]", "tab点", "T_max[K]", "D_max", "PCC损", "状态")
    for s in summaries
        @printf("  %-4d %-4s %-5d %-22s %-6.1f %-7.4f %-6s %-10.2f %-10.3e %-8d %-8s\n",
            s.id, s.group, s.n_pos, s.theta_str, s.W_mm, s.L_mm, s.tab_nodes,
            s.T_max, s.D_max, s.pcc_dmg, s.status)
    end

    open(joinpath(outdir, "tab_case_matrix.csv"), "w") do io
        println(io, "case,group,n_pos,theta_pos_over_pi,W_mm,L_mm,pos_tab_nodes,neg_tab_nodes,",
            "T_max_K,T_mean_end_K,D_max,D_mean,pcc_D_gt0,ncc_D_gt0,n_fractured,",
            "capacity_discharge_Ah,status,wall_s,note")
        for s in summaries
            println(io, join((
                s.id, s.group, s.n_pos, s.theta_str, s.W_mm, s.L_mm,
                s.pos_tab_nodes, s.neg_tab_nodes,
                s.T_max, s.T_mean_end, s.D_max, s.D_mean, s.pcc_dmg, s.ncc_dmg,
                s.n_fractured, s.capacity_Ah, s.status, s.wall_s,
                replace(s.note, "," => "，")), ","))
        end
    end
    println("\n  tab_case_matrix.csv 已导出")
    println("\n" * "="^80)
    println("全部完成")
    println("="^80)
end

struct TabCaseSummary
    id::Int; group::String; n_pos::Int; theta_str::String
    W_mm::Float64; L_mm::Float64; tab_nodes::String
    pos_tab_nodes::Int; neg_tab_nodes::Int
    T_max::Float64; T_mean_end::Float64
    D_max::Float64; D_mean::Float64; pcc_dmg::Int; ncc_dmg::Int; n_fractured::Int
    capacity_Ah::Float64; status::String; wall_s::Float64; note::String
end

function build_tab_case(tc::TabCase, thetas::Vector{Float64}, h_tab::Float64,
        outdir::String, mesh_data; n_theta, n_cycles, iter_method, tau_seconds,
        phase_times=(discharge=1800.0, rest1=600.0, charge=1800.0, rest2=600.0), dt_cycle=[1.0, 10.0])
    casedir = joinpath(outdir, "case$(tc.id)")
    mkpath(casedir)
    println("\n" * "-"^78)
    @printf("[Case %d / %s] n_pos=%d, θ_pos/π=[%s], W=%.1fmm, L=%.4fmm（%s）\n",
        tc.id, tc.group, tc.n_pos, join(round.(thetas ./ π, digits=2), ", "),
        tc.W * 1e3, tc.L * 1e3, tc.note)

    # 1. 参数：电化学/热/CZM 同调试口径 + 极耳工况设计
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    # PE-PCC 界面：强度 [Pa]、刚度 [Pa/m]、断裂能 [J/m²]
    param_dim.PCC.σ_max = 6.97e6
    param_dim.PCC.τ_max = 6.97e6
    param_dim.PCC.K_n = 2.4e16
    param_dim.PCC.K_t = 2.4e16
    param_dim.PCC.G_c = 50.6
    param_dim.PCC.G_c_t = 50.6
    # NE-NCC 界面：强度 [Pa]、刚度 [Pa/m]、断裂能 [J/m²]
    param_dim.NCC.σ_max = 7.82e6
    param_dim.NCC.τ_max = 7.82e6
    param_dim.NCC.K_n = 1.2e16
    param_dim.NCC.K_t = 1.2e16
    param_dim.NCC.G_c = 12.4
    param_dim.NCC.G_c_t = 12.4
    for ip in (param_dim.PCC, param_dim.NCC)
        ip.δ_0 = ip.σ_max / ip.K_n
        ip.δ_0_t = ip.τ_max / ip.K_t
        ip.δ_c = 2.0 * ip.G_c / ip.σ_max
        ip.δ_c_t = 2.0 * ip.G_c_t / ip.τ_max
    end
    param_dim.scale.σ_czm = param_dim.PCC.σ_max
    param_dim.scale.δ_czm = 2.0 * param_dim.PCC.G_c / param_dim.PCC.σ_max
    param_dim.scale.G_czm = param_dim.scale.σ_czm * param_dim.scale.δ_czm
    param_dim.scale.K_czm = param_dim.scale.σ_czm / param_dim.scale.δ_czm
    for (label, ip) in (("PE-PCC", param_dim.PCC), ("NE-NCC", param_dim.NCC))
        @printf("  %s: σ_max/τ_max=%.2f/%.2f MPa, K_n/K_t=%.1e/%.1e Pa/m, G_c/G_c_t=%.1f/%.1f J/m²\n",
            label, ip.σ_max * 1e-6, ip.τ_max * 1e-6, ip.K_n, ip.K_t, ip.G_c, ip.G_c_t)
    end
    apply_tab_design!(param_dim, tc, thetas, h_tab)

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
    opt.cool_method = "tab"          # 极耳体冷却 + I²R 产热（本研究的物理通道）
    opt.per_element_spme = true

    opt.debug_coupling = false
    opt.debug_log_path = joinpath(casedir, "simple_coupling_debug.log")
    opt.czm.enabled = true
    opt.czm.model = "mix"
    opt.czm.fix_inner = false
    opt.czm.iter_method = iter_method
    opt.czm.max_iter = parse(Int, get(ENV, "JUBAT_CZM_MAX_ITER", string(opt.czm.max_iter)))
    opt.czm.load_steps = 10
    opt.czm.tol = 1e-3
    opt.czm.geo_nonlinear = false
    opt.czm.j2_plasticity = false
    opt.czm.update_interval = 1
    if !isempty(tau_seconds)
        opt.czm.viscous_enabled = true
        opt.czm.viscous_tau = parse(Float64, tau_seconds)
    end

    cycle_opt = JuBat.CycleOption(
        n_cycles = n_cycles,
        t_charge = phase_times.charge,
        t_discharge = phase_times.discharge,
        t_rest1 = phase_times.rest1,
        t_rest2 = phase_times.rest2,
        SOC_init = 0.65,
        dt_cycle = dt_cycle,
    )

    # 2. 案例 + 网格（几何与极耳无关 → 首工况构建后复用）
    case = JuBat.SetCase(param_dim, opt)
    if mesh_data === nothing
        mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, czm_enabled=true, gsorder=2)
    end
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    case.mesh["czm"] = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_th)
    mechanical_seed=JuBat.MechState(case.mesh["czm"])
    return (; case, mesh_data, mechanical_seed, cycle_opt, casedir)
end

function run_tab_case(tc::TabCase, thetas::Vector{Float64}, h_tab::Float64,
        outdir::String, mesh_data; n_theta, n_cycles, iter_method, tau_seconds,
        phase_times=(discharge=1800.0, rest1=600.0, charge=1800.0, rest2=600.0), dt_cycle=[1.0, 10.0])
    (; case, mesh_data, mechanical_seed, cycle_opt, casedir) =
        build_tab_case(tc, thetas, h_tab, outdir, mesh_data;
            n_theta, n_cycles,
            iter_method, tau_seconds, phase_times, dt_cycle)
    mesh_th = case.mesh["thermal2D"]
    runtime=Dict{String,Any}()


    # 极耳窗口节点组核对（空组=极耳物理静默，必须在结果中可见）
    pos_groups, neg_groups = JuBat.jellyroll_tab_node_groups(mesh_th, case.param)
    n_pos_nodes = sum(length, pos_groups; init=0)
    n_neg_nodes = sum(length, neg_groups; init=0)
    @printf("  极耳窗口节点: pos %s / neg %s（空组=该耳静默）\n",
        join(length.(pos_groups), "|"), join(length.(neg_groups), "|"))
    (n_pos_nodes > 0 && n_neg_nodes > 0) ||
        @warn "Case $(tc.id) 极耳窗口存在空节点组" pos_groups neg_groups

    # 3. 求解（失败不中断矩阵，状态记录后继续）
    t0 = time_ns()
    local result = nothing
    status = "OK"
    err_msg = ""
    try
        result = JuBat.solve_cycling(case, cycle_opt, mechanical_seed; variables=runtime, save_detailed = true)
    catch err
        status = "失败"
        err_msg = first(sprint(showerror, err), 160)
        println("  [诊断] Case $(tc.id) solve_cycling 失败: $err_msg")
        mech = runtime["accepted"].mech
        if mech.gp_damage_states !== nothing
            gp = mech.gp_damage_states
            @printf("  GP: D_max=%.4f, D>=0.99=%d / %d\n",
                maximum(s.D for s in gp), count(s -> s.D >= 0.99, gp), length(gp))
        end
    end
    wall_s = (time_ns() - t0) * 1e-9

    # 4. 结果汇总 + 导出
    if status == "OK"
        k = length(result.cycle_idx)
        fm = result.final_mech
        idx_pcc = findall(e -> e.interface_type == :PE_PCC, case.mesh["czm"].cohesive_elements)
        idx_ncc = findall(e -> e.interface_type == :NE_NCC, case.mesh["czm"].cohesive_elements)
        pcc_dmg = count(s -> s.D > 1e-8, fm.damage_states[idx_pcc])
        ncc_dmg = count(s -> s.D > 1e-8, fm.damage_states[idx_ncc])
        @printf("  OK: 墙钟 %.1f s ｜ T_max=%.2f K ｜ D_max=%.3e（PCC %d / NCC %d 损伤）｜ 放电 %.3f Ah\n",
            wall_s, result.T_max[k], result.D_max[k], pcc_dmg, ncc_dmg, result.capacity_discharge[k])

        open(joinpath(casedir, "cycle_summary.csv"), "w") do io
            println(io, "cycle_idx,capacity_charge_Ah,capacity_discharge_Ah,coulombic_efficiency,D_max,D_mean,n_fractured,T_max_K,soh")
            for c in eachindex(result.cycle_idx)
                println(io, join((
                    result.cycle_idx[c], result.capacity_charge[c], result.capacity_discharge[c],
                    result.coulombic_efficiency[c], result.D_max[c], result.D_mean[c],
                    result.n_fractured[c], result.T_max[c], result.soh[c]), ","))
            end
        end
        open(joinpath(casedir, "phase_summary.csv"), "w") do io
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
        snapshot_rows = write_czm_snapshot_csv(result,case,casedir)
        if !isempty(result.czm_snapshots)
            snaps = [row.snapshot for row in snapshot_rows]
            t_abs = [row.time_s for row in snapshot_rows]
            phase_of_snap(i) = snapshot_rows[i].phase
        end
        # —— 云图数据导出（Python 后处理契约，仿 run_boundary_current_clouds.jl）：
        #    损伤最大时刻瞬态（损伤/温度/电流）+ 网格 + 极耳窗口节点 + 起始信息 manifest
        if !isempty(result.czm_snapshots)
            L_scale = case.param.scale.L
            open(joinpath(casedir, "mesh_nodes.csv"), "w") do io
                println(io, "node_id,x_m,y_m")
                for n in 1:size(mesh_th.node, 1)
                    println(io, "$n,$(mesh_th.node[n, 1] * L_scale),$(mesh_th.node[n, 2] * L_scale)")
                end
            end
            open(joinpath(casedir, "mesh_elements.csv"), "w") do io
                println(io, "elem_id,n1,n2,n3,n4")
                for e in 1:size(mesh_th.element, 1)
                    el = mesh_th.element[e, :]
                    println(io, "$e,$(el[1]),$(el[2]),$(el[3]),$(el[4])")
                end
            end
            open(joinpath(casedir, "tab_nodes.csv"), "w") do io
                println(io, "polarity,group_idx,node_id")
                for (g, grp) in enumerate(pos_groups), n in grp
                    println(io, "pos,$g,$n")
                end
                for (g, grp) in enumerate(neg_groups), n in grp
                    println(io, "neg,$g,$n")
                end
            end
            # 选定时刻 = 全局损伤最大首次达到（≥0.999×终值，按绝对时间顺序）
            Dmat = reduce(hcat, [s.damage for s in snaps])
            gmax = vec(maximum(Dmat, dims = 1))
            i_sel = findfirst(>=(0.999 * maximum(gmax)), gmax)
            t_sel = t_abs[i_sel]
            phase_sel = phase_of_snap(i_sel)
            D_sel = Dmat[:, i_sel]
            # 损伤起始 = 首个 max D>0.01 的快照；该时刻最深单元为起始位置
            i_init = findfirst(>(0.01), gmax)
            t_init = i_init === nothing ? NaN : t_abs[i_init]
            e_init = i_init === nothing ? 0 : argmax(Dmat[:, i_init])
            init_xy = if e_init == 0
                (NaN, NaN)
            else
                ns = vcat(case.mesh["czm"].cohesive_elements[e_init].nodes_bottom,
                          case.mesh["czm"].cohesive_elements[e_init].nodes_top)
                (sum(case.mesh["czm"].node[ns, 1]) / length(ns) * L_scale,
                 sum(case.mesh["czm"].node[ns, 2]) / length(ns) * L_scale)
            end
            # Sample only the selected snapshot's real cycle/phase, including at phase seams.
            selected = snaps[i_sel]
            ph = cycling_csv_phases(result)[(selected.cycle, selected.phase)]
            sr = ph.solve_result
            raw_T = sr["thermal2D temperature at nodes [K]"]
            raw_I = sr["thermal2D element current"]
            size(raw_T, 2) == size(raw_I, 2) || throw(DimensionMismatch("cloud field times differ"))
            tvec = cycling_csv_times(ph, size(raw_T, 2))
            k_sel = argmin(abs.(tvec .- t_sel))
            node_T = raw_T[:, k_sel]
            elem_I = raw_I[:, k_sel]
            open(joinpath(casedir, "cloud_damage.csv"), "w") do io
                println(io, "coh_id,D")
                for j in eachindex(D_sel)
                    println(io, "$j,$(D_sel[j])")
                end
            end
            open(joinpath(casedir, "cloud_temperature.csv"), "w") do io
                println(io, "node_id,T_K")
                for n in eachindex(node_T)
                    println(io, "$n,$(node_T[n])")
                end
            end
            open(joinpath(casedir, "cloud_current.csv"), "w") do io
                println(io, "elem_id,I_e")
                for e in eachindex(elem_I)
                    println(io, "$e,$(elem_I[e])")
                end
            end
            open(joinpath(casedir, "case_manifest.csv"), "w") do io
                println(io, "key,value")
                for (mk, mv) in (("case_id", tc.id), ("group", tc.group), ("n_pos", tc.n_pos),
                                 ("theta_pos_over_pi", join(round.(thetas ./ π, digits = 3), ";")),
                                 ("W_mm", tc.W * 1e3), ("L_mm", tc.L * 1e3),
                                 ("t_sel_s", t_sel), ("phase_sel", phase_sel),
                                 ("D_max", maximum(gmax)),
                                 ("t_init_s", t_init),
                                 ("phase_init", i_init === nothing ? "NA" : phase_of_snap(i_init)),
                                 ("init_coh_id", e_init), ("init_x_m", init_xy[1]),
                                 ("init_y_m", init_xy[2]),
                                 ("T_max_K", result.T_max[k]), ("status", "OK"))
                    println(io, "$mk,$mv")
                end
            end
            @printf("  云图导出: t_sel=%.1f s（%s, 损伤最大）｜ t_init=%.1f s @ coh#%d (%.2f, %.2f) mm\n",
                t_sel, phase_sel, t_init, e_init, init_xy[1] * 1e3, init_xy[2] * 1e3)
        end
        last_phase = result.cycle_results[end].rest2
        summary = TabCaseSummary(tc.id, tc.group, tc.n_pos,
            join(round.(thetas ./ π, digits=3), ";"),
            tc.W * 1e3, tc.L * 1e3, "$(n_pos_nodes)/$(n_neg_nodes)",
            n_pos_nodes, n_neg_nodes,
            result.T_max[k], last_phase === nothing ? NaN : last_phase.T_mean_end,
            result.D_max[k], result.D_mean[k], pcc_dmg, ncc_dmg,
            result.n_fractured[k], result.capacity_discharge[k], status, wall_s,
            status == "OK" ? tc.note : err_msg)
    else
        summary = TabCaseSummary(tc.id, tc.group, tc.n_pos,
            join(round.(thetas ./ π, digits=3), ";"),
            tc.W * 1e3, tc.L * 1e3, "$(n_pos_nodes)/$(n_neg_nodes)",
            n_pos_nodes, n_neg_nodes,
            NaN, NaN, NaN, NaN, 0, 0, 0, NaN, status, wall_s, err_msg)
    end
    return (; summary, mesh_data, result, case, casedir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
