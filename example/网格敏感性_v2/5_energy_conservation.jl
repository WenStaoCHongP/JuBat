"""
Script 5: 全场耦合能量守恒检查

在全耦合算例（SPMe + distributed2D thermal + CZM）上验证：
  Q_generated = ΔE_thermal + Q_loss + ΔE_elastic + E_fracture

输出：R(t)、ε_R(t) 图、归一化 RMS 残余
"""

using Printf, Plots, LinearAlgebra, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

function main()
    println("=" ^ 70)
    println("Script 5: 全场耦合能量守恒检查")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    nθ = 80

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    println("\n[1] 创建案例 (nθ=$nθ, 全耦合)...")
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, nθ_czm=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_data.thermal2D, case.param)
    case.czm_mesh = czm_mesh
    n_coh = czm_mesh.n_cohesive

    @printf("  热单元: %d,  节点: %d,  CZM单元: %d\n", ne, mesh_th.nlen, n_coh)

    println("\n[2] 运行全耦合求解...")
    t0 = time_ns()
    result = JuBat.Solve(case)
    dt_wall = (time_ns() - t0) * 1e-9
    @printf("  耗时 %.2f s\n", dt_wall)

    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I_arr = result["cell current [A]"]
    nt = length(t)

    println("\n[3] 计算能量平衡...")

    # 电功 W_elec = ∫ V·I dt
    W_elec = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        W_elec[k] = W_elec[k-1] + 0.5 * (V[k]*I_arr[k] + V[k-1]*I_arr[k-1]) * dt_k
    end
    @printf("  W_elec(终) = %.4e J\n", W_elec[end])

    # 热能变化 ΔE_th
    scale = param_dim.scale
    ρ  = param_dim.cell.rho
    cp = param_dim.cell.heat_Q
    Vol = param_dim.cell.volume
    T_mean = result["temperature [K]"]
    T0 = T_mean[1]
    C_total = ρ * cp * Vol

    ΔE_th = C_total .* (T_mean .- T0)
    @printf("  ΔE_th(终)  = %.4e J  (T: %.2f → %.2f K)\n", ΔE_th[end], T0, T_mean[end])

    # 断裂能 E_frac
    # G_c = param_dim.cohesive.G_c_n                                        # TODO Chunk 2 Task 2.1
    G_c = NaN  # TODO Chunk 2 Task 2.1
    coh_lengths = [elem.length * scale.L for elem in czm_mesh.cohesive_elements]

    E_frac = zeros(Float64, nt)
    if haskey(result, "czm damage [0-1]")
        D_all = result["czm damage [0-1]"]
        if ndims(D_all) == 2
            for k in 1:nt
                for e in 1:min(n_coh, size(D_all, 1))
                    E_frac[k] += G_c * coh_lengths[e] * D_all[e, k]
                end
            end
        end
    end
    @printf("  E_frac(终) = %.4e J\n", E_frac[end])

    # 边界热损失 Q_loss
    h_cool   = param_dim.cell.h
    T_amb    = param_dim.cell.T_amb
    A_surface = param_dim.cell.cooling_surface

    Q_loss_inst = h_cool .* A_surface .* (T_mean .- T_amb)
    Q_loss = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        Q_loss[k] = Q_loss[k-1] + 0.5 * (Q_loss_inst[k] + Q_loss_inst[k-1]) * dt_k
    end
    @printf("  Q_loss(终) = %.4e J\n", Q_loss[end])

    # 热源积分 Q_gen
    Q_gen = zeros(Float64, nt)
    if haskey(result, "total heat source [W]")
        Qs = result["total heat source [W]"]
        for k in 2:nt
            dt_k = t[k] - t[k-1]
            Q_gen[k] = Q_gen[k-1] + 0.5 * (Qs[k] + Qs[k-1]) * dt_k
        end
    end
    @printf("  Q_gen(终)  = %.4e J\n", Q_gen[end])

    println("\n[4] 能量残余计算...")
    R = Q_gen .- Q_loss .- ΔE_th .- E_frac

    ε_R = zeros(Float64, nt)
    for k in 1:nt
        denom = max(abs(Q_gen[k]), 1e-12)
        ε_R[k] = abs(R[k]) / denom * 100
    end

    ε_R_rms = sqrt(mean(R[2:end].^2)) / abs(W_elec[end]) * 100

    @printf("  R(终)      = %.4e J\n", R[end])
    @printf("  ε_R(终)    = %.4f %%\n", ε_R[end])
    @printf("  max ε_R    = %.4f %%\n", maximum(ε_R[2:end]))
    @printf("  ε_R,rms   = %.4f %%\n", ε_R_rms)

    println("\n" * "=" ^ 70)
    println("能量平衡汇总")
    println("=" ^ 70)
    @printf("  %-25s %15.4e J\n", "W_elec (电功)", W_elec[end])
    @printf("  %-25s %15.4e J\n", "Q_gen  (热源积分)", Q_gen[end])
    @printf("  %-25s %15.4e J\n", "ΔE_th  (热能变化)", ΔE_th[end])
    @printf("  %-25s %15.4e J\n", "Q_loss (边界热损失)", Q_loss[end])
    @printf("  %-25s %15.4e J\n", "E_frac (断裂能)", E_frac[end])
    @printf("  %-25s %15.4e J\n", "R      (残余)", R[end])
    @printf("  %-25s %15.4f %%\n", "ε_R    (相对误差)", ε_R[end])
    @printf("  %-25s %15.4f %%\n", "ε_R,rms (归一化RMS)", ε_R_rms)

    # ── 绘图 ──
    println("\n[5] 绘图...")
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    p1 = plot(xlabel="Time [s]", ylabel="Energy [J]",
              title="Energy Components", legend=:topleft)
    plot!(p1, t, Q_gen,     label="Q_gen (heat source)",    lw=2, color=:blue)
    plot!(p1, t, ΔE_th,    label="ΔE_th (thermal)",        lw=2, color=:red)
    plot!(p1, t, Q_loss,   label="Q_loss (boundary)",      lw=2, color=:green)
    plot!(p1, t, E_frac,   label="E_frac (fracture)",      lw=2, color=:purple)
    plot!(p1, t, W_elec,   label="W_elec (electrical)",     lw=2, color=:orange, ls=:dash)
    savefig(p1, joinpath(out_dir, "energy_components.png"))

    p2 = plot(xlabel="Time [s]", ylabel="Residual R [J]",
              title="Energy Balance Residual", legend=:topright)
    plot!(p2, t, R, label="R = Q_gen - Q_loss - ΔE_th - E_frac", lw=2, color=:black)
    hline!(p2, [0.0], label="", color=:gray, ls=:dot)
    savefig(p2, joinpath(out_dir, "energy_residual.png"))

    # ── 归一化 RMS 残余报告 ──
    println("\n[6] 生成归一化 RMS 残余报告...")
    report = """
    能量守恒检查 - 归一化 RMS 残余报告
    =================================

    网格分辨率: nθ = $nθ
    计算时间: $dt_wall s

    能量平衡残余统计:
    ----------------
    最终残余 R(终)    = $(@sprintf("%.4e", R[end])) J
    最终相对误差 ε_R = $(@sprintf("%.4f", ε_R[end])) %
    最大相对误差    = $(@sprintf("%.4f", maximum(ε_R[2:end]))) %
    归一化 RMS 残余  = $(@sprintf("%.4f", ε_R_rms)) %

    能量分量汇总:
    ------------
    电功 W_elec      = $(@sprintf("%.4e", W_elec[end])) J
    热源积分 Q_gen   = $(@sprintf("%.4e", Q_gen[end])) J
    热能变化 ΔE_th   = $(@sprintf("%.4e", ΔE_th[end])) J
    边界热损失 Q_loss = $(@sprintf("%.4e", Q_loss[end])) J
    断裂能 E_frac    = $(@sprintf("%.4e", E_frac[end])) J

    能量平衡方程:
    ------------
    Q_gen = ΔE_th + Q_loss + E_frac + R

    检查结果: $(if ε_R_rms < 1.0 "通过" else "警告" end)
    """

    open(joinpath(out_dir, "energy_conservation_report.txt"), "w") do f
        write(f, report)
    end

    println(report)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end