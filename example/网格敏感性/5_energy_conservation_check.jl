"""
Script 5: 全场耦合能量守恒检查

在全耦合算例（SPMe + distributed2D thermal + CZM）上验证：
  W_elec = ΔE_thermal + ΔE_elastic + E_fracture + Q_loss + ΔE_chem

简化方案（首选）：
  Q_generated = ΔE_thermal + Q_loss + ΔE_elastic + E_fracture

其中 Q_generated 由 `total heat source [W]` 积分得到（已可用）。

输出：R(t)、ε_R(t) 图
"""

using Printf, Plots, LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "0_rmspe_utils.jl"))

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

    # ── 配置 ──
    nθ = 80  # 推荐网格（可在收敛分析后调整）

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
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_th, case.param)
    case.czm_mesh = czm_mesh
    n_coh = czm_mesh.n_cohesive

    @printf("  热单元: %d,  节点: %d,  CZM单元: %d\n", ne, mesh_th.nlen, n_coh)

    # ── 求解 ──
    println("\n[2] 运行全耦合求解...")
    t0 = time_ns()
    result = JuBat.Solve(case)
    dt_wall = (time_ns() - t0) * 1e-9
    @printf("  耗时 %.2f s\n", dt_wall)

    # ── 提取时间序列 ──
    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I = result["cell current [A]"]
    nt = length(t)

    # ================================================================
    # [3] 计算各项能量
    # ================================================================
    println("\n[3] 计算能量平衡...")

    # ── 3.1 电功 W_elec = ∫ V·I dt ──
    W_elec = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        W_elec[k] = W_elec[k-1] + 0.5 * (V[k]*I[k] + V[k-1]*I[k-1]) * dt_k
    end
    @printf("  W_elec(终) = %.4e J\n", W_elec[end])

    # ── 3.2 热能变化 ΔE_th ──
    # 从节点温度场积分: ΔE_th = Σ_n M_n * (T_n - T_0)
    # 简化：用体积平均温度 × 总热容
    scale = param_dim.scale
    ρ  = param_dim.cell.rho      # kg/m³ (等效)
    cp = param_dim.cell.heat_Q   # J/kg/K (等效)
    Vol = param_dim.cell.volume  # m³

    T_mean = result["temperature [K]"]
    T0 = T_mean[1]
    C_total = ρ * cp * Vol  # J/K

    ΔE_th = C_total .* (T_mean .- T0)
    @printf("  ΔE_th(终)  = %.4e J  (T: %.2f → %.2f K)\n", ΔE_th[end], T0, T_mean[end])

    # ── 3.3 断裂能 E_frac = Σ_e G_c · l_e · D_e ──
    # G_c = param_dim.cohesive.G_c_n  # J/m²                                # TODO Chunk 2 Task 2.1
    G_c = NaN  # TODO Chunk 2 Task 2.1

    # 获取每个 CZM 单元的长度（使用 create_czm_mesh 预计算值）
    coh_lengths = Float64[]
    for elem in czm_mesh.cohesive_elements
        push!(coh_lengths, elem.length * scale.L)  # 归一化长度 → 物理长度 [m]
    end

    E_frac = zeros(Float64, nt)
    if haskey(result, "czm damage [0-1]")
        D_all = result["czm damage [0-1]"]  # n_coh × nt 或类似
        if ndims(D_all) == 2
            for k in 1:nt
                for e in 1:min(n_coh, size(D_all, 1))
                    E_frac[k] += G_c * coh_lengths[e] * D_all[e, k]
                end
            end
        end
    end
    @printf("  E_frac(终) = %.4e J\n", E_frac[end])

    # ── 3.4 边界热损失 Q_loss = ∫ h·A_s·(T_s - T_amb) dt ──
    # 简化：用体积平均温度近似表面温度
    h_cool   = param_dim.cell.h      # W/m²/K
    T_amb    = param_dim.cell.T_amb  # K
    A_surface = param_dim.cell.cooling_surface  # m²

    Q_loss_inst = h_cool .* A_surface .* (T_mean .- T_amb)  # W
    Q_loss = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        Q_loss[k] = Q_loss[k-1] + 0.5 * (Q_loss_inst[k] + Q_loss_inst[k-1]) * dt_k
    end
    @printf("  Q_loss(终) = %.4e J\n", Q_loss[end])

    # ── 3.5 热源积分 Q_gen = ∫ total_heat_source dt ──
    Q_gen = zeros(Float64, nt)
    if haskey(result, "total heat source [W]")
        Qs = result["total heat source [W]"]
        for k in 2:nt
            dt_k = t[k] - t[k-1]
            Q_gen[k] = Q_gen[k-1] + 0.5 * (Qs[k] + Qs[k-1]) * dt_k
        end
    end
    @printf("  Q_gen(终)  = %.4e J\n", Q_gen[end])

    # ================================================================
    # [4] 能量残余
    # ================================================================
    println("\n[4] 能量残余计算...")

    # 完整方案：R = W_elec - Q_loss - ΔE_th - E_frac - ΔE_chem
    #   其中 ΔE_chem ≈ W_elec - Q_gen（由热源分解隐含）
    #   代入得：R = Q_gen - Q_loss - ΔE_th - E_frac（简化方案）
    R = Q_gen .- Q_loss .- ΔE_th .- E_frac

    # 相对误差
    ε_R = zeros(Float64, nt)
    for k in 1:nt
        denom = max(abs(Q_gen[k]), 1e-12)
        ε_R[k] = abs(R[k]) / denom * 100
    end

    @printf("  R(终)      = %.4e J\n", R[end])
    @printf("  ε_R(终)    = %.4f %%\n", ε_R[end])
    @printf("  max ε_R    = %.4f %%\n", maximum(ε_R[2:end]))

    # ── 归一化 RMS 残余 (spec §3.4) ──
    ε_R_rms = sqrt(mean(R[2:end].^2)) / abs(W_elec[end]) * 100
    @printf("  ε_R,rms   = %.4f %%\n", ε_R_rms)

    # 汇总表
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

    # ================================================================
    # [5] 绘图
    # ================================================================
    println("\n[5] 绘图...")
    out_dir = joinpath(root_dir, "output", "mesh_sensitivity")
    mkpath(out_dir)

    # 图1: 各能量分量随时间变化
    p1 = plot(xlabel="Time [s]", ylabel="Energy [J]",
              title="Energy Components", legend=:topleft)
    plot!(p1, t, Q_gen,     label="Q_gen (heat source)",    lw=2, color=:blue)
    plot!(p1, t, ΔE_th,    label="ΔE_th (thermal)",        lw=2, color=:red)
    plot!(p1, t, Q_loss,   label="Q_loss (boundary)",      lw=2, color=:green)
    plot!(p1, t, E_frac,   label="E_frac (fracture)",      lw=2, color=:purple)
    plot!(p1, t, W_elec,   label="W_elec (electrical)",     lw=2, color=:orange, ls=:dash)
    savefig(p1, joinpath(out_dir, "energy_components.png"))

    # 图2: 残余 R(t)
    p2 = plot(xlabel="Time [s]", ylabel="Residual R [J]",
              title="Energy Balance Residual", legend=:topright)
    plot!(p2, t, R, label="R = Q_gen - Q_loss - ΔE_th - E_frac",
          lw=2, color=:black)
    hline!(p2, [0.0], label="", color=:gray, ls=:dot)
    savefig(p2, joinpath(out_dir, "energy_residual.png"))

    # 图3: 相对误差 ε_R(t)
    p3 = plot(xlabel="Time [s]", ylabel="Relative Error ε_R [%]",
              title="Energy Conservation Error", legend=:topright,
              yscale=:log10, ylims=(1e-6, max(100.0, maximum(ε_R[2:end]) * 2)))
    plot!(p3, t[2:end], ε_R[2:end], label="ε_R", lw=2, color=:blue)
    hline!(p3, [1.0],  label="1%",  color=:green, ls=:dash)
    hline!(p3, [5.0],  label="5%",  color=:red,   ls=:dash)
    savefig(p3, joinpath(out_dir, "energy_relative_error.png"))

    # 图4: 分解验证 Q_gen = ΔE_th + Q_loss + E_frac + R
    p4 = plot(xlabel="Time [s]", ylabel="Energy [J]",
              title="Energy Balance: Q_gen Decomposition", legend=:topleft)
    plot!(p4, t, Q_gen, label="Q_gen", lw=3, color=:black)
    plot!(p4, t, ΔE_th .+ Q_loss .+ E_frac, label="ΔE_th + Q_loss + E_frac",
          lw=2, color=:blue, ls=:dash)
    plot!(p4, t, ΔE_th .+ Q_loss .+ E_frac .+ R, label="+ R (should ≈ Q_gen)",
          lw=2, color=:red, ls=:dot)
    savefig(p4, joinpath(out_dir, "energy_decomposition.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
    println("能量守恒检查完成")
    println("=" ^ 70)
end

main()
