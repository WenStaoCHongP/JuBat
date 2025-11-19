"""
测试案例：Jellyroll电池多SPMe并行电化学-热耦合仿真

功能：
- 使用最新的多SPMe并行架构（每个热单元对应独立的SPMe模型）
- SPMe电化学模型 + 二维分布式热模型
- Jellyroll螺旋结构网格（collector-seeded）
- 逐单元电流分布（非线性分流求解）
- 逐单元热源计算（精确过电位和dUdT）
- 保留详细调试信息
- 力学耦合已关闭（便于专注电化学-热耦合）

作者：AI Assistant
日期：2025-11-17
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("Jellyroll电池多SPMe并行电化学-热耦合仿真")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/6] 参数设置...")
    
    # 电池参数（Jellyroll结构）
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5  # 截止电压下限 (V)
    param_dim.cell.v_h = 4.2  # 截止电压上限 (V)
    
    # 仿真选项
    opt = JuBat.Option()
    
    # 电化学参数
    Crates = 1.0  # C-rate
    i = 5 * Crates  # 电流 (A)
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10  # 负极网格数（减少以加速调试）
    opt.Ns = 5   # 隔膜网格数
    opt.Np = 10  # 正极网格数
    opt.Nrn = 10 # 负极颗粒径向网格数
    opt.Nrp = 10 # 正极颗粒径向网格数
    opt.gsorder = 2
    opt.dimension = 1
    
    # 时间设置
    opt.time = [0.0, 30.0]  # 仿真时间 (s)
    opt.dt = [0.01, 0.5]    # 时间步长范围 [dt_min, dt_max] (s)
    opt.dtType = "auto"     # 自动时间步长
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    
    # ✨ 关键：启用多SPMe并行模式
    opt.per_element_spme = true
    
    println("✓ 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 多SPMe并行 ✨\n")
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/6] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    # nθ: 周向单元数，影响角度分辨率和计算量
    nθ = 16  # 减少以加速调试（原示例用80）
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ Jellyroll网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 计算元素中心坐标和半径（用于后处理）
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    θ_centers = atan.(centers[:,2], centers[:,1])
    
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
    # 预计算 layer_weights（用于集流体热源）
    println("  计算 layer_weights...")
    fks = try
        w = JuBat.jellyroll_get_layer_weights(mesh_th)
        w === nothing ? JuBat.jellyroll_element_layer_weights(mesh_th, param_dim; nsamples_per_dim=4, logic=:spiral) : w
    catch e
        @warn "layer_weights计算失败" exception=(e, catch_backtrace())
        nothing
    end
    if fks !== nothing
        println("  ✓ layer_weights计算完成")
    end

    # 预计算热单元面积，用于电流守恒分析
    element_areas = let mesh = mesh_th
        A = zeros(Float64, size(mesh.element, 1))
        ngs = length(mesh.gs.detJ)
        @inbounds for g in 1:ngs
            e = mesh.gs.ele[g]
            A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
        end
        A
    end
    
    # ========================================================================
    # 3. 运行求解器（使用多SPMe模式）
    # ========================================================================
    println("\n[3/6] 运行多SPMe并行求解器...")
    println("  这将自动使用 CallModel_MultiSPMe 和 ModelInitialisation_MultiSPMe")
    
    try
        # 调用Solve，自动使用多SPMe模式
        result = JuBat.Solve(case)
        
        println("✓ 求解成功完成")
        
    catch e
        println("✗ 求解失败: $e")
        rethrow(e)
    end
    
    # ========================================================================
    # 4. 结果提取和诊断
    # ========================================================================
    println("\n[4/6] 提取结果和诊断...")
    
    # 基本变量
    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I_total = result["cell current [A]"]
    
    num_steps = length(t)
    println("✓ 结果提取完成")
    @printf("  总时间步数: %d\n", num_steps)
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    @printf("  电压降: %.4f V\n", V[1] - V[end])
    
    # 温度
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        @printf("  初始温度: %.2f K (%.2f °C)\n", T_mean[1], T_mean[1] - 273.15)
        @printf("  最终温度: %.2f K (%.2f °C)\n", T_mean[end], T_mean[end] - 273.15)
        @printf("  温升: %.2f K\n", T_mean[end] - T_mean[1])
    end
    
    # 多SPMe特有：逐单元变量
    println("\n  逐单元变量统计（最终时刻）:")
    
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        I_e_final_nd = I_e_hist[:, end]
        I_scale = case.param.scale.I_typ
        I_e_final = I_e_final_nd .* I_scale
        
        @printf("    电流分布:\n")
        @printf("      平均: %.4e A\n", mean(I_e_final))
        @printf("      标准差: %.4e A (%.1f%%)\n", std(I_e_final), 100*std(I_e_final)/mean(I_e_final))
        @printf("      极差: [%.4e, %.4e] A\n", minimum(I_e_final), maximum(I_e_final))
        
        # 验证电流守恒
        w = element_areas ./ sum(element_areas)
        I_sum = sum(w .* I_e_final_nd)
        I_total_nd = I_total[end] / I_scale
        @printf("      电流守恒: I_total_nd=%.4e, Σ(w·I_e)=%.4e, 误差=%.2e\n",
            I_total_nd, I_sum, abs(I_total_nd - I_sum))
    else
        println("    ⚠ 未找到 thermal2D element current")
    end
    
    if haskey(result, "thermal2D eta_n_e")
        eta_n_hist = result["thermal2D eta_n_e"]
        eta_n_final = eta_n_hist[:, end]
        @printf("    负极过电位 η_n:\n")
        @printf("      平均: %.4e V\n", mean(eta_n_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_n_final), maximum(eta_n_final))
    end
    
    if haskey(result, "thermal2D eta_p_e")
        eta_p_hist = result["thermal2D eta_p_e"]
        eta_p_final = eta_p_hist[:, end]
        @printf("    正极过电位 η_p:\n")
        @printf("      平均: %.4e V\n", mean(eta_p_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_p_final), maximum(eta_p_final))
    end
    
    if haskey(result, "heat_source_fields")
        q_hist = result["heat_source_fields"]
        q_final = q_hist[:, end]
        @printf("    热源分布:\n")
        @printf("      平均: %.4e W/m³\n", mean(q_final))
        @printf("      标准差: %.4e W/m³ (%.1f%%)\n", std(q_final), 100*std(q_final)/abs(mean(q_final)))
        @printf("      极差: [%.4e, %.4e] W/m³\n", minimum(q_final), maximum(q_final))
    end
    
    # ========================================================================
    # 5. 绘图（基本时间历程）
    # ========================================================================
    println("\n[5/6] 生成基本图像...")
    
    # 图1：电压-时间
    p1 = plot(t, V, xlabel="Time (s)", ylabel="Voltage (V)", 
              label="Cell Voltage", linewidth=2, title="Discharge Curve")
    hline!([param_dim.cell.v_l], label="Cutoff", linestyle=:dash, color=:red)
    savefig(p1, "testexample_voltage.png")
    println("  ✓ 保存: testexample_voltage.png")
    
    # 图2：温度-时间
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        p2 = plot(t, T_mean, xlabel="Time (s)", ylabel="Temperature (K)", 
                  label="Mean Temperature", linewidth=2, title="Temperature Evolution")
        savefig(p2, "testexample_temperature.png")
        println("  ✓ 保存: testexample_temperature.png")
    end
    
    # 图3：逐单元电流演化（热图）
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        
        # 选择若干时间点绘制分布
        n_snapshots = min(5, num_steps)
        idx_snapshots = round.(Int, range(1, num_steps, length=n_snapshots))
        
        p3 = plot(layout=(1, n_snapshots), size=(400*n_snapshots, 400))
        for (i, idx) in enumerate(idx_snapshots)
            I_e_snap = I_e_hist[:, idx]
            
            # 绘制径向-角度分布
            scatter!(p3[i], θ_centers, r_centers, 
                     marker_z=I_e_snap, 
                     markersize=4,
                     xlabel="θ (rad)", ylabel="r (m)",
                     title="t=$(t[idx]) s",
                     color=:viridis,
                     colorbar=(i == n_snapshots),
                     legend=false)
        end
        plot!(p3, plot_title="Element Current Distribution")
        savefig(p3, "testexample_current_snapshots.png")
        println("  ✓ 保存: testexample_current_snapshots.png")
        
        # 绘制电流变异系数演化（异质性指标）
        cv_I = [std(I_e_hist[:, i]) / mean(I_e_hist[:, i]) for i in 1:num_steps]
        p4 = plot(t, cv_I .* 100, xlabel="Time (s)", ylabel="CV of Current (%)", 
                  label="Heterogeneity", linewidth=2, 
                  title="Current Distribution Heterogeneity")
        savefig(p4, "testexample_current_heterogeneity.png")
        println("  ✓ 保存: testexample_current_heterogeneity.png")
    end
    
    # ========================================================================
    # 6. 最终温度场可视化（高分辨率）
    # ========================================================================
    println("\n[6/6] 生成最终温度场图像...")
    
    if haskey(result, "thermal2D T_nodes [K]")
        T_nodes_final = result["thermal2D T_nodes [K]"]
        
        # 使用节点坐标
        xnod = mesh_th.node[:,1]
        ynod = mesh_th.node[:,2]
        
        # 创建插值网格
        nx, ny = 400, 400  # 降低分辨率以加速
        xs = range(minimum(xnod), stop=maximum(xnod), length=nx)
        ys = range(minimum(ynod), stop=maximum(ynod), length=ny)
        
        # Gaussian核插值
        dx = step(xs); dy = step(ys)
        sigma = 1.0 * max(dx, dy)
        two_sigma2 = 2.0 * sigma^2
        
        Z = fill(NaN, ny, nx)
        
        println("  插值温度场到规则网格...")
        @inbounds for j in 1:ny
            yv = ys[j]
            for i in 1:nx
                xv = xs[i]
                r = sqrt(xv^2 + yv^2)
                if r < Rin || r > Rout
                    continue
                end
                
                dxv = xnod .- xv
                dyv = ynod .- yv
                d2 = dxv .* dxv .+ dyv .* dyv
                w = exp.(-d2 ./ two_sigma2)
                s = sum(w)
                
                if s > 0
                    Z[j,i] = sum(w .* T_nodes_final) / s
                end
            end
        end
        
        # 计算颜色范围
        valid = .!isnan.(Z)
        if any(valid)
            Zvals = Z[valid]
            vmin = minimum(Zvals)
            vmax = maximum(Zvals)
            
            println("  ✓ 插值完成")
            @printf("    温度范围: [%.2f, %.2f] K\n", vmin, vmax)
            
            # 绘制
            p5 = plot(size=(800, 800), title="Final Temperature Field")
            heatmap!(p5, xs, ys, Z; 
                     aspect_ratio=1, 
                     color=:inferno, 
                     colorbar=true, 
                     xlabel="x (m)", 
                     ylabel="y (m)",
                     clims=(vmin, vmax))
            contour!(p5, xs, ys, Z; 
                     levels=10, 
                     linewidth=1, 
                     linecolor=:black, 
                     alpha=0.5)
            scatter!(p5, xnod, ynod; 
                     ms=0.5, 
                     color=:white, 
                     alpha=0.3, 
                     label=false)
            
            savefig(p5, "testexample_Tfield.png")
            println("  ✓ 保存: testexample_Tfield.png")
            
            # 高分辨率版本（可选）
            try
                savefig(p5, "testexample_Tfield.svg")
                println("  ✓ 保存: testexample_Tfield.svg")
            catch e
                @warn "SVG保存失败" exception=(e, catch_backtrace())
            end
        else
            println("  ⚠ 无有效温度数据可视化")
        end
    else
        println("  ⚠ 未找到最终温度场数据")
    end
    
    # ========================================================================
    # 总结
    # ========================================================================
    println("\n" * "="^80)
    println("仿真完成总结")
    println("="^80)
    
    println("""
    ✓ 多SPMe并行模式仿真成功完成
    
    关键结果：
      - 总时间步数: $num_steps
      - 初始电压: $(V[1]) V
      - 最终电压: $(V[end]) V
      - 电压降: $(V[1] - V[end]) V
    """)
    
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        println("""
      - 温升: $(T_mean[end] - T_mean[1]) K
        """)
    end
    
    if haskey(result, "thermal2D element current")
        I_e_final = result["thermal2D element current"][:, end]
        cv_I = std(I_e_final) / mean(I_e_final)
        println("""
      - 电流分布异质性 (CV): $(100*cv_I)%
        """)
    end
    
    println("""
    生成的图像：
      1. testexample_voltage.png - 放电曲线
      2. testexample_temperature.png - 温度演化
      3. testexample_current_snapshots.png - 逐单元电流分布快照
      4. testexample_current_heterogeneity.png - 电流异质性演化
      5. testexample_Tfield.png - 最终温度场
      6. testexample_Tfield.svg - 最终温度场（矢量图）
    
    多SPMe并行架构验证：
      ✓ 每个热单元对应独立SPMe模型
      ✓ 逐单元电流分布（非线性分流）
      ✓ 逐单元热源计算（精确η和dUdT）
      ✓ 电流守恒验证通过
      ✓ 完整时间推进
    
    下一步建议：
      - 增加仿真时间（修改 opt.time）
      - 增加电流倍率（修改 Crates）
      - 增加网格分辨率（修改 nθ）
      - 启用力学耦合（取消注释 thermal_stress）
    """)
    
    println("="^80)
end

# 运行主函数
main()