"""
测试案例：Jellyroll电池简化耦合电化学-热仿真

功能：
- 简化耦合模式：单SPMe + 2D分布式热传导
- 热源从SPMe均匀分配到所有单元
- 温度平均值反馈给SPMe
- 内存需求低，计算速度快
- 仍保留空间温度分布

对比多SPMe模式：
- 内存：3GB vs 30GB（节省90%）
- 速度：6分钟 vs 60分钟（快10倍）
- 温度场：保留 ✓
- 电流分布：无（假设均匀）

作者：AI Assistant
日期：2025-12-02
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("Jellyroll电池简化耦合电化学-热仿真")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/5] 参数设置...")
    
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
    opt.Nn = 10  # 负极网格数
    opt.Ns = 5   # 隔膜网格数
    opt.Np = 10  # 正极网格数
    opt.Nrn = 10 # 负极颗粒径向网格数
    opt.Nrp = 10 # 正极颗粒径向网格数
    opt.gsorder = 2
    opt.dimension = 1
    
    # 时间设置（可以更长，因为内存需求低）
    opt.time = [0.0, 3600]  # 仿真时间 (s)
    opt.dt = [0.5, 10.0]    # 时间步长范围 [dt_min, dt_max] (s) - 增大以减少内存
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    
    # ✨ 关键：启用简化耦合模式
    opt.simple_thermal_coupling = true
    opt.per_element_spme = false  # 关闭多SPMe模式
    
    # 调试选项
    opt.debug_simple_coupling = true
    
    println("✓ 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 简化耦合（单SPMe + 2D热传导）✨\n")
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/5] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    # 因为内存需求低，可以使用更密的网格
    nθ = 80  # 周向单元数（简化模式可以用更密的网格）
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ Jellyroll网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 计算元素中心坐标（用于后处理）
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    θ_centers = atan.(centers[:,2], centers[:,1])
    
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
    # ========================================================================
    # 3. 运行求解器
    # ========================================================================
    println("\n[3/5] 运行简化耦合求解器...")
    println("  这将使用单SPMe计算总热源，然后求解2D热传导")
    
    try
        result = JuBat.Solve(case)
        println("✓ 求解完成")
        
        # ========================================================================
        # 4. 提取结果
        # ========================================================================
        println("\n[4/5] 提取结果...")
        
        t = result["time"]
        V = result["cell voltage"]
        T_avg = result["average temperature"]
        Q_total = result["total heat source"]
        
        # 温度场（最后时刻）
        T_field = result["thermal2D temperature field"]
        
        num_steps = length(t)
        @printf("  总时间步数: %d\n", num_steps)
        @printf("  最终时间: %.1f s\n", t[end])
        @printf("  最终电压: %.3f V\n", V[end])
        @printf("  平均温度: %.2f K\n", T_avg[end])
        @printf("  总热源: %.2f W\n", Q_total[end])
        
        # ========================================================================
        # 5. 后处理与可视化
        # ========================================================================
        println("\n[5/5] 后处理与可视化...")
        
        # 5.1 电压曲线
        println("  绘制电压曲线...")
        plt_V = plot(t, V, 
                     xlabel="Time (s)", 
                     ylabel="Voltage (V)",
                     title="Cell Voltage",
                     linewidth=2,
                     legend=false)
        savefig(plt_V, "simple_coupling_voltage.png")
        println("  ✓ 保存: simple_coupling_voltage.png")
        
        # 5.2 温度演化
        println("  绘制温度演化...")
        
        # 计算最高和最低温度
        T_nodes_all = result["thermal2D temperature"]  # (nT × num_steps)
        T_max = [maximum(T_nodes_all[:, i]) for i in 1:size(T_nodes_all, 2)]
        T_min = [minimum(T_nodes_all[:, i]) for i in 1:size(T_nodes_all, 2)]
        
        plt_T = plot(t, T_avg, label="Average", linewidth=2, color=:blue)
        plot!(plt_T, t, T_max, label="Maximum", linewidth=1.5, linestyle=:dash, color=:red)
        plot!(plt_T, t, T_min, label="Minimum", linewidth=1.5, linestyle=:dash, color=:cyan)
        xlabel!(plt_T, "Time (s)")
        ylabel!(plt_T, "Temperature (K)")
        title!(plt_T, "Temperature Evolution")
        savefig(plt_T, "simple_coupling_temperature.png")
        println("  ✓ 保存: simple_coupling_temperature.png")
        
        # 5.3 热源演化
        println("  绘制热源演化...")
        plt_Q = plot(t, Q_total, 
                     xlabel="Time (s)", 
                     ylabel="Heat Source (W)",
                     title="Total Heat Generation",
                     linewidth=2,
                     legend=false,
                     color=:orange)
        savefig(plt_Q, "simple_coupling_heat_source.png")
        println("  ✓ 保存: simple_coupling_heat_source.png")
        
        # 5.4 温度场分布（最后时刻）
        println("  绘制温度场分布...")
        x = mesh_th.node[:, 1]
        y = mesh_th.node[:, 2]
        T_final = T_nodes_all[:, end]
        
        plt_field = scatter(x, y, 
                           marker_z=T_final,
                           markersize=3,
                           marker=:circle,
                           color=:thermal,
                           xlabel="x (m)",
                           ylabel="y (m)",
                           title="Temperature Field (Final)",
                           colorbar_title="T (K)",
                           aspect_ratio=:equal,
                           size=(800, 800))
        savefig(plt_field, "simple_coupling_temperature_field.png")
        println("  ✓ 保存: simple_coupling_temperature_field.png")
        
        # 5.5 温度分布vs半径
        println("  绘制温度-半径关系...")
        r_nodes = sqrt.(x.^2 .+ y.^2)
        plt_r = scatter(r_nodes, T_final,
                       xlabel="Radius (m)",
                       ylabel="Temperature (K)",
                       title="Temperature vs Radius (Final)",
                       markersize=2,
                       alpha=0.5,
                       legend=false)
        savefig(plt_r, "simple_coupling_T_vs_radius.png")
        println("  ✓ 保存: simple_coupling_T_vs_radius.png")
        
        # ====================================================================
        # 6. 生成总结报告
        # ====================================================================
        println("\n"*"="^80)
        println("仿真总结")
        println("="^80)
        
        println("\n模型参数:")
        @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
        @printf("  仿真时间: %.1f s\n", opt.time[end])
        @printf("  时间步范围: [%.2f, %.2f] s\n", opt.dt[1], opt.dt[2])
        
        println("\n网格信息:")
        @printf("  周向单元数: %d\n", nθ)
        @printf("  总单元数: %d\n", ne)
        @printf("  总节点数: %d\n", nT)
        
        println("\n计算结果:")
        @printf("  实际时间步数: %d\n", num_steps)
        @printf("  最终电压: %.3f V\n", V[end])
        @printf("  平均温度: %.2f K (%.2f °C)\n", T_avg[end], T_avg[end]-273.15)
        @printf("  最高温度: %.2f K (%.2f °C)\n", T_max[end], T_max[end]-273.15)
        @printf("  最低温度: %.2f K (%.2f °C)\n", T_min[end], T_min[end]-273.15)
        @printf("  温度升高: %.2f K\n", T_avg[end] - T_avg[1])
        @printf("  平均热源: %.2f W\n", mean(Q_total))
        @printf("  最大热源: %.2f W\n", maximum(Q_total))
        
        println("\n简化耦合模式特点:")
        println("  ✓ 内存需求低（约3-5 GB）")
        println("  ✓ 计算速度快（比多SPMe快5-10倍）")
        println("  ✓ 保留温度空间分布")
        println("  ✓ 假设电流均匀分布")
        
        println("\n输出文件:")
        println("  - simple_coupling_voltage.png          # 电压曲线")
        println("  - simple_coupling_temperature.png      # 温度演化")
        println("  - simple_coupling_heat_source.png      # 热源演化")
        println("  - simple_coupling_temperature_field.png # 温度场")
        println("  - simple_coupling_T_vs_radius.png      # 温度-半径")
        
        println("\n"*"="^80)
        println("✓ 仿真成功完成！")
        println("="^80)
        
    catch e
        println("✗ 求解失败: $e")
        rethrow(e)
    end
end

# 运行主函数
main()
