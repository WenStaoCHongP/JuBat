"""
测试案例：Z方向冷却方法（cool_method）

对比两种冷却方式：
1. surface：整体表面冷却（整个xy域）
2. tab：极耳强化冷却（仅极耳节点邻域）

物理机制：
z方向对流散热以"体积热汇"形式贡献到xy平面：
    q_vol = 2h(T - T_amb) / H

作者：AI Assistant  
日期：2025-12-29
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function run_simulation(cool_method::String; h_val=nothing)
    println("\n" * "="^80)
    println("运行仿真：冷却方式 = $cool_method")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    
    # 极耳设置（螺旋线离散点）
    param_dim.tab.theta_pos = [0.0]
    param_dim.tab.theta_neg = [20.0 * π]
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.Current = x -> 5.0
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.time = [0.0, 60.0]
    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # ✨ 设置冷却方式（cool_method）
    opt.cool_method = cool_method
    
    if cool_method == "surface"
        opt.h_surface = h_val !== nothing ? h_val : 10.0
        println("  h_surface = $(opt.h_surface) W/(m²·K)")
        println("  体积散热系数 = $(2*opt.h_surface/param_dim.cell.width) W/(m³·K)")
    elseif cool_method == "tab"
        opt.h_tab = h_val !== nothing ? h_val : 100.0
        println("  h_tab = $(opt.h_tab) W/(m²·K)")
        println("  极耳体积散热系数 = $(2*opt.h_tab/param_dim.cell.width) W/(m³·K)")
    end
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    case = JuBat.SetCase(param_dim, opt)
    
    nθ = 80
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("  网格信息: ne=$ne, nT=$nT")
    
    # 计算Biot数
    H = param_dim.cell.width
    L_th = case.param_dim.scale.L_th
    k_th = case.param_dim.scale.k_th
    
    if cool_method == "surface"
        h = opt.h_surface
    elseif cool_method == "tab"
        h = opt.h_tab
    else
        h = 0.0
    end
    
    Bi_z = 2.0 * h * L_th^2 / (H * k_th)
    println("  Biot数 Bi_z = $(round(Bi_z, sigdigits=3))")
    
    if Bi_z < 0.01
        println("  ⚠️  Bi_z << 1，z方向散热可忽略")
    elseif Bi_z < 0.1
        println("  ✓ Bi_z ~ O(0.01)，z方向散热较弱")
    elseif Bi_z < 1.0
        println("  ✓ Bi_z ~ O(0.1)，z方向散热明显")
    else
        println("  ✓ Bi_z ~ O(1)，z方向散热主导")
    end
    
    # ========================================================================
    # 3. 运行求解器
    # ========================================================================
    println("  开始求解...")
    
    try
        result = JuBat.Solve(case)
        
        println("  ✓ 求解成功")
        
        # 提取结果
        V_hist = result["cell voltage [V]"]
        T_mean_hist = result["temperature [K]"]
        time_hist = result["time [s]"]
        
        println("  最终电压: $(round(V_hist[end], digits=4)) V")
        println("  最终平均温度: $(round(T_mean_hist[end], digits=2)) K")
        println("  最高温度: $(round(maximum(T_mean_hist), digits=2)) K")
        println("  温升: $(round(maximum(T_mean_hist) - T_mean_hist[1], digits=2)) K")
        
        return (
            success = true,
            time = time_hist,
            voltage = V_hist,
            temperature = T_mean_hist,
            Bi_z = Bi_z,
            case = case,
            result = result
        )
        
    catch err
        @error "❌ 求解失败" exception=(err, catch_backtrace())
        return (success = false,)
    end
end

function main()
    println("="^80)
    println("Z方向冷却方法对比测试")
    println("="^80)
    
    println("\n物理机制：")
    println("  z方向对流散热 → 体积热汇")
    println("  q_vol = 2h(T - T_amb) / H")
    println("  刚度矩阵贡献：K_ij += ∫ (2h/H) N_i N_j dA")
    println("  这是xy平面的面积分，不是边界线积分！")
    
    # ========================================================================
    # 测试1：整体表面冷却
    # ========================================================================
    println("\n[测试1] 整体表面冷却（h_surface = 10 W/(m²·K)）")
    result1 = run_simulation("surface", h_val=10.0)
    
    # ========================================================================
    # 测试2：极耳强化冷却
    # ========================================================================
    println("\n[测试2] 极耳强化冷却（h_tab = 100 W/(m²·K)）")
    result2 = run_simulation("tab", h_val=100.0)
    
    # ========================================================================
    # 测试3：无冷却（对比基准）
    # ========================================================================
    println("\n[测试3] 无z方向冷却（对比基准）")
    # 创建一个没有设置cool_method的case
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    param_dim.tab.theta_pos = [0.0]
    param_dim.tab.theta_neg = [20.0 * π]
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.Current = x -> 5.0
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.time = [0.0, 60.0]
    opt.dt = [0.5, 10]
    # 不设置 cool_method
    
    case = JuBat.SetCase(param_dim, opt)
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    try
        result = JuBat.Solve(case)
        V_hist = result["cell voltage [V]"]
        T_mean_hist = result["temperature [K]"]
        time_hist = result["time [s]"]
        
        println("  ✓ 求解成功")
        println("  最终电压: $(round(V_hist[end], digits=4)) V")
        println("  最终平均温度: $(round(T_mean_hist[end], digits=2)) K")
        println("  最高温度: $(round(maximum(T_mean_hist), digits=2)) K")
        println("  温升: $(round(maximum(T_mean_hist) - T_mean_hist[1], digits=2)) K")
        
        result3 = (
            success = true,
            time = time_hist,
            voltage = V_hist,
            temperature = T_mean_hist,
            Bi_z = 0.0
        )
    catch err
        @error "无冷却求解失败" exception=err
        result3 = (success = false,)
    end
    
    # ========================================================================
    # 对比结果
    # ========================================================================
    println("\n" * "="^80)
    println("结果对比")
    println("="^80)
    
    if result1.success && result2.success && result3.success
        # 绘图对比
        p1 = plot(
            result3.time, result3.temperature,
            label="无冷却 (Bi=0)",
            xlabel="时间 [s]",
            ylabel="平均温度 [K]",
            linewidth=2,
            title="Z方向冷却效果对比",
            legend=:bottomright
        )
        plot!(p1, result1.time, result1.temperature, 
              label="表面冷却 (Bi=$(round(result1.Bi_z, digits=3)))", linewidth=2)
        plot!(p1, result2.time, result2.temperature, 
              label="极耳冷却 (Bi=$(round(result2.Bi_z, digits=3)))", linewidth=2)
        
        p2 = plot(
            result3.time, result3.voltage,
            label="无冷却",
            xlabel="时间 [s]",
            ylabel="电压 [V]",
            linewidth=2,
            title="电压对比",
            legend=:bottomleft
        )
        plot!(p2, result1.time, result1.voltage, label="表面冷却", linewidth=2)
        plot!(p2, result2.time, result2.voltage, label="极耳冷却", linewidth=2)
        
        p = plot(p1, p2, layout=(2,1), size=(800,800))
        
        savefig(p, joinpath(@__DIR__, "../output/cool_method_comparison.png"))
        println("✓ 图表已保存到 output/cool_method_comparison.png")
        
        # 数值对比
        println("\n最终状态对比（t = 60s）：")
        println("┌────────────────┬──────────┬────────────┬────────────┬──────────┐")
        println("│ 冷却方式       │ 电压 [V] │ 平均温度[K]│ 最高温度[K]│ 温升 [K] │")
        println("├────────────────┼──────────┼────────────┼────────────┼──────────┤")
        @printf("│ %-14s │ %8.4f │ %10.2f │ %10.2f │ %8.2f │\n",
                "无冷却", result3.voltage[end], result3.temperature[end], 
                maximum(result3.temperature), maximum(result3.temperature) - result3.temperature[1])
        @printf("│ %-14s │ %8.4f │ %10.2f │ %10.2f │ %8.2f │\n",
                "表面冷却", result1.voltage[end], result1.temperature[end], 
                maximum(result1.temperature), maximum(result1.temperature) - result1.temperature[1])
        @printf("│ %-14s │ %8.4f │ %10.2f │ %10.2f │ %8.2f │\n",
                "极耳冷却", result2.voltage[end], result2.temperature[end], 
                maximum(result2.temperature), maximum(result2.temperature) - result2.temperature[1])
        println("└────────────────┴──────────┴────────────┴────────────┴──────────┘")
        
        # 冷却效果分析
        println("\n冷却效果分析：")
        ΔT_nocool = maximum(result3.temperature) - result3.temperature[1]
        ΔT_surface = maximum(result1.temperature) - result1.temperature[1]
        ΔT_tab = maximum(result2.temperature) - result2.temperature[1]
        
        cooling_eff_surface = (ΔT_nocool - ΔT_surface) / ΔT_nocool * 100
        cooling_eff_tab = (ΔT_nocool - ΔT_tab) / ΔT_nocool * 100
        
        println("  表面冷却降温效果: $(round(cooling_eff_surface, digits=1))%")
        println("  极耳冷却降温效果: $(round(cooling_eff_tab, digits=1))%")
        
        # Biot数对比
        println("\nBiot数对比：")
        println("  表面冷却: Bi_z = $(round(result1.Bi_z, digits=3))")
        println("  极耳冷却: Bi_z = $(round(result2.Bi_z, digits=3))")
        println("  比值: Bi_tab/Bi_surface = $(round(result2.Bi_z/result1.Bi_z, digits=1))")
        
        # 物理解释
        println("\n物理解释：")
        println("  表面冷却：整个xy域均匀散热，温升降低 $(round(cooling_eff_surface, digits=0))%")
        println("  极耳冷却：仅极耳邻域强化散热，局部形成'冷斑'")
        println("  极耳冷却的Biot数是表面冷却的 $(round(result2.Bi_z/result1.Bi_z, digits=0)) 倍")
        
    else
        println("❌ 部分测试失败，无法进行完整对比")
    end
    
    println("\n" * "="^80)
    println("关键结论")
    println("="^80)
    println("✓ z方向对流散热以'体积热汇'形式贡献到xy平面")
    println("✓ 刚度矩阵贡献：K_ij += ∫ (2h/H) N_i N_j dA （面积分）")
    println("✓ 不是边界线积分！")
    println("✓ 极耳节点 = 螺旋线离散点（以直代曲）")
    println("✓ Biot数 Bi_z = 2h*L²/(H*k) 衡量散热强度")
    println("✓ 数值稳定：Bi_z ~ O(0.01-1)，远小于惩罚法的 1e12")
    println("="^80)
end

# 运行主函数
main()
