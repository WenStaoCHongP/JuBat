"""
测试案例：极耳对流散热边界条件

对比三种边界处理方式：
1. 惩罚法（penalty）
2. 一般表面散热（surface_convection）
3. 极耳强化散热（tab_convection）

作者：AI Assistant
日期：2025-12-29
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function run_simulation(bc_type::String; h_val=nothing, penalty_val=nothing)
    println("\n" * "="^80)
    println("运行仿真：边界类型 = $bc_type")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    
    # 极耳设置
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
    
    # ✨ 设置边界条件类型
    opt.tab_bc_type = bc_type
    
    if bc_type == "surface_convection"
        opt.h_surface = h_val !== nothing ? h_val : 10.0
        println("  h_surface = $(opt.h_surface) W/(m²·K)")
    elseif bc_type == "tab_convection"
        opt.h_tab = h_val !== nothing ? h_val : 100.0
        println("  h_tab = $(opt.h_tab) W/(m²·K)")
    elseif bc_type == "penalty"
        opt.tab_penalty = penalty_val !== nothing ? penalty_val : 1e6
        println("  penalty = $(opt.tab_penalty)")
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
        
        return (
            success = true,
            time = time_hist,
            voltage = V_hist,
            temperature = T_mean_hist,
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
    println("极耳边界条件对比测试")
    println("="^80)
    
    # ========================================================================
    # 测试1：惩罚法（降低惩罚值）
    # ========================================================================
    println("\n[测试1] 惩罚法（penalty = 1e6）")
    result1 = run_simulation("penalty", penalty_val=1e6)
    
    # ========================================================================
    # 测试2：一般表面散热
    # ========================================================================
    println("\n[测试2] 一般表面散热（h_surface = 10 W/(m²·K)）")
    result2 = run_simulation("surface_convection", h_val=10.0)
    
    # ========================================================================
    # 测试3：极耳强化散热
    # ========================================================================
    println("\n[测试3] 极耳强化散热（h_tab = 100 W/(m²·K)）")
    result3 = run_simulation("tab_convection", h_val=100.0)
    
    # ========================================================================
    # 对比结果
    # ========================================================================
    println("\n" * "="^80)
    println("结果对比")
    println("="^80)
    
    if result1.success && result2.success && result3.success
        # 绘图对比
        p1 = plot(
            result1.time, result1.voltage,
            label="惩罚法 (1e6)",
            xlabel="时间 [s]",
            ylabel="电压 [V]",
            linewidth=2,
            title="电压对比"
        )
        plot!(p1, result2.time, result2.voltage, label="表面散热", linewidth=2)
        plot!(p1, result3.time, result3.voltage, label="极耳散热", linewidth=2)
        
        p2 = plot(
            result1.time, result1.temperature,
            label="惩罚法 (1e6)",
            xlabel="时间 [s]",
            ylabel="平均温度 [K]",
            linewidth=2,
            title="温度对比"
        )
        plot!(p2, result2.time, result2.temperature, label="表面散热", linewidth=2)
        plot!(p2, result3.time, result3.temperature, label="极耳散热", linewidth=2)
        
        p = plot(p1, p2, layout=(2,1), size=(800,800))
        
        savefig(p, joinpath(@__DIR__, "../output/tab_bc_comparison.png"))
        println("✓ 图表已保存到 output/tab_bc_comparison.png")
        
        # 数值对比
        println("\n最终状态对比（t = 60s）：")
        println("┌─────────────────┬──────────┬────────────┬────────────┐")
        println("│ 边界类型        │ 电压 [V] │ 平均温度[K]│ 最高温度[K]│")
        println("├─────────────────┼──────────┼────────────┼────────────┤")
        @printf("│ %-15s │ %8.4f │ %10.2f │ %10.2f │\n", 
                "惩罚法", result1.voltage[end], result1.temperature[end], maximum(result1.temperature))
        @printf("│ %-15s │ %8.4f │ %10.2f │ %10.2f │\n",
                "表面散热", result2.voltage[end], result2.temperature[end], maximum(result2.temperature))
        @printf("│ %-15s │ %8.4f │ %10.2f │ %10.2f │\n",
                "极耳散热", result3.voltage[end], result3.temperature[end], maximum(result3.temperature))
        println("└─────────────────┴──────────┴────────────┴────────────┘")
        
        # 分析差异
        println("\n差异分析：")
        ΔT_12 = result1.temperature[end] - result2.temperature[end]
        ΔT_13 = result1.temperature[end] - result3.temperature[end]
        ΔT_23 = result2.temperature[end] - result3.temperature[end]
        
        println("  惩罚法 vs 表面散热：ΔT = $(round(ΔT_12, digits=2)) K")
        println("  惩罚法 vs 极耳散热：ΔT = $(round(ΔT_13, digits=2)) K")
        println("  表面散热 vs 极耳散热：ΔT = $(round(ΔT_23, digits=2)) K")
        
        if abs(ΔT_23) > 1.0
            println("\n✓ 极耳强化散热效果明显（温度降低 > 1K）")
        else
            println("\n⚠️  极耳散热效果不明显，可能需要增大 h_tab")
        end
    else
        println("❌ 部分测试失败，无法进行完整对比")
    end
    
    println("\n" * "="^80)
    println("测试完成")
    println("="^80)
end

# 运行主函数
main()
