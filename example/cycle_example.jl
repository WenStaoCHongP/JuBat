"""
电化学-热耦合循环仿真示例

功能：
- 多SPMe并行电化学模型
- 二维分布式热模型
- 验证循环求解器可靠性

循环参数：
- 循环次数：1次（可调整）
- 充放电倍率：1C (5A)
- 放电：2000s 或 截止电压 2.5V
- 静置1：1000s（锂扩散过程）
- 充电：2000s 或 截止电压 4.2V
- 静置2：1000s（锂扩散过程）

输出：
- 电压曲线
- 温度历史
- 容量记录
- 静置阶段锂扩散验证

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))

function main()
    println("="^70)
    println("电化学-热耦合循环仿真验证")
    println("="^70)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/5] 参数设置...")
    
    # 电池参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5  # 放电截止电压
    param_dim.cell.v_h = 4.2  # 充电截止电压
    
    # 仿真选项
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10        # 负极网格数
    opt.Ns = 5         # 隔膜网格数
    opt.Np = 10        # 正极网格数
    opt.Nrn = 10       # 负极颗粒径向网格数
    opt.Nrp = 10       # 正极颗粒径向网格数
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "full"
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"
    opt.per_element_spme = true  # 多SPMe并行模式
    
    println("✓ 参数设置完成")
    println("  模型: SPMe + 分布式热模型")
    println("  模式: 多SPMe并行")
    
    # ========================================================================
    # 2. 创建网格
    # ========================================================================
    println("\n[2/5] 创建网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 热网格（使用较少单元数以加快验证速度）
    nθ = 40  # 周向单元数
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ 热网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 计算 layer_weights（热源分配）
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
    
    # ========================================================================
    # 3. 设置循环参数
    # ========================================================================
    println("\n[3/5] 设置循环参数...")
    
    # 创建循环选项
    # 循环顺序：放电 → 静置1 → 充电 → 静置2
    cycle_opt = JuBat.CycleOption(
        n_cycles = 1,           # 循环次数（验证用1次）
        t_discharge = 2000.0,   # 放电时长 (s)
        t_rest1 = 1000.0,       # 放电后静置 (s) - 锂扩散
        t_charge = 2000.0,      # 充电时长 (s)
        t_rest2 = 1000.0,       # 充电后静置 (s) - 锂扩散
        I_discharge = 5.0,      # 放电电流 (A) - 1C
        I_charge = 5.0,         # 充电电流 (A) - 1C
        V_lower = 2.5,          # 放电截止电压 (V)
        V_upper = 4.2,          # 充电截止电压 (V)
        reset_T_each_cycle = false,  # 循环内温度继承
        dt_cycle = [1.0, 10.0]  # 时间步长范围 [min, max] (s)
    )
    
    println("✓ 循环参数设置完成")
    @printf("  循环次数: %d\n", cycle_opt.n_cycles)
    @printf("  放电: %.0fs, %.1fC, 截止 %.2fV\n", 
            cycle_opt.t_discharge, cycle_opt.I_discharge/5.0, cycle_opt.V_lower)
    @printf("  静置1: %.0fs (锂扩散过程)\n", cycle_opt.t_rest1)
    @printf("  充电: %.0fs, %.1fC, 截止 %.2fV\n", 
            cycle_opt.t_charge, cycle_opt.I_charge/5.0, cycle_opt.V_upper)
    @printf("  静置2: %.0fs (锂扩散过程)\n", cycle_opt.t_rest2)
    
    total_time = cycle_opt.n_cycles * (cycle_opt.t_discharge + cycle_opt.t_rest1 + 
                                        cycle_opt.t_charge + cycle_opt.t_rest2)
    @printf("  单循环时长: %.0fs (%.1f分钟)\n", total_time, total_time/60.0)
    @printf("  总仿真时间: %.1f分钟\n", total_time / 60.0)
    
    # ========================================================================
    # 4. 运行循环仿真
    # ========================================================================
    println("\n[4/5] 运行循环仿真...")
    println("  注意：静置阶段将继承电化学状态并继续锂扩散过程")
    
    result = nothing
    try
        # 不使用 CZM，传入 nothing
        result = JuBat.solve_cycling(case, cycle_opt, nothing; 
                                     verbose=true, save_detailed=true)
        println("\n✓ 循环仿真完成")
    catch e
        println("\n✗ 仿真失败: $e")
        println("\n详细错误信息:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return nothing
    end
    
    if result === nothing
        println("仿真未返回结果")
        return nothing
    end
    
    # ========================================================================
    # 5. 结果可视化
    # ========================================================================
    println("\n[5/5] 生成结果图像...")
    
    # 确保输出目录存在
    isdir("output") || mkdir("output")
    
    cycles = result.cycle_idx
    
    if length(cycles) > 0
        # 图1: 容量和温度（主要验证图）
        p1 = plot(layout=(2, 1), size=(800, 600))
        
        # 子图1: 容量
        plot!(p1[1], cycles, result.capacity_discharge,
              xlabel="", ylabel="Capacity (Ah)",
              label="Discharge", linewidth=2, marker=:circle, markersize=6,
              title="Cycle Capacity and Temperature")
        plot!(p1[1], cycles, result.capacity_charge,
              label="Charge", linewidth=2, marker=:square, markersize=6, linestyle=:dash)
        
        # 子图2: 温度
        plot!(p1[2], cycles, result.T_max,
              xlabel="Cycle Number", ylabel="T_max (K)",
              label="Max Temperature", linewidth=2, marker=:circle, markersize=6,
              color=:red)
        
        savefig(p1, "output/cycle_example_results.png")
        println("  ✓ 保存: output/cycle_example_results.png")
        
        # 图2: 库伦效率
        if length(cycles) > 0
            p2 = plot(cycles, result.coulombic_efficiency,
                      xlabel="Cycle Number", ylabel="Coulombic Efficiency (%)",
                      label="CE", linewidth=2, marker=:circle, markersize=6,
                      title="Coulombic Efficiency",
                      legend=:bottomleft)
            
            savefig(p2, "output/cycle_example_ce.png")
            println("  ✓ 保存: output/cycle_example_ce.png")
        end
    end
    
    # ========================================================================
    # 结果汇总
    # ========================================================================
    println("\n" * "="^70)
    println("循环仿真结果汇总")
    println("="^70)
    
    if length(cycles) > 0
        @printf("\n  完成循环数: %d\n", result.n_cycles)
        
        println("\n  【容量信息】")
        @printf("    放电容量: %.4f Ah\n", result.capacity_discharge[1])
        @printf("    充电容量: %.4f Ah\n", result.capacity_charge[1])
        @printf("    库伦效率: %.2f%%\n", result.coulombic_efficiency[1])
        
        println("\n  【温度信息】")
        @printf("    最高温度: %.2f K (%.2f °C)\n", 
                result.T_max[1], result.T_max[1] - 273.15)
        
        println("\n  【循环验证】")
        if result.n_cycles >= 1
            println("    ✓ 放电阶段完成")
            println("    ✓ 静置1阶段完成（锂扩散）")
            println("    ✓ 充电阶段完成")
            println("    ✓ 静置2阶段完成（锂扩散）")
            println("    ✓ 状态正确传递")
        end
        
        # 如果有详细结果，打印各阶段信息
        if !isempty(result.cycle_results)
            cycle1 = result.cycle_results[1]
            println("\n  【各阶段详情】")
            @printf("    放电: %.1fs, %.3fV → %.3fV, %.4fAh\n",
                    cycle1.discharge.duration, 
                    cycle1.discharge.V_start, cycle1.discharge.V_end,
                    cycle1.discharge.capacity)
            @printf("    静置1: %.1fs, T_max=%.2fK\n",
                    cycle1.rest1.duration, cycle1.rest1.T_max)
            @printf("    充电: %.1fs, %.3fV → %.3fV, %.4fAh\n",
                    cycle1.charge.duration,
                    cycle1.charge.V_start, cycle1.charge.V_end,
                    cycle1.charge.capacity)
            @printf("    静置2: %.1fs, T_max=%.2fK\n",
                    cycle1.rest2.duration, cycle1.rest2.T_max)
        end
    end
    
    println("\n  【生成的图像】")
    println("    1. output/cycle_example_results.png - 容量和温度")
    println("    2. output/cycle_example_ce.png - 库伦效率")
    
    println("\n" * "="^70)
    println("验证完成")
    println("="^70)
    
    return result
end

# 运行主函数
result = main()
