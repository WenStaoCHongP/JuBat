"""
充放电循环CZM损伤累积仿真示例

功能：
- 多SPMe并行电化学模型
- 二维分布式热模型
- 内聚力损伤模型（层间脱粘）
- 充放电循环：放电3600s → 静置600s → 充电3600s → 静置600s
- 损伤跨循环累积

循环参数：
- 循环次数：50次（可调整）
- 充放电倍率：1C (5A)
- 放电截止：2.5V 或 3600s
- 充电截止：4.2V 或 3600s
- 初始SOC：由Jellyroll参数决定（高SOC状态）
- 温度：每循环重置，循环内继承

输出：
- 容量衰减曲线
- 损伤演化曲线
- 库伦效率曲线
- 温度历史

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^70)
    println("充放电循环CZM损伤累积仿真")
    println("="^70)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/5] 参数设置...")
    
    # 电池参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    
    # 打印内聚力参数
    println("\n  内聚力参数：")
    @printf("    σ_max_n = %.1f MPa, δ_c_n = %.1f μm\n", 
            param_dim.cohesive.σ_max_n / 1e6, param_dim.cohesive.δ_c_n * 1e6)
    @printf("    τ_max_t = %.1f MPa, δ_c_t = %.1f μm\n",
            param_dim.cohesive.τ_max_t / 1e6, param_dim.cohesive.δ_c_t * 1e6)
    
    # 仿真选项
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "full"
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"
    opt.per_element_spme = true
    
    println("✓ 参数设置完成")
    
    # ========================================================================
    # 2. 创建网格
    # ========================================================================
    println("\n[2/5] 创建网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 热网格
    nθ = 60
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ 热网格创建完成")
    @printf("  单元数: %d, 节点数: %d\n", ne, nT)
    
    # CZM网格
    czm_mesh = JuBat.create_czm_mesh(mesh_th, param_dim; tol=1e-8)
    
    println("✓ CZM网格创建完成")
    @printf("  内聚力单元数: %d\n", czm_mesh.n_cohesive)
    
    # ========================================================================
    # 3. 设置循环参数
    # ========================================================================
    println("\n[3/5] 设置循环参数...")
    
    # 创建循环选项
    # 注意：为了演示，这里设置较少的循环次数
    # 实际仿真可以设置 n_cycles = 50
    # 循环顺序：放电 → 静置 → 充电 → 静置（适配Jellyroll初始高SOC状态）
    cycle_opt = JuBat.CycleOption(
        n_cycles = 5,           # 演示用5次循环，实际可设50次
        t_discharge = 3600.0,   # 放电时长 (s)
        t_rest1 = 600.0,        # 放电后静置 (s)
        t_charge = 3600.0,      # 充电时长 (s)
        t_rest2 = 600.0,        # 充电后静置 (s)
        I_discharge = 5.0,      # 放电电流 (A) - 1C
        I_charge = 5.0,         # 充电电流 (A) - 1C
        V_lower = 2.5,          # 放电截止电压
        V_upper = 4.2,          # 充电截止电压
        SOC_init = 0.8,         # 初始SOC（高SOC开始放电）
        reset_T_each_cycle = true,  # 每循环重置温度
        dt_cycle = [1.0, 10.0]  # 时间步长范围
    )
    
    println("✓ 循环参数设置完成")
    @printf("  循环次数: %d\n", cycle_opt.n_cycles)
    @printf("  充电: %.0fs, %.1fC\n", cycle_opt.t_charge, cycle_opt.I_charge/5.0)
    @printf("  放电: %.0fs, %.1fC\n", cycle_opt.t_discharge, cycle_opt.I_discharge/5.0)
    @printf("  静置: %.0fs + %.0fs\n", cycle_opt.t_rest1, cycle_opt.t_rest2)
    
    total_time = cycle_opt.n_cycles * (cycle_opt.t_charge + cycle_opt.t_rest1 + 
                                        cycle_opt.t_discharge + cycle_opt.t_rest2)
    @printf("  总仿真时间: %.1f小时\n", total_time / 3600.0)
    
    # ========================================================================
    # 4. 运行循环仿真
    # ========================================================================
    println("\n[4/5] 运行循环仿真...")
    
    result = nothing
    try
        result = JuBat.solve_cycling(case, cycle_opt, czm_mesh; 
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
        # 图1: 容量衰减和损伤演化（主要图）
        p1 = plot(layout=(2, 1), size=(800, 700))
        
        # 子图1: 容量
        plot!(p1[1], cycles, result.capacity_discharge,
              xlabel="", ylabel="Discharge Capacity (Ah)",
              label="Discharge", linewidth=2, marker=:circle, markersize=4,
              title="Capacity Fade and Damage Evolution")
        plot!(p1[1], cycles, result.capacity_charge,
              label="Charge", linewidth=2, marker=:square, markersize=4, linestyle=:dash)
        
        # 子图2: 损伤
        plot!(p1[2], cycles, result.D_max .* 100,
              xlabel="Cycle Number", ylabel="Damage (%)",
              label="D_max", linewidth=2.5, color=:red)
        plot!(p1[2], cycles, result.D_mean .* 100,
              label="D_mean", linewidth=2.5, color=:blue, linestyle=:dash)
        
        savefig(p1, "output/cycle_capacity_damage.png")
        println("  ✓ 保存: output/cycle_capacity_damage.png")
        
        # 图2: 温度和库伦效率
        p2 = plot(layout=(2, 1), size=(800, 600))
        
        plot!(p2[1], cycles, result.T_max,
              xlabel="", ylabel="T_max (K)",
              label="Max Temperature", linewidth=2, marker=:circle,
              title="Temperature and Coulombic Efficiency")
        
        plot!(p2[2], cycles, result.coulombic_efficiency,
              xlabel="Cycle Number", ylabel="CE (%)",
              label="Coulombic Efficiency", linewidth=2, marker=:circle)
        
        savefig(p2, "output/cycle_temperature_ce.png")
        println("  ✓ 保存: output/cycle_temperature_ce.png")
        
        # 图3: 综合汇总
        try
            p_all = JuBat.plot_cycling_results(result; save_path="output/")
        catch
            println("  ⚠ 综合图生成失败")
        end
    end
    
    # ========================================================================
    # 结果汇总
    # ========================================================================
    println("\n" * "="^70)
    println("仿真结果汇总")
    println("="^70)
    
    if length(cycles) > 0
        @printf("  完成循环数: %d\n", result.n_cycles)
        @printf("  初始放电容量: %.4f Ah\n", result.capacity_discharge[1])
        @printf("  最终放电容量: %.4f Ah\n", result.capacity_discharge[end])
        @printf("  容量保持率: %.1f%%\n", 
                100 * result.capacity_discharge[end] / result.capacity_discharge[1])
        @printf("  最终最大损伤: %.2f%%\n", result.D_max[end] * 100)
        @printf("  最终平均损伤: %.2f%%\n", result.D_mean[end] * 100)
        @printf("  最终断裂单元数: %d / %d\n", result.n_fractured[end], czm_mesh.n_cohesive)
    end
    
    println("\n生成的图像:")
    println("  1. output/cycle_capacity_damage.png - 容量衰减和损伤演化")
    println("  2. output/cycle_temperature_ce.png - 温度和库伦效率")
    println("  3. output/cycling_summary.png - 综合汇总图")
    
    println("\n" * "="^70)
    
    return result, czm_mesh
end

# ========================================================================
# 快速测试函数（单循环）
# ========================================================================
"""
    quick_test()

快速测试：只运行1个循环，用于验证代码正确性。
"""
function quick_test()
    println("="^50)
    println("快速测试模式（1循环）")
    println("="^50)
    
    param_dim = JuBat.ChooseCell("Jellyroll")
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 5; opt.Ns = 3; opt.Np = 5
    opt.Nrn = 5; opt.Nrp = 5
    opt.gsorder = 2
    opt.dimension = 1
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.per_element_spme = true
    opt.dtType = "auto"
    opt.solveType = "Crank-Nicolson"
    
    case = JuBat.SetCase(param_dim, opt)
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=30, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    czm_mesh = JuBat.create_czm_mesh(mesh_th, param_dim; tol=1e-8)
    
    # 短循环测试（放电 → 静置 → 充电 → 静置）
    cycle_opt = JuBat.CycleOption(
        n_cycles = 1,
        t_discharge = 60.0,   # 只放电60s
        t_rest1 = 10.0,
        t_charge = 60.0,      # 只充电60s
        t_rest2 = 10.0,
        I_discharge = 5.0,
        I_charge = 5.0,
        V_lower = 2.5,
        V_upper = 4.2,
        SOC_init = 0.8,       # 从80%开始放电
        dt_cycle = [1.0, 5.0]
    )
    
    println("开始快速测试...")
    result = JuBat.solve_cycling(case, cycle_opt, czm_mesh; verbose=true)
    
    println("\n✓ 快速测试完成")
    return result
end

# 运行主函数
if abspath(PROGRAM_FILE) == @__FILE__
    result, czm_mesh = main()
end
