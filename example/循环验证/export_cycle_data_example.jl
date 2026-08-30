"""
电化学-热耦合数据导出示例

功能：
- 运行一次电化学-热耦合仿真（不含CZM）
- 导出每个时间步的节点温度和SOC数据到CSV文件
- 输出的数据可供后续CZM仿真使用，加快全耦合仿真速度

循环参数：
- 循环次数：1次
- 充放电倍率：1C (5A)
- 放电：1800s 或 截止电压 2.5V
- 静置1：600s
- 充电：1800s 或 截止电压 4.2V
- 静置2：600s

输出文件：
- cycle_timesteps.csv: 时间步汇总（时间、电压、电流、温度、SOC）
- cycle_T_nodes.csv: 节点温度场历史
- cycle_soc_n.csv: 负极SOC场历史
- cycle_soc_p.csv: 正极SOC场历史
- cycle_mesh_nodes.csv: 网格节点坐标
- cycle_mesh_elements.csv: 网格单元连接

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^70)
    println("电化学-热耦合循环数据导出")
    println("="^70)
    println("\n此示例运行电化学-热仿真，并导出数据到CSV文件，")
    println("以便后续CZM仿真可以直接读取，加快仿真速度。")
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/4] 参数设置...")
    
    # 电池参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    
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
    opt.czm.enabled = false
    
    println("  ✓ 参数设置完成")
    
    # ========================================================================
    # 2. 创建网格
    # ========================================================================
    println("\n[2/4] 创建网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 热网格
    nθ = 60
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("  ✓ 热网格创建完成")
    @printf("  单元数: %d, 节点数: %d\n", ne, nT)
    
    # ========================================================================
    # 3. 设置循环参数
    # ========================================================================
    println("\n[3/4] 设置循环参数...")
    
    # 创建循环选项（单次循环）
    cycle_opt = JuBat.CycleOption(
        n_cycles = 1,            # 单次循环
        t_discharge = 1800.0,    # 放电时长 (s)
        t_rest1 = 600.0,         # 放电后静置 (s)
        t_charge = 1800.0,       # 充电时长 (s)
        t_rest2 = 600.0,         # 充电后静置 (s)
        I_discharge = 5.0,       # 放电电流 (A) - 1C
        I_charge = 5.0,          # 充电电流 (A) - 1C
        V_lower = 2.5,           # 放电截止电压
        V_upper = 4.2,           # 充电截止电压
        SOC_init = 0.65,         # 初始SOC
        reset_T_each_cycle = false,
        dt_cycle = [1.0, 10.0]   # 时间步长范围 [min, max] (s)
    )
    
    println("  ✓ 循环参数设置完成")
    @printf("  放电: %.0fs, %.1fC\n", cycle_opt.t_discharge, cycle_opt.I_discharge/5.0)
    @printf("  静置1: %.0fs\n", cycle_opt.t_rest1)
    @printf("  充电: %.0fs, %.1fC\n", cycle_opt.t_charge, cycle_opt.I_charge/5.0)
    @printf("  静置2: %.0fs\n", cycle_opt.t_rest2)
    
    total_time = cycle_opt.t_discharge + cycle_opt.t_rest1 + cycle_opt.t_charge + cycle_opt.t_rest2
    @printf("  总仿真时间: %.0fs (%.1f分钟)\n", total_time, total_time/60.0)
    
    # ========================================================================
    # 4. 运行仿真并导出数据
    # ========================================================================
    println("\n[4/4] 运行仿真并导出数据...")
    
    # 设置导出间隔（每N个时间步导出一次）
    # 设为1表示每个时间步都导出，设为10表示每10步导出一次
    export_interval = 5  # 每5步导出一次，平衡精度和文件大小
    
    export_data = nothing
    try
        export_data = JuBat.solve_cycle_with_export(
            case, cycle_opt;
            verbose=true,
            export_interval=export_interval
        )
        println("\n  ✓ 仿真完成")
    catch e
        println("\n  ✗ 仿真失败: $e")
        println("\n详细错误信息:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return nothing
    end
    
    if export_data === nothing
        println("仿真未返回结果")
        return nothing
    end
    
    # 导出到CSV
    println("\n导出数据到CSV文件...")
    output_dir = joinpath(@__DIR__, "..", "output", "export_cycle_data_example")
    
    files = JuBat.export_cycle_data_to_csv(export_data, output_dir; prefix="cycle")
    
    # ========================================================================
    # 结果汇总
    # ========================================================================
    println("\n" * "="^70)
    println("导出完成")
    println("="^70)
    
    n_steps = length(export_data.timesteps)
    @printf("\n  总时间步数: %d\n", n_steps)
    @printf("  网格信息: %d 单元, %d 节点\n", export_data.ne, export_data.nT)
    
    # 统计温度范围
    T_max_all = maximum(ts.T_max for ts in export_data.timesteps)
    T_min_all = minimum(ts.T_mean for ts in export_data.timesteps)
    @printf("  温度范围: %.2f K ~ %.2f K\n", T_min_all, T_max_all)
    
    # 统计SOC范围
    if !isempty(export_data.timesteps[1].soc_n)
        soc_min = minimum(minimum(ts.soc_n) for ts in export_data.timesteps if !isempty(ts.soc_n))
        soc_max = maximum(maximum(ts.soc_n) for ts in export_data.timesteps if !isempty(ts.soc_n))
        @printf("  负极SOC范围: %.4f ~ %.4f\n", soc_min, soc_max)
    end
    
    println("\n  输出文件目录: $output_dir")
    println("  文件列表:")
    for f in files
        if f !== nothing
            println("    - $(basename(f))")
        end
    end
    
    println("\n" * "="^70)
    println("使用说明")
    println("="^70)
    println("\n  导出的 CSV 用于外部可视化、统计分析和结果归档。")
    println("  JuBat 不再读取这些 CSV 来驱动后续求解。")
    
    return export_data
end

# 运行主函数
if abspath(PROGRAM_FILE) == @__FILE__
    export_data = main()
end
