"""
测试案例：单元放电截止电压逻辑验证

功能：
- 测试当有单元放电达到截止电压时的系统行为
- 验证是整个系统break还是单个单元电流归零
- 使用多SPMe并行架构（每个热单元对应独立的SPMe模型）
- 延长仿真时间以观察截止现象

测试目标：
1. 观察各单元SOC和电压的演化
2. 检测何时有单元达到截止电压
3. 验证截止后的电流分布变化
4. 对比整体系统break vs 单元电流归零的行为

日期：2025-01-23
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))

function main()
    println("="^80)
    println("截止电压逻辑测试：单元放电达到截止时的系统行为")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/7] 参数设置...")
    
    # 电池参数（Jellyroll结构）
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5  # 截止电压下限 (V)
    param_dim.cell.v_h = 4.2  # 截止电压上限 (V)
    
    # 仿真选项
    opt = JuBat.Option()
    
    # 电化学参数 - 使用较高的C-rate加速放电
    Crates = 2.0  # 2C放电，加速达到截止电压
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
    opt.mechanicalmodel = "full"
    
    # 时间设置 - 延长仿真时间以达到截止电压
    # 2C放电约需要30分钟（1800秒）完全放完
    opt.time = [0.0, 2000]  # 仿真时间 (s) - 足够长以观察截止
    opt.dt = [0.5, 10]      # 时间步长范围 [dt_min, dt_max] (s)
    opt.dtType = "auto"     # 自动时间步长
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"
    
    # ✨ 关键：启用多SPMe并行模式
    opt.per_element_spme = true
    
    println("✓ 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒 (约 %.1f 分钟)\n", opt.time[end], opt.time[end]/60)
    @printf("  截止电压: %.2f V\n", param_dim.cell.v_l)
    @printf("  模式: 多SPMe并行 ✨\n")
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/7] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    # 使用较少的单元数以便更清晰地观察单元级别的行为
    nθ = 40  # 周向单元数（减少以便观察）
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ Jellyroll网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 计算元素中心坐标和半径
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    θ_centers = atan.(centers[:,2], centers[:,1])
    
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
    # 预计算 layer_weights
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

    # 预计算热单元面积
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
    # 3. 运行求解器
    # ========================================================================
    println("\n[3/7] 运行多SPMe并行求解器...")
    println("  观察截止电压行为...")
    
    # 获取初始化的状态向量
    if case.opt.per_element_spme
        y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    else
        y0 = JuBat.ModelInitialisation_SimpleCoupling(case)
    end

    # 提取初始颗粒浓度
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen

    csn_init = y0[1:Nrn]
    csp_init = y0[(Nrn+1):(Nrn+Nrp)]

    theta_n = mean(csn_init)
    theta_p = mean(csp_init)

    @printf("  负极 theta = %.4f\n", theta_n)
    @printf("  正极 theta = %.4f\n", theta_p)

    # 计算理论 OCV
    U_p = Base.invokelatest(case.param_dim.PE.U, theta_p)
    U_n = Base.invokelatest(case.param_dim.NE.U, theta_n)
    OCV = U_p - U_n
    @printf("  初始 OCV = %.4f V\n", OCV)
    @printf("  截止电压 = %.4f V\n", param_dim.cell.v_l)
    
    result = nothing
    termination_reason = "unknown"
    
    try
        result = JuBat.Solve(case)
        
        # 检查终止原因
        if result !== nothing
            t = result["time [s]"]
            V = result["cell voltage [V]"]
            
            if V[end] <= param_dim.cell.v_l + 0.01
                termination_reason = "voltage_cutoff"
                println("✓ 求解因达到截止电压而终止")
            elseif t[end] >= opt.time[end] - 0.1
                termination_reason = "time_limit"
                println("✓ 求解因达到时间限制而终止")
            else
                termination_reason = "other"
                println("✓ 求解因其他原因终止")
            end
        end
        
    catch e
        println("✗ 求解失败: $e")
        rethrow(e)
    end
    
    if result === nothing
        error("Solve(case) did not return a result; aborting post-processing")
    end
    
    # ========================================================================
    # 4. 截止电压行为分析
    # ========================================================================
    println("\n[4/7] 截止电压行为分析...")
    println("="^60)
    
    # 基本变量
    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I_total = result["cell current [A]"]
    
    num_steps = length(t)
    println("  总时间步数: $num_steps")
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    @printf("  截止电压: %.4f V\n", param_dim.cell.v_l)
    @printf("  仿真结束时间: %.1f s\n", t[end])
    @printf("  设定结束时间: %.1f s\n", opt.time[end])
    
    # 判断终止原因
    println("\n  【终止原因分析】")
    if termination_reason == "voltage_cutoff"
        println("  ★ 系统因整体电压达到截止电压 ($(param_dim.cell.v_l) V) 而终止")
        println("    这表明当前实现是：整个系统 break")
    elseif termination_reason == "time_limit"
        println("  ★ 系统完成了全部仿真时间，未触发截止电压")
        println("    需要更长的仿真时间或更高的C-rate来触发截止")
    else
        println("  ★ 系统因其他原因终止")
    end
    
    # ========================================================================
    # 5. 单元级别分析
    # ========================================================================
    println("\n[5/7] 单元级别截止状态分析...")
    println("="^60)
    
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        I_scale = case.param.scale.I_typ
        
        println("\n  【逐单元电流演化分析】")
        
        # 找出电流变化最大的时刻
        n_snapshots = min(10, num_steps)
        idx_snapshots = round.(Int, range(1, num_steps, length=n_snapshots))
        
        println("\n  时间点采样分析：")
        println("  " * "-"^70)
        @printf("  %10s | %12s | %12s | %12s | %10s\n", 
                "时间(s)", "电流均值(A)", "电流标准差", "零电流单元", "电压(V)")
        println("  " * "-"^70)
        
        zero_current_threshold = 1e-6  # 判定为零电流的阈值
        
        for idx in idx_snapshots
            I_e_snap = I_e_hist[:, idx] .* I_scale
            n_zero = sum(abs.(I_e_snap) .< zero_current_threshold)
            
            @printf("  %10.1f | %12.4e | %12.4e | %12d | %10.4f\n",
                    t[idx], mean(I_e_snap), std(I_e_snap), n_zero, V[idx])
        end
        println("  " * "-"^70)
        
        # 检查是否有单元电流归零但系统继续运行
        println("\n  【单元电流归零检测】")
        
        # 在每个时间步检查是否有电流接近零的单元
        cutoff_detected_step = -1
        cutoff_element_count = 0
        
        for step in 1:num_steps
            I_e_step = I_e_hist[:, step] .* I_scale
            n_near_zero = sum(abs.(I_e_step) .< zero_current_threshold)
            
            if n_near_zero > 0 && cutoff_detected_step < 0
                cutoff_detected_step = step
                cutoff_element_count = n_near_zero
            end
        end
        
        if cutoff_detected_step > 0
            println("  ★ 首次检测到单元电流归零:")
            @printf("    时间步: %d (t = %.1f s)\n", cutoff_detected_step, t[cutoff_detected_step])
            @printf("    零电流单元数: %d / %d\n", cutoff_element_count, ne)
            @printf("    此时系统电压: %.4f V\n", V[cutoff_detected_step])
            
            # 检查系统是否继续运行
            if cutoff_detected_step < num_steps
                remaining_steps = num_steps - cutoff_detected_step
                @printf("    系统在检测到单元截止后继续运行了 %d 步\n", remaining_steps)
                println("    → 这表明实现支持：单元电流归零而系统继续运行")
            end
        else
            println("  未检测到任何单元电流归零的情况")
            println("  可能原因：")
            println("    1. 所有单元同时达到截止电压")
            println("    2. 系统在单元截止前就因整体电压截止而终止")
        end
        
        # 检查截止相关的变量
        if haskey(result, "thermal2D n_cutoff_elements")
            n_cutoff_hist = result["thermal2D n_cutoff_elements"]
            println("\n  【截止单元数历史】")
            println("  注：如果此数据存在，表明系统支持单元级别截止追踪")
            
            max_cutoff = maximum(n_cutoff_hist)
            @printf("  最大截止单元数: %d\n", max_cutoff)
            
            if max_cutoff > 0
                first_cutoff_idx = findfirst(n_cutoff_hist .> 0)
                @printf("  首次出现截止: t = %.1f s\n", t[first_cutoff_idx])
            end
        end
        
        if haskey(result, "thermal2D active_mask")
            active_mask_hist = result["thermal2D active_mask"]
            n_active_final = sum(active_mask_hist[:, end])
            @printf("\n  最终活跃单元数: %d / %d\n", n_active_final, ne)
        end
        
    else
        println("  ⚠ 未找到 thermal2D element current 数据")
    end
    
    # SOC分析
    if haskey(result, "thermal2D element soc_n") && haskey(result, "thermal2D element soc_p")
        soc_n_hist = result["thermal2D element soc_n"]
        soc_p_hist = result["thermal2D element soc_p"]
        
        println("\n  【单元SOC分布分析】（最终时刻）")
        soc_n_final = soc_n_hist[:, end]
        soc_p_final = soc_p_hist[:, end]
        
        @printf("  负极SOC: 均值=%.4f, 范围=[%.4f, %.4f]\n", 
                mean(soc_n_final), minimum(soc_n_final), maximum(soc_n_final))
        @printf("  正极SOC: 均值=%.4f, 范围=[%.4f, %.4f]\n", 
                mean(soc_p_final), minimum(soc_p_final), maximum(soc_p_final))
        
        # SOC非均匀性分析
        soc_spread_n = maximum(soc_n_final) - minimum(soc_n_final)
        soc_spread_p = maximum(soc_p_final) - minimum(soc_p_final)
        
        @printf("  SOC分布差异: 负极=%.4f, 正极=%.4f\n", soc_spread_n, soc_spread_p)
        
        if soc_spread_n > 0.01 || soc_spread_p > 0.01
            println("  → SOC存在明显的空间不均匀性，这可能导致某些单元先达到截止")
        end
    end
    
    # ========================================================================
    # 6. 绘图
    # ========================================================================
    println("\n[6/7] 生成诊断图像...")
    
    # 图1：电压-时间曲线（标注截止电压）
    p1 = plot(t, V, xlabel="Time (s)", ylabel="Voltage (V)", 
              label="Cell Voltage", linewidth=2, title="Discharge Curve with Cutoff Analysis",
              legend=:topright)
    hline!([param_dim.cell.v_l], label="Cutoff Voltage ($(param_dim.cell.v_l) V)", 
           linestyle=:dash, color=:red, linewidth=2)
    
    # 标注终止点
    scatter!([t[end]], [V[end]], label="Termination Point", 
             markersize=8, color=:red, markershape=:star5)
    
    savefig(p1, "testexample_cutoff_voltage.png")
    println("  ✓ 保存: testexample_cutoff_voltage.png")
    
    # 图2：逐单元电流演化
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        I_scale = case.param.scale.I_typ
        
        # 选择几个代表性单元
        sample_elements = [1, round(Int, ne/4), round(Int, ne/2), round(Int, 3*ne/4), ne]
        
        p2 = plot(xlabel="Time (s)", ylabel="Element Current (A)", 
                  title="Element Current Evolution",
                  legend=:topright)
        
        for e in sample_elements
            plot!(p2, t, I_e_hist[e, :] .* I_scale, 
                  label="Element $e", linewidth=1.5)
        end
        
        hline!([0.0], label="Zero Current", linestyle=:dash, color=:black)
        
        savefig(p2, "testexample_element_currents.png")
        println("  ✓ 保存: testexample_element_currents.png")
        
        # 图3：电流分布热图
        p3 = heatmap(t, 1:ne, I_e_hist .* I_scale,
                     xlabel="Time (s)", ylabel="Element Index",
                     title="Element Current Distribution (A)",
                     color=:viridis)
        
        savefig(p3, "testexample_current_heatmap.png")
        println("  ✓ 保存: testexample_current_heatmap.png")
    end
    
    # 图4：SOC演化
    if haskey(result, "thermal2D element soc_n") && haskey(result, "thermal2D element soc_p")
        soc_n_hist = result["thermal2D element soc_n"]
        soc_p_hist = result["thermal2D element soc_p"]
        
        # SOC均值和范围
        soc_n_mean = [mean(soc_n_hist[:, i]) for i in 1:num_steps]
        soc_n_min = [minimum(soc_n_hist[:, i]) for i in 1:num_steps]
        soc_n_max = [maximum(soc_n_hist[:, i]) for i in 1:num_steps]
        
        soc_p_mean = [mean(soc_p_hist[:, i]) for i in 1:num_steps]
        soc_p_min = [minimum(soc_p_hist[:, i]) for i in 1:num_steps]
        soc_p_max = [maximum(soc_p_hist[:, i]) for i in 1:num_steps]
        
        p4 = plot(layout=(2,1), size=(800, 600))
        
        # 负极SOC
        plot!(p4[1], t, soc_n_mean, ribbon=(soc_n_mean .- soc_n_min, soc_n_max .- soc_n_mean),
              label="Negative Electrode", fillalpha=0.3, linewidth=2,
              xlabel="", ylabel="SOC (θ)", title="Negative Electrode SOC")
        
        # 正极SOC
        plot!(p4[2], t, soc_p_mean, ribbon=(soc_p_mean .- soc_p_min, soc_p_max .- soc_p_mean),
              label="Positive Electrode", fillalpha=0.3, linewidth=2, color=:red,
              xlabel="Time (s)", ylabel="SOC (θ)", title="Positive Electrode SOC")
        
        savefig(p4, "testexample_soc_evolution.png")
        println("  ✓ 保存: testexample_soc_evolution.png")
    end
    
    # 图5：电流异质性演化
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        
        cv_I = [std(I_e_hist[:, i]) / max(abs(mean(I_e_hist[:, i])), 1e-10) for i in 1:num_steps]
        
        p5 = plot(t, cv_I .* 100, xlabel="Time (s)", ylabel="CV of Current (%)", 
                  label="Current Heterogeneity", linewidth=2, 
                  title="Current Distribution Heterogeneity")
        
        savefig(p5, "testexample_current_heterogeneity.png")
        println("  ✓ 保存: testexample_current_heterogeneity.png")
    end
    
    # ========================================================================
    # 7. 总结
    # ========================================================================
    println("\n[7/7] 测试总结")
    println("="^80)
    
    println("""
    【截止电压逻辑测试结果】
    
    1. 仿真概况：
       - 总时间步数: $num_steps
       - 初始电压: $(V[1]) V
       - 最终电压: $(V[end]) V
       - 截止电压: $(param_dim.cell.v_l) V
       - 终止原因: $termination_reason
    """)
    
    # 关键结论
    println("    2. 关键发现：")
    
    if termination_reason == "voltage_cutoff"
        println("       ★ 系统在整体电压达到截止电压时终止 (break)")
        
        # 检查是否有单元级别的差异
        if haskey(result, "thermal2D element current")
            I_e_final = result["thermal2D element current"][:, end]
            cv_final = std(I_e_final) / max(abs(mean(I_e_final)), 1e-10)
            
            if cv_final > 0.1
                println("       ★ 最终时刻单元电流存在明显差异 (CV=$(round(cv_final*100, digits=1))%)")
                println("         这表明不同单元的放电深度不同")
            else
                println("       ★ 最终时刻单元电流分布较均匀 (CV=$(round(cv_final*100, digits=1))%)")
            end
        end
    end
    
    println("""
    
    3. 当前实现行为：
       - 当整体电池电压 < v_l 时：整个系统 break（Solve.jl:264-266）
       - 单元级截止检测：由 _detect_cutoff_elements 函数实现
       - 达到截止的单元：电流设为0，但需系统继续运行才能观察到
    
    4. 生成的图像：
       - testexample_cutoff_voltage.png   - 放电曲线与截止电压标注
       - testexample_element_currents.png - 代表性单元电流演化
       - testexample_current_heatmap.png  - 全部单元电流热图
       - testexample_soc_evolution.png    - SOC演化（带范围）
       - testexample_current_heterogeneity.png - 电流异质性演化
    """)
    
    println("="^80)
    println("测试完成")
end

# 运行主函数
main()
