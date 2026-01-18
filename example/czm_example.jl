"""
内聚力模型(CZM)示例：电化学-热-力-损伤耦合仿真

功能：
- 多SPMe并行电化学模型
- 二维分布式热模型
- 内聚力损伤模型（层间脱粘）
- 热应变和化学应变耦合
- 损伤累积与断裂判据

输出图像：
- 图1：最大/平均损伤变量随时间变化曲线
- 图2：温度场和损伤场分布

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("内聚力模型(CZM)电化学-热-力-损伤耦合仿真")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/8] 参数设置...")
    
    # 电池参数（Jellyroll结构，包含内聚力参数）
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5  # 截止电压下限 (V)
    param_dim.cell.v_h = 4.2  # 截止电压上限 (V)
    
    # 打印内聚力参数
    println("\n  内聚力参数（有量纲）：")
    @printf("    法向最大牵引力 σ_max_n = %.1f MPa\n", param_dim.cohesive.σ_max_n / 1e6)
    @printf("    法向临界位移 δ_c_n = %.1f μm\n", param_dim.cohesive.δ_c_n * 1e6)
    @printf("    法向断裂能 G_c_n = %.1f J/m²\n", param_dim.cohesive.G_c_n)
    @printf("    切向最大牵引力 τ_max_t = %.1f MPa\n", param_dim.cohesive.τ_max_t / 1e6)
    @printf("    BK指数 η = %.2f\n", param_dim.cohesive.eta)
    
    # 仿真选项
    opt = JuBat.Option()
    
    # 电化学参数
    Crates = 2.0  # 高倍率以产生更大应力
    i = 5 * Crates  # 电流 (A)
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "full"
    
    # 时间设置（较短时间用于演示）
    opt.time = [0.0, 120]  # 仿真时间 (s)
    opt.dt = [0.5, 5]      # 时间步长范围 (s)
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"
    
    # 启用多SPMe并行模式
    opt.per_element_spme = true
    
    println("✓ 参数设置完成")
    @printf("  电流: %.2f A (%.1f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/8] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    nθ = 60  # 周向单元数（适中分辨率）
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ 热网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 几何参数
    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
    # ========================================================================
    # 3. 创建内聚力网格
    # ========================================================================
    println("\n[3/8] 创建内聚力网格...")
    
    czm_mesh = JuBat.create_czm_mesh(mesh_th, param_dim; tol=1e-8)
    
    println("✓ 内聚力网格创建完成")
    @printf("  扩展后节点数: %d\n", czm_mesh.nnode)
    @printf("  内聚力单元数: %d\n", czm_mesh.n_cohesive)
    @printf("  层数: %d\n", czm_mesh.n_layers)
    @printf("  层间界面数: %d\n", czm_mesh.n_layers - 1)
    
    if czm_mesh.n_cohesive == 0
        println("\n⚠ 警告：没有检测到内聚力单元。")
        println("  这可能是因为热网格只覆盖一圈，没有层间界面。")
        println("  继续运行电化学-热仿真，但跳过损伤计算。")
    end
    
    # ========================================================================
    # 4. 运行电化学-热求解器
    # ========================================================================
    println("\n[4/8] 运行多SPMe并行求解器...")
    
    result = nothing
    try
        result = JuBat.Solve(case)
        println("✓ 求解成功完成")
    catch e
        println("✗ 求解失败: $e")
        rethrow(e)
    end
    
    if result === nothing
        error("求解未返回结果")
    end
    
    # 提取基本结果
    t = result["time [s]"]
    V = result["cell voltage [V]"]
    num_steps = length(t)
    
    println("✓ 结果提取完成")
    @printf("  总时间步数: %d\n", num_steps)
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    
    # ========================================================================
    # 5. 计算每个时间步的CZM损伤
    # ========================================================================
    println("\n[5/8] 计算时间历程损伤场...")
    
    # 初始化损伤历史数组
    damage_max_hist = zeros(Float64, num_steps)
    damage_mean_hist = zeros(Float64, num_steps)
    damage_fractured_hist = zeros(Int64, num_steps)
    
    # 材料参数
    E_eff = 0.5 * (param_dim.NE.E + param_dim.PE.E)
    ν_eff = 0.5 * (param_dim.NE.nu + param_dim.PE.nu)
    α_eff = 0.5 * (param_dim.NE.alphaT + param_dim.PE.alphaT)
    β_n = param_dim.NE.Omega / 3.0
    β_p = param_dim.PE.Omega / 3.0
    
    @printf("  有效杨氏模量 E = %.2f GPa\n", E_eff / 1e9)
    @printf("  有效泊松比 ν = %.2f\n", ν_eff)
    @printf("  有效热膨胀系数 α = %.2e /K\n", α_eff)
    
    # 只有存在内聚力单元时才计算损伤
    if czm_mesh.n_cohesive > 0
        # 获取SOC和温度历史
        # 尝试两个可能的键名
        T_key = haskey(result, "thermal2D temperature [K]") ? "thermal2D temperature [K]" :
                (haskey(result, "thermal2D T_nodes [K]") ? "thermal2D T_nodes [K]" : nothing)
        
        has_data = haskey(result, "thermal2D element soc_n") && 
                   haskey(result, "thermal2D element soc_p") &&
                   T_key !== nothing
        
        if has_data
            soc_n_hist = result["thermal2D element soc_n"]
            soc_p_hist = result["thermal2D element soc_p"]
            T_nodes_hist_K = result[T_key]
            
            # 参考值
            soc_ref_n = case.param.NE.cs0
            soc_ref_p = case.param.PE.cs0
            T_ref = param_dim.scale.T_ref
            T0 = param_dim.cell.T0
            
            # 初始化位移
            ndof = 2 * czm_mesh.nnode
            u_prev = zeros(Float64, ndof)
            F_ext = zeros(Float64, ndof)  # 无外部机械载荷
            
            # 重置损伤状态
            JuBat.reset_damage_states!(czm_mesh)
            
            println("  计算 $(num_steps) 个时间步的损伤...")
            
            for step in 1:num_steps
                try
                    # 计算温度变化 - 处理不同的数据维度
                    if ndims(T_nodes_hist_K) == 1
                        T_nodes_K = T_nodes_hist_K  # 单一时刻数据
                    else
                        T_nodes_K = T_nodes_hist_K[:, step]
                    end
                    
                    # 扩展温度到CZM网格节点
                    # CZM网格的前nT个节点对应原热网格节点
                    T_czm = zeros(Float64, czm_mesh.nnode)
                    nT_actual = min(nT, length(T_nodes_K))
                    T_czm[1:nT_actual] = T_nodes_K[1:nT_actual]
                    # 复制的节点使用相同温度
                    for (orig, new_nodes) in czm_mesh.node_map
                        if orig <= nT_actual
                            for new_n in new_nodes
                                if new_n > nT_actual && new_n <= czm_mesh.nnode
                                    T_czm[new_n] = T_nodes_K[orig]
                                end
                            end
                        end
                    end
                    
                    # 计算单元温度变化
                    dT_elem = zeros(Float64, ne)
                    for e in 1:ne
                        nodes = czm_mesh.bulk_element[e, :]
                        T_avg = mean(T_czm[nodes])
                        dT_elem[e] = T_avg - T0
                    end
                    
                    # SOC变化 - 处理不同的数据维度
                    if ndims(soc_n_hist) == 1
                        Δsoc_n_elem = soc_n_hist .- soc_ref_n
                        Δsoc_p_elem = soc_p_hist .- soc_ref_p
                    else
                        Δsoc_n_elem = soc_n_hist[:, step] .- soc_ref_n
                        Δsoc_p_elem = soc_p_hist[:, step] .- soc_ref_p
                    end
                    
                    # 求解CZM系统
                    czm_result = JuBat.solve_czm_step(
                        czm_mesh, F_ext, E_eff, ν_eff,
                        param_dim.cohesive, param_dim, u_prev;
                        α_eff=α_eff, β_n=β_n, β_p=β_p,
                        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
                        max_iter=30, tol=1e-6)
                    
                    if czm_result.converged
                        u_prev = czm_result.displacement
                    end
                    
                    # 记录损伤统计
                    stats = JuBat.get_damage_statistics(czm_mesh)
                    damage_max_hist[step] = stats.max_D
                    damage_mean_hist[step] = stats.mean_D
                    damage_fractured_hist[step] = stats.n_fractured
                    
                    # 进度输出
                    if step % 10 == 0 || step == num_steps
                        @printf("    步 %d/%d: t=%.1fs, D_max=%.2f%%, D_mean=%.2f%%\n",
                                step, num_steps, t[step], 
                                stats.max_D * 100, stats.mean_D * 100)
                    end
                    
                catch e
                    @warn "时间步 $step 损伤计算失败" exception=(e, catch_backtrace())
                    damage_max_hist[step] = step > 1 ? damage_max_hist[step-1] : 0.0
                    damage_mean_hist[step] = step > 1 ? damage_mean_hist[step-1] : 0.0
                end
            end
            
            println("\n✓ 损伤历史计算完成")
            @printf("  最终最大损伤: %.2f%%\n", damage_max_hist[end] * 100)
            @printf("  最终平均损伤: %.2f%%\n", damage_mean_hist[end] * 100)
            @printf("  断裂单元数: %d / %d\n", damage_fractured_hist[end], czm_mesh.n_cohesive)
        else
            println("  ⚠ 未找到SOC或温度历史数据，跳过损伤计算")
        end
    else
        println("  ⚠ 无内聚力单元，跳过损伤计算")
    end
    
    # ========================================================================
    # 6. 图1：损伤变量随时间变化曲线
    # ========================================================================
    println("\n[6/8] 生成图1：损伤变量演化曲线...")
    
    # 确保输出目录存在
    isdir("output") || mkdir("output")
    
    if czm_mesh.n_cohesive > 0 && any(damage_max_hist .> 0)
        p1 = plot(size=(800, 500), dpi=150)
        
        # 最大损伤曲线
        plot!(p1, t, damage_max_hist .* 100,
              label="Maximum Damage D_max",
              linewidth=2.5, color=:red,
              xlabel="Time (s)", ylabel="Damage (%)",
              title="Damage Evolution During Discharge ($(Crates)C)")
        
        # 平均损伤曲线
        plot!(p1, t, damage_mean_hist .* 100,
              label="Average Damage D_mean",
              linewidth=2.5, color=:blue, linestyle=:dash)
        
        # 添加断裂阈值参考线
        hline!(p1, [99.0], label="Fracture threshold (99%)", 
               linestyle=:dot, color=:black, linewidth=1.5)
        
        # 添加网格
        plot!(p1, grid=true, gridalpha=0.3)
        
        # 添加图例
        plot!(p1, legend=:topleft)
        
        savefig(p1, "output/czm_damage_evolution.png")
        println("  ✓ 保存: output/czm_damage_evolution.png")
        
        # 保存SVG版本
        try
            savefig(p1, "output/czm_damage_evolution.svg")
            println("  ✓ 保存: output/czm_damage_evolution.svg")
        catch
        end
    else
        println("  ⚠ 无损伤数据，跳过图1")
    end
    
    # ========================================================================
    # 7. 图2：温度场和损伤场分布
    # ========================================================================
    println("\n[7/8] 生成图2：温度场和损伤场分布...")
    
    # 计算单元中心坐标（热网格）
    x_elem = zeros(Float64, ne)
    y_elem = zeros(Float64, ne)
    for e in 1:ne
        nodes = mesh_th.element[e, :]
        x_elem[e] = mean(mesh_th.node[nodes, 1])
        y_elem[e] = mean(mesh_th.node[nodes, 2])
    end
    
    # 温度场 - 处理不同的数据格式
    T_final = nothing
    T_nodes_K = nothing
    
    # 优先使用 T_nodes 数据（最终时刻）
    if haskey(result, "thermal2D T_nodes [K]")
        T_data = result["thermal2D T_nodes [K]"]
        if ndims(T_data) == 1
            T_nodes_K = T_data
        else
            T_nodes_K = T_data[:, end]
        end
    elseif haskey(result, "thermal2D temperature [K]")
        T_data = result["thermal2D temperature [K]"]
        if ndims(T_data) == 1
            T_nodes_K = T_data
        else
            T_nodes_K = T_data[:, end]
        end
    end
    
    if T_nodes_K !== nothing
        T_elem = zeros(Float64, ne)
        for e in 1:ne
            nodes = mesh_th.element[e, :]
            # 确保节点索引有效
            valid_nodes = [n for n in nodes if n <= length(T_nodes_K)]
            if !isempty(valid_nodes)
                T_elem[e] = mean(T_nodes_K[valid_nodes])
            end
        end
        T_final = T_elem
    end
    
    # 损伤场（映射到单元）
    D_elem = zeros(Float64, ne)
    if czm_mesh.n_cohesive > 0
        # 计算每个热单元关联的内聚力单元的平均损伤
        coh_count = zeros(Int64, ne)
        for coh_elem in czm_mesh.cohesive_elements
            # 找到关联的热单元（通过节点位置）
            x_coh = mean([czm_mesh.node[n, 1] for n in coh_elem.nodes])
            y_coh = mean([czm_mesh.node[n, 2] for n in coh_elem.nodes])
            
            # 找最近的热单元
            min_dist = Inf
            nearest_e = 1
            for e in 1:ne
                dist = hypot(x_elem[e] - x_coh, y_elem[e] - y_coh)
                if dist < min_dist
                    min_dist = dist
                    nearest_e = e
                end
            end
            
            D_elem[nearest_e] += czm_mesh.damage_states[coh_elem.id].D
            coh_count[nearest_e] += 1
        end
        
        # 计算平均
        for e in 1:ne
            if coh_count[e] > 0
                D_elem[e] /= coh_count[e]
            end
        end
    end
    
    # 创建组合图
    p2 = plot(layout=(1, 2), size=(1400, 600), dpi=150)
    
    # 子图1：温度场
    if T_final !== nothing
        T_min, T_max = extrema(T_final)
        scatter!(p2[1], x_elem .* 1000, y_elem .* 1000,
                marker_z=T_final,
                color=:inferno, markersize=4, markerstrokewidth=0,
                xlabel="x (mm)", ylabel="y (mm)",
                title="Temperature Field (t=$(t[end])s)",
                colorbar=true, colorbar_title="T (K)",
                aspect_ratio=:equal,
                clims=(T_min, T_max))
        
        @printf("  温度范围: [%.2f, %.2f] K\n", T_min, T_max)
    else
        annotate!(p2[1], 0.5, 0.5, text("No temperature data", 12))
    end
    
    # 子图2：损伤场
    if czm_mesh.n_cohesive > 0 && maximum(D_elem) > 0
        scatter!(p2[2], x_elem .* 1000, y_elem .* 1000,
                marker_z=D_elem .* 100,
                color=:hot, markersize=4, markerstrokewidth=0,
                xlabel="x (mm)", ylabel="y (mm)",
                title="Damage Field (t=$(t[end])s)",
                colorbar=true, colorbar_title="D (%)",
                aspect_ratio=:equal,
                clims=(0, max(maximum(D_elem) * 100, 1)))
        
        @printf("  损伤范围: [%.2f, %.2f]%%\n", minimum(D_elem) * 100, maximum(D_elem) * 100)
    else
        scatter!(p2[2], x_elem .* 1000, y_elem .* 1000,
                marker_z=zeros(ne),
                color=:hot, markersize=4, markerstrokewidth=0,
                xlabel="x (mm)", ylabel="y (mm)",
                title="Damage Field (No damage)",
                colorbar=true, colorbar_title="D (%)",
                aspect_ratio=:equal,
                clims=(0, 1))
    end
    
    savefig(p2, "output/czm_temperature_damage_field.png")
    println("  ✓ 保存: output/czm_temperature_damage_field.png")
    
    try
        savefig(p2, "output/czm_temperature_damage_field.svg")
        println("  ✓ 保存: output/czm_temperature_damage_field.svg")
    catch
    end
    
    # ========================================================================
    # 8. 额外图像：电压和温度演化
    # ========================================================================
    println("\n[8/8] 生成额外图像...")
    
    # 电压曲线
    p_v = plot(t, V, xlabel="Time (s)", ylabel="Voltage (V)",
              label="Cell Voltage", linewidth=2, 
              title="Discharge Curve ($(Crates)C)")
    hline!([param_dim.cell.v_l], label="Cutoff", linestyle=:dash, color=:red)
    savefig(p_v, "output/czm_voltage.png")
    println("  ✓ 保存: output/czm_voltage.png")
    
    # 温度演化
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        p_t = plot(t, T_mean, xlabel="Time (s)", ylabel="Temperature (K)",
                  label="Mean Temperature", linewidth=2, 
                  title="Temperature Evolution")
        savefig(p_t, "output/czm_temperature.png")
        println("  ✓ 保存: output/czm_temperature.png")
        
        @printf("  温升: %.2f K\n", T_mean[end] - T_mean[1])
    end
    
    # ========================================================================
    # 总结
    # ========================================================================
    println("\n" * "="^80)
    println("CZM耦合仿真完成总结")
    println("="^80)
    
    println("""
    ✓ 电化学-热-力-损伤耦合仿真完成
    
    关键结果：
      - 总时间步数: $num_steps
      - 初始电压: $(round(V[1], digits=4)) V
      - 最终电压: $(round(V[end], digits=4)) V
      - 电压降: $(round(V[1] - V[end], digits=4)) V
    """)
    
    if czm_mesh.n_cohesive > 0
        println("""
      - 内聚力单元数: $(czm_mesh.n_cohesive)
      - 最终最大损伤: $(round(damage_max_hist[end] * 100, digits=2))%
      - 最终平均损伤: $(round(damage_mean_hist[end] * 100, digits=2))%
      - 断裂单元数: $(damage_fractured_hist[end])
        """)
    end
    
    println("""
    生成的图像：
      1. output/czm_damage_evolution.png - 损伤变量随时间变化曲线
      2. output/czm_temperature_damage_field.png - 温度场和损伤场分布
      3. output/czm_voltage.png - 放电曲线
      4. output/czm_temperature.png - 温度演化
    
    物理耦合：
      ✓ 电化学 → 热源（焦耳热 + 反应热）
      ✓ 温度场 → 热应变
      ✓ SOC变化 → 化学应变
      ✓ 应变 → 层间应力 → 损伤演化
    
    下一步建议：
      - 增加仿真时间观察损伤累积
      - 实现充放电循环以模拟疲劳损伤
      - 调整内聚力参数研究参数敏感性
    """)
    
    println("="^80)
    
    return result, czm_mesh
end

# 运行主函数
result, czm_mesh = main()
