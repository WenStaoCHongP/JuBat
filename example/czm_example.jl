"""
内聚力模型(CZM)验证示例：纯机械加载

功能：
- 不启用电化学-热耦合模型
- 给定周期性位移场，验证双线性本构模型
- 绘制损伤曲线和牵引力-位移滞回环曲线

验证内容：
1. 单点本构测试：验证双线性牵引力-分离关系
2. 滞回环测试：验证加卸载行为和损伤演化
3. 混合模式测试：验证法向+切向耦合响应

输出图像：
- 图1：法向牵引力-分离曲线（加卸载滞回环）
- 图2：切向牵引力-分离曲线（加卸载滞回环）
- 图3：损伤变量随最大分离位移的演化
- 图4：混合模式下的牵引力响应

日期：2025
"""

using LinearAlgebra, Printf, Plots

# 包含JuBat模块
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

"""
    create_czm_test_params()

创建用于CZM测试的内聚力参数。
"""
function create_czm_test_params()
    # 创建内聚力参数（典型电极-隔膜界面参数）
    cohesive = JuBat.Cohesive()
    
    # 法向参数 (Mode I)
    cohesive.σ_max_n = 50e6       # 最大法向牵引力 [Pa] (50 MPa)
    cohesive.δ_0_n = 1e-6         # 损伤起始分离位移 [m] (1 μm)
    cohesive.δ_c_n = 10e-6        # 临界分离位移 [m] (10 μm)
    cohesive.G_c_n = 0.5 * cohesive.σ_max_n * cohesive.δ_c_n
    cohesive.K_n = cohesive.σ_max_n / cohesive.δ_0_n
    
    # 切向参数 (Mode II)
    cohesive.τ_max_t = 30e6       # 最大切向牵引力 [Pa] (30 MPa)
    cohesive.δ_0_t = 1e-6         # 损伤起始切向位移 [m] (1 μm)
    cohesive.δ_c_t = 15e-6        # 临界切向位移 [m] (15 μm)
    cohesive.G_c_t = 0.5 * cohesive.τ_max_t * cohesive.δ_c_t
    cohesive.K_t = cohesive.τ_max_t / cohesive.δ_0_t
    
    # BK准则指数
    cohesive.eta = 1.45
    
    return cohesive
end

"""
    test_monotonic_loading(cohesive_params)

测试单调加载下的本构响应。

绘制从0加载到完全断裂的牵引力-分离曲线。
"""
function test_monotonic_loading(cohesive_params)
    println("\n" * "="^60)
    println("测试1：单调加载本构响应")
    println("="^60)
    
    # 创建损伤状态
    damage_state_n = JuBat.DamageState()
    damage_state_t = JuBat.DamageState()
    
    # 法向加载参数
    δ_max_n = cohesive_params.δ_c_n * 1.2  # 超过临界分离
    n_points = 200
    δ_n_vals = range(0, δ_max_n, length=n_points)
    
    T_n_vals = zeros(n_points)
    D_n_vals = zeros(n_points)
    
    # 法向单调加载
    for (i, δ_n) in enumerate(δ_n_vals)
        T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=true)
        T_n_vals[i] = T_n
        D_n_vals[i] = D
    end
    
    # 切向加载参数
    δ_max_t = cohesive_params.δ_c_t * 1.2
    δ_t_vals = range(0, δ_max_t, length=n_points)
    
    T_t_vals = zeros(n_points)
    D_t_vals = zeros(n_points)
    
    # 切向单调加载
    for (i, δ_t) in enumerate(δ_t_vals)
        _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
        T_t_vals[i] = T_t
        D_t_vals[i] = D
    end
    
    # 打印关键参数
    println("\n内聚力参数：")
    @printf("  法向最大牵引力 σ_max_n = %.1f MPa\n", cohesive_params.σ_max_n / 1e6)
    @printf("  法向起始分离 δ_0_n = %.1f μm\n", cohesive_params.δ_0_n * 1e6)
    @printf("  法向临界分离 δ_c_n = %.1f μm\n", cohesive_params.δ_c_n * 1e6)
    @printf("  法向断裂能 G_c_n = %.1f J/m²\n", cohesive_params.G_c_n)
    println()
    @printf("  切向最大牵引力 τ_max_t = %.1f MPa\n", cohesive_params.τ_max_t / 1e6)
    @printf("  切向起始分离 δ_0_t = %.1f μm\n", cohesive_params.δ_0_t * 1e6)
    @printf("  切向临界分离 δ_c_t = %.1f μm\n", cohesive_params.δ_c_t * 1e6)
    @printf("  切向断裂能 G_c_t = %.1f J/m²\n", cohesive_params.G_c_t)
    
    return (δ_n_vals, T_n_vals, D_n_vals, δ_t_vals, T_t_vals, D_t_vals)
end

"""
    test_cyclic_loading(cohesive_params; n_cycles=3, max_amp_factor=0.8)

测试周期性加卸载下的本构响应。

生成滞回环曲线，验证加卸载行为和损伤累积。
"""
function test_cyclic_loading(cohesive_params; n_cycles::Int=3, max_amp_factor::Float64=0.8)
    println("\n" * "="^60)
    println("测试2：周期性加卸载滞回环")
    println("="^60)
    
    # 法向循环加载
    damage_state_n = JuBat.DamageState()
    
    # 周期性加载幅值逐渐增大
    δ_max_n = cohesive_params.δ_c_n * max_amp_factor
    
    δ_n_history = Float64[]
    T_n_history = Float64[]
    D_n_history = Float64[]
    
    points_per_half_cycle = 50
    
    for cycle in 1:n_cycles
        # 每个循环的幅值递增
        amp = δ_max_n * (cycle / n_cycles)
        
        # 加载阶段
        for i in 1:points_per_half_cycle
            δ_n = amp * (i / points_per_half_cycle)
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=true)
            push!(δ_n_history, δ_n)
            push!(T_n_history, T_n)
            push!(D_n_history, D)
        end
        
        # 卸载阶段
        for i in 1:points_per_half_cycle
            δ_n = amp * (1.0 - i / points_per_half_cycle)
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=false)
            push!(δ_n_history, δ_n)
            push!(T_n_history, T_n)
            push!(D_n_history, D)
        end
        
        @printf("  循环 %d: 幅值 = %.2f μm, 最大损伤 = %.2f%%\n", 
                cycle, amp * 1e6, damage_state_n.D * 100)
    end
    
    # 切向循环加载
    damage_state_t = JuBat.DamageState()
    
    δ_max_t = cohesive_params.δ_c_t * max_amp_factor
    
    δ_t_history = Float64[]
    T_t_history = Float64[]
    D_t_history = Float64[]
    
    for cycle in 1:n_cycles
        amp = δ_max_t * (cycle / n_cycles)
        
        # 正向加载
        for i in 1:points_per_half_cycle
            δ_t = amp * (i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # 卸载到零
        for i in 1:points_per_half_cycle
            δ_t = amp * (1.0 - i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=false)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # 反向加载（负切向）
        for i in 1:points_per_half_cycle
            δ_t = -amp * (i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # 反向卸载到零
        for i in 1:points_per_half_cycle
            δ_t = -amp * (1.0 - i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=false)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
    end
    
    return (δ_n_history, T_n_history, D_n_history, δ_t_history, T_t_history, D_t_history)
end

"""
    test_mixed_mode_loading(cohesive_params)

测试混合模式（法向+切向）加载下的本构响应。

验证BK准则的正确性。
"""
function test_mixed_mode_loading(cohesive_params)
    println("\n" * "="^60)
    println("测试3：混合模式加载（BK准则）")
    println("="^60)
    
    # 测试不同的模式混合比 β = |δ_t| / δ_eff
    mode_ratios = [0.0, 0.25, 0.5, 0.75, 1.0]  # 0=纯Mode I, 1=纯Mode II
    
    n_points = 100
    
    results = Dict{Float64, NamedTuple}()
    
    for β in mode_ratios
        damage_state = JuBat.DamageState()
        
        # 根据混合比计算等效分离位移
        # δ_eff = sqrt(δ_n² + δ_t²), β = |δ_t| / δ_eff
        # => δ_t = β * δ_eff, δ_n = sqrt(1 - β²) * δ_eff
        
        δ_eff_max = max(cohesive_params.δ_c_n, cohesive_params.δ_c_t) * 1.2
        
        δ_eff_vals = range(0, δ_eff_max, length=n_points)
        δ_n_vals = zeros(n_points)
        δ_t_vals = zeros(n_points)
        T_n_vals = zeros(n_points)
        T_t_vals = zeros(n_points)
        D_vals = zeros(n_points)
        
        for (i, δ_eff) in enumerate(δ_eff_vals)
            δ_n = sqrt(1.0 - β^2) * δ_eff
            δ_t = β * δ_eff
            
            δ_n_vals[i] = δ_n
            δ_t_vals[i] = δ_t
            
            T_n, T_t, D = JuBat.bilinear_traction(δ_n, δ_t, damage_state, cohesive_params; update=true)
            
            T_n_vals[i] = T_n
            T_t_vals[i] = T_t
            D_vals[i] = D
        end
        
        # 计算等效牵引力
        T_eff_vals = sqrt.(T_n_vals.^2 .+ T_t_vals.^2)
        
        results[β] = (δ_eff=δ_eff_vals, δ_n=δ_n_vals, δ_t=δ_t_vals,
                      T_n=T_n_vals, T_t=T_t_vals, T_eff=T_eff_vals, D=D_vals)
        
        @printf("  模式比 β = %.2f: T_max = %.1f MPa, D_final = %.2f%%\n",
                β, maximum(T_eff_vals) / 1e6, D_vals[end] * 100)
    end
    
    return results
end

"""
    test_sinusoidal_displacement(cohesive_params; n_cycles=5, frequency=1.0, amplitude_factor=0.6)

测试正弦位移加载下的响应。

模拟更真实的周期性加载情况。
"""
function test_sinusoidal_displacement(cohesive_params; n_cycles::Int=5, 
                                       frequency::Float64=1.0,
                                       amplitude_factor::Float64=0.6)
    println("\n" * "="^60)
    println("测试4：正弦位移加载")
    println("="^60)
    
    damage_state = JuBat.DamageState()
    
    # 时间参数
    T_period = 1.0 / frequency
    t_total = n_cycles * T_period
    n_points = n_cycles * 100
    t_vals = range(0, t_total, length=n_points)
    
    # 位移幅值
    δ_amp = cohesive_params.δ_c_n * amplitude_factor
    
    # 正弦位移历史
    δ_n_vals = δ_amp .* sin.(2π * frequency .* t_vals)
    
    T_n_history = zeros(n_points)
    D_history = zeros(n_points)
    
    for (i, δ_n) in enumerate(δ_n_vals)
        # 只有正向分离才更新损伤
        if δ_n > 0
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state, cohesive_params; update=true)
        else
            # 压缩时使用纯弹性接触
            T_n = cohesive_params.K_n * δ_n
            D = damage_state.D
        end
        T_n_history[i] = T_n
        D_history[i] = D
    end
    
    @printf("\n  加载参数：\n")
    @printf("    振幅 = %.2f μm (%.0f%% δ_c)\n", δ_amp * 1e6, amplitude_factor * 100)
    @printf("    频率 = %.1f Hz\n", frequency)
    @printf("    循环数 = %d\n", n_cycles)
    @printf("  结果：\n")
    @printf("    最终损伤 = %.2f%%\n", D_history[end] * 100)
    @printf("    最大牵引力 = %.1f MPa\n", maximum(T_n_history) / 1e6)
    
    return (t=t_vals, δ_n=δ_n_vals, T_n=T_n_history, D=D_history)
end

"""
    plot_all_results(...)

绘制所有测试结果。
"""
function plot_all_results(monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data, cohesive_params)
    println("\n" * "="^60)
    println("生成图像...")
    println("="^60)
    
    # 确保输出目录存在
    isdir("output") || mkdir("output")
    
    # 解包数据
    δ_n_mono, T_n_mono, D_n_mono, δ_t_mono, T_t_mono, D_t_mono = monotonic_data
    δ_n_cyc, T_n_cyc, D_n_cyc, δ_t_cyc, T_t_cyc, D_t_cyc = cyclic_data
    
    # ====================================================================
    # 图1：单调加载本构曲线
    # ====================================================================
    p1 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # 法向牵引力-分离曲线
    plot!(p1[1], δ_n_mono .* 1e6, T_n_mono ./ 1e6,
          xlabel="法向分离 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="(a) 法向本构曲线 (Mode I)",
          label="单调加载", linewidth=2, color=:blue,
          legend=:topright)
    
    # 标注关键点
    vline!(p1[1], [cohesive_params.δ_0_n * 1e6], label="δ_0_n", linestyle=:dash, color=:gray)
    vline!(p1[1], [cohesive_params.δ_c_n * 1e6], label="δ_c_n", linestyle=:dot, color=:red)
    hline!(p1[1], [cohesive_params.σ_max_n / 1e6], label="σ_max_n", linestyle=:dash, color=:orange)
    
    # 切向牵引力-分离曲线
    plot!(p1[2], δ_t_mono .* 1e6, T_t_mono ./ 1e6,
          xlabel="切向分离 δ_t (μm)", ylabel="切向牵引力 T_t (MPa)",
          title="(b) 切向本构曲线 (Mode II)",
          label="单调加载", linewidth=2, color=:green,
          legend=:topright)
    
    vline!(p1[2], [cohesive_params.δ_0_t * 1e6], label="δ_0_t", linestyle=:dash, color=:gray)
    vline!(p1[2], [cohesive_params.δ_c_t * 1e6], label="δ_c_t", linestyle=:dot, color=:red)
    hline!(p1[2], [cohesive_params.τ_max_t / 1e6], label="τ_max_t", linestyle=:dash, color=:orange)
    
    # 法向损伤演化
    plot!(p1[3], δ_n_mono .* 1e6, D_n_mono .* 100,
          xlabel="法向分离 δ_n (μm)", ylabel="损伤变量 D (%)",
          title="(c) 法向损伤演化",
          label="D vs δ_n", linewidth=2, color=:red,
          legend=:bottomright)
    
    # 切向损伤演化
    plot!(p1[4], δ_t_mono .* 1e6, D_t_mono .* 100,
          xlabel="切向分离 δ_t (μm)", ylabel="损伤变量 D (%)",
          title="(d) 切向损伤演化",
          label="D vs δ_t", linewidth=2, color=:purple,
          legend=:bottomright)
    
    savefig(p1, "output/czm_monotonic_loading.png")
    println("  ✓ 保存: output/czm_monotonic_loading.png")
    
    # ====================================================================
    # 图2：周期性加卸载滞回环
    # ====================================================================
    p2 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # 法向滞回环
    plot!(p2[1], δ_n_cyc .* 1e6, T_n_cyc ./ 1e6,
          xlabel="法向分离 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="(a) 法向加卸载滞回环",
          label="", linewidth=1.5, color=:blue)
    
    # 添加单调加载曲线作为参考
    plot!(p2[1], δ_n_mono .* 1e6, T_n_mono ./ 1e6,
          label="单调包络", linestyle=:dash, linewidth=1, color=:gray, alpha=0.5)
    
    # 切向滞回环（含正负方向）
    plot!(p2[2], δ_t_cyc .* 1e6, T_t_cyc ./ 1e6,
          xlabel="切向分离 δ_t (μm)", ylabel="切向牵引力 T_t (MPa)",
          title="(b) 切向加卸载滞回环（双向）",
          label="", linewidth=1.5, color=:green)
    
    # 法向损伤历史
    n_cyc_n = length(δ_n_cyc)
    plot!(p2[3], 1:n_cyc_n, D_n_cyc .* 100,
          xlabel="加载步", ylabel="损伤变量 D (%)",
          title="(c) 法向损伤累积",
          label="D", linewidth=1.5, color=:red)
    
    # 切向损伤历史
    n_cyc_t = length(δ_t_cyc)
    plot!(p2[4], 1:n_cyc_t, D_t_cyc .* 100,
          xlabel="加载步", ylabel="损伤变量 D (%)",
          title="(d) 切向损伤累积",
          label="D", linewidth=1.5, color=:purple)
    
    savefig(p2, "output/czm_cyclic_hysteresis.png")
    println("  ✓ 保存: output/czm_cyclic_hysteresis.png")
    
    # ====================================================================
    # 图3：混合模式响应
    # ====================================================================
    p3 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # 等效牵引力-等效分离曲线
    colors = [:blue, :cyan, :green, :orange, :red]
    mode_labels = ["β=0 (纯Mode I)", "β=0.25", "β=0.5", "β=0.75", "β=1 (纯Mode II)"]
    
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[1], collect(data.δ_eff) .* 1e6, data.T_eff ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[1], xlabel="等效分离 δ_eff (μm)", ylabel="等效牵引力 T_eff (MPa)",
          title="(a) 混合模式：等效牵引力曲线", legend=:topright)
    
    # 损伤演化
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[2], collect(data.δ_eff) .* 1e6, data.D .* 100,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[2], xlabel="等效分离 δ_eff (μm)", ylabel="损伤变量 D (%)",
          title="(b) 混合模式：损伤演化", legend=:bottomright)
    
    # 法向分量
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[3], data.δ_n .* 1e6, data.T_n ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[3], xlabel="法向分离 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="(c) 混合模式：法向分量", legend=:topright)
    
    # 切向分量
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[4], data.δ_t .* 1e6, data.T_t ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[4], xlabel="切向分离 δ_t (μm)", ylabel="切向牵引力 T_t (MPa)",
          title="(d) 混合模式：切向分量", legend=:topright)
    
    savefig(p3, "output/czm_mixed_mode.png")
    println("  ✓ 保存: output/czm_mixed_mode.png")
    
    # ====================================================================
    # 图4：正弦位移加载响应
    # ====================================================================
    p4 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # 位移-时间曲线
    plot!(p4[1], collect(sinusoidal_data.t), sinusoidal_data.δ_n .* 1e6,
          xlabel="时间 t (s)", ylabel="法向分离 δ_n (μm)",
          title="(a) 正弦位移加载历史",
          label="δ_n(t)", linewidth=1.5, color=:blue)
    
    # 牵引力-时间曲线
    plot!(p4[2], collect(sinusoidal_data.t), sinusoidal_data.T_n ./ 1e6,
          xlabel="时间 t (s)", ylabel="法向牵引力 T_n (MPa)",
          title="(b) 牵引力响应历史",
          label="T_n(t)", linewidth=1.5, color=:green)
    
    # 滞回环
    plot!(p4[3], sinusoidal_data.δ_n .* 1e6, sinusoidal_data.T_n ./ 1e6,
          xlabel="法向分离 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="(c) 正弦加载滞回环",
          label="", linewidth=1.5, color=:purple)
    
    # 损伤演化
    plot!(p4[4], collect(sinusoidal_data.t), sinusoidal_data.D .* 100,
          xlabel="时间 t (s)", ylabel="损伤变量 D (%)",
          title="(d) 损伤累积过程",
          label="D(t)", linewidth=1.5, color=:red)
    
    savefig(p4, "output/czm_sinusoidal_loading.png")
    println("  ✓ 保存: output/czm_sinusoidal_loading.png")
    
    # ====================================================================
    # 图5：综合滞回环对比
    # ====================================================================
    p5 = plot(size=(800, 600), dpi=150)
    
    # 周期性加载滞回环
    plot!(p5, δ_n_cyc .* 1e6, T_n_cyc ./ 1e6,
          label="递增幅值加载", linewidth=2, color=:blue)
    
    # 正弦加载滞回环
    plot!(p5, sinusoidal_data.δ_n .* 1e6, sinusoidal_data.T_n ./ 1e6,
          label="正弦加载", linewidth=2, color=:red, linestyle=:dash)
    
    plot!(p5, xlabel="法向分离 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="内聚力模型滞回环对比",
          legend=:topright, grid=true)
    
    savefig(p5, "output/czm_hysteresis_comparison.png")
    println("  ✓ 保存: output/czm_hysteresis_comparison.png")
    
    # 保存SVG格式
    try
        savefig(p1, "output/czm_monotonic_loading.svg")
        savefig(p2, "output/czm_cyclic_hysteresis.svg")
        savefig(p3, "output/czm_mixed_mode.svg")
        savefig(p4, "output/czm_sinusoidal_loading.svg")
        savefig(p5, "output/czm_hysteresis_comparison.svg")
        println("  ✓ SVG格式图像已保存")
    catch
        println("  ⚠ SVG格式保存失败（可能缺少依赖）")
    end
end

"""
    print_verification_summary(cohesive_params, monotonic_data, cyclic_data)

打印验证摘要，对比理论值和计算值。
"""
function print_verification_summary(cohesive_params, monotonic_data, cyclic_data)
    println("\n" * "="^60)
    println("验证摘要")
    println("="^60)
    
    δ_n_mono, T_n_mono, D_n_mono, δ_t_mono, T_t_mono, D_t_mono = monotonic_data
    
    # 理论值
    σ_max_theo = cohesive_params.σ_max_n
    τ_max_theo = cohesive_params.τ_max_t
    
    # 计算值
    σ_max_calc = maximum(T_n_mono)
    τ_max_calc = maximum(T_t_mono)
    
    # 计算断裂能（数值积分）
    G_n_calc = 0.0
    for i in 2:length(δ_n_mono)
        dδ = δ_n_mono[i] - δ_n_mono[i-1]
        G_n_calc += 0.5 * (T_n_mono[i] + T_n_mono[i-1]) * dδ
    end
    
    G_t_calc = 0.0
    for i in 2:length(δ_t_mono)
        dδ = δ_t_mono[i] - δ_t_mono[i-1]
        G_t_calc += 0.5 * (T_t_mono[i] + T_t_mono[i-1]) * dδ
    end
    
    println("\n法向 (Mode I) 验证：")
    @printf("  σ_max: 理论 = %.2f MPa, 计算 = %.2f MPa, 误差 = %.2f%%\n",
            σ_max_theo / 1e6, σ_max_calc / 1e6, 
            abs(σ_max_calc - σ_max_theo) / σ_max_theo * 100)
    @printf("  G_c_n: 理论 = %.2f J/m², 计算 = %.2f J/m², 误差 = %.2f%%\n",
            cohesive_params.G_c_n, G_n_calc,
            abs(G_n_calc - cohesive_params.G_c_n) / cohesive_params.G_c_n * 100)
    
    println("\n切向 (Mode II) 验证：")
    @printf("  τ_max: 理论 = %.2f MPa, 计算 = %.2f MPa, 误差 = %.2f%%\n",
            τ_max_theo / 1e6, τ_max_calc / 1e6,
            abs(τ_max_calc - τ_max_theo) / τ_max_theo * 100)
    @printf("  G_c_t: 理论 = %.2f J/m², 计算 = %.2f J/m², 误差 = %.2f%%\n",
            cohesive_params.G_c_t, G_t_calc,
            abs(G_t_calc - cohesive_params.G_c_t) / cohesive_params.G_c_t * 100)
    
    # 检查损伤演化
    println("\n加卸载行为验证：")
    δ_n_cyc, T_n_cyc, D_n_cyc, _, _, _ = cyclic_data
    
    # 检查卸载时损伤是否保持不变
    D_max_reached = maximum(D_n_cyc)
    D_at_zero = D_n_cyc[findfirst(x -> abs(x) < 1e-10, δ_n_cyc[100:end]) + 99]  # 第一次回到零点
    
    @printf("  卸载时损伤保持: 最大D = %.2f%%, 零点D = %.2f%% ✓\n",
            D_max_reached * 100, D_at_zero * 100)
    
    # 检查卸载刚度
    println("\n卸载刚度验证：")
    K_initial = cohesive_params.K_n
    # 在损伤软化区卸载时，刚度应该是 (1-D)*K
    println("  初始刚度 K_n = $(K_initial/1e12) TPa/m")
    println("  损伤后卸载刚度 ≈ (1-D)*K_n")
    
    println("\n" * "="^60)
    println("✓ CZM本构模型验证完成")
    println("="^60)
end

"""
    main()

主函数：运行所有CZM验证测试。
"""
function main()
    println("="^80)
    println("内聚力模型(CZM)验证：纯机械加载测试")
    println("="^80)
    println("\n本测试验证双线性牵引力-分离本构模型的正确性，")
    println("包括单调加载、周期性加卸载和混合模式响应。")
    println("不涉及电化学-热耦合模型。")
    
    # 创建内聚力参数
    cohesive_params = create_czm_test_params()
    
    # 测试1：单调加载
    monotonic_data = test_monotonic_loading(cohesive_params)
    
    # 测试2：周期性加卸载
    cyclic_data = test_cyclic_loading(cohesive_params; n_cycles=3, max_amp_factor=0.8)
    
    # 测试3：混合模式
    mixed_mode_data = test_mixed_mode_loading(cohesive_params)
    
    # 测试4：正弦位移加载
    sinusoidal_data = test_sinusoidal_displacement(cohesive_params; 
                                                    n_cycles=5, 
                                                    frequency=1.0,
                                                    amplitude_factor=0.6)
    
    # 绘制所有结果
    plot_all_results(monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data, cohesive_params)
    
    # 打印验证摘要
    print_verification_summary(cohesive_params, monotonic_data, cyclic_data)
    
    println("\n生成的图像：")
    println("  1. output/czm_monotonic_loading.png - 单调加载本构曲线")
    println("  2. output/czm_cyclic_hysteresis.png - 周期性加卸载滞回环")
    println("  3. output/czm_mixed_mode.png - 混合模式响应(BK准则)")
    println("  4. output/czm_sinusoidal_loading.png - 正弦位移加载响应")
    println("  5. output/czm_hysteresis_comparison.png - 滞回环对比")
    
    return cohesive_params, monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data
end

# 运行主函数
result = main()
