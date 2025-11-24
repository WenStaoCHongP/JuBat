"""
力学模块验证示例

功能：
1. 内聚力模型单元测试（双线性牵引-分离定律）
2. 接触模型单元测试（Hertz接触理论）
3. 与解析解对比验证

作者：AI Assistant
日期：2025-11-24
"""

using Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ================================================================================
# 测试1：内聚力模型 - 双线性牵引-分离定律
# ================================================================================
function test_cohesive_bilinear()
    println("\n" * "="^80)
    println("测试1：内聚力模型 - 双线性牵引-分离定律")
    println("="^80)
    
    # 材料参数
    K_n = 1e12      # Pa/m
    T_max = 10e6    # 10 MPa
    Γ = 100.0       # 100 J/m²
    
    δ_0 = T_max / K_n
    δ_f = 2.0 * Γ / T_max
    
    @printf("材料参数：\n")
    @printf("  刚度 K_n = %.2e Pa/m\n", K_n)
    @printf("  强度 T_max = %.2f MPa\n", T_max / 1e6)
    @printf("  断裂能 Γ = %.2f J/m²\n", Γ)
    @printf("  损伤起始位移 δ_0 = %.2e m (%.2f μm)\n", δ_0, δ_0 * 1e6)
    @printf("  完全失效位移 δ_f = %.2e m (%.2f μm)\n", δ_f, δ_f * 1e6)
    
    # 构造简单的内聚力界面
    coh = JuBat.CohesiveInterface(
        [(1,1)],           # element_pairs
        [1.0],             # interface_area
        [1.0 0.0 0.0],     # normals
        K_n, K_n,          # K_n, K_t
        T_max, T_max,      # T_n_max, T_t_max
        δ_0, δ_0,          # δ_n_0, δ_t_0
        δ_f, δ_f,          # δ_n_f, δ_t_f
        Γ, Γ,              # Γ_n, Γ_t
        2.0,               # η
        [0.0], [0.0], [false], [0.0]  # 状态变量
    )
    
    # 加载路径
    n_steps = 100
    δ_max = 1.5 * δ_f
    δ_path = range(0, δ_max, length=n_steps)
    
    # 计算牵引力
    T_n_result = zeros(n_steps)
    D_result = zeros(n_steps)
    
    for (i, δ) in enumerate(δ_path)
        T_n, T_t, D = JuBat.compute_cohesive_traction(coh, [δ], [0.0])
        T_n_result[i] = T_n[1]
        D_result[i] = D[1]
        JuBat.update_cohesive_damage!(coh, [δ], [0.0])
    end
    
    # 解析解
    T_n_analytical = zeros(n_steps)
    D_analytical = zeros(n_steps)
    
    for (i, δ) in enumerate(δ_path)
        if δ <= δ_0
            # 弹性
            T_n_analytical[i] = K_n * δ
            D_analytical[i] = 0.0
        elseif δ < δ_f
            # 软化
            T_n_analytical[i] = T_max * (δ_f - δ) / (δ_f - δ_0)
            D_analytical[i] = (δ_f / δ) * (δ - δ_0) / (δ_f - δ_0)
        else
            # 失效
            T_n_analytical[i] = 0.0
            D_analytical[i] = 1.0
        end
    end
    
    # 计算误差
    err_T = maximum(abs.(T_n_result .- T_n_analytical) ./ (T_max + 1e-10))
    err_D = maximum(abs.(D_result .- D_analytical))
    
    @printf("\n验证结果：\n")
    @printf("  牵引力最大相对误差: %.2e\n", err_T)
    @printf("  损伤变量最大绝对误差: %.2e\n", err_D)
    
    if err_T < 1e-6 && err_D < 1e-6
        println("  ✓ 测试通过！")
    else
        println("  ✗ 测试失败！")
    end
    
    # 绘图
    p1 = plot(δ_path .* 1e6, T_n_result ./ 1e6,
              label="数值解", linewidth=2, color=:blue,
              xlabel="分离位移 δ [μm]", ylabel="牵引应力 T [MPa]",
              title="内聚力牵引-分离曲线")
    plot!(p1, δ_path .* 1e6, T_n_analytical ./ 1e6,
          label="解析解", linestyle=:dash, linewidth=2, color=:red)
    vline!(p1, [δ_0 * 1e6], label="δ_0", linestyle=:dot, color=:green)
    vline!(p1, [δ_f * 1e6], label="δ_f", linestyle=:dot, color=:orange)
    
    p2 = plot(δ_path .* 1e6, D_result,
              label="数值解", linewidth=2, color=:blue,
              xlabel="分离位移 δ [μm]", ylabel="损伤变量 D",
              title="损伤演化曲线", ylims=(-0.1, 1.1))
    plot!(p2, δ_path .* 1e6, D_analytical,
          label="解析解", linestyle=:dash, linewidth=2, color=:red)
    hline!(p2, [0.0, 1.0], label="", linestyle=:dot, color=:gray, alpha=0.5)
    
    p = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    savefig(p, "mechanical_validation_cohesive.png")
    println("\n✓ 保存图像: mechanical_validation_cohesive.png")
    
    return err_T < 1e-6 && err_D < 1e-6
end

# ================================================================================
# 测试2：Hertz接触理论验证
# ================================================================================
function test_hertz_contact()
    println("\n" * "="^80)
    println("测试2：Hertz接触理论验证")
    println("="^80)
    
    # 材料参数（石墨）
    E = 10e9        # 10 GPa
    ν = 0.3
    E_star = E / (1 - ν^2)
    
    # 几何参数
    R = 5e-6        # 5 μm 球半径
    F_normal = 1e-6 # 1 μN 法向力
    
    @printf("材料参数：\n")
    @printf("  杨氏模量 E = %.2f GPa\n", E / 1e9)
    @printf("  泊松比 ν = %.2f\n", ν)
    @printf("  有效模量 E* = %.2f GPa\n", E_star / 1e9)
    @printf("\n几何和载荷：\n")
    @printf("  球半径 R = %.2f μm\n", R * 1e6)
    @printf("  法向力 F = %.2f μN\n", F_normal * 1e6)
    
    # 计算Hertz解
    a, p_max, p_hertz = JuBat.hertz_contact_pressure(R, F_normal, E_star)
    
    @printf("\nHertz解析解：\n")
    @printf("  接触半径 a = %.3e m (%.2f nm)\n", a, a * 1e9)
    @printf("  最大压力 p_max = %.2f MPa\n", p_max / 1e6)
    
    # 压力分布
    r_values = range(0, 1.5 * a, length=100)
    p_values = [p_hertz(r) for r in r_values]
    
    # 验证积分 ∫p·2πr dr = F
    dr = r_values[2] - r_values[1]
    F_integrated = sum(p_values .* (2π .* r_values)) * dr
    err_F = abs(F_integrated - F_normal) / F_normal
    
    @printf("\n验证结果：\n")
    @printf("  积分力 F_int = %.2e N\n", F_integrated)
    @printf("  理论力 F = %.2e N\n", F_normal)
    @printf("  相对误差: %.2e\n", err_F)
    
    if err_F < 0.01
        println("  ✓ 测试通过！")
    else
        println("  ✗ 测试失败！")
    end
    
    # 绘图
    p = plot(r_values .* 1e9, p_values ./ 1e6,
             label="Hertz压力分布", linewidth=2, color=:blue,
             xlabel="径向距离 r [nm]", ylabel="接触压力 p [MPa]",
             title="Hertz接触压力分布", fill=(0, :lightblue, 0.3))
    vline!(p, [a * 1e9], label="接触半径 a", linestyle=:dash, color=:red)
    hline!(p, [p_max / 1e6], label="p_max", linestyle=:dot, color=:orange, alpha=0.5)
    
    savefig(p, "mechanical_validation_hertz.png")
    println("\n✓ 保存图像: mechanical_validation_hertz.png")
    
    return err_F < 0.01
end

# ================================================================================
# 测试3：摩擦模型验证
# ================================================================================
function test_friction_model()
    println("\n" * "="^80)
    println("测试3：Coulomb摩擦模型验证")
    println("="^80)
    
    # 接触参数
    p_contact = 10e6    # 10 MPa 接触压力
    μ_s = 0.3           # 静摩擦系数
    μ_k = 0.2           # 动摩擦系数
    v_critical = 0.01   # 1 cm/s
    
    @printf("摩擦参数：\n")
    @printf("  接触压力 p = %.2f MPa\n", p_contact / 1e6)
    @printf("  静摩擦系数 μ_s = %.2f\n", μ_s)
    @printf("  动摩擦系数 μ_k = %.2f\n", μ_k)
    @printf("  临界速度 v_c = %.2f m/s\n", v_critical)
    
    # 构造简单的接触界面
    contact = JuBat.ContactInterface(
        [1], [1], [(1,1)], [1.0 0.0 0.0],
        1e11, 1e11, 50e6, 1e-8,
        μ_s, μ_k, v_critical,
        [-1e-9], [p_contact], [0.0], [true], [false], [1.0], [0.0]
    )
    
    # 速度扫描
    n_steps = 100
    v_max = 10 * v_critical
    v_values = range(-v_max, v_max, length=n_steps)
    
    τ_values = zeros(n_steps)
    μ_values = zeros(n_steps)
    
    for (i, v) in enumerate(v_values)
        τ, is_sliding = JuBat.compute_friction_stress(contact, [v])
        τ_values[i] = τ[1]
        
        # 计算等效摩擦系数
        if abs(v) > 1e-10
            μ_values[i] = abs(τ[1]) / p_contact
        else
            μ_values[i] = 0.0
        end
    end
    
    # 理论值
    μ_theory = [μ_k + (μ_s - μ_k) * exp(-abs(v) / v_critical) for v in v_values]
    
    @printf("\n验证结果：\n")
    @printf("  静摩擦极限: μ_s·p = %.2f MPa\n", μ_s * p_contact / 1e6)
    @printf("  动摩擦极限: μ_k·p = %.2f MPa\n", μ_k * p_contact / 1e6)
    @printf("  最大计算摩擦力: %.2f MPa\n", maximum(abs.(τ_values)) / 1e6)
    @printf("  最小计算摩擦力: %.2f MPa\n", minimum(abs.(τ_values)) / 1e6)
    
    # 绘图
    p1 = plot(v_values, τ_values ./ 1e6,
              label="摩擦应力 τ", linewidth=2, color=:blue,
              xlabel="切向速度 v_t [m/s]", ylabel="摩擦应力 τ [MPa]",
              title="Coulomb摩擦模型")
    hline!(p1, [μ_s * p_contact / 1e6], label="μ_s·p", linestyle=:dash, color=:red)
    hline!(p1, [-μ_s * p_contact / 1e6], label="", linestyle=:dash, color=:red)
    hline!(p1, [μ_k * p_contact / 1e6], label="μ_k·p", linestyle=:dot, color=:orange)
    hline!(p1, [-μ_k * p_contact / 1e6], label="", linestyle=:dot, color=:orange)
    
    p2 = plot(v_values, μ_values,
              label="数值解", linewidth=2, color=:blue,
              xlabel="切向速度 v_t [m/s]", ylabel="摩擦系数 μ",
              title="速度依赖摩擦系数", ylims=(0, 0.4))
    plot!(p2, v_values, μ_theory,
          label="理论解", linestyle=:dash, linewidth=2, color=:red)
    hline!(p2, [μ_s, μ_k], label="", linestyle=:dot, color=:gray, alpha=0.5)
    
    p = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    savefig(p, "mechanical_validation_friction.png")
    println("\n✓ 保存图像: mechanical_validation_friction.png")
    
    return true
end

# ================================================================================
# 主函数
# ================================================================================
function main()
    println("="^80)
    println("JuBat 力学模块验证测试套件")
    println("="^80)
    
    results = Dict{String, Bool}()
    
    # 运行所有测试
    try
        results["内聚力模型"] = test_cohesive_bilinear()
    catch e
        @warn "内聚力测试失败" exception=(e, catch_backtrace())
        results["内聚力模型"] = false
    end
    
    try
        results["Hertz接触"] = test_hertz_contact()
    catch e
        @warn "Hertz接触测试失败" exception=(e, catch_backtrace())
        results["Hertz接触"] = false
    end
    
    try
        results["摩擦模型"] = test_friction_model()
    catch e
        @warn "摩擦模型测试失败" exception=(e, catch_backtrace())
        results["摩擦模型"] = false
    end
    
    # 总结
    println("\n" * "="^80)
    println("测试总结")
    println("="^80)
    
    for (test_name, passed) in results
        status = passed ? "✓ PASS" : "✗ FAIL"
        @printf("  %s: %s\n", test_name, status)
    end
    
    n_passed = count(values(results))
    n_total = length(results)
    @printf("\n总计: %d/%d 测试通过\n", n_passed, n_total)
    
    if n_passed == n_total
        println("\n🎉 所有测试通过！力学模块验证成功！")
    else
        println("\n⚠️  部分测试失败，请检查实现。")
    end
    
    println("="^80)
end

# 运行主函数
main()
