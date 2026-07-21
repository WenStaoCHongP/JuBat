"""
验证CZM参数设置 - 简化测试

目的：验证调整后的内聚力参数能够实现：
1. 损伤能够持续演化（发散）
2. 循环次数 N > 1000 时达到断裂

损伤演化机制（双线性本构，无疲劳模型）：
- D = δ_c × (δ_eff - δ_0) / (δ_eff × (δ_c - δ_0))
- 损伤只随分离位移增加而增加
- 分离位移通过热-化学应变的循环累积增长

日期：2025
"""

using LinearAlgebra, Printf

# 包含JuBat模块
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ========================================================================
# 1. 参数加载和显示
# ========================================================================

function load_and_display_parameters()
    println("="^70)
    println("CZM参数验证")
    println("="^70)
    
    # 加载Jellyroll参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive
    
    println("\n[当前内聚力参数]")
    println("-"^50)
    @printf("法向 (Mode I):\n")
    # @printf("  σ_max_n = %.2f MPa\n", coh.σ_max_n / 1e6)      # TODO Chunk 2 Task 2.1
    # @printf("  δ_0_n = %.1f nm\n", coh.δ_0_n * 1e9)            # TODO Chunk 2 Task 2.1
    # @printf("  δ_c_n = %.1f nm\n", coh.δ_c_n * 1e9)            # TODO Chunk 2 Task 2.1
    # @printf("  K_n = %.2e Pa/m\n", coh.K_n)                    # TODO Chunk 2 Task 2.1
    # @printf("  G_c_n = %.4f J/m²\n", coh.G_c_n)                # TODO Chunk 2 Task 2.1

    @printf("\n切向 (Mode II):\n")
    # @printf("  τ_max_t = %.2f MPa\n", coh.τ_max_t / 1e6)       # TODO Chunk 2 Task 2.1
    # @printf("  δ_0_t = %.1f nm\n", coh.δ_0_t * 1e9)            # TODO Chunk 2 Task 2.1
    # @printf("  δ_c_t = %.1f nm\n", coh.δ_c_t * 1e9)            # TODO Chunk 2 Task 2.1
    # @printf("  K_t = %.2e Pa/m\n", coh.K_t)                    # TODO Chunk 2 Task 2.1
    # @printf("  G_c_t = %.4f J/m²\n", coh.G_c_t)                # TODO Chunk 2 Task 2.1
    
    @printf("\n混合模式:\n")
    @printf("  eta (BK指数) = %.2f\n", coh.eta)
    
    return param_dim
end

# ========================================================================
# 2. 双线性本构模型分析
# ========================================================================

"""
计算给定分离位移下的损伤值（双线性本构）
"""
function compute_damage(δ_eff, δ_0, δ_c)
    if δ_eff <= δ_0
        return 0.0
    elseif δ_eff >= δ_c
        return 1.0
    else
        return δ_c * (δ_eff - δ_0) / (δ_eff * (δ_c - δ_0))
    end
end

"""
模拟损伤演化（纯双线性本构，无疲劳模型）

假设分离位移随循环次数线性累积：δ_N = δ_init + N × Δδ
"""
function simulate_damage_evolution(n_cycles, δ_init, Δδ_per_cycle, δ_0, δ_c)
    δ_history = zeros(n_cycles + 1)
    D_history = zeros(n_cycles + 1)
    
    δ_history[1] = δ_init
    D_history[1] = compute_damage(δ_init, δ_0, δ_c)
    
    for i in 1:n_cycles
        δ_history[i+1] = δ_init + i * Δδ_per_cycle
        D_history[i+1] = compute_damage(δ_history[i+1], δ_0, δ_c)
    end
    
    return D_history, δ_history
end

"""
分析参数：计算损伤演化和断裂循环数
"""
function analyze_parameters()
    println("\n" * "="^70)
    println("损伤演化分析")
    println("="^70)
    
    # 加载参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive
    
    # δ_0 = coh.δ_0_n  # TODO Chunk 2 Task 2.1
    # δ_c = coh.δ_c_n  # TODO Chunk 2 Task 2.1
    
    # 假设的初始分离位移和每循环增量
    # 基于原始问题：损伤收敛于 6.83%，反推初始分离位移
    # 使用原参数：δ_0 = 15 nm, δ_c = 100 nm, D = 6.83%
    # 0.0683 = 100 × (δ - 15) / (δ × 85) => δ ≈ 15.93 nm
    δ_init = 16e-9  # 16 nm
    
    # 每循环分离位移增量（需要根据实际模型调整）
    Δδ_scenarios = [0.005e-9, 0.01e-9, 0.02e-9]  # 0.005, 0.01, 0.02 nm/cycle
    
    println("\n[参数设置]")
    @printf("  δ_0 = %.1f nm\n", δ_0 * 1e9)
    @printf("  δ_c = %.1f nm\n", δ_c * 1e9)
    @printf("  δ_init = %.1f nm (估算)\n", δ_init * 1e9)
    
    println("\n[初始损伤]")
    D_init = compute_damage(δ_init, δ_0, δ_c)
    @printf("  D(δ=%.1f nm) = %.1f%%\n", δ_init * 1e9, D_init * 100)
    
    println("\n[损伤演化预测]")
    println("-"^60)
    @printf("%-15s %-12s %-12s %-12s %-12s\n", 
            "Δδ(nm/cyc)", "D@100cyc", "D@500cyc", "D@1000cyc", "断裂循环")
    println("-"^60)
    
    results = Dict()
    
    for Δδ in Δδ_scenarios
        n_cycles = 3000
        D_hist, δ_hist = simulate_damage_evolution(n_cycles, δ_init, Δδ, δ_0, δ_c)
        
        # 找断裂循环
        fracture_cycle = findfirst(D_hist .>= 0.99)
        fracture_str = fracture_cycle === nothing ? ">$(n_cycles)" : string(fracture_cycle - 1)
        
        @printf("%-15.3f %-12.1f %-12.1f %-12.1f %-12s\n",
                Δδ * 1e9, 
                D_hist[101] * 100, 
                D_hist[501] * 100, 
                D_hist[1001] * 100,
                fracture_str)
        
        results[Δδ] = (D=D_hist, δ=δ_hist, fracture=fracture_cycle)
    end
    
    return results, param_dim
end

# ========================================================================
# 3. 与原参数对比
# ========================================================================

function compare_with_original()
    println("\n" * "="^70)
    println("新旧参数对比")
    println("="^70)
    
    # 原参数
    δ_0_old = 15e-9
    δ_c_old = 100e-9
    
    # 新参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    # δ_0_new = param_dim.cohesive.δ_0_n  # TODO Chunk 2 Task 2.1
    # δ_c_new = param_dim.cohesive.δ_c_n  # TODO Chunk 2 Task 2.1
    
    # 典型分离位移范围
    δ_test = [15e-9, 16e-9, 18e-9, 20e-9, 25e-9, 30e-9]
    
    println("\n[不同分离位移下的损伤值]")
    println("-"^60)
    @printf("%-12s %-15s %-15s\n", "δ (nm)", "D_old (%)", "D_new (%)")
    println("-"^60)
    
    for δ in δ_test
        D_old = compute_damage(δ, δ_0_old, δ_c_old)
        D_new = compute_damage(δ, δ_0_new, δ_c_new)
        @printf("%-12.1f %-15.1f %-15.1f\n", δ * 1e9, D_old * 100, D_new * 100)
    end
    
    println("\n[关键差异]")
    println("  原参数: δ_0=15nm, δ_c=100nm → 损伤演化缓慢，容易收敛")
    println("  新参数: δ_0=$(round(δ_0_new*1e9, digits=1))nm, δ_c=$(round(δ_c_new*1e9, digits=1))nm → 损伤演化加快，能够发散")
end

# ========================================================================
# 4. 推荐和总结
# ========================================================================

function print_summary()
    println("\n" * "="^70)
    println("参数调整总结")
    println("="^70)
    
    println("""
    
    [原参数问题]
    - δ_c_n = 100 nm 远大于实际分离位移（~16 nm）
    - 损伤公式: D = δ_c × (δ - δ_0) / (δ × (δ_c - δ_0))
    - 当 δ = 16 nm, δ_0 = 15 nm, δ_c = 100 nm 时:
      D = 100 × (16-15) / (16 × 85) = 100/1360 = 7.4% ≈ 6.83%
    - 损伤收敛是因为 δ_c 太大，δ 难以达到
    
    [新参数设计]
    - 降低 δ_c 到 30 nm，使其在循环累积分离后可达
    - 设置 δ_0 = 14 nm，使初始损伤 ~25%
    - 设计目标：每循环累积 0.01 nm 时，约 1400 循环断裂
    
    [损伤演化机制]
    - 双线性本构：损伤随分离位移单调增加
    - 分离位移通过热-化学应变循环累积
    - 无额外疲劳模型：损伤完全由分离位移决定
    
    [验证方法]
    1. 运行 czm_cycle_example.jl 进行完整耦合仿真
    2. 观察损伤是否持续增长而非收敛
    3. 确认断裂循环数 N > 1000
    
    [进一步调整建议]
    - 如断裂太快：增大 δ_c 或降低 δ_0
    - 如断裂太慢：降低 δ_c 或增大热膨胀系数
    - 可通过观察实际分离位移增量来精调参数
    """)
end

# ========================================================================
# 主函数
# ========================================================================

function main()
    # 1. 显示参数
    param_dim = load_and_display_parameters()
    
    # 2. 分析损伤演化
    results, _ = analyze_parameters()
    
    # 3. 新旧对比
    compare_with_original()
    
    # 4. 总结
    print_summary()
    
    println("\n" * "="^70)
    println("验证完成！")
    println("="^70)
    
    return results
end

# 运行
if abspath(PROGRAM_FILE) == @__FILE__
    results = main()
end
