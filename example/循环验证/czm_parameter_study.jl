"""
CZM参数敏感性研究

目标：找到合适的内聚力参数使得：
1. 损伤能够持续演化（发散）而不是收敛到固定值
2. 循环次数 N > 1000 时达到完全断裂

分析方法：
1. 估算热-化学应变产生的实际分离位移
2. 参数敏感性扫描
3. 确定合理的参数范围

损伤模型：双线性本构（无疲劳模型）
D = δ_c × (δ_eff - δ_0) / (δ_eff × (δ_c - δ_0))

日期：2025
"""

using LinearAlgebra, Printf

# 包含JuBat模块
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ========================================================================
# 1. 应变和分离位移估算
# ========================================================================

"""
估算热-化学应变产生的分离位移
"""
function estimate_separation_displacement()
    println("="^70)
    println("应变和分离位移估算")
    println("="^70)
    
    # 材料参数 (from Jellyroll.jl)
    α_NE = 3.0e-6    # 负极热膨胀系数 [1/K]
    α_PE = 1.0e-5    # 正极热膨胀系数 [1/K]
    Ω_NE = 3.1e-6    # 负极部分摩尔体积 [m³/mol]
    Ω_PE = -7.28e-7  # 正极部分摩尔体积 [m³/mol]
    
    # 有效参数（厚度加权）
    t_NE = 8.52e-5
    t_PE = 7.56e-5
    α_eff = (α_NE * t_NE + α_PE * t_PE) / (t_NE + t_PE)
    β_n = Ω_NE / 3.0
    β_p = Ω_PE / 3.0
    
    println("\n[材料参数]")
    @printf("  α_eff = %.2e 1/K\n", α_eff)
    @printf("  β_n (负极) = %.2e\n", β_n)
    @printf("  β_p (正极) = %.2e\n", β_p)
    
    # 典型工况
    ΔT_typical = 10.0     # 典型温升 [K]
    ΔSOC_typical = 0.5    # 典型SOC变化
    
    # 应变计算
    ε_thermal = α_eff * ΔT_typical
    ε_chem_n = β_n * ΔSOC_typical
    ε_chem_p = β_p * ΔSOC_typical
    
    println("\n[应变估算] (ΔT=$ΔT_typical K, ΔSOC=$ΔSOC_typical)")
    @printf("  热应变 ε_th = %.2e (%.1f με)\n", ε_thermal, ε_thermal * 1e6)
    @printf("  化学应变 (负极) ε_chem_n = %.2e (%.1f με)\n", ε_chem_n, ε_chem_n * 1e6)
    @printf("  化学应变 (正极) ε_chem_p = %.2e (%.1f με)\n", ε_chem_p, ε_chem_p * 1e6)
    @printf("  总应变 ε_total ≈ %.2e (%.1f με)\n", ε_thermal + abs(ε_chem_n), (ε_thermal + abs(ε_chem_n)) * 1e6)
    
    # 分离位移估算
    L_char = 1e-3  # 特征长度 ~ 1 mm (典型单元尺寸)
    t_layer = 2e-4 # 单层厚度 ~ 200 μm
    
    # 情景1：均匀应变，由于几何约束产生的分离
    Δε_gradient = 1e-5  # 假设相邻单元间的应变差 ~ 10 με
    δ_from_gradient = Δε_gradient * L_char
    
    # 情景2：弯曲效应（内外圈曲率不同）
    R_mean = 8e-3  # 平均半径 ~ 8 mm
    ΔR = t_layer   # 层厚
    ε_bending = ε_thermal * ΔR / R_mean
    δ_from_bending = ε_bending * L_char
    
    println("\n[分离位移估算]")
    @printf("  梯度效应: δ ≈ %.1f nm (Δε=%.0f με, L=%.0f mm)\n", 
            δ_from_gradient * 1e9, Δε_gradient * 1e6, L_char * 1e3)
    @printf("  弯曲效应: δ ≈ %.1f nm\n", δ_from_bending * 1e9)
    
    # 综合估算
    δ_eff_estimate = sqrt(δ_from_gradient^2 + δ_from_bending^2)
    @printf("  综合估算: δ_eff ≈ %.1f nm\n", δ_eff_estimate * 1e9)
    
    return δ_eff_estimate
end

# ========================================================================
# 2. 双线性本构参数研究
# ========================================================================

"""
计算给定分离位移下的损伤值
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
参数敏感性研究：扫描不同的 δ_c 和 δ_0 组合
"""
function parameter_sensitivity_study()
    println("\n" * "="^70)
    println("CZM参数敏感性研究")
    println("="^70)
    
    # 估算的实际分离位移（初始）
    δ_init = 16e-9  # 16 nm (基于原损伤值 6.83% 反推)
    
    # 参数扫描范围
    δ_0_range = [5e-9, 10e-9, 12e-9, 14e-9, 15e-9]  # nm
    δ_c_range = [20e-9, 25e-9, 30e-9, 40e-9, 50e-9]  # nm
    
    println("\n[参数扫描] δ_init = $(δ_init*1e9) nm")
    println("-"^70)
    @printf("%-10s %-10s %-12s %-12s %-15s\n", 
            "δ_0(nm)", "δ_c(nm)", "D_init(%)", "δ_c/δ_0", "评估")
    println("-"^70)
    
    for δ_0 in δ_0_range
        for δ_c in δ_c_range
            if δ_c <= δ_0 || δ_c <= δ_init
                continue
            end
            
            D_init = compute_damage(δ_init, δ_0, δ_c)
            ratio = δ_c / δ_0
            
            # 评估
            # 目标：初始损伤 20-40%，分离位移累积后能达到 δ_c
            if D_init < 0.1
                evaluation = "损伤起始太晚"
            elseif D_init < 0.2
                evaluation = "偏保守"
            elseif D_init < 0.4
                evaluation = "✓ 推荐"
            elseif D_init < 0.6
                evaluation = "偏激进"
            else
                evaluation = "断裂过快"
            end
            
            @printf("%-10.0f %-10.0f %-12.1f %-12.1f %-15s\n", 
                    δ_0*1e9, δ_c*1e9, D_init*100, ratio, evaluation)
        end
    end
end

# ========================================================================
# 3. 循环断裂预测
# ========================================================================

"""
预测达到断裂所需的循环次数
"""
function predict_fracture_cycles()
    println("\n" * "="^70)
    println("断裂循环数预测")
    println("="^70)
    
    # 加载当前参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive
    
    δ_0 = coh.δ_0_n
    δ_c = coh.δ_c_n
    δ_init = 16e-9  # 初始分离位移
    
    println("\n[当前参数]")
    @printf("  δ_0 = %.1f nm\n", δ_0 * 1e9)
    @printf("  δ_c = %.1f nm\n", δ_c * 1e9)
    @printf("  δ_init = %.1f nm\n", δ_init * 1e9)
    
    # 不同分离位移增量下的断裂循环数
    Δδ_range = [0.002e-9, 0.005e-9, 0.01e-9, 0.02e-9, 0.05e-9]  # nm/cycle
    
    println("\n[断裂循环数预测]")
    println("-"^50)
    @printf("%-15s %-15s %-15s\n", "Δδ(nm/cyc)", "断裂循环N", "D_init(%)")
    println("-"^50)
    
    D_init = compute_damage(δ_init, δ_0, δ_c)
    
    for Δδ in Δδ_range
        # 计算达到 δ_c 需要的循环数
        if δ_init >= δ_c
            N_fracture = 0
        else
            N_fracture = ceil(Int, (δ_c - δ_init) / Δδ)
        end
        
        @printf("%-15.3f %-15d %-15.1f\n", Δδ * 1e9, N_fracture, D_init * 100)
    end
    
    println("\n[说明]")
    println("  - 断裂循环 N = (δ_c - δ_init) / Δδ")
    println("  - Δδ 是每循环分离位移累积增量")
    println("  - 实际 Δδ 取决于热-化学应变累积效应")
    println("  - 典型值范围：0.005 - 0.02 nm/cycle")
end

# ========================================================================
# 4. 推荐参数
# ========================================================================

"""
生成推荐的CZM参数
"""
function recommend_parameters()
    println("\n" * "="^70)
    println("推荐的CZM参数")
    println("="^70)
    
    println("""
    
    [分析结论]
    原参数问题：
    1. δ_c_n = 100 nm 太大，实际分离位移约 16 nm
    2. 损伤 D = 6.83% 收敛，因为 δ 无法达到 δ_c
    3. 比值 δ_c/δ_0 = 100/15 ≈ 6.7
    
    [推荐参数]
    法向参数：
      σ_max_n = 1.5 MPa (降低界面强度)
      δ_0_n = 14 nm    (初始分离 ~16nm 时开始损伤)
      δ_c_n = 30 nm    (约 1400 循环达到，假设 Δδ=0.01nm/cyc)
      → K_n = 1.07e14 Pa/m
      → G_c_n = 0.0225 J/m²
    
    切向参数：
      τ_max_t = 0.8 MPa
      δ_0_t = 40 nm
      δ_c_t = 90 nm
      → K_t = 2.0e13 Pa/m
      → G_c_t = 0.036 J/m²
    
    混合模式：
      eta = 1.45 (BK指数，保持不变)
    
    [预期效果]
    - 初始损伤（δ=16nm）: D ≈ 25%
    - 损伤持续演化，不会收敛
    - 预计断裂循环: 1000-2000 (取决于实际 Δδ)
    
    [调整建议]
    - 如需加快断裂: 降低 δ_c (如 25nm → ~1100循环)
    - 如需减慢断裂: 增大 δ_c (如 40nm → ~2400循环)
    """)
end

# ========================================================================
# 5. 生成新的参数代码
# ========================================================================

"""
生成可以直接复制到 Jellyroll.jl 的参数代码
"""
function generate_parameter_code()
    println("\n" * "="^70)
    println("可复制的参数代码 (已更新到 Jellyroll.jl)")
    println("="^70)
    
    # 加载当前参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive
    
    code = """
# Cohesive zone model parameters
cohesive = Cohesive()

# 法向参数 (Mode I)
cohesive.σ_max_n = $(coh.σ_max_n)       # [Pa] ($(coh.σ_max_n/1e6) MPa)
cohesive.δ_0_n = $(coh.δ_0_n)          # [m] ($(coh.δ_0_n*1e9) nm)
cohesive.δ_c_n = $(coh.δ_c_n)          # [m] ($(coh.δ_c_n*1e9) nm)
cohesive.G_c_n = 0.5 * cohesive.σ_max_n * cohesive.δ_c_n  # ≈ $(round(coh.G_c_n, digits=4)) J/m²
cohesive.K_n = cohesive.σ_max_n / cohesive.δ_0_n          # ≈ $(round(coh.K_n, sigdigits=3)) Pa/m

# 切向参数 (Mode II)
cohesive.τ_max_t = $(coh.τ_max_t)       # [Pa] ($(coh.τ_max_t/1e6) MPa)
cohesive.δ_0_t = $(coh.δ_0_t)          # [m] ($(coh.δ_0_t*1e9) nm)
cohesive.δ_c_t = $(coh.δ_c_t)          # [m] ($(coh.δ_c_t*1e9) nm)
cohesive.G_c_t = 0.5 * cohesive.τ_max_t * cohesive.δ_c_t  # ≈ $(round(coh.G_c_t, digits=4)) J/m²
cohesive.K_t = cohesive.τ_max_t / cohesive.δ_0_t          # ≈ $(round(coh.K_t, sigdigits=3)) Pa/m

# 混合模式参数
cohesive.eta = $(coh.eta)
"""
    
    println(code)
    return code
end

# ========================================================================
# 主函数
# ========================================================================

function main()
    println("="^70)
    println("CZM参数敏感性研究 - 21700果冻卷电池界面脱粘模型")
    println("="^70)
    
    # 1. 应变和分离位移估算
    δ_estimate = estimate_separation_displacement()
    
    # 2. 参数敏感性扫描
    parameter_sensitivity_study()
    
    # 3. 断裂循环预测
    predict_fracture_cycles()
    
    # 4. 推荐参数
    recommend_parameters()
    
    # 5. 生成代码
    generate_parameter_code()
    
    println("\n" * "="^70)
    println("研究完成！")
    println("="^70)
    println("\n下一步：")
    println("1. 运行 czm_cycle_example.jl 验证损伤演化")
    println("2. 根据实际结果微调 δ_c 参数")
end

# 运行
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
