"""
调试热应力为零的问题

此脚本用于验证温度场和热应力计算是否正确
"""

using LinearAlgebra, Statistics, Printf

println("="^80)
println("热应力调试工具")
println("="^80)

# 测试参数
α = 5e-6        # 1/K (热膨胀系数)
Ω = 3.1e-6      # m³/mol
cs_max = 33133  # mol/m³
β_eff = Ω * cs_max / 3.0  # 有效扩散应变系数
T_ref = 298.15  # K
E = 3e10        # Pa
ν = 0.29

println("\n材料参数:")
@printf("  α = %.2e 1/K\n", α)
@printf("  β_eff = %.4f (无量纲)\n", β_eff)
@printf("  T_ref = %.2f K\n", T_ref)
@printf("  E = %.1f GPa\n", E/1e9)
@printf("  ν = %.2f\n", ν)

println("\n" * "="^80)
println("场景1：仅热应力（无SOC变化）")
println("="^80)

ΔT = 30.0  # K，温升30度
Δsoc = 0.0  # 无SOC变化

ε_thermal = α * ΔT
ε_diffusion = β_eff * Δsoc
ε_total = ε_thermal + ε_diffusion

σ_approx = E / (1 - ν^2) * ε_total

println("\n输入:")
@printf("  ΔT = %.1f K\n", ΔT)
@printf("  Δsoc = %.4f\n", Δsoc)

println("\n计算的应变:")
@printf("  ε_thermal = %.4e\n", ε_thermal)
@printf("  ε_diffusion = %.4e\n", ε_diffusion)
@printf("  ε_total = %.4e\n", ε_total)

println("\n计算的应力:")
@printf("  σ ≈ %.2f MPa\n", σ_approx/1e6)

println("\n应力分离:")
if abs(ε_thermal) + abs(ε_diffusion) > 1e-15
    ratio_thermal = abs(ε_thermal) / (abs(ε_thermal) + abs(ε_diffusion))
    ratio_diffusion = abs(ε_diffusion) / (abs(ε_thermal) + abs(ε_diffusion))
    σ_thermal = ratio_thermal * σ_approx
    σ_diffusion = ratio_diffusion * σ_approx
    
    @printf("  σ_thermal = %.2f MPa (%.1f%%)\n", σ_thermal/1e6, ratio_thermal*100)
    @printf("  σ_diffusion = %.2f MPa (%.1f%%)\n", σ_diffusion/1e6, ratio_diffusion*100)
else
    println("  应变太小，应力为0")
end

println("\n✅ 预期：热应力 ≈ $(round(σ_approx/1e6, digits=1)) MPa")

println("\n" * "="^80)
println("场景2：仅扩散应力（无温度变化）")
println("="^80)

ΔT = 0.0    # 无温升
Δsoc = 0.1  # 10% SOC变化

ε_thermal = α * ΔT
ε_diffusion = β_eff * Δsoc
ε_total = ε_thermal + ε_diffusion

σ_approx = E / (1 - ν^2) * ε_total

println("\n输入:")
@printf("  ΔT = %.1f K\n", ΔT)
@printf("  Δsoc = %.4f\n", Δsoc)

println("\n计算的应变:")
@printf("  ε_thermal = %.4e\n", ε_thermal)
@printf("  ε_diffusion = %.4e\n", ε_diffusion)
@printf("  ε_total = %.4e\n", ε_total)

println("\n计算的应力:")
@printf("  σ ≈ %.2f MPa\n", σ_approx/1e6)

println("\n应力分离:")
if abs(ε_thermal) + abs(ε_diffusion) > 1e-15
    ratio_thermal = abs(ε_thermal) / (abs(ε_thermal) + abs(ε_diffusion))
    ratio_diffusion = abs(ε_diffusion) / (abs(ε_thermal) + abs(ε_diffusion))
    σ_thermal = ratio_thermal * σ_approx
    σ_diffusion = ratio_diffusion * σ_approx
    
    @printf("  σ_thermal = %.2f MPa (%.1f%%)\n", σ_thermal/1e6, ratio_thermal*100)
    @printf("  σ_diffusion = %.2f MPa (%.1f%%)\n", σ_diffusion/1e6, ratio_diffusion*100)
else
    println("  应变太小，应力为0")
end

println("\n✅ 预期：扩散应力 ≈ $(round(σ_approx/1e6, digits=1)) MPa")

println("\n" * "="^80)
println("场景3：热应力 + 扩散应力（典型放电情况）")
println("="^80)

ΔT = 30.0   # 温升30度
Δsoc = 0.1  # 10% SOC变化

ε_thermal = α * ΔT
ε_diffusion = β_eff * Δsoc
ε_total = ε_thermal + ε_diffusion

σ_approx = E / (1 - ν^2) * ε_total

println("\n输入:")
@printf("  ΔT = %.1f K\n", ΔT)
@printf("  Δsoc = %.4f\n", Δsoc)

println("\n计算的应变:")
@printf("  ε_thermal = %.4e (%.1f%%)\n", ε_thermal, 100*ε_thermal/ε_total)
@printf("  ε_diffusion = %.4e (%.1f%%)\n", ε_diffusion, 100*ε_diffusion/ε_total)
@printf("  ε_total = %.4e\n", ε_total)

println("\n计算的应力:")
@printf("  σ ≈ %.2f MPa\n", σ_approx/1e6)

println("\n应力分离:")
if abs(ε_thermal) + abs(ε_diffusion) > 1e-15
    ratio_thermal = abs(ε_thermal) / (abs(ε_thermal) + abs(ε_diffusion))
    ratio_diffusion = abs(ε_diffusion) / (abs(ε_thermal) + abs(ε_diffusion))
    σ_thermal = ratio_thermal * σ_approx
    σ_diffusion = ratio_diffusion * σ_approx
    
    @printf("  σ_thermal = %.2f MPa (%.1f%%)\n", σ_thermal/1e6, ratio_thermal*100)
    @printf("  σ_diffusion = %.2f MPa (%.1f%%)\n", σ_diffusion/1e6, ratio_diffusion*100)
    @printf("  σ_total = %.2f MPa\n", (σ_thermal+σ_diffusion)/1e6)
else
    println("  应变太小，应力为0")
end

println("\n✅ 预期：总应力 ≈ $(round(σ_approx/1e6, digits=1)) MPa")

println("\n" * "="^80)
println("诊断检查")
println("="^80)

println("""
如果testexample.jl中热应力为0，可能的原因：

1. ❌ 温度没有变化
   检查：T_avg_hist是否所有值相同
   解决：确保温度场有时间演化

2. ❌ 温度数据没有正确传递
   检查：variables_step["T_nodes"]是否正确设置
   解决：使用插值方法重建温度场历史

3. ❌ 热膨胀系数为0
   检查：α_eff的值
   解决：确保使用有量纲参数 param_dim

4. ❌ 应变分离逻辑错误
   检查：epsilon_thermal_elem的值
   解决：使用绝对值比例分配

当前修复：
✅ 使用平均温度历史插值重建节点温度场
✅ 改进应变分离逻辑（绝对值比例）
✅ 添加详细的调试输出
""")

println("="^80)
