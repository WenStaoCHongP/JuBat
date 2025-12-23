"""
对比颗粒应力和宏观应力的数值和物理意义

此脚本用于理解SPMe中两个尺度的力学模型
"""

using Printf

println("="^80)
println("SPMe力学模型：颗粒应力 vs 宏观应力")
println("="^80)

# ============================================================================
# 1. 颗粒级应力（微观）
# ============================================================================
println("\n" * "="^80)
println("1. 颗粒级应力（Calstressdisp）")
println("="^80)

# 典型参数
Ω_n = 3.1e-6        # m³/mol
E_n = 2.0e10        # Pa = 20 GPa
ν_n = 0.28          # -
r_s = 5.86e-6       # m

# 典型浓度变化
c_max = 33133       # mol/m³
Δc = 0.3 * c_max    # 30% SOC变化
c_avg = 0.6 * c_max
c_center = 0.7 * c_max
c_surf = 0.5 * c_max

# 计算颗粒应力
σ_r_center = (2 * Ω_n * E_n * (c_avg - c_center)) / (9 * (1 - ν_n))
σ_θ_surf = (Ω_n * E_n * (c_avg - c_surf)) / (3 * (1 - ν_n))
u_surf = (Ω_n * r_s * c_avg) / 3

println("\n输入参数:")
println("  Ω = $(Ω_n*1e6) × 10⁻⁶ m³/mol")
println("  E = $(E_n/1e9) GPa")
println("  ν = $(ν_n)")
println("  r_s = $(r_s*1e6) μm")
println("  Δc = $(Δc) mol/m³ ($(round(100*Δc/c_max))% SOC变化)")

println("\n颗粒内浓度分布:")
println("  c_avg = $(round(c_avg)) mol/m³")
println("  c_center = $(round(c_center)) mol/m³")
println("  c_surf = $(round(c_surf)) mol/m³")

println("\n计算的颗粒应力:")
@printf("  σ_r(中心) = %.1f MPa\n", σ_r_center/1e6)
@printf("  σ_θ(表面) = %.1f MPa\n", σ_θ_surf/1e6)
@printf("  u(表面) = %.3f nm\n", u_surf*1e9)

println("\n物理意义:")
println("  ✓ 描述单个颗粒内部的应力状态")
println("  ✓ 球形对称，解析解")
println("  ✓ 浓度梯度驱动")
println("  ✓ 影响电化学反应（修正过电位）")

println("\n代表性:")
println("  SPMe: 1个颗粒代表整个负极")
println("  假设: 所有颗粒行为相同")
println("  适用: 均匀电极")

# ============================================================================
# 2. 宏观应力（介观/宏观）
# ============================================================================
println("\n" * "="^80)
println("2. 宏观应力（thermal_diffusion_stress_2D）")
println("="^80)

# 典型参数
α_eff = 5e-6        # 1/K（混合热膨胀系数）
β_n_eff = Ω_n * c_max / 3  # 无量纲（有效扩散应变系数）
E_eff = 3.0e10      # Pa = 30 GPa（混合模量）
ν_eff = 0.29        # -

# 典型变化
ΔT_nd = 0.1         # 无量纲
T_ref = 298.15      # K
ΔT = ΔT_nd * T_ref  # K
Δsoc = 0.1          # 无量纲（10% SOC变化）

# 计算应变
ε_thermal = α_eff * ΔT
ε_diffusion = β_n_eff * Δsoc
ε_total = ε_thermal + ε_diffusion

# 计算应力（平面应力近似）
σ_macro = E_eff / (1 - ν_eff^2) * ε_total

println("\n输入参数:")
println("  α_eff = $(α_eff*1e6) × 10⁻⁶ 1/K")
println("  β_n_eff = $(round(β_n_eff, digits=4))")
println("  E_eff = $(E_eff/1e9) GPa")
println("  ν_eff = $(ν_eff)")

println("\n空间分布:")
println("  ΔT = $(round(ΔT, digits=1)) K (温度变化)")
println("  Δsoc = $(round(Δsoc*100))% (SOC变化)")

println("\n计算的应变:")
@printf("  ε_thermal = %.2e (热应变)\n", ε_thermal)
@printf("  ε_diffusion = %.2e (扩散应变)\n", ε_diffusion)
@printf("  ε_total = %.2e (总应变)\n", ε_total)

println("\n计算的宏观应力:")
@printf("  σ_macro ≈ %.1f MPa\n", σ_macro/1e6)

println("\n物理意义:")
println("  ✓ 描述电极/电池整体的应力场")
println("  ✓ 2D/3D复杂几何，有限元求解")
println("  ✓ 温度和SOC空间分布驱动")
println("  ✓ 用于机械失效分析")

println("\n代表性:")
println("  多SPMe: 每个热单元独立的SOC和温度")
println("  考虑: 空间非均匀性")
println("  适用: 复杂边界和载荷")

# ============================================================================
# 3. 对比分析
# ============================================================================
println("\n" * "="^80)
println("3. 对比分析")
println("="^80)

println("\n数值对比:")
@printf("  颗粒应力: %.1f MPa (颗粒中心/表面)\n", max(abs(σ_r_center), abs(σ_θ_surf))/1e6)
@printf("  宏观应力: %.1f MPa (电极整体)\n", σ_macro/1e6)
@printf("  比值: %.2f\n", max(abs(σ_r_center), abs(σ_θ_surf)) / σ_macro)

println("\n")
println("┌─────────────────┬──────────────────────┬──────────────────────┐")
println("│ 特性            │ 颗粒应力             │ 宏观应力             │")
println("├─────────────────┼──────────────────────┼──────────────────────┤")
println("│ 尺度            │ 微米（颗粒内）       │ 厘米（电池整体）     │")
println("│ 物理机制        │ 浓度梯度             │ 温度+SOC分布         │")
println("│ 几何            │ 球形对称             │ 复杂2D/3D            │")
println("│ 求解方法        │ 解析解               │ 有限元               │")
println("│ 边界条件        │ 自由表面             │ 固定边界             │")
@printf("│ 典型值          │ %-20s │ %-20s │\n", "$(round(Int,max(abs(σ_r_center), abs(σ_θ_surf))/1e6)) MPa", "$(round(Int,σ_macro/1e6)) MPa")
println("│ 影响            │ 电化学反应           │ 机械失效             │")
println("│ SPMe表示        │ ✅ 单颗粒可以        │ ❌ 需要空间分布      │")
println("└─────────────────┴──────────────────────┴──────────────────────┘")

# ============================================================================
# 4. 能否用单颗粒表示？
# ============================================================================
println("\n" * "="^80)
println("4. SPMe能否用单颗粒表示电极力学场？")
println("="^80)

println("\n✅ 可以表示：")
println("  1. 颗粒内部应力状态")
println("  2. 电化学-力学耦合效应（过电位修正）")
println("  3. 均匀电极的平均应力（近似）")

println("\n❌ 无法表示：")
println("  1. 电极的空间非均匀应力场")
println("  2. 边界效应和应力集中")
println("  3. 宏观几何的影响")
println("  4. 集流体约束的作用")

println("\n💡 当前实现（两个独立模型）：")
println("  • 颗粒应力 → 用于电化学计算（修正η）")
println("  • 宏观应力 → 用于力学分析（失效预测）")
println("  • 两者目前未直接联系")

println("\n🔧 改进方向：")
println("  1. 短期：保持现状，明确文档")
println("  2. 中期：建立均匀化联系")
println("  3. 长期：完整多尺度耦合")

println("\n" * "="^80)
println("结论")
println("="^80)
println("""
SPMe的单颗粒：
• ✅ 可以表示颗粒内应力
• ✅ 可以近似表示均匀电极
• ❌ 无法完整描述非均匀宏观应力场

建议：
• 保持当前的双尺度独立模型
• 根据需求选择合适的模型
• 未来可考虑多尺度耦合
""")

println("="^80)
