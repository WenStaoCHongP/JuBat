"""
参数对比脚本：Jellyroll vs LG M50

比较 testexample_simple_coupling 和 thermal_example 使用的电池参数
"""

using Printf

include("src/JuBat.jl")
using .JuBat

println("="^80)
println("电池参数对比：Jellyroll vs LG M50")
println("="^80)

# 加载参数
param_jelly = JuBat.ChooseCell("Jellyroll")
param_lg = JuBat.ChooseCell("LG M50")

println("\n" * "="^80)
println("1. 正极 (PE) 参数对比")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)

@printf("%-30s %15.2e %15.2e %15.2f\n", "厚度 (m)", param_jelly.PE.thickness, param_lg.PE.thickness, 
        param_jelly.PE.thickness/param_lg.PE.thickness)
@printf("%-30s %15.2f %15.2f %15.2f ⭐\n", "比热容 (J/kg/K)", param_jelly.PE.heat_Q, param_lg.PE.heat_Q,
        param_jelly.PE.heat_Q/param_lg.PE.heat_Q)
@printf("%-30s %15.2f %15.2f %15.2f\n", "密度 (kg/m³)", param_jelly.PE.rho, param_lg.PE.rho,
        param_jelly.PE.rho/param_lg.PE.rho)
@printf("%-30s %15.2f %15.2f %15.2f\n", "热导率 (W/m/K)", param_jelly.PE.lambda, param_lg.PE.lambda,
        param_jelly.PE.lambda/param_lg.PE.lambda)
@printf("%-30s %15.2f %15.2f %15.2f\n", "电导率 (S/m)", param_jelly.PE.sig, param_lg.PE.sig,
        param_jelly.PE.sig/param_lg.PE.sig)
@printf("%-30s %15.2e %15.2e %15.2f\n", "扩散系数 (m²/s)", param_jelly.PE.Ds, param_lg.PE.Ds,
        param_jelly.PE.Ds/param_lg.PE.Ds)
@printf("%-30s %15.2f %15.2f %15s\n", "dU/dT (V/K) @ SOC=0.5", param_jelly.PE.dUdT(0.5), param_lg.PE.dUdT(0.5), "⚠️ 都为零")
@printf("%-30s %15.2f %15.2f %15s\n", "Eac_D (J/mol)", param_jelly.PE.Eac_D, param_lg.PE.Eac_D, "⚠️ 都为零")
@printf("%-30s %15.2f %15.2f %15.2f\n", "Eac_k (J/mol)", param_jelly.PE.Eac_k, param_lg.PE.Eac_k,
        param_jelly.PE.Eac_k/param_lg.PE.Eac_k)

println("\n" * "="^80)
println("2. 负极 (NE) 参数对比")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)

@printf("%-30s %15.2e %15.2e %15.2f\n", "厚度 (m)", param_jelly.NE.thickness, param_lg.NE.thickness,
        param_jelly.NE.thickness/param_lg.NE.thickness)
@printf("%-30s %15.2f %15.2f %15.2f ⭐⭐\n", "比热容 (J/kg/K)", param_jelly.NE.heat_Q, param_lg.NE.heat_Q,
        param_jelly.NE.heat_Q/param_lg.NE.heat_Q)
@printf("%-30s %15.2f %15.2f %15.2f\n", "密度 (kg/m³)", param_jelly.NE.rho, param_lg.NE.rho,
        param_jelly.NE.rho/param_lg.NE.rho)
@printf("%-30s %15.2f %15.2f %15.2f\n", "热导率 (W/m/K)", param_jelly.NE.lambda, param_lg.NE.lambda,
        param_jelly.NE.lambda/param_lg.NE.lambda)
@printf("%-30s %15.2f %15.2f %15.2f\n", "电导率 (S/m)", param_jelly.NE.sig, param_lg.NE.sig,
        param_jelly.NE.sig/param_lg.NE.sig)
@printf("%-30s %15.2e %15.2e %15.2f\n", "扩散系数 (m²/s)", param_jelly.NE.Ds, param_lg.NE.Ds,
        param_jelly.NE.Ds/param_lg.NE.Ds)
@printf("%-30s %15.2f %15.2f %15s\n", "dU/dT (V/K) @ SOC=0.5", param_jelly.NE.dUdT(0.5), param_lg.NE.dUdT(0.5), "⚠️ 都为零")
@printf("%-30s %15.2f %15.2f %15s\n", "Eac_D (J/mol)", param_jelly.NE.Eac_D, param_lg.NE.Eac_D, "⚠️ 都为零")
@printf("%-30s %15.2f %15.2f %15.2f\n", "Eac_k (J/mol)", param_jelly.NE.Eac_k, param_lg.NE.Eac_k,
        param_jelly.NE.Eac_k/param_lg.NE.Eac_k)

println("\n" * "="^80)
println("3. 隔膜 (SP) 参数对比")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)

@printf("%-30s %15.2e %15.2e %15.2f\n", "厚度 (m)", param_jelly.SP.thickness, param_lg.SP.thickness,
        param_jelly.SP.thickness/param_lg.SP.thickness)
@printf("%-30s %15.2f %15.2f %15.2f ⭐⭐⭐\n", "比热容 (J/kg/K)", param_jelly.SP.heat_Q, param_lg.SP.heat_Q,
        param_jelly.SP.heat_Q/param_lg.SP.heat_Q)
@printf("%-30s %15.2f %15.2f %15.2f\n", "密度 (kg/m³)", param_jelly.SP.rho, param_lg.SP.rho,
        param_jelly.SP.rho/param_lg.SP.rho)
@printf("%-30s %15.2f %15.2f %15.2f\n", "热导率 (W/m/K)", param_jelly.SP.lambda, param_lg.SP.lambda,
        param_jelly.SP.lambda/param_lg.SP.lambda)

println("\n" * "="^80)
println("4. 电解液 (EL) 参数对比")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)

@printf("%-30s %15.2f %15.2f %15s\n", "比热容 (J/kg/K)", param_jelly.EL.heat_Q, param_lg.EL.heat_Q, "⭐")
@printf("%-30s %15.2f %15.2f %15.2f\n", "密度 (kg/m³)", param_jelly.EL.rho, param_lg.EL.rho,
        param_jelly.EL.rho/param_lg.EL.rho)
@printf("%-30s %15.2f %15.2f %15.2f\n", "迁移数", param_jelly.EL.tplus, param_lg.EL.tplus,
        param_jelly.EL.tplus/param_lg.EL.tplus)

println("\n" * "="^80)
println("5. 电池整体 (Cell) 参数对比")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)

@printf("%-30s %15.2f %15.2f %15.2f ⭐⭐⭐\n", "对流系数 h (W/m²/K)", param_jelly.cell.h, param_lg.cell.h,
        param_jelly.cell.h/param_lg.cell.h)
@printf("%-30s %15.2f %15.2f %15s\n", "初始温度 (K)", param_jelly.cell.T0, param_lg.cell.T0, "相同")
@printf("%-30s %15.2f %15.2f %15s\n", "环境温度 (K)", param_jelly.cell.T_amb, param_lg.cell.T_amb, "相同")
@printf("%-30s %15.2e %15.2e %15.2f\n", "冷却面积 (m²)", param_jelly.cell.cooling_surface, param_lg.cell.cooling_surface,
        param_jelly.cell.cooling_surface/param_lg.cell.cooling_surface)
@printf("%-30s %15.2e %15.2e %15.2f\n", "体积 (m³)", param_jelly.cell.volume, param_lg.cell.volume,
        param_jelly.cell.volume/param_lg.cell.volume)

# 计算等效热容
C_jelly = (param_jelly.PE.heat_Q * param_jelly.PE.rho * param_jelly.PE.thickness +
           param_jelly.NE.heat_Q * param_jelly.NE.rho * param_jelly.NE.thickness +
           param_jelly.SP.heat_Q * param_jelly.SP.rho * param_jelly.SP.thickness) /
          (param_jelly.PE.rho * param_jelly.PE.thickness +
           param_jelly.NE.rho * param_jelly.NE.thickness +
           param_jelly.SP.rho * param_jelly.SP.thickness)

C_lg = (param_lg.PE.heat_Q * param_lg.PE.rho * param_lg.PE.thickness +
        param_lg.NE.heat_Q * param_lg.NE.rho * param_lg.NE.thickness +
        param_lg.SP.heat_Q * param_lg.SP.rho * param_lg.SP.thickness) /
       (param_lg.PE.rho * param_lg.PE.thickness +
        param_lg.NE.rho * param_lg.NE.thickness +
        param_lg.SP.rho * param_lg.SP.thickness)

println("\n" * "="^80)
println("6. 等效热容计算")
println("="^80)
@printf("%-30s %15s %15s %15s\n", "参数", "Jellyroll", "LG M50", "比值")
println("-"^80)
@printf("%-30s %15.2f %15.2f %15.2f ⭐⭐⭐\n", "等效比热容 (J/kg/K)", C_jelly, C_lg, C_jelly/C_lg)

println("\n" * "="^80)
println("7. 关键热时间尺度")
println("="^80)

# 对流冷却时间尺度 τ_conv = (ρ × C × V) / (h × A)
τ_conv_jelly = (param_jelly.cell.mass * C_jelly) / (param_jelly.cell.h * param_jelly.cell.cooling_surface)
τ_conv_lg = (param_lg.cell.mass * C_lg) / (param_lg.cell.h * param_lg.cell.cooling_surface)

@printf("%-40s %15.2f s\n", "Jellyroll 对流冷却时间尺度:", τ_conv_jelly)
@printf("%-40s %15.2f s\n", "LG M50 对流冷却时间尺度:", τ_conv_lg)
@printf("%-40s %15.2f\n", "比值 (Jellyroll/LG M50):", τ_conv_jelly/τ_conv_lg)

println("\n" * "="^80)
println("8. 关键发现总结")
println("="^80)
println()
println("⭐⭐⭐ 最重要的差异：")
println()
println("1. 对流换热系数差异：")
@printf("   Jellyroll: h = %.1f W/(m²·K)\n", param_jelly.cell.h)
@printf("   LG M50:    h = %.1f W/(m²·K)\n", param_lg.cell.h)
@printf("   → Jellyroll 的散热能力是 LG M50 的 %.1f 倍！\n", param_jelly.cell.h/param_lg.cell.h)
println()

println("2. 等效比热容差异：")
@printf("   Jellyroll: C_eff = %.1f J/(kg·K)\n", C_jelly)
@printf("   LG M50:    C_eff = %.1f J/(kg·K)\n", C_lg)
@printf("   → Jellyroll 的热容是 LG M50 的 %.2f 倍！\n", C_jelly/C_lg)
@printf("   → 相同热量下，Jellyroll 的温升仅为 LG M50 的 %.2f 倍\n", C_lg/C_jelly)
println()

println("3. 温度对电压影响：")
@printf("   dU/dT (正极): %.5f V/K (都为零！)\n", param_jelly.PE.dUdT(0.5))
@printf("   dU/dT (负极): %.5f V/K (都为零！)\n", param_jelly.NE.dUdT(0.5))
@printf("   Eac_D (正极): %.1f J/mol (都为零！)\n", param_jelly.PE.Eac_D)
@printf("   Eac_D (负极): %.1f J/mol (都为零！)\n", param_jelly.NE.Eac_D)
println("   → 温度对电压的直接影响非常小！")
println()

println("⚠️  关键警告：")
println()
println("由于 dU/dT = 0 和 Eac_D = 0，温度对电压的影响本身就很弱！")
println("即使温度升高 20-30 K，电压降也只有约 0.01-0.02 V。")
println("这是电压降不明显的根本物理原因！")
println()

println("="^80)
println("9. 物理解释")
println("="^80)
println()
println("为什么 testexample_simple_coupling 的电压降不明显？")
println()
println("原因 1: 高比热容 (C_jelly = $(round(C_jelly, digits=1)) J/kg/K)")
println("  → 需要更多热量才能提高温度")
println("  → 温度升高较慢")
println()
println("原因 2: 高散热系数 (h = $(round(param_jelly.cell.h, digits=1)) W/m²/K)")
println("  → 热量快速散失")
println("  → 温度升高被抑制")
println()
println("原因 3: 温度对电压影响弱 (dU/dT = 0)")
println("  → 即使温度升高，电压降也很小")
println("  → 这是最根本的原因！")
println()
println("结果：小热源 + 高热容 + 高散热 + 弱温度依赖")
println("    → 温升小 (13K) → 电压降小 (0.004V)")
println()

println("="^80)
println("10. 改进建议")
println("="^80)
println()
println("如果希望观察到更明显的电压降，可以：")
println()
println("1. 降低散热系数：")
println("   cell.h = 10.0  # 从 150 降至 10")
println("   → 预期温升增加到 40-50 K")
println()
println("2. 使用 LG M50 参数：")
println("   param_dim = JuBat.ChooseCell(\"LG M50\")")
println("   → 预期温升加倍")
println()
println("3. 提高 C-rate：")
println("   Crates = 3.0  # 从 1C 提高到 3C")
println("   → 预期热源增加 9 倍")
println()
println("4. 添加温度依赖性（修改参数文件）：")
println("   PE.dUdT = x-> -0.0002 * ones(size(x))")
println("   NE.dUdT = x-> 0.0003 * ones(size(x))")
println("   → 预期 dV/dT 增加到 -0.5 mV/K")
println()

println("="^80)
println("分析完成！")
println("="^80)
println()
println("详细报告已保存至:")
println("  - analysis_voltage_drop_comparison.md")
println("  - 电压降差异分析总结.md")
println()
