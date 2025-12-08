"""
分析 Jellyroll 参数的归一化
"""

using Printf

include("src/JuBat.jl")
using .JuBat

println("="^80)
println("Jellyroll 参数归一化分析")
println("="^80)

# 加载参数
param_dim = JuBat.ChooseCell("Jellyroll")

println("\n原始参数 (有量纲):")
println("-"^80)
@printf("  I1C (标称电流):          %.2f A\n", param_dim.cell.I1C)
@printf("  capacity (容量):         %.2f Ah\n", param_dim.cell.capacity)
@printf("  cell.area (面积):        %.6f m²\n", param_dim.cell.area)
@printf("  cell.length (长度):      %.6f m\n", param_dim.cell.length)
@printf("  cell.width (宽度):       %.6f m\n", param_dim.cell.width)
@printf("  cell.volume (体积):      %.9f m³\n", param_dim.cell.volume)

println("\n电极厚度:")
@printf("  PE.thickness:            %.2e m\n", param_dim.PE.thickness)
@printf("  NE.thickness:            %.2e m\n", param_dim.NE.thickness)
@printf("  SP.thickness:            %.2e m\n", param_dim.SP.thickness)
@printf("  L (总厚度):              %.2e m\n", 
        param_dim.PE.thickness + param_dim.NE.thickness + param_dim.SP.thickness)

println("\n特征尺度:")
@printf("  scale.L:                 %.2e m\n", param_dim.scale.L)
@printf("  scale.r0:                %.2e m\n", param_dim.scale.r0)
@printf("  scale.a0 (1/r0):         %.2e m⁻¹\n", param_dim.scale.a0)
@printf("  scale.I_typ:             %.2f A\n", param_dim.scale.I_typ)
@printf("  scale.t0:                %.0f s\n", param_dim.scale.t0)
@printf("  scale.T_ref:             %.0f K\n", param_dim.scale.T_ref)

println("\n归一化电流密度尺度:")
@printf("  scale.j = I_typ/(a0×L×A)\n")
@printf("  scale.j = %.2f/(%.2e × %.2e × %.6f)\n", 
        param_dim.scale.I_typ, param_dim.scale.a0, param_dim.scale.L, param_dim.cell.area)
@printf("  scale.j = %.6e A/(m²·m⁻¹) = %.6e A/m\n", param_dim.scale.j, param_dim.scale.j)

println("\n容量时间尺度:")
@printf("  ts_p = F×cs_max×A×L/I_typ\n")
@printf("  ts_p = %.0f × %.0f × %.6f × %.2e / %.2f\n",
        param_dim.scale.F, param_dim.PE.cs_max, param_dim.cell.area, param_dim.scale.L, param_dim.scale.I_typ)
@printf("  ts_p = %.2e s = %.2f h\n", param_dim.scale.ts_p, param_dim.scale.ts_p/3600)

@printf("\n  ts_n = %.2e s = %.2f h\n", param_dim.scale.ts_n, param_dim.scale.ts_n/3600)

println("\n实际容量计算:")
Q_p = param_dim.scale.F * param_dim.PE.cs_max * param_dim.cell.area * param_dim.scale.L / 3600
Q_n = param_dim.scale.F * param_dim.NE.cs_max * param_dim.cell.area * param_dim.scale.L / 3600
@printf("  Q_p = F×cs_max×A×L/3600\n")
@printf("  Q_p = %.2f Ah (正极容量)\n", Q_p)
@printf("  Q_n = %.2f Ah (负极容量)\n", Q_n)
@printf("  Q_min = %.2f Ah (限制容量)\n", min(Q_p, Q_n))
@printf("  capacity (设定) = %.2f Ah\n", param_dim.cell.capacity)
@printf("  比值 = %.1f 倍\n", min(Q_p, Q_n) / param_dim.cell.capacity)

println("\n理论放电时间 (1C):")
t_discharge = min(Q_p, Q_n) / param_dim.cell.I1C
@printf("  t = Q/I1C = %.2f Ah / %.2f A\n", min(Q_p, Q_n), param_dim.cell.I1C)
@printf("  t = %.2f h = %.0f s\n", t_discharge, t_discharge * 3600)

println("\n"*"="^80)
println("关键发现:")
println("="^80)

if abs(min(Q_p, Q_n) - param_dim.cell.capacity) > 0.1
    println("\n⚠️  计算容量与设定容量不匹配！")
    @printf("  计算: %.2f Ah\n", min(Q_p, Q_n))
    @printf("  设定: %.2f Ah\n", param_dim.cell.capacity)
    @printf("  差异: %.1f 倍\n", min(Q_p, Q_n) / param_dim.cell.capacity)
    println("\n这意味着 cell.area 的定义可能确实是物理意义上的")
    println("'展开总面积'，而不是'单层等效面积'！")
end

println("\n"*"="^80)
println("物理解释:")
println("="^80)
println()
println("在归一化方程中，cell.area 的作用:")
println("1. 定义容量: Q = F × cs_max × A × L")
println("2. 定义时间尺度: ts = Q / I = (F × cs_max × A × L) / I")
println("3. 定义电流密度尺度: j = I / (a0 × L × A)")
println()
println("如果 A 很大 (48.95 m²)，则:")
println("- 容量很大 ($(round(min(Q_p, Q_n), digits=1)) Ah)")
println("- 时间尺度很大 ($(round(param_dim.scale.ts_p/3600, digits=1)) h)")
println("- 电流密度很小 ($(round(param_dim.scale.j, sigdigits=3)) A/m)")
println()
println("这与观察到的现象一致:")
println("- 1 小时只消耗了很小比例的容量")
println("- 电压几乎不变")
