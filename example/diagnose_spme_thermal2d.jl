# 诊断脚本：分析 spme_thermal2d_example.jl 的温度异常
# 问题：给定恒定热源，但温度出现不正常（可能下降）

using LinearAlgebra, SparseArrays, Statistics
include("../src/JuBat.jl")

println("="^80)
println("spme_thermal2d_example.jl 温度异常诊断")
println("="^80)

# 1) 复现问题设置
param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
opt = JuBat.Option()
Crates = 1
i = 5*Crates
opt.Current = x-> i 
opt.model = "SPMe" 
opt.Nn = 20; opt.Ns = 10; opt.Np = 20
opt.Nrn = 15; opt.Nrp = 15
opt.gsorder = 2
opt.dimension = 1
opt.time = [0.0 50]
opt.dt = [0.5, 0.5]
opt.jacobi = "update"
opt.solveType = "backward" 
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"

electrochem_enabled = false  # 纯热模拟
case = JuBat.SetCase(param_dim, opt)

# 2) 创建热网格
mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case.mesh["thermal2D"] = mesh_th
ne = size(mesh_th.element, 1)
nnode = mesh_th.nlen

println("\n网格信息:")
println("  单元数: $ne")
println("  节点数: $nnode")

# 3) 几何参数
Rin = getfield(param_dim.cell, :Rin)
Rout = getfield(param_dim.cell, :Rout)
width = param_dim.cell.width
h_conv = param_dim.cell.h
T_amb = param_dim.cell.T_amb

println("\n几何与物性:")
println("  内径 Rin: $(Rin*1000) mm")
println("  外径 Rout: $(Rout*1000) mm")
println("  宽度 width: $(width*1000) mm")
println("  对流换热系数 h: $(h_conv) W/(m²·K)")
println("  环境温度 T_amb: $(T_amb) K")

# 4) 计算对流面积
A_conv = 2π * Rout * width + 2 * π * Rout^2
A_conv_outer_circumference = 2π * Rout * width  # 仅外周
A_conv_end_caps = 2 * π * Rout^2  # 两端盖

println("\n对流散热面积:")
println("  外周面积: $(A_conv_outer_circumference*1e4) cm²")
println("  端盖面积: $(A_conv_end_caps*1e4) cm²")
println("  总面积: $(A_conv*1e4) cm²")

# 5) 计算体积
V_total = π * (Rout^2 - Rin^2) * width
println("\n电池体积:")
println("  总体积: $(V_total*1e6) cm³")

# 6) 热源设置
q_user = 1000.0  # W/m³（默认值）
println("\n热源设置:")
println("  体积热源密度 q: $(q_user) W/m³")

# 7) 计算总产热功率
Q_gen_total = q_user * V_total
println("\n总产热功率:")
println("  Q_gen = q × V = $(round(Q_gen_total, digits=4)) W")

# 8) 估算稳态温升
# 稳态时：Q_gen = h * A_conv * (T_ss - T_amb)
# => T_ss = T_amb + Q_gen / (h * A_conv)
ΔT_ss = Q_gen_total / (h_conv * A_conv)
T_ss = T_amb + ΔT_ss

println("\n稳态温升估算（忽略热传导）:")
println("  ΔT_ss = Q_gen / (h × A) = $(round(ΔT_ss, digits=4)) K")
println("  T_ss = T_amb + ΔT_ss = $(round(T_ss, digits=4)) K")
println("  相对温升: $(round(ΔT_ss, digits=4)) K")

if ΔT_ss < 0.1
    println("\n⚠️  警告：稳态温升仅 $(round(ΔT_ss, digits=4)) K！")
    println("   散热能力远大于产热，温度几乎不会升高。")
end

# 9) 无量纲Biot数
scale = case.param_dim.scale
L_th = scale.L_th
k_ref = scale.k_th
Bi = scale.h_th

println("\n无量纲参数:")
println("  特征长度 L_th: $(round(L_th*1000, digits=4)) mm")
println("  参考热导率 k_ref: $(round(k_ref, digits=4)) W/(m·K)")
println("  Biot数 Bi = h×L_th/k_ref: $(round(Bi, digits=6))")

if Bi < 0.1
    println("  → Bi < 0.1：热传导主导，温度场均匀")
elseif Bi > 10
    println("  → Bi > 10：对流主导，边界温度接近环境温度")
else
    println("  → 0.1 < Bi < 10：传导与对流竞争")
end

# 10) 时间尺度分析
t_th = scale.t_th
t0 = scale.t0
println("\n时间尺度:")
println("  热扩散时间尺度 t_th: $(round(t_th, digits=4)) s")
println("  电化学时间尺度 t0: $(round(t0, digits=4)) s")
println("  比值 t0/t_th: $(round(t0/t_th, digits=4))")

# 11) 计算对流散热功率随温度的变化
T_range = T_amb .+ (0:0.1:5)  # T_amb + 0 到 5K
Q_conv_range = h_conv * A_conv .* (T_range .- T_amb)

println("\n对流散热功率（在不同温升下）:")
for (i, T) in enumerate(T_range)
    if (i-1) % 10 == 0  # 每隔1K打印一次
        ΔT = T - T_amb
        Q_conv = Q_conv_range[i]
        println("  ΔT = $(round(ΔT, digits=1)) K → Q_conv = $(round(Q_conv, digits=4)) W")
    end
end

# 12) 能量平衡判断
println("\n能量平衡分析:")
println("  产热功率 Q_gen: $(round(Q_gen_total, digits=4)) W")
println("  在 ΔT = 1K 时的散热: $(round(h_conv * A_conv * 1.0, digits=4)) W")
println("  在 ΔT = 0.1K 时的散热: $(round(h_conv * A_conv * 0.1, digits=4)) W")

if Q_gen_total < h_conv * A_conv * 0.1
    println("\n❌ 问题诊断：")
    println("   Q_gen ($(round(Q_gen_total, digits=4)) W) << Q_conv@ΔT=0.1K ($(round(h_conv * A_conv * 0.1, digits=4)) W)")
    println("   散热能力远超产热，温度升高极其微小，可能在数值噪声范围内！")
    println("\n   可能的\"温度不正常\"表现：")
    println("   1. 温度几乎不变（预期升高但实际不变）")
    println("   2. 温度轻微波动（数值误差）")
    println("   3. 温度略微下降（如果初始温度略高于环境温度）")
elseif Q_gen_total < h_conv * A_conv * 1.0
    println("\n⚠️  注意：")
    println("   Q_gen ($(round(Q_gen_total, digits=4)) W) < Q_conv@ΔT=1K ($(round(h_conv * A_conv * 1.0, digits=4)) W)")
    println("   稳态温升将小于1K，温度变化很小。")
else
    println("\n✓ 能量平衡合理：产热功率足够大，预期有明显温升。")
end

# 13) 建议
println("\n"*"="^80)
println("诊断建议:")
println("="^80)

if Q_gen_total < h_conv * A_conv * 0.1
    println("\n1. 问题根源：**散热能力远超产热能力**")
    println("   - 对流换热系数 h = $(h_conv) W/(m²·K) 过大")
    println("   - 热源密度 q = $(q_user) W/m³ 过小")
    println("   - 稳态温升仅 $(round(ΔT_ss, digits=4)) K")
    
    println("\n2. 解决方案：")
    println("   方案A：减小对流换热系数")
    println("      param_dim.cell.h = 10  # 从150降至10 W/(m²·K)")
    println("      预期温升：$(round(Q_gen_total / (10 * A_conv), digits=2)) K")
    
    println("\n   方案B：增大热源密度")
    q_new = h_conv * A_conv * 10 / V_total  # 使ΔT_ss = 10K
    println("      q_user = $(round(q_new, digits=0))  # 从1000增至 $(round(q_new, digits=0)) W/m³")
    println("      预期温升：10 K")
    
    println("\n   方案C：使用自然对流边界")
    println("      param_dim.cell.h = 5  # 自然对流")
    println("      预期温升：$(round(Q_gen_total / (5 * A_conv), digits=2)) K")
else
    println("\n能量平衡合理。如果仍出现温度异常，检查：")
    println("  1. 初始温度是否高于环境温度")
    println("  2. 边界条件是否施加了额外的强制冷却")
    println("  3. 热源单位是否正确转换（SI vs 无量纲）")
end

println("\n"*"="^80)
println("诊断完成")
println("="^80)
