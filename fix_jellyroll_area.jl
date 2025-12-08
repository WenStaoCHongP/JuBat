"""
诊断和修正 Jellyroll 面积参数

问题: cell.area 使用了螺旋展开后的总长度，导致面积过大 476 倍
结果: 电流密度过低，放电速率 0.0021C 而非 1C
"""

using Printf

include("src/JuBat.jl")
using .JuBat

println("="^80)
println("Jellyroll 面积参数诊断")
println("="^80)

# 加载当前参数
param = JuBat.ChooseCell("Jellyroll")

println("\n当前参数:")
println("-"^80)
@printf("  Rout (外半径):     %.6f m\n", param.cell.Rout)
@printf("  Rin (内半径):      %.6f m\n", param.cell.Rin)
@printf("  width (高度):      %.6f m\n", param.cell.width)
@printf("  length (定义值):   %.6f m  ⚠️ 螺旋展开长度\n", param.cell.length)
@printf("  area (当前):       %.6f m²  ⚠️ 过大\n", param.cell.area)
@printf("  I1C (标称电流):    %.2f A\n", param.cell.I1C)
@printf("  capacity (容量):   %.2f Ah\n", param.cell.capacity)

# 计算平均半径
R_avg = (param.cell.Rout + param.cell.Rin) / 2

# 计算单层周长（应该用这个作为 length）
L_single = 2 * π * R_avg

# 计算单层面积
A_single = L_single * param.cell.width

# 计算螺旋圈数
N_turns = param.cell.length / L_single

# 计算当前电流密度
j_current = param.cell.I1C / param.cell.area

# 计算修正后的电流密度
j_corrected = param.cell.I1C / A_single

println("\n分析结果:")
println("-"^80)
@printf("  平均半径 R_avg:         %.6f m\n", R_avg)
@printf("  单层周长:               %.6f m\n", L_single)
@printf("  单层面积:               %.6f m²  ✅ 合理\n", A_single)
@printf("  螺旋圈数 (估算):        %.0f 圈\n", N_turns)
@printf("  面积比 (当前/单层):     %.1f 倍\n", param.cell.area / A_single)

println("\n电流密度对比:")
println("-"^80)
@printf("  当前电流密度:           %.3f A/m²  ⚠️ 极低\n", j_current)
@printf("  修正后电流密度:         %.1f A/m²  ✅ 正常\n", j_corrected)
@printf("  比值:                   %.0f 倍\n", j_corrected / j_current)

println("\n放电速率对比:")
println("-"^80)
@printf("  名义 C-rate:            1C\n")
@printf("  当前实际 C-rate:        %.4fC  ⚠️ 极慢\n", j_current / j_corrected)
@printf("  修正后实际 C-rate:      1C  ✅ 正确\n")

println("\n1 小时放电后 SOC 变化:")
println("-"^80)
@printf("  当前:  SOC 变化 %.2f%%  ⚠️ 几乎不变\n", 100 * j_current / j_corrected)
@printf("  修正:  SOC 变化 100%%  ✅ 完全放电\n")

println("\n修正方案:")
println("="^80)
println()
println("方案 1: 修正 area (推荐)")
println("-"^40)
println("修改 src/parameters/Jellyroll.jl:")
println()
println("# 原来:")
println("cell.length = 938.7693138")
println("cell.area = cell.width * cell.length * cell.no_layers  # $(round(param.cell.area, digits=4)) m²")
println()
println("# 改为:")
println("R_avg = (cell.Rout + cell.Rin) / 2")
println("cell.length = 2 * π * R_avg  # 单层周长 = $(round(L_single, digits=6)) m")
println("cell.area = cell.width * cell.length  # 单层面积 = $(round(A_single, digits=6)) m²")
println()

println("\n方案 2: 修正 I1C")
println("-"^40)
@printf("cell.I1C = %.1f  # 从 %.1f 改为 %.1f\n", 
        param.cell.I1C * (param.cell.area / A_single),
        param.cell.I1C,
        param.cell.I1C * (param.cell.area / A_single))
println()

println("\n预期改进效果:")
println("="^80)
println()
println("修正前:")
println("  - 电流密度: $(round(j_current, digits=3)) A/m²")
println("  - 实际 C-rate: $(round(j_current/j_corrected, digits=4))C")
println("  - 1小时后 SOC: 27% → 26.79%")
println("  - 电压降: 0.004V")
println("  - 热源: 0.7W")
println("  - 温升: 13K")
println()
println("修正后:")
println("  - 电流密度: $(round(j_corrected, digits=1)) A/m²")
println("  - 实际 C-rate: 1C")
println("  - 1小时后 SOC: 27% → 0%")
println("  - 电压降: ~1.7V  ✅ 与 thermal_example 可比")
println("  - 热源: ~50-100W  (增加 $(round(Int, j_corrected/j_current))² 倍)")
println("  - 温升: ~100-200K")
println()

println("="^80)
println("建议:")
println("="^80)
println()
println("1. 备份原文件:")
println("   cp src/parameters/Jellyroll.jl src/parameters/Jellyroll.jl.backup")
println()
println("2. 应用修正 (方案 1):")
println("   修改 Jellyroll.jl 中的 cell.length 和 cell.area 计算")
println()
println("3. 重新运行:")
println("   julia example/testexample_simple_coupling.jl")
println()
println("4. 验证结果:")
println("   - 检查 1 小时后电压是否降到接近 2.5V")
println("   - 检查温升是否显著 (50-100K)")
println("   - 与 thermal_example 对比")
println()

println("="^80)
println("诊断完成！")
println("="^80)
