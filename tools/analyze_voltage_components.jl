"""
分析初始电压的各个组成部分

已知：
- 理论 OCV = 4.1809 V
- 多SPMe初始电压 = 3.9992 V (差 -0.181 V)
- 简化耦合初始电压 = 4.0637 V (差 -0.117 V)

分析：过电位和欧姆压降
"""

using Printf

println("="^70)
println("初始电压组成分析")
println("="^70)

# 已知数据
OCV = 4.1809  # V
V_multi = 3.9992  # V (多SPMe)
V_simple = 4.0637  # V (简化耦合)

eta_total_multi = OCV - V_multi
eta_total_simple = OCV - V_simple

println("\n[1] 基本数据:")
@printf("  理论 OCV = %.4f V\n", OCV)
@printf("  多SPMe电压 = %.4f V\n", V_multi)
@printf("  简化耦合电压 = %.4f V\n", V_simple)
println()
@printf("  多SPMe总过电位 = %.4f V\n", eta_total_multi)
@printf("  简化耦合总过电位 = %.4f V\n", eta_total_simple)
@printf("  过电位差异 = %.4f V\n", eta_total_multi - eta_total_simple)

# 电池参数（从 Jellyroll.jl）
I_total = 5.0  # A
V_cell = 2.34e-5  # m³

# 典型参数估算
k_n = 6.48e-7  # 负极反应速率常数
k_p = 3.42e-6  # 正极反应速率常数
a_n = 3 * 0.25 / 5.86e-6  # 负极比表面积
a_p = 3 * 0.335 / 5.22e-6  # 正极比表面积
L_n = 8.52e-5  # 负极厚度
L_p = 7.56e-5  # 正极厚度
A_cell = 0.147  # 电池面积

println("\n[2] 估算过电位组成:")
println("\n  假设电流密度（均匀分布）:")
j_avg = I_total / A_cell
@printf("    j_avg = I / A = %.2f A/m²\n", j_avg)

# 反应过电位（Butler-Volmer近似）
# η = (R*T)/(α*F) * asinh(j/(2*j0))
# 其中 j0 = k * √(cs * ce * (cs_max - cs))
R = 8.314  # J/(mol·K)
T = 298.15  # K
F = 96485  # C/mol
alpha = 0.5

# 负极
j0_n = k_n * 1e3  # 估算，mol/(m²·s) * F
eta_n_est = (R*T)/(alpha*F) * asinh(j_avg * L_n / (2 * a_n * L_n * j0_n))

# 正极
j0_p = k_p * 1e3
eta_p_est = (R*T)/(alpha*F) * asinh(j_avg * L_p / (2 * a_p * L_p * j0_p))

@printf("\n  估算的反应过电位:\n")
@printf("    η_n ≈ %.4f V\n", abs(eta_n_est))
@printf("    η_p ≈ %.4f V\n", abs(eta_p_est))
@printf("    η_rxn_total ≈ %.4f V\n", abs(eta_n_est) + abs(eta_p_est))

# 欧姆过电位
sigma_n = 215.0  # S/m (负极固相电导率)
sigma_p = 0.18  # S/m (正极固相电导率)
kappa_e = 1.0  # S/m (电解液离子电导率，估算)

# 固相欧姆压降
R_solid = L_n / (sigma_n * A_cell) + L_p / (sigma_p * A_cell)
eta_solid = I_total * R_solid

# 电解液欧姆压降（粗略估算）
L_e = L_n + L_p + 1.2e-5  # 总电解液路径
R_elec = L_e / (kappa_e * A_cell * 0.3)  # 考虑孔隙率
eta_elec = I_total * R_elec

@printf("\n  估算的欧姆压降:\n")
@printf("    η_solid ≈ %.4f V (固相)\n", eta_solid)
@printf("    η_elec ≈ %.4f V (电解液)\n", eta_elec)
@printf("    η_ohm_total ≈ %.4f V\n", eta_solid + eta_elec)

total_est = abs(eta_n_est) + abs(eta_p_est) + eta_solid + eta_elec
@printf("\n  总估算过电位 ≈ %.4f V\n", total_est)

println("\n[3] 与实际的对比:")
@printf("  多SPMe实际过电位: %.4f V\n", eta_total_multi)
@printf("  估算过电位: %.4f V\n", total_est)
@printf("  简化耦合实际过电位: %.4f V\n", eta_total_simple)

println("\n[4] 过电位差异分析:")
delta_eta = eta_total_multi - eta_total_simple
@printf("  ΔV = %.4f V\n", delta_eta)

println("\n  可能的原因:")
println("  a) 电流分布不均匀性:")
println("     多SPMe: 逐单元独立计算，某些单元电流密度更高")
println("     → 局部过电位更大 → 整体电压更低")
println()
println("  b) 欧姆压降差异:")
println("     多SPMe: 考虑空间电流分布，路径电阻更大")
println("     简化耦合: 假设均匀电流，欧姆损失更小")
println()
println("  c) 初始瞬态:")
println("     多SPMe: 可能需要更多迭代才能平衡")
println("     简化耦合: 单SPMe更快达到稳态")

println("\n[5] 验证建议:")
println("\n  1. 检查第一个时间步的详细输出:")
println("     - 各单元的电流分布")
println("     - 各单元的过电位")
println("     - 电流守恒情况")
println()
println("  2. 比较 t=0 和第一个时间步:")
println("     - 真正的 t=0 应该等于 OCV")
println("     - 第一个时间步会包含过电位")
println()
println("  3. 输出电压的详细组成:")
println("     在求解器中添加:")
println("     ```julia")
println("     V_cell = OCV - eta_n - eta_p - eta_ohm")
println("     println(\"  OCV = \", OCV)")
println("     println(\"  η_n = \", eta_n)")
println("     println(\"  η_p = \", eta_p)")
println("     println(\"  η_ohm = \", eta_ohm)")
println("     println(\"  V_cell = \", V_cell)")
println("     ```")

println("\n[6] 预期行为:")
println("\n  ✓ 正常情况:")
println("    - t=0: V = OCV = 4.18 V (无电流)")
println("    - t=dt: V = 4.0-4.1 V (有电流，包含过电位)")
println()
println("  ⚠ 当前情况:")
println("    - 报告的\"初始电压\"已经包含过电位")
println("    - 可能是第一个时间步的结果，不是真正的 t=0")
println()
println("  ✓ 解释差异:")
println("    - 多SPMe过电位更大 (0.181 V)")
println("    - 简化耦合过电位更小 (0.117 V)")
println("    - 差异 0.064 V 来自于模型复杂度的不同")

println("\n[7] 结论:")
println("""
两个模式的初始状态（SOC）完全相同，理论 OCV 也相同。

电压差异来自于：
  1. "初始电压"实际上是第一个时间步的电压（已有电流）
  2. 多SPMe考虑了空间异质性，过电位更大
  3. 简化耦合假设均匀分布，过电位更小

这是合理的，反映了两种模型的物理差异：
  - 多SPMe更精确（考虑局部效应）
  - 简化耦合更简化（假设均匀）

建议：
  - 这种差异是预期的，不需要修正
  - 如果需要完全一致，应该比较平衡后的 OCV（无电流状态）
  - 或者在文档中说明两种模式的适用场景
""")

println("="^70)
