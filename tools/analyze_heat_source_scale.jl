"""
分析thermal_example和testexample_simple_coupling的内热源尺度差异
"""

include("../src/JuBat.jl")

println("="^80)
println("内热源尺度参数分析")
println("="^80)

# ========================================================================
# 1. LG M50 电池（thermal_example使用）
# ========================================================================
println("\n[1] LG M50 电池参数:")
println("-"^80)

param_lg = JuBat.ChooseCell("LG M50")
println("基本参数:")
println("  cell.volume = $(param_lg.cell.volume) m³")
println("  cell.area = $(param_lg.cell.area) m²")
println("  cell.I1C = $(param_lg.cell.I1C) A")
println("  cell.Rout = $(param_lg.cell.Rout) m")
println("  cell.Rin = $(param_lg.cell.Rin) m")

println("\n尺度参数 (scale):")
println("  I_typ = $(param_lg.scale.I_typ) A")
println("  phi = $(param_lg.scale.phi) V")
println("  T_ref = $(param_lg.scale.T_ref) K")
println("  L = $(param_lg.scale.L) m (电化学特征长度)")
println("  L_th = $(param_lg.scale.L_th) m (热特征长度)")
println("  k_th = $(param_lg.scale.k_th) W/(m·K) (参考热导率)")
println("  rho_c_th = $(param_lg.scale.rho_c_th) J/(m³·K) (参考体积热容)")
println("  q_th = $(param_lg.scale.q_th) W/m³ (参考体积热源密度)")
println("  t_th = $(param_lg.scale.t_th) s (热扩散时间尺度)")

# 集总模型热源计算
# 根据 Thermal.jl 第62行: heat_internal_W = heat_internal_nd * I_typ * phi
heat_scale_lumped = param_lg.scale.I_typ * param_lg.scale.phi
println("\n集总模型热源尺度:")
println("  heat_scale = I_typ × phi = $(heat_scale_lumped) W")
println("  (热源单位: W, 总热源)")

# ========================================================================
# 2. Jellyroll 电池（testexample_simple_coupling使用）
# ========================================================================
println("\n\n[2] Jellyroll 电池参数:")
println("-"^80)

param_jr = JuBat.ChooseCell("Jellyroll")
println("基本参数:")
println("  cell.volume = $(param_jr.cell.volume) m³")
println("  cell.area = $(param_jr.cell.area) m²")
println("  cell.I1C = $(param_jr.cell.I1C) A")
println("  cell.Rout = $(param_jr.cell.Rout) m")
println("  cell.Rin = $(param_jr.cell.Rin) m")

println("\n尺度参数 (scale):")
println("  I_typ = $(param_jr.scale.I_typ) A")
println("  phi = $(param_jr.scale.phi) V")
println("  T_ref = $(param_jr.scale.T_ref) K")
println("  L = $(param_jr.scale.L) m (电化学特征长度)")
println("  L_th = $(param_jr.scale.L_th) m (热特征长度)")
println("  k_th = $(param_jr.scale.k_th) W/(m·K) (参考热导率)")
println("  rho_c_th = $(param_jr.scale.rho_c_th) J/(m³·K) (参考体积热容)")
println("  q_th = $(param_jr.scale.q_th) W/m³ (参考体积热源密度)")
println("  t_th = $(param_jr.scale.t_th) s (热扩散时间尺度)")

# 分布式模型热源计算
# 根据 CallModel_SimpleCoupling.jl 第93行: total_heat_W = q_ref * sum(q_elem .* elem_volumes)
println("\n简化耦合模型热源尺度:")
println("  q_ref = q_th = $(param_jr.scale.q_th) W/m³")
println("  total_heat_W = q_ref × ∑(q_elem × V_elem)")
println("  (热源单位: W/m³ → W, 需要乘以单元体积和进行积分)")

# ========================================================================
# 3. 关键差异分析
# ========================================================================
println("\n\n[3] 关键差异分析:")
println("="^80)

println("\n3.1 热源计算公式差异:")
println("  集总模型 (thermal_example):")
println("    heat_internal_W = heat_internal_nd × I_typ × phi")
println("    其中 heat_internal_nd 是无量纲热源")
println("    尺度因子 = $(heat_scale_lumped) W")
println()
println("  简化耦合模型 (testexample_simple_coupling):")
println("    total_heat_W = q_th × ∑(q_elem × V_elem)")
println("    其中 q_elem 是无量纲体积热源密度")
println("    尺度因子 = $(param_jr.scale.q_th) W/m³")

println("\n3.2 尺度因子比值:")
ratio1 = heat_scale_lumped / param_jr.scale.q_th
println("  (I_typ × phi) / q_th = $(heat_scale_lumped) / $(param_jr.scale.q_th)")
println("  = $(ratio1) m³")
println("  注意: 单位不匹配！ [W] / [W/m³] = [m³]")

println("\n3.3 体积比较:")
println("  LG M50 体积: $(param_lg.cell.volume) m³")
println("  Jellyroll 体积: $(param_jr.cell.volume) m³")
println("  体积比: $(param_lg.cell.volume / param_jr.cell.volume)")

println("\n3.4 q_th 的计算公式:")
println("  q_th = k_th × T_ref / L_th²")
println("  其中:")
println("    - L_th 是热特征长度 (通常取 Rout)")
println("    - k_th 是参考热导率")
println("    - T_ref 是参考温度")

println("\n  对于 LG M50:")
if param_lg.scale.L_th == 0.0
    println("    ⚠️ L_th = 0 (Rout未设置)")
    println("    这会导致 q_th 异常！")
else
    println("    L_th = $(param_lg.scale.L_th) m")
    println("    q_th = $(param_lg.scale.k_th) × $(param_lg.scale.T_ref) / $(param_lg.scale.L_th)²")
    println("        = $(param_lg.scale.q_th) W/m³")
end

println("\n  对于 Jellyroll:")
println("    L_th = $(param_jr.scale.L_th) m (= Rout)")
println("    q_th = $(param_jr.scale.k_th) × $(param_jr.scale.T_ref) / $(param_jr.scale.L_th)²")
println("        = $(param_jr.scale.q_th) W/m³")

println("\n3.5 根本原因分析:")
println("  ✓ LG M50使用集总模型:")
println("    - 热源是总功率 (W)")
println("    - 无量纲化尺度: I_typ × phi = $(heat_scale_lumped) W")
println("    - 计算: Q_total = (Q_rxn + Q_ohm + Q_rev) × I_typ × phi")
println()
println("  ✓ Jellyroll使用简化耦合模型:")
println("    - 热源是体积热源密度 (W/m³)")
println("    - 无量纲化尺度: q_th = $(param_jr.scale.q_th) W/m³")
println("    - 计算: Q_total = q_th × ∑(q_elem × V_elem)")
println()
println("  ⚠️ 关键问题:")
println("    1. 两个模型使用了不同的热源尺度")
println("    2. 集总模型: [W] = [无量纲] × [W]")
println("    3. 分布式模型: [W] = [无量纲] × [W/m³] × [m³]")
println("    4. 如果分布式模型的无量纲热源 q_elem 和集总模型的 heat_internal_nd")
println("       使用相同的方式计算，那么它们的物理意义不同:")
println("       - heat_internal_nd 应该是总热源的无量纲形式")
println("       - q_elem 应该是体积热源密度的无量纲形式")

# ========================================================================
# 4. 数值示例
# ========================================================================
println("\n\n[4] 数值示例 (假设相同的电化学条件):")
println("="^80)

# 假设无量纲热源值为 0.1 (典型值)
Q_nd = 0.1

println("假设无量纲热源 = $(Q_nd)")
println()
println("集总模型 (thermal_example):")
Q_lumped = Q_nd * heat_scale_lumped
println("  Q_total = $(Q_nd) × $(heat_scale_lumped) = $(Q_lumped) W")
println()
println("简化耦合模型 (testexample_simple_coupling):")
V_jr = param_jr.cell.volume
Q_distributed = Q_nd * param_jr.scale.q_th * V_jr
println("  假设热源均匀分布在整个体积")
println("  Q_total = $(Q_nd) × $(param_jr.scale.q_th) × $(V_jr)")
println("         = $(Q_distributed) W")
println()
println("热源比值:")
ratio_Q = Q_lumped / Q_distributed
println("  Q_lumped / Q_distributed = $(ratio_Q)")

println("\n" * "="^80)
println("结论:")
println("="^80)
println("两个例子计算得到的内热源数量级相差很大的根本原因是:")
println()
println("1. 使用了不同的热模型:")
println("   - thermal_example: 集总模型 (lumped)")
println("   - testexample_simple_coupling: 2D分布式模型 (distributed2D)")
println()
println("2. 不同的热源无量纲化方式:")
println("   - 集总模型: 以总功率为尺度 (I_typ × phi ≈ $(heat_scale_lumped) W)")
println("   - 分布式模型: 以体积热源密度为尺度 (q_th ≈ $(param_jr.scale.q_th) W/m³)")
println()
println("3. 尺度因子的巨大差异:")
println("   - q_th 的值严重依赖于 L_th (热特征长度)")
println("   - Jellyroll 的 L_th = Rout = $(param_jr.scale.L_th) m")
println("   - q_th ∝ 1/L_th² → L_th 越小，q_th 越大")
println()
println("4. 物理解释:")
println("   - 集总模型假设温度均匀，只关心总热源 (W)")
println("   - 分布式模型需要空间分辨的体积热源密度 (W/m³)")
println("   - 即使电化学条件相同，两者的无量纲热源物理意义不同")
println()
println("5. 建议修正:")
println("   - 如果要比较，应该将分布式模型的热源对体积积分")
println("   - 或者将集总模型的热源除以体积得到体积热源密度")
println("   - 确保比较的是相同物理量")
println("="^80)
