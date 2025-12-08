"""
诊断热边界条件的无量纲化

检查：
1. 质量矩阵、刚度矩阵、载荷向量的无量纲化
2. 边界条件的无量纲化
3. h（Biot数）的使用是否正确
"""

using Printf

println("="^70)
println("热边界条件无量纲化诊断")
println("="^70)

# 参数
h = 150.0  # W/(m²·K) - 对流换热系数
k_ref = 2.1  # W/(m·K) - 参考热导率
L_th = 0.0105  # m - 特征长度
T_ref = 298.0  # K - 参考温度
rho_c_ref = 1.5e6  # J/(m³·K) - 体积热容

println("\n[1] 参数:")
@printf("  h = %.1f W/(m²·K)\n", h)
@printf("  k_ref = %.2f W/(m·K)\n", k_ref)
@printf("  L_th = %.4f m\n", L_th)
@printf("  T_ref = %.1f K\n", T_ref)

# 计算 Biot 数
Bi = h * L_th / k_ref
@printf("\n[2] Biot 数:\n")
@printf("  Bi = h × L_th / k_ref = %.3f\n", Bi)

println("\n[3] 无量纲化分析:")

println("\n  === 热传导方程 ===")
println("  有量纲: ρc ∂T/∂t = ∇·(k∇T) + q")
println("  无量纲: ∂T*/∂t* = ∇*·(k*∇*T*) + q*")
println()
println("  其中:")
println("    T* = T / T_ref")
println("    x* = x / L_th")
println("    t* = t × (k/(ρc×L²))")
println("    k* = k / k_ref")
println("    q* = q / (k_ref×T_ref/L_th²)")

println("\n  === 对流边界条件 ===")
println("  有量纲: -k ∂T/∂n = h (T - T_amb)")
println("  推导:")
println("    -k (T_ref/L_th) ∂T*/∂n* = h T_ref (T* - T_amb*)")
println("    -∂T*/∂n* = (h×L_th/k) (T* - T_amb*)")
println("    -∂T*/∂n* = Bi (T* - T_amb*)")
@printf("    -∂T*/∂n* = %.3f (T* - T_amb*)\n", Bi)

println("\n[4] 弱形式分析:")

println("\n  === 刚度矩阵（体积积分）===")
println("  有量纲: KT = ∬_Ω k ∇N^T ∇N dΩ")
println("  无量纲化:")
println("    设 x* = x/L, dΩ* = dΩ/L²")
println("    ∇* = L∇ (链式法则)")
println()
println("    KT* = ∬_Ω* (k/k_ref) (L∇N)^T (L∇N) dΩ*")
println("        = ∬_Ω (k/k_ref) L² (∇N)^T (∇N) (dΩ/L²)")
println("        = ∬_Ω (k/k_ref) (∇N)^T (∇N) dΩ")
println()
println("  离散形式:")
println("    KT* = Σ (k/k_ref) (∇N)^T (∇N) |J|")
println("    其中 |J| = w × detJ (有量纲, m²)")
println()
println("  代码中（ThermalDistributed.jl:291）:")
println("    weights = -λ_iso_nd .* wJ")
println("    ⚠️ 注意：wJ 是有量纲的（m²）")
println("    ⚠️ 但没有除以 L_th²！")

println("\n  === 质量矩阵（体积积分）===")
println("  有量纲: MT = ∬_Ω ρc N^T N dΩ")
println("  无量纲化:")
println("    MT* = ∬_Ω* (ρc/ρc_ref) N^T N dΩ*")
println("        = ∬_Ω (ρc/ρc_ref) N^T N (dΩ/L²)")
println()
println("  代码中（ThermalDistributed.jl:179）:")
println("    ρc_weights = ρc_e[ele_of_gp] .* (wJ ./ L_th^2)")
println("    ✓ 正确：wJ / L_th² 是无量纲的体积元")

println("\n  === 边界条件（线积分）===")
println("  弱形式: ∫_∂Ω h (T - T_amb) φ dS")
println("  无量纲化:")
println("    ∫_∂Ω* Bi (T* - T_amb*) φ dS*")
println("    其中 dS* = dS / L_th")
println()
println("  离散形式:")
println("    ∫ Bi (T* - T_amb*) φ (J/L_th) dξ")
println("    其中 J = |∂s/∂ξ| (有量纲边长/2, m)")
println()
println("  代码中（ThermalDistributed.jl:394）:")
println("    wt = Bi * w * (J / L_th)")
println("    ✓ 看起来正确：J/L_th 是无量纲的线元")

println("\n[5] 问题分析:")

println("\n  🔴 发现不一致！")
println()
println("  质量矩阵: wJ / L_th² ← 有无量纲化因子")
println("  刚度矩阵: wJ        ← 缺少无量纲化因子！")
println("  边界条件: J / L_th  ← 有无量纲化因子")
println()
println("  推测：")
println("  1. 如果刚度矩阵没有除以 L_th²，那么它不是真正的无量纲矩阵")
println("  2. 边界条件除以 L_th 是对的")
println("  3. 但如果刚度矩阵的尺度不对，边界条件的权重也会不匹配")

println("\n[6] 数值验证:")

# 假设一个简单的单元
elem_area = 1e-6  # m² - 单元面积
edge_length = 1e-3  # m - 边长

wJ_example = elem_area  # 体积积分权重
J_example = edge_length / 2  # 边界积分雅可比

println("\n  假设单元:")
@printf("    单元面积: %.2e m²\n", elem_area)
@printf("    边长: %.2e m\n", edge_length)

println("\n  刚度矩阵项（当前代码）:")
@printf("    weight = k* × wJ = 1.0 × %.2e = %.2e\n", wJ_example, wJ_example)
println("    单位: [?] (不清楚)")

println("\n  刚度矩阵项（如果正确无量纲化）:")
k_nd_weight_correct = wJ_example / L_th^2
@printf("    weight = k* × (wJ/L²) = 1.0 × %.2e = %.2e\n", k_nd_weight_correct, k_nd_weight_correct)
println("    单位: [无量纲]")

println("\n  边界条件项:")
bc_weight = Bi * J_example / L_th
@printf("    weight = Bi × (J/L) = %.3f × %.2e = %.2e\n", Bi, J_example/L_th, bc_weight)
println("    单位: [无量纲]")

println("\n  比例:")
ratio_current = wJ_example / bc_weight
ratio_correct = k_nd_weight_correct / bc_weight
@printf("    当前: K_weight / BC_weight = %.2e\n", ratio_current)
@printf("    正确: K_weight / BC_weight = %.2e\n", ratio_correct)
@printf("    差异因子: L_th² = %.2e\n", L_th^2)

println("\n[7] 推测的bug:")

println("""
如果刚度矩阵没有正确无量纲化，可能导致：
  1. 边界条件的相对权重不正确
  2. Biot 数的实际效果被缩放了 L_th² 倍
  3. 有效 Biot 数 = Bi × L_th² = %.3f × %.2e = %.2e
  4. 这会导致对流换热的效果被大幅放大！

数值影响：
  - 理论 Bi = %.3f (适中的对流)
  - 实际 Bi_eff ≈ %.3f × %.2e = %.2e (非常强的对流！)
  - 这可能解释为什么外圈温度被"锁定"在 T_amb
""" % (Bi, L_th^2, Bi * L_th^2, Bi, Bi, L_th^2, Bi * L_th^2))

println("\n[8] 验证方法:")

println("""
方法1: 修正刚度矩阵无量纲化
  在 ThermalDistributed.jl 第291行：
  ```julia
  # 原来
  weights = -λ_iso_nd .* wJ
  
  # 修正为
  weights = -λ_iso_nd .* (wJ ./ (L_th^2))
  ```

方法2: 修正边界条件权重
  在 ThermalDistributed.jl 第394行：
  ```julia
  # 原来
  wt = Bi * w * (J / L_th)
  
  # 修正为（补偿刚度矩阵的尺度）
  wt = Bi * w * J * L_th  # = Bi * w * J * L_th（增大L_th²倍）
  ```

推荐: 修正刚度矩阵（方法1），因为它应该与质量矩阵一致。
""")

println("\n[9] 检查清单:")
println("""
需要检查的代码位置：
  1. ThermalDistributed.jl:291 - 各向同性刚度矩阵
  2. ThermalDistributed.jl:242 - 各向异性刚度矩阵
  3. ThermalDistributed.jl:394 - 对流边界条件
  4. 确保所有矩阵的无量纲化一致
""")

println("="^70)
