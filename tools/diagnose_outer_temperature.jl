"""
诊断最外圈温度始终为298K的问题

分析可能的原因：
1. 对流换热系数太大
2. 边界条件设置错误
3. 热源分布问题
4. 仿真时间太短
"""

using Printf

println("="^70)
println("最外圈温度诊断")
println("="^70)

# 从 Jellyroll 参数
println("\n[1] Jellyroll 参数:")
h = 150.0  # W/(m²·K) - 对流换热系数
T_amb = 298.0  # K - 环境温度
T0 = 298.0  # K - 初始温度
L_th = 0.0105  # m - 特征长度
k_ref = 2.1  # W/(m·K) - 参考热导率

@printf("  对流换热系数: h = %.1f W/(m²·K)\n", h)
@printf("  环境温度: T_amb = %.1f K\n", T_amb)
@printf("  初始温度: T0 = %.1f K\n", T0)
@printf("  特征长度: L_th = %.4f m\n", L_th)
@printf("  热导率: k_ref = %.2f W/(m·K)\n", k_ref)

# 计算 Biot 数
Bi = h * L_th / k_ref
@printf("\n[2] Biot 数分析:\n")
@printf("  Bi = h × L_th / k_ref = %.3f\n", Bi)

if Bi < 0.1
    println("  → Bi < 0.1: 对流换热弱，内部温差小")
elseif Bi < 1.0
    println("  → 0.1 < Bi < 1.0: 对流和传导都重要")
else
    println("  → Bi > 1.0: 对流换热强，外表面温度接近 T_amb")
end

# 估算温度响应时间
rho_c = 1.5e6  # J/(m³·K) - 体积热容（估算）
alpha = k_ref / rho_c  # m²/s - 热扩散率
tau_diffusion = L_th^2 / alpha  # s - 扩散时间尺度
tau_convection = rho_c * L_th / h  # s - 对流时间尺度

@printf("\n[3] 时间尺度分析:\n")
@printf("  热扩散率: α = %.2e m²/s\n", alpha)
@printf("  扩散时间尺度: τ_diff = L²/α = %.1f s\n", tau_diffusion)
@printf("  对流时间尺度: τ_conv = ρc×L/h = %.1f s\n", tau_convection)

if tau_convection < tau_diffusion / 10
    println("  ⚠️ 对流时间尺度远小于扩散时间尺度")
    println("     → 外表面温度会被快速"锁定"在 T_amb 附近")
else
    println("  ✓ 对流和扩散时间尺度相当")
end

# 估算稳态温升
Q_gen = 0.7  # W - 典型生成功率
V_cell = 2.34e-5  # m³ - 电池体积
R_out = 0.0105  # m - 外半径
A_surf = 2 * π * R_out * 0.07  # m² - 外表面积（估算）

q_vol = Q_gen / V_cell  # W/m³ - 体积热源密度

# 稳态能量平衡: Q_gen = h × A × (T_surf - T_amb)
delta_T_steady = Q_gen / (h * A_surf)

@printf("\n[4] 稳态温升估算:\n")
@printf("  生成功率: Q_gen = %.2f W\n", Q_gen)
@printf("  外表面积: A_surf ≈ %.4f m²\n", A_surf)
@printf("  体积热源: q_vol = %.1f W/m³\n", q_vol)
@printf("\n  稳态温升: ΔT = Q_gen / (h × A)\n")
@printf("           = %.2f / (%.1f × %.4f)\n", Q_gen, h, A_surf)
@printf("           = %.2f K\n", delta_T_steady)
@printf("\n  → 外表面温度 ≈ %.2f K\n", T_amb + delta_T_steady)

println("\n[5] 可能原因分析:")

println("\n  原因 1: 换热系数过大 ⭐⭐⭐")
println("    当前 h = 150 W/(m²·K) 是典型的强制对流水平")
println("    如果实际应该是自然对流（h = 5-25 W/(m²·K)），")
println("    那么当前设置会导致外表面温度过于接近 T_amb")
println()
@printf("    如果 h = 10 W/(m²·K): ΔT ≈ %.2f K\n", Q_gen / (10.0 * A_surf))
@printf("    如果 h = 150 W/(m²·K): ΔT ≈ %.2f K\n", delta_T_steady)

println("\n  原因 2: T_amb = T0 = 298K")
println("    环境温度恰好等于初始温度")
println("    如果对流很强，外圈温度会始终接近 298K")
println("    这在物理上是合理的（热平衡）")

println("\n  原因 3: 仿真时间太短")
t_sim = 360.0  # s - 仿真时间（testexample.jl）
@printf("    仿真时间: %d s\n", Int(t_sim))
@printf("    扩散时间: %.1f s\n", tau_diffusion)
@printf("    对流时间: %.1f s\n", tau_convection)
if t_sim < tau_diffusion
    println("    ⚠️ 仿真时间短于扩散时间，热量可能还未到达外圈")
else
    println("    ✓ 仿真时间足够长")
end

println("\n  原因 4: 热源分布不均")
println("    如果热源主要在内圈，外圈温升会很小")
println("    特别是在短时间仿真中")

println("\n  原因 5: 边界条件实施错误（需检查）")
println("    可能外圈节点被错误地强制固定在 T_amb")
println("    例如：极耳边界条件误应用到外圈")

println("\n[6] 验证方法:")
println("""
1. 检查外圈节点是否被固定:
   - 查看 KT 矩阵对角元素（外圈节点）
   - 如果有超大值（如 1e12），说明被惩罚法固定了

2. 修改对流换热系数:
   在测试脚本中添加：
   ```julia
   # 降低对流换热系数（模拟绝热）
   param_dim.cell.h = 10.0  # 或 0.0（完全绝热）
   ```

3. 延长仿真时间:
   ```julia
   opt.time = [0.0, 3600]  # 1小时
   ```

4. 输出外圈节点温度:
   在结果分析中：
   ```julia
   T_final = result["thermal2D T_nodes [K]"]
   r_nodes = sqrt.(mesh.node[:,1].^2 .+ mesh.node[:,2].^2)
   idx_outer = findall(r_nodes .> 0.95 * maximum(r_nodes))
   println("外圈节点温度:")
   @printf("  最小: %.2f K\\n", minimum(T_final[idx_outer]))
   @printf("  平均: %.2f K\\n", mean(T_final[idx_outer]))
   @printf("  最大: %.2f K\\n", maximum(T_final[idx_outer]))
   ```

5. 检查能量守恒:
   ```julia
   Q_gen = result["generation power"]
   Q_conv = result["convection power"]
   Q_storage = result["storage rate"]
   
   # 应该满足: Q_gen ≈ Q_conv + Q_storage
   ```
""")

println("\n[7] 推荐修改:")

println("""
方案 A: 调整换热系数（如果实际是绝热/自然对流）
  ```julia
  # 在 testexample.jl 或 testexample_simple_coupling.jl 中
  param_dim = JuBat.ChooseCell("Jellyroll")
  param_dim.cell.h = 10.0  # W/(m²·K) - 自然对流
  # 或
  param_dim.cell.h = 0.0  # 完全绝热
  ```

方案 B: 检查边界条件实施
  在 ThermalDistributed.jl 中添加调试输出：
  ```julia
  # 在 _apply_convection_bc! 函数中
  println("外边界节点数: ", sum(is_outer))
  println("对流换热系数: Bi = ", Bi)
  ```

方案 C: 延长仿真时间并观察
  ```julia
  opt.time = [0.0, 3600]  # 足够长以达到稳态
  ```
""")

println("\n[8] 预期结果:")

println("""
如果设置正确，外圈温度应该：
  ✓ 初始时刻: T = 298 K (T0)
  ✓ 短时间后: T = 298 + 几K (开始升温)
  ✓ 稳态时: T ≈ 298 + ΔT_steady (取决于 h)

如果 h = 150 (强对流):
  → ΔT_steady ≈ 1-2 K (外圈温度接近 T_amb)

如果 h = 10 (自然对流):
  → ΔT_steady ≈ 10-20 K (外圈温度明显升高)

如果 h = 0 (绝热):
  → ΔT 持续增加（无散热）
""")

println("\n" * "="^70)
println("诊断完成")
println("="^70)

println("""
最可能的原因：
  ⭐⭐⭐ 对流换热系数设置为强制对流（h=150），
       导致外表面温度被"锁定"在 T_amb = 298K 附近

建议：
  1. 确认实际换热条件（强制/自然对流/绝热）
  2. 相应调整 param_dim.cell.h
  3. 如果是绝热条件，设置 h = 0.0
""")
