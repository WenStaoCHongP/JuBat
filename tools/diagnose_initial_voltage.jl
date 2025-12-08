"""
诊断初始电压差异

分析 testexample.jl (多SPMe) 和 testexample_simple_coupling.jl (简化耦合)
的初始电压为什么不同。

已知：
- 多SPMe: 初始电压 3.9992 V
- 简化耦合: 初始电压 4.0637 V
- 差距: 0.064 V
"""

using Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

println("="^70)
println("初始电压差异诊断")
println("="^70)

# ============================================================================
# 1. 加载 Jellyroll 电池参数
# ============================================================================
println("\n[1] 加载 Jellyroll 电池参数...")

param_dim = JuBat.ChooseCell("Jellyroll")

println("  ✓ 电池参数加载完成")
println("\n  负极参数:")
@printf("    cs0 = %.2f mol/m³\n", param_dim.NE.cs0)
@printf("    cs_max = %.2f mol/m³\n", param_dim.NE.cs_max)
@printf("    theta_init = cs0/cs_max = %.4f\n", param_dim.NE.cs0 / param_dim.NE.cs_max)
@printf("    theta_100 = %.4f (SOC=100%%)\n", param_dim.NE.theta_100)
@printf("    theta_0 = %.4f (SOC=0%%)\n", param_dim.NE.theta_0)

println("\n  正极参数:")
@printf("    cs0 = %.2f mol/m³\n", param_dim.PE.cs0)
@printf("    cs_max = %.2f mol/m³\n", param_dim.PE.cs_max)
@printf("    theta_init = cs0/cs_max = %.4f\n", param_dim.PE.cs0 / param_dim.PE.cs_max)
@printf("    theta_100 = %.4f (SOC=100%%)\n", param_dim.PE.theta_100)
@printf("    theta_0 = %.4f (SOC=0%%)\n", param_dim.PE.theta_0)

# ============================================================================
# 2. 计算初始 SOC
# ============================================================================
println("\n[2] 估算初始 SOC...")

# 负极：SOC越高，theta越接近theta_100
theta_n = param_dim.NE.cs0 / param_dim.NE.cs_max
soc_n = (theta_n - param_dim.NE.theta_0) / (param_dim.NE.theta_100 - param_dim.NE.theta_0)

# 正极：SOC越高，theta越接近theta_0（锂转移到负极）
theta_p = param_dim.PE.cs0 / param_dim.PE.cs_max
soc_p = (param_dim.PE.theta_0 - theta_p) / (param_dim.PE.theta_0 - param_dim.PE.theta_100)

@printf("  负极 SOC (基于theta): %.2f%%\n", soc_n * 100)
@printf("  正极 SOC (基于theta): %.2f%%\n", soc_p * 100)

# 平衡 SOC（理论上应该一致）
soc_avg = (soc_n + soc_p) / 2
@printf("  平均 SOC: %.2f%%\n", soc_avg * 100)

# ============================================================================
# 3. 计算初始 OCV
# ============================================================================
println("\n[3] 计算初始开路电压 (OCV)...")

# 使用 theta 计算 OCV
U_p = param_dim.PE.U(theta_p)
U_n = param_dim.NE.U(theta_n)
OCV_init = U_p - U_n

@printf("  正极电势 U_p(%.4f) = %.4f V\n", theta_p, U_p)
@printf("  负极电势 U_n(%.4f) = %.4f V\n", theta_n, U_n)
@printf("  初始 OCV = U_p - U_n = %.4f V\n", OCV_init)

# ============================================================================
# 4. 分析可能的差异原因
# ============================================================================
println("\n[4] 分析可能的初始电压差异原因...")

println("\n  可能原因 1: 初始 SOC 实际不同")
println("    如果多SPMe使用不同的初始SOC（例如均匀分布但略低）")
println("    即使参数文件相同，初始化代码可能有差异")

println("\n  可能原因 2: 初始浓度分布不同")
println("    - 简化耦合: 单个SPMe，颗粒内均匀初始化")
println("    - 多SPMe: 每个单元独立SPMe，可能有不同的初始化逻辑")

println("\n  可能原因 3: 电解质浓度初始化")
println("    - 电解质浓度分布影响离子电导，从而影响欧姆压降")
@printf("    - 初始电解质浓度: %.2f mol/m³\n", param_dim.EL.ce0)

println("\n  可能原因 4: 初始化时的瞬态过程")
println("    - 某个模式可能在初始化时已经施加电流")
println("    - 导致初始电压包含过电位和欧姆压降")

# ============================================================================
# 5. 估算可能的 SOC 差异
# ============================================================================
println("\n[5] 反推 SOC 差异...")

V1 = 3.9992  # 多SPMe
V2 = 4.0637  # 简化耦合
delta_V = V2 - V1

@printf("  多SPMe初始电压: %.4f V\n", V1)
@printf("  简化耦合初始电压: %.4f V\n", V2)
@printf("  电压差: %.4f V\n", delta_V)

# 估算对应的 SOC 差异（假设 dV/dSOC ≈ 0.5 V 在这个SOC范围）
# 典型锂电池在高SOC区域，dV/dSOC ≈ 0.3-0.8 V
dV_dSOC_typical = 0.5  # V
delta_SOC_estimate = delta_V / dV_dSOC_typical

@printf("\n  如果是 SOC 差异造成 (假设 dV/dSOC ≈ %.1f V):\n", dV_dSOC_typical)
@printf("    估算 SOC 差异: %.2f%%\n", delta_SOC_estimate * 100)

# ============================================================================
# 6. 测试不同 SOC 的 OCV
# ============================================================================
println("\n[6] 测试不同 SOC 的 OCV 曲线...")

println("\n  SOC -> OCV 关系:")
for soc_test in [0.85, 0.90, 0.95, 0.98, 1.00]
    # 计算对应的 theta
    theta_n_test = param_dim.NE.theta_0 + soc_test * (param_dim.NE.theta_100 - param_dim.NE.theta_0)
    theta_p_test = param_dim.PE.theta_0 - soc_test * (param_dim.PE.theta_0 - param_dim.PE.theta_100)
    
    # 计算 OCV
    U_p_test = param_dim.PE.U(theta_p_test)
    U_n_test = param_dim.NE.U(theta_n_test)
    OCV_test = U_p_test - U_n_test
    
    @printf("    SOC = %.0f%% → OCV = %.4f V\n", soc_test * 100, OCV_test)
end

# ============================================================================
# 7. 建议的检查步骤
# ============================================================================
println("\n" * "="^70)
println("建议的检查步骤")
println("="^70)

println("""
1. 检查初始化代码:
   - src/Initialisation.jl (标准初始化)
   - src/ModelInitialisation_SimpleCoupling.jl (简化耦合)
   - 查看是否使用相同的 cs0 值

2. 检查初始 SOC 设置:
   - testexample.jl 中是否有 initial_soc_distribution 参数
   - testexample_simple_coupling.jl 的初始 SOC

3. 运行诊断测试:
   a) 在两个脚本中添加以下代码（初始化后）:
      ```julia
      println("初始状态诊断:")
      y0 = result["initial state"]
      Nrn = case.mesh["negative particle"].nlen
      Nrp = case.mesh["positive particle"].nlen
      
      csn_init = y0[1:Nrn]
      csp_init = y0[(Nrn+1):(Nrn+Nrp)]
      
      theta_n_init = mean(csn_init) * case.param_dim.NE.cs_max / case.param_dim.NE.cs_max
      theta_p_init = mean(csp_init) * case.param_dim.PE.cs_max / case.param_dim.PE.cs_max
      
      @printf("  负极 theta: %.4f\\n", theta_n_init)
      @printf("  正极 theta: %.4f\\n", theta_p_init)
      ```

4. 比较两种模式的初始状态:
   - 颗粒浓度分布
   - 电解质浓度
   - 温度场

5. 可能需要修正的地方:
   - 如果是初始化代码不一致，应该统一
   - 如果是有意设计（不同的初始SOC），应该在文档中说明
""")

println("="^70)
