"""
分析Jellyroll参数导致电压降过小和温升不明显的原因

目标：计算关键尺度因子，分析热源计算公式
"""

using Printf

println("="^80)
println("Jellyroll参数分析")
println("="^80)

# 从Jellyroll.jl参数文件提取关键值
println("\n[1] 几何参数")
Rout = 0.065/2  # 外半径
Rin = 1.92e-3   # 内半径
width = 5.212683e-2  # 高度
cell_length = 938.7693138  # 展开长度
area = width * cell_length  # 有效电化学面积

@printf("  外半径 Rout: %.4f m\n", Rout)
@printf("  内半径 Rin: %.4f m\n", Rin)
@printf("  高度 width: %.4f m\n", width)
@printf("  展开长度: %.2f m\n", cell_length)
@printf("  电化学面积 area: %.4f m²\n", area)

# 电极参数
println("\n[2] 电极几何参数")
PE_thickness = 7.56e-5   # 正极厚度
NE_thickness = 8.52e-5   # 负极厚度
SP_thickness = 1.2e-5    # 隔膜厚度
PCC_thickness = 10e-6    # 正极集流体厚度
NCC_thickness = 10e-6    # 负极集流体厚度

L = PE_thickness + NE_thickness + SP_thickness  # 特征电化学长度

@printf("  正极厚度: %.2e m\n", PE_thickness)
@printf("  负极厚度: %.2e m\n", NE_thickness)
@printf("  隔膜厚度: %.2e m\n", SP_thickness)
@printf("  总厚度 L: %.2e m\n", L)

# 材料参数
println("\n[3] 材料参数")
PE_lambda = 1.04  # 正极热导率
NE_lambda = 1.58  # 负极热导率
SP_lambda = 0.16  # 隔膜热导率
PCC_lambda = 237.  # 正极集流体热导率
NCC_lambda = 401.  # 负极集流体热导率

PE_rho = 2895  # 正极密度
NE_rho = 1555  # 负极密度
SP_rho = 1100  # 隔膜密度
PCC_rho = 2700  # 正极集流体密度
NCC_rho = 8940  # 负极集流体密度

PE_cp = 1.27e3  # 正极比热容
NE_cp = 1.437e3  # 负极比热容
SP_cp = 1.978e3  # 隔膜比热容
PCC_cp = 9.03e2  # 正极集流体比热容
NCC_cp = 3.85e2  # 负极集流体比热容

@printf("  正极热导率: %.2f W/m/K\n", PE_lambda)
@printf("  负极热导率: %.2f W/m/K\n", NE_lambda)

# 计算有效热导率（径向和切向）
lambda_r = (NE_thickness + SP_thickness + PE_thickness + PCC_thickness + NCC_thickness) /
           (NE_thickness/NE_lambda + SP_thickness/SP_lambda + PE_thickness/PE_lambda + 
            PCC_thickness/PCC_lambda + NCC_thickness/NCC_lambda)

lambda_t = (NE_thickness*NE_lambda + SP_thickness*SP_lambda + PE_thickness*PE_lambda + 
            PCC_thickness*PCC_lambda + NCC_thickness*NCC_lambda) /
           (NE_thickness + SP_thickness + PE_thickness + PCC_thickness + NCC_thickness)

@printf("  径向有效热导率 λ_r: %.4f W/m/K\n", lambda_r)
@printf("  切向有效热导率 λ_t: %.4f W/m/K\n", lambda_t)

# 计算整体热容（用于热模型）
volume = pi * (Rout^2 - Rin^2) * width
mass = 8.521789028721335e-01

# 计算有效比热容
cell_heat_Q = (PE_cp * PE_rho * PE_thickness + NE_cp * NE_rho * NE_thickness +
               SP_cp * SP_rho * SP_thickness + NCC_cp * NCC_rho * NCC_thickness +
               PCC_cp * PCC_rho * PCC_thickness) /
              (PE_rho * PE_thickness + NE_rho * NE_thickness + SP_rho * SP_thickness + 
               NCC_rho * NCC_thickness + PCC_rho * PCC_thickness)

cell_rho = mass / volume

@printf("  电池质量: %.4f kg\n", mass)
@printf("  电池体积: %.6e m³\n", volume)
@printf("  电池密度: %.2f kg/m³\n", cell_rho)
@printf("  有效比热容: %.2f J/kg/K\n", cell_heat_Q)

# 电化学参数
println("\n[4] 电化学参数")
I1C = 60  # 1C电流
capacity = 60  # 容量

PE_eps_s = 1 - 0.335 - 0.025  # 正极固相体积分数
NE_eps_s = 1 - 0.25 - 0.0326  # 负极固相体积分数
PE_rs = 5.22e-6  # 正极颗粒半径
NE_rs = 5.86e-6  # 负极颗粒半径
PE_as = 3 * PE_eps_s / PE_rs  # 正极比表面积
NE_as = 3 * NE_eps_s / NE_rs  # 负极比表面积

PE_sig = 0.18  # 正极电导率
NE_sig = 100   # 负极电导率

@printf("  1C电流 I1C: %.1f A\n", I1C)
@printf("  容量: %.1f Ah\n", capacity)
@printf("  正极比表面积 as_p: %.2e 1/m\n", PE_as)
@printf("  负极比表面积 as_n: %.2e 1/m\n", NE_as)
@printf("  正极电导率 σ_p: %.2f S/m\n", PE_sig)
@printf("  负极电导率 σ_n: %.1f S/m\n", NE_sig)

# 尺度因子
println("\n[5] 尺度因子计算")
T_ref = 298  # 参考温度
F = 96485.33289  # 法拉第常数
R = 8.314  # 气体常数

I_typ = I1C
phi_scale = T_ref * R / F

@printf("  典型电流 I_typ: %.1f A\n", I_typ)
@printf("  电势尺度 φ: %.6f V\n", phi_scale)

# 热尺度因子
L_th = Rout  # 特征热长度（使用外半径）
k_ref = PE_lambda  # 参考热导率
rho_c_ref = cell_rho * cell_heat_Q  # 体积热容
q_ref = k_ref * T_ref / L_th^2  # 参考体积热源
t_th = rho_c_ref * L_th^2 / k_ref  # 热扩散时间尺度

@printf("  特征热长度 L_th: %.4f m\n", L_th)
@printf("  参考热导率 k_ref: %.2f W/m/K\n", k_ref)
@printf("  体积热容 ρc: %.2e J/m³/K\n", rho_c_ref)
@printf("  参考体积热源 q_ref: %.2e W/m³\n", q_ref)
@printf("  热扩散时间尺度 t_th: %.2f s\n", t_th)

# 关键发现1：电化学尺度与热尺度的匹配
println("\n[6] 尺度匹配分析")
q_typical_chem = I_typ * phi_scale / volume  # 典型电化学功率密度

@printf("  典型电化学功率: %.2f W\n", I_typ * phi_scale)
@printf("  典型电化学功率密度: %.2e W/m³\n", q_typical_chem)
@printf("  热源尺度 q_ref: %.2e W/m³\n", q_ref)
@printf("  比值 q_chem/q_ref: %.4f\n", q_typical_chem / q_ref)

# 关键发现2：电流密度和归一化
println("\n[7] 电流密度分析")
j_typ = I_typ / area  # 典型电流密度（A/m²）
i_vol = I_typ / volume  # 典型体积电流（A/m³）

@printf("  典型面电流密度 j: %.2f A/m²\n", j_typ)
@printf("  典型体积电流: %.2e A/m³\n", i_vol)

# 关键发现3：热源公式检查
println("\n[8] 热源公式单位检查")
println("  多SPMe模式中的热源计算（第510-529行）：")

# 假设典型值
I_e_nd = 1.0  # 无量纲单元电流
eta_nd = 0.05  # 无量纲过电位（典型值~0.05V/0.025V）
j_nd = 1.0  # 无量纲界面电流密度
T_nd = 1.0  # 无量纲温度
as_nd = NE_as * 1e-6  # 归一化比表面积（as * r0）
sig_nd = NE_sig * L / (I_typ / area) / phi_scale  # 归一化电导率

println("\n  反应热公式: Q_rxn = as * |j| * |η|")
@printf("    无量纲值: %.2e * %.2e * %.2e = %.2e\n", as_nd, j_nd, eta_nd, as_nd * j_nd * eta_nd)
println("    单位检查: [1/m] * [A/m²] * [V] = [W/m³] ✓")

println("\n  欧姆热公式: Q_ohm = I²/(3σ)")
@printf("    无量纲值: %.2e² / (3 * %.2e) = %.2e\n", I_e_nd, sig_nd, I_e_nd^2 / (3 * sig_nd))
println("    单位检查: [A]² / [S/m] = [A² m / (A/V)] = [V·A/m] = [W/m³] ✗")
println("    ❌ 问题：公式应该是 I²/(σ·A²·L) 或类似形式，而不是 I²/(3σ)")

# 关键发现4：面积缺失问题
println("\n[9] 关键问题诊断")
println("  ❌ 问题1：热源计算缺少面积因子")
println("     热源公式中使用了电流I，但应该使用电流密度j = I/A")
println("     导致热源被高估了 A² 倍（约 $(area^2)）")
println()
println("  ❌ 问题2：归一化不一致")
println("     电流I使用 I_typ=$(I_typ)A 归一化")
println("     但热源使用 q_ref=$(q_ref) W/m³ 归一化")
println("     两者之间的转换需要 I_typ*phi/volume = $(I_typ*phi_scale/volume) W/m³")
println("     而 q_ref = $(q_ref) W/m³")
println("     比值: $(I_typ*phi_scale/volume/q_ref)")
println()
println("  ⚠️  问题3：dUdT = 0")
println("     Jellyroll参数中 PE.dUdT = x-> 0*x 和 NE.dUdT = x-> 0*x")
println("     这导致可逆热源为0，但这不是主要问题")

# 关键发现5：电压降过小的原因
println("\n[10] 电压降过小的原因")
println("  正常情况（testexample_simple_coupling）:")
println("    单SPMe模型，电流全部集中在一个电化学单元")
println("    内阻正常，电压降正常（~1.6V）")
println()
println("  多SPMe情况（testexample）:")
println("    如果电流分配不当，每个单元分到的电流过小")
println("    例如：总电流5A分配到80个单元，每单元仅0.0625A")
println("    但这不应该导致总电压降过小，除非...")
println()
println("  ❌ 可能原因：电流分配算法问题")
println("     如果分流求解器认为所有单元并联，")
println("     那么每个单元看到相同电压，但分担不同电流")
println("     总电压降 = 单元电压降（而不是累加）")
println("     这会导致电压降比串联情况小很多！")

println("\n" * "="^80)
println("结论")
println("="^80)
println()
println("1. 电压降过小（~0.001V）的主要原因：")
println("   - 多SPMe模式将80个热单元视为并联")
println("   - 每个单元独立计算SPMe，共享相同电压")
println("   - 总电压降 = 单元电压降（~0.001V），而不是串联累加")
println("   - 这与物理实际不符：电池应该是单层结构，不是80层并联")
println()
println("2. 温升不明显的主要原因：")
println("   - 热源计算公式中电流-热源转换可能有单位问题")
println("   - 欧姆热公式 Q = I²/(3σ) 缺少几何因子")
println("   - 应该是 Q = I² * R_ohm / Volume，其中 R_ohm = L/(σ·A)")
println()
println("3. 建议修改：")
println("   - 检查热源计算公式（CallModel_MultiSPMe第510-529行）")
println("   - 添加正确的面积/体积因子")
println("   - 或者：考虑使用简化耦合模式（simple_thermal_coupling）")
println("   - 简化模式假设电流均匀，使用单个SPMe，更符合实际")
println()
println("="^80)
