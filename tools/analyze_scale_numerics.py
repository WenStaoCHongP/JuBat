#!/usr/bin/env python3
"""
数值分析：内热源尺度差异
"""

import math

print("=" * 80)
print("内热源数量级差异分析 - 数值计算")
print("=" * 80)

# ============================================================================
# 参数设置
# ============================================================================

print("\n[1] 参数设置")
print("-" * 80)

# LG M50 参数
I_typ_lg = 5.0  # A
phi_lg = 0.0257  # V (R×T/F)
V_lg = 2.42e-5  # m³
area_lg = 0.1027  # m²

print("LG M50 电池:")
print(f"  I_typ = {I_typ_lg} A")
print(f"  phi = {phi_lg} V")
print(f"  体积 = {V_lg:.2e} m³")
print(f"  面积 = {area_lg:.4f} m²")

# Jellyroll 参数
I_typ_jr = 5.0  # A
phi_jr = 0.0257  # V
Rout_jr = 0.0105  # m
Rin_jr = 1.92e-3  # m
width_jr = 0.07  # m (高度)
V_jr = math.pi * (Rout_jr**2 - Rin_jr**2) * width_jr
area_jr = 0.147  # m²
L_chem = 1.728e-4  # m (电化学层厚度)
k_th = 2.1  # W/(m·K)
T_ref = 298  # K

print("\nJellyroll 电池:")
print(f"  I_typ = {I_typ_jr} A")
print(f"  phi = {phi_jr} V")
print(f"  Rout = {Rout_jr*1000:.2f} mm")
print(f"  Rin = {Rin_jr*1000:.2f} mm")
print(f"  高度 = {width_jr*1000:.1f} mm")
print(f"  体积 = {V_jr:.2e} m³")
print(f"  面积 = {area_jr:.4f} m²")
print(f"  电化学层厚度 = {L_chem*1e6:.1f} μm")

# ============================================================================
# 尺度因子计算
# ============================================================================

print("\n[2] 尺度因子计算")
print("-" * 80)

# 集总模型尺度
scale_lumped = I_typ_lg * phi_lg
print(f"\n集总模型:")
print(f"  尺度 = I_typ × phi")
print(f"       = {I_typ_lg} × {phi_lg}")
print(f"       = {scale_lumped:.6f} W")

# 分布式模型尺度 - 方案1：当前实现 (用Rout)
q_th_current = k_th * T_ref / Rout_jr**2
print(f"\n分布式模型 - 方案1 (当前实现):")
print(f"  尺度 = k_th × T_ref / L_th²")
print(f"  L_th = Rout = {Rout_jr*1000:.2f} mm")
print(f"  q_th = {k_th} × {T_ref} / ({Rout_jr})²")
print(f"       = {q_th_current:.3e} W/m³")

# 分布式模型尺度 - 方案2：用电化学长度
q_th_chem = k_th * T_ref / L_chem**2
print(f"\n分布式模型 - 方案2 (用电化学长度):")
print(f"  L_th = L_chem = {L_chem*1e6:.1f} μm")
print(f"  q_th = {k_th} × {T_ref} / ({L_chem})²")
print(f"       = {q_th_chem:.3e} W/m³")

# 分布式模型尺度 - 方案3：用电功率基准
q_th_power = I_typ_jr * phi_jr / V_jr
print(f"\n分布式模型 - 方案3 (基于电功率):")
print(f"  尺度 = I_typ × phi / V_cell")
print(f"       = {I_typ_jr} × {phi_jr} / {V_jr:.3e}")
print(f"       = {q_th_power:.3e} W/m³")

# ============================================================================
# 尺度比较
# ============================================================================

print("\n[3] 尺度因子比较")
print("-" * 80)

ratio1 = q_th_current * V_jr / scale_lumped
ratio2 = q_th_chem * V_jr / scale_lumped
ratio3 = q_th_power * V_jr / scale_lumped

print(f"\n将分布式尺度转换为等效功率尺度:")
print(f"  方案1: q_th × V = {q_th_current:.3e} × {V_jr:.3e}")
print(f"                  = {q_th_current * V_jr:.3f} W")
print(f"         与集总模型比值 = {ratio1:.1f}")
print(f"\n  方案2: q_th × V = {q_th_chem:.3e} × {V_jr:.3e}")
print(f"                  = {q_th_chem * V_jr:.3e} W")
print(f"         与集总模型比值 = {ratio2:.2e}")
print(f"\n  方案3: q_th × V = {q_th_power:.3e} × {V_jr:.3e}")
print(f"                  = {q_th_power * V_jr:.3f} W")
print(f"         与集总模型比值 = {ratio3:.1f}")

# ============================================================================
# 总热源计算示例
# ============================================================================

print("\n[4] 总热源计算示例")
print("-" * 80)

# 假设无量纲热源为 6.0 (从日志文件观察到的典型值)
Q_nd = 6.0

print(f"\n假设无量纲热源值 = {Q_nd}")
print(f"(这是从 thermal_example 运行日志中观察到的典型值)")

# 集总模型
Q_lumped = Q_nd * scale_lumped
print(f"\n集总模型:")
print(f"  Q_total = Q_nd × (I_typ × phi)")
print(f"          = {Q_nd} × {scale_lumped:.6f}")
print(f"          = {Q_lumped:.4f} W")
print(f"  状态: 合理 ✓ (与实测的 0.73-0.80 W 一致)")

# 简化耦合模型 - 方案1
Q_dist_current = Q_nd * q_th_current * V_jr
print(f"\n简化耦合模型 - 方案1 (当前实现):")
print(f"  Q_total = Q_nd × q_th × V")
print(f"          = {Q_nd} × {q_th_current:.3e} × {V_jr:.3e}")
print(f"          = {Q_dist_current:.2f} W")
print(f"  状态: 偏大约 {Q_dist_current/Q_lumped:.0f} 倍 ✗")

# 简化耦合模型 - 方案2
Q_dist_chem = Q_nd * q_th_chem * V_jr
print(f"\n简化耦合模型 - 方案2 (用电化学长度):")
print(f"  Q_total = {Q_nd} × {q_th_chem:.3e} × {V_jr:.3e}")
print(f"          = {Q_dist_chem:.3e} W")
print(f"  状态: 偏大约 {Q_dist_chem/Q_lumped:.2e} 倍 ✗✗")

# 简化耦合模型 - 方案3
Q_dist_power = Q_nd * q_th_power * V_jr
print(f"\n简化耦合模型 - 方案3 (基于电功率):")
print(f"  Q_total = {Q_nd} × {q_th_power:.3e} × {V_jr:.3e}")
print(f"          = {Q_dist_power:.4f} W")
print(f"  状态: 与集总模型比值 {Q_dist_power/Q_lumped:.2f} ≈ 1 ✓")

# ============================================================================
# 合理性分析
# ============================================================================

print("\n[5] 合理性分析")
print("-" * 80)

print("\n典型锂离子电池的热源范围:")
Q_typical_min = 0.5  # W
Q_typical_max = 2.0  # W
print(f"  - 1C 放电: {Q_typical_min} - {Q_typical_max} W")
print(f"  - 平均体积热源密度: {Q_typical_min/V_jr:.2e} - {Q_typical_max/V_jr:.2e} W/m³")

print("\n各方案的状态:")
scenarios = [
    ("集总模型 (实测)", Q_lumped),
    ("简化耦合 - 方案1 (当前)", Q_dist_current),
    ("简化耦合 - 方案2 (用L_chem)", Q_dist_chem),
    ("简化耦合 - 方案3 (用I×φ/V)", Q_dist_power),
]

for name, Q in scenarios:
    in_range = Q_typical_min <= Q <= Q_typical_max
    status = "✓ 合理" if in_range else "✗ 不合理"
    print(f"  {name:30s}: {Q:10.3e} W  {status}")

# ============================================================================
# 根本原因分析
# ============================================================================

print("\n[6] 根本原因分析")
print("-" * 80)

print("\nq_th = k×T/L² 公式的来源:")
print("  这是傅里叶热传导方程无量纲化的标准尺度")
print("  适用于纯热传导问题")
print("")
print("但在电池产热问题中:")
print("  - 热源来自电化学反应和欧姆热")
print("  - 其尺度应与电功率 I×V 相关")
print("  - 而不是傅里叶传导尺度 k×T/L²")

print("\n数值证明:")
print(f"  电功率尺度:      I×φ = {scale_lumped:.6f} W")
print(f"  傅里叶尺度×体积: k×T/L²×V = {q_th_current * V_jr:.2f} W")
print(f"  相差约 {q_th_current * V_jr / scale_lumped:.0f} 倍")

print("\n为什么会这样?")
print("  因为 L_th (热特征长度) 和电化学产热无关")
print("  L_th 反映的是热传导的空间尺度")
print("  而产热强度取决于电流密度和电阻")

# ============================================================================
# L_th 敏感性分析
# ============================================================================

print("\n[7] L_th 敏感性分析")
print("-" * 80)

print("\nq_th 与 L_th 的关系: q_th ∝ 1/L_th²")
print("L_th 的不同选择:")

L_th_options = [
    ("Rout (当前)", Rout_jr),
    ("Rout - Rin (径向)", Rout_jr - Rin_jr),
    ("电化学层厚度", L_chem),
    ("高度", width_jr),
]

print(f"\n{'选择':<20s} {'L_th (m)':<15s} {'q_th (W/m³)':<20s} {'q_th×V (W)':<15s}")
print("-" * 70)

for name, L in L_th_options:
    q = k_th * T_ref / L**2
    Q_equiv = q * V_jr
    print(f"{name:<20s} {L:<15.3e} {q:<20.3e} {Q_equiv:<15.3f}")

print("\n观察:")
print("  - L_th 减小 → q_th 增大 (平方反比)")
print("  - 用电化学层厚度最糟糕 (太小)")
print("  - 用 Rout 已经导致 q_th 过大")
print("  - 用高度 (70mm) 会更合理，但仍不如用电功率基准")

# ============================================================================
# 推荐方案
# ============================================================================

print("\n[8] 推荐方案")
print("=" * 80)

print("\n问题症结:")
print("  集总模型和分布式模型使用了不同的无量纲化哲学")
print("  - 集总: 基于电功率 (I×V)")
print("  - 分布式: 基于傅里叶传导 (k×T/L²)")

print("\n短期解决方案:")
print("  修改 SetParams.jl 中 q_th 的计算:")
print("")
print("  # 当前实现 (第262行)")
print("  q_ref = k_ref * param_dim.scale.T_ref / L_th^2")
print("")
print("  # 建议改为")
print("  q_ref = param_dim.scale.I_typ * param_dim.scale.phi / param_dim.cell.volume")
print("")
print("  这样可以确保:")
print("  1. 与集总模型的尺度一致")
print("  2. 总热源数值合理")
print("  3. 两种模型可以直接比较")

print("\n长期解决方案:")
print("  统一整个代码库的无量纲化方案")
print("  建议使用电功率相关的尺度，因为:")
print("  1. 物理意义清晰 (热源 = 功耗)")
print("  2. 与电化学模型自然衔接")
print("  3. 数值稳定性好")

# ============================================================================
# 总结
# ============================================================================

print("\n[9] 总结")
print("=" * 80)

print("\n造成内热源数量级差异的原因:")
print("  1. 不同的无量纲化方式")
print("     - 集总模型: I×φ ≈ 0.13 W")
print("     - 分布式模型: k×T/L² ≈ 5.7×10⁶ W/m³")
print("")
print("  2. 尺度差异巨大")
print(f"     - (k×T/L²) × V / (I×φ) ≈ {q_th_current * V_jr / scale_lumped:.0f}")
print("")
print("  3. 物理意义不同")
print("     - 电功率尺度 vs 热传导尺度")
print("")
print("  4. L_th 选择不当")
print("     - 用 Rout 导致 q_th 过大")
print("     - L_th ∝ q_th⁻¹/²，影响平方级")

print("\n验证:")
print(f"  - 集总模型实测: 0.73-0.80 W ✓")
print(f"  - 理论计算: {Q_lumped:.3f} W ✓")
print(f"  - 分布式(当前): {Q_dist_current:.1f} W (偏大{Q_dist_current/Q_lumped:.0f}倍) ✗")
print(f"  - 分布式(建议): {Q_dist_power:.3f} W (合理) ✓")

print("\n" + "=" * 80)

# 保存结果
output_file = "/workspace/docs/内热源尺度差异_数值分析结果.txt"
with open(output_file, "w", encoding="utf-8") as f:
    # 这里可以把所有输出重定向到文件，但为了简单起见暂不实现
    pass

print(f"\n分析完成！")
print(f"建议查看:")
print(f"  - /workspace/docs/内热源差异分析_简要版.md")
print(f"  - /workspace/docs/内热源数量级差异分析.md")
