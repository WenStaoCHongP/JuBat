#!/usr/bin/env python3
"""
验证尺度转换的正确性
"""

# 参数值（从jellyroll参数）
I_typ = 5.0  # A
phi = 0.0257  # V
V_cell = 2.34e-5  # m³
k_th = 2.1  # W/(m·K)
T_ref = 298.0  # K
L_th = 0.0105  # m

# 计算两个尺度
q_ec_scale = I_typ * phi / V_cell
q_fourier_scale = k_th * T_ref / L_th**2

print("=" * 60)
print("尺度参数")
print("=" * 60)
print(f"电化学尺度: q_ec = I×φ/V = {q_ec_scale:.2f} W/m³")
print(f"傅里叶尺度: q_fourier = k×T/L² = {q_fourier_scale:.2e} W/m³")
print(f"尺度比: q_ec / q_fourier = {q_ec_scale / q_fourier_scale:.6f}")
print()

# 假设一个物理热源值
Q_phys = 30000.0  # W/m³（典型值）

print("=" * 60)
print("无量纲化示例")
print("=" * 60)
print(f"物理热源: Q_phys = {Q_phys:.0f} W/m³")
print()

# 电化学尺度归一化
q_elem_ec = Q_phys / q_ec_scale
print(f"电化学无量纲: q_elem = Q_phys / q_ec")
print(f"                     = {Q_phys:.0f} / {q_ec_scale:.2f}")
print(f"                     = {q_elem_ec:.4f}")
print()

# 傅里叶尺度归一化  
q_elem_fourier = Q_phys / q_fourier_scale
print(f"傅里叶无量纲: q_elem = Q_phys / q_fourier")
print(f"                     = {Q_phys:.0f} / {q_fourier_scale:.2e}")
print(f"                     = {q_elem_fourier:.6f}")
print()

# 转换验证
ratio = q_ec_scale / q_fourier_scale
q_converted = q_elem_ec * ratio
print("=" * 60)
print("转换验证")
print("=" * 60)
print(f"转换公式: q_fourier_nd = q_ec_nd × (q_ec / q_fourier)")
print(f"                       = {q_elem_ec:.4f} × {ratio:.6f}")
print(f"                       = {q_converted:.6f}")
print()
print(f"直接计算的傅里叶无量纲: {q_elem_fourier:.6f}")
print(f"误差: {abs(q_converted - q_elem_fourier):.2e}")
print()

# 恢复物理值验证
Q_recovered_ec = q_elem_ec * q_ec_scale
Q_recovered_fourier = q_elem_fourier * q_fourier_scale
print("=" * 60)
print("物理值恢复验证")
print("=" * 60)
print(f"从电化学无量纲恢复: {q_elem_ec:.4f} × {q_ec_scale:.2f}")
print(f"                    = {Q_recovered_ec:.0f} W/m³")
print()
print(f"从傅里叶无量纲恢复: {q_elem_fourier:.6f} × {q_fourier_scale:.2e}")
print(f"                    = {Q_recovered_fourier:.0f} W/m³")
print()

# 总功率计算
print("=" * 60)
print("总功率计算")
print("=" * 60)
Q_total = Q_phys * V_cell
print(f"物理热源密度: {Q_phys:.0f} W/m³")
print(f"电池体积: {V_cell:.2e} m³")
print(f"总功率: {Q_total:.3f} W")
print()

# 通过傅里叶无量纲计算总功率
Q_total_from_fourier = q_elem_fourier * q_fourier_scale * V_cell
print(f"从傅里叶无量纲计算: {q_elem_fourier:.6f} × {q_fourier_scale:.2e} × {V_cell:.2e}")
print(f"                    = {Q_total_from_fourier:.3f} W")
print()

# 关键结论
print("=" * 60)
print("关键结论")
print("=" * 60)
print("✓ 转换公式正确: q_fourier_nd = q_ec_nd × (q_ec / q_fourier)")
print(f"✓ 虽然 q_ec / q_fourier ≈ {ratio:.3f} < 1，")
print("  这是正常的，因为傅里叶尺度比电化学尺度大约1000倍")
print("✓ 物理值在两种尺度下保持一致")
print()
print("预期效果:")
print("- 如果原始代码 q_elem = 6.0（电化学尺度）")
print(f"- 转换后 heat_fields = 6.0 × {ratio:.6f} = {6.0 * ratio:.6f}（傅里叶尺度）")
print(f"- 物理热源 = 6.0 × {q_ec_scale:.2f} = {6.0 * q_ec_scale:.0f} W/m³")
print(f"- 总功率 = {6.0 * q_ec_scale * V_cell:.3f} W")
