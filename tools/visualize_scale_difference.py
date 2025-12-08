#!/usr/bin/env python3
"""
可视化两种热模型的尺度差异
"""

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================================
# 参数设置
# ============================================================================

# LG M50 参数
I_typ_lg = 5.0  # A
phi_lg = 0.0257  # V (R×T/F)
V_lg = 2.42e-5  # m³
area_lg = 0.1027  # m²

# Jellyroll 参数
I_typ_jr = 5.0  # A
phi_jr = 0.0257  # V
V_jr = 2.28e-5  # m³
area_jr = 0.147  # m²
Rout_jr = 0.0105  # m
Rin_jr = 1.92e-3  # m
L_chem = 1.728e-4  # m (电化学层厚度)
k_th = 2.1  # W/(m·K)
T_ref = 298  # K

# ============================================================================
# 计算尺度参数
# ============================================================================

# 集总模型尺度
scale_lumped = I_typ_lg * phi_lg
print(f"集总模型尺度: {scale_lumped:.4f} W")

# 分布式模型尺度 (当前实现：用Rout)
q_th_current = k_th * T_ref / Rout_jr**2
print(f"分布式模型尺度 (用Rout): {q_th_current:.2e} W/m³")

# 分布式模型尺度 (建议：用电化学长度)
q_th_suggested = k_th * T_ref / L_chem**2
print(f"分布式模型尺度 (用L_chem): {q_th_suggested:.2e} W/m³")

# 分布式模型尺度 (另一建议：用电功率/体积)
q_th_power_based = I_typ_jr * phi_jr / V_jr
print(f"分布式模型尺度 (用I×φ/V): {q_th_power_based:.2e} W/m³")

# ============================================================================
# 数值示例
# ============================================================================

# 假设无量纲热源为 6.0
Q_nd = 6.0

# 集总模型
Q_lumped = Q_nd * scale_lumped
print(f"\n集总模型总热源: {Q_lumped:.3f} W")

# 简化耦合模型 (当前实现)
Q_distributed_current = Q_nd * q_th_current * V_jr
print(f"简化耦合总热源 (当前): {Q_distributed_current:.1f} W")

# 简化耦合模型 (建议实现 - 用L_chem)
Q_distributed_suggested = Q_nd * q_th_suggested * V_jr
print(f"简化耦合总热源 (建议L_chem): {Q_distributed_suggested:.1e} W")

# 简化耦合模型 (建议实现 - 用电功率基准)
Q_distributed_power = Q_nd * q_th_power_based * V_jr
print(f"简化耦合总热源 (建议电功率): {Q_distributed_power:.3f} W")

# ============================================================================
# 可视化
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# -------------------- 图1: 尺度因子对比 --------------------
ax1 = axes[0, 0]
scales = [
    ('集总模型\n(I×φ)', scale_lumped, 'W'),
    ('分布式(当前)\n(k×T/R²)', q_th_current, 'W/m³'),
    ('分布式(建议1)\n(k×T/L²)', q_th_suggested, 'W/m³'),
    ('分布式(建议2)\n(I×φ/V)', q_th_power_based, 'W/m³')
]

labels = [s[0] for s in scales]
values = [s[1] for s in scales]
units = [s[2] for s in scales]

# 使用对数坐标
ax1.bar(range(len(labels)), values, color=['blue', 'red', 'orange', 'green'], alpha=0.7)
ax1.set_yscale('log')
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, fontsize=9)
ax1.set_ylabel('数值 (对数坐标)', fontsize=10)
ax1.set_title('尺度因子对比', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

# 添加数值标签
for i, (v, u) in enumerate(zip(values, units)):
    ax1.text(i, v, f'{v:.2e}\n{u}', ha='center', va='bottom', fontsize=8)

# -------------------- 图2: 总热源对比 --------------------
ax2 = axes[0, 1]
Q_values = [
    ('集总模型', Q_lumped, 'blue'),
    ('简化耦合\n(当前)', Q_distributed_current, 'red'),
    ('简化耦合\n(建议1)', Q_distributed_suggested, 'orange'),
    ('简化耦合\n(建议2)', Q_distributed_power, 'green')
]

Q_labels = [q[0] for q in Q_values]
Q_vals = [q[1] for q in Q_values]
Q_colors = [q[2] for q in Q_values]

bars = ax2.bar(range(len(Q_labels)), Q_vals, color=Q_colors, alpha=0.7)
ax2.set_yscale('log')
ax2.set_xticks(range(len(Q_labels)))
ax2.set_xticklabels(Q_labels, fontsize=9)
ax2.set_ylabel('总热源 (W, 对数坐标)', fontsize=10)
ax2.set_title(f'总热源对比 (无量纲值={Q_nd})', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')
ax2.axhline(y=0.77, color='black', linestyle='--', linewidth=1, label='实测值 (~0.77 W)')
ax2.legend(fontsize=9)

# 添加数值标签
for i, v in enumerate(Q_vals):
    if v < 1:
        label = f'{v:.3f} W'
    elif v < 1000:
        label = f'{v:.1f} W'
    else:
        label = f'{v:.2e} W'
    ax2.text(i, v, label, ha='center', va='bottom', fontsize=8)

# -------------------- 图3: L_th 对 q_th 的影响 --------------------
ax3 = axes[1, 0]
L_th_range = np.logspace(-5, -1, 100)  # 10^-5 到 10^-1 m
q_th_range = k_th * T_ref / L_th_range**2

ax3.plot(L_th_range, q_th_range, linewidth=2, color='purple')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('热特征长度 L_th (m)', fontsize=10)
ax3.set_ylabel('q_th (W/m³)', fontsize=10)
ax3.set_title('q_th 与 L_th 的关系 (q ∝ 1/L²)', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, which='both')

# 标记关键点
ax3.scatter([Rout_jr], [q_th_current], s=100, color='red', marker='o', 
           label=f'当前 (Rout={Rout_jr*1000:.1f}mm)')
ax3.scatter([L_chem], [q_th_suggested], s=100, color='orange', marker='s', 
           label=f'建议1 (L_chem={L_chem*1e6:.0f}μm)')
ax3.legend(fontsize=9)

# -------------------- 图4: 合理性分析 --------------------
ax4 = axes[1, 1]

# 典型值范围
Q_typical_min = 0.5  # W
Q_typical_max = 2.0  # W
q_typical_min = Q_typical_min / V_jr  # W/m³
q_typical_max = Q_typical_max / V_jr  # W/m³

# 比较各种方案得到的热源
scenarios = [
    '集总模型\n(实测)', 
    '简化耦合\n(当前实现)',
    '简化耦合\n(建议1)',
    '简化耦合\n(建议2)'
]
Q_scenarios = [Q_lumped, Q_distributed_current, Q_distributed_suggested, Q_distributed_power]
colors_scenarios = ['blue', 'red', 'orange', 'green']
markers = ['✓' if Q_typical_min <= q <= Q_typical_max else '✗' for q in Q_scenarios]

# 绘制合理范围
ax4.axhspan(Q_typical_min, Q_typical_max, alpha=0.2, color='green', label='合理范围')

# 绘制各方案的值
for i, (s, Q, c, m) in enumerate(zip(scenarios, Q_scenarios, colors_scenarios, markers)):
    # 限制显示范围
    Q_plot = min(Q, 100)  # 超过100的截断显示
    marker_style = 'o' if m == '✓' else 'x'
    marker_size = 150 if m == '✓' else 200
    ax4.scatter(i, Q_plot, s=marker_size, color=c, marker=marker_style, alpha=0.7)
    
    # 添加标签
    if Q <= 100:
        label = f'{Q:.2f} W\n{m}'
    else:
        label = f'{Q:.2e} W\n{m}'
    ax4.text(i, Q_plot, label, ha='center', va='bottom', fontsize=8)

ax4.set_xticks(range(len(scenarios)))
ax4.set_xticklabels(scenarios, fontsize=9)
ax4.set_ylabel('总热源 (W)', fontsize=10)
ax4.set_title('合理性检查', fontsize=12, fontweight='bold')
ax4.set_ylim(0, 10)
ax4.grid(True, alpha=0.3, axis='y')
ax4.legend(fontsize=9, loc='upper right')

# 添加说明文本
textstr = '✓ = 在合理范围内\n✗ = 超出合理范围'
ax4.text(0.02, 0.98, textstr, transform=ax4.transAxes, fontsize=9,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# ============================================================================
# 保存图像
# ============================================================================

plt.tight_layout()
plt.savefig('/workspace/docs/内热源尺度差异分析.png', dpi=150, bbox_inches='tight')
print("\n图像已保存到: /workspace/docs/内热源尺度差异分析.png")

# ============================================================================
# 生成文本总结
# ============================================================================

summary = f"""
{'='*80}
内热源数量级差异分析总结
{'='*80}

1. 尺度因子对比:
   集总模型 (I×φ):               {scale_lumped:.4f} W
   分布式模型 (当前, k×T/R²):     {q_th_current:.2e} W/m³
   分布式模型 (建议1, k×T/L²):    {q_th_suggested:.2e} W/m³
   分布式模型 (建议2, I×φ/V):     {q_th_power_based:.2e} W/m³

2. 总热源对比 (假设无量纲值 = {Q_nd}):
   集总模型:                      {Q_lumped:.3f} W      ← 合理 ✓
   简化耦合 (当前实现):           {Q_distributed_current:.1f} W      ← 偏大1000倍 ✗
   简化耦合 (建议1, 用L_chem):    {Q_distributed_suggested:.2e} W  ← 偏大更多 ✗
   简化耦合 (建议2, 用I×φ/V):     {Q_distributed_power:.3f} W      ← 合理 ✓

3. 关键发现:
   - 当前实现用 Rout 作为热特征长度，导致 q_th 过大
   - 用电化学长度会更糟糕 (因为 L_chem << Rout)
   - 建议使用基于电功率的尺度: q_th = I×φ/V

4. 根本问题:
   q_th = k×T/L² 这个公式来自傅里叶热传导的无量纲化
   但它与电化学产热的尺度不匹配

5. 推荐解决方案:
   修改 SetParams.jl 中 q_th 的计算:
   
   # 当前实现
   q_th = k_th × T_ref / L_th²
   
   # 建议改为
   q_th = I_typ × phi / V_cell
   
   这样可以确保两种模型的尺度一致。

{'='*80}
"""

print(summary)

# 保存到文件
with open('/workspace/docs/内热源尺度差异分析_数值总结.txt', 'w', encoding='utf-8') as f:
    f.write(summary)

print("数值总结已保存到: /workspace/docs/内热源尺度差异分析_数值总结.txt")
