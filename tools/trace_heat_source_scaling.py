#!/usr/bin/env python3
"""
追踪热源的归一化流程
"""

print("=" * 80)
print("热源归一化流程追踪")
print("=" * 80)

print("""
[步骤1] _compute_layer_heat_sources 计算热源
────────────────────────────────────────────────────────────────────────────
文件: ThermalDistributed.jl 第525行

Q_NE = as_n * abs(j_n) * abs(eta_n) + ...

其中:
- as_n = 归一化的比表面积 (无量纲, as/a0)
- j_n = 归一化的电流密度 (无量纲, j/j0)  
- eta_n = 归一化的过电位 (无量纲, eta/phi)

乘积 Q_NE 是无量纲的，但它的尺度是什么？

有量纲形式:
    q = as * j * eta
      = [1/m] × [A/m²] × [V]
      = [W/m³]

归一化:
    q* = (as/a0) * (j/j0) * (eta/phi)
    
恢复有量纲:
    q = q* × (a0 × j0 × phi)
      = q* × (a0 × I_typ/(a0×L×A) × phi)
      = q* × (I_typ × phi / (L×A))
      = q* × (I_typ × phi / V_cell)

结论: Q_layers 的尺度是 I_typ × phi / V_cell ← 这是电功率尺度!
""")

print("""
[步骤2] CallModel_SimpleCoupling 使用热源
────────────────────────────────────────────────────────────────────────────
文件: CallModel_SimpleCoupling.jl 第93行

total_heat_W = q_ref * sum(q_elem .* elem_volumes)

其中:
- q_elem = 无量纲热源密度 (来自 Q_layers)
- elem_volumes = 有量纲体积 [m³]
- q_ref = case.param_dim.scale.q_th

计算:
    total_heat_W = q_ref × Σ(q_elem × V_elem)
                 = [W/m³] × [无量纲] × [m³]
                 = [W]

问题: q_ref 应该与 q_elem 的尺度匹配!

如果 q_elem 的尺度是 I×phi/V，那么 q_ref 也应该是 I×phi/V。
""")

print("""
[步骤3] _assemble_force_vector 构建载荷向量
────────────────────────────────────────────────────────────────────────────
文件: ThermalDistributed.jl 第317行

coeff_f = q_gs .* (wJ ./ L_th^2)

其中:
- q_gs = 高斯点的无量纲热源
- wJ = 雅可比 × 权重 [m²]
- L_th = 特征长度 [m]

单位分析:
    coeff_f = [无量纲] × [m²] / [m²]
            = [无量纲]

这是正确的，因为载荷向量 F_T 应该是无量纲的。

但关键是: q_gs 是相对于什么尺度无量纲化的?
""")

print("""
[问题诊断]
════════════════════════════════════════════════════════════════════════════

情况A: 如果 q_ref 使用傅里叶尺度 (k×T/L²)
────────────────────────────────────────────────────────────────────────────
q_ref = 5.68×10⁶ W/m³

但 q_elem 的实际尺度是 I×phi/V = 5483 W/m³

尺度不匹配! 
    total_heat_W = 5.68×10⁶ × q_elem × V
                 = 1000 × (5483 × q_elem × V)
                 = 1000 × Q_actual

结果: 热源偏大1000倍! ✗


情况B: 如果 q_ref 使用电功率尺度 (I×phi/V)
────────────────────────────────────────────────────────────────────────────
q_ref = 5483 W/m³

q_elem 的实际尺度也是 5483 W/m³

尺度匹配! 
    total_heat_W = 5483 × q_elem × V
                 = Q_actual

结果: 热源正常! ✓


但温度呢?
════════════════════════════════════════════════════════════════════════════

在热传导方程中:
    ρc ∂T/∂t = ∇·(k∇T) + q

无量纲化后需要平衡:
    (ρc/t_th) × ∂T*/∂t* = (k/L²) × ∇*²T* + (q_ref/T_ref) × q*

标准平衡条件:
    ρc/t_th = k/L² = q_ref/T_ref

但是! 在耦合系统中，通过 t_ratio 调整后:
    (ρc/t_th) × (t0/t_th) = ρc/t0

实际平衡条件变成:
    ρc/t0 = q_ref/T_ref

检查:
    ρc/t0 = 1.0×10⁶ / 3600 = 278 [J/(m³·K·s)] = [W/(m³·K)]
    q_ref/T_ref = 5483 / 298 = 18.4 [W/(m³·K)]

不匹配! 差了 278/18.4 ≈ 15 倍!

这就是为什么温度有问题的原因!
""")

print("""
[根本问题]
════════════════════════════════════════════════════════════════════════════

在电化学模型中，热源是相对于 I×phi 归一化的（总功率）。
在热模型中，体积热源密度应该相对于 q_th 归一化。

但 q_th 的定义方式导致了不一致:
- 电化学热源尺度: I×phi/V (基于电功率)
- 热传导方程尺度: k×T/L² (基于傅里叶传导)

两者不匹配!

解决方案1: 在热源计算中转换尺度
────────────────────────────────────────────────────────────────────────────
在 _compute_layer_heat_sources 返回后，需要转换:

    Q_layers_fourier = Q_layers × (I×phi/V) / (k×T/L²)
                     = Q_layers × scale_factor

其中 scale_factor ≈ 0.001


解决方案2: 使用一致的尺度系统
────────────────────────────────────────────────────────────────────────────
重新设计无量纲化，使电化学和热模型使用相同的尺度基准。

这需要大量重构，不太现实。


解决方案3: 在简化耦合模式中特殊处理 (推荐!)
────────────────────────────────────────────────────────────────────────────
在 CallModel_SimpleCoupling 中:

1. 识别 Q_layers 使用的是电化学尺度
2. 计算 scale_conversion = (I×phi/V) / (k×T/L²)
3. 转换: Q_layers_thermal = Q_layers × scale_conversion
4. 然后使用傅里叶的 q_th 继续计算
""")

import math

# 参数
I_typ = 5
phi = 0.0257
V = 2.34e-5
k = 2.1
T_ref = 298
L_th = 0.0105

q_electrochemical = I_typ * phi / V
q_fourier = k * T_ref / L_th**2
scale_factor = q_electrochemical / q_fourier

print("\n[数值计算]")
print("=" * 80)
print(f"电化学热源尺度: I×φ/V = {q_electrochemical:.2e} W/m³")
print(f"傅里叶热源尺度: k×T/L² = {q_fourier:.2e} W/m³")
print(f"尺度转换因子: {scale_factor:.6e}")
print(f"即: Q_thermal = Q_electrochemical × {scale_factor:.6e}")
