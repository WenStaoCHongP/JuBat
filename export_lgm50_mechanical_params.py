#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 PyBaMM 导出 LG M50 (Chen2020) 的力学参数

包括：
1. 偏摩尔体积 (Partial Molar Volume) -> 用于计算化学膨胀系数
2. 弹性模量、泊松比等力学性质
3. 热膨胀系数
4. 计算化学膨胀系数 beta_c

理论背景：
化学膨胀系数 beta_c 与偏摩尔体积 Ω (Omega) 的关系：
    beta_c ≈ 3 * Ω / V_m
其中 V_m 是摩尔体积

或者更直接地：
    ε_vol = Ω * Δc_s / (1 + Ω * c_s)
简化为线性关系（小应变）：
    ε_vol ≈ Ω * Δc_s
    beta_c = Ω / c_s_max

依赖：pybamm>=23.4
"""

import os
import csv
import json
from typing import Dict, Any

try:
    import pybamm
    HAS_PYBAMM = True
except ImportError:
    pybamm = None
    HAS_PYBAMM = False


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def get_parameter_safely(param, key: str, default=None):
    """安全地从参数集中获取参数值"""
    try:
        val = param[key]
        return val if val is not None else default
    except Exception:
        return default


def export_mechanical_parameters():
    """导出力学相关参数"""
    
    if not HAS_PYBAMM:
        print("❌ 错误: PyBaMM 未安装")
        print("请安装: pip install pybamm")
        return
    
    print("="*70)
    print("从 PyBaMM 导出 LG M50 (Chen2020) 力学参数")
    print("="*70)
    
    # 加载 Chen2020 参数集
    try:
        param = pybamm.ParameterValues("Chen2020")
        print("✓ 成功加载 Chen2020 参数集")
    except Exception as e:
        print(f"❌ 加载参数集失败: {e}")
        return
    
    # 输出目录
    repo_root = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(repo_root, "JuBat", "output")
    ensure_dir(out_dir)
    
    # ====================================================================
    # 1. 偏摩尔体积 (Partial Molar Volume)
    # ====================================================================
    
    print("\n" + "="*70)
    print("1. 偏摩尔体积 (Partial Molar Volume)")
    print("="*70)
    
    # PyBaMM 中的参数名称
    pmv_keys = [
        "Lithium metal partial molar volume [m3.mol-1]",
        "Positive electrode partial molar volume [m3.mol-1]",
        "Negative electrode partial molar volume [m3.mol-1]",
    ]
    
    pmv_data = {}
    for key in pmv_keys:
        val = get_parameter_safely(param, key)
        if val is not None:
            pmv_data[key] = float(val)
            print(f"  {key}: {val:.6e} m³/mol")
    
    # ====================================================================
    # 2. 最大浓度 (用于计算 beta_c)
    # ====================================================================
    
    print("\n" + "="*70)
    print("2. 最大锂浓度")
    print("="*70)
    
    cs_max_n = get_parameter_safely(
        param, 
        "Maximum concentration in negative electrode [mol.m-3]",
        33133.0
    )
    cs_max_p = get_parameter_safely(
        param,
        "Maximum concentration in positive electrode [mol.m-3]",
        63104.0
    )
    
    print(f"  负极最大浓度: {cs_max_n:.1f} mol/m³")
    print(f"  正极最大浓度: {cs_max_p:.1f} mol/m³")
    
    # ====================================================================
    # 3. 计算化学膨胀系数 beta_c
    # ====================================================================
    
    print("\n" + "="*70)
    print("3. 计算化学膨胀系数 beta_c")
    print("="*70)
    
    # 方法1: 从偏摩尔体积计算
    # beta_c = Omega / V_m ≈ Omega * c_s_max
    # 其中 V_m 是单位体积，取 1 m³ 时，beta_c ≈ Omega * c_s_max
    
    # 获取 Omega（PyBaMM 中通常存储为 partial molar volume）
    Omega_n = get_parameter_safely(
        param,
        "Negative electrode partial molar volume [m3.mol-1]",
        3.1e-6  # 石墨的典型值
    )
    
    Omega_p = get_parameter_safely(
        param,
        "Positive electrode partial molar volume [m3.mol-1]",
        -7.28e-7  # NMC 的典型值（负值表示收缩）
    )
    
    print(f"\n偏摩尔体积:")
    print(f"  负极 Ω_n: {Omega_n:.6e} m³/mol")
    print(f"  正极 Ω_p: {Omega_p:.6e} m³/mol")
    
    # 计算 beta_c
    # 公式: ε_vol = Ω * Δc_s
    # 归一化: ε_vol / (c_s_max * ΔSOC) = Ω * c_s_max
    # 因此: beta_c = Ω * c_s_max
    
    beta_c_n = abs(Omega_n * cs_max_n)  # 取绝对值，因为我们关心体积变化大小
    beta_c_p = abs(Omega_p * cs_max_p)
    
    print(f"\n化学膨胀系数 beta_c:")
    print(f"  负极: {beta_c_n:.6f} [-]")
    print(f"  正极: {beta_c_p:.6f} [-]")
    print(f"\n解释: 这表示当 SOC 从 0→1 变化时，体积应变约为 {beta_c_n*100:.2f}% (负极) 和 {beta_c_p*100:.2f}% (正极)")
    
    # 与文献值对比
    print(f"\n文献典型值对比:")
    print(f"  石墨负极: 1-4% (文献)  vs  {beta_c_n*100:.2f}% (计算)")
    print(f"  NMC正极:  1-2% (文献)  vs  {beta_c_p*100:.2f}% (计算)")
    
    # ====================================================================
    # 4. 其他力学参数（如果 PyBaMM 有的话）
    # ====================================================================
    
    print("\n" + "="*70)
    print("4. 其他力学参数")
    print("="*70)
    
    # 尝试获取弹性模量、泊松比等
    mechanical_keys = [
        "Negative electrode Young's modulus [Pa]",
        "Positive electrode Young's modulus [Pa]",
        "Negative electrode Poisson's ratio",
        "Positive electrode Poisson's ratio",
        "Negative electrode thermal expansion coefficient [K-1]",
        "Positive electrode thermal expansion coefficient [K-1]",
    ]
    
    mechanical_data = {}
    found_any = False
    for key in mechanical_keys:
        val = get_parameter_safely(param, key)
        if val is not None:
            mechanical_data[key] = float(val)
            print(f"  {key}: {val:.6e}")
            found_any = True
    
    if not found_any:
        print("  ⚠️  PyBaMM Chen2020 参数集中未找到显式的力学参数")
        print("  建议使用文献值:")
        print("    负极 (石墨):")
        print("      E = 10-15 GPa")
        print("      ν = 0.3")
        print("      α = 1-2×10⁻⁵ K⁻¹")
        print("    正极 (NMC811):")
        print("      E = 100-200 GPa")
        print("      ν = 0.2-0.3")
        print("      α = 1-1.5×10⁻⁵ K⁻¹")
    
    # ====================================================================
    # 5. 导出为 CSV 和 JSON
    # ====================================================================
    
    print("\n" + "="*70)
    print("5. 导出结果")
    print("="*70)
    
    # CSV 格式
    csv_file = os.path.join(out_dir, "lgm50_mechanical_params.csv")
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value", "unit", "description"])
        
        # 偏摩尔体积
        writer.writerow(["Omega_n", f"{Omega_n:.6e}", "m³/mol", "负极偏摩尔体积"])
        writer.writerow(["Omega_p", f"{Omega_p:.6e}", "m³/mol", "正极偏摩尔体积"])
        
        # 最大浓度
        writer.writerow(["cs_max_n", f"{cs_max_n:.1f}", "mol/m³", "负极最大锂浓度"])
        writer.writerow(["cs_max_p", f"{cs_max_p:.1f}", "mol/m³", "正极最大锂浓度"])
        
        # 化学膨胀系数
        writer.writerow(["beta_c_n", f"{beta_c_n:.6f}", "-", "负极化学膨胀系数"])
        writer.writerow(["beta_c_p", f"{beta_c_p:.6f}", "-", "正极化学膨胀系数"])
        writer.writerow(["volume_change_n", f"{beta_c_n*100:.2f}", "%", "负极体积变化 (0→100% SOC)"])
        writer.writerow(["volume_change_p", f"{beta_c_p*100:.2f}", "%", "正极体积变化 (0→100% SOC)"])
        
        # 其他力学参数（如果有）
        for key, val in mechanical_data.items():
            writer.writerow([key, f"{val:.6e}", "", ""])
    
    print(f"✓ CSV 导出至: {csv_file}")
    
    # JSON 格式
    json_file = os.path.join(out_dir, "lgm50_mechanical_params.json")
    output_data = {
        "source": "PyBaMM Chen2020 parameter set",
        "partial_molar_volumes": {
            "negative": {
                "value": Omega_n,
                "unit": "m³/mol"
            },
            "positive": {
                "value": Omega_p,
                "unit": "m³/mol"
            }
        },
        "max_concentrations": {
            "negative": {
                "value": cs_max_n,
                "unit": "mol/m³"
            },
            "positive": {
                "value": cs_max_p,
                "unit": "mol/m³"
            }
        },
        "chemical_expansion_coefficients": {
            "negative": {
                "beta_c": beta_c_n,
                "volume_change_percent": beta_c_n * 100,
                "description": "体积应变 / SOC变化"
            },
            "positive": {
                "beta_c": beta_c_p,
                "volume_change_percent": beta_c_p * 100,
                "description": "体积应变 / SOC变化"
            }
        },
        "other_mechanical_properties": mechanical_data,
        "notes": {
            "calculation": "beta_c = |Omega * cs_max|",
            "interpretation": "beta_c 表示 SOC 从 0→1 时的体积应变",
            "literature_comparison": {
                "graphite": "1-4% (typical)",
                "NMC": "1-2% (typical)"
            }
        }
    }
    
    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    print(f"✓ JSON 导出至: {json_file}")
    
    # ====================================================================
    # 6. 生成 Julia 代码片段
    # ====================================================================
    
    print("\n" + "="*70)
    print("6. Julia 代码片段")
    print("="*70)
    
    julia_code = f"""
# LG M50 (Chen2020) 力学参数
# 来源: PyBaMM Chen2020 参数集

# 负极 (石墨)
NE.Omega = {Omega_n:.6e}     # 偏摩尔体积 [m³/mol]
NE.cs_max = {cs_max_n:.1f}    # 最大锂浓度 [mol/m³]
NE.beta_c = {beta_c_n:.6f}    # 化学膨胀系数 [-]
NE.E = 15e9                   # 弹性模量 [Pa] (文献值)
NE.nu = 0.3                   # 泊松比 [-] (文献值)
NE.alphaT = 1.5e-5            # 热膨胀系数 [1/K] (文献值)

# 正极 (NMC811)
PE.Omega = {Omega_p:.6e}     # 偏摩尔体积 [m³/mol]
PE.cs_max = {cs_max_p:.1f}    # 最大锂浓度 [mol/m³]
PE.beta_c = {beta_c_p:.6f}    # 化学膨胀系数 [-]
PE.E = 150e9                  # 弹性模量 [Pa] (文献值)
PE.nu = 0.3                   # 泊松比 [-] (文献值)
PE.alphaT = 1.0e-5            # 热膨胀系数 [1/K] (文献值)

# 说明:
# - beta_c 从 PyBaMM 的偏摩尔体积计算得到
# - E, nu, alphaT 使用文献典型值（PyBaMM 中未提供）
# - 负极体积变化约 {beta_c_n*100:.1f}%，正极约 {beta_c_p*100:.1f}%
"""
    
    julia_file = os.path.join(out_dir, "lgm50_mechanical_params.jl")
    with open(julia_file, "w", encoding="utf-8") as f:
        f.write(julia_code.strip())
    
    print(f"✓ Julia 代码导出至: {julia_file}")
    print("\n可直接复制到 src/parameters/LGM50.jl 中使用")
    
    print("\n" + "="*70)
    print("完成! ✅")
    print("="*70)
    
    # ====================================================================
    # 7. 检查所有可用的参数（调试用）
    # ====================================================================
    
    print("\n提示: 如需查看 PyBaMM 中所有可用的参数，运行:")
    print("  python -c \"import pybamm; p=pybamm.ParameterValues('Chen2020'); print('\\n'.join(sorted(p.keys())))\"")


def list_all_mechanical_keys():
    """列出 PyBaMM 中所有可能与力学相关的参数键"""
    
    if not HAS_PYBAMM:
        print("PyBaMM 未安装")
        return
    
    param = pybamm.ParameterValues("Chen2020")
    
    keywords = [
        "molar", "volume", "expansion", "modulus", "poisson", 
        "thermal", "mechanical", "stress", "strain", "elastic"
    ]
    
    print("\nPyBaMM Chen2020 中可能的力学相关参数:")
    print("="*70)
    
    for key in sorted(param.keys()):
        if any(kw in key.lower() for kw in keywords):
            try:
                val = param[key]
                print(f"{key}: {val}")
            except Exception:
                print(f"{key}: <无法获取>")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == "--list":
        list_all_mechanical_keys()
    else:
        export_mechanical_parameters()
