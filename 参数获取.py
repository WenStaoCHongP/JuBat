#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
使用 PyBaMM 与 Jellyroll 参数函数，导出用于交叉验证的关键数据：
- OCV_n(x)、OCV_p(x)、OCV_total(x)（无量纲与物理电压）
- 电解液 κ(c, T)、D_e(c, T)
- 参考标度 phi = R*T_ref/F

输出 CSV 到 JuBat/output/ 下，供与 Julia 侧 PostProcessing 曲线对比。

注意：这里直接复刻 Jellyroll.jl 中的函数形式，确保一致性；
如你切换了参数文件，请相应修改下方函数/常数。
"""

import os
import math
from typing import Iterable

# 可选：使用 pybamm 的常数，若不可用则退化到内置常数
try:
	import pybamm
	R = float(pybamm.constants.R.value)  # 8.314...
	F = float(pybamm.constants.F.value)  # 96485...
except Exception:
	R = 8.314
	F = 96485.33289

# ---------------- Jellyroll 参数（与 src/parameters/Jellyroll.jl 对齐） ----------------
T_ref = 298.0  # K
phi = R * T_ref / F

# 电压上下限（物理）
V_L = 2.5
V_H = 4.3

# 正极/负极 cs_max 与初始 cs0（用于给出一个 OCV@初始点做 sanity check）
PE_cs_max = 63104.0
PE_cs0 = 17038.0
NE_cs_max = 33133.0
NE_cs0 = 29866.0

def x_clip(x):
	# 与 Julia 侧保持一致，避免在 x->0 或 x->1 的奇点
	eps = 1e-8
	return max(eps, min(1.0 - eps, x))

# OCV_p(x) 与 OCV_n(x) – 物理电压（V），复制自 Jellyroll.jl（新版 tanh 形式）
def U_p_V(x: float) -> float:
	x = x_clip(x)
	# -0.8090*x + 4.4875 - 0.0428*tanh(18.5138*(x-0.5542))
	# - 17.7326*tanh(15.7890*(x-0.3117)) + 17.5842*tanh(15.9308*(x-0.3120))
	return (
		-0.8090 * x
		+ 4.4875
		- 0.0428 * math.tanh(18.5138 * (x - 0.5542))
		- 17.7326 * math.tanh(15.7890 * (x - 0.3117))
		+ 17.5842 * math.tanh(15.9308 * (x - 0.3120))
	)

def U_n_V(x: float) -> float:
	x = x_clip(x)
	# 1.97938*exp(-39.3631*x) + 0.2482 - 0.0909*tanh(29.8538*(x-0.1234))
	# - 0.04478*tanh(14.9159*(x-0.2769)) - 0.0205*tanh(30.4444*(x-0.6103))
	return (
		1.97938 * math.exp(-39.3631 * x)
		+ 0.2482
		- 0.0909 * math.tanh(29.8538 * (x - 0.1234))
		- 0.04478 * math.tanh(14.9159 * (x - 0.2769))
		- 0.0205 * math.tanh(30.4444 * (x - 0.6103))
	)

# dUdT 这里在 Jellyroll.jl 被设为 0
def dUdT_p(x: float) -> float:
	return 0.0

def dUdT_n(x: float) -> float:
	return 0.0

# 电解液传输性质（物理单位）：κ 与 D_e（复制新版 Jellyroll.jl 表达式）
def kappa_SI(c_e: float, T: float = 298.0) -> float:
	# c_e [mol/m^3]; 原式用 (x/1000) 即单位转成 mol/L
	y = c_e / 1000.0
	return 0.1297 * y**3 - 2.51 * y ** 1.5 + 3.329 * y

def De_SI(c_e: float, T: float = 298.0) -> float:
	y = c_e / 1000.0
	return 8.794e-11 * y**2 - 3.972e-10 * y + 4.862e-10


# ---------------- 导出 CSV ----------------
def linspace(a: float, b: float, n: int) -> Iterable[float]:
	if n <= 1:
		yield a
		return
	step = (b - a) / (n - 1)
	for i in range(n):
		yield a + i * step


def ensure_dir(path: str):
	os.makedirs(path, exist_ok=True)


def export_ocv_curves(out_dir: str, npts: int = 2001):
	fp = os.path.join(out_dir, "ocv_curves.csv")
	with open(fp, "w", encoding="utf-8") as f:
		f.write(
			"x,Up_V,Un_V,OCV_V,Up_nd,Un_nd,OCV_nd\n"
		)
		for x in linspace(0.001, 0.999, npts):
			up = U_p_V(x)
			un = U_n_V(x)
			ocv = up - un
			f.write(
				f"{x:.8f},{up:.8f},{un:.8f},{ocv:.8f},{up/phi:.8f},{un/phi:.8f},{ocv/phi:.8f}\n"
			)
	# 初始点诊断
	x_p0 = PE_cs0 / PE_cs_max
	x_n0 = NE_cs0 / NE_cs_max
	up0 = U_p_V(x_p0)
	un0 = U_n_V(x_n0)
	ocv0 = up0 - un0
	print(
		f"[ocv] phi={phi:.6f} V; x_p0={x_p0:.5f}, x_n0={x_n0:.5f}; OCV0={ocv0:.6f} V (nd={ocv0/phi:.6f})"
	)
	print(f"[ocv] bounds: [{V_L}, {V_H}] V -> nd=[{V_L/phi:.3f}, {V_H/phi:.3f}]")
	print(f"[ocv] CSV written: {fp}")


def export_electrolyte_props(out_dir: str, c_min=100.0, c_max=2000.0, npts: int = 100):
	fp = os.path.join(out_dir, "electrolyte_props.csv")
	with open(fp, "w", encoding="utf-8") as f:
		f.write("c_e_molm3,T_K,kappa_S_per_m,De_m2_per_s\n")
		for T in (273.15, 298.0, 323.15):
			for c in linspace(c_min, c_max, npts):
				f.write(f"{c:.3f},{T:.2f},{kappa_SI(c,T):.6e},{De_SI(c,T):.6e}\n")
	print(f"[elec] CSV written: {fp}")


def main():
	repo_root = os.path.dirname(os.path.abspath(__file__))
	out_dir = os.path.join(repo_root, "JuBat", "output")
	ensure_dir(out_dir)
	print(f"phi = {phi:.6f} V (R={R}, F={F}, T_ref={T_ref})")
	print(f"输出目录: {out_dir}")
	export_ocv_curves(out_dir)
	export_electrolyte_props(out_dir)


if __name__ == "__main__":
	main()

