#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
导出 PyBaMM 中 LGM50T（Chen2020）用于 SPMe 的关键电化学参数。

输出到 JuBat/output/ 下，包含：
- 标量参数（厚度、孔隙率、半径、最大浓度、传递数等） -> lgm50t_spme_scalars.csv
- OCP 曲线（Un/Up/OCV vs x） -> lgm50t_ocp_curves.csv
- 电解液性质（κ、D_e vs c,T） -> lgm50t_electrolyte_props.csv
- 固相扩散系数（D_s^n, D_s^p vs x,T） -> lgm50t_solid_diffusivity.csv

实现说明：
- 优先直接从 PyBaMM 的参数集中取值/调用；
- 若函数型参数在当前 PyBaMM 版本不可直接调用，则回退到常见的 Chen2020 显式近似式（与工程中广泛使用的拟合一致）。

依赖：pybamm>=23（建议最新稳定版）。
"""

import os
import csv
import math
from typing import Iterable, Tuple, Callable, Optional

# 兼容无 pybamm 环境：允许脚本在未安装时给出友好提示
try:
    import pybamm  # type: ignore
    HAS_PYBAMM = True
except Exception:
    pybamm = None  # type: ignore
    HAS_PYBAMM = False

# ---------------- 通用工具 ----------------

def linspace(a: float, b: float, n: int):
    if n <= 1:
        yield a
        return
    step = (b - a) / (n - 1)
    for i in range(n):
        yield a + i * step


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


# ---------------- Chen2020（LGM50T）常见显式近似（回退用） ----------------
# 注：这些式子与 PyBaMM Chen2020 参数集等价/等效；在某些版本中直接从参数对象提取函数
# 可能较繁琐，因此作为回退保持脚本可运行性。

# OCP_p(x) 与 OCP_n(x)
# 参考：常见的 Chen2020 NMC811/Graphite OCP 拟合
_DEF_EPS = 1e-8

def _x_clip(x: float) -> float:
    return max(_DEF_EPS, min(1.0 - _DEF_EPS, x))


def U_p_V_fallback(x: float) -> float:
    x = _x_clip(x)
    return (
        -0.8090 * x
        + 4.4875
        - 0.0428 * math.tanh(18.5138 * (x - 0.5542))
        - 17.7326 * math.tanh(15.7890 * (x - 0.3117))
        + 17.5842 * math.tanh(15.9308 * (x - 0.3120))
    )


def U_n_V_fallback(x: float) -> float:
    x = _x_clip(x)
    return (
        1.97938 * math.exp(-39.3631 * x)
        + 0.2482
        - 0.0909 * math.tanh(29.8538 * (x - 0.1234))
        - 0.04478 * math.tanh(14.9159 * (x - 0.2769))
        - 0.0205 * math.tanh(30.4444 * (x - 0.6103))
    )


def kappa_SI_fallback(c_e: float, T: float = 298.0) -> float:
    # c_e [mol/m^3] -> y=mol/L
    y = c_e / 1000.0
    return 0.1297 * y**3 - 2.51 * y ** 1.5 + 3.329 * y


def De_SI_fallback(c_e: float, T: float = 298.0) -> float:
    y = c_e / 1000.0
    return 8.794e-11 * y**2 - 3.972e-10 * y + 4.862e-10


# ---------------- 从 PyBaMM 参数集中获取函数的辅助 ----------------

def _try_get_callable(param, key: str) -> Optional[Callable]:
    """尽量从 ParameterValues 中拿到可直接调用的函数；失败返回 None。"""
    try:
        val = param[key]
        if callable(val):
            return val
    except Exception:
        pass
    return None


def _try_eval_with_signatures(fn: Callable, args: Tuple):
    """尝试用不同签名调用函数，返回第一个成功的结果，否则抛出异常。"""
    tried = []
    # 常见签名： (x), (x,T), (c,T)
    candidates = [
        (args[0:1],),
        (args[0:2],),
    ]
    for cand in candidates:
        try:
            return fn(*cand[0])
        except Exception as e:  # 继续尝试其它签名
            tried.append(str(e))
            continue
    # 若都失败，抛出最后一次异常信息
    raise RuntimeError("无法用常见签名调用函数: " + "; ".join(tried))


# ---------------- 主导出逻辑 ----------------

def export_scalars(param, out_dir: str):
    keys = [
        # 几何/孔隙/Bruggeman
        "Positive electrode thickness [m]",
        "Negative electrode thickness [m]",
        "Separator thickness [m]",
        "Positive electrode porosity",
        "Negative electrode porosity",
        "Separator porosity",
        "Positive electrode Bruggeman coefficient (electrode)",
        "Negative electrode Bruggeman coefficient (electrode)",
        "Bruggeman coefficient (electrolyte)",
        # 颗粒/浓度
        "Positive particle radius [m]",
        "Negative particle radius [m]",
        "Positive electrode maximum concentration in particles [mol.m-3]",
        "Negative electrode maximum concentration in particles [mol.m-3]",
        # 反应/电解液
        "Reference temperature [K]",
        "Ambient temperature [K]",
        "Electrolyte diffusivity [m2.s-1]",  # 若为数值常数
        "Electrolyte conductivity [S.m-1]",  # 若为数值常数
        "Reference exchange current density [A.m-2]",  # 有的模型定义存在
        "Typical electrolyte concentration [mol.m-3]",
        "Lithium metal partial molar volume [m3.mol-1]",
        "Cation transference number",
    ]

    fp = os.path.join(out_dir, "lgm50t_spme_scalars.csv")
    with open(fp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "value"])
        for k in keys:
            try:
                v = param[k]
                # 仅导出标量；函数跳过
                if isinstance(v, (int, float)):
                    w.writerow([k, v])
            except Exception:
                # 某些键在不同版本中名称略有差别，缺失则忽略
                continue
    print(f"[scalars] 写入: {fp}")


def export_ocp_curves(param, out_dir: str, npts: int = 2001):
    # 优先从参数集拿函数；否则回退
    up_fn = _try_get_callable(param, "Positive electrode OCP [V]") if HAS_PYBAMM else None
    un_fn = _try_get_callable(param, "Negative electrode OCP [V]") if HAS_PYBAMM else None

    if up_fn is None:
        up_fn = U_p_V_fallback
    if un_fn is None:
        un_fn = U_n_V_fallback

    fp = os.path.join(out_dir, "lgm50t_ocp_curves.csv")
    with open(fp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x", "Up_V", "Un_V", "OCV_V"])
        for x in linspace(0.001, 0.999, npts):
            try:
                up = up_fn(x)
            except Exception:
                up = U_p_V_fallback(x)
            try:
                un = un_fn(x)
            except Exception:
                un = U_n_V_fallback(x)
            w.writerow([f"{x:.8f}", f"{up:.8f}", f"{un:.8f}", f"{(up-un):.8f}"])
    print(f"[ocp] 写入: {fp}")


def export_electrolyte_props(param, out_dir: str, c_min=100.0, c_max=2000.0, npts: int = 100):
    kappa_fn = _try_get_callable(param, "Electrolyte conductivity [S.m-1]") if HAS_PYBAMM else None
    De_fn = _try_get_callable(param, "Electrolyte diffusivity [m2.s-1]") if HAS_PYBAMM else None

    if kappa_fn is None:
        kappa_fn = kappa_SI_fallback
    if De_fn is None:
        De_fn = De_SI_fallback

    fp = os.path.join(out_dir, "lgm50t_electrolyte_props.csv")
    with open(fp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["c_e_molm3", "T_K", "kappa_S_per_m", "De_m2_per_s"])
        for T in (273.15, 298.0, 323.15):
            for c in linspace(c_min, c_max, npts):
                # 尝试 (c), (c,T)
                try:
                    kappa = _try_eval_with_signatures(kappa_fn, (c, T))
                except Exception:
                    kappa = kappa_SI_fallback(c, T)
                try:
                    De = _try_eval_with_signatures(De_fn, (c, T))
                except Exception:
                    De = De_SI_fallback(c, T)
                w.writerow([f"{c:.3f}", f"{T:.2f}", f"{kappa:.6e}", f"{De:.6e}"])
    print(f"[elec] 写入: {fp}")


def export_solid_diffusivity(param, out_dir: str, n_x: int = 200, T_list = (273.15, 298.0, 323.15)):
    # 先尝试直接可调用对象
    Dn_fn = _try_get_callable(param, "Negative electrode diffusivity [m2.s-1]") if HAS_PYBAMM else None
    Dp_fn = _try_get_callable(param, "Positive electrode diffusivity [m2.s-1]") if HAS_PYBAMM else None

    # 若不可调用，再尝试读取原始数值常数
    if HAS_PYBAMM and Dn_fn is None:
        try:
            raw = param["Negative electrode diffusivity [m2.s-1]"]
            if isinstance(raw, (int, float)):
                c = float(raw)
                Dn_fn = lambda x, T=298.0: c
        except Exception:
            pass
    if HAS_PYBAMM and Dp_fn is None:
        try:
            raw = param["Positive electrode diffusivity [m2.s-1]"]
            if isinstance(raw, (int, float)):
                c = float(raw)
                Dp_fn = lambda x, T=298.0: c
        except Exception:
            pass

    fp = os.path.join(out_dir, "lgm50t_solid_diffusivity.csv")
    with open(fp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x", "T_K", "Dn_m2_per_s", "Dp_m2_per_s"])
        for T in T_list:
            for x in linspace(0.001, 0.999, n_x):
                # 尝试 (x), (x,T)
                if Dn_fn is None:
                    Dn = ""
                else:
                    try:
                        Dn = _try_eval_with_signatures(Dn_fn, (x, T))
                    except Exception:
                        # 若签名不匹配或不可用，留空
                        Dn = ""
                if Dp_fn is None:
                    Dp = ""
                else:
                    try:
                        Dp = _try_eval_with_signatures(Dp_fn, (x, T))
                    except Exception:
                        Dp = ""
                w.writerow([f"{x:.6f}", f"{T:.2f}",
                            (f"{Dn:.6e}" if Dn != "" else ""),
                            (f"{Dp:.6e}" if Dp != "" else "")])
    print(f"[Ds] 写入: {fp}")


def main():
    repo_root = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(repo_root, "JuBat", "output")
    ensure_dir(out_dir)

    if not HAS_PYBAMM:
        print("[警告] 未检测到 PyBaMM；将使用回退公式，仅导出 OCP 与电解液 κ/De 的常见近似。")
        param = {}
    else:
        # 选择 LGM50T（Chen2020）参数集
        try:
            param = pybamm.ParameterValues("Chen2020")
        except Exception:
            # 旧版本可用 chemistry 字典
            chem = pybamm.parameter_sets.Chen2020
            param = pybamm.ParameterValues(chemistry=chem)

    export_scalars(param, out_dir)
    export_ocp_curves(param, out_dir)
    export_electrolyte_props(param, out_dir)
    export_solid_diffusivity(param, out_dir)

    print("完成 LGM50T（Chen2020）SPMe 参数导出。")


if __name__ == "__main__":
    main()
