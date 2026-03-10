#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制电化学-热耦合（SPMe + 热Q4）代表性单元示意图：
- 上方：电化学 L2 元素（两节点+电流方向）
- 下方：热 Q4 元素（四节点），在单元内部绘制热容-热阻(RC)等效示意
- 中间：耦合箭头：电化学 -> 热 (q_gen 发热)，热 -> 电化学 (温度 T 影响参数)

输出：figure_echem_thermal_coupling.png 与 .svg

依赖：matplotlib、numpy（标准科学栈）
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, FancyBboxPatch
from matplotlib.lines import Line2D

# 字体设置：Windows 优先使用微软雅黑/黑体，避免中文缺字；并允许负号正常显示
plt.rcParams.setdefault('font.sans-serif', ['Microsoft YaHei', 'SimHei', 'DejaVu Sans'])
plt.rcParams.setdefault('axes.unicode_minus', False)


def draw_echem_l2(ax, start=(0.16, 0.74), length=0.68):
    """绘制电化学 L2 元素（两节点 + 电流箭头）。

    返回用于耦合箭头定位的包围盒 (x, y, w, h)。
    """
    x0, y0 = start
    x1 = x0 + length
    # 连接线
    ax.add_line(Line2D([x0, x1], [y0, y0], color="#1f77b4", lw=2))
    # 两端节点
    ax.add_patch(Circle((x0, y0), 0.01, facecolor="#1f77b4", edgecolor="none"))
    ax.add_patch(Circle((x1, y0), 0.01, facecolor="#1f77b4", edgecolor="none"))
    ax.text(x0, y0 + 0.03, "node i", ha="center", va="bottom", fontsize=10, color="#1f77b4")
    ax.text(x1, y0 + 0.03, "node i+1", ha="center", va="bottom", fontsize=10, color="#1f77b4")

    # L2 标签
    ax.text((x0 + x1) / 2, y0 + 0.085, "L2 element (electrochemical)",
        ha="center", va="bottom", fontsize=11, color="#1f77b4")

    # 电流方向箭头
    arr = FancyArrowPatch((x0 - 0.05, y0), (x1 + 0.05, y0),
                          arrowstyle="->", mutation_scale=15, lw=2, color="#d62728")
    ax.add_patch(arr)
    ax.text(x1 + 0.07, y0, "Current I", ha="left", va="center", fontsize=10.5, color="#d62728")

    return (x0, y0 - 0.02, length, 0.1)


def draw_thermal_q4(ax, origin=(0.32, 0.29), width=0.36, height=0.24):
    """绘制单个热 Q4 元素（四节点），并返回包围盒 (x, y, w, h)。"""
    x0, y0 = origin
    # 单元外框
    ax.add_patch(Rectangle((x0, y0), width, height, fill=False, lw=2, ec="#2ca02c"))
    # 四个角节点
    corners = [
        (x0, y0), (x0 + width, y0), (x0 + width, y0 + height), (x0, y0 + height)
    ]
    for (cx, cy) in corners:
        ax.add_patch(Circle((cx, cy), 0.01, facecolor="#2ca02c", edgecolor="none"))
    # 标签
    ax.text(x0 + width/2, y0 - 0.04, "Q4 element (thermal)",
        ha="center", va="top", fontsize=11, color="#2ca02c")

    # 在单元内部绘制 RC 示意（使用单元内部 70% 大小的内框作为布局参考）
    inner = (x0 + 0.15 * width, y0 + 0.2 * height, 0.7 * width, 0.6 * height)
    draw_rc_inside(ax, inner, label="Equivalent RC")
    return (x0, y0, width, height)


def draw_rc_inside(ax, bbox=(0.25, 0.2, 0.18, 0.12), label: str | None = None):
    """在给定矩形区域内部绘制 RC 等效示意（不画外框）。

    bbox: (x, y, w, h) 画布归一化坐标
    label: 可选的文字标签
    """
    x, y, w, h = bbox
    if label:
        ax.text(x + w/2, y + h + 0.02, label, ha="center", va="bottom", fontsize=10)

    # 节点温度 T_i
    node_x = x + 0.2 * w
    node_y = y + 0.6 * h
    ax.add_patch(Circle((node_x, node_y), 0.006, color="black"))
    ax.text(node_x - 0.01, node_y + 0.035, "T_i", ha="right", va="bottom", fontsize=10)

    # 向下线到电容
    ax.add_line(Line2D([node_x, node_x], [node_y, y + 0.28 * h], color="black", lw=1))

    # 电容 C（两条平行线）
    cap_y = y + 0.25 * h
    ax.add_line(Line2D([node_x - 0.03, node_x + 0.03], [cap_y, cap_y], color="black", lw=2))
    ax.add_line(Line2D([node_x - 0.03, node_x + 0.03], [cap_y - 0.02, cap_y - 0.02], color="black", lw=2))
    ax.text(node_x + 0.04, cap_y - 0.01, "C", ha="left", va="center", fontsize=10)

    # 下方“地”
    gnd_y = y + 0.14 * h
    ax.add_line(Line2D([node_x - 0.04, node_x + 0.04], [gnd_y, gnd_y], color="black", lw=1))
    ax.add_line(Line2D([node_x - 0.03, node_x + 0.03], [gnd_y - 0.01, gnd_y - 0.01], color="black", lw=1))
    ax.add_line(Line2D([node_x - 0.02, node_x + 0.02], [gnd_y - 0.02, gnd_y - 0.02], color="black", lw=1))

    # 向右通过热阻 R
    r_x1 = node_x + 0.03
    r_x2 = x + 0.75 * w
    ax.add_line(Line2D([node_x, r_x1], [node_y, node_y], color="black", lw=1))
    # 锯齿表示电阻
    zigs = 6
    xs = np.linspace(r_x1, r_x2, zigs * 2 + 1)
    ys = np.array([node_y + (0.01 if i % 2 == 0 else -0.01) for i in range(len(xs))])
    ys[0] = ys[-1] = node_y
    ax.plot(xs, ys, color="black", lw=1)
    ax.text((r_x1 + r_x2)/2, node_y + 0.03, "R", ha="center", va="bottom", fontsize=10)

    # 右端边界（与相邻单元/对流边界）
    ax.add_line(Line2D([r_x2, x + 0.93 * w], [node_y, node_y], color="black", lw=1))
    ax.add_line(Line2D([x + 0.93 * w, x + 0.93 * w], [y + 0.22 * h, y + 0.78 * h], color="black", lw=1))


def draw_coupling_side_arrows(ax, echem_box, thermal_box):
    """Draw coupling arrows on the left and right sides between models.

    Left: electrochem -> thermal (q_gen) downward.
    Right: thermal -> electrochem (T feedback) upward.
    """
    ex, ey, ew, eh = echem_box
    tx, ty, tw, th = thermal_box

    left_x = max(0.05, min(ex, tx) - 0.05)
    right_x = min(0.95, max(ex + ew, tx + tw) + 0.05)

    # Left downward arrow: q_gen
    start = (left_x, ey - 0.005)
    end = (left_x, ty + th + 0.005)
    arr_q = FancyArrowPatch(start, end, arrowstyle="->", mutation_scale=15, lw=2, color="#ff7f0e")
    ax.add_patch(arr_q)
    ax.text(left_x - 0.015, (start[1] + end[1]) / 2, "q_gen = q_ohmic + q_rxn + q_entropic",
        ha="right", va="center", fontsize=9.8, color="#ff7f0e", rotation=90)

    # Right upward arrow: T feedback
    start2 = (right_x, ty + th + 0.005)
    end2 = (right_x, ey - 0.005)
    arr_t = FancyArrowPatch(start2, end2, arrowstyle="->", mutation_scale=15, lw=2, color="#9467bd")
    ax.add_patch(arr_t)
    ax.text(right_x + 0.015, (start2[1] + end2[1]) / 2, "T affects: D_s(T), \u03BA_e(T), \u03C3_s(T), k(T), \u2202U/\u2202T",
        ha="left", va="center", fontsize=9.8, color="#9467bd", rotation=270)


def add_equation_box(ax, x, y, lines, color="#333333", align="left", where="top"):
    """Add a rounded text box with model equations.

    x, y in figure normalized coords (0..1). 'where' controls text anchor.
    """
    ha = "left" if align == "left" else "center" if align == "center" else "right"
    va = "top" if where == "top" else "bottom"
    text = "\n".join(lines)
    ax.text(x, y, text,
        ha=ha, va=va, fontsize=10, color=color,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, alpha=0.9))


def main():
    fig = plt.figure(figsize=(10, 7), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # 绘制代表性单元
    l2_box = draw_echem_l2(ax, start=(0.16, 0.74), length=0.68)
    q4_box = draw_thermal_q4(ax, origin=(0.32, 0.29), width=0.36, height=0.24)

    # 耦合箭头（左右两侧）
    draw_coupling_side_arrows(ax, echem_box=l2_box, thermal_box=q4_box)

    # 模型控制方程文本框（上：电化学；下：热）
    ec_lines = [
        "Electrochem (SPMe, L2):",
        "\u2202I/\u2202x = a_s F j",
        "I = -\u03C3_s \u2202\u03C6_s/\u2202x",
        "j = i0(T) sinh(\u03B1 F \u03B7 / (R T))",
        "\u2202c_s/\u2202t = D_s \u2202\u00B2c_s/\u2202x\u00B2",
    ]
    th_lines = [
        "Thermal (Q4, 2D):",
        "\u03C1 c_p \u2202T/\u2202t = \u2207\u00B7(k \u2207T) + q_gen",
        "-k \u2207T \u00B7 n = h (T - T_\u221E) + q_rad",
    ]
    add_equation_box(ax, x=0.08, y=0.92, lines=ec_lines, color="#1f77b4", align="left", where="top")
    add_equation_box(ax, x=0.08, y=0.20, lines=th_lines, color="#2ca02c", align="left", where="top")

    # 输出
    out_png = "figure_echem_thermal_coupling.png"
    out_svg = "figure_echem_thermal_coupling.svg"
    plt.savefig(out_png, bbox_inches='tight')
    plt.savefig(out_svg, bbox_inches='tight')
    print(f"Saved: {os.path.abspath(out_png)}")
    print(f"Saved: {os.path.abspath(out_svg)}")


if __name__ == "__main__":
    main()
