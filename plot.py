#!/usr/bin/env python3
"""
Plot a paper-style illustration showing coupling between a 1D SPMe
electrochemical model (through-thickness) and a 2D distributed thermal
conduction model (top view jellyroll mesh).

- Left: Cylindrical jellyroll top view with discretized thermal elements.
- Middle: SPMe 1D domains (NE | SP | PE) and key fields (j, η, φ_s, φ_e).
- Right: 2D thermal element with lumped mass (ρc) and anisotropic R(λ) to neighbors.

Outputs:
- figure_model_coupling.png
- figure_model_coupling.svg

Dependencies:
- matplotlib, numpy

Run:
    python plot.py
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Arc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import csv
from matplotlib import patheffects as pe


def add_label(ax, xy, text, ha="center", va="center", fontsize=11, weight="bold"):
    ax.text(
        *xy, text, ha=ha, va=va, fontsize=fontsize, weight=weight,
        path_effects=[pe.withStroke(linewidth=3, foreground="white")],
    )


def _polar_point(r: float, theta: float, center=(0.0, 0.0)):
    """Utility: polar to Cartesian on given center."""
    return (center[0] + r * np.cos(theta), center[1] + r * np.sin(theta))


def draw_thermal_bc_on_annulus(
    ax,
    r_in: float,
    r_out_vis: float,
    center=(0.0, 0.0),
    show_legend: bool = True,
    *,
    legend_loc: str = "upper left",
    legend_bbox_to_anchor=(1.02, 1.0, 1.0, 1.0),
    legend_width: float = 2.0,
    legend_height: float = 1.4,
):
    """
    Draw thermal boundary annotations on an annulus-like sketch:
    - Outer boundary (radius ~ r_out_vis): convection arrows and label
    - Inner boundary (radius ~ r_in): adiabatic ticks and label
    - Optional legend box (Convection / Adiabatic / Dirichlet / Radiation)
    """
    # Outer convection arrows
    thetas = np.deg2rad(np.array([ -60, -25, 10, 45, 80 ], dtype=float))
    L = 0.28  # arrow length outward
    for th in thetas:
        x0, y0 = _polar_point(r_out_vis, th, center)
        x1, y1 = _polar_point(r_out_vis + L, th, center)
        ax.annotate(
            "",
            xy=(x1, y1), xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#1f77b4"),
            zorder=5,
        )
    # Outer label
    ax.text(
        center[0] + r_out_vis + L + 0.1,
        center[1] + 0.0,
        r"对流: $-k\,\partial T/\partial n = h\,(T-T_{amb})$",
        fontsize=9,
        color="#5b9bd5",
        va="center",
        alpha=0.9,
    )

    # Inner adiabatic ticks (short radial ticks)
    tick_thetas = np.deg2rad(np.array([0, 90, 180, 270]))
    tL = 0.12
    for th in tick_thetas:
        ux, uy = np.cos(th), np.sin(th)
        x0, y0 = center[0] + (r_in - tL) * ux, center[1] + (r_in - tL) * uy
        x1, y1 = center[0] + (r_in + tL) * ux, center[1] + (r_in + tL) * uy
        ax.plot([x0, x1], [y0, y1], color="#2ca02c", lw=2.0, zorder=5)
    # Inner label
    ax.text(
        center[0] - r_out_vis - 0.1,
        center[1] - 0.0,
        r"内侧绝热: $-k\,\partial T/\partial n = 0$",
        fontsize=9,
        color="#2ca02c",
        ha="right",
        va="center",
    )

    # Draw a faint inner boundary circle to clarify the location
    ax.add_patch(Circle(center, r_in, fill=False, ec="#7f7f7f", lw=1.0, ls=(0, (4, 3)), alpha=0.9))

    if show_legend:
        # Small legend panel placed outside upper-right by default to avoid overlap
        # Convert 4-tuple to Bbox for robust handling of percentage sizes
        try:
            from matplotlib.transforms import Bbox as _Bbox
            _bbox = _Bbox.from_bounds(*legend_bbox_to_anchor) if isinstance(legend_bbox_to_anchor, (tuple, list)) and len(legend_bbox_to_anchor) == 4 else legend_bbox_to_anchor
        except Exception:
            _bbox = legend_bbox_to_anchor
        legend_ax = inset_axes(
            ax,
            width=legend_width,
            height=legend_height,
            loc=legend_loc,
            bbox_to_anchor=_bbox,
            bbox_transform=ax.transAxes,
            borderpad=0.4,
        )
        legend_ax.axis("off")
        legend_ax.set_xlim(0, 100)
        legend_ax.set_ylim(0, 100)
        # Panel frame
        legend_ax.add_patch(Rectangle((2, 2), 96, 96, fill=True, ec="#e0e0e0", fc="#f8f9fa", lw=0.8))
        # Convection item
        legend_ax.annotate(
            "", xy=(32, 82), xytext=(8, 82),
            arrowprops=dict(arrowstyle="-|>", lw=2.0, color="#1f77b4"),
        )
        legend_ax.text(40, 84, "对流", fontsize=9, weight="bold")
        legend_ax.text(40, 72, r"q = h (T - T_{amb})", fontsize=8, family="monospace")
        # Adiabatic item
        legend_ax.plot([8, 32], [52, 52], color="#2ca02c", lw=2.0)
        legend_ax.plot([14, 14], [44, 60], color="#2ca02c", lw=2.0)
        legend_ax.plot([26, 26], [44, 60], color="#2ca02c", lw=2.0)
        legend_ax.text(40, 54, "绝热", fontsize=9, weight="bold")
        legend_ax.text(40, 42, r"q = 0  ⇔  ∂T/∂n = 0", fontsize=8, family="monospace")
        # Dirichlet item (isothermal)
        legend_ax.add_patch(Rectangle((8, 18), 24, 16, fill=True, ec="#64b5f6", fc="#e3f2fd", lw=1.4))
        legend_ax.text(40, 30, "恒温", fontsize=9, weight="bold")
        legend_ax.text(40, 18, r"T = T_{spec}", fontsize=8, family="monospace")
        # Radiation (optional)
        circ = Circle((20, 8), 6, fill=False, ec="#fb8c00", lw=1.6)
        legend_ax.add_patch(circ)
        legend_ax.text(40, 10, "辐射(可选)", fontsize=9, weight="bold")
        legend_ax.text(40, -2, r"q = εσ (T^4 - T_{amb}^4)", fontsize=8, family="monospace")


def _arch_spiral(a: float, b: float, theta: np.ndarray, center=(0, 0)):
    r = a + b * theta
    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    return x, y, r


def draw_collector_seeded_band_mesh(
    ax,
    center=(0, 0),
    r_in=0.25,
    r_out=1.9,
    turns=3,
    nbands=6,
    seg_per_turn=24,
    phase=0.0,
    rail_lw=1.0,
    edge_lw=0.6,
    rail_color="#A45A2A",
    edge_color="#444444",
    node_stride=6,
    show_layers=False,
    layer_fracs=(0.12, 0.30, 0.16, 0.30, 0.12),
    layer_names=("PCC", "PE", "SP", "NE", "NCC"),
    layer_colors=("#8c564b", "#ff7f0e", "#7f7f7f", "#1f77b4", "#5B5B5B"),
    layer_tick_stride=14,
    fill_layers=False,
    fill_alpha=0.35,
    fill_homogenized=False,
    overlay_layers_band_index=None,
    overlay_tick_stride=2,
    overlay_tick_lw=2.4,
    only_band_index=None,
    force_layered=False,
):
    """
    基于集流体"导轨"的条带网格生成器（以直代曲）
    
    网格生成原理（仿照 Jellyrollmodel.jl）:
    - 阿基米德螺旋方程: r(θ) = a + bθ，其中 a=Rin, b=t_repeat/(2π)
    - 沿内外螺旋线等角度采样：
      * 内螺旋: r_in(θ) = a + bθ + s_in（PCC内侧）
      * 外螺旋: r_out(θ) = a + bθ + s_out（NCC外侧）
    - 连接对应采样点形成Q4条带单元: [(in_i, out_i, out_{i+1}, in_{i+1})]
    - 通过径向堆叠多个条带填充整个环形区域
    
    层序定义（从内到外）:
    - PCC (正极集流体) → PE (正极) → SP (隔膜) → NE (负极) → NCC (负极集流体)
    
    参数:
    - r_in, r_out: 内外半径 [m]
    - turns: 螺旋圈数
    - nbands: 径向条带数量
    - seg_per_turn: 每圈的周向分段数（建议≥24）
    - phase: 相位对齐角度 [rad]
    - show_layers: 是否显示层序标记
    - fill_layers: 是否填充各层颜色
    - fill_homogenized: 是否使用均质化颜色
    """
    theta_max = 2 * np.pi * turns
    
    # 螺旋参数（与Julia代码保持一致）
    # a: 螺旋起始半径（Rin）
    # b: 螺旋节距系数 = t_repeat / (2π)
    a = r_in
    t_repeat = (r_out - r_in) / max(1, nbands)  # 估算单层厚度
    b = t_repeat / (2.0 * np.pi)

    # 将环形区域划分为nbands个径向条带（每个条带厚度相等）
    dr = (r_out - r_in) / nbands
    dtheta = 2 * np.pi / float(max(1, seg_per_turn))  # 周向角度步长

    # 用于颜色填充的多边形集合
    poly_verts = []
    poly_fc = []

    # 圈间径向节距 Δr_turn = 2πb（螺旋每转一圈的径向增量）
    # 当条带厚度 dr > Δr_turn 时，存在视觉重叠风险
    delta_r_turn = 2 * np.pi * b
    do_homog = bool(fill_homogenized)
    
    # 自动切换到均质化填充以避免视觉重叠
    if fill_layers and (dr > delta_r_turn) and (not force_layered):
        do_homog = True
        print(f"[注意] dr={dr:.3f} > Δr_turn={delta_r_turn:.3f}: 切换到均质化填充以避免视觉重叠")

    # 预计算均质化颜色（按层体积分数加权混合）
    def _mix_color_hex(hex_colors, weights):
        import matplotlib.colors as mcolors
        w = np.asarray(weights, dtype=float)
        w = w / np.sum(w)
        cols = np.array([mcolors.to_rgb(c) for c in hex_colors], dtype=float)
        rgb = np.clip(np.sum(cols * w[:, None], axis=0), 0, 1)
        return rgb
    homog_color = _mix_color_hex(layer_colors, layer_fracs)

    # 遍历每个径向条带
    for k in range(nbands):
        # 仅绘制单条带模式（用于演示）
        if (only_band_index is not None) and (k != int(only_band_index)):
            continue
        
        # 当前条带的径向偏移范围 [s_in, s_out]
        s_in = k * dr
        s_out = (k + 1) * dr
        
        # θ 范围裁剪（确保内外螺旋都在 [r_in, r_out] 内）
        # θ0: 起始角度，确保 r_in(θ0) ≥ r_in
        # θ1: 终止角度，确保 r_out(θ1) ≤ r_out
        theta0 = max(0.0, (r_in - a - s_in) / b)
        theta1_lim = (r_out - a - s_out) / b
        if theta1_lim <= theta0:
            continue
        theta1 = min(theta_max, theta1_lim)
        
        # 等角度网格相位对齐（与Julia代码保持一致）
        k0 = int(np.ceil((theta0 - phase) / dtheta))
        k1 = int(np.floor((theta1 - phase) / dtheta))
        if k1 <= k0:
            # 备选方案：至少两个采样点以绘制短条带段
            theta = np.array([theta0, min(theta1, theta0 + dtheta*0.75)])
        else:
            theta = phase + dtheta * np.arange(k0, k1 + 1)

        # 生成内外螺旋线上的采样点
        xin, yin, _ = _arch_spiral(a + s_in, b, theta, center)
        xout, yout, _ = _arch_spiral(a + s_out, b, theta, center)

        # 绘制螺旋"导轨"（避免相邻条带间重复绘制共享导轨）
        if only_band_index is not None:
            # 单条带模式：绘制选中条带的内外两条导轨
            ax.plot(xin, yin, color=rail_color, lw=rail_lw, alpha=0.9)
            ax.plot(xout, yout, color=rail_color, lw=rail_lw, alpha=0.9)
        else:
            # 多条带模式：仅绘制最内侧导轨和每个条带的外侧导轨
            if k == 0:
                ax.plot(xin, yin, color=rail_color, lw=rail_lw, alpha=0.9)
            ax.plot(xout, yout, color=rail_color, lw=rail_lw, alpha=0.9)

        # 生成Q4条带单元（连接内外导轨形成四边形网格）
        # 单元节点顺序: [n1(内侧当前), n2(外侧当前), n3(外侧下一个), n4(内侧下一个)]
        for i in range(len(theta) - 1):
            # 横向边（连接内外导轨）
            ax.plot([xin[i], xout[i]], [yin[i], yout[i]], color=edge_color, lw=edge_lw)
            ax.plot([xin[i + 1], xout[i + 1]], [yin[i + 1], yout[i + 1]], color=edge_color, lw=edge_lw)
            # 纵向边（沿螺旋线方向连接相邻采样点）
            ax.plot([xout[i], xout[i + 1]], [yout[i], yout[i + 1]], color=edge_color, lw=edge_lw)
            ax.plot([xin[i], xin[i + 1]], [yin[i], yin[i + 1]], color=edge_color, lw=edge_lw)

            # 可选：分层填充（每个子层使用单独的颜色）
            # 将每个条带按层序体积分数细分，生成多个子Q4多边形
            if fill_layers and (not do_homog):
                fracs = np.asarray(layer_fracs, dtype=float)
                fracs = fracs / fracs.sum()
                cum = np.concatenate(([0.0], np.cumsum(fracs)))
                th_i = theta[i]
                th_ip1 = theta[i + 1]
                for j in range(len(fracs)):
                    # 子层径向范围 [s_lo, s_hi]
                    s_lo = s_in + cum[j] * dr
                    s_hi = s_in + cum[j + 1] * dr
                    # 四个角点坐标
                    x_lo_i, y_lo_i, _ = _arch_spiral(a + s_lo, b, np.array([th_i]))
                    x_hi_i, y_hi_i, _ = _arch_spiral(a + s_hi, b, np.array([th_i]))
                    x_lo_ip1, y_lo_ip1, _ = _arch_spiral(a + s_lo, b, np.array([th_ip1]))
                    x_hi_ip1, y_hi_ip1, _ = _arch_spiral(a + s_hi, b, np.array([th_ip1]))
                    poly_verts.append([
                        (x_lo_i[0], y_lo_i[0]),
                        (x_hi_i[0], y_hi_i[0]),
                        (x_hi_ip1[0], y_hi_ip1[0]),
                        (x_lo_ip1[0], y_lo_ip1[0]),
                    ])
                    poly_fc.append(layer_colors[j % len(layer_colors)])

            # 均质化填充：每个单元一个Q4多边形（全厚度），使用混合颜色
            if do_homog:
                th_i = theta[i]
                th_ip1 = theta[i + 1]
                x_lo_i, y_lo_i, _ = _arch_spiral(a + s_in, b, np.array([th_i]))
                x_hi_i, y_hi_i, _ = _arch_spiral(a + s_out, b, np.array([th_i]))
                x_lo_ip1, y_lo_ip1, _ = _arch_spiral(a + s_in, b, np.array([th_ip1]))
                x_hi_ip1, y_hi_ip1, _ = _arch_spiral(a + s_out, b, np.array([th_ip1]))
                poly_verts.append([
                    (x_lo_i[0], y_lo_i[0]),
                    (x_hi_i[0], y_hi_i[0]),
                    (x_hi_ip1[0], y_hi_ip1[0]),
                    (x_lo_ip1[0], y_lo_ip1[0]),
                ])
                poly_fc.append(homog_color)

        # 稀疏节点标记（用于视觉清晰度）
        if node_stride > 0 and len(theta) > 0:
            sel = np.arange(0, len(theta), max(1, node_stride))
            ax.scatter(xin[sel], yin[sel], s=8, color="black", zorder=3)
            ax.scatter(xout[sel], yout[sel], s=8, color="black", zorder=3)

    # 可选：使用彩色径向刻度线显示层序卷绕结构
    if show_layers and (not fill_layers) and len(theta) > 0:
            # 从条带内侧开始的累积偏移
            fracs = np.asarray(layer_fracs, dtype=float)
            fracs = fracs / fracs.sum()
            cum = np.concatenate(([0.0], np.cumsum(fracs)))  # 长度 = 层数+1
            # 选择角度采样点的子集来绘制刻度
            step = max(1, layer_tick_stride)
            for idx in range(0, len(theta), step):
                th = theta[idx]
                # 沿径向为每一层绘制短彩色线段
                for j in range(len(fracs)):
                    s_lo = s_in + cum[j] * dr
                    s_hi = s_in + cum[j + 1] * dr
                    # 该角度处层的内外边界点
                    x_lo, y_lo, _ = _arch_spiral(a + s_lo, b, np.array([th]), center)
                    x_hi, y_hi, _ = _arch_spiral(a + s_hi, b, np.array([th]), center)
                    ax.plot([x_lo[0], x_hi[0]], [y_lo[0], y_hi[0]],
                            color=layer_colors[j % len(layer_colors)], lw=2.2, alpha=0.95, zorder=2)

    # 仅在指定条带索引上叠加层序刻度（即使启用填充也适用）
    if overlay_layers_band_index is not None and k == int(overlay_layers_band_index) and len(theta) > 0:
            fracs = np.asarray(layer_fracs, dtype=float)
            fracs = fracs / fracs.sum()
            cum = np.concatenate(([0.0], np.cumsum(fracs)))
            step = max(1, overlay_tick_stride)
            for idx in range(0, len(theta), step):
                th = theta[idx]
                for j in range(len(fracs)):
                    s_lo = s_in + cum[j] * dr
                    s_hi = s_in + cum[j + 1] * dr
                    x_lo, y_lo, _ = _arch_spiral(a + s_lo, b, np.array([th]), center)
                    x_hi, y_hi, _ = _arch_spiral(a + s_hi, b, np.array([th]), center)
                    ax.plot([x_lo[0], x_hi[0]], [y_lo[0], y_hi[0]],
                            color=layer_colors[j % len(layer_colors)], lw=overlay_tick_lw, alpha=0.98, zorder=3)

    # 添加填充多边形（如果请求）
    if (fill_layers or do_homog) and len(poly_verts) > 0:
        from matplotlib.collections import PolyCollection
        coll = PolyCollection(poly_verts, facecolors=poly_fc, edgecolors='none', alpha=fill_alpha, zorder=1)
        ax.add_collection(coll)


def draw_cylinder_panel(ax):
    """
    绘制圆柱体俯视图面板
    
    展示集流体导轨条带Q4网格的俯视图，用于论文插图。
    """
    # 顶视图外边界
    R = 2.3
    ax.add_patch(Circle((0, 0), R, fill=False, ec="black", lw=1.5))

    # 集流体导轨条带Q4网格（以直代曲）
    draw_collector_seeded_band_mesh(
        ax,
        r_in=0.25,
        r_out=1.9,
        turns=3.2,
        nbands=6,
        seg_per_turn=24,
        phase=0.0,
        rail_lw=1.0,
        edge_lw=0.6,
        rail_color="#B05300",
        edge_color="#555555",
        node_stride=8,
    )

    # 平面坐标轴三联组 (x,y)
    ax.arrow(0, 0, 0.9, 0.0, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.arrow(0, 0, 0.0, 0.9, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.text(1.05, 0.0, "x", fontsize=10)
    ax.text(0.0, 1.05, "y", fontsize=10)

    ax.set_aspect("equal")
    ax.set_xlim(-2.7, 2.7)
    ax.set_ylim(-2.7, 2.7)
    ax.axis("off")
    add_label(ax, (0, -2.95), "集流体导轨条带Q4网格（俯视图）", fontsize=11)


def draw_resistor(ax, x, y, w=0.5, h=0.18, color="#1f77b4"):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fill=False, ec=color, lw=1.8))


def draw_capacitor(ax, x, y, gap=0.06, height=0.36, color="#1f77b4"):
    ax.plot([x - gap, x - gap], [y - height/2, y + height/2], color=color, lw=1.6)
    ax.plot([x + gap, x + gap], [y - height/2, y + height/2], color=color, lw=1.6)


def draw_spme_panel(ax):
    """Draw a compact SPMe 1D domains panel: NE | SP | PE with main symbols."""
    y0 = 0.0
    # Domain rectangles
    widths = [1.2, 0.6, 1.2]
    labels = ["NE", "SP", "PE"]
    colors = ["#1f77b4", "#7f7f7f", "#ff7f0e"]
    x_left = 0.0
    rects = []
    for w, lab, col in zip(widths, labels, colors):
        ax.add_patch(Rectangle((x_left, -0.6), w, 1.2, fill=False, ec=col, lw=1.8))
        add_label(ax, (x_left + w/2, -0.9), lab, fontsize=10, weight="bold")
        rects.append((x_left, w))
        x_left += w
    x_right = x_left

    # Through-thickness coordinate x
    ax.arrow(-0.5, -0.85, x_right + 0.8, 0, head_width=0.07, head_length=0.1, fc="k", ec="k")
    ax.text(x_right + 0.4, -0.95, "x", fontsize=10)

    # Applied current and interfacial current density j
    ax.arrow(-0.35, 0.0, 0.3, 0.0, head_width=0.07, head_length=0.1, fc="#1f77b4", ec="#1f77b4")
    ax.text(-0.45, 0.08, r"$I_{app}$", color="#1f77b4")
    # interface arrows j at NE/SP and SP/PE
    x_ne_sp = rects[0][0] + rects[0][1]
    x_sp_pe = rects[1][0] + rects[1][1]
    for xi in [x_ne_sp, x_sp_pe]:
        ax.arrow(xi, -0.1, 0.0, 0.2, head_width=0.06, head_length=0.08, fc="#2ca02c", ec="#2ca02c")
        ax.text(xi+0.03, 0.18, r"$j$", color="#2ca02c")

    # potentials φ_s (top line) and φ_e (bottom line)
    ax.plot([0, x_right], [0.45, 0.45], color="#1f77b4", lw=1.2)
    ax.plot([0, x_right], [-0.45, -0.45], color="#9467bd", lw=1.2)
    ax.text(0.02, 0.53, r"$\phi_s$", color="#1f77b4")
    ax.text(0.02, -0.60, r"$\phi_e$", color="#9467bd")

    # Overpotentials and heat terms
    ax.text(0.15, 0.05, r"$\eta = 2 T\,\mathrm{asinh}\left(\frac{j}{2 j_0}\right)$", fontsize=9)
    ax.text(0.15, -0.12, r"$q= q_{rxn} + q_{rev} + q_{ohm}$", fontsize=9)

    ax.set_xlim(-0.6, x_right + 0.9)
    ax.set_ylim(-1.2, 1.1)
    ax.axis("off")
    add_label(ax, (x_right/2, -1.05), "SPMe 1D electrochemical model", fontsize=11)


def _load_spiral_csv(dirname: str):
    """Load spiral example CSVs: nodes, elements, T_nodes (K), stress per element.
    Returns: nodes(Nx2), elements(Ex4 int-1based), T_nodes(N), stress_e(E)
    """
    def _read_csv(path):
        with open(path, newline='') as f:
            reader = csv.reader(f)
            return [list(map(float, row)) for row in reader]
    def _read_csv_int(path):
        with open(path, newline='') as f:
            reader = csv.reader(f)
            return [list(map(int, row)) for row in reader]

    nodes_rows = _read_csv(os.path.join(dirname, 'example', 'spiral_nodes.csv'))
    elems_rows = _read_csv_int(os.path.join(dirname, 'example', 'spiral_elements.csv'))
    T_rows = _read_csv(os.path.join(dirname, 'example', 'spiral_T_nodes_K.csv'))
    s_rows = _read_csv(os.path.join(dirname, 'example', 'spiral_stress_elem.csv'))

    nodes = np.array([[r[1], r[2]] for r in nodes_rows], dtype=float)
    elems = np.array([r[1:5] for r in elems_rows], dtype=int)  # 1-based indices
    Tn = np.array([r[1] for r in T_rows], dtype=float)
    se = np.array([r[1] for r in s_rows], dtype=float)
    return nodes, elems, Tn, se


def _apply_zoom(ax, nodes_xy: np.ndarray, frac: float = 0.65):
    """Zoom the axes to a fraction of the mesh bounding box around its center.
    frac in (0,1]; smaller means tighter zoom (enlarged cells)."""
    x = nodes_xy[:, 0]; y = nodes_xy[:, 1]
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    wx = (xmax - xmin) * frac
    wy = (ymax - ymin) * frac
    # keep aspect 1:1 by using the larger of wx, wy
    w = max(wx, wy)
    pad = 0.02 * w
    ax.set_xlim(cx - 0.5*w - pad, cx + 0.5*w + pad)
    ax.set_ylim(cy - 0.5*w - pad, cy + 0.5*w + pad)
    ax.set_aspect('equal')


def draw_thermal_panel(ax):
    """2D Thermal (Q4) panel — single illustrative element with colorbar (Kelvin)."""
    from matplotlib.collections import PolyCollection
    # One Q4 element centered at origin
    square = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    polys = [square]
    Tval = np.array([300.0])  # example temperature [K]
    coll = PolyCollection(polys, array=Tval, cmap='turbo', edgecolors='k', linewidths=1.0)
    # Stable color scale so the colorbar is meaningful
    coll.set_clim(293.0, 313.0)
    ax.add_collection(coll)
    ax.set_aspect('equal')
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    cax = inset_axes(ax, width="6%", height="70%", loc="right", borderpad=0.8)
    cb = plt.colorbar(coll, cax=cax)
    cb.set_label("T [K]")
    # Legend (proxy artist)
    leg_handle = Rectangle((0, 0), 0.2, 0.2, fill=False, ec='k', lw=1.0, label='Q4 单元')
    ax.legend(handles=[leg_handle], loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=9, frameon=True, framealpha=0.9)
    # Governing equation annotation
    ax.text(0.0, 0.72, r"$\rho c\,\partial T/\partial t=\nabla\cdot(\lambda\nabla T)+q$", ha='center', fontsize=10)
    ax.text(0.0, 0.58, r"BC: $-\mathbf{n}\cdot\lambda\nabla T=h\,(T-T_{amb})$", ha='center', fontsize=9, color='gray')
    ax.axis("off")
    add_label(ax, (0.0, -0.9), "2D Thermal (Q4)", fontsize=11)


def draw_mechanical_panel(ax):
    """2D Mechanical (Q4) panel — single illustrative element with colorbar (thermal stress)."""
    from matplotlib.collections import PolyCollection
    square = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    polys = [square]
    sval = np.array([1.0e6])  # example stress [Pa]
    coll = PolyCollection(polys, array=sval, cmap='coolwarm', edgecolors='k', linewidths=1.0)
    # Stable symmetric color scale for stress
    coll.set_clim(-5.0e6, 5.0e6)
    ax.add_collection(coll)
    ax.set_aspect('equal')
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    cax = inset_axes(ax, width="6%", height="70%", loc="right", borderpad=0.8)
    cb = plt.colorbar(coll, cax=cax)
    cb.set_label(r"$\sigma_{th}$ [Pa]")
    # Legend (proxy artist)
    leg_handle = Rectangle((0, 0), 0.2, 0.2, fill=False, ec='k', lw=1.0, label='Q4 单元')
    ax.legend(handles=[leg_handle], loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=9, frameon=True, framealpha=0.9)
    # Thermal stress equations
    ax.text(0.0, 0.72, r"$\varepsilon_{th}=\alpha\,\Delta T$", ha='center', fontsize=10)
    ax.text(0.0, 0.58, r"$\sigma=\mathbf{C}:(\varepsilon-\varepsilon_{th})$", ha='center', fontsize=10)
    ax.axis("off")
    add_label(ax, (0.0, -0.9), "2D Mechanical (Q4)", fontsize=11)


def draw_coupling_arrows_bidirectional(ax):
    """Bidirectional coupling arrows (echem ↔ thermal)."""
    arrow1 = FancyArrowPatch(
        posA=(0.0, 0.4), posB=(1.0, 0.4),
        connectionstyle="arc3,rad=0.25", arrowstyle="-|>", mutation_scale=12,
        lw=1.4, fc="gray", ec="gray", alpha=0.9,
    )
    arrow2 = FancyArrowPatch(
        posA=(1.0, -0.4), posB=(0.0, -0.4),
        connectionstyle="arc3,rad=-0.25", arrowstyle="-|>", mutation_scale=12,
        lw=1.4, fc="gray", ec="gray", alpha=0.9,
    )
    ax.add_patch(arrow1)
    ax.add_patch(arrow2)
    ax.text(0.5, 0.62, "q (heat)", ha="center", fontsize=9, color="gray")
    ax.text(0.5, -0.62, "T → kinetics", ha="center", fontsize=9, color="gray")


def draw_coupling_arrows_unidirectional(ax):
    """Unidirectional coupling arrow (thermal → mechanical)."""
    arrow = FancyArrowPatch(
        posA=(0.0, 0.0), posB=(1.0, 0.0),
        connectionstyle="arc3,rad=0.0", arrowstyle="-|>", mutation_scale=12,
        lw=1.4, fc="gray", ec="gray", alpha=0.9,
    )
    ax.add_patch(arrow)
    ax.text(0.5, 0.18, "T → σ_th", ha="center", fontsize=9, color="gray")


# === New: Individual figures per theory section ===
def figure_topview_thermal_mesh(
    r_in=0.25, r_out=1.9, turns=3.2, nbands=6, seg_per_turn=24, phase=0.0,
    save=True
):
    """
    生成热网格俯视图（集流体导轨条带Q4网格，带层填充）
    
    展示完整的螺旋网格结构，包含层序颜色填充和图例。
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    # 外边界圆
    R = 2.3
    ax.add_patch(Circle((0, 0), R, fill=False, ec="black", lw=1.5))
    draw_collector_seeded_band_mesh(
        ax,
        r_in=r_in, r_out=r_out, turns=turns, nbands=nbands,
        seg_per_turn=seg_per_turn, phase=phase,
        rail_lw=1.0, edge_lw=0.6, rail_color="#B05300", edge_color="#555555",
        node_stride=8, show_layers=False, fill_layers=True, fill_alpha=0.5,
        overlay_layers_band_index=max(0, nbands//2 - 1), overlay_tick_stride=2,
    )
    # 坐标轴三联组
    ax.arrow(0, 0, 0.9, 0.0, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.arrow(0, 0, 0.0, 0.9, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.text(1.05, 0.0, "x", fontsize=10)
    ax.text(0.0, 1.05, "y", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-2.7, 2.7)
    ax.set_ylim(-2.7, 2.7)
    ax.axis("off")
    add_label(ax, (0, -2.95), "热网格 (集流体导轨条带Q4网格；层填充)", fontsize=11)
    
    # 层序图例（颜色含义：各层颜色或均质化混合色）
    from matplotlib.lines import Line2D
    layer_names = ("PCC", "PE", "SP", "NE", "NCC")
    layer_colors = ("#8c564b", "#ff7f0e", "#7f7f7f", "#1f77b4", "#5B5B5B")
    handles = [Line2D([0], [0], color=c, lw=2.2, label=n) for n, c in zip(layer_names, layer_colors)]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(-0.02, 1.02), fontsize=8, frameon=True, framealpha=0.9, ncol=5)
    if save:
        for fn in ("figure_mesh_topview.png", "figure_mesh_topview.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def figure_single_spiral_layered(
    r_in=0.25, r_out=1.9, turns=3.2, nbands=6, seg_per_turn=48, phase=0.0,
    band_index=None, save=True
):
    """
    生成单条螺旋带示意图（热单元）
    
    突出显示单个螺旋条带的结构，配合热边界条件标注。
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    R = 2.3
    ax.add_patch(Circle((0, 0), R, fill=False, ec="black", lw=1.5))
    if band_index is None:
        band_index = max(0, nbands//2 - 1)
    draw_collector_seeded_band_mesh(
        ax,
        r_in=r_in, r_out=r_out, turns=turns, nbands=nbands,
        seg_per_turn=seg_per_turn, phase=phase,
        rail_lw=1.2, edge_lw=0.8, rail_color="#B05300", edge_color="#BBBBBB",
        node_stride=0,
        show_layers=False, fill_layers=False, fill_alpha=0.0, force_layered=False,
        only_band_index=band_index,
    )
    # 坐标轴三联组
    ax.arrow(0, 0, 0.9, 0.0, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.arrow(0, 0, 0.0, 0.9, head_width=0.08, head_length=0.12, fc="k", ec="k")
    ax.text(1.05, 0.0, "x", fontsize=10)
    ax.text(0.0, 1.05, "y", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-2.7, 2.7)
    ax.set_ylim(-2.7, 2.7)
    ax.axis("off")
    add_label(ax, (0, -2.95), "单条螺旋带（热单元）", fontsize=11)

    # Add thermal BC annotations (outer convection, inner adiabatic) + legend
    try:
        draw_thermal_bc_on_annulus(
            ax,
            r_in=r_in,
            r_out_vis=R,
            center=(0.0, 0.0),
            show_legend=True,
            legend_loc="upper left",
            legend_bbox_to_anchor=(1.02, 1.0, 1.0, 1.0),  # outside upper-right
            legend_width=2.0,
            legend_height=1.4,
        )
    except Exception as e:
        # non-fatal: ensure plotting continues even if legend drawing fails
        print("[warn] draw_thermal_bc_on_annulus failed:", e)
    # Mark one thermal element region with an "SPMe" label to indicate heat source comes from SPMe inside this thermal cell
    try:
        r_mid = 0.5 * (r_in + r_out * 0.9)
        theta_deg = 35.0
        px, py = _polar_point(r_mid, np.deg2rad(theta_deg))
        w, h = 0.28, 0.18  # illustrative footprint of one Q4 thermal element
        rect = Rectangle((px - w/2, py - h/2), w, h,
                         fill=False, ec="#d62728", lw=1.6, ls=(0, (4, 2)))
        ax.add_patch(rect)
        ax.text(px, py, "SPMe", ha="center", va="center", color="#d62728", fontsize=10,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    except Exception:
        pass
    if save:
        for fn in ("figure_single_spiral_layered.png", "figure_single_spiral_layered.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def figure_echem_1d_mesh(n_ne=8, n_sp=4, n_pe=8, save=True):
    fig, ax = plt.subplots(1, 1, figsize=(7.2, 3.2))
    # Domain spans
    widths = [1.2, 0.6, 1.2]
    labels = ["NE", "SP", "PE"]
    colors = ["#1f77b4", "#7f7f7f", "#ff7f0e"]
    x_left = 0.0
    x_marks = [x_left]
    rects = []
    for w, lab, col in zip(widths, labels, colors):
        ax.add_patch(Rectangle((x_left, -0.6), w, 1.2, fill=False, ec=col, lw=1.8))
        add_label(ax, (x_left + w/2, -0.9), lab, fontsize=10, weight="bold")
        rects.append((x_left, w))
        x_left += w
        x_marks.append(x_left)
    x_right = x_left
    # 1D mesh nodes/elements per region
    Ns = [n_ne, n_sp, n_pe]
    x0 = 0.0
    node_y = -0.02
    for (xl, w), N in zip(rects, Ns):
        xs = np.linspace(xl, xl + w, N + 1)
        # element vertical ticks
        for x in xs:
            ax.plot([x, x], [-0.6, 0.6], color="#CCCCCC", lw=0.8, alpha=0.9)
        # node markers on center line
        ax.scatter(xs, np.full_like(xs, node_y), s=14, color="black", zorder=3)
    # NCC/PCC and current
    ax.text(-0.08, 0.68, "NCC", fontsize=9)
    ax.text(x_right + 0.02, 0.68, "PCC", fontsize=9)
    ax.arrow(-0.25, 0.0, 0.22, 0.0, head_width=0.07, head_length=0.1, fc="#1f77b4", ec="#1f77b4")
    ax.text(-0.35, 0.08, r"$I_{app}$", color="#1f77b4")
    # potentials φ_s and φ_e
    ax.plot([0, x_right], [0.45, 0.45], color="#1f77b4", lw=1.2)
    ax.plot([0, x_right], [-0.45, -0.45], color="#9467bd", lw=1.2)
    ax.text(0.02, 0.53, r"$\phi_s$", color="#1f77b4")
    ax.text(0.02, -0.60, r"$\phi_e$", color="#9467bd")
    # j arrows at interfaces
    for xi in [rects[0][0] + rects[0][1], rects[1][0] + rects[1][1]]:
        ax.arrow(xi, -0.1, 0.0, 0.2, head_width=0.06, head_length=0.08, fc="#2ca02c", ec="#2ca02c")
        ax.text(xi + 0.03, 0.18, r"$j$", color="#2ca02c")
    # x-axis
    ax.arrow(-0.5, -0.85, x_right + 0.8, 0, head_width=0.07, head_length=0.1, fc="k", ec="k")
    ax.text(x_right + 0.4, -0.95, "x", fontsize=10)
    ax.set_xlim(-0.6, x_right + 0.9)
    ax.set_ylim(-1.2, 1.1)
    ax.axis("off")
    add_label(ax, (x_right/2, -1.05), "Electrochemical 1D mesh (NE|SP|PE)", fontsize=11)
    if save:
        for fn in ("figure_echem_mesh.png", "figure_echem_mesh.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def figure_thermal_element(save=True):
    fig, ax = plt.subplots(1, 1, figsize=(4.2, 4.2))
    draw_thermal_panel(ax)
    if save:
        for fn in ("figure_thermal_element.png", "figure_thermal_element.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def figure_mechanical_element(save=True):
    fig, ax = plt.subplots(1, 1, figsize=(4.2, 4.2))
    draw_mechanical_panel(ax)
    if save:
        for fn in ("figure_mechanical_element.png", "figure_mechanical_element.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def figure_coupling_flow(save=True):
    fig = plt.figure(figsize=(10, 3.6))
    gs = fig.add_gridspec(1, 5, width_ratios=[1.6, 0.6, 1.6, 0.6, 1.6])
    ax_spme = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs[0, 2])
    ax_mech = fig.add_subplot(gs[0, 4])
    for ax in (ax_spme, ax_mid, ax_mech):
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    # boxes
    ax_spme.add_patch(Rectangle((0.1, 0.2), 0.8, 0.6, fill=False, lw=1.8))
    ax_mid.add_patch(Rectangle((0.1, 0.2), 0.8, 0.6, fill=False, lw=1.8))
    ax_mech.add_patch(Rectangle((0.1, 0.2), 0.8, 0.6, fill=False, lw=1.8))
    ax_spme.text(0.5, 0.75, "SPMe (1D)", ha="center", weight="bold")
    ax_mid.text(0.5, 0.75, "2D Thermal (Q4)", ha="center", weight="bold")
    ax_mech.text(0.5, 0.75, "2D Mechanical (Q4)", ha="center", weight="bold")
    # labels inside
    ax_spme.text(0.5, 0.5, r"$q_{rxn}, q_{rev}, q_{ohm}$", ha="center", color="#2ca02c")
    ax_mid.text(0.5, 0.5, r"$T(x,y)$,  $T_e$", ha="center", color="#1f77b4")
    ax_mech.text(0.5, 0.5, r"$\varepsilon_{th}=\alpha\,\Delta T$", ha="center")
    # arrows
    ax_a1 = fig.add_subplot(gs[0, 1]); ax_a1.axis("off"); ax_a1.set_xlim(0,1); ax_a1.set_ylim(0,1)
    draw_coupling_arrows_bidirectional(ax_a1)
    ax_a1.text(0.5, 0.74, r"$q = q_{rxn}+q_{rev}+q_{ohm}$", ha="center", fontsize=9)
    ax_a1.text(0.5, 0.26, r"$T_e$", ha="center", fontsize=9)
    ax_a2 = fig.add_subplot(gs[0, 3]); ax_a2.axis("off"); ax_a2.set_xlim(0,1); ax_a2.set_ylim(0,1)
    draw_coupling_arrows_unidirectional(ax_a2)
    # footnote on aggregation/parallel split
    ax_mid.text(0.5, 0.28, "f_k 聚合, 并联分流 I_e", ha="center", fontsize=8, color="gray")
    fig.suptitle("Echem–Thermal–Mechanical coupling flow", y=0.98)
    if save:
        for fn in ("figure_coupling_flow.png", "figure_coupling_flow.svg"):
            fig.savefig(fn, dpi=300, bbox_inches="tight")
    return fig


def main():
    # Fonts: prefer Windows CJK-capable fonts, fallback to DejaVu Sans
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams.update({
        "font.size": 11,
        "figure.figsize": (15, 4.2),
        "axes.titlesize": 12,
    })

    # Combined figure (legacy)
    fig = plt.figure(figsize=(15, 4.2), constrained_layout=True)
    # Layout: [Jellyroll] - [space] - [SPMe] - [arrows] - [Thermal Q4] - [arrows] - [Mechanical Q4]
    gs = fig.add_gridspec(1, 7, width_ratios=[2.6, 0.2, 1.6, 0.25, 1.8, 0.25, 1.8])

    ax_left = fig.add_subplot(gs[0, 0])
    draw_cylinder_panel(ax_left)

    ax_echem = fig.add_subplot(gs[0, 2])
    draw_spme_panel(ax_echem)

    ax_th = fig.add_subplot(gs[0, 4])
    draw_thermal_panel(ax_th)

    ax_mech = fig.add_subplot(gs[0, 6])
    draw_mechanical_panel(ax_mech)

    # Coupling arrows: SPMe ↔ Thermal, Thermal → Mechanical
    ax_ar1 = fig.add_subplot(gs[0, 3])
    ax_ar1.axis("off")
    ax_ar1.set_xlim(0, 1)
    ax_ar1.set_ylim(-1, 1)
    draw_coupling_arrows_bidirectional(ax_ar1)

    ax_ar2 = fig.add_subplot(gs[0, 5])
    ax_ar2.axis("off")
    ax_ar2.set_xlim(0, 1)
    ax_ar2.set_ylim(-1, 1)
    draw_coupling_arrows_unidirectional(ax_ar2)

    fig.suptitle("SPMe (1D)  ↔  2D Thermal (Q4)  →  2D Mechanical (Q4)", y=1.02)

    for fn in ("figure_model_coupling.png", "figure_model_coupling.svg"):
        fig.savefig(fn, dpi=300, bbox_inches="tight")
    print("Saved:", "figure_model_coupling.png, figure_model_coupling.svg")

    # New: five separate figures
    figure_topview_thermal_mesh()
    figure_echem_1d_mesh()
    figure_thermal_element()
    figure_mechanical_element()
    figure_coupling_flow()
    figure_single_spiral_layered()
    print("Saved:", 
          "figure_mesh_topview.png, figure_echem_mesh.png, figure_thermal_element.png, figure_mechanical_element.png, figure_coupling_flow.png, figure_single_spiral_layered.png")


if __name__ == "__main__":
    main()