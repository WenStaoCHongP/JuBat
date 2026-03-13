#!/usr/bin/env python3
"""
Paper figures for JuBat: Jellyroll battery SPMe-thermal-CZM coupled model.

Outputs:
- figure_1_mesh_discretization.png/svg (Mesh discretization diagram)
- figure_2_coupling_schematic.png/svg (Coupling schematic)

Run:
    python plot.py
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import os

# Output directory: photo folder under workspace root
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "photo")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Journal-compliant settings (Elsevier double-column width: 7.2 in)
FIG_WIDTH_DOUBLE = 7.2  # inches
DPI = 300

# Color scheme (colorblind-friendly)
COLORS = {
    "spme": "#2ca02c",      # green - electrochemical
    "thermal": "#1f77b4",   # blue - thermal
    "czm": "#d62728",       # red - mechanical
    "coupling": "#7f7f7f",  # gray - coupling
    "branch": "#9467bd",    # purple - branch solver
    "itr": "#ff7f0e",       # orange - interface thermal resistance
}

# Layer colors for mesh
LAYER_COLORS = ("#8c564b", "#ff7f0e", "#7f7f7f", "#1f77b4", "#5B5B5B")
LAYER_NAMES = ("PCC", "PE", "SP", "NE", "NCC")
LAYER_FRACS = (0.12, 0.30, 0.16, 0.30, 0.12)


def _polar_point(r: float, theta: float, center=(0.0, 0.0)):
    """Utility: polar to Cartesian on given center."""
    return (center[0] + r * np.cos(theta), center[1] + r * np.sin(theta))


def _arch_spiral(a: float, b: float, theta: np.ndarray, center=(0, 0)):
    """Archimedean spiral: r = a + b*theta"""
    r = a + b * theta
    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    return x, y, r


def _draw_collector_seeded_band_mesh(
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
    fill_layers=False,
    fill_alpha=0.35,
    layer_colors=LAYER_COLORS,
    layer_fracs=LAYER_FRACS,
    overlay_layers_band_index=None,
    overlay_tick_stride=2,
    overlay_tick_lw=2.4,
    only_band_index=None,
):
    """
    Draw collector-seeded Q4 band mesh for jellyroll battery.
    
    Based on Jellyrollmodel.jl:
    - Archimedean spiral: r(θ) = a + bθ
    - Q4 elements connecting inner and outer spiral rails
    """
    theta_max = 2 * np.pi * turns
    a = r_in
    t_repeat = (r_out - r_in) / max(1, nbands)
    b = t_repeat / (2.0 * np.pi)
    dr = (r_out - r_in) / nbands
    dtheta = 2 * np.pi / float(max(1, seg_per_turn))

    poly_verts = []
    poly_fc = []

    for k in range(nbands):
        if (only_band_index is not None) and (k != int(only_band_index)):
            continue

        s_in = k * dr
        s_out = (k + 1) * dr

        theta0 = max(0.0, (r_in - a - s_in) / b)
        theta1_lim = (r_out - a - s_out) / b
        if theta1_lim <= theta0:
            continue
        theta1 = min(theta_max, theta1_lim)

        k0 = int(np.ceil((theta0 - phase) / dtheta))
        k1 = int(np.floor((theta1 - phase) / dtheta))
        if k1 <= k0:
            theta = np.array([theta0, min(theta1, theta0 + dtheta*0.75)])
        else:
            theta = phase + dtheta * np.arange(k0, k1 + 1)

        xin, yin, _ = _arch_spiral(a + s_in, b, theta, center)
        xout, yout, _ = _arch_spiral(a + s_out, b, theta, center)

        # Draw spiral rails
        if only_band_index is not None:
            ax.plot(xin, yin, color=rail_color, lw=rail_lw, alpha=0.9)
            ax.plot(xout, yout, color=rail_color, lw=rail_lw, alpha=0.9)
        else:
            if k == 0:
                ax.plot(xin, yin, color=rail_color, lw=rail_lw, alpha=0.9)
            ax.plot(xout, yout, color=rail_color, lw=rail_lw, alpha=0.9)

        # Generate Q4 elements
        for i in range(len(theta) - 1):
            ax.plot([xin[i], xout[i]], [yin[i], yout[i]], color=edge_color, lw=edge_lw)
            ax.plot([xin[i + 1], xout[i + 1]], [yin[i + 1], yout[i + 1]], color=edge_color, lw=edge_lw)
            ax.plot([xout[i], xout[i + 1]], [yout[i], yout[i + 1]], color=edge_color, lw=edge_lw)
            ax.plot([xin[i], xin[i + 1]], [yin[i], yin[i + 1]], color=edge_color, lw=edge_lw)

            # Layer filling
            if fill_layers:
                fracs = np.asarray(layer_fracs, dtype=float)
                fracs = fracs / fracs.sum()
                cum = np.concatenate(([0.0], np.cumsum(fracs)))
                th_i, th_ip1 = theta[i], theta[i + 1]
                for j in range(len(fracs)):
                    s_lo = s_in + cum[j] * dr
                    s_hi = s_in + cum[j + 1] * dr
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

        # Sparse node markers
        if node_stride > 0 and len(theta) > 0:
            sel = np.arange(0, len(theta), max(1, node_stride))
            ax.scatter(xin[sel], yin[sel], s=8, color="black", zorder=3)
            ax.scatter(xout[sel], yout[sel], s=8, color="black", zorder=3)

    # Add filled polygons
    if fill_layers and len(poly_verts) > 0:
        from matplotlib.collections import PolyCollection
        coll = PolyCollection(poly_verts, facecolors=poly_fc, edgecolors='none', alpha=fill_alpha, zorder=1)
        ax.add_collection(coll)


def figure_1_mesh_discretization(save=True):
    """
    Figure 1: Mesh discretization diagram for jellyroll battery.
    
    Shows:
    - 3-turn spiral mesh with layer color filling
    - Thermal boundary conditions (convection outer, adiabatic inner)
    - Tab positions with bold markers
    - Layer legend
    """
    fig, ax = plt.subplots(1, 1, figsize=(FIG_WIDTH_DOUBLE, FIG_WIDTH_DOUBLE * 0.85))
    
    r_in, r_out, turns = 0.25, 1.9, 3.0
    nbands, seg_per_turn = 6, 24
    R_outer = 2.3  # Visual outer boundary
    
    # Outer boundary circle
    ax.add_patch(Circle((0, 0), R_outer, fill=False, ec="black", lw=1.5))
    
    # Draw spiral mesh with layer filling
    _draw_collector_seeded_band_mesh(
        ax,
        r_in=r_in, r_out=r_out, turns=turns, nbands=nbands,
        seg_per_turn=seg_per_turn, phase=0.0,
        rail_lw=1.0, edge_lw=0.6, rail_color="#B05300", edge_color="#555555",
        node_stride=8, fill_layers=True, fill_alpha=0.5,
        layer_colors=LAYER_COLORS, layer_fracs=LAYER_FRACS,
        overlay_layers_band_index=max(0, nbands//2 - 1),
    )
    
    # Coordinate axes
    ax.arrow(0, 0, 0.8, 0.0, head_width=0.06, head_length=0.1, fc="k", ec="k")
    ax.arrow(0, 0, 0.0, 0.8, head_width=0.06, head_length=0.1, fc="k", ec="k")
    ax.text(0.95, 0.0, "x", fontsize=10, fontname="Times New Roman")
    ax.text(0.0, 0.95, "y", fontsize=10, fontname="Times New Roman")
    
    # === Thermal boundary conditions ===
    # Outer convection arrows
    conv_thetas = np.deg2rad(np.array([-60, -25, 10, 45, 80], dtype=float))
    L_arrow = 0.25
    for th in conv_thetas:
        x0, y0 = _polar_point(R_outer, th)
        x1, y1 = _polar_point(R_outer + L_arrow, th)
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                   arrowprops=dict(arrowstyle="-|>", lw=1.4, color=COLORS["thermal"]))
    ax.text(R_outer + L_arrow + 0.15, 0.0,
            r"Convection: $-k\,\partial T/\partial n = h(T-T_{amb})$",
            fontsize=8, color=COLORS["thermal"], va="center", fontname="Times New Roman")
    
    # Inner adiabatic ticks
    adiab_thetas = np.deg2rad(np.array([0, 90, 180, 270]))
    tL = 0.1
    for th in adiab_thetas:
        ux, uy = np.cos(th), np.sin(th)
        x0, y0 = (r_in - tL) * ux, (r_in - tL) * uy
        x1, y1 = (r_in + tL) * ux, (r_in + tL) * uy
        ax.plot([x0, x1], [y0, y1], color="#2ca02c", lw=2.0, zorder=5)
    ax.text(-R_outer - 0.1, 0.0, r"Adiabatic: $\partial T/\partial n = 0$",
            fontsize=8, color="#2ca02c", ha="right", va="center", fontname="Times New Roman")
    
    # Inner boundary circle (dashed)
    ax.add_patch(Circle((0, 0), r_in, fill=False, ec="#7f7f7f", lw=1.0, ls=(0, (4, 3)), alpha=0.9))
    
    # === Tab positions (bold markers) ===
    # Positive tab (near θ=0, inner spiral start)
    pos_tab_theta = np.deg2rad(10)
    pos_tab_r = r_in + 0.15
    px, py = _polar_point(pos_tab_r, pos_tab_theta)
    ax.plot(px, py, "s", markersize=12, markerfacecolor="none", markeredgecolor="#d62728", markeredgewidth=2.5, zorder=10)
    ax.text(px + 0.15, py + 0.15, "+Tab", fontsize=9, color="#d62728", fontweight="bold", fontname="Times New Roman")
    
    # Negative tab (near outer boundary)
    neg_tab_theta = np.deg2rad(-170)
    neg_tab_r = r_out - 0.15
    nx, ny = _polar_point(neg_tab_r, neg_tab_theta)
    ax.plot(nx, ny, "s", markersize=12, markerfacecolor="none", markeredgecolor="#1f77b4", markeredgewidth=2.5, zorder=10)
    ax.text(nx - 0.35, ny - 0.15, "-Tab", fontsize=9, color="#1f77b4", fontweight="bold", fontname="Times New Roman")
    
    # Layer legend
    handles = [Line2D([0], [0], color=c, lw=2.5, label=n) for n, c in zip(LAYER_NAMES, LAYER_COLORS)]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(-0.02, 1.02),
              fontsize=8, frameon=True, framealpha=0.9, ncol=5, prop={'family': 'Times New Roman'})
    
    ax.set_aspect("equal")
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-2.6, 2.6)
    ax.axis("off")
    
    if save:
        for ext in ("png", "svg"):
            fn = os.path.join(OUTPUT_DIR, f"figure_1_mesh_discretization.{ext}")
            fig.savefig(fn, dpi=DPI, bbox_inches="tight", facecolor="white")
        print(f"Saved: figure_1_mesh_discretization.png/svg")
    
    return fig


def _draw_module_box(ax, xy, width, height, label, color, sublabel=None, fontsize=9):
    """Draw a module box with label and optional sublabel."""
    x, y = xy
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.05",
                         fill=False, ec=color, lw=1.8, zorder=2)
    ax.add_patch(box)
    ax.text(x, y + height/2 - 0.12, label, ha="center", va="top",
            fontsize=fontsize, fontweight="bold", color=color, fontname="Times New Roman")
    if sublabel:
        ax.text(x, y - 0.05, sublabel, ha="center", va="center",
                fontsize=7, color="#555555", fontname="Times New Roman")


def _draw_arrow(ax, start, end, color="gray", style="-|>", lw=1.2, curved=False, label=None, label_pos=0.5):
    """Draw an arrow between two points with optional label."""
    if curved:
        arrow = FancyArrowPatch(start, end, connectionstyle="arc3,rad=0.2",
                               arrowstyle=style, mutation_scale=10,
                               lw=lw, fc=color, ec=color, alpha=0.85, zorder=1)
    else:
        arrow = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=10,
                               lw=lw, fc=color, ec=color, alpha=0.85, zorder=1)
    ax.add_patch(arrow)
    if label:
        mx = start[0] + label_pos * (end[0] - start[0])
        my = start[1] + label_pos * (end[1] - start[1])
        ax.text(mx, my + 0.08, label, ha="center", va="bottom",
                fontsize=7, color=color, fontname="Times New Roman")


def _draw_bidirectional_arrow(ax, start, end, color="gray", lw=1.2, label_top=None, label_bot=None):
    """Draw bidirectional arrows with labels."""
    # Top arrow (forward)
    arrow1 = FancyArrowPatch(start, end, connectionstyle="arc3,rad=0.15",
                            arrowstyle="-|>", mutation_scale=10,
                            lw=lw, fc=color, ec=color, alpha=0.85, zorder=1)
    ax.add_patch(arrow1)
    # Bottom arrow (backward)
    arrow2 = FancyArrowPatch(end, start, connectionstyle="arc3,rad=0.15",
                            arrowstyle="-|>", mutation_scale=10,
                            lw=lw, fc=color, ec=color, alpha=0.85, zorder=1)
    ax.add_patch(arrow2)
    if label_top:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2 + 0.18
        ax.text(mx, my, label_top, ha="center", va="bottom",
                fontsize=7, color=color, fontname="Times New Roman")
    if label_bot:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2 - 0.18
        ax.text(mx, my, label_bot, ha="center", va="top",
                fontsize=7, color=color, fontname="Times New Roman")


def figure_2_coupling_schematic(save=True):
    """
    Figure 2: SPMe-Thermal-CZM coupling schematic.
    
    Shows complete coupling relationships:
    - SPMe: I_e, T_e → q, σ_diff
    - Thermal: q, h_eff → T, σ_th
    - CZM: σ_diff + σ_th → D, δ
    - Branch Solver: I_total, T_e, D → I_e, V
    - Interface Thermal Resistance: D, δ → h_eff
    """
    fig, ax = plt.subplots(1, 1, figsize=(FIG_WIDTH_DOUBLE, FIG_WIDTH_DOUBLE * 0.7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    
    # === Module positions ===
    # Row 1 (top): Branch Solver
    bs_pos = (5.0, 5.2)
    
    # Row 2 (middle): SPMe, Thermal, ITR
    spme_pos = (1.5, 3.2)
    thermal_pos = (5.0, 3.2)
    itr_pos = (8.5, 3.2)
    
    # Row 3 (bottom): CZM
    czm_pos = (5.0, 1.0)
    
    box_w, box_h = 1.8, 1.0
    
    # === Draw modules ===
    # Branch Solver (top center)
    _draw_module_box(ax, bs_pos, box_w, box_h, "Branch Solver", COLORS["branch"],
                    "Newton-Raphson\n$I_{total}$, $T_e$, $D$ → $I_e$, $V$")
    
    # SPMe (middle left)
    _draw_module_box(ax, spme_pos, box_w, box_h, "SPMe (per elem)", COLORS["spme"],
                    "$I_e$, $T_e$ → $q$, $\\sigma_{diff}$")
    
    # Thermal (middle center)
    _draw_module_box(ax, thermal_pos, box_w, box_h, "2D Thermal", COLORS["thermal"],
                    "$q$, $h_{eff}$ → $T$, $\\sigma_{th}$")
    
    # ITR (middle right)
    _draw_module_box(ax, itr_pos, box_w * 0.9, box_h, "Interface\nThermal Res.", COLORS["itr"],
                    "$D$, $\\delta$ → $h_{eff}$")
    
    # CZM (bottom center)
    _draw_module_box(ax, czm_pos, box_w, box_h, "CZM (COH2D4)", COLORS["czm"],
                    "$\\sigma = \\sigma_{diff} + \\sigma_{th}$\n→ $D$, $\\delta$")
    
    # === Draw coupling arrows ===
    # Branch Solver ↔ SPMe (bidirectional: I_e down, T_e up)
    _draw_bidirectional_arrow(ax, (bs_pos[0] - 0.5, bs_pos[1] - box_h/2),
                              (spme_pos[0] + 0.3, spme_pos[1] + box_h/2),
                              color=COLORS["branch"], label_top="$I_e$", label_bot="$T_e$")
    
    # SPMe ↔ Thermal (bidirectional: q right, T_e left)
    _draw_bidirectional_arrow(ax, (spme_pos[0] + box_w/2, spme_pos[1]),
                              (thermal_pos[0] - box_w/2, thermal_pos[1]),
                              color=COLORS["coupling"], label_top="$q$", label_bot="$T_e$")
    
    # Thermal → CZM (σ_th down)
    _draw_arrow(ax, (thermal_pos[0] - 0.3, thermal_pos[1] - box_h/2),
               (czm_pos[0] - 0.3, czm_pos[1] + box_h/2),
               color=COLORS["thermal"], label="$\\sigma_{th}$")
    
    # SPMe → CZM (σ_diff diagonal)
    _draw_arrow(ax, (spme_pos[0] + box_w/2, spme_pos[1] - box_h/2),
               (czm_pos[0] - box_w/2, czm_pos[1] + box_h/2),
               color=COLORS["spme"], curved=True, label="$\\sigma_{diff}$")
    
    # CZM → ITR (D, δ right-up)
    _draw_arrow(ax, (czm_pos[0] + box_w/2, czm_pos[1] + box_h/3),
               (itr_pos[0] - box_w*0.45, itr_pos[1] - box_h/3),
               color=COLORS["czm"], label="$D$, $\\delta$")
    
    # ITR → Thermal (h_eff left)
    _draw_arrow(ax, (itr_pos[0] - box_w*0.45, itr_pos[1]),
               (thermal_pos[0] + box_w/2, thermal_pos[1]),
               color=COLORS["itr"], label="$h_{eff}$")
    
    # CZM → Branch Solver (D up, dashed for feedback)
    arrow_d = FancyArrowPatch((czm_pos[0] + 0.5, czm_pos[1] + box_h/2),
                              (bs_pos[0] + 0.5, bs_pos[1] - box_h/2),
                              arrowstyle="-|>", mutation_scale=10,
                              lw=1.2, fc=COLORS["czm"], ec=COLORS["czm"],
                              alpha=0.85, linestyle="--", zorder=1)
    ax.add_patch(arrow_d)
    ax.text(czm_pos[0] + 0.8, (czm_pos[1] + bs_pos[1])/2, "$D$ (deactivated)",
            fontsize=7, color=COLORS["czm"], fontname="Times New Roman")
    
    # === Legend for coupling types ===
    legend_items = [
        ("Bidirectional coupling", COLORS["coupling"], "-"),
        ("Unidirectional coupling", COLORS["coupling"], "-"),
        ("Feedback (CZM failure)", COLORS["czm"], "--"),
    ]
    legend_y = 0.3
    for i, (text, color, ls) in enumerate(legend_items):
        x = 0.5 + i * 3.2
        ax.plot([x, x + 0.5], [legend_y, legend_y], color=color, lw=1.5, linestyle=ls)
        ax.text(x + 0.6, legend_y, text, fontsize=7, va="center", fontname="Times New Roman")
    
    # Title
    ax.text(5.0, 5.9, "Multi-physics Coupling Schematic", ha="center", va="bottom",
            fontsize=11, fontweight="bold", fontname="Times New Roman")
    
    if save:
        for ext in ("png", "svg"):
            fn = os.path.join(OUTPUT_DIR, f"figure_2_coupling_schematic.{ext}")
            fig.savefig(fn, dpi=DPI, bbox_inches="tight", facecolor="white")
        print(f"Saved: figure_2_coupling_schematic.png/svg")
    
    return fig


def main():
    """Generate both paper figures."""
    # Journal-compliant font settings
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"]
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9
    
    print("Generating paper figures...")
    print(f"Output directory: {OUTPUT_DIR}")
    
    figure_1_mesh_discretization(save=True)
    figure_2_coupling_schematic(save=True)
    
    print("\nDone. Generated:")
    print("  - figure_1_mesh_discretization.png/svg")
    print("  - figure_2_coupling_schematic.png/svg")


if __name__ == "__main__":
    main()