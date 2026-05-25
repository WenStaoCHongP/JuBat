#!/usr/bin/env python3
"""Figure 9: triple-tab debonding evolution (tabs at 0, 1/2, 1 of spiral).

CZM damage lives on the interface (the Archimedean spiral).  The spiral is
rendered as a coloured line whose colour encodes D ∈ [0, 1].

Three tabs: start (θ=0), middle (θ=θ_total/2), end (θ=θ_total).
Damage hierarchy: mid-tab > start-tab >> end-tab.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


def _spiral_xy(theta: float | np.ndarray, a: float, b: float):
    """Cartesian coordinates on the Archimedean spiral r = a + bθ."""
    r = a + b * theta
    return r * np.cos(theta), r * np.sin(theta)


def _spiral_damage_profile_triple(
    theta: np.ndarray,
    cycle: int,
    tab_thetas: tuple[float, float, float],
) -> np.ndarray:
    """Damage D(θ) with three tabs."""
    if cycle == 0:
        return np.zeros_like(theta)

    total_angle = theta[-1]
    sigma = np.pi * 0.25

    intensities = [0.55, 1.0, 0.30]
    delays = [0.45, 0.0, 0.70]

    damage = np.zeros_like(theta)
    for tab_th, intensity, delay in zip(tab_thetas, intensities, delays):
        d_th = np.abs(theta - tab_th)
        d_th = np.minimum(d_th, total_angle - d_th)
        hotspot = np.exp(-d_th**2 / (2.0 * sigma**2))
        ramp = 1.0 / (1.0 + np.exp(0.06 * (delay * 350 - cycle)))
        damage += intensity * ramp * hotspot

    rate = 1.0 - np.exp(-cycle / 200.0)
    background = 0.01
    damage = np.clip(rate * (damage + background), 0.0, 1.0)
    return damage


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    n_turns = 3.5
    r_inner = 0.30
    growth = 0.22 / (2.0 * np.pi)
    repeat_thickness = 0.22
    r_outer = r_inner + growth * 2.0 * np.pi * n_turns + repeat_thickness
    a, b = r_inner, growth

    total_angle = 2.0 * np.pi * n_turns
    tab_thetas = (0.0, total_angle * 0.5, total_angle)

    theta_spiral = np.linspace(0.0, total_angle, 3000)
    x_sp, y_sp = _spiral_xy(theta_spiral, a, b)

    cycles = [0, 10, 50, 150, 350]

    cmap_damage = LinearSegmentedColormap.from_list(
        "debond",
        [
            (0.00, "#ffffff"),
            (0.15, "#fff7bc"),
            (0.35, "#fec44f"),
            (0.55, "#f16913"),
            (0.75, "#d62728"),
            (1.00, "#67000d"),
        ],
        N=512,
    )
    norm_damage = Normalize(vmin=0.0, vmax=1.0)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(
        1, 5, figsize=(7.4, 2.4), dpi=300, constrained_layout=True
    )

    span = r_outer + 0.12
    tab_positions = [_spiral_xy(th, a, b) for th in tab_thetas]
    tab_labels = ["Tab 0 (start)", r"Tab 1/2 (mid)", "Tab 1 (end)"]
    tab_colors = ["#1f77b4", "#d62728", "#2ca02c"]

    for col, cyc in enumerate(cycles):
        ax = axes[col]

        # Annular fill
        th_fill = np.linspace(0, 2 * np.pi, 300)
        ax.fill(r_outer * np.cos(th_fill), r_outer * np.sin(th_fill), color="white", zorder=0)
        ax.fill(r_inner * np.cos(th_fill), r_inner * np.sin(th_fill), color="white", zorder=1)

        # Boundary circles
        ax.add_patch(plt.Circle((0, 0), r_inner, fill=False, ec="#999999", lw=0.4, zorder=3))
        ax.add_patch(plt.Circle((0, 0), r_outer, fill=False, ec="#666666", lw=0.6, zorder=3))

        # Spiral damage line
        damage_profile = _spiral_damage_profile_triple(theta_spiral, cyc, tab_thetas)

        points = np.column_stack([x_sp, y_sp]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        seg_damage = 0.5 * (damage_profile[:-1] + damage_profile[1:])

        lc = LineCollection(segments, cmap=cmap_damage, norm=norm_damage, linewidths=0.8, zorder=5)
        lc.set_array(seg_damage)
        ax.add_collection(lc)

        # Tab markers
        for (xt, yt), color in zip(tab_positions, tab_colors):
            ax.plot(
                xt, yt, "D", markersize=4,
                markerfacecolor=color, markeredgecolor="white",
                markeredgewidth=0.5, zorder=10,
            )

        ax.set_aspect("equal")
        margin = span + 0.02
        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.set_facecolor("white")
        ax.axis("off")
        ax.set_title(f"N = {cyc}", fontsize=9, pad=3)

    # Colour bar
    sm_d = plt.cm.ScalarMappable(cmap=cmap_damage, norm=norm_damage)
    sm_d.set_array([])
    cbar_d = fig.colorbar(
        sm_d, ax=axes.ravel().tolist(),
        fraction=0.015, pad=0.015, aspect=25,
    )
    cbar_d.set_label("Damage D", fontsize=8)
    cbar_d.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar_d.ax.tick_params(labelsize=7)

    # Legend
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor=c,
               markersize=5, label=l, markeredgecolor="white", markeredgewidth=0.4)
        for c, l in zip(tab_colors, tab_labels)
    ]
    fig.legend(
        handles=legend_handles, loc="lower center", ncol=3,
        fontsize=6.5, frameon=True, framealpha=0.85,
        prop={"family": "Times New Roman"},
        bbox_to_anchor=(0.5, -0.02),
    )

    out_path = os.path.join(output_dir, "fig9_triple_tab_debonding.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
