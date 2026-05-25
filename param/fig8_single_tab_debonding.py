#!/usr/bin/env python3
"""Figure 8: single-tab debonding evolution with tab at 1/3 and 1/2 of spiral.

CZM is an *interface* model – damage lives on the cohesive surface, i.e. along
the Archimedean spiral.  Each sub-plot renders the spiral as a coloured line
whose colour encodes the local damage value D ∈ [0, 1].

Two rows compare tab positions at 1/3 and 1/2 of the total spiral angle.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


# ── Spiral helpers ───────────────────────────────────────────────────────

def _spiral_xy(theta: float | np.ndarray, a: float, b: float):
    """Cartesian coordinates on the Archimedean spiral r = a + bθ."""
    r = a + b * theta
    return r * np.cos(theta), r * np.sin(theta)


# ── Damage along the spiral ──────────────────────────────────────────────

def _spiral_damage_profile(
    theta: np.ndarray,
    cycle: int,
    tab_theta: float,
) -> np.ndarray:
    """Return damage D(θ) for every point on the spiral."""
    if cycle == 0:
        return np.zeros_like(theta)

    rate = 1.0 - np.exp(-cycle / 200.0)
    d_theta = np.abs(theta - tab_theta)
    total_angle = theta[-1]
    d_theta = np.minimum(d_theta, total_angle - d_theta)

    # Narrow hotspot: D=1 only right at the tab after ~350 cycles
    sigma = np.pi * 0.25
    hotspot = np.exp(-d_theta**2 / (2.0 * sigma**2))

    opp_theta = tab_theta + np.pi
    d_opp = np.abs(theta - opp_theta)
    d_opp = np.minimum(d_opp, total_angle - d_opp)
    secondary = 0.05 * np.exp(-d_opp**2 / (2.0 * (sigma * 0.8)**2))

    background = 0.01
    damage = np.clip(rate * (hotspot + secondary + background), 0.0, 1.0)
    return damage


# ── Figure entry point ───────────────────────────────────────────────────

def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Spiral geometry parameters ──
    n_turns = 3.5
    r_inner = 0.30
    growth = 0.22 / (2.0 * np.pi)
    repeat_thickness = 0.22
    r_outer = r_inner + growth * 2.0 * np.pi * n_turns + repeat_thickness
    a, b = r_inner, growth

    total_angle = 2.0 * np.pi * n_turns

    theta_spiral = np.linspace(0.0, total_angle, 3000)
    x_sp, y_sp = _spiral_xy(theta_spiral, a, b)

    tab_theta_third = total_angle * (1.0 / 3.0)
    tab_theta_half = total_angle * (1.0 / 2.0)

    cycles = [0, 10, 50, 150, 350]

    # ── Colourmap for damage ──
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
        2, 5, figsize=(7.4, 4.2), dpi=300, constrained_layout=True
    )

    tab_configs = [
        (tab_theta_third, r"Tab at 1/3 spiral"),
        (tab_theta_half, r"Tab at 1/2 spiral"),
    ]

    span = r_outer + 0.12

    for row, (tab_theta, row_title) in enumerate(tab_configs):
        x_tab, y_tab = _spiral_xy(tab_theta, a, b)

        for col, cyc in enumerate(cycles):
            ax = axes[row, col]

            # ── Annular fill (light grey for context) ──
            th_fill = np.linspace(0, 2 * np.pi, 300)
            ax.fill(
                r_outer * np.cos(th_fill),
                r_outer * np.sin(th_fill),
                color="white", zorder=0,
            )
            ax.fill(
                r_inner * np.cos(th_fill),
                r_inner * np.sin(th_fill),
                color="white", zorder=1,
            )

            # ── Boundary circles ──
            ax.add_patch(plt.Circle((0, 0), r_inner, fill=False, ec="#999999", lw=0.4, zorder=3))
            ax.add_patch(plt.Circle((0, 0), r_outer, fill=False, ec="#666666", lw=0.6, zorder=3))

            # ── Spiral damage line ──
            damage_profile = _spiral_damage_profile(theta_spiral, cyc, tab_theta)

            points = np.column_stack([x_sp, y_sp]).reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            seg_damage = 0.5 * (damage_profile[:-1] + damage_profile[1:])

            lc = LineCollection(
                segments, cmap=cmap_damage, norm=norm_damage,
                linewidths=0.8, zorder=5,
            )
            lc.set_array(seg_damage)
            ax.add_collection(lc)

            # Tab marker
            ax.plot(
                x_tab, y_tab, "D",
                markersize=4.0,
                markerfacecolor="#d62728",
                markeredgecolor="white",
                markeredgewidth=0.5,
                zorder=10,
            )

            ax.set_aspect("equal")
            margin = span + 0.02
            ax.set_xlim(-margin, margin)
            ax.set_ylim(-margin, margin)
            ax.set_facecolor("white")
            ax.axis("off")

            if row == 0:
                ax.set_title(f"N = {cyc}", fontsize=9, pad=3)

        # Row label
        axes[row, 0].text(
            -0.15, 0.5, row_title,
            transform=axes[row, 0].transAxes,
            fontsize=7.5, va="center", ha="right", rotation=90,
            fontname="Times New Roman",
        )

    # ── Colour bar ──
    sm_d = plt.cm.ScalarMappable(cmap=cmap_damage, norm=norm_damage)
    sm_d.set_array([])
    cbar_d = fig.colorbar(
        sm_d, ax=axes.ravel().tolist(),
        fraction=0.012, pad=0.015, aspect=30,
    )
    cbar_d.set_label("Damage D", fontsize=8)
    cbar_d.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar_d.ax.tick_params(labelsize=7)

    out_path = os.path.join(output_dir, "fig8_single_tab_debonding.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
