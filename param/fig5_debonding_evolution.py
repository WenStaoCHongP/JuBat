#!/usr/bin/env python3
"""Figure 5: debonding / delamination evolution at different cycle counts."""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def _damage_field_at_cycle(
    n: int,
    n_turns: float,
    cycle: int,
    r_inner: float,
    growth: float,
    repeat_thickness: float,
) -> tuple[np.ndarray, float]:
    """Generate a circular damage field for a given cycle number.

    Returns (damage_field, r_outer).
    """
    r_outer = r_inner + growth * 2.0 * np.pi * n_turns + repeat_thickness
    span = r_outer + 0.08  # small padding

    x = np.linspace(-span, span, n)
    y = np.linspace(-span, span, n)
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    tt = np.arctan2(yy, xx)

    r_norm = np.clip((rr - r_inner) / (r_outer - r_inner), 0.0, 1.0)

    if cycle == 0:
        damage = np.zeros_like(rr)
    else:
        rate = 1.0 - np.exp(-cycle / 120.0)
        radial_profile = np.exp(-3.5 * r_norm)
        angular_mod = 1.0 + 0.35 * np.cos(2.0 * tt - 0.3) + 0.2 * np.sin(3.0 * tt + 0.7)
        damage = np.clip(rate * radial_profile * angular_mod, 0.0, 1.0)

    mask = (rr < r_inner) | (rr > r_outer)
    damage[mask] = np.nan
    return damage, span


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    cycles = [0, 10, 50, 150, 350]
    n_turns = 3.5
    r_inner = 0.30
    growth = 0.22 / (2.0 * np.pi)
    repeat_thickness = 0.22
    r_outer = r_inner + growth * 2.0 * np.pi * n_turns + repeat_thickness
    n_px = 500

    cmap_damage = LinearSegmentedColormap.from_list(
        "debond",
        [
            (0.00, "#e0e0e0"),
            (0.15, "#fff7bc"),
            (0.35, "#fec44f"),
            (0.55, "#f16913"),
            (0.75, "#d62728"),
            (1.00, "#67000d"),
        ],
        N=512,
    )
    cmap_damage.set_bad(color="white", alpha=0)

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 5, figsize=(7.4, 2.4), dpi=300, constrained_layout=True)

    for ax, cyc in zip(axes, cycles):
        damage, span = _damage_field_at_cycle(
            n=n_px, n_turns=n_turns, cycle=cyc,
            r_inner=r_inner, growth=growth,
            repeat_thickness=repeat_thickness,
        )

        ax.imshow(
            damage,
            extent=[-span, span, -span, span],
            origin="lower",
            cmap=cmap_damage,
            vmin=0.0,
            vmax=1.0,
            interpolation="bilinear",
        )

        circle_inner = plt.Circle((0, 0), r_inner, fill=False, ec="#444444", lw=0.6)
        circle_outer = plt.Circle((0, 0), r_outer, fill=False, ec="#444444", lw=0.6)
        ax.add_patch(circle_inner)
        ax.add_patch(circle_outer)

        ax.set_aspect("equal")
        margin = span + 0.02
        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.axis("off")
        ax.set_title(f"N = {cyc}", fontsize=9, pad=3)

    sm = plt.cm.ScalarMappable(cmap=cmap_damage, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.012, pad=0.015, aspect=25)
    cbar.set_label("Damage D", fontsize=8)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.tick_params(labelsize=7)

    out_path = os.path.join(output_dir, "fig5_debonding_evolution.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
