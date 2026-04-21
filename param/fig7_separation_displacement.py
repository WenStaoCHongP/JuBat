#!/usr/bin/env python3
"""Figure 7: interface separation / gap displacement field at different cycles."""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch


def _separation_field(
    n: int,
    n_turns: float,
    cycle: int,
    r_inner: float,
    growth: float,
    repeat_thickness: float,
    n_interfaces: int = 4,
) -> tuple[np.ndarray, float]:
    """Generate separation displacement field on jellyroll annulus.

    The separation δ is concentrated along discrete interface bands
    (modelling cohesive zone locations) and grows with cycle count.
    """
    r_outer = r_inner + growth * 2.0 * np.pi * n_turns + repeat_thickness
    span = r_outer + 0.08

    x = np.linspace(-span, span, n)
    y = np.linspace(-span, span, n)
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    tt = np.arctan2(yy, xx)

    r_norm = np.clip((rr - r_inner) / (r_outer - r_inner), 0.0, 1.0)

    if cycle == 0:
        sep = np.zeros_like(rr)
    else:
        # Overall separation growth
        rate = 1.0 - np.exp(-cycle / 100.0)

        # Radial profile: inner interfaces open first
        radial = np.exp(-3.0 * r_norm)

        # Discrete interface bands (Gaussian peaks at interface locations)
        interface_signal = np.zeros_like(rr)
        for k in range(n_interfaces):
            r_k = (k + 1.0) / (n_interfaces + 1.0)
            sigma = 0.04 + 0.02 * k / n_interfaces
            interface_signal += np.exp(-((r_norm - r_k) ** 2) / (2.0 * sigma**2))

        # Angular modulation: non-uniform opening
        angular = 1.0 + 0.4 * np.cos(2.0 * tt - 0.5) + 0.25 * np.sin(3.0 * tt + 1.0)

        sep = rate * radial * interface_signal * angular
        sep = np.clip(sep, 0.0, None)

    mask = (rr < r_inner) | (rr > r_outer)
    sep[mask] = np.nan
    return sep, span


def _draw_interface_lines(ax, r_inner, r_outer, n_interfaces, growth, n_turns, n_pts=300):
    """Draw thin arcs at each interface location."""
    for k in range(n_interfaces):
        r_k = r_inner + (k + 1.0) / (n_interfaces + 1.0) * (r_outer - r_inner)
        theta = np.linspace(0, 2.0 * np.pi, n_pts)
        # Spiral: slight angular offset with radius
        r_line = r_k + growth * theta * 0.15
        x = r_line * np.cos(theta)
        y = r_line * np.sin(theta)
        ax.plot(x, y, color="#333333", lw=0.3, alpha=0.35, zorder=3)


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
    n_interfaces = 4

    # Colormap: white (no separation) → light blue → cyan → purple → dark purple
    cmap_sep = LinearSegmentedColormap.from_list(
        "separation",
        [
            (0.00, "#f7f7f7"),   # no separation: near white
            (0.10, "#d4e6f1"),   # barely open: pale blue
            (0.25, "#85c1e9"),   # slight: sky blue
            (0.45, "#2e86c1"),   # moderate: blue
            (0.65, "#8e44ad"),   # significant: purple
            (0.85, "#6c3483"),   # large: dark purple
            (1.00, "#1a0a2e"),   # fully open: near black
        ],
        N=512,
    )
    cmap_sep.set_bad(color="white", alpha=0)

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 5, figsize=(7.4, 2.4), dpi=300, constrained_layout=True)

    for ax, cyc in zip(axes, cycles):
        sep, span = _separation_field(
            n=n_px, n_turns=n_turns, cycle=cyc,
            r_inner=r_inner, growth=growth,
            repeat_thickness=repeat_thickness,
            n_interfaces=n_interfaces,
        )

        ax.imshow(
            sep,
            extent=[-span, span, -span, span],
            origin="lower",
            cmap=cmap_sep,
            vmin=0.0,
            vmax=1.0,
            interpolation="bilinear",
        )

        # Interface lines
        _draw_interface_lines(ax, r_inner, r_outer, n_interfaces, growth, n_turns)

        # Boundary circles
        ax.add_patch(plt.Circle((0, 0), r_inner, fill=False, ec="#444444", lw=0.6))
        ax.add_patch(plt.Circle((0, 0), r_outer, fill=False, ec="#444444", lw=0.6))

        ax.set_aspect("equal")
        margin = span + 0.02
        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.axis("off")
        ax.set_title(f"N = {cyc}", fontsize=9, pad=3)

    # Shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_sep, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.012, pad=0.015, aspect=25)
    cbar.set_label(r"Separation $\delta$ (nm)", fontsize=8)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0", "50", "100", "150", "200"])
    cbar.ax.tick_params(labelsize=7)

    out_path = os.path.join(output_dir, "fig7_separation_displacement.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
