#!/usr/bin/env python3
"""Figure 10: average damage D vs cycle number for different tab configurations.

Three curves comparing:
  - Single tab:   earliest onset, fastest growth
  - Double tabs:  moderate onset, moderate growth
  - Triple tabs:  latest onset, slowest growth

Fewer tabs → higher local current density → earlier and faster damage.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _damage_curve(
    N: np.ndarray,
    N_onset: float,
    D_sat: float,
    rate: float,
    smooth: float = 12.0,
    accel: float = 2.2,
) -> np.ndarray:
    """Saturating curve with non-linear slow-start, fast-growth profile."""
    delta = N - N_onset
    sigmoid = 1.0 / (1.0 + np.exp(-delta / smooth))
    base = np.where(delta > 0, 1.0 - np.exp(-rate * np.clip(delta, 0, None)), 0.0)
    growth = np.power(base, accel)
    return D_sat * sigmoid * growth


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Style ──
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    N = np.linspace(0, 500, 2000)

    # ── Curve parameters ──
    # Fewer tabs → earlier onset, higher saturation, faster rate
    #              onset  D_sat   rate    description
    configs = [
        (55,  0.185, 0.0080, "Single tab"),    # earliest onset, fastest growth
        (76,  0.122, 0.0052, "Double tabs"),    # moderate
        (97,  0.085, 0.0038, "Triple tabs"),    # latest onset, slowest growth
    ]
    colors = ["#d62728", "#ff7f0e", "#1f77b4"]
    markers = ["D", "s", "o"]

    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=300)

    for (N_onset, D_sat, rate, label), color, marker in zip(configs, colors, markers):
        D = _damage_curve(N, N_onset, D_sat, rate)

        # Continuous curve
        ax.plot(N, D, color=color, lw=1.4, alpha=0.9)

        # Sampled markers at integer cycle milestones
        milestone_N = np.arange(0, 501, 50)
        milestone_D = _damage_curve(milestone_N.astype(float), N_onset, D_sat, rate)
        ax.plot(
            milestone_N, milestone_D,
            marker=marker, markersize=4,
            markerfacecolor="none", markeredgecolor=color,
            markeredgewidth=1.2, linestyle="none",
        )

        # Vertical dashed line at onset
        ax.axvline(N_onset, color=color, lw=0.7, ls=":", alpha=0.5)
        ax.text(
            N_onset + 3, 0.03, f"N={N_onset}",
            fontsize=6, color=color, fontname="Times New Roman",
            rotation=90, va="bottom",
        )

    # ── Axes ──
    ax.set_xlabel("Cycle number N", fontsize=9)
    ax.set_ylabel(r"Average damage $\bar{D}$", fontsize=9)
    ax.set_xlim(0, 500)
    ax.set_ylim(0, 0.20)
    ax.set_xticks(np.arange(0, 501, 50))
    ax.tick_params(labelsize=7.5)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Legend
    handles = [
        Line2D([0], [0], color=c, lw=1.4, marker=m, markersize=4,
               markerfacecolor="none", markeredgecolor=c, markeredgewidth=1.2, label=l)
        for c, m, l in zip(colors, markers, [c[3] for c in configs])
    ]
    ax.legend(
        handles=handles, fontsize=7.5, frameon=True, framealpha=0.9,
        loc="upper left", prop={"family": "Times New Roman"},
    )

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig10_damage_cycle_comparison.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
