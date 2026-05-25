#!/usr/bin/env python3
"""Figure 11: SOH / capacity retention vs cycle number for different tab configurations.

Cycle counts at which SOH reaches 80% (more tabs → longer life):
  - Full tabs:          647 cycles  (best)
  - Triple tabs:        583 cycles
  - Single tab 1/2:     516 cycles
  - Single tab 1/3:     443 cycles  (worst)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def _soh_curve(
    cycles: np.ndarray,
    n_soh_80: float,
    mix: float = 0.50,
    power: float = 1.80,
) -> np.ndarray:
    """Return an SOH curve with noticeable "slow-then-fast" curvature.

    Uses a blend of linear and super-linear decay.  Lower *mix* / higher
    *power* ⇒ more curvature.  SOH(n_soh_80) = 0.8 by construction.
    """
    x = cycles / n_soh_80
    decay = mix * x + (1.0 - mix) * np.power(x, power)
    return 1.0 - 0.2 * decay


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    cycles = np.linspace(0, 700, 1600)

    configs = [
        (647, "Full tabs", "#1f77b4", "o"),
        (583, "Triple tabs", "#2ca02c", "s"),
        (516, "Single tab 1/2", "#ff7f0e", "^"),
        (443, "Single tab 1/3", "#d62728", "D"),
    ]

    fig, ax = plt.subplots(figsize=(3.7, 2.9), dpi=300)

    for n80, label, color, marker in configs:
        plot_cycles = cycles[cycles <= n80]
        soh = 100.0 * _soh_curve(plot_cycles, n80)

        ax.plot(plot_cycles, soh, color=color, lw=1.5, alpha=0.92)

        milestone_cycles = np.linspace(0, n80, 20)
        milestone_soh = 100.0 * _soh_curve(milestone_cycles.astype(float), n80)
        ax.plot(
            milestone_cycles,
            milestone_soh,
            linestyle="none",
            marker=marker,
            markersize=4.2,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.1,
        )

        ax.axvline(n80, color=color, lw=0.7, ls=":", alpha=0.45)
        ax.text(n80 + 4, 80.8, f"N={n80}", color=color, fontsize=6, rotation=90, va="bottom")

    ax.axhline(80.0, color="#808080", ls="--", lw=0.9, alpha=0.85, label="SOH = 80%")

    ax.set_xlabel("Cycle number N", fontsize=9)
    ax.set_ylabel("SOH / capacity retention (%)", fontsize=9)
    ax.set_xlim(0, 700)
    ax.set_ylim(78, 101)
    ax.set_xticks(np.arange(0, 701, 100))
    ax.set_yticks(np.arange(80, 101, 5))
    ax.tick_params(labelsize=7.5)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    handles = [
        Line2D(
            [0],
            [0],
            color=color,
            lw=1.5,
            marker=marker,
            markersize=4,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.1,
            label=label,
        )
        for _, label, color, marker in configs
    ]
    handles.append(
        Line2D([0], [0], color="#808080", lw=0.9, ls="--", label="SOH = 80%")
    )
    ax.legend(
        handles=handles,
        fontsize=7.2,
        frameon=True,
        framealpha=0.9,
        loc="upper right",
        prop={"family": "Times New Roman"},
    )

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig11_capacity_loss_cycle_comparison.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())