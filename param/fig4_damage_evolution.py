#!/usr/bin/env python3
"""Figure 4: expected damage evolution and fractured elements."""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    cycles = np.arange(0, 401)
    d_max = 1.0 - np.exp(-cycles / 170.0)
    d_mean = 0.65 * (1.0 - np.exp(-cycles / 230.0))
    fractured = np.round(320.0 * np.clip((cycles - 60.0) / 300.0, 0.0, 1.0) ** 1.25)

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300)

    axes[0].plot(cycles, d_max * 100.0, lw=2.0, label="D_max")
    axes[0].plot(cycles, d_mean * 100.0, lw=2.0, ls="--", label="D_mean")
    axes[0].set_title("(a) Damage index")
    axes[0].set_xlabel("Cycle number")
    axes[0].set_ylabel("Damage (%)")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)

    axes[1].plot(cycles, fractured, color="#d62728", lw=2.0)
    axes[1].set_title("(b) Fractured cohesive elements")
    axes[1].set_xlabel("Cycle number")
    axes[1].set_ylabel("Count")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig4_damage_evolution.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
