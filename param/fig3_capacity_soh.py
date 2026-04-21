#!/usr/bin/env python3
"""Figure 3: expected capacity fade and SOH evolution."""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    cycles = np.arange(0, 401)
    cap0 = 5.0
    capacity = cap0 * (1.0 - 0.00028 * cycles - 0.00000055 * cycles**2)
    capacity = np.clip(capacity, 3.6, None)
    soh = 100.0 * capacity / cap0

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300)

    axes[0].plot(cycles, capacity, color="#1f77b4", lw=2.0)
    axes[0].axhline(4.0, color="#808080", ls="--", lw=1.0, label="4 Ah reference")
    axes[0].set_xlabel("Cycle number")
    axes[0].set_ylabel("Discharge capacity (Ah)")
    axes[0].set_title("(a) Capacity fade")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)

    axes[1].plot(cycles, soh, color="#d62728", lw=2.0)
    axes[1].axhline(80.0, color="#808080", ls="--", lw=1.0, label="SOH = 80%")
    axes[1].set_xlabel("Cycle number")
    axes[1].set_ylabel("SOH (%)")
    axes[1].set_title("(b) SOH evolution")
    axes[1].set_ylim(70, 102)
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig3_capacity_soh.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
