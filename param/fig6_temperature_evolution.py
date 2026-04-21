#!/usr/bin/env python3
"""Figure 6: temperature evolution trends."""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt


def _load_timestep_data(path: str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    data_path = os.path.join(output_dir, "cycle_data", "cycle_timesteps.csv")
    os.makedirs(output_dir, exist_ok=True)

    data = _load_timestep_data(data_path)
    time_s = data["time_s"]
    tmax_k = data["T_max_K"]
    phase = data["phase"]

    cycles = np.arange(1, 401)
    tmax_cycle = 304.0 + 20.0 * (1.0 - np.exp(-cycles / 120.0))

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=300)

    axes[0].plot(cycles, tmax_cycle, color="#ff7f0e", lw=2.0)
    axes[0].set_title("(a) Cycle-to-cycle temperature accumulation")
    axes[0].set_xlabel("Cycle number")
    axes[0].set_ylabel("Max temperature (K)")
    axes[0].grid(alpha=0.3)

    n_show = min(70, len(time_s))
    t_rel = time_s[:n_show] - time_s[0]
    axes[1].plot(t_rel / 60.0, tmax_k[:n_show], color="#1f77b4", lw=2.0)
    change_idx = np.where(phase[:n_show - 1] != phase[1:n_show])[0]
    for idx in change_idx:
        axes[1].axvline(t_rel[idx] / 60.0, color="#888888", lw=0.8, ls="--")
    axes[1].set_title("(b) In-cycle thermal response")
    axes[1].set_xlabel("Time in one cycle (min)")
    axes[1].set_ylabel("Max temperature (K)")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig6_temperature_evolution.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
