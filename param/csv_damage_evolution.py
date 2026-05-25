#!/usr/lib/python3
"""Damage and separation displacement evolution from CSV data.

Plots:
  1. D_max(t) over time (per phase)
  2. D(theta) profile at specified time
  3. delta_n(theta) profile at specified time
  4. D(t) for the most damaged element
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def make_figure(
    csv_dir: str | None = None,
    target_time: float = 1200.0,
    cycle: int = 1,
) -> str:
    """Plot CZM damage and separation evolution. Returns path to saved PNG."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if csv_dir is None:
        csv_dir = os.path.join(root, "output", "csv", "czm_study_1")
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Load cohesive damage ──
    df = pd.read_csv(os.path.join(csv_dir, "cohesive_damage.csv"))
    df = df[df["cycle"] == cycle].copy()

    if len(df) == 0:
        print(f"No CZM data for cycle {cycle}")
        return ""

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=300, constrained_layout=True)

    # ── (a) D_max over time per phase ──
    ax = axes[0, 0]
    for phase, color, ls in [
        ("PHASE_DISCHARGE", "#1f77b4", "-"),
        ("PHASE_REST", "#2ca02c", "--"),
        ("PHASE_CHARGE", "#d62728", "-."),
    ]:
        sub = df[df["phase"] == phase]
        if len(sub) == 0:
            continue
        max_d = sub.groupby("time_s")["D"].max()
        ax.plot(max_d.index, max_d.values, color=color, ls=ls, lw=1.2, label=phase.replace("PHASE_", ""))
    ax.set_xlabel("Time (s)", fontsize=9)
    ax.set_ylabel(r"$D_{\max}$", fontsize=10)
    ax.set_title("(a) Maximum damage evolution", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # ── (b) D(theta) profile at target_time ──
    ax = axes[0, 1]
    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - target_time))]
    sub_t = df[df["time_s"] == t_closest]
    if len(sub_t) > 0:
        # Sort by theta
        sub_sorted = sub_t.sort_values("theta_deg")
        ax.plot(sub_sorted["theta_deg"], sub_sorted["D"], "b-", lw=0.6, alpha=0.8)
        ax.fill_between(sub_sorted["theta_deg"], sub_sorted["D"], alpha=0.15, color="blue")

        # Also plot at final time
        t_final = times[-1]
        sub_f = df[df["time_s"] == t_final].sort_values("theta_deg")
        ax.plot(sub_f["theta_deg"], sub_f["D"], "r-", lw=0.8, alpha=0.7, label=f"t={t_final:.0f}s")

    ax.set_xlabel(r"$\theta$ (deg)", fontsize=9)
    ax.set_ylabel("Damage D", fontsize=9)
    ax.set_title(f"(b) Damage profile at t = {t_closest:.0f} s", fontsize=10)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # ── (c) delta_n(theta) at target_time ──
    ax = axes[1, 0]
    if len(sub_t) > 0:
        sub_sorted = sub_t.sort_values("theta_deg")
        dn_um = sub_sorted["delta_n"].values * 1e6  # um
        ax.plot(sub_sorted["theta_deg"], dn_um, "b-", lw=0.6, alpha=0.8)
        ax.fill_between(sub_sorted["theta_deg"], dn_um, alpha=0.15, color="blue")

        # Final time
        if len(sub_f) > 0:
            dn_f = sub_f.sort_values("theta_deg")["delta_n"].values * 1e6
            th_f = sub_f.sort_values("theta_deg")["theta_deg"].values
            ax.plot(th_f, dn_f, "r-", lw=0.8, alpha=0.7, label=f"t={t_final:.0f}s")

    ax.set_xlabel(r"$\theta$ (deg)", fontsize=9)
    ax.set_ylabel(r"Normal separation $\delta_n$ ($\mu$m)", fontsize=9)
    ax.set_title(f"(c) Separation profile at t = {t_closest:.0f} s", fontsize=10)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # ── (d) D(t) for the most damaged element ──
    ax = axes[1, 1]
    # Find element with max D at final time
    final_time = df["time_s"].max()
    sub_final = df[df["time_s"] == final_time]
    if len(sub_final) > 0:
        max_elem = sub_final.loc[sub_final["D"].idxmax(), "coh_id"]

        elem_data = df[df["coh_id"] == max_elem].sort_values("time_s")
        ax.plot(elem_data["time_s"], elem_data["D"], "r-", lw=1.2)
        ax.fill_between(elem_data["time_s"], elem_data["D"], alpha=0.15, color="red")

        # Also plot delta_n on secondary axis
        ax2 = ax.twinx()
        dn_elem = elem_data["delta_n"].values * 1e6  # um
        ax2.plot(elem_data["time_s"], dn_elem, "b--", lw=0.8, alpha=0.7)
        ax2.set_ylabel(r"$\delta_n$ ($\mu$m)", fontsize=9, color="blue")
        ax2.tick_params(axis="y", labelsize=8, labelcolor="blue")

        ax.set_xlabel("Time (s)", fontsize=9)
        ax.set_ylabel("Damage D", fontsize=9, color="red")
        ax.tick_params(axis="y", labelsize=8, labelcolor="red")
        ax.set_title(f"(d) Most damaged element (id={max_elem})", fontsize=10)
    else:
        ax.set_title("(d) No damage data", fontsize=10)

    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"CZM Damage & Separation Evolution  |  Cycle {cycle}", fontsize=12, fontweight="bold")

    out_path = os.path.join(output_dir, "csv_damage_evolution.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


if __name__ == "__main__":
    make_figure()
