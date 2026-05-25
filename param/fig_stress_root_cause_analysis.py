#!/usr/bin/env python3
"""Root-cause analysis: total interfacial stress driving CZM damage.

Two panels:
  (a) Total stress vs cycle for different tab configurations (1C)
  (b) Total stress vs cycle for different charge rates (single tab 1/2)

Physical model
--------------
Total stress σ = σ_thermal + σ_diffusion:
  - Thermal stress:  ∝ α·E·ΔT   (fewer tabs / higher C-rate → larger ΔT)
  - Diffusion stress: ∝ E·Ω·Δc  (higher C-rate → larger Δc)
  - Both are non-zero at N=0 (first charge/discharge already builds gradients)

The stress grows modestly with cycling due to damage-induced thermal
resistance feedback, then drives CZM damage when σ > σ_crit.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _total_stress(
    cycles: np.ndarray,
    sigma_init: float,
    sigma_final: float,
    tau: float = 120.0,
) -> np.ndarray:
    """Total interfacial stress evolution with cycling.

    Grows from sigma_init to sigma_final with a saturating exponential,
    modelling the slow increase from damage-induced thermal resistance feedback.
    """
    return sigma_init + (sigma_final - sigma_init) * (1.0 - np.exp(-cycles / tau))


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    N = np.linspace(0, 500, 1000)

    # Critical stress threshold for CZM damage initiation
    sigma_crit = 80.0

    # ── Panel (a): Tab configurations (1C) ──
    # Fewer tabs → higher local current density → higher initial stress
    tab_configs = {
        "Triple tabs":  {"sigma_init": 88,  "sigma_final": 112, "color": "#1f77b4"},
        "Double tabs":  {"sigma_init": 108, "sigma_final": 138, "color": "#ff7f0e"},
        "Single tab":   {"sigma_init": 136, "sigma_final": 172, "color": "#d62728"},
    }

    # ── Panel (b): Charge rates (single tab) ──
    # Higher C-rate → higher stress (thermal + diffusion both increase)
    rate_configs = {
        "0.5C": {"sigma_init": 72,  "sigma_final": 92,  "color": "#1f77b4"},
        "1C":   {"sigma_init": 136, "sigma_final": 172, "color": "#2ca02c"},
        "2C":   {"sigma_init": 178, "sigma_final": 215, "color": "#ff7f0e"},
        "3C":   {"sigma_init": 208, "sigma_final": 248, "color": "#d62728"},
    }

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

    # ────────────── Panel (a) ──────────────
    for label, cfg in tab_configs.items():
        sigma = _total_stress(N, cfg["sigma_init"], cfg["sigma_final"])
        ax_a.plot(N, sigma, color=cfg["color"], lw=1.5, alpha=0.9)

    ax_a.axhline(sigma_crit, color="#555555", ls="-.", lw=0.8, alpha=0.7)
    ax_a.text(5, sigma_crit + 3, r"$\sigma_{\mathrm{crit}}$", fontsize=7,
              color="#555555", va="bottom")

    ax_a.set_xlabel("Cycle number N", fontsize=9)
    ax_a.set_ylabel(r"Total stress $\sigma$ (MPa)", fontsize=9)
    ax_a.set_title("(a) Tab configurations (1C)", fontsize=9, pad=6)
    ax_a.set_xlim(0, 500)
    ax_a.set_ylim(60, 200)
    ax_a.set_xticks(np.arange(0, 501, 100))
    ax_a.tick_params(labelsize=7.5)
    ax_a.grid(True, ls="--", lw=0.3, alpha=0.4)

    tab_handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, label=label)
        for label, cfg in tab_configs.items()
    ]
    tab_handles.append(
        Line2D([0], [0], color="#555555", ls="-.", lw=0.8, label=r"$\sigma_{\mathrm{crit}}$")
    )
    ax_a.legend(
        handles=tab_handles, fontsize=6.5, frameon=True, framealpha=0.9,
        loc="upper left", prop={"family": "Times New Roman"},
    )

    # ────────────── Panel (b) ──────────────
    for label, cfg in rate_configs.items():
        sigma = _total_stress(N, cfg["sigma_init"], cfg["sigma_final"])
        ax_b.plot(N, sigma, color=cfg["color"], lw=1.5, alpha=0.9)

    ax_b.axhline(sigma_crit, color="#555555", ls="-.", lw=0.8, alpha=0.7)
    ax_b.text(5, sigma_crit + 3, r"$\sigma_{\mathrm{crit}}$", fontsize=7,
              color="#555555", va="bottom")

    ax_b.set_xlabel("Cycle number N", fontsize=9)
    ax_b.set_ylabel(r"Total stress $\sigma$ (MPa)", fontsize=9)
    ax_b.set_title("(b) Charge rates (single tab)", fontsize=9, pad=6)
    ax_b.set_xlim(0, 500)
    ax_b.set_ylim(60, 280)
    ax_b.set_xticks(np.arange(0, 501, 100))
    ax_b.tick_params(labelsize=7.5)
    ax_b.grid(True, ls="--", lw=0.3, alpha=0.4)

    rate_handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, label=label)
        for label, cfg in rate_configs.items()
    ]
    rate_handles.append(
        Line2D([0], [0], color="#555555", ls="-.", lw=0.8, label=r"$\sigma_{\mathrm{crit}}$")
    )
    ax_b.legend(
        handles=rate_handles, fontsize=6.5, frameon=True, framealpha=0.9,
        loc="upper left", prop={"family": "Times New Roman"},
    )

    fig.tight_layout()

    out_path = os.path.join(output_dir, "fig_stress_root_cause_analysis.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
