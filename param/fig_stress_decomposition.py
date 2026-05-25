#!/usr/bin/env python3
"""Figure: thermal stress vs diffusion stress decomposition.

Two panels:
  (a) Stress components vs cycle for different tab configurations (1C)
  (b) Stress components vs cycle for different charge rates (single tab)

Physical model
--------------
Thermal stress σ_th ∝ α·E·ΔT:
  - Fewer tabs → higher local current density → more Joule heating → larger ΔT
  - Higher C-rate → larger reaction heat + Joule heat → larger ΔT
  - Non-zero at N=0, grows modestly with damage-induced thermal resistance feedback

Diffusion stress σ_diff ∝ E·Ω·Δc:
  - Higher C-rate → larger concentration gradient → larger σ_diff
  - Temperature partially relaxes it: higher T → faster diffusion → lower Δc
  - Non-zero at N=0, grows then slightly saturates

σ_total = σ_th + σ_diff drives CZM damage when σ_total > σ_crit
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _thermal_stress(
    cycles: np.ndarray,
    sigma_init: float,
    sigma_final: float,
    tau: float = 120.0,
) -> np.ndarray:
    """Thermal stress: non-zero at N=0, rises with damage feedback."""
    return sigma_init + (sigma_final - sigma_init) * (1.0 - np.exp(-cycles / tau))


def _diffusion_stress(
    cycles: np.ndarray,
    sigma_init: float,
    sigma_peak: float,
    sigma_final: float,
    tau_rise: float = 80.0,
    tau_relax: float = 300.0,
) -> np.ndarray:
    """Diffusion stress: rises from init to peak, then relaxes toward final.

    At N=0:   sigma_init
    At peak:  ~sigma_peak
    At N→∞:   sigma_final
    """
    rise = (sigma_peak - sigma_init) * (1.0 - np.exp(-cycles / tau_rise))
    relax = (sigma_peak - sigma_final) * (1.0 - np.exp(-cycles / tau_relax))
    return sigma_init + rise - relax


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    N = np.linspace(0, 500, 1000)
    sigma_crit = 80.0

    # ── Panel (a): Tab configurations (1C) ──
    # σ_th dominates the difference; σ_diff similar across tabs at same C-rate
    tab_configs = {
        "Triple tabs": {
            "th_init": 22, "th_final": 28,
            "diff_init": 66, "diff_peak": 74, "diff_final": 72,
            "color": "#1f77b4",
        },
        "Double tabs": {
            "th_init": 28, "th_final": 36,
            "diff_init": 82, "diff_peak": 92, "diff_final": 88,
            "color": "#ff7f0e",
        },
        "Single tab": {
            "th_init": 34, "th_final": 44,
            "diff_init": 102, "diff_peak": 114, "diff_final": 110,
            "color": "#d62728",
        },
    }

    # ── Panel (b): Charge rates (single tab) ──
    rate_configs = {
        "0.5C": {
            "th_init": 18, "th_final": 22,
            "diff_init": 54, "diff_peak": 60, "diff_final": 58,
            "color": "#1f77b4",
        },
        "1C": {
            "th_init": 34, "th_final": 44,
            "diff_init": 102, "diff_peak": 114, "diff_final": 110,
            "color": "#2ca02c",
        },
        "2C": {
            "th_init": 42, "th_final": 52,
            "diff_init": 138, "diff_peak": 154, "diff_final": 148,
            "color": "#ff7f0e",
        },
        "3C": {
            "th_init": 48, "th_final": 60,
            "diff_init": 168, "diff_peak": 188, "diff_final": 180,
            "color": "#d62728",
        },
    }

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)

    # ────────────── Panel (a) ──────────────
    for label, cfg in tab_configs.items():
        sth = _thermal_stress(N, cfg["th_init"], cfg["th_final"])
        sdiff = _diffusion_stress(N, cfg["diff_init"], cfg["diff_peak"], cfg["diff_final"])

        ax_a.plot(N, sth,   color=cfg["color"], lw=1.2, ls="--", alpha=0.85)
        ax_a.plot(N, sdiff, color=cfg["color"], lw=1.2, ls=":",  alpha=0.85)

    ax_a.axhline(sigma_crit, color="#555555", ls="-.", lw=0.8, alpha=0.7)
    ax_a.text(5, sigma_crit + 2, r"$\sigma_{\mathrm{crit}}$", fontsize=7,
              color="#555555", va="bottom")

    ax_a.set_xlabel("Cycle number N", fontsize=9)
    ax_a.set_ylabel(r"Stress $\sigma$ (MPa)", fontsize=9)
    ax_a.set_title("(a) Tab configurations (1C)", fontsize=9, pad=6)
    ax_a.set_xlim(0, 500)
    ax_a.set_ylim(10, 200)
    ax_a.set_xticks(np.arange(0, 501, 100))
    ax_a.tick_params(labelsize=7.5)
    ax_a.grid(True, ls="--", lw=0.3, alpha=0.4)

    # ────────────── Panel (b) ──────────────
    for label, cfg in rate_configs.items():
        sth = _thermal_stress(N, cfg["th_init"], cfg["th_final"])
        sdiff = _diffusion_stress(N, cfg["diff_init"], cfg["diff_peak"], cfg["diff_final"])

        ax_b.plot(N, sth,   color=cfg["color"], lw=1.2, ls="--", alpha=0.85)
        ax_b.plot(N, sdiff, color=cfg["color"], lw=1.2, ls=":",  alpha=0.85)

    ax_b.axhline(sigma_crit, color="#555555", ls="-.", lw=0.8, alpha=0.7)
    ax_b.text(5, sigma_crit + 2, r"$\sigma_{\mathrm{crit}}$", fontsize=7,
              color="#555555", va="bottom")

    ax_b.set_xlabel("Cycle number N", fontsize=9)
    ax_b.set_ylabel(r"Stress $\sigma$ (MPa)", fontsize=9)
    ax_b.set_title("(b) Charge rates (single tab)", fontsize=9, pad=6)
    ax_b.set_xlim(0, 500)
    ax_b.set_ylim(10, 300)
    ax_b.set_xticks(np.arange(0, 501, 100))
    ax_b.tick_params(labelsize=7.5)
    ax_b.grid(True, ls="--", lw=0.3, alpha=0.4)

    # ── Legends ──
    # Line type legend (shared, bottom center)
    type_handles = [
        Line2D([0], [0], color="#333333", lw=1.2, ls="--", alpha=0.7,
               label=r"$\sigma_{\mathrm{th}}$ (thermal)"),
        Line2D([0], [0], color="#333333", lw=1.2, ls=":", alpha=0.7,
               label=r"$\sigma_{\mathrm{diff}}$ (diffusion)"),
        Line2D([0], [0], color="#555555", ls="-.", lw=0.8,
               label=r"$\sigma_{\mathrm{crit}}$"),
    ]
    fig.legend(
        handles=type_handles, loc="lower center",
        fontsize=6.5, frameon=True, framealpha=0.9,
        prop={"family": "Times New Roman"},
        ncol=3, bbox_to_anchor=(0.5, -0.04),
    )

    # Panel (a) config legend
    tab_handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, label=label)
        for label, cfg in tab_configs.items()
    ]
    ax_a.legend(
        handles=tab_handles, fontsize=6.5, frameon=True, framealpha=0.9,
        loc="upper left", prop={"family": "Times New Roman"},
    )

    # Panel (b) config legend
    rate_handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, label=label)
        for label, cfg in rate_configs.items()
    ]
    ax_b.legend(
        handles=rate_handles, fontsize=6.5, frameon=True, framealpha=0.9,
        loc="upper left", prop={"family": "Times New Roman"},
    )

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)

    out_path = os.path.join(output_dir, "fig_stress_decomposition.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
