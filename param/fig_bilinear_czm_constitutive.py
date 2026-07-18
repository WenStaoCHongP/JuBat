#!/usr/bin/env python3
"""Bilinear traction-separation constitutive law for CZM.

Plots the Mode I bilinear cohesive zone model with annotated key
parameters (stiffness, peak traction, fracture energy, critical
separation).
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np


# =====================================================================
# CZM parameters (from src/parameters/Jellyroll.jl)
# =====================================================================
SIGMA_MAX = 82e6        # Pa  – maximum normal traction
K_N       = 2.4e17      # Pa/m – initial stiffness
G_C       = 25.3        # J/m² – Mode I fracture energy

DELTA_0   = SIGMA_MAX / K_N               # damage onset separation [m]
DELTA_C   = 2.0 * G_C / SIGMA_MAX         # complete fracture separation [m]


def _bilinear(delta: np.ndarray) -> np.ndarray:
    """Evaluate bilinear traction-separation law."""
    sigma = np.zeros_like(delta)
    mask1 = delta <= DELTA_0
    sigma[mask1] = K_N * delta[mask1]
    mask2 = (delta > DELTA_0) & (delta <= DELTA_C)
    sigma[mask2] = SIGMA_MAX * (DELTA_C - delta[mask2]) / (DELTA_C - DELTA_0)
    return sigma


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Style
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    # Data
    delta = np.linspace(0, DELTA_C * 1.15, 600)
    sigma = _bilinear(delta)

    delta_um   = delta * 1e6       # m -> μm
    sigma_MPa  = sigma * 1e-6      # Pa -> MPa
    s_max      = SIGMA_MAX * 1e-6  # MPa
    delta0_um  = DELTA_0 * 1e6
    deltaC_um  = DELTA_C * 1e6

    fig, ax = plt.subplots(figsize=(5.0, 3.8), dpi=300)

    # ---- Main curve ----
    mask_el    = delta <= DELTA_0
    mask_soft  = (delta > DELTA_0) & (delta <= DELTA_C)
    mask_frac  = delta > DELTA_C

    ax.plot(delta_um[mask_el],   sigma_MPa[mask_el],   "b-",  lw=1.8, label="Elastic")
    ax.plot(delta_um[mask_soft], sigma_MPa[mask_soft],  "r-",  lw=1.8, label="Softening")
    ax.plot(delta_um[mask_frac], sigma_MPa[mask_frac],  "k--", lw=1.0, alpha=0.5, label="Fully fractured")

    # ---- Key points ----
    ax.plot(delta0_um, s_max, "ko", ms=6, zorder=5)
    ax.plot(deltaC_um, 0,     "kx", ms=8, mew=2, zorder=5)
    ax.plot(0, 0, "ko", ms=5, zorder=5)

    # ---- G_c shaded area (draw first so annotations sit on top) ----
    tri_x = [0, delta0_um, deltaC_um, 0]
    tri_y = [0, s_max, 0, 0]
    ax.fill(tri_x, tri_y, color="green", alpha=0.08)

    # ---- Dashed guide lines ----
    ax.plot([delta0_um, delta0_um], [0, s_max], "k:", lw=0.6, alpha=0.4)
    ax.plot([0, delta0_um], [s_max, s_max], "k:", lw=0.6, alpha=0.4)

    # ========== Annotations (spread into 4 corners) ==========

    # 1. σ_max — upper-left, arrow points at peak
    ax.annotate(
        r"$\sigma_{\max}$" + f" = {s_max:.0f} MPa",
        xy=(delta0_um, s_max),
        xytext=(deltaC_um * 0.05, s_max * 1.15),
        fontsize=8.5, ha="left",
        arrowprops=dict(arrowstyle="->", lw=0.8, color="black"),
    )

    # 2. δ_0 — below x-axis, left
    ax.annotate(
        r"$\delta_0$",
        xy=(delta0_um, 0),
        xytext=(deltaC_um * 0.08, -s_max * 0.12),
        fontsize=9, ha="center",
        arrowprops=dict(arrowstyle="->", lw=0.6, color="gray"),
    )

    # 3. δ_c — below x-axis, right
    ax.annotate(
        r"$\delta_c$",
        xy=(deltaC_um, 0),
        xytext=(deltaC_um * 0.92, -s_max * 0.12),
        fontsize=9, ha="center",
        arrowprops=dict(arrowstyle="->", lw=0.6, color="gray"),
    )

    # 4. G_c — centre of triangle (no arrow, plain text)
    ax.text(
        deltaC_um * 0.45, s_max * 0.35,
        r"$G_c = \frac{1}{2}\sigma_{\max}\delta_c$" + f" = {G_C:.1f} J/m²",
        fontsize=8.5, color="green", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"),
    )

    # 5. Parameter box — lower-right (clearest zone)
    param_text = (
        f"$K_n = {K_N:.1e}$ Pa/m\n"
        f"$\\sigma_{{\\max}} = {s_max:.0f}$ MPa\n"
        f"$\\delta_0 = {delta0_um:.4f}\\;\\mu$m\n"
        f"$G_c = {G_C:.1f}$ J/m²\n"
        f"$\\delta_c = {deltaC_um:.3f}\\;\\mu$m"
    )
    ax.text(
        0.97, 0.40, param_text,
        transform=ax.transAxes, fontsize=7.5,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.85, edgecolor="gray"),
    )

    # ---- Axes ----
    ax.set_xlabel(r"Normal separation $\delta_n$ (μm)", fontsize=10)
    ax.set_ylabel(r"Normal traction $\sigma_n$ (MPa)", fontsize=10)
    ax.set_xlim(-deltaC_um * 0.08, deltaC_um * 1.20)
    ax.set_ylim(-s_max * 0.18, s_max * 1.30)
    ax.tick_params(labelsize=8)
    ax.grid(True, ls="--", lw=0.3, alpha=0.4)

    # Legend — upper-right
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "fig_bilinear_czm_constitutive.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    path = make_figure()
    print(f"Saved: {path}")
