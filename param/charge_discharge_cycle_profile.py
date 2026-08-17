#!/usr/bin/env python3
"""Plot the four-step constant-current charge-discharge cycle profile.

The JuBat sign convention is used: discharge current is positive and charge
current is negative.  Run this file directly to create the schematic in the
repository's ``output`` directory.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt


DISCHARGE_END_S = 2400
REST_1_END_S = 3000
CHARGE_END_S = 5400
CYCLE_END_S = 6000


def make_figure() -> Path:
    """Create and save the 1C charge-discharge cycle schematic."""
    root = Path(__file__).resolve().parent.parent
    output_dir = root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 9,
            "axes.unicode_minus": False,
        }
    )

    # Explicit vertices retain the vertical current changes at phase boundaries.
    time_s = [
        0,
        0,
        DISCHARGE_END_S,
        DISCHARGE_END_S,
        REST_1_END_S,
        REST_1_END_S,
        CHARGE_END_S,
        CHARGE_END_S,
        CYCLE_END_S,
    ]
    c_rate = [0, 1, 1, 0, 0, -1, -1, 0, 0]

    phases = [
        (0, DISCHARGE_END_S, "Step 1\nDischarge\n1C, 2400 s", "#dbeafe"),
        (DISCHARGE_END_S, REST_1_END_S, "Step 2\nRest\n600 s", "#f3f4f6"),
        (REST_1_END_S, CHARGE_END_S, "Step 3\nCharge\n1C, 2400 s", "#fee2e2"),
        (CHARGE_END_S, CYCLE_END_S, "Step 4\nRest\n600 s", "#f3f4f6"),
    ]

    fig, ax = plt.subplots(figsize=(5.6, 4.0), dpi=300)

    for start_s, end_s, label, color in phases:
        ax.axvspan(start_s, end_s, color=color, alpha=0.75, zorder=0)
        label_y = 1.24 if end_s - start_s > 600 else 0.48
        ax.text(
            (start_s + end_s) / 2,
            label_y,
            label,
            ha="center",
            va="center",
            fontsize=8,
        )

    ax.plot(time_s, c_rate, color="#1f4e79", linewidth=2.2, zorder=2)
    ax.axhline(0, color="#4b5563", linewidth=0.8, zorder=1)

    for boundary_s in (DISCHARGE_END_S, REST_1_END_S, CHARGE_END_S):
        ax.axvline(boundary_s, color="#6b7280", linestyle="--", linewidth=0.7, alpha=0.7)

    ax.set_title("Four-step 1C charge-discharge cycle profile")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Current rate (C)")
    ax.set_xlim(0, CYCLE_END_S)
    ax.set_ylim(-1.45, 1.55)
    ax.set_xticks([0, 2400, 3000, 5400, 6000])
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(["-1C  Charge", "0  Rest", "+1C  Discharge"])
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.45)

    fig.tight_layout()
    output_path = output_dir / "charge_discharge_cycle_profile.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_path


if __name__ == "__main__":
    print(f"Saved: {make_figure()}")
