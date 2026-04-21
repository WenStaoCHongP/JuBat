#!/usr/bin/env python3
"""Figure 2: multiphysics coupling framework."""

from __future__ import annotations

import os
from plot import figure_2_coupling_schematic


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    fig = figure_2_coupling_schematic(save=False)
    out_path = os.path.join(output_dir, "fig2_coupling_framework.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    return out_path


if __name__ == "__main__":
    print(make_figure())
