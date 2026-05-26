#!/usr/bin/env python3
"""Unified CZM post-processing plotter.

Produces cloud maps, evolution plots, and geometry visualizations from
JuBat cycling simulation CSV output.

Usage: Edit CONFIG dict at top, then run:
    python param/plot_czm.py
"""

from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from matplotlib.colors import Normalize

# =====================================================================
# User Configuration
# =====================================================================

CONFIG = {
    "data_dir": "output/csv/czm_study_1",   # simulation result CSVs
    "mesh_dir": "output/csv/20",             # mesh geometry CSVs
    "output_dir": "output",                  # PNG output directory

    "cycle": 1,                              # which cycle to plot
    "target_time": 1200.0,                   # target time (s)
    "phase": "discharge",                    # phase name (for thermal/current)

    "plots": [
        "temperature",
        "displacement",
        "current",
        "damage_evolution",
        "capacity_soh",
        "mesh_geometry",
        "traction_separation",
        "separation_cloud",
    ],
}

# =====================================================================
# Shared Utilities
# =====================================================================

def _resolve_path(cfg: dict, name: str) -> str:
    """Find a CSV file in data_dir or mesh_dir."""
    for d in (cfg["data_dir"], cfg["mesh_dir"]):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(f"{name} not found in {cfg['data_dir']} or {cfg['mesh_dir']}")


def _load_mesh_geometry(cfg: dict):
    """Load thermal mesh geometry, deduplicate, return Triangulation.

    Returns (node_x_orig, node_y_orig, tri, orig_to_unique)
    """
    df_nodes = pd.read_csv(_resolve_path(cfg, "mesh_nodes.csv"))
    df_elems = pd.read_csv(_resolve_path(cfg, "mesh_elements.csv"))

    node_x_orig = df_nodes["x"].values * 1000  # m -> mm
    node_y_orig = df_nodes["y"].values * 1000

    xy_rounded = np.round(np.column_stack([node_x_orig, node_y_orig]), decimals=12)
    _, unique_idx, orig_to_unique = np.unique(
        xy_rounded, axis=0, return_index=True, return_inverse=True
    )

    triangles = []
    for _, row in df_elems.iterrows():
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        u1, u2, u3, u4 = orig_to_unique[n1], orig_to_unique[n2], orig_to_unique[n3], orig_to_unique[n4]
        triangles.append([u1, u2, u3])
        triangles.append([u1, u3, u4])
    triangles = np.array(triangles)

    node_x = node_x_orig[unique_idx]
    node_y = node_y_orig[unique_idx]
    tri = Triangulation(node_x, node_y, triangles)
    return node_x_orig, node_y_orig, tri, orig_to_unique


def _load_czm_geometry(cfg: dict):
    """Load CZM mesh geometry, deduplicate, return Triangulation.

    Returns (node_x_orig, node_y_orig, tri, orig_to_unique)
    """
    df_nodes = pd.read_csv(_resolve_path(cfg, "czm_nodes.csv"))
    df_elems = pd.read_csv(_resolve_path(cfg, "czm_bulk_elements.csv"))

    node_x_orig = df_nodes["x"].values * 1000
    node_y_orig = df_nodes["y"].values * 1000

    xy_rounded = np.round(np.column_stack([node_x_orig, node_y_orig]), decimals=12)
    _, unique_idx, orig_to_unique = np.unique(
        xy_rounded, axis=0, return_index=True, return_inverse=True
    )

    triangles = []
    for _, row in df_elems.iterrows():
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        u1, u2, u3, u4 = orig_to_unique[n1], orig_to_unique[n2], orig_to_unique[n3], orig_to_unique[n4]
        triangles.append([u1, u2, u3])
        triangles.append([u1, u3, u4])
    triangles = np.array(triangles)

    node_x = node_x_orig[unique_idx]
    node_y = node_y_orig[unique_idx]
    tri = Triangulation(node_x, node_y, triangles)
    return node_x_orig, node_y_orig, tri, orig_to_unique


def _setup_rc():
    """Set matplotlib style."""
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False


def _save_fig(fig, cfg, name):
    """Save figure to output_dir."""
    os.makedirs(cfg["output_dir"], exist_ok=True)
    out = os.path.join(cfg["output_dir"], name)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}")
    return out


# =====================================================================
# Dispatch
# =====================================================================

PLOTTERS: dict[str, callable] = {}  # filled as functions are defined below

def main(cfg: dict | None = None):
    _setup_rc()
    c = cfg or CONFIG
    for plot_type in c["plots"]:
        fn = PLOTTERS.get(plot_type)
        if fn:
            try:
                fn(c)
            except Exception as e:
                print(f"  [!] {plot_type} failed: {e}")
        else:
            print(f"  [!] Unknown plot type: {plot_type}")

if __name__ == "__main__":
    main()
