#!/usr/bin/env python3
"""Current distribution cloud map from CSV data.

Renders the element-level branch current as a scatter-based cloud map on the
jellyroll cross-section at a specified time step.
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation


def make_figure(
    csv_dir: str | None = None,
    target_time: float = 1200.0,
    cycle: int = 1,
    phase: str = "discharge",
) -> str:
    """Plot element current distribution cloud map.

    Returns path to saved PNG.
    """
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if csv_dir is None:
        csv_dir = os.path.join(root, "output", "csv", "czm_study_1")
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Load data ──
    df_elem = pd.read_csv(os.path.join(csv_dir, "element_currents.csv"))
    df_node = pd.read_csv(os.path.join(csv_dir, "node_temperature.csv"))

    # Filter cycle and phase
    df_e = df_elem[(df_elem["cycle"] == cycle) & (df_elem["phase"] == phase)].copy()
    df_n = df_node[(df_node["cycle"] == cycle) & (df_node["phase"] == phase)].copy()

    # Find closest timestep
    times = df_e["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - target_time))]

    df_e_t = df_e[df_e["time_s"] == t_closest].copy()
    df_n_t = df_n[df_n["time_s"] == t_closest].copy()

    if len(df_e_t) == 0:
        print(f"No data at t={t_closest:.1f}s for cycle={cycle} phase={phase}")
        return ""

    # ── Compute element centers from node coordinates ──
    # Node coordinates (unique per node_id)
    node_coords = df_n_t[["node_id", "x", "y"]].drop_duplicates("node_id")
    node_coords = node_coords.set_index("node_id")

    # Approximate element centers by mapping element index to node positions
    # Since we don't have element-node connectivity in CSV,
    # we use the thermal mesh structure: Q4 elements with sequential numbering.
    # For each element e (1-based), its 4 nodes are approximately:
    #   n1 = 2*e-1, n2 = 2*e, n3 = 2*e+2, n4 = 2*e+1  (strip mesh pattern)
    # But the jellyroll mesh is more complex. Instead, use scatter of node T
    # as background, and overlay element data at estimated centers.

    # Better approach: use node-level data to create triangulation,
    # then map element current to a cell-centered quantity using nearest-neighbor.
    # Simplest: scatter plot of elements colored by current, positioned via
    # interpolation from node coords using element index → center estimate.

    # For jellyroll, each element has a center that lies on the spiral.
    # We can approximate: for element e, its center is at the average position
    # of the 4 surrounding nodes. Since mesh is structured in a strip,
    # we estimate using the midpoint of two adjacent nodes.

    # Use the approach: compute node-based triangulation, then average per element
    # For simplicity, use scatter on node coordinates with T data, and
    # add element current as a separate scatter with interpolated positions.

    # Actually, for the current cloud map, we can use a different strategy:
    # Create a triangulation from node coords and plot temperature-like field,
    # then overlay element current as colored markers.

    # Most robust: use node triangulation and interpolate element centers
    x_n = df_n_t["x"].values * 1000  # mm
    y_n = df_n_t["y"].values * 1000
    T_n = df_n_t["T_K"].values

    # Element data
    I_e = df_e_t["I_e"].values
    T_e = df_e_t["T_e"].values
    elem_ids = df_e_t["elem_id"].values

    # Compute approximate element centers:
    # For Q4 strip mesh, element e spans nodes (2e-1, 2e, 2e+1, 2e+2)
    # (1-based indexing). Center = average of 4 node coords.
    x_ec = np.full(len(elem_ids), np.nan)
    y_ec = np.full(len(elem_ids), np.nan)
    for i, eid in enumerate(elem_ids):
        n1, n2 = int(2 * eid - 1), int(2 * eid)
        if n1 in node_coords.index and n2 in node_coords.index:
            x_ec[i] = 0.5 * (node_coords.loc[n1, "x"] + node_coords.loc[n2, "x"]) * 1000
            y_ec[i] = 0.5 * (node_coords.loc[n1, "y"] + node_coords.loc[n2, "y"]) * 1000

    # Fallback: use triangulation-based interpolation for missing centers
    valid = ~np.isnan(x_ec)
    if valid.sum() < len(elem_ids) * 0.5:
        # If node-based estimate fails, compute from theta (angle from origin)
        # This works for spiral mesh where elements are ordered angularly
        for i, eid in enumerate(elem_ids):
            n1 = int(2 * eid - 1)
            if n1 in node_coords.index:
                x_ec[i] = node_coords.loc[n1, "x"] * 1000
                y_ec[i] = node_coords.loc[n1, "y"] * 1000
            elif eid <= len(node_coords):
                x_ec[i] = node_coords.iloc[min(eid - 1, len(node_coords) - 1)]["x"] * 1000
                y_ec[i] = node_coords.iloc[min(eid - 1, len(node_coords) - 1)]["y"] * 1000

    # ── Plot ──
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=300)

    # Scatter plot of element currents
    valid = ~np.isnan(x_ec)
    if valid.sum() > 0:
        sc = ax.scatter(
            x_ec[valid], y_ec[valid],
            c=I_e[valid],
            cmap="coolwarm",
            s=3.0,
            marker="s",
            edgecolors="none",
            alpha=0.9,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, format="%.3f")
        cbar.set_label("Branch current $I_e$ (A)", fontsize=9)
        cbar.ax.tick_params(labelsize=8)
    else:
        ax.text(
            0.5, 0.5, "No valid element positions",
            transform=ax.transAxes, ha="center", va="center",
        )

    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=9)
    ax.set_ylabel("y (mm)", fontsize=9)
    ax.set_title(
        f"Branch current distribution\nt = {t_closest:.0f} s, {phase}, cycle {cycle}",
        fontsize=10,
    )
    ax.tick_params(labelsize=8)

    out_path = os.path.join(output_dir, "csv_current_cloud.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


if __name__ == "__main__":
    make_figure()
