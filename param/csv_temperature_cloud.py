#!/usr/bin/env python3
"""Temperature cloud map from CSV data.

Renders node-level temperature as a mesh-based contour map on the
jellyroll cross-section at a specified time step.

Uses mesh_nodes.csv + mesh_elements.csv for geometry (Q4 → triangles)
and node_temperature.csv for field data.
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator


def _find_csv(csv_dir: str, name: str) -> str:
    """Find CSV in csv_dir or its parent directory."""
    p = os.path.join(csv_dir, name)
    if os.path.isfile(p):
        return p
    parent = os.path.dirname(csv_dir)
    p2 = os.path.join(parent, name)
    if os.path.isfile(p2):
        return p2
    raise FileNotFoundError(f"{name} not found in {csv_dir} or its parent")


def _load_mesh_geometry(csv_dir: str):
    """Load mesh geometry, deduplicate overlapping spiral nodes, return Triangulation.

    The jellyroll mesh has duplicate coordinates at spiral wrap-around points
    (inner/outer spiral overlap). We merge them before building triangulation.

    Returns
    -------
    node_x, node_y : original coordinate arrays (length N_original)
    tri : Triangulation on deduplicated nodes
    orig_to_unique : mapping from original node index → unique node index
    """
    df_nodes = pd.read_csv(_find_csv(csv_dir, "mesh_nodes.csv"))
    df_elems = pd.read_csv(_find_csv(csv_dir, "mesh_elements.csv"))

    node_x_orig = df_nodes["x"].values * 1000  # m -> mm
    node_y_orig = df_nodes["y"].values * 1000

    # Deduplicate: round to ~12 significant digits, find unique positions
    xy_rounded = np.round(np.column_stack([node_x_orig, node_y_orig]), decimals=12)
    _, unique_idx, orig_to_unique = np.unique(
        xy_rounded, axis=0, return_index=True, return_inverse=True
    )

    # Remap element connectivity: split Q4 → 2 triangles, use unique indices
    triangles = []
    for _, row in df_elems.iterrows():
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        u1, u2, u3, u4 = orig_to_unique[n1], orig_to_unique[n2], orig_to_unique[n3], orig_to_unique[n4]
        triangles.append([u1, u2, u3])
        triangles.append([u1, u3, u4])
    triangles = np.array(triangles)

    # Deduplicated coordinates
    node_x = node_x_orig[unique_idx]
    node_y = node_y_orig[unique_idx]

    tri = Triangulation(node_x, node_y, triangles)
    return node_x_orig, node_y_orig, tri, orig_to_unique


def make_figure(
    csv_dir: str | None = None,
    target_time: float = 1200.0,
    cycle: int = 1,
    phase: str = "discharge",
) -> str:
    """Plot node temperature cloud map. Returns path to saved PNG."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if csv_dir is None:
        csv_dir = os.path.join(root, "output", "csv", "czm_study_1")
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Load mesh geometry ──
    node_x_orig, node_y_orig, tri, orig_to_unique = _load_mesh_geometry(csv_dir)
    n_unique = tri.x.size

    # ── Load node temperature ──
    df = pd.read_csv(os.path.join(csv_dir, "node_temperature.csv"))
    df = df[(df["cycle"] == cycle) & (df["phase"] == phase)].copy()

    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - target_time))]
    sub = df[df["time_s"] == t_closest].copy()

    if len(sub) == 0:
        print(f"No data at t={t_closest:.1f}s")
        return ""

    # Map temperature to original node IDs, then to deduplicated nodes
    T_orig = np.full(len(node_x_orig), np.nan)
    for _, row in sub.iterrows():
        nid = int(row["node_id"]) - 1
        T_orig[nid] = row["T_K"]

    T_unique = np.full(n_unique, np.nan)
    for i in range(len(T_orig)):
        u = orig_to_unique[i]
        if np.isnan(T_unique[u]):
            T_unique[u] = T_orig[i]
        elif not np.isnan(T_orig[i]):
            T_unique[u] = 0.5 * (T_unique[u] + T_orig[i])

    if np.all(np.isnan(T_unique)):
        print("No valid temperature data")
        return ""

    vmin, vmax = np.nanmin(T_unique), np.nanmax(T_unique)
    if vmax - vmin < 0.01:
        vmax = vmin + 0.1

    # ── Interpolation grid ──
    interp = LinearTriInterpolator(tri, T_unique)
    xg = np.linspace(tri.x.min(), tri.x.max(), 400)
    yg = np.linspace(tri.y.min(), tri.y.max(), 400)
    Xg, Yg = np.meshgrid(xg, yg)
    Zg = interp(Xg, Yg)

    # ── Plot ──
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=300)
    cf = ax.contourf(Xg, Yg, Zg, levels=50, cmap="hot", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, format="%.1f")
    cbar.set_label("Temperature (K)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=9)
    ax.set_ylabel("y (mm)", fontsize=9)
    ax.set_title(
        f"Temperature distribution\nt = {t_closest:.0f} s, {phase}, cycle {cycle}",
        fontsize=10,
    )
    ax.tick_params(labelsize=8)

    out_path = os.path.join(output_dir, "csv_temperature_cloud.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


if __name__ == "__main__":
    make_figure()
