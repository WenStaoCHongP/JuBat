#!/usr/bin/env python3
"""Displacement cloud map from CSV data.

Renders node-level displacement magnitude as a mesh-based contour map on the
jellyroll cross-section at a specified time step.

Uses czm_nodes.csv + czm_bulk_elements.csv for geometry (Q4 → triangles)
and node_displacement.csv for field data.
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


def _load_czm_geometry(csv_dir: str):
    """Load CZM mesh geometry, deduplicate overlapping spiral nodes.

    Returns
    -------
    node_x_orig, node_y_orig : original coordinate arrays
    tri : Triangulation on deduplicated nodes
    orig_to_unique : mapping from original node index → unique node index
    """
    df_nodes = pd.read_csv(_find_csv(csv_dir, "czm_nodes.csv"))
    df_elems = pd.read_csv(_find_csv(csv_dir, "czm_bulk_elements.csv"))

    node_x_orig = df_nodes["x"].values * 1000  # m -> mm
    node_y_orig = df_nodes["y"].values * 1000

    # Deduplicate overlapping spiral coordinates
    xy_rounded = np.round(np.column_stack([node_x_orig, node_y_orig]), decimals=12)
    _, unique_idx, orig_to_unique = np.unique(
        xy_rounded, axis=0, return_index=True, return_inverse=True
    )

    # Split Q4 bulk elements → 2 triangles, remap to unique indices
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


def make_figure(
    csv_dir: str | None = None,
    target_time: float = 1200.0,
    cycle: int = 1,
    phase: str = "PHASE_DISCHARGE",
) -> str:
    """Plot node displacement cloud map. Returns path to saved PNG."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if csv_dir is None:
        csv_dir = os.path.join(root, "output", "csv", "czm_study_1")
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ── Load CZM mesh geometry ──
    node_x_orig, node_y_orig, tri, orig_to_unique = _load_czm_geometry(csv_dir)
    n_unique = tri.x.size

    # ── Load node displacement ──
    df = pd.read_csv(os.path.join(csv_dir, "node_displacement.csv"))
    df = df[(df["cycle"] == cycle) & (df["phase"] == phase)].copy()

    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - target_time))]
    sub = df[df["time_s"] == t_closest].copy()

    if len(sub) == 0:
        print(f"No data at t={t_closest:.1f}s for phase={phase}")
        return ""

    # Map displacement to original node IDs, then to deduplicated nodes
    ux_orig = np.full(len(node_x_orig), np.nan)
    uy_orig = np.full(len(node_x_orig), np.nan)
    for _, row in sub.iterrows():
        nid = int(row["node_id"]) - 1
        ux_orig[nid] = row["ux"] * 1e6  # m -> um
        uy_orig[nid] = row["uy"] * 1e6

    ux_unique = np.full(n_unique, np.nan)
    uy_unique = np.full(n_unique, np.nan)
    for i in range(len(ux_orig)):
        u = orig_to_unique[i]
        if np.isnan(ux_unique[u]):
            ux_unique[u] = ux_orig[i]
            uy_unique[u] = uy_orig[i]
        elif not np.isnan(ux_orig[i]):
            ux_unique[u] = 0.5 * (ux_unique[u] + ux_orig[i])
            uy_unique[u] = 0.5 * (uy_unique[u] + uy_orig[i])

    if np.all(np.isnan(ux_unique)):
        print("No valid displacement data")
        return ""

    disp_mag = np.sqrt(np.nan_to_num(ux_unique) ** 2 + np.nan_to_num(uy_unique) ** 2)

    # ── Plot ──
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), dpi=300, constrained_layout=True)

    titles = [r"$|u|$ ($\mu$m)", r"$u_x$ ($\mu$m)", r"$u_y$ ($\mu$m)"]
    fields = [disp_mag, ux_unique, uy_unique]
    cmaps = ["inferno", "coolwarm", "coolwarm"]

    for ax, field, title, cmap in zip(axes, fields, titles, cmaps):
        vmin, vmax = np.nanmin(field), np.nanmax(field)
        if vmax - vmin < 1e-12:
            vmax = vmin + 1e-3

        interp = LinearTriInterpolator(tri, field)
        xg = np.linspace(tri.x.min(), tri.x.max(), 400)
        yg = np.linspace(tri.y.min(), tri.y.max(), 400)
        Xg, Yg = np.meshgrid(xg, yg)
        Zg = interp(Xg, Yg)

        cf = ax.contourf(Xg, Yg, Zg, levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, format="%.2f")

        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)", fontsize=9)
        ax.set_ylabel("y (mm)", fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.tick_params(labelsize=8)

    fig.suptitle(
        f"Displacement field  |  t = {t_closest:.0f} s, cycle {cycle}",
        fontsize=11, fontweight="bold",
    )

    out_path = os.path.join(output_dir, "csv_displacement_cloud.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


if __name__ == "__main__":
    make_figure()
