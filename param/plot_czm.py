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
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

# =====================================================================
# User Configuration
# =====================================================================

CONFIG = {
    "data_dir": "output/csv/synthetic_core_collapse",  # synthetic collapse CSV
    "mesh_dir": "output/csv/80",                       # mesh geometry CSVs
    "output_dir": "output",                  # PNG output directory

    "cycle": 200,                            # match the 200-cycle reference image
    "target_time": 1800.0,                   # target time (s)
    "phase": "discharge",                    # phase name (for thermal/current)
    "deformation_amplification": 1.0,        # show the constructed displacement at true scale

    "plots": [
        "deformed_shape",
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
# Plotter Registry
# =====================================================================

PLOTTERS: dict[str, callable] = {}


# =====================================================================
# 1. Temperature cloud map
# =====================================================================

def make_temperature(cfg: dict):
    node_x_orig, node_y_orig, tri, orig_to_unique = _load_mesh_geometry(cfg)
    n_unique = tri.x.size

    df = pd.read_csv(os.path.join(cfg["data_dir"], "node_temperature.csv"))
    df = df[(df["cycle"] == cfg["cycle"]) & (df["phase"] == cfg["phase"])].copy()

    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    sub = df[df["time_s"] == t_closest].copy()

    if len(sub) == 0:
        print(f"  [temperature] No data at t={t_closest:.1f}s")
        return

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

    vmin, vmax = np.nanmin(T_unique), np.nanmax(T_unique)
    if vmax - vmin < 0.01:
        vmax = vmin + 0.1

    interp = LinearTriInterpolator(tri, T_unique)
    xg = np.linspace(tri.x.min(), tri.x.max(), 400)
    yg = np.linspace(tri.y.min(), tri.y.max(), 400)
    Xg, Yg = np.meshgrid(xg, yg)
    Zg = interp(Xg, Yg)

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=300)
    cf = ax.contourf(Xg, Yg, Zg, levels=50, cmap="hot", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, format="%.1f")
    cbar.set_label("Temperature (K)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=9)
    ax.set_ylabel("y (mm)", fontsize=9)
    ax.set_title(f"Temperature distribution\nt = {t_closest:.0f} s, {cfg['phase']}, cycle {cfg['cycle']}", fontsize=10)
    ax.tick_params(labelsize=8)

    _save_fig(fig, cfg, "csv_temperature_cloud.png")

PLOTTERS["temperature"] = make_temperature


# =====================================================================
# 2. Displacement cloud map
# =====================================================================

def make_displacement(cfg: dict):
    node_x_orig, node_y_orig, tri, orig_to_unique = _load_czm_geometry(cfg)
    n_unique = tri.x.size

    df = pd.read_csv(os.path.join(cfg["data_dir"], "node_displacement.csv"))
    czm_phase = "PHASE_" + cfg["phase"].upper()
    df = df[(df["cycle"] == cfg["cycle"]) & (df["phase"] == czm_phase)].copy()

    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    sub = df[df["time_s"] == t_closest].copy()

    if len(sub) == 0:
        print(f"  [displacement] No data at t={t_closest:.1f}s for phase={czm_phase}")
        return

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

    disp_mag = np.sqrt(np.nan_to_num(ux_unique) ** 2 + np.nan_to_num(uy_unique) ** 2)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), dpi=300, constrained_layout=True)
    titles = [r"$|u|$ ($\mu$m)", r"$u_x$ ($\mu$m)", r"$u_y$ ($\mu$m)"]
    fields = [disp_mag, ux_unique, uy_unique]
    cmaps = ["coolwarm", "coolwarm", "coolwarm"]

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

    fig.suptitle(f"Displacement field  |  t = {t_closest:.0f} s, cycle {cfg['cycle']}", fontsize=11, fontweight="bold")
    _save_fig(fig, cfg, "csv_displacement_cloud.png")

PLOTTERS["displacement"] = make_displacement


# =====================================================================
# 3. Current distribution cloud map
# =====================================================================

def make_current(cfg: dict):
    df_elem = pd.read_csv(os.path.join(cfg["data_dir"], "element_currents.csv"))
    df_mesh_nodes = pd.read_csv(_resolve_path(cfg, "mesh_nodes.csv"))
    df_mesh_elems = pd.read_csv(_resolve_path(cfg, "mesh_elements.csv"))

    df_e = df_elem[(df_elem["cycle"] == cfg["cycle"]) & (df_elem["phase"] == cfg["phase"])].copy()

    times = df_e["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]

    df_e_t = df_e[df_e["time_s"] == t_closest].copy()

    if len(df_e_t) == 0:
        print(f"  [current] No data at t={t_closest:.1f}s")
        return

    # Node coordinates (1-indexed -> 0-indexed)
    x_n = df_mesh_nodes["x"].values * 1000  # m -> mm
    y_n = df_mesh_nodes["y"].values * 1000

    # Build elem_id -> I_e lookup
    I_lookup = dict(zip(df_e_t["elem_id"].values, df_e_t["I_e"].values))

    # Build Q4 polygon patches colored by branch current
    patches = []
    colors = []
    for _, row in df_mesh_elems.iterrows():
        eid = int(row["elem_id"]) if "elem_id" in df_mesh_elems.columns else int(row.name) + 1
        if eid not in I_lookup:
            continue
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        quad = [[x_n[n1], y_n[n1]], [x_n[n2], y_n[n2]], [x_n[n3], y_n[n3]], [x_n[n4], y_n[n4]]]
        patches.append(MplPolygon(quad, closed=True))
        colors.append(I_lookup[eid])

    if not patches:
        print("  [current] No matching elements between mesh and current data")
        return

    colors = np.array(colors)

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=300)
    pc = PatchCollection(patches, cmap="coolwarm", edgecolors="face", linewidths=0.1)
    pc.set_array(colors)
    pc.set_clim(colors.min(), colors.max())
    ax.add_collection(pc)

    cbar = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, format="%.3f")
    cbar.set_label("Branch current $I_e$ (A)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_xlim(x_n.min(), x_n.max())
    ax.set_ylim(y_n.min(), y_n.max())
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=9)
    ax.set_ylabel("y (mm)", fontsize=9)
    ax.set_title(f"Branch current distribution\nt = {t_closest:.0f} s, {cfg['phase']}, cycle {cfg['cycle']}", fontsize=10)
    ax.tick_params(labelsize=8)
    _save_fig(fig, cfg, "csv_current_cloud.png")

PLOTTERS["current"] = make_current


# =====================================================================
# 4. Damage & separation evolution
# =====================================================================

def make_damage_evolution(cfg: dict):
    df = pd.read_csv(os.path.join(cfg["data_dir"], "cohesive_damage.csv"))
    df = df[df["cycle"] == cfg["cycle"]].copy()

    if len(df) == 0:
        print(f"  [damage_evolution] No CZM data for cycle {cfg['cycle']}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=300, constrained_layout=True)

    # (a) D_max over time per phase
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

    # (b) D(theta) profile
    ax = axes[0, 1]
    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    sub_t = df[df["time_s"] == t_closest]
    t_final = times[-1]
    sub_f = df[df["time_s"] == t_final]

    if len(sub_t) > 0:
        sub_sorted = sub_t.sort_values("theta_deg")
        ax.plot(sub_sorted["theta_deg"], sub_sorted["D"], "b-", lw=0.6, alpha=0.8)
        ax.fill_between(sub_sorted["theta_deg"], sub_sorted["D"], alpha=0.15, color="blue")
    if len(sub_f) > 0:
        sub_f_sorted = sub_f.sort_values("theta_deg")
        ax.plot(sub_f_sorted["theta_deg"], sub_f_sorted["D"], "r-", lw=0.8, alpha=0.7, label=f"t={t_final:.0f}s")
    ax.set_xlabel(r"$\theta$ (deg)", fontsize=9)
    ax.set_ylabel("Damage D", fontsize=9)
    ax.set_title(f"(b) Damage profile at t = {t_closest:.0f} s", fontsize=10)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # (c) delta_n(theta)
    ax = axes[1, 0]
    if len(sub_t) > 0:
        sub_sorted = sub_t.sort_values("theta_deg")
        dn_um = sub_sorted["delta_n"].values * 1e6
        ax.plot(sub_sorted["theta_deg"], dn_um, "b-", lw=0.6, alpha=0.8)
        ax.fill_between(sub_sorted["theta_deg"], dn_um, alpha=0.15, color="blue")
    if len(sub_f) > 0:
        sub_f_sorted = sub_f.sort_values("theta_deg")
        dn_f = sub_f_sorted["delta_n"].values * 1e6
        ax.plot(sub_f_sorted["theta_deg"].values, dn_f, "r-", lw=0.8, alpha=0.7, label=f"t={t_final:.0f}s")
    ax.set_xlabel(r"$\theta$ (deg)", fontsize=9)
    ax.set_ylabel(r"Normal separation $\delta_n$ ($\mu$m)", fontsize=9)
    ax.set_title(f"(c) Separation profile at t = {t_closest:.0f} s", fontsize=10)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)

    # (d) D(t) for most damaged element
    ax = axes[1, 1]
    final_time = df["time_s"].max()
    sub_final = df[df["time_s"] == final_time]
    if len(sub_final) > 0:
        max_elem = sub_final.loc[sub_final["D"].idxmax(), "coh_id"]
        elem_data = df[df["coh_id"] == max_elem].sort_values("time_s")
        ax.plot(elem_data["time_s"], elem_data["D"], "r-", lw=1.2)
        ax.fill_between(elem_data["time_s"], elem_data["D"], alpha=0.15, color="red")

        ax2 = ax.twinx()
        dn_elem = elem_data["delta_n"].values * 1e6
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

    fig.suptitle(f"CZM Damage & Separation Evolution  |  Cycle {cfg['cycle']}", fontsize=12, fontweight="bold")
    _save_fig(fig, cfg, "csv_damage_evolution.png")

PLOTTERS["damage_evolution"] = make_damage_evolution


# =====================================================================
# 5. Capacity loss & SOH curve
# =====================================================================

def make_capacity_soh(cfg: dict):
    df = pd.read_csv(os.path.join(cfg["data_dir"], "cycle_summary.csv"))
    df_discharge = df[df["phase"] == "discharge"].copy()

    if len(df_discharge) == 0:
        print("  [capacity_soh] No discharge data in cycle_summary")
        return

    cycles = df_discharge["cycle"].values
    cap = df_discharge["capacity_ah"].values
    soh = df_discharge["soh"].values * 100

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=300, constrained_layout=True)

    # Capacity
    ax1.bar(cycles, cap, color="#1f77b4", alpha=0.8, width=0.6)
    ax1.plot(cycles, cap, "o-", color="#1f77b4", ms=5, lw=1.2)
    ax1.set_xlabel("Cycle", fontsize=9)
    ax1.set_ylabel("Discharge Capacity (Ah)", fontsize=9)
    ax1.set_title("(a) Capacity vs Cycle", fontsize=10)
    ax1.tick_params(labelsize=8)
    ax1.grid(True, alpha=0.3, axis="y")

    # SOH
    ax2.plot(cycles, soh, "s-", color="#d62728", ms=5, lw=1.2)
    ax2.fill_between(cycles, soh, alpha=0.1, color="red")
    ax2.set_xlabel("Cycle", fontsize=9)
    ax2.set_ylabel("SOH (%)", fontsize=9)
    ax2.set_title("(b) State of Health", fontsize=10)
    ax2.tick_params(labelsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(min(soh) - 2, 100.5)

    fig.suptitle("Capacity & SOH Evolution", fontsize=12, fontweight="bold")
    _save_fig(fig, cfg, "csv_capacity_soh.png")

PLOTTERS["capacity_soh"] = make_capacity_soh


# =====================================================================
# 6. Mesh geometry overview
# =====================================================================

def make_mesh_geometry(cfg: dict):
    df_mesh_nodes = pd.read_csv(_resolve_path(cfg, "mesh_nodes.csv"))
    df_mesh_elems = pd.read_csv(_resolve_path(cfg, "mesh_elements.csv"))
    df_czm_nodes = pd.read_csv(_resolve_path(cfg, "czm_nodes.csv"))
    df_czm_elems = pd.read_csv(_resolve_path(cfg, "czm_elements.csv"))
    df_czm_bulk = pd.read_csv(_resolve_path(cfg, "czm_bulk_elements.csv"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=300, constrained_layout=True)

    # (a) Thermal mesh Q4 element outlines
    x_m = df_mesh_nodes["x"].values
    y_m = df_mesh_nodes["y"].values
    for _, row in df_mesh_elems.iterrows():
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        xs = [x_m[n1], x_m[n2], x_m[n3], x_m[n4], x_m[n1]]
        ys = [y_m[n1], y_m[n2], y_m[n3], y_m[n4], y_m[n1]]
        ax1.plot(xs, ys, "k-", lw=0.15)

    ax1.set_aspect("equal")
    ax1.set_xlabel("x (m)", fontsize=9)
    ax1.set_ylabel("y (m)", fontsize=9)
    ax1.set_title("(a) Thermal mesh outline", fontsize=10)
    ax1.tick_params(labelsize=8)

    # (b) CZM bulk mesh outline
    x_c = df_czm_nodes["x"].values * 1000
    y_c = df_czm_nodes["y"].values * 1000
    ax2.scatter(x_c, y_c, s=0.3, c="blue", alpha=0.5, label="CZM nodes")
    for _, row in df_czm_bulk.iterrows():
        n1, n2, n3, n4 = int(row["n1"]) - 1, int(row["n2"]) - 1, int(row["n3"]) - 1, int(row["n4"]) - 1
        xs = [x_c[n1], x_c[n2], x_c[n3], x_c[n4], x_c[n1]]
        ys = [y_c[n1], y_c[n2], y_c[n3], y_c[n4], y_c[n1]]
        ax2.plot(xs, ys, "b-", lw=0.15, alpha=0.3)

    ax2.set_aspect("equal")
    ax2.set_xlabel("x (mm)", fontsize=9)
    ax2.set_ylabel("y (mm)", fontsize=9)
    ax2.set_title("(b) CZM bulk mesh", fontsize=10)
    ax2.tick_params(labelsize=8)
    ax2.legend(fontsize=7, markerscale=10)

    fig.suptitle("Mesh Geometry Overview", fontsize=12, fontweight="bold")
    _save_fig(fig, cfg, "csv_mesh_geometry.png")

PLOTTERS["mesh_geometry"] = make_mesh_geometry


# =====================================================================
# 7. Traction-separation phase diagram
# =====================================================================

def make_traction_separation(cfg: dict):
    df = pd.read_csv(os.path.join(cfg["data_dir"], "cohesive_damage.csv"))
    df = df[df["cycle"] == cfg["cycle"]].copy()

    if len(df) == 0:
        print(f"  [traction_separation] No CZM data for cycle {cfg['cycle']}")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300, constrained_layout=True)

    # (a) T_n vs delta_n for selected elements
    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    t_final = times[-1]

    # Pick a few representative elements (evenly spaced)
    all_ids = sorted(df["coh_id"].unique())
    sample_ids = all_ids[::max(1, len(all_ids) // 8)]

    for eid in sample_ids:
        elem_data = df[df["coh_id"] == eid].sort_values("delta_n")
        dn = elem_data["delta_n"].values * 1e6  # um
        tn = elem_data["T_n"].values * 1e-6     # MPa (from Pa)
        ax1.plot(dn, tn, "-", lw=0.6, alpha=0.7)
        # Mark final state
        if len(dn) > 0:
            ax1.plot(dn[-1], tn[-1], "o", ms=3, alpha=0.8)

    ax1.set_xlabel(r"$\delta_n$ ($\mu$m)", fontsize=9)
    ax1.set_ylabel(r"$T_n$ (MPa)", fontsize=9)
    ax1.set_title("(a) Traction-separation curves", fontsize=10)
    ax1.tick_params(labelsize=8)
    ax1.grid(True, alpha=0.3)

    # (b) All elements at final time: delta_n vs T_n colored by D
    sub_f = df[df["time_s"] == t_final]

    if len(sub_f) > 0:
        sc = ax2.scatter(sub_f["delta_n"] * 1e6, sub_f["T_n"] * 1e-6,
                         c=sub_f["D"], cmap="coolwarm", s=10, edgecolors="none", alpha=0.8)
        fig.colorbar(sc, ax=ax2, fraction=0.046, pad=0.04, label="Damage D")

    ax2.set_xlabel(r"$\delta_n$ ($\mu$m)", fontsize=9)
    ax2.set_ylabel(r"$T_n$ (MPa)", fontsize=9)
    ax2.set_title(f"(b) All elements at t = {t_final:.0f} s", fontsize=10)
    ax2.tick_params(labelsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(f"Traction-Separation Phase Diagram  |  Cycle {cfg['cycle']}", fontsize=12, fontweight="bold")
    _save_fig(fig, cfg, "csv_traction_separation.png")

PLOTTERS["traction_separation"] = make_traction_separation


# =====================================================================
# 8. Separation displacement cloud map
# =====================================================================

def make_separation_cloud(cfg: dict):
    df = pd.read_csv(os.path.join(cfg["data_dir"], "cohesive_damage.csv"))
    df = df[df["cycle"] == cfg["cycle"]].copy()

    if len(df) == 0:
        print(f"  [separation_cloud] No CZM data for cycle {cfg['cycle']}")
        return

    times = df["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    t_final = times[-1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), dpi=300, constrained_layout=True)

    sc = None
    for ax, t_val, label in [(ax1, t_closest, f"t = {t_closest:.0f} s"), (ax2, t_final, f"t = {t_final:.0f} s")]:
        sub = df[df["time_s"] == t_val].sort_values("theta_deg")
        if len(sub) == 0:
            ax.set_title(f"{label} (no data)")
            continue

        theta = sub["theta_deg"].values
        dn = sub["delta_n"].values * 1e6  # um
        D = sub["D"].values

        sc = ax.scatter(theta, dn, c=D, cmap="coolwarm", s=5, edgecolors="none", alpha=0.8,
                        vmin=0, vmax=1)
        ax.fill_between(theta, dn, alpha=0.1, color="blue")

        ax.set_xlabel(r"$\theta$ (deg)", fontsize=9)
        ax.set_ylabel(r"$\delta_n$ ($\mu$m)", fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)

    if sc is not None:
        fig.colorbar(sc, ax=ax2, fraction=0.046, pad=0.04, label="Damage D")
    fig.suptitle(f"Separation Displacement Profile  |  Cycle {cfg['cycle']}", fontsize=12, fontweight="bold")
    _save_fig(fig, cfg, "csv_separation_cloud.png")

PLOTTERS["separation_cloud"] = make_separation_cloud


# =====================================================================
# 9. Deformed shape visualization
# =====================================================================

def make_deformed_shape(cfg: dict):
    """Plot full and zoomed deformed Q4 meshes without a scalar cloud map."""
    df_nodes = pd.read_csv(_resolve_path(cfg, "czm_nodes.csv"))
    df_elems = pd.read_csv(_resolve_path(cfg, "czm_bulk_elements.csv"))

    df_disp = pd.read_csv(os.path.join(cfg["data_dir"], "node_displacement.csv"))
    czm_phase = "PHASE_" + cfg["phase"].upper()
    df_disp = df_disp[(df_disp["cycle"] == cfg["cycle"]) & (df_disp["phase"] == czm_phase)].copy()

    if len(df_disp) == 0:
        print(f"  [deformed_shape] No displacement data for cycle {cfg['cycle']}, phase {czm_phase}")
        return

    times = df_disp["time_s"].unique()
    t_closest = times[np.argmin(np.abs(times - cfg["target_time"]))]
    sub = df_disp[df_disp["time_s"] == t_closest].copy()

    if len(sub) == 0:
        print(f"  [deformed_shape] No data at t={t_closest:.1f}s")
        return

    # Original coordinates (m) -> mm
    x0 = df_nodes["x"].values
    y0 = df_nodes["y"].values
    s = 1000
    x0_mm, y0_mm = x0 * s, y0 * s

    # Build displacement arrays (m)
    ux = np.zeros(len(x0))
    uy = np.zeros(len(y0))
    for _, row in sub.iterrows():
        nid = int(row["node_id"]) - 1
        ux[nid] = row["ux"]
        uy[nid] = row["uy"]

    amp = float(cfg.get("deformation_amplification", 1.0))
    if not np.isfinite(amp) or amp <= 0:
        raise ValueError("deformation_amplification must be finite and positive")

    # Deformed coordinates (mm)
    x_def_mm = (x0 + amp * ux) * s
    y_def_mm = (y0 + amp * uy) * s

    nn = df_elems[["n1", "n2", "n3", "n4"]].values.astype(int) - 1

    fig, (ax_full, ax_zoom) = plt.subplots(
        1, 2, figsize=(11, 5.2), dpi=300, constrained_layout=True
    )

    def draw_mesh(ax, x_nodes, y_nodes, *, color, linewidth, alpha):
        for n1, n2, n3, n4 in nn:
            xs = [x_nodes[n1], x_nodes[n2], x_nodes[n3], x_nodes[n4], x_nodes[n1]]
            ys = [y_nodes[n1], y_nodes[n2], y_nodes[n3], y_nodes[n4], y_nodes[n1]]
            ax.plot(xs, ys, color=color, lw=linewidth, alpha=alpha)

    for ax in (ax_full, ax_zoom):
        draw_mesh(ax, x0_mm, y0_mm, color="#b8b8b8", linewidth=0.12, alpha=0.55)
        draw_mesh(ax, x_def_mm, y_def_mm, color="#171717", linewidth=0.22, alpha=0.95)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)", fontsize=9)
        ax.set_ylabel("y (mm)", fontsize=9)
        ax.tick_params(labelsize=8)

    all_x = np.concatenate([x0_mm, x_def_mm])
    all_y = np.concatenate([y0_mm, y_def_mm])
    x_pad = 0.04 * np.ptp(all_x)
    y_pad = 0.04 * np.ptp(all_y)
    ax_full.set_xlim(all_x.min() - x_pad, all_x.max() + x_pad)
    ax_full.set_ylim(all_y.min() - y_pad, all_y.max() + y_pad)
    ax_full.set_title("(a) Full deformed jellyroll mesh", fontsize=10)
    ax_full.plot([], [], color="#b8b8b8", lw=1.0, label="Original mesh")
    ax_full.plot([], [], color="#171717", lw=1.0, label="Deformed mesh")
    ax_full.legend(loc="lower left", fontsize=8, frameon=False)

    # Center the detail view on the node with the largest displacement and
    # shift slightly outward so the alternating folded layers remain visible.
    displacement_magnitude = np.hypot(ux, uy)
    radius = np.hypot(x0, y0)
    radial_fraction = (radius - radius.min()) / (radius.max() - radius.min())
    fold_candidates = np.where((radial_fraction >= 0.08) & (radial_fraction <= 0.25))[0]
    fold_node = int(fold_candidates[np.argmax(displacement_magnitude[fold_candidates])])
    model_span_mm = max(np.ptp(x0_mm), np.ptp(y0_mm))
    zoom_half_width = 0.16 * model_span_mm
    zoom_center_x = x0_mm[fold_node] + 0.15 * zoom_half_width
    zoom_center_y = y0_mm[fold_node]
    ax_zoom.set_xlim(zoom_center_x - zoom_half_width, zoom_center_x + zoom_half_width)
    ax_zoom.set_ylim(
        zoom_center_y - 0.82 * zoom_half_width,
        zoom_center_y + 0.82 * zoom_half_width,
    )
    ax_zoom.set_title("(b) Inner-core collapse folds", fontsize=10)
    ax_zoom.annotate(
        "Fold zone",
        xy=(x_def_mm[fold_node], y_def_mm[fold_node]),
        xytext=(
            zoom_center_x + 0.48 * zoom_half_width,
            zoom_center_y + 0.58 * zoom_half_width,
        ),
        color="#d62728",
        fontsize=8,
        ha="center",
        arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2),
    )

    fig.suptitle(
        f"Synthetic core-collapse deformation  (amplification ×{amp:.1f})\n"
        f"t = {t_closest:.0f} s, {cfg['phase']}, cycle {cfg['cycle']}",
        fontsize=11,
        fontweight="bold",
    )

    _save_fig(fig, cfg, "csv_deformed_shape.png")

PLOTTERS["deformed_shape"] = make_deformed_shape


# =====================================================================
# Dispatch
# =====================================================================

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
