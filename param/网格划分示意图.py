#!/usr/bin/env python3
"""Generate a compact top-view jellyroll thermal-mesh schematic for publication."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import PolyCollection


LAYER_SEQUENCE = [
    ("PE", 0.16, "#FF7A00"),
    ("PCC", 0.08, "#B64A1F"),
    ("PE", 0.16, "#FF7A00"),
    ("SP", 0.12, "#B3B3B3"),
    ("NE", 0.18, "#0077CC"),
    ("NCC", 0.08, "#3E3E3E"),
    ("NE", 0.18, "#0077CC"),
    ("SP", 0.12, "#B3B3B3"),
]

MATERIAL_COLORS = {
    "PCC": "#B64A1F",
    "PE": "#FF7A00",
    "SP": "#B3B3B3",
    "NE": "#0077CC",
    "NCC": "#3E3E3E",
}


def polar_point(radius: float, theta: float, center: tuple[float, float] = (0.0, 0.0)) -> tuple[float, float]:
    return center[0] + radius * np.cos(theta), center[1] + radius * np.sin(theta)


def spiral_curve(base_radius: float, growth: float, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radius = base_radius + growth * theta
    return radius * np.cos(theta), radius * np.sin(theta)


def draw_boundary_annotations(ax, inner_radius: float, outer_radius: float) -> None:
    convection_color = "#56B4E9"
    outer_thetas = np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False)
    arrow_length = 0.11
    for theta in outer_thetas:
        x0, y0 = polar_point(outer_radius, theta)
        x1, y1 = polar_point(outer_radius + arrow_length, theta)
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=convection_color, shrinkA=0.0, shrinkB=0.0),
            zorder=6,
        )


def add_legends(ax_leg) -> None:
    convection_color = "#56B4E9"
    bc_handles = [
        Line2D([0], [0], color=convection_color, lw=1.2, marker=r"$\rightarrow$", markersize=10, markevery=[1]),
        Line2D([0], [0], color="#2A9D5B", lw=2.6),
    ]
    bc_labels = [
        "Outer boundary: convection",
        "Inner boundary: adiabatic",
    ]
    boundary_legend = ax_leg.legend(
        bc_handles,
        bc_labels,
        loc="upper left",
        bbox_to_anchor=(0.03, 0.98),
        ncol=1,
        frameon=True,
        framealpha=0.96,
        facecolor="white",
        edgecolor="#D9D9D9",
        fontsize=8.2,
        handlelength=1.8,
        borderpad=0.45,
        labelspacing=0.35,
    )
    boundary_legend.get_frame().set_linewidth(0.8)
    ax_leg.add_artist(boundary_legend)

    handles = [
        Rectangle((0.0, 0.0), 1.0, 1.0, facecolor=color, edgecolor="none")
        for _, color in MATERIAL_COLORS.items()
    ]
    labels = list(MATERIAL_COLORS.keys())
    material_legend = ax_leg.legend(
        handles,
        labels,
        loc="lower left",
        bbox_to_anchor=(0.06, 0.03),
        ncol=1,
        frameon=True,
        framealpha=0.96,
        facecolor="white",
        edgecolor="#D9D9D9",
        fontsize=8.5,
        borderpad=0.45,
        labelspacing=0.3,
        handlelength=1.0,
    )
    material_legend.get_frame().set_linewidth(0.8)


def add_material_legend(ax) -> None:
    # kept for backward compatibility; legend rendering is now handled by add_legends
    return None


def build_layer_polygons(r_inner: float, growth: float, turns: float, segments_per_turn: int) -> tuple[list[list[tuple[float, float]]], list[str], float]:
    theta = np.linspace(0.0, 2.0 * np.pi * turns, int(turns * segments_per_turn) + 1)
    repeat_thickness = 0.22
    polygons: list[list[tuple[float, float]]] = []
    colors: list[str] = []
    cumulative = np.concatenate(([0.0], np.cumsum([item[1] for item in LAYER_SEQUENCE])))
    total_fraction = cumulative[-1]

    for idx in range(len(theta) - 1):
        theta_a = theta[idx]
        theta_b = theta[idx + 1]
        for layer_index, (_, _, color) in enumerate(LAYER_SEQUENCE):
            s0 = repeat_thickness * cumulative[layer_index] / total_fraction
            s1 = repeat_thickness * cumulative[layer_index + 1] / total_fraction
            x00, y00 = spiral_curve(r_inner + s0, growth, np.array([theta_a]))
            x10, y10 = spiral_curve(r_inner + s1, growth, np.array([theta_a]))
            x11, y11 = spiral_curve(r_inner + s1, growth, np.array([theta_b]))
            x01, y01 = spiral_curve(r_inner + s0, growth, np.array([theta_b]))
            polygons.append([
                (x00[0], y00[0]),
                (x10[0], y10[0]),
                (x11[0], y11[0]),
                (x01[0], y01[0]),
            ])
            colors.append(color)

    return polygons, colors, repeat_thickness


def draw_spiral_mesh(ax, r_inner: float = 0.34, turns: float = 3.0, segments_per_turn: int = 14) -> tuple[float, float]:
    repeat_thickness = 0.22
    growth = repeat_thickness / (2.0 * np.pi)
    polygons, colors, repeat_thickness = build_layer_polygons(r_inner, growth, turns, segments_per_turn)
    collection = PolyCollection(
        polygons,
        facecolors=colors,
        edgecolors="none",
        linewidths=0.0,
        antialiaseds=False,
        zorder=2,
    )
    ax.add_collection(collection)

    theta_dense = np.linspace(0.0, 2.0 * np.pi * turns, 900)
    cumulative = np.concatenate(([0.0], np.cumsum([item[1] for item in LAYER_SEQUENCE])))
    total_fraction = cumulative[-1]

    # Draw cross-band mesh lines once per theta station.
    theta_nodes = np.linspace(0.0, 2.0 * np.pi * turns, int(turns * segments_per_turn) + 1)
    for theta in theta_nodes:
        x0, y0 = polar_point(r_inner + growth * theta, theta)
        x1, y1 = polar_point(r_inner + growth * theta + repeat_thickness, theta)
        ax.plot([x0, x1], [y0, y1], color="#5A5A5A", lw=0.42, alpha=0.9, zorder=3)

    outer_radius = r_inner + growth * theta_dense[-1] + repeat_thickness
    outer_bc_radius = outer_radius
    inner_bc_radius = r_inner
    ax.add_patch(Circle((0.0, 0.0), outer_bc_radius, fill=False, ec="#1F1F1F", lw=1.0, zorder=1))
    ax.add_patch(Circle((0.0, 0.0), inner_bc_radius, fill=False, ec="#2A9D5B", lw=2.4, zorder=5))
    return inner_bc_radius, outer_bc_radius


def figure_topview_thermal_mesh(save: bool = True):
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 9,
            "axes.linewidth": 0.8,
        }
    )

    fig = plt.figure(figsize=(7.8, 5.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[4.6, 1.8], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])

    ax_leg.axis("off")
    ax_leg.set_xlim(0.0, 1.0)
    ax_leg.set_ylim(0.0, 1.0)

    inner_radius, visible_outer_radius = draw_spiral_mesh(ax)
    draw_boundary_annotations(ax, inner_radius=inner_radius, outer_radius=visible_outer_radius)
    add_legends(ax_leg)

    # Include arrow length margin to avoid clipping any convection arrows.
    limit = visible_outer_radius + 0.12
    ax.set_aspect("equal")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axis("off")

    if save:
        for file_name in ("figure_mesh_topview.png", "figure_mesh_topview.svg"):
            fig.savefig(file_name, dpi=600, bbox_inches="tight", pad_inches=0.02)
    return fig


def main() -> None:
    figure_topview_thermal_mesh(save=True)
    print("Saved: figure_mesh_topview.png, figure_mesh_topview.svg")


if __name__ == "__main__":
    main()