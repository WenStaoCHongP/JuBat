#!/usr/bin/env python3
"""Generate all proposal expected-result figures."""

from __future__ import annotations

from fig1_mesh_schematic import make_figure as make_fig1
from fig2_coupling_framework import make_figure as make_fig2
from fig3_capacity_soh import make_figure as make_fig3
from fig4_damage_evolution import make_figure as make_fig4
from fig5_debonding_evolution import make_figure as make_fig5
from fig6_temperature_evolution import make_figure as make_fig6
from fig7_separation_displacement import make_figure as make_fig7


def main() -> None:
    makers = [
        ("fig1_mesh_schematic", make_fig1),
        ("fig2_coupling_framework", make_fig2),
        ("fig3_capacity_soh", make_fig3),
        ("fig4_damage_evolution", make_fig4),
        ("fig5_debonding_evolution", make_fig5),
        ("fig6_temperature_evolution", make_fig6),
        ("fig7_separation_displacement", make_fig7),
    ]
    for name, maker in makers:
        path = maker()
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
