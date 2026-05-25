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
from fig8_single_tab_debonding import make_figure as make_fig8
from fig9_triple_tab_debonding import make_figure as make_fig9
from fig10_damage_cycle_comparison import make_figure as make_fig10
from fig11_capacity_loss_cycle_comparison import make_figure as make_fig11
from fig12_charge_rate_soh_comparison import make_figure as make_fig12


def main() -> None:
    makers = [
        ("fig1_mesh_schematic", make_fig1),
        ("fig2_coupling_framework", make_fig2),
        ("fig3_capacity_soh", make_fig3),
        ("fig4_damage_evolution", make_fig4),
        ("fig5_debonding_evolution", make_fig5),
        ("fig6_temperature_evolution", make_fig6),
        ("fig7_separation_displacement", make_fig7),
        ("fig8_single_tab_debonding", make_fig8),
        ("fig9_triple_tab_debonding", make_fig9),
        ("fig10_damage_cycle_comparison", make_fig10),
        ("fig11_capacity_loss_cycle_comparison", make_fig11),
        ("fig12_charge_rate_soh_comparison", make_fig12),
    ]
    for name, maker in makers:
        path = maker()
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
