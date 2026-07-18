# CZM Plot Post-Processing Design

**Date**: 2026-05-26
**Status**: Approved
**Scope**: `param/plot_czm.py` (new unified script)

## Problem

Four independent Python scripts (`csv_temperature_cloud.py`, `csv_displacement_cloud.py`, `csv_current_cloud.py`, `csv_damage_evolution.py`) duplicate code for mesh loading, style setup, and CSV finding. Each script produces one type of plot. Users need to run multiple scripts and manually edit each one for different cycles/times.

## Design

### 1. Single entry point with config dict

Replace 4 scripts with 1 `param/plot_czm.py`. All parameters controlled by a `CONFIG` dict at script top:

```python
CONFIG = {
    "data_dir": "output/csv/czm_study_1",
    "mesh_dir": "output/csv/20",
    "output_dir": "output",

    "cycle": 1,
    "target_time": 1200.0,
    "phase": "discharge",

    "plots": [
        "temperature",
        "displacement",
        "current",
        "damage_evolution",
        "capacity_soh",
        "mesh_geometry",
        "traction_separation",
        "separation_cloud",
    ]
}
```

### 2. Shared internal functions

All extracted into the same script (no external module dependency):

| Function | Purpose |
|----------|---------|
| `_find_csv(csv_dir, name)` | Find CSV in data_dir or mesh_dir |
| `_load_mesh_geometry(csv_dir)` | Load thermal mesh, deduplicate, triangulate |
| `_load_czm_geometry(csv_dir)` | Load CZM mesh, deduplicate, triangulate |
| `_setup_rc()` | Set matplotlib rcParams (Times New Roman, STIX) |

### 3. Eight plot types

Each plot type has a `make_<type>()` function. All take `CONFIG` as parameter.

| Plot | Data Source | Description |
|------|-------------|-------------|
| `temperature` | node_temperature.csv + mesh | Node temperature contourf cloud map |
| `displacement` | node_displacement.csv + czm_mesh | 3-panel: \|u\|, ux, uy contourf |
| `current` | element_currents.csv + mesh | Element branch current scatter cloud map |
| `damage_evolution` | cohesive_damage.csv | 4-panel: D_max(t), D(θ), δ_n(θ), D(t) most damaged element |
| `capacity_soh` | cycle_summary.csv | Capacity and SOH vs cycle number |
| `mesh_geometry` | mesh + czm_mesh CSVs | Thermal nodes, CZM elements, layer topology |
| `traction_separation` | cohesive_damage.csv | T_n vs δ_n phase-space trajectories |
| `separation_cloud` | cohesive_damage.csv + czm_mesh | δ_n(θ) distribution cloud map at target time |

### 4. Main dispatch logic

```python
PLOTTERS = {
    "temperature": make_temperature,
    "displacement": make_displacement,
    "current": make_current,
    "damage_evolution": make_damage_evolution,
    "capacity_soh": make_capacity_soh,
    "mesh_geometry": make_mesh_geometry,
    "traction_separation": make_traction_separation,
    "separation_cloud": make_separation_cloud,
}

if __name__ == "__main__":
    for plot_type in CONFIG["plots"]:
        if plot_type in PLOTTERS:
            PLOTTERS[plot_type](CONFIG)
```

### 5. File structure

| File | Action |
|------|--------|
| `param/plot_czm.py` | Create (new unified script) |

The existing 4 scripts remain untouched (not deleted) for backward compatibility.

### 6. Phase name handling

Note: CZM snapshot data uses `"PHASE_DISCHARGE"`, `"PHASE_CHARGE"`, `"PHASE_REST"` while thermal data uses `"discharge"`, `"charge"`, `"rest1"`, `"rest2"`. The plot functions must use the correct phase name for each data source.

### 7. Output files

Each plot saves to `output/csv_<type>.png`. Example:
- `output/csv_temperature_cloud.png`
- `output/csv_displacement_cloud.png`
- `output/csv_damage_evolution.png`
- etc.

### 8. Dependencies

Only standard scientific Python: numpy, pandas, matplotlib. No additional packages.
