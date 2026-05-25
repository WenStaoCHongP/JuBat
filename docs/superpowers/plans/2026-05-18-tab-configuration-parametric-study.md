# Tab Configuration Parametric Study — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create the technical documentation (`md/15`), Julia simulation script (7 tab cases), and Python visualization scripts (fig_tab1–4) for orthogonal tab configuration parametric study.

**Architecture:** Single Julia script runs 7 independent cases (3 groups: count effect, position effect, spacing effect), collects results into a Dict, and saves scalar metrics as CSV. Four separate Python scripts (one per figure) currently use placeholder data (to be replaced with CSV-loaded simulation results after Julia run). The `md/15` document describes the methodology and serves as the canonical reference.

> **Note on single-discharge limitation:** The initial script runs single 1C discharge (`opt.time = [0, 3600]`). In this regime, `n_fractured` may be zero and `D_max` differences will be subtle (3rd–4th decimal). The SOH metric is approximated as `1 - D_mean` (proxy), since `"soh"` is only available via `CycleSolver`. For meaningful `n_fractured` comparison, a multi-cycle variant should be added later.

**Tech Stack:** Julia (JuBat framework, Plots.jl, CSV output), Python 3 (matplotlib, numpy, pandas for CSV)

**Spec:** `docs/superpowers/specs/2026-05-18-tab-configuration-parametric-study-design.md`

---

## File Structure

| Action | Path | Responsibility |
|--------|------|----------------|
| Create | `md/15_参数研究_极耳工况.md` | Technical documentation for tab parametric study |
| Create | `example/极耳参数研究/tab_parametric_study.jl` | Main Julia script: angle calculation + 7-case loop + result collection |
| Create | `param/fig_tab1_count_effect.py` | Python plot: tab count effect (Group 1, Cases 1–3) |
| Create | `param/fig_tab2_position_effect.py` | Python plot: tab position effect (Group 2, Cases 4–6) |
| Create | `param/fig_tab3_spacing_effect.py` | Python plot: tab spacing effect (Group 3, Cases 7–9) |
| Create | `param/fig_tab4_summary_bar.py` | Python plot: cross-group summary bar chart |
| Modify | `param/generate_all_figures.py` | Register fig_tab1–4 in the master figure generator |
| Modify | `CLAUDE.md` | Add md/15 to documentation index and example table |

---

## Chunk 1: Technical Documentation

### Task 1: Create `md/15_参数研究_极耳工况.md`

**Files:**
- Create: `md/15_参数研究_极耳工况.md`
- Reference: `md/06_内聚力模型_CZM.md` (document style), `md/14_粘性正则化.md` (latest style), spec document

- [ ] **Step 1: Write the technical document**

Create `md/15_参数研究_极耳工况.md` following the established `md/` conventions (numbered sections, markdown tables, LaTeX math where needed). The document structure:

```markdown
# 15 参数研究：极耳工况设计

## 1. 研究背景
极耳是电池与外部电路连接的关键部件，其数量和位置直接影响：
- 局部热分布（极耳冷却效果）
- 电流密度分布
- 热-力学耦合应力
- CZM 界面损伤起始与传播

## 2. 控制变量法
### 2.1 固定参数
[Table from spec §2.1 + §2.2, all fixed parameters]

### 2.2 变量
- 正极耳数量 n_pos: 1, 2, 3
- 正极耳周向位置 theta_pos: 向量 [rad]

## 3. 工况矩阵
### 3.1 第一组：极耳数量效应
[Table: Case 1, 2, 3 — fixed equal-spacing, varying count]
对比分析：Case 1 vs 2 vs 3

### 3.2 第二组：极耳位置效应
[Table: Case 4, 5, 6 — fixed count=1, varying position]
对比分析：Case 4 vs 5 vs 6

### 3.3 第三组：双极耳间距效应
[Table: Case 7, 8, 9 — fixed count=2, varying spacing]
对比分析：Case 7 vs 8 vs 9

### 3.4 工况交叉关系
独立工况 7 个，重复 Case 1=5, Case 2=8

## 4. 角度参数计算
[Angle computation formulas from spec §4, including normalization note]
[Julia code snippet for angle extraction from mesh]

## 5. 输出指标
### 5.1 CZM 损伤指标
### 5.2 热场指标
### 5.3 电化学性能指标
[Tables from spec §5]

## 6. 分析方法
### 6.1 组内对比
### 6.2 跨组效应分析
### 6.3 可视化规范
[From spec §6]
```

Write each section in full, following the spec document content. Use the same tone and format as `md/14_粘性正则化.md` (problem-oriented, practical). All tables carry exact values from the spec.

- [ ] **Step 2: Verify document renders correctly**

Read back the file and verify:
- All section numbers are sequential (1–6)
- Tables have correct column alignment
- Julia code blocks have correct syntax highlighting markers
- No broken markdown

- [ ] **Step 3: Commit**

```bash
git add md/15_参数研究_极耳工况.md
git commit -m "docs: add tab parametric study technical document (md/15)"
```

---

## Chunk 2: Julia Simulation Script

### Task 2: Create the main Julia parametric study script

**Files:**
- Create: `example/极耳参数研究/tab_parametric_study.jl`
- Reference: `example/网格敏感性/4_czm_mesh_sensitivity.jl` (loop pattern), spec §4.1 (angle computation), spec §7 (config template)

- [ ] **Step 1: Create directory and write the script**

Create `example/极耳参数研究/tab_parametric_study.jl`. The script structure follows `4_czm_mesh_sensitivity.jl`:

```julia
"""
Script: Tab Configuration Parametric Study

Orthogonal experimental design with 3 groups (7 independent cases):
  Group 1 (count effect):    Cases 1, 2, 3
  Group 2 (position effect): Cases 4, 5, 6   (5 = 1)
  Group 3 (spacing effect):  Cases 7, 8, 9   (8 = 2)

Fixed: theta_neg = [44π], tab dimensions, solver config.
Variable: theta_pos (positive tab count and angular position).

Output: JLD2 file with all case results + Julia Plots summary.
"""

using Printf, Plots, Statistics

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

# =====================================================================
# Angle computation (from spec §4.1 — mesh-based)
# =====================================================================

"""Compute tab angle parameters from mesh node coordinates.

Builds a reference mesh once to extract the inner spiral cumulative angle range,
then derives case-specific angle markers. This matches the approach used by
`jellyroll_tab_node_indices` in `src/Jellyrollmodel.jl`.
"""
function compute_tab_angles(param_dim; nθ::Int=80, gsorder::Int=2)
    # Build a reference mesh to get actual node positions
    param = JuBat.NormaliseParam(param_dim)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=nθ, gsorder=gsorder)
    mesh = mesh_data.mesh
    nn = size(mesh.node, 1)

    a = param.cell.Rin
    b = param.cell.layer / (2 * pi)

    # Compute inner spiral cumulative angle (same formula as jellyroll_tab_node_indices)
    theta_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
    theta_min_in, theta_max_in = extrema(theta_cum_in)

    range_in = theta_max_in - theta_min_in
    θ_mid   = (theta_min_in + theta_max_in) / 2
    Δθ      = range_in / 4
    θ_start = theta_min_in + range_in * 0.1
    θ_end   = theta_min_in + range_in * 0.9
    θ_1     = theta_min_in + range_in / 6
    θ_3     = theta_max_in - range_in / 6

    return (θ_mid=θ_mid, Δθ=Δθ, θ_start=θ_start, θ_end=θ_end, θ_1=θ_1, θ_3=θ_3,
            theta_min=theta_min_in, theta_max=theta_max_in)
end

# =====================================================================
# Case definitions
# =====================================================================

"""Define all 7 independent cases. Returns Vector of (case_id, label, group, theta_pos)."""
function define_cases(angles)
    θ = angles
    cases = [
        # (id, label, group, theta_pos)
        (1, "1-tab mid (baseline)", "count", [θ.θ_mid]),
        (2, "2-tab equal-spaced",   "count", [θ.θ_mid - θ.Δθ, θ.θ_mid + θ.Δθ]),
        (3, "3-tab equal-spaced",   "count", [θ.θ_1, θ.θ_mid, θ.θ_3]),
        (4, "1-tab start",          "pos",   [θ.θ_start]),
        (6, "1-tab end",            "pos",   [θ.θ_end]),
        (7, "2-tab tight",          "space", [θ.θ_mid - θ.Δθ/2, θ.θ_mid + θ.Δθ/2]),
        (9, "2-tab wide",           "space", [θ.θ_start, θ.θ_end]),
    ]
    return cases
end

# =====================================================================
# Run single case
# =====================================================================

function run_tab_case(param_dim, theta_pos; nθ::Int=80)
    # Override tab configuration
    param_dim.tab.theta_pos = theta_pos
    param_dim.tab.theta_neg = [44π]

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "tab"
    opt.per_element_spme = true

    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
    case.czm_mesh = czm_mesh

    result = JuBat.Solve(case)
    return result, mesh_data
end

# =====================================================================
# Result extraction
# =====================================================================

function extract_metrics(result)
    m = Dict{String, Float64}()
    m["D_max_end"]  = get(result, "czm D_max",  [NaN])[end]
    m["D_mean_end"] = get(result, "czm D_mean", [NaN])[end]
    m["n_frac_end"] = get(result, "czm n_fractured", [NaN])[end]
    # SOH proxy: not available in single-discharge result dict (only in CycleSolver path)
    # Use 1 - D_mean as structural integrity proxy
    D_mean = m["D_mean_end"]
    m["soh_proxy"] = isnan(D_mean) ? NaN : 1.0 - D_mean

    T_field = get(result, "thermal2D temperature [K]", nothing)
    if T_field !== nothing
        T_end = T_field[:, end]
        m["T_max_end"] = maximum(T_end)
        m["T_min_end"] = minimum(T_end)
        m["T_mean_end"] = mean(T_end)
        m["ΔT_end"] = maximum(T_end) - minimum(T_end)
    else
        m["T_max_end"] = m["T_min_end"] = m["T_mean_end"] = m["ΔT_end"] = NaN
    end

    m["V_end"] = get(result, "cell voltage [V]", [NaN])[end]
    return m
end

# =====================================================================
# Main
# =====================================================================

function main()
    println("=" ^ 70)
    println("Tab Configuration Parametric Study")
    println("3 groups, 7 independent cases, controlled variable method")
    println("=" ^ 70)

    # Load base parameters
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    # Compute angles from reference mesh
    angles = compute_tab_angles(param_dim; nθ=80, gsorder=2)
    @printf("Inner spiral angle range: [%.2f, %.2f] rad\n", angles.theta_min, angles.theta_max)
    @printf("θ_mid=%.2f, Δθ=%.2f, θ_start=%.2f, θ_end=%.2f\n",
            angles.θ_mid, angles.Δθ, angles.θ_start, angles.θ_end)

    # Define cases
    cases = define_cases(angles)

    # Run all cases
    all_results = Dict{Int, Dict{String, Float64}}()
    all_raw     = Dict{Int, Dict{String, Any}}()

    for (case_id, label, group, theta_pos) in cases
        @printf("\n--- Case %d: %s ---\n", case_id, label)
        @printf("    theta_pos = %s\n", string(theta_pos))
        t0 = time_ns()
        result, mesh_data = run_tab_case(deepcopy(param_dim), theta_pos)
        dt_wall = (time_ns() - t0) * 1e-9

        metrics = extract_metrics(result)
        all_results[case_id] = metrics
        all_raw[case_id] = result

        @printf("  %.2f s | D_max=%.4f D_mean=%.6f T_max=%.1fK ΔT=%.2fK\n",
                dt_wall, metrics["D_max_end"], metrics["D_mean_end"],
                metrics["T_max_end"], metrics["ΔT_end"])
    end

    # Fill repeated cases
    all_results[5] = all_results[1]  # Case 5 = Case 1
    all_results[8] = all_results[2]  # Case 8 = Case 2
    all_raw[5] = all_raw[1]
    all_raw[8] = all_raw[2]

    # Print summary
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    @printf("%-6s %-25s %8s %8s %10s %8s %8s %8s\n",
            "Case", "Label", "D_max", "D_mean", "T_max[K]", "ΔT[K]", "1-D_mean", "V_end")
    println("-" ^ 95)

    all_labels = Dict(
        1 => "1-tab mid", 2 => "2-tab equal", 3 => "3-tab equal",
        4 => "1-tab start", 5 => "1-tab mid", 6 => "1-tab end",
        7 => "2-tab tight", 8 => "2-tab equal", 9 => "2-tab wide")

    for id in [1,2,3,4,5,6,7,8,9]
        m = all_results[id]
        @printf("%-6d %-25s %8.4f %8.6f %10.1f %8.2f %8.4f %8.3f\n",
                id, all_labels[id], m["D_max_end"], m["D_mean_end"],
                m["T_max_end"], m["ΔT_end"], m["soh_proxy"], m["V_end"])
    end

    # ── Save results ──
    out_dir = joinpath(root_dir, "output", "tab_parametric")
    mkpath(out_dir)

    # Save scalar metrics as CSV for Python post-processing
    # (Serialization format is Julia-specific; CSV is cross-language)
    metric_keys = ["D_max_end", "D_mean_end", "n_frac_end", "soh_proxy",
                   "T_max_end", "T_min_end", "T_mean_end", "ΔT_end", "V_end"]
    open(joinpath(out_dir, "tab_metrics.csv"), "w") do f
        println(f, "case_id,label,", join(metric_keys, ","))
        for id in [1,2,3,4,5,6,7,8,9]
            m = all_results[id]
            vals = [get(m, k, "NaN") for k in metric_keys]
            println(f, id, ",", all_labels[id], ",", join(vals, ","))
        end
    end

    # ── Julia summary plots ──
    groups = [
        ("Count Effect (Group 1)", [1, 2, 3]),
        ("Position Effect (Group 2)", [4, 5, 6]),
        ("Spacing Effect (Group 3)", [7, 8, 9]),
    ]
    metrics_keys = ["D_max_end", "ΔT_end", "soh_proxy"]
    titles = ["D_max", "ΔT [K]", "1-D_mean"]

    for (gi, (gname, case_ids)) in enumerate(groups)
        p = plot(title=gname, layout=(1, 3), size=(900, 300))
        for (mi, (key, ttl)) in enumerate(zip(metrics_keys, titles))
            vals = [all_results[c][key] for c in case_ids]
            bar!(p[mi], string.(case_ids), vals,
                 xlabel="Case", ylabel=ttl, label=false,
                 color=:steelblue, alpha=0.8)
        end
        savefig(p, joinpath(out_dir, "group$(gi)_summary.png"))
    end

    @printf("\nResults saved to %s\n", out_dir)
    println("=" ^ 70)
end

main()
```

- [ ] **Step 2: Verify script structure**

Read back the file and verify:
- `compute_tab_angles` uses normalized `param` (not `param_dim`)
- `run_tab_case` creates a `deepcopy` of `param_dim` before modifying `theta_pos`
- All 7 independent cases are defined in `define_cases`
- Repeated cases (5=1, 8=2) are filled after the loop
- Output directory is `output/tab_parametric/`

- [ ] **Step 3: Commit**

```bash
mkdir -p example/极耳参数研究
git add example/极耳参数研究/tab_parametric_study.jl
git commit -m "feat: add tab parametric study Julia simulation script"
```

---

## Chunk 3: Python Visualization Scripts

### Task 3: Create `param/fig_tab1_count_effect.py`

**Files:**
- Create: `param/fig_tab1_count_effect.py`
- Reference: `param/fig10_damage_cycle_comparison.py` (plot style), `param/fig8_single_tab_debonding.py` (spiral polar plot)

- [ ] **Step 1: Write the plot script**

```python
#!/usr/bin/env python3
"""Figure tab1: Tab count effect (Group 1, Cases 1-3).

Compares 1-tab, 2-tab, and 3-tab configurations at equal spacing.
Panels: (a) D_max(t), (b) T_max(t), (c) ΔT(t).
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    # ── Placeholder data (replace with simulation results) ──
    # TODO: Load from output/tab_parametric/tab_results.bin after Julia run
    t = np.linspace(0, 3600, 200)

    # Representative damage curves for different tab counts
    cases = {
        "1-tab (mid)": {"color": "#d62728", "marker": "o", "D_rate": 0.003, "T_max": 310.2, "dT": 4.8},
        "2-tab (equal)": {"color": "#2ca02c", "marker": "s", "D_rate": 0.002, "T_max": 307.5, "dT": 3.2},
        "3-tab (equal)": {"color": "#1f77b4", "marker": "D", "D_rate": 0.001, "T_max": 305.8, "dT": 2.1},
    }

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.6), dpi=300)

    # Panel (a): D_max(t)
    ax = axes[0]
    for label, cfg in cases.items():
        D = cfg["D_rate"] * t / 3600
        ax.plot(t, D, color=cfg["color"], lw=1.4, alpha=0.9)
        milestone_t = np.arange(0, 3601, 600)
        ax.plot(milestone_t, cfg["D_rate"] * milestone_t / 3600,
                marker=cfg["marker"], markersize=4,
                markerfacecolor="none", markeredgecolor=cfg["color"],
                markeredgewidth=1.2, linestyle="none")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$D_{\max}$")
    ax.set_title("(a) Damage")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Panel (b): T_max(t)
    ax = axes[1]
    for label, cfg in cases.items():
        T = 298.15 + (cfg["T_max"] - 298.15) * (1 - np.exp(-t / 1200))
        ax.plot(t, T, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$T_{\max}$ [K]")
    ax.set_title("(b) Temperature")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Panel (c): ΔT(t)
    ax = axes[2]
    for label, cfg in cases.items():
        dT = cfg["dT"] * (1 - np.exp(-t / 1000))
        ax.plot(t, dT, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\Delta T$ [K]")
    ax.set_title("(c) Temp. gradient")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Shared legend
    handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, marker=cfg["marker"],
               markersize=4, markerfacecolor="none", markeredgecolor=cfg["color"],
               markeredgewidth=1.2, label=label)
        for label, cfg in cases.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=7.5, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    out_path = os.path.join(output_dir, "fig_tab1_count_effect.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
```

- [ ] **Step 2: Verify script runs**

```bash
cd param && python fig_tab1_count_effect.py
```

Expected: prints `output/fig_tab1_count_effect.png` path without errors.

- [ ] **Step 3: Commit**

```bash
git add param/fig_tab1_count_effect.py
git commit -m "feat: add fig_tab1 tab count effect plot script"
```

### Task 4: Create `param/fig_tab2_position_effect.py`

**Files:**
- Create: `param/fig_tab2_position_effect.py`

- [ ] **Step 1: Write the plot script**

Follow the same structure as `fig_tab1_count_effect.py` but with:
- 3 cases: "1-tab start" (#d62728), "1-tab mid" (#2ca02c), "1-tab end" (#1f77b4)
- Same 3-panel layout: (a) D_max(t), (b) T_max(t), (c) ΔT(t)
- Title: "Tab Position Effect (Group 2, Cases 4-6)"
- Output: `output/fig_tab2_position_effect.png`

```python
#!/usr/bin/env python3
"""Figure tab2: Tab position effect (Group 2, Cases 4-6).

Compares single tab at start, mid, and end of positive spiral region.
Panels: (a) D_max(t), (b) T_max(t), (c) ΔT(t).
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    t = np.linspace(0, 3600, 200)

    cases = {
        "1-tab start":  {"color": "#d62728", "marker": "o", "D_rate": 0.0035, "T_max": 311.0, "dT": 5.2},
        "1-tab mid":    {"color": "#2ca02c", "marker": "s", "D_rate": 0.0030, "T_max": 310.2, "dT": 4.8},
        "1-tab end":    {"color": "#1f77b4", "marker": "D", "D_rate": 0.0025, "T_max": 309.5, "dT": 4.3},
    }

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.6), dpi=300)

    # Panel (a): D_max(t)
    ax = axes[0]
    for label, cfg in cases.items():
        D = cfg["D_rate"] * t / 3600
        ax.plot(t, D, color=cfg["color"], lw=1.4, alpha=0.9)
        milestone_t = np.arange(0, 3601, 600)
        ax.plot(milestone_t, cfg["D_rate"] * milestone_t / 3600,
                marker=cfg["marker"], markersize=4,
                markerfacecolor="none", markeredgecolor=cfg["color"],
                markeredgewidth=1.2, linestyle="none")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$D_{\max}$")
    ax.set_title("(a) Damage")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Panel (b): T_max(t)
    ax = axes[1]
    for label, cfg in cases.items():
        T = 298.15 + (cfg["T_max"] - 298.15) * (1 - np.exp(-t / 1200))
        ax.plot(t, T, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$T_{\max}$ [K]")
    ax.set_title("(b) Temperature")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    # Panel (c): ΔT(t)
    ax = axes[2]
    for label, cfg in cases.items():
        dT = cfg["dT"] * (1 - np.exp(-t / 1000))
        ax.plot(t, dT, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\Delta T$ [K]")
    ax.set_title("(c) Temp. gradient")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, marker=cfg["marker"],
               markersize=4, markerfacecolor="none", markeredgecolor=cfg["color"],
               markeredgewidth=1.2, label=label)
        for label, cfg in cases.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=7.5, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    out_path = os.path.join(output_dir, "fig_tab2_position_effect.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
```

- [ ] **Step 2: Verify script runs**

```bash
cd param && python fig_tab2_position_effect.py
```

Expected: prints output path without errors.

- [ ] **Step 3: Commit**

```bash
git add param/fig_tab2_position_effect.py
git commit -m "feat: add fig_tab2 tab position effect plot script"
```

### Task 5: Create `param/fig_tab3_spacing_effect.py`

**Files:**
- Create: `param/fig_tab3_spacing_effect.py`

- [ ] **Step 1: Write the plot script**

Same 3-panel structure as fig_tab1/tab2, but with:
- 3 cases: "2-tab tight" (#d62728), "2-tab equal" (#2ca02c), "2-tab wide" (#1f77b4)
- Title: "Tab Spacing Effect (Group 3, Cases 7-9)"
- Output: `output/fig_tab3_spacing_effect.png`

```python
#!/usr/bin/env python3
"""Figure tab3: Tab spacing effect (Group 3, Cases 7-9).

Compares 2-tab configurations with tight, moderate, and wide spacing.
Panels: (a) D_max(t), (b) T_max(t), (c) ΔT(t).
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    t = np.linspace(0, 3600, 200)

    cases = {
        "2-tab tight":  {"color": "#d62728", "marker": "o", "D_rate": 0.0025, "T_max": 308.8, "dT": 3.8},
        "2-tab equal":  {"color": "#2ca02c", "marker": "s", "D_rate": 0.0020, "T_max": 307.5, "dT": 3.2},
        "2-tab wide":   {"color": "#1f77b4", "marker": "D", "D_rate": 0.0018, "T_max": 306.8, "dT": 2.5},
    }

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.6), dpi=300)

    ax = axes[0]
    for label, cfg in cases.items():
        D = cfg["D_rate"] * t / 3600
        ax.plot(t, D, color=cfg["color"], lw=1.4, alpha=0.9)
        milestone_t = np.arange(0, 3601, 600)
        ax.plot(milestone_t, cfg["D_rate"] * milestone_t / 3600,
                marker=cfg["marker"], markersize=4,
                markerfacecolor="none", markeredgecolor=cfg["color"],
                markeredgewidth=1.2, linestyle="none")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$D_{\max}$")
    ax.set_title("(a) Damage")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    ax = axes[1]
    for label, cfg in cases.items():
        T = 298.15 + (cfg["T_max"] - 298.15) * (1 - np.exp(-t / 1200))
        ax.plot(t, T, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$T_{\max}$ [K]")
    ax.set_title("(b) Temperature")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    ax = axes[2]
    for label, cfg in cases.items():
        dT = cfg["dT"] * (1 - np.exp(-t / 1000))
        ax.plot(t, dT, color=cfg["color"], lw=1.4, alpha=0.9)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\Delta T$ [K]")
    ax.set_title("(c) Temp. gradient")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)

    handles = [
        Line2D([0], [0], color=cfg["color"], lw=1.4, marker=cfg["marker"],
               markersize=4, markerfacecolor="none", markeredgecolor=cfg["color"],
               markeredgewidth=1.2, label=label)
        for label, cfg in cases.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=7.5, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    out_path = os.path.join(output_dir, "fig_tab3_spacing_effect.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
```

- [ ] **Step 2: Verify script runs**

```bash
cd param && python fig_tab3_spacing_effect.py
```

Expected: prints output path without errors.

- [ ] **Step 3: Commit**

```bash
git add param/fig_tab3_spacing_effect.py
git commit -m "feat: add fig_tab3 tab spacing effect plot script"
```

### Task 6: Create `param/fig_tab4_summary_bar.py`

**Files:**
- Create: `param/fig_tab4_summary_bar.py`
- Reference: `param/fig12_charge_rate_soh_comparison.py` (bar chart style)

- [ ] **Step 1: Write the summary bar chart script**

```python
#!/usr/bin/env python3
"""Figure tab4: Cross-group summary bar chart.

Grouped bar chart comparing D_max, T_max, ΔT, SOH across all 9 cases
(3 groups × 3 cases each).
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt


def make_figure() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(root, "output")
    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 9

    # ── Placeholder data (replace with simulation results) ──
    group_labels = ["Count\n(1/2/3-tab)", "Position\n(start/mid/end)", "Spacing\n(tight/equal/wide)"]
    # Each group has 3 values: [case_a, case_b, case_c]
    data = {
        r"$D_{\max}$": {
            "vals": [[0.003, 0.002, 0.001], [0.0035, 0.003, 0.0025], [0.0025, 0.002, 0.0018]],
            "colors": [["#d62728", "#2ca02c", "#1f77b4"]] * 3,
        },
        r"$\Delta T$ [K]": {
            "vals": [[4.8, 3.2, 2.1], [5.2, 4.8, 4.3], [3.8, 3.2, 2.5]],
            "colors": [["#d62728", "#2ca02c", "#1f77b4"]] * 3,
        },
        "SOH": {
            "vals": [[0.997, 0.998, 0.999], [0.996, 0.997, 0.998], [0.998, 0.998, 0.999]],
            "colors": [["#d62728", "#2ca02c", "#1f77b4"]] * 3,
        },
    }

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.8), dpi=300)

    bar_labels = ["a", "b", "c"]

    for ax, (metric, d) in zip(axes, data.items()):
        x = np.arange(len(group_labels))
        width = 0.25

        for j in range(3):
            vals = [d["vals"][g][j] for g in range(3)]
            color = d["colors"][0][j]
            offset = (j - 1) * width
            bars = ax.bar(x + offset, vals, width, label=bar_labels[j],
                          color=color, alpha=0.85, edgecolor="white", linewidth=0.5)

        ax.set_xlabel("Group")
        ax.set_ylabel(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(group_labels, fontsize=7)
        ax.grid(True, axis="y", ls="--", lw=0.4, alpha=0.5)

    # Legend with case descriptions
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d62728", alpha=0.85, label="Case a (1st)"),
        Patch(facecolor="#2ca02c", alpha=0.85, label="Case b (2nd)"),
        Patch(facecolor="#1f77b4", alpha=0.85, label="Case c (3rd)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=7.5, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, -0.04))

    fig.tight_layout(rect=[0, 0.06, 1, 1])
    out_path = os.path.join(output_dir, "fig_tab4_summary_bar.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(make_figure())
```

- [ ] **Step 2: Verify script runs**

```bash
cd param && python fig_tab4_summary_bar.py
```

Expected: prints output path without errors.

- [ ] **Step 3: Commit**

```bash
git add param/fig_tab4_summary_bar.py
git commit -m "feat: add fig_tab4 cross-group summary bar chart"
```

---

## Chunk 4: Integration

### Task 7: Update `param/generate_all_figures.py`

**Files:**
- Modify: `param/generate_all_figures.py`

- [ ] **Step 1: Add fig_tab1–4 imports and entries**

In `param/generate_all_figures.py`, make two precise edits:

**Edit A:** After the line `from fig12_charge_rate_soh_comparison import make_figure as make_fig12` (line 17), add 4 new import lines:

```python
from fig_tab1_count_effect import make_figure as make_fig_tab1
from fig_tab2_position_effect import make_figure as make_fig_tab2
from fig_tab3_spacing_effect import make_figure as make_fig_tab3
from fig_tab4_summary_bar import make_figure as make_fig_tab4
```

**Edit B:** Inside the `main()` function, after the last existing entry in the `makers` list (`("fig12_charge_rate_soh_comparison", make_fig12),` on line 33), **append** (do not replace) 4 new entries:

```python
        ("fig_tab1_count_effect", make_fig_tab1),
        ("fig_tab2_position_effect", make_fig_tab2),
        ("fig_tab3_spacing_effect", make_fig_tab3),
        ("fig_tab4_summary_bar", make_fig_tab4),
```

- [ ] **Step 2: Run generate_all_figures.py to verify all plots render**

```bash
cd param && python generate_all_figures.py
```

Expected: All 16 figures (fig1–fig12, fig_tab1–4) print their output paths without errors.

- [ ] **Step 3: Commit**

```bash
git add param/generate_all_figures.py
git commit -m "feat: register fig_tab1-4 in generate_all_figures.py"
```

### Task 8: Update CLAUDE.md documentation index

**Files:**
- Modify: `CLAUDE.md` (if needed — add md/15 reference to the documentation table)

- [ ] **Step 1: Add md/15 to the documentation index in CLAUDE.md**

In the "第二层：模型实现 (04-07)" section or a new "第五层：参数研究 (15)" section, add:

```markdown
### 第五层：参数研究 (15)
| 编号 | 文档 | 内容 |
|------|------|------|
| 15 | `15_参数研究_极耳工况.md` | 正交分组法极耳工况设计、控制变量法、输出指标 |
```

Also add to the example files table:

```markdown
| `example/极耳参数研究/tab_parametric_study.jl` | 极耳工况参数研究 (7 cases) |
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add md/15 tab parametric study to CLAUDE.md index"
```
