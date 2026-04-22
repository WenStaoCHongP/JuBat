# 网格敏感性分析实施计划

> **Spec Document:** [2026-04-22-grid-sensitivity-analysis-design.md](../specs/2026-04-22-grid-sensitivity-analysis-design.md)
>
> **For agentic workers:** REQUIRED: Keep this plan synchronized with findings/progress. Do not cross-derive thermal and cohesive criteria unless the spec explicitly says so.
>
> **Context:** Based on the discussion recorded in [findings.md](../../planning-with-files/网格敏感性分析/findings.md).
>
> **Updated:** 2026-04-22 (rev 2) — thermal criterion replaced by Biot number + systematic refinement; cohesive $l_c$ confirmed with $E_{eff}$; pure thermal/mechanical paths clarified.

**Goal:** For the Jellyroll coupled case, establish three independent mesh sensitivity tracks:

1. Electrochemical mesh
2. Thermal mesh
3. Cohesive mesh

Each track must:

- compare candidates against its own finest mesh
- keep the governing metric deviation within 5%
- use a fixed time-step policy within the track
- record conclusions in a form that can be reused by later implementation work

**Locked Decisions:**

- Thermal mesh criterion: Biot number analysis ($Bi_t \approx 0.004 \ll 0.1$) shows near-axisymmetric temperature field; use systematic refinement $n\theta \in \{20, 40, 80, 160\}$
- Cohesive length: $l_c = G_c \cdot E_{eff} / \sigma_{max}^2 \approx 96\ \mu\text{m}$, where $E_{eff} \approx 25.6\ \text{GPa}$ from `Mechanical.jl:181`
- Thermal benchmark: uniform volumetric heat source, pure thermal path (`opt.model = "thermal"`), following `example/热模块验证/thermal_verify.jl`
- Lumped thermal gradient metric: max $|dT/dt|$, not spatial gradient
- CZM load-displacement curve: pure mechanical model with displacement boundary conditions triggering damage
- nθ definition: per-revolution circumferential resolution

**Working assumptions:**

- `gsorder` is fixed at 2 unless a later review changes it
- `Nrn/Nrp` remain fixed at their current default values in the first pass; if electrochemical error stays above threshold, particle refinement becomes a second-stage task

---

## 1. Electrochemical Track

### Scope

Candidate sets stay fixed as the four existing `opt.Nn / opt.Ns / opt.Np` combinations:

- `(40, 20, 40)`
- `(20, 10, 20)`
- `(20, 5, 20)`
- `(10, 5, 10)`

### Metrics

- Cell voltage curve
- Temperature peak
- Maximum lumped thermal gradient $\max |dT/dt|$

### Execution notes

- Use the same thermal benchmark and the same time-step policy for all four cases
- If the error budget is not met after fixing the electrochemical mesh, revisit the particle mesh as a separate extension rather than mixing criteria in the first pass

---

## 2. Thermal Track

### Scope

Thermal mesh is determined by Biot number analysis + systematic refinement.

**Physical basis:**

The Jellyroll cell has strongly anisotropic thermal conductivity:

$$\lambda_r \approx 1.32\ \text{W/m/K},\quad \lambda_t \approx 25.3\ \text{W/m/K}\quad (\lambda_t/\lambda_r \approx 19\times)$$

The tangential Biot number is:

$$Bi_t = \frac{h \cdot R_{out}}{\lambda_t} = \frac{10 \times 0.01015}{25.3} \approx 0.004 \ll 0.1$$

This means the temperature field is nearly axisymmetric — angular resolution has minimal effect on thermal accuracy. The thermal penetration depth $\delta_T = \sqrt{\alpha_t \cdot 3600} \approx 194\ \text{mm} \gg R_{out} = 10.15\ \text{mm}$ also confirms full thermal penetration.

**Conclusion:** Thermal nθ sensitivity is expected to be very weak. Use systematic refinement to confirm.

**Candidate values (per revolution):**

| nθ | Outer arc length | Inner arc length | Total elements ≈ |
|----|-----------------|-----------------|------------------|
| 20 | 3.19 mm | 0.60 mm | 440 |
| 40 | 1.59 mm | 0.30 mm | 880 |
| 80 | 0.80 mm | 0.15 mm | 1760 |
| 160 | 0.40 mm | 0.08 mm | 3520 |

### Model basis

- Pure thermal model: `opt.model = "thermal"`, `opt.thermalmodel = "distributed2D"`
- Uniform volumetric heat source: $q_0 = 2 \times 10^5\ \text{W/m}^3$ (following `example/热模块验证/thermal_verify.jl`)
- Surface cooling boundary condition
- No SPMe coupling inside the thermal sensitivity sweep

### Metrics

- Temperature peak $T_{max}$
- Temperature range $T_{max} - T_{min}$ (spatial uniformity indicator)
- Angular variation convergence (reference: `thermal_verify.jl` angular_variation function)

### Execution notes

- Keep the uniform heat source $q_0$ fixed across all four nθ levels
- Monitor both absolute T_max convergence and angular symmetry convergence
- Expected result: rapid convergence due to near-axisymmetry; likely nθ ≥ 40 sufficient

---

## 3. Cohesive Track

### Scope

Cohesive mesh is determined by the cohesive characteristic length.

**Physical basis:**

$$l_c = \frac{G_c \cdot E_{eff}}{\sigma_{max}^2}$$

All parameters confirmed from code:

| Parameter | Value | Source |
|-----------|-------|--------|
| $G_c$ | 25.3 J/m² | `Jellyroll.jl` Cohesive 法向断裂能 |
| $\sigma_{max}$ | 82 MPa | `Jellyroll.jl` Cohesive 法向峰值强度 |
| $E_{eff}$ | 25.6 GPa | `Mechanical.jl:181` 厚度加权平均 |

$$l_c = \frac{25.3 \times 25.6 \times 10^9}{(82 \times 10^6)^2} \approx 96\ \mu\text{m}$$

**nθ interval (per revolution, arc length < $l_c$):**

- Inner constraint: $n\theta > 2\pi R_{in}/l_c \approx 126$
- Outer constraint: $n\theta > 2\pi R_{out}/l_c \approx 664$

Split this interval into four levels (round to convenient integers).

### Model basis

- Pure mechanical model with displacement boundary conditions
- Load-displacement curve obtained by incrementally applying prescribed displacement to trigger CZM damage evolution
- No thermal or electrochemical coupling in this track

### Metrics

- Stress peak (von Mises)
- Load-displacement curve (global response under prescribed displacement)
- Interface traction-separation curve (local CZM constitutive response)
- Damage onset time ($D_{max}$ reaching threshold)
- $D_{max}$ evolution
- `n_fractured(t)` growth rate

### Execution notes

- Keep the local traction-separation response and the global load-displacement response as separate outputs
- Use the same geometric mesh family only if the final comparison step explicitly needs a shared geometry; do not force a shared rule before the independent sweeps are done

---

## 4. Energy Conservation Check (Full Coupling)

### Purpose

Verify that the fully coupled simulation (electrochemical + thermal + elastic + CZM) satisfies the first law of thermodynamics. This is a global verification step, not a mesh sensitivity track.

### Energy balance equation

$$R(t) = W_{elec}(t) - Q_{loss}(t) - \Delta E_{thermal}(t) - \Delta E_{elastic}(t) - E_{fracture}(t) - \Delta E_{chem}(t)$$

$$\epsilon_R(t) = \frac{|R(t)|}{|W_{elec}(t)|}$$

### Energy terms

| Term | Formula | Data source | Status |
|------|---------|-------------|--------|
| $W_{elec}$ | $\int_0^t V \cdot I \, d\tau$ | `cell voltage [V]`, `cell current [A]` | ✅ Direct |
| $\Delta E_{th}$ | $\sum_n M_n [T_n(t) - T_n(0)]$ | Node temperatures + thermal mass matrix | ⚠️ Needs implementation |
| $\Delta E_{el}$ | $\frac{1}{2} u^T K_{bulk} u \big\vert_0^t$ | CZM displacement + bulk stiffness (`czm.jl:329`) | ⚠️ Needs implementation |
| $E_{frac}$ | $\sum_e G_c \cdot l_e \cdot D_e(t)$ | Cohesive damage states + element lengths | ⚠️ Needs implementation |
| $Q_{loss}$ | $\int_0^t \int_\Gamma h(T - T_{amb}) \, d\Gamma \, d\tau$ | Boundary node temperatures + areas | ⚠️ Needs implementation |
| $\Delta E_{chem}$ | $\int_V [U_{OCV}(SOC(t)) - U_{OCV}(SOC_0)] \cdot F \cdot c_{s,max} \, dV$ | `element soc_n/soc_p` + OCV functions | ⚠️ Needs implementation |

### Simplified fallback

If full chemical energy proves difficult to compute, use the thermal-elastic-fracture subsystem:

$$R_{sub}(t) = Q_{generated}(t) - \Delta E_{th}(t) - Q_{loss}(t) - \Delta E_{el}(t) - E_{frac}(t)$$

where $Q_{generated}$ is obtained by integrating `total heat source [W]` over time (already available).

### Acceptance criteria

- $\epsilon_R(t) < 1\%$ throughout the simulation
- Residual $R(t)$ should not exhibit systematic growth (indicates energy leak)
- Finer mesh should yield smaller residual

---

## 5. Planned Execution Order

- [x] Confirm the final formulas and assumptions (thermal: Biot + systematic; cohesive: $l_c$ with $E_{eff}$)
- [x] Confirm thermal nθ candidates: $\{20, 40, 80, 160\}$ (per revolution)
- [x] Define energy conservation check formulation
- [ ] Derive the four cohesive `nθ` levels from $l_c \approx 96\ \mu\text{m}$ (interval 126–664)
- [ ] Implement energy balance post-processing functions
- [ ] Run the electrochemical four-case sweep
- [ ] Run the thermal four-case sweep (pure thermal, uniform $q_0$)
- [ ] Run the cohesive four-case sweep (pure mechanical, displacement BCs)
- [ ] Run the fully coupled case at recommended mesh and compute energy balance
- [ ] Compare every candidate against the finest case in the same family
- [ ] Compute error tables, curves, and the final recommendation
- [ ] Sync the conclusions back into findings and progress notes

---

## 6. Deliverables

- Mesh configuration table
- Metric comparison table
- Error curves and convergence plots
- Final recommended mesh choices for each family
- Energy balance residual plot $R(t)$ and $\epsilon_R(t)$ for the fully coupled case
- Residual issue log for any metric that still needs a follow-up definition

---

## 7. Risks

- If a reference value is near zero, relative error becomes misleading; switch to absolute error or event-time difference
- The thermal nθ sensitivity is expected to be weak ($Bi_t \approx 0.004$); if convergence is too fast to produce meaningful differentiation, consider adding a non-uniform heat source as a secondary test
- If the cohesive nθ range (126–664) leads to excessive computational cost, the upper bound can be relaxed (e.g., accepting $l_c / 2$ instead of $l_c$ as the minimum element size)
- If the global load-displacement curve and the local traction-separation curve are mixed together, the CZM interpretation becomes ambiguous
- Keep `gsorder` and the time-step policy fixed during the sweep so spatial and temporal errors are not conflated
- Energy balance residual may be dominated by time integration error rather than spatial error; differentiate between the two by also varying `dt` if needed
- If the full chemical energy term $\Delta E_{chem}$ is too complex to implement, fall back to the simplified subsystem check (§4 Simplified fallback)

---

## 8. Reference Files

- [findings.md](../../planning-with-files/网格敏感性分析/findings.md) — upstream discussion and constraints
- [progress.md](../../planning-with-files/网格敏感性分析/progress.md) — session log and status
- [2026-04-22-grid-sensitivity-analysis-design.md](../specs/2026-04-22-grid-sensitivity-analysis-design.md) — formal definitions and assumptions
