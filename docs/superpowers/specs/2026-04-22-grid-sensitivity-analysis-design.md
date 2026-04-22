# 网格敏感性分析设计规格

> **Date:** 2026-04-22
>
> **Status:** Draft (rev 2 — thermal criterion updated)
>
> **Plan:** [2026-04-22-grid-sensitivity-analysis-plan.md](../plans/2026-04-22-grid-sensitivity-analysis-plan.md)
>
> **Source discussion:** [findings.md](../../planning-with-files/网格敏感性分析/findings.md)

---

## 1. Purpose

This specification defines the mesh-sensitivity policy for the Jellyroll coupled case in three separate families:

1. Electrochemical mesh
2. Thermal mesh
3. Cohesive mesh

The key requirement is that each family must be evaluated against its own finest candidate, with the same-family deviation kept within 5% unless a near-zero reference requires absolute or event-time error instead.

---

## 2. Confirmed Physical Criteria

### 2.1 Thermal criterion

The thermal family uses systematic refinement based on Biot number analysis, not derived from the cohesive process-zone length.

**Physical basis:** The Jellyroll cell has strongly anisotropic thermal conductivity ($\lambda_t / \lambda_r \approx 19\times$). The tangential Biot number $Bi_t = h \cdot R_{out} / \lambda_t \approx 0.004 \ll 0.1$, indicating the temperature field is nearly axisymmetric. The thermal penetration depth $\delta_T \approx 194\ \text{mm} \gg R_{out} = 10.15\ \text{mm}$, confirming full thermal penetration and no angular boundary layers.

**Consequence:** Thermal accuracy is insensitive to angular resolution (nθ). Use systematic refinement to confirm convergence.

**nθ definition:** per-revolution circumferential resolution.

**Candidate values:** $n\theta \in \{20, 40, 80, 160\}$.

### 2.2 Cohesive criterion

The cohesive family is determined by the cohesive characteristic length.

Use the confirmed estimate:

$$
l_c = \frac{G_c \cdot E_{eff}}{\sigma_{max}^2} \approx 96\ \mu\text{m}
$$

where:

- $G_c = 25.3\ \text{J/m}^2$ — fracture energy (`Jellyroll.jl`)
- $E_{eff} = 25.6\ \text{GPa}$ — thickness-weighted bulk Young's modulus (`Mechanical.jl:181`)
- $\sigma_{max} = 82\ \text{MPa}$ — peak strength (`Jellyroll.jl`)

This is a bulk-material length scale computed as $E_{eff} = (NE.E \cdot t_{NE} + PE.E \cdot t_{PE}) / (t_{NE} + t_{PE})$. It is not the cohesive penalty stiffness and must not be substituted by `K_n`.

**nθ interval (per revolution):** 126 (inner) to 664 (outer), split into four levels.

### 2.3 Heat-source basis

The thermal sweep is based on a **uniform volumetric heat source** conduction benchmark.

Implementation follows `example/热模块验证/thermal_verify.jl`:
- `opt.model = "thermal"`, `opt.thermalmodel = "distributed2D"`
- $q_0 = 2 \times 10^5\ \text{W/m}^3$, applied as `q_func = (r, theta, t) -> q0`
- Surface cooling boundary condition
- No SPMe coupling inside the thermal sensitivity sweep

### 2.4 Lumped-gradient rule

For the lumped thermal model, the gradient metric is the time derivative of temperature, not a spatial gradient.

Use:

$$
\max |dT/dt|
$$

as the lumped gradient comparison metric.

### 2.5 CZM response policy

CZM outputs keep both of the following response families:

- **load-displacement curve**: obtained from pure mechanical model with displacement boundary conditions, incrementally applying prescribed displacement to trigger CZM damage evolution
- **traction-separation curve**: local CZM constitutive response at interface elements

The two curves are distinct and should be compared separately.

---

## 3. Metric Definitions

### 3.1 Electrochemical metrics

- Cell voltage curve
- Temperature peak
- Maximum lumped thermal gradient $\max |dT/dt|$

### 3.2 Thermal metrics

- Temperature peak $T_{max}$
- Temperature range $T_{max} - T_{min}$ (spatial uniformity)
- Angular variation convergence

### 3.3 Cohesive metrics

- Stress peak
- Load-displacement curve deviation
- Traction-separation curve deviation
- Damage onset time
- `D_max`
- `n_fractured(t)` growth rate

### 3.4 Energy conservation (full coupling)

Energy balance residual:

$$R(t) = W_{elec}(t) - Q_{loss}(t) - \Delta E_{thermal}(t) - \Delta E_{elastic}(t) - E_{fracture}(t) - \Delta E_{chem}(t)$$

Relative error: $\epsilon_R = |R| / |W_{elec}| < 1\%$.

Simplified fallback (thermal-elastic-fracture subsystem):

$$R_{sub}(t) = Q_{generated}(t) - \Delta E_{th}(t) - Q_{loss}(t) - \Delta E_{el}(t) - E_{frac}(t)$$

If the metric is an event time or a near-zero amplitude, use absolute error, time difference, or curve-area difference rather than raw relative error.

---

## 4. Working Assumptions

- `gsorder` remains fixed at 2 during the first-pass sweep
- The time-step strategy remains fixed within each family so that temporal error does not contaminate mesh-error interpretation
- `Nrn/Nrp` remain fixed in the first pass; if the electrochemical family does not meet the target, particle refinement can be added as a second-stage scope
- The heat-source history remains fixed for the thermal family
- The final thermal and cohesive grid families may be kept independent even if they are later aligned geometrically for a joint comparison

---

## 5. Open Questions

- ~~What exact thermal candidate values should be chosen?~~ **Resolved:** $n\theta \in \{20, 40, 80, 160\}$
- Should the thermal family and cohesive family share a final geometry after independent sweeps, or remain fully separate?
- If the electrochemical family fails to meet the target with fixed particle mesh, should the particle sweep be promoted into the main scope or kept as a follow-up?
- What are the exact four cohesive nθ levels after splitting the 126–664 interval?

---

## 6. Acceptance Criteria

A candidate family is considered acceptable when:

- Every candidate is compared against the finest candidate in the same family
- The target metric deviation is less than 5%, or the alternative error definition is used for near-zero references
- The definitions, assumptions, and unresolved questions are recorded in the plan and progress notes
- The chosen mesh family can be reused without rewriting the physical criterion from scratch

---

## 7. Reference Material

- [findings.md](../../planning-with-files/网格敏感性分析/findings.md) — review conclusions and discussion record
- [2026-04-22-grid-sensitivity-analysis-plan.md](../plans/2026-04-22-grid-sensitivity-analysis-plan.md) — execution plan
