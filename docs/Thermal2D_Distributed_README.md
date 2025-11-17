# 2D Distributed Thermal Model for Jellyroll Batteries

This document describes the implementation of a 2D distributed thermal model for simulating temperature field evolution in jellyroll batteries with electrochemical-thermal coupling.

## Overview

The model simulates the electrochemical behavior and thermal response of jellyroll (wound cylindrical) lithium-ion batteries. Key features include:

- **2D Distributed Thermal Model**: Captures spatial temperature variations across the battery cross-section
- **Electrochemical-Thermal Coupling**: Strong coupling between electrochemical reactions and heat generation
- **Per-Element SPMe**: Each thermal element runs its own SPMe model for accurate local behavior
- **Parallel Current Distribution**: Non-linear current splitting based on voltage consistency
- **Collector-Seeded Mesh**: Optimized mesh generation that respects layer boundaries

## Mathematical Model

### Electrochemical Model (SPMe)
- Single Particle Model with electrolyte (SPMe)
- Includes lithium diffusion in particles, electrolyte transport, and electrochemical kinetics
- Butler-Volmer kinetics with overpotentials

### Thermal Model
- 2D heat conduction equation: `(ρc)∂T/∂t = ∇·(K∇T) + q`
- Anisotropic thermal conductivity for jellyroll geometry
- Heat sources: reaction heat, reversible heat, ohmic heat
- Boundary conditions: natural convection on outer surface, insulation on inner surface

### Coupling Strategy
- **Strong Coupling**: Iterative solution within each time step
  1. Estimate temperature field
  2. Solve parallel current distribution
  3. Run per-element SPMe for heat sources
  4. Solve 2D thermal equation
  5. Check convergence and iterate

## Usage

### Basic Example

```julia
using JuBat

# Set up parameters
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()

# Enable thermal modeling
opt.thermal_enabled = true
opt.collector_seeded = true
opt.per_element_spme = true
opt.coupling_mode = "strong"
opt.thermalmodel = "distributed2D"
opt.units_thermal = "SI"

# Time and current settings
opt.dt = [1.0, 5.0]
opt.time = [0.0, 300.0]  # 5 minutes
opt.Current = t -> 2.0 * param_dim.cell.I1C  # 2C discharge

# Create case and mesh
case = JuBat.SetCase(param_dim, opt)
mesh_th = JuBat.jellyroll_Q4_mesh(param_dim; nx=180, gsorder=2, crop_mode=:collector_seeded)
case.mesh["thermal2D"] = mesh_th

# Run simulation
result = JuBat.Solve(case)
```

### Key Options

- `opt.coupling_mode`: `"strong"` or `"weak"` coupling
- `opt.per_element_spme`: Enable per-element electrochemical calculations
- `opt.collector_seeded`: Use optimized mesh generation
- `opt.parallel_solve_V`: Enable voltage-based current distribution
- `opt.units_thermal`: `"SI"` or `"nd"` (dimensionless) units

## Output Analysis

### Temperature Field
- `result["thermal2D T_nodes [K]"]`: Nodal temperatures
- `result["thermal2D nodes xy [m]"]`: Node coordinates

### Current Distribution
- `variables["thermal2D element current"]`: Current through each element
- `variables["thermal2D common voltage"]`: Common terminal voltage

### Heat Sources
- `variables["heat_source_fields"]`: Volumetric heat generation per element

## Validation

The implementation includes several validation checks:

- **Energy Conservation**: Heat generation vs. heat dissipation vs. energy storage
- **Current Conservation**: Sum of element currents equals total current
- **Temperature Bounds**: Physically reasonable temperature ranges
- **Convergence**: Iterative coupling convergence monitoring

## Mesh Generation

The collector-seeded mesh generation creates elements that span complete layer sequences:

```julia
mesh = jellyroll_Q4_mesh(param_dim; nx=180, crop_mode=:collector_seeded)
```

This ensures each thermal element contains all material layers (PCC, PE, SP, NE, NCC) with proper weighting.

## Performance Considerations

- **Strong Coupling**: More accurate but computationally expensive
- **Mesh Resolution**: Higher `nx` improves accuracy but increases computation time
- **Time Stepping**: Adaptive time stepping based on convergence

## Troubleshooting

### Common Issues

1. **Convergence Failure**: Reduce time step or switch to weak coupling
2. **Memory Issues**: Reduce mesh resolution or disable debugging
3. **Unphysical Results**: Check parameter values and boundary conditions

### Debugging Options

```julia
opt.debug_coupling = true
opt.debug_sample_elems = 5  # Monitor first 5 elements
```

## References

- Full theoretical derivation in `docs/SPMe_Thermal2D_Theory_vs_Code.md`
- Example implementation in `example/jellyroll_thermal_2d_example.jl`
- Validation cases in `example/jellyroll_coupled_example.jl`