# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

JuBat is a Julia-based framework for battery modeling based on Newman's battery models, including P2D, SPM, and SPMe models. The framework uses 2nd order finite element method (FEM) for solving equations and supports electro-thermal-mechanical coupling with cohesive zone modeling (CZM) for degradation analysis.

Reference: W. Ai, Y. Liu, Improving the convergence rate of Newman's battery model using 2nd order finite element method, J. Energy Storage. 67 (2023) 107512.

## Development Workflow

### Running Examples
```bash
# Basic usage - include JuBat module and run examples
julia example/minimal_example.jl
julia example/cycle_example.jl
julia example/czm_cycle_example.jl
julia example/SPMe_Thermal_example.jl
julia example/testexample.jl  # Full coupled simulation
```

### Module Structure
The main module is defined in `src/JuBat.jl` with all submodules included. To use in scripts:
```julia
include("src/JuBat.jl")
using .JuBat
```

## Core Architecture

### Model Types
- **P2D** (Pseudo-Two-Dimensional): Full electrochemical model with particle diffusion and electrolyte transport
- **SPM** (Single Particle Model): Simplified model assuming electrolyte concentration uniformity
- **SPMe** (Single Particle Model extended): Includes electrolyte dynamics

### Coupling Capabilities
- **Electro-thermal coupling**: Per-element SPMe with distributed 2D thermal model
- **Electro-mechanical coupling**: Thermal and diffusion stress calculations
- **CZM (Cohesive Zone Model)**: Interlayer degradation and fracture modeling
- **Cycle solver**: Multi-cycle charge/discharge simulation with state management

### Key Data Structures

#### Case Configuration
```julia
param_dim = JuBat.ChooseCell("LG M50" or "Jellyroll")
opt = JuBat.Option()  # Configure model, time stepping, thermal/mechanical options
case = JuBat.SetCase(param_dim, opt)
```

#### Mesh Creation
```julia
# Standard 1D mesh (automatically created by SetCase)
# Jellyroll collector-seeded mesh for thermal2D
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
```

#### Main Solver
```julia
result = JuBat.Solve(case)
# Returns Dict with time, voltage, temperature, concentrations, stresses, etc.
```

## Configuration Options

### Option Structure (`src/Option.jl`)
Key fields:
- `model`: "P2D", "SPM", or "SPMe"
- `Np`, `Ns`, `Nn`: Grid points for positive electrode, separator, negative electrode
- `Nrp`, `Nrn`: Radial grid points for positive/negative particles
- `gsorder`: Order of Gaussian quadrature (affects accuracy)
- `Current`: Function defining current profile I(t)
- `time`: [t_start, t_end] simulation time range
- `dt`: [dt_min, dt_max] adaptive time stepping range
- `solveType`: "Crank-Nicolson", "forward", or "backward"

### Thermal Options
- `thermal_enabled`: Enable thermal coupling
- `thermalmodel`: "none", "lumped", "distributed1D", "distributed2D"
- `thermal_dim`: "1D" or "2D"
- `cool_method`: "tab" or "surface" cooling
- `per_element_spme`: Enable per-element SPMe for thermal coupling

### Mechanical/CZM Options
- `mechanicalmodel`: "none" or "full"
- `czm_enabled`: Enable cohesive zone modeling
- `czm_update_interval`: How often to update damage (1=every step)
- `czm_iter_method`: "basic", "load_substep", or "arc_length"

## Advanced Features

### Multi-SPMe Parallel Architecture
When `opt.per_element_spme=true`, each thermal element has its own SPMe model:
- Current distribution solved via Newton-Raphson (`solve_branch_currents_newton`)
- Per-element heat sources computed from overpotentials
- Thermal variables passed back to electrochemical solver
- See `src/Parallelsolution.jl` for implementation

### Jellyroll Geometry
Specialized mesh for spiral-wound batteries:
- `jellyroll_collector_seed_mesh`: Creates Q4 mesh with collector-seeded semantics
- `jellyroll_element_properties`: Returns layer weights, radii, angles for each element
- Layer weights (`fks`) determine material properties per element
- See `src/Jellyrollmodel.jl`

### Cohesive Zone Model
Models interlayer degradation:
- `create_czm_mesh`: Creates cohesive elements at layer interfaces
- `bilinear_traction`: Traction-separation law
- `newton_raphson_czm`: Solve coupled electro-thermal-CZM system
- Damage state tracked per element, affects thermal conductivity
- See `src/czm.jl` and `src/CzmSolve.jl`

### Cycle Solver
Multi-cycle simulation with state management:
```julia
cycle_opt = JuBat.CycleOption(
    n_cycles=50,
    I_charge=5.0, I_discharge=5.0,
    t_charge=3600, t_discharge=3600,
    V_upper=4.2, V_lower=2.5,
    SOC_init=0.05
)
result = JuBat.solve_cycling(case, cycle_opt)
```
- Supports charge, rest, discharge phases
- State preservation between cycles
- Damage accumulation across cycles
- See `src/CycleSolver.jl`

## File Organization

### Core Model Files
- `SetCase.jl`: Setup case configuration and mesh indices
- `SetParams.jl`: Parameter normalization and scaling
- `SetMesh.jl`: Mesh generation (1D L2/L3, 2D Q4 elements)
- `Option.jl`: Option structure definitions

### Solver Files
- `Solve.jl`: Main time-stepping solver
- `SPM.jl`, `SPMe.jl`, `P2D.jl`: Model-specific implementations
- `Parallelsolution.jl`: Multi-SPMe parallel solver with current distribution
- `CzmSolve.jl`: CZM-specific solver

### Coupling Files
- `Thermal.jl`, `ThermalDistributed.jl`: Thermal model implementations
- `ThermalPolar2D.jl`: Polar coordinate thermal solver for ring geometry
- `mechanical.jl`: Stress calculations (thermal and diffusion stress)
- `czm.jl`: Cohesive zone model implementation

### Supporting Files
- `Variables.jl`: Variable initialization and storage
- `PostProcessing.jl`: Result processing and visualization
- `Tools.jl`: Utility functions
- `Materialmatrix.jl`: Material property calculations
- `Assemble.jl`: FEM assembly functions

### Parameter Files
`src/parameters/`: Cell-specific parameter sets
- `LG M50.jl`, `Enertech.jl`, `Northrop.jl`: Standard pouch cells
- `Jellyroll.jl`, `Ring.jl`: Specialized geometries

## Debugging and Verification

### Debug Options
- `opt.debug_coupling=true`: Print detailed coupling logs
- `opt.debug_log_path="output/debug.log"`: Write debug logs to file
- Debug logs include prefactors, coefficients, voltages for each element

### Verification Scripts
`tools/`: Verification and debugging scripts
- `verify_czm_system.jl`: Verify CZM implementation
- `check_branch_currents.jl`: Verify current distribution
- `check_thermal_kernels.jl`: Verify thermal model
- Various mesh checking scripts

### Documentation
`docs/`: Technical documentation and verification notes
- `SPMe_Thermal2D_Theory_vs_Code.md`: Theory-code comparison
- `Jellyroll_Thermal_Nondim.md`: Nondimensionalization details
- `CZM_Thermal_Coupling_Plan.md`: CZM coupling architecture

## Common Patterns

### Setting Up a Simulation
1. Choose cell parameters: `param_dim = JuBat.ChooseCell("Jellyroll")`
2. Configure options: `opt = JuBat.Option(); opt.model="SPMe"; opt.thermal_enabled=true`
3. Create case: `case = JuBat.SetCase(param_dim, opt)`
4. Create mesh (if thermal2D): `mesh_data = jellyroll_collector_seed_mesh(...)`
5. Solve: `result = JuBat.Solve(case)`
6. Plot: `plot(result["time [s]"], result["cell voltage [V]"])`

### Accessing Results
Results are stored as Dict with keys like:
- `time [s]`: Time vector
- `cell voltage [V]`: Cell voltage over time
- `temperature`: Temperature field (if thermal enabled)
- `negative particle lithium concentration`: Concentration profiles
- Mechanical stress fields (if mechanical enabled)
- Damage variables (if CZM enabled)

### Modifying Current Profiles
```julia
# Constant current
opt.Current = x -> 5.0  # 5A constant

# Time-dependent current
opt.Current = x -> (x < 1800) ? 5.0 : -5.0  # Charge then discharge

# Multi-step profile
opt.Current = x -> begin
    if x < 3600; return 5.0
    elseif x < 7200; return 0.0
    else; return -5.0
    end
end
```

## Important Notes

### Dimensionless vs Dimensional
- Internal calculations use dimensionless variables
- Parameters are normalized in `NormaliseParam`
- Results may need dimensionalization for interpretation
- Reference values in `param_dim.scale` (T_ref, I_ref, etc.)

### Mesh Resolution
- Higher `gsorder` (Gaussian quadrature order) = better accuracy, slower
- Thermal mesh resolution controlled by `n_theta` in jellyroll mesh
- Balance between accuracy and computational cost

### Time Stepping
- Use adaptive time stepping for stiff problems: `dtType="auto"`
- Start with small `dt_min`, allow larger `dt_max`
- Convergence issues: reduce `dt_max` or check parameters

### State Management
- For cycle simulations, use `final_state` from one phase as `initial_state` for next
- Damage state accumulates across cycles
- Temperature can be reset between cycles with `reset_T_each_cycle=true`
