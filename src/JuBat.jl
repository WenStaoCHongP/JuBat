module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator, Statistics, Printf

include("Option.jl") 
include("SetMesh.jl") 
include("SetParams.jl") 
include("SetCase.jl")
include("czm.jl")  # 内聚力区域模型
include("CzmSolve.jl")
include("Assemble.jl") 
include("ElectrodeDiffusion.jl")
include("ElectrolyteDiffusion.jl")
include("ElectrodePotential.jl")
include("ElectrolytePotential.jl")
include("SPM.jl") 
include("SPMe.jl") 
include("P2D.jl") 
include("Parallelsolution.jl")
include("Tools.jl")
include("Solve.jl")
include("PostProcessing.jl")
include("Materialmatrix.jl")
include("Thermal.jl")
include("ThermalDistributed.jl")
include("ThermalPolar2D.jl")
include("Variables.jl")
include("Initialisation.jl")
include("mechanical.jl")
include("Jellyrollmodel.jl")
include("ring.jl")
include("CycleSolver.jl")  # 充放电循环求解器
include("CycleData.jl")

export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, Mesh2D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe, SPMe_element, ModelInitialisation
export ModelInitialisation_MultiSPMe, MultiSPMe_extract_element_state, MultiSPMe_get_thermal_dofs
export MultiSPMe_update_state
export Arrhenius, IntV, IntQ4
export jellyroll_collector_seed_mesh, jellyroll_element_properties
export jellyroll_tab_node_indices, edge_boundary
export jellyroll_element_centers
export ring_mesh
export setup_thermal2D_mesh
export thermal2D_volume_average_temperature
export ThermalDistributed2D, ThermalDistributed2D_BC
export ThermalDistributed2D_Ring, ThermalRing2D_BC
export ThermalPolar2D_Ring
export identify_boundary_nodes, apply_convection_bc, apply_cool_method
export compute_heat_sources, compute_heat_sources_with_czm, solve_branch_currents_newton
export ThermalModel, ThermalLumpedModel, ThermalDistributed2DModel
export thermal_diffusion_stress_2D
# CZM exports
export CohesiveElement, CohesiveMesh, DamageState, CZMResult
export create_czm_mesh, compute_separation
export bilinear_traction, bilinear_tangent, update_damage
export assemble_czm_system, assemble_coupled_system, assemble_bulk_stiffness
export assemble_thermal_chemical_load, assemble_coupled_system_full
export apply_bc_czm, identify_bc_nodes_czm
export newton_raphson_czm, solve_czm_step
export get_damage_statistics, check_fracture_criterion, reset_damage_states, accumulate_cycle_damage
export czm_output_to_variables
# CZM gap conductance exports
export compute_gap_conductance, compute_element_gap_conductance
export get_fractured_elements, get_active_elements, compute_all_gap_conductances
# Cycle solver exports
export CycleOption, PhaseType, PHASE_CHARGE, PHASE_REST, PHASE_DISCHARGE
export PhaseResult, CycleResult, CyclingResult
export solve_phase, solve_cycling, plot_cycling_results
export compute_cs0_from_soc, apply_initial_soc!
# Cycle data export/import
export TimeStepData, CycleExportData
export solve_phase_with_export, solve_cycle_with_export
export export_cycle_data_to_csv, load_cycle_data_from_csv
end