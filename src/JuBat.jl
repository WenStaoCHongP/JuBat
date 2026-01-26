module JuBat
using LinearAlgebra, SparseArrays, Plots, Parameters, CSV, Infiltrator, Statistics

include("Option.jl") 
include("SetMesh.jl") 
include("SetParams.jl") 
include("SetCase.jl")
include("Assemble.jl") 
include("ElectrodeDiffusion.jl")
include("ElectrolyteDiffusion.jl")
include("ElectrodePotential.jl")
include("ElectrolytePotential.jl")
include("SPM.jl") 
include("SPMe.jl") 
include("P2D.jl") 
include("Solve.jl")
include("PostProcessing.jl")
include("Tools.jl")
include("Thermal.jl")
include("ThermalDistributed.jl")
include("Variables.jl")
include("Initialisation.jl")
include("mechanical.jl")
include("Jellyrollmodel.jl")
include("czm.jl")  # 内聚力区域模型
include("CycleSolver.jl")  # 充放电循环求解器
include("ModelInitialisation_SimpleCoupling.jl")
include("CallModel_SimpleCoupling.jl")

export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, Mesh2D, GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe, SPMe_element, ModelInitialisation
export ModelInitialisation_MultiSPMe, MultiSPMe_extract_element_state, MultiSPMe_get_thermal_dofs
export MultiSPMe_update_element_state!, MultiSPMe_update_thermal_dofs!
export ModelInitialisation_SimpleCoupling, CallModel_SimpleCoupling
export extract_states_simple_coupling, compute_average_temperature
export Arrhenius, IntV
export jellyroll_spiral_params, cart2pol
export jellyroll_collector_seed_mesh, jellyroll_get_layer_weights, jellyroll_get_element_areas
export jellyroll_tab_node_indices, edge_boundary, jellyroll_element_properties
export jellyroll_element_centers
export ThermalDistributed1D, ThermalDistributed2D, ThermalDistributed2D_BC
export heatQ_Source, compute_element_heat_sources, solve_branch_currents_newton
export ThermalModel, ThermalLumpedModel, ThermalDistributed1DModel, ThermalDistributed2DModel
export thermal_diffusion_stress_2D
# CZM exports
export CohesiveElement, CohesiveMesh, DamageState, CZMResult
export create_czm_mesh, identify_interface_node_pairs
export bilinear_traction, bilinear_tangent, update_damage!
export cohesive_element_matrices, compute_separation
export assemble_czm_system, assemble_coupled_system, assemble_bulk_stiffness
export assemble_thermal_chemical_load, assemble_coupled_system_full
export apply_bc_czm!, identify_bc_nodes_czm
export newton_raphson_czm, solve_czm_step
export get_damage_statistics, check_fracture_criterion, reset_damage_states!, accumulate_cycle_damage!
export czm_output_to_variables
# Cycle solver exports
export CycleOption, PhaseType, PHASE_CHARGE, PHASE_REST, PHASE_DISCHARGE
export PhaseResult, CycleResult, CyclingResult
export solve_phase, solve_cycling, plot_cycling_results
export compute_cs0_from_soc, apply_initial_soc!
end