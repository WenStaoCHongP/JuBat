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
include("cohesive_contact.jl")
include("Jellyrollmodel.jl")

export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, Mesh2D,GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe, SPMe_element, ModelInitialisation
export ModelInitialisation_MultiSPMe, MultiSPMe_extract_element_state, MultiSPMe_get_thermal_dofs
export MultiSPMe_update_element_state!, MultiSPMe_update_thermal_dofs!
export Arrhenius, IntV
export jellyroll_spiral_params, cart2pol, material_at
export jellyroll_collector_seed_mesh, jellyroll_get_layer_weights
export jellyroll_tab_node_indices, edge_boundary
export jellyroll_element_centers, jellyroll_effective_K_at
export ThermalDistributed1D, ThermalDistributed2D, ThermalDistributed_BC, heatQ_Source, solve_branch_currents_newton
export ThermalModel, ThermalLumpedModel, ThermalDistributed1DModel, ThermalDistributed2DModel
export thermal_stress, homogenize_particle_stress_to_2D
export CohesiveInterface, ContactInterface
export init_cohesive_interface, compute_cohesive_traction, update_cohesive_damage!, cohesive_stiffness_matrix
export init_contact_interface, detect_contact, compute_contact_pressure, compute_friction_stress
export cohesive_output, contact_output, hertz_contact_pressure, effective_modulus
end