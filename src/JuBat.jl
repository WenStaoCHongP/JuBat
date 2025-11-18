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

export Assemble, ElectrodeDiffusion, ElectrolyteDiffusion, Postprocessing, SetCase, SetMesh, ChooseCell
export Mesh1D, Mesh2D,GetGS, LagrangeBasis, GSweight, ShapeFunction1D, NormaliseParam, StandardVariables
export SPM, Solve, SPMe, SPMe_element, ModelInitialisation
export Arrhenius, IntV
export jellyroll_spiral_params, cart2pol, pol2cart, material_at
export jellyroll_element_layer_weights, get_element_layer_weights
export ThermalDistributed1D, ThermalDistributed2D, ThermalDistributed_BC, heatQ_Source, solve_branch_currents_newton
export ThermalModel, ThermalLumpedModel, ThermalDistributed1DModel, ThermalDistributed2DModel
export thermal_stress
end