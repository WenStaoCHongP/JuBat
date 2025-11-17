using Parameters

@with_kw mutable struct Option
#   option for a lithium-ion battery model
    Np::Int64 = 10
    Ns::Int64 = 10
    Nn::Int64 = 10
    Nrp::Int64 = 10
    Nrn::Int64 = 10
    model::String  = "SPM"
    time::Array{Float64} = [0 3600]
    meshType::String  = "L2"
    gsorder::Int64 = 4
    dimension::Int64 = 1
    #opt.load = {"constant discharge 1C for 1h"}
    Current::Function = x-> 0
    coupleMethod:: String  = "fully coupled"
    coupleOrder::Int64 = 0
    y0::Array{Float64} = []
    dt::Array{Float64} = [1, 100]
    dtType::String  = "constant" # auto or manual
    dtThreshold::Float64 = 0.01
    solveType::String  = "Crank-Nicolson" # forward, backward or Crank-Nicolson
    outputType::String  = "auto" # auto or manual
    jacobi::String = "constant" # constant or update
    thermalmodel::String  = "none" # none, lumped, distributed1D, distributed2D
    mechanicalmodel::String = "none" #none or full
    cite::Vector{String} = String[]
    # Thermal module flags (Phase A)
    thermal_enabled::Bool = false      # whether thermal module is active
    thermal_dim::String = "1D"        # "1D" or "2D" for distributed models
    thermalmeshType::String = "L2"    # 1D: L2/L3; 2D: Q4 (default)
    # --- New coupling/parallel options (default keep legacy behavior) ---
    collector_seeded::Bool = false     # use collector-seeded band mesh semantics (layer_weights)
    per_element_spme::Bool = false     # allow passing per-element I_app and T to SPMe
    units_thermal::String = "nd"      # "nd" (dimensionless Scheme B) or "SI"
    # --- Strong coupling (per-element SPMe) controls ---
    strong_max_iter::Int = 5           # inner iterations per thermal time step
    strong_tol_T::Float64 = 5e-3       # temperature change tol in nd (≈0.005*T_ref)
    strong_tol_Irel::Float64 = 1e-3    # relative change of I_e tolerance
    strong_relaxation::Float64 = 0.5   # relaxation factor for T (0<ω<=1)
    # --- Debug/Diagnostics ---
    debug_coupling::Bool = false       # print detailed logs for electro-thermal coupling
    debug_sample_elems::Int = 3        # number of sample elements to print per iteration
end