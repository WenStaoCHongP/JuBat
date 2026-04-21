include(joinpath(@__DIR__, "..", "..", "..", "src", "JuBat.jl"))
using .JuBat
using Printf

param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
param_dim.cell.v_h = 4.2

opt = JuBat.Option()
opt.Current = x -> 5.0
opt.model = "SPMe"
opt.Nn = 10
opt.Ns = 5
opt.Np = 10
opt.Nrn = 10
opt.Nrp = 10
opt.gsorder = 2
opt.dimension = 1
opt.mechanicalmodel = "none"
opt.time = [0.0, 3600.0]
opt.dt = [0.5, 10.0]
opt.dtType = "auto"
opt.jacobi = "update"
opt.solveType = "Crank-Nicolson"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.cool_method = "surface"
opt.per_element_spme = true
opt.debug_coupling = false
opt.czm_enabled = true

case = JuBat.SetCase(param_dim, opt)

function print_row(label, δ_n, δ_t, state, coh)
    T_n, T_t, D, new_state = JuBat.bilinear_traction_state(δ_n, δ_t, state, coh)
    K = JuBat.bilinear_tangent(δ_n, δ_t, new_state, coh)
    @printf("%-18s δ_n=% .6e  D=% .6e  fractured=%s  T_n=% .6e  K_nn=% .6e  K_tt=% .6e\n",
        label, δ_n, D, string(new_state.fractured), T_n, K[1,1], K[2,2])
    return new_state
end

function main()
    coh = case.param.cohesive
    state = JuBat.DamageState()
    prefracture_state = nothing

    println("=== Loading sweep ===")
    δ0 = coh.δ_0_n
    δc = coh.δ_c_n
    for δ in [0.0, 0.25*δ0, 0.5*δ0, 0.9*δ0, 0.99*δ0, 1.01*δ0, 1.1*δ0, 0.5*δc, 0.9*δc, 0.99*δc, 1.01*δc, 1.1*δc]
        state = print_row("load", δ, 0.0, state, coh)
        if abs(δ - 0.99*δc) < 1e-12
            prefracture_state = state
        end
    end

    println("=== Damaged but not fractured zero-crossing ===")
    if prefracture_state !== nothing
        state2 = prefracture_state
        for δ in [-1e-6, -1e-9, 0.0, 1e-9, 1e-6]
            state2 = print_row("near-zero", δ, 0.0, state2, coh)
        end
    end

    println("=== Unloading / closing sweep ===")
    for δ in [0.5*δc, 0.2*δc, 0.0, -1e-9, -1e-6]
        state = print_row("unload", δ, 0.0, state, coh)
    end
end

main()
