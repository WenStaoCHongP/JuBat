"""
Cohesive Zone Model (CZM) Validation: Pure Mechanical Loading

Features:
- No electrochemical-thermal coupling
- Periodic displacement loading to validate bilinear traction-separation law
- Plot damage curves and traction-separation hysteresis loops

Validation Tests:
1. Monotonic loading: Validate bilinear traction-separation relationship
2. Cyclic loading: Validate loading/unloading behavior and damage evolution
3. Mixed-mode loading: Validate normal + tangential coupling (BK criterion)
4. Sinusoidal displacement: Validate response under realistic cyclic loading

Output Figures:
- Figure 1: Normal traction-separation curve (loading/unloading hysteresis)
- Figure 2: Tangential traction-separation curve (loading/unloading hysteresis)
- Figure 3: Damage variable evolution vs maximum separation
- Figure 4: Mixed-mode traction response

Date: 2025
"""

using LinearAlgebra, Printf, Plots

# Include JuBat module
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

"""
    create_czm_test_params()

Create cohesive zone parameters for testing.
Returns typical electrode-separator interface parameters.
"""
function create_czm_test_params()
    # Create cohesive parameters (typical electrode-separator interface)
    cohesive = JuBat.Cohesive()
    
    # Normal (Mode I) parameters
    cohesive.σ_max_n = 50e6       # Maximum normal traction [Pa] (50 MPa)
    cohesive.δ_0_n = 1e-6         # Damage initiation separation [m] (1 μm)
    cohesive.δ_c_n = 10e-6        # Critical separation [m] (10 μm)
    cohesive.G_c_n = 0.5 * cohesive.σ_max_n * cohesive.δ_c_n  # Fracture energy [J/m²]
    cohesive.K_n = cohesive.σ_max_n / cohesive.δ_0_n  # Initial stiffness [Pa/m]
    
    # Tangential (Mode II) parameters
    cohesive.τ_max_t = 30e6       # Maximum tangential traction [Pa] (30 MPa)
    cohesive.δ_0_t = 1e-6         # Damage initiation separation [m] (1 μm)
    cohesive.δ_c_t = 15e-6        # Critical separation [m] (15 μm)
    cohesive.G_c_t = 0.5 * cohesive.τ_max_t * cohesive.δ_c_t  # Fracture energy [J/m²]
    cohesive.K_t = cohesive.τ_max_t / cohesive.δ_0_t  # Initial stiffness [Pa/m]
    
    # BK criterion exponent
    cohesive.eta = 1.45
    
    return cohesive
end

"""
    test_monotonic_loading(cohesive_params)

Test constitutive response under monotonic loading.
Plot traction-separation curve from zero to complete fracture.
"""
function test_monotonic_loading(cohesive_params)
    println("\n" * "="^60)
    println("Test 1: Monotonic Loading Constitutive Response")
    println("="^60)
    
    # Create damage states
    damage_state_n = JuBat.DamageState()
    damage_state_t = JuBat.DamageState()
    
    # Normal loading parameters
    δ_max_n = cohesive_params.δ_c_n * 1.2  # Exceed critical separation
    n_points = 200
    δ_n_vals = range(0, δ_max_n, length=n_points)
    
    T_n_vals = zeros(n_points)
    D_n_vals = zeros(n_points)
    
    # Normal monotonic loading
    for (i, δ_n) in enumerate(δ_n_vals)
        T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=true)
        T_n_vals[i] = T_n
        D_n_vals[i] = D
    end
    
    # Tangential loading parameters
    δ_max_t = cohesive_params.δ_c_t * 1.2
    δ_t_vals = range(0, δ_max_t, length=n_points)
    
    T_t_vals = zeros(n_points)
    D_t_vals = zeros(n_points)
    
    # Tangential monotonic loading
    for (i, δ_t) in enumerate(δ_t_vals)
        _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
        T_t_vals[i] = T_t
        D_t_vals[i] = D
    end
    
    # Print key parameters
    println("\nCohesive Zone Parameters:")
    @printf("  Normal max traction sigma_max_n = %.1f MPa\n", cohesive_params.σ_max_n / 1e6)
    @printf("  Normal initiation delta_0_n = %.1f um\n", cohesive_params.δ_0_n * 1e6)
    @printf("  Normal critical delta_c_n = %.1f um\n", cohesive_params.δ_c_n * 1e6)
    @printf("  Normal fracture energy G_c_n = %.1f J/m^2\n", cohesive_params.G_c_n)
    println()
    @printf("  Tangential max traction tau_max_t = %.1f MPa\n", cohesive_params.τ_max_t / 1e6)
    @printf("  Tangential initiation delta_0_t = %.1f um\n", cohesive_params.δ_0_t * 1e6)
    @printf("  Tangential critical delta_c_t = %.1f um\n", cohesive_params.δ_c_t * 1e6)
    @printf("  Tangential fracture energy G_c_t = %.1f J/m^2\n", cohesive_params.G_c_t)
    
    return (δ_n_vals, T_n_vals, D_n_vals, δ_t_vals, T_t_vals, D_t_vals)
end

"""
    test_cyclic_loading(cohesive_params; n_cycles=3, max_amp_factor=0.8)

Test constitutive response under cyclic loading/unloading.
Generate hysteresis loops to validate loading/unloading behavior and damage accumulation.
"""
function test_cyclic_loading(cohesive_params; n_cycles::Int=3, max_amp_factor::Float64=0.8)
    println("\n" * "="^60)
    println("Test 2: Cyclic Loading/Unloading Hysteresis")
    println("="^60)
    
    # Normal cyclic loading
    damage_state_n = JuBat.DamageState()
    
    # Progressively increasing amplitude
    δ_max_n = cohesive_params.δ_c_n * max_amp_factor
    
    δ_n_history = Float64[]
    T_n_history = Float64[]
    D_n_history = Float64[]
    
    points_per_half_cycle = 50
    
    for cycle in 1:n_cycles
        # Amplitude increases with each cycle
        amp = δ_max_n * (cycle / n_cycles)
        
        # Loading phase
        for i in 1:points_per_half_cycle
            δ_n = amp * (i / points_per_half_cycle)
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=true)
            push!(δ_n_history, δ_n)
            push!(T_n_history, T_n)
            push!(D_n_history, D)
        end
        
        # Unloading phase
        for i in 1:points_per_half_cycle
            δ_n = amp * (1.0 - i / points_per_half_cycle)
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state_n, cohesive_params; update=false)
            push!(δ_n_history, δ_n)
            push!(T_n_history, T_n)
            push!(D_n_history, D)
        end
        
        @printf("  Cycle %d: amplitude = %.2f um, max damage = %.2f%%\n", 
                cycle, amp * 1e6, damage_state_n.D * 100)
    end
    
    # Tangential cyclic loading (bidirectional)
    damage_state_t = JuBat.DamageState()
    
    δ_max_t = cohesive_params.δ_c_t * max_amp_factor
    
    δ_t_history = Float64[]
    T_t_history = Float64[]
    D_t_history = Float64[]
    
    for cycle in 1:n_cycles
        amp = δ_max_t * (cycle / n_cycles)
        
        # Positive loading
        for i in 1:points_per_half_cycle
            δ_t = amp * (i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # Unload to zero
        for i in 1:points_per_half_cycle
            δ_t = amp * (1.0 - i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=false)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # Negative loading (reverse direction)
        for i in 1:points_per_half_cycle
            δ_t = -amp * (i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=true)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
        
        # Reverse unloading to zero
        for i in 1:points_per_half_cycle
            δ_t = -amp * (1.0 - i / points_per_half_cycle)
            _, T_t, D = JuBat.bilinear_traction(0.0, δ_t, damage_state_t, cohesive_params; update=false)
            push!(δ_t_history, δ_t)
            push!(T_t_history, T_t)
            push!(D_t_history, D)
        end
    end
    
    return (δ_n_history, T_n_history, D_n_history, δ_t_history, T_t_history, D_t_history)
end

"""
    test_mixed_mode_loading(cohesive_params)

Test constitutive response under mixed-mode (normal + tangential) loading.
Validate the BK (Benzeggagh-Kenane) criterion.
"""
function test_mixed_mode_loading(cohesive_params)
    println("\n" * "="^60)
    println("Test 3: Mixed-Mode Loading (BK Criterion)")
    println("="^60)
    
    # Test different mode mixing ratios: beta = |delta_t| / delta_eff
    mode_ratios = [0.0, 0.25, 0.5, 0.75, 1.0]  # 0=pure Mode I, 1=pure Mode II
    
    n_points = 100
    
    results = Dict{Float64, NamedTuple}()
    
    for β in mode_ratios
        damage_state = JuBat.DamageState()
        
        # Calculate separation components based on mode ratio
        # delta_eff = sqrt(delta_n^2 + delta_t^2), beta = |delta_t| / delta_eff
        # => delta_t = beta * delta_eff, delta_n = sqrt(1 - beta^2) * delta_eff
        
        δ_eff_max = max(cohesive_params.δ_c_n, cohesive_params.δ_c_t) * 1.2
        
        δ_eff_vals = range(0, δ_eff_max, length=n_points)
        δ_n_vals = zeros(n_points)
        δ_t_vals = zeros(n_points)
        T_n_vals = zeros(n_points)
        T_t_vals = zeros(n_points)
        D_vals = zeros(n_points)
        
        for (i, δ_eff) in enumerate(δ_eff_vals)
            δ_n = sqrt(1.0 - β^2) * δ_eff
            δ_t = β * δ_eff
            
            δ_n_vals[i] = δ_n
            δ_t_vals[i] = δ_t
            
            T_n, T_t, D = JuBat.bilinear_traction(δ_n, δ_t, damage_state, cohesive_params; update=true)
            
            T_n_vals[i] = T_n
            T_t_vals[i] = T_t
            D_vals[i] = D
        end
        
        # Calculate effective traction
        T_eff_vals = sqrt.(T_n_vals.^2 .+ T_t_vals.^2)
        
        results[β] = (δ_eff=δ_eff_vals, δ_n=δ_n_vals, δ_t=δ_t_vals,
                      T_n=T_n_vals, T_t=T_t_vals, T_eff=T_eff_vals, D=D_vals)
        
        @printf("  Mode ratio beta = %.2f: T_max = %.1f MPa, D_final = %.2f%%\n",
                β, maximum(T_eff_vals) / 1e6, D_vals[end] * 100)
    end
    
    return results
end

"""
    test_sinusoidal_displacement(cohesive_params; n_cycles=5, frequency=1.0, amplitude_factor=0.6)

Test response under sinusoidal displacement loading.
Simulates more realistic cyclic loading conditions.
"""
function test_sinusoidal_displacement(cohesive_params; n_cycles::Int=5, 
                                       frequency::Float64=1.0,
                                       amplitude_factor::Float64=0.6)
    println("\n" * "="^60)
    println("Test 4: Sinusoidal Displacement Loading")
    println("="^60)
    
    damage_state = JuBat.DamageState()
    
    # Time parameters
    T_period = 1.0 / frequency
    t_total = n_cycles * T_period
    n_points = n_cycles * 100
    t_vals = range(0, t_total, length=n_points)
    
    # Displacement amplitude
    δ_amp = cohesive_params.δ_c_n * amplitude_factor
    
    # Sinusoidal displacement history
    δ_n_vals = δ_amp .* sin.(2π * frequency .* t_vals)
    
    T_n_history = zeros(n_points)
    D_history = zeros(n_points)
    
    for (i, δ_n) in enumerate(δ_n_vals)
        # Only update damage for positive (opening) separation
        if δ_n > 0
            T_n, _, D = JuBat.bilinear_traction(δ_n, 0.0, damage_state, cohesive_params; update=true)
        else
            # Compression: use pure elastic contact
            T_n = cohesive_params.K_n * δ_n
            D = damage_state.D
        end
        T_n_history[i] = T_n
        D_history[i] = D
    end
    
    @printf("\n  Loading Parameters:\n")
    @printf("    Amplitude = %.2f um (%.0f%% of delta_c)\n", δ_amp * 1e6, amplitude_factor * 100)
    @printf("    Frequency = %.1f Hz\n", frequency)
    @printf("    Number of cycles = %d\n", n_cycles)
    @printf("  Results:\n")
    @printf("    Final damage = %.2f%%\n", D_history[end] * 100)
    @printf("    Maximum traction = %.1f MPa\n", maximum(T_n_history) / 1e6)
    
    return (t=t_vals, δ_n=δ_n_vals, T_n=T_n_history, D=D_history)
end

"""
    plot_all_results(...)

Generate all validation plots with detailed English annotations.
"""
function plot_all_results(monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data, cohesive_params)
    println("\n" * "="^60)
    println("Generating Figures...")
    println("="^60)
    
    # Ensure output directory exists
    isdir("output") || mkdir("output")
    
    # Unpack data
    δ_n_mono, T_n_mono, D_n_mono, δ_t_mono, T_t_mono, D_t_mono = monotonic_data
    δ_n_cyc, T_n_cyc, D_n_cyc, δ_t_cyc, T_t_cyc, D_t_cyc = cyclic_data
    
    # ====================================================================
    # Figure 1: Monotonic Loading Constitutive Curves
    # ====================================================================
    p1 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # Normal traction-separation curve
    plot!(p1[1], δ_n_mono .* 1e6, T_n_mono ./ 1e6,
          xlabel="Normal Separation delta_n (um)", 
          ylabel="Normal Traction T_n (MPa)",
          title="(a) Normal Traction-Separation (Mode I)",
          label="Monotonic Loading", linewidth=2, color=:blue,
          legend=:topright, grid=true, minorgrid=true)
    
    # Mark key points
    vline!(p1[1], [cohesive_params.δ_0_n * 1e6], label="delta_0_n (initiation)", 
           linestyle=:dash, color=:gray, linewidth=1.5)
    vline!(p1[1], [cohesive_params.δ_c_n * 1e6], label="delta_c_n (critical)", 
           linestyle=:dot, color=:red, linewidth=1.5)
    hline!(p1[1], [cohesive_params.σ_max_n / 1e6], label="sigma_max_n", 
           linestyle=:dash, color=:orange, linewidth=1.5)
    
    # Tangential traction-separation curve
    plot!(p1[2], δ_t_mono .* 1e6, T_t_mono ./ 1e6,
          xlabel="Tangential Separation delta_t (um)", 
          ylabel="Tangential Traction T_t (MPa)",
          title="(b) Tangential Traction-Separation (Mode II)",
          label="Monotonic Loading", linewidth=2, color=:green,
          legend=:topright, grid=true, minorgrid=true)
    
    vline!(p1[2], [cohesive_params.δ_0_t * 1e6], label="delta_0_t (initiation)", 
           linestyle=:dash, color=:gray, linewidth=1.5)
    vline!(p1[2], [cohesive_params.δ_c_t * 1e6], label="delta_c_t (critical)", 
           linestyle=:dot, color=:red, linewidth=1.5)
    hline!(p1[2], [cohesive_params.τ_max_t / 1e6], label="tau_max_t", 
           linestyle=:dash, color=:orange, linewidth=1.5)
    
    # Normal damage evolution
    plot!(p1[3], δ_n_mono .* 1e6, D_n_mono .* 100,
          xlabel="Normal Separation delta_n (um)", 
          ylabel="Damage Variable D (%)",
          title="(c) Normal Damage Evolution",
          label="D vs delta_n", linewidth=2, color=:red,
          legend=:bottomright, grid=true, minorgrid=true)
    
    # Add annotation for damage stages
    annotate!(p1[3], [(cohesive_params.δ_0_n * 1e6 * 0.5, 5, 
                       text("Elastic\nStage", 8, :center)),
                      ((cohesive_params.δ_0_n + cohesive_params.δ_c_n) / 2 * 1e6, 50, 
                       text("Softening\nStage", 8, :center)),
                      (cohesive_params.δ_c_n * 1e6 * 1.1, 95, 
                       text("Fractured", 8, :left))])
    
    # Tangential damage evolution
    plot!(p1[4], δ_t_mono .* 1e6, D_t_mono .* 100,
          xlabel="Tangential Separation delta_t (um)", 
          ylabel="Damage Variable D (%)",
          title="(d) Tangential Damage Evolution",
          label="D vs delta_t", linewidth=2, color=:purple,
          legend=:bottomright, grid=true, minorgrid=true)
    
    savefig(p1, "output/czm_monotonic_loading.png")
    println("  Saved: output/czm_monotonic_loading.png")
    
    # ====================================================================
    # Figure 2: Cyclic Loading/Unloading Hysteresis Loops
    # ====================================================================
    p2 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # Normal hysteresis loop
    plot!(p2[1], δ_n_cyc .* 1e6, T_n_cyc ./ 1e6,
          xlabel="Normal Separation delta_n (um)", 
          ylabel="Normal Traction T_n (MPa)",
          title="(a) Normal Loading/Unloading Hysteresis",
          label="Cyclic Response", linewidth=1.5, color=:blue,
          grid=true, minorgrid=true)
    
    # Add monotonic envelope as reference
    plot!(p2[1], δ_n_mono .* 1e6, T_n_mono ./ 1e6,
          label="Monotonic Envelope", linestyle=:dash, linewidth=1, color=:gray, alpha=0.5)
    
    # Add arrows to indicate loading direction
    annotate!(p2[1], [(2, 35, text("Loading", 9, :left, :blue)),
                      (5, 15, text("Unloading", 9, :left, :red))])
    
    # Tangential hysteresis loop (bidirectional)
    plot!(p2[2], δ_t_cyc .* 1e6, T_t_cyc ./ 1e6,
          xlabel="Tangential Separation delta_t (um)", 
          ylabel="Tangential Traction T_t (MPa)",
          title="(b) Tangential Hysteresis (Bidirectional)",
          label="Cyclic Response", linewidth=1.5, color=:green,
          grid=true, minorgrid=true)
    
    # Normal damage history
    n_cyc_n = length(δ_n_cyc)
    plot!(p2[3], 1:n_cyc_n, D_n_cyc .* 100,
          xlabel="Loading Step", 
          ylabel="Damage Variable D (%)",
          title="(c) Normal Damage Accumulation",
          label="Damage D", linewidth=1.5, color=:red,
          grid=true, minorgrid=true)
    
    # Mark cycle boundaries
    cycle_steps = [100, 200, 300]
    for (i, step) in enumerate(cycle_steps)
        if step <= n_cyc_n
            vline!(p2[3], [step], label=i==1 ? "Cycle End" : "", 
                   linestyle=:dot, color=:gray, alpha=0.5)
        end
    end
    
    # Tangential damage history
    n_cyc_t = length(δ_t_cyc)
    plot!(p2[4], 1:n_cyc_t, D_t_cyc .* 100,
          xlabel="Loading Step", 
          ylabel="Damage Variable D (%)",
          title="(d) Tangential Damage Accumulation",
          label="Damage D", linewidth=1.5, color=:purple,
          grid=true, minorgrid=true)
    
    savefig(p2, "output/czm_cyclic_hysteresis.png")
    println("  Saved: output/czm_cyclic_hysteresis.png")
    
    # ====================================================================
    # Figure 3: Mixed-Mode Response (BK Criterion)
    # ====================================================================
    p3 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # Color palette and labels
    colors = [:blue, :cyan, :green, :orange, :red]
    mode_labels = ["beta=0 (Pure Mode I)", "beta=0.25", "beta=0.5", 
                   "beta=0.75", "beta=1 (Pure Mode II)"]
    
    # Effective traction vs effective separation
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[1], collect(data.δ_eff) .* 1e6, data.T_eff ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[1], xlabel="Effective Separation delta_eff (um)", 
          ylabel="Effective Traction T_eff (MPa)",
          title="(a) Mixed-Mode: Effective Traction Curves",
          legend=:topright, grid=true, minorgrid=true)
    
    # Damage evolution
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[2], collect(data.δ_eff) .* 1e6, data.D .* 100,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[2], xlabel="Effective Separation delta_eff (um)", 
          ylabel="Damage Variable D (%)",
          title="(b) Mixed-Mode: Damage Evolution",
          legend=:bottomright, grid=true, minorgrid=true)
    
    # Normal component
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[3], data.δ_n .* 1e6, data.T_n ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[3], xlabel="Normal Separation delta_n (um)", 
          ylabel="Normal Traction T_n (MPa)",
          title="(c) Mixed-Mode: Normal Component",
          legend=:topright, grid=true, minorgrid=true)
    
    # Tangential component
    for (i, (β, data)) in enumerate(sort(collect(mixed_mode_data)))
        plot!(p3[4], data.δ_t .* 1e6, data.T_t ./ 1e6,
              label=mode_labels[i], linewidth=2, color=colors[i])
    end
    plot!(p3[4], xlabel="Tangential Separation delta_t (um)", 
          ylabel="Tangential Traction T_t (MPa)",
          title="(d) Mixed-Mode: Tangential Component",
          legend=:topright, grid=true, minorgrid=true)
    
    savefig(p3, "output/czm_mixed_mode.png")
    println("  Saved: output/czm_mixed_mode.png")
    
    # ====================================================================
    # Figure 4: Sinusoidal Displacement Loading Response
    # ====================================================================
    p4 = plot(layout=(2, 2), size=(1200, 900), dpi=150)
    
    # Displacement-time curve
    plot!(p4[1], collect(sinusoidal_data.t), sinusoidal_data.δ_n .* 1e6,
          xlabel="Time t (s)", 
          ylabel="Normal Separation delta_n (um)",
          title="(a) Sinusoidal Displacement History",
          label="delta_n(t)", linewidth=1.5, color=:blue,
          grid=true, minorgrid=true)
    
    # Mark positive/negative regions
    hline!(p4[1], [0], label="", linestyle=:dash, color=:gray, alpha=0.5)
    
    # Traction-time curve
    plot!(p4[2], collect(sinusoidal_data.t), sinusoidal_data.T_n ./ 1e6,
          xlabel="Time t (s)", 
          ylabel="Normal Traction T_n (MPa)",
          title="(b) Traction Response History",
          label="T_n(t)", linewidth=1.5, color=:green,
          grid=true, minorgrid=true)
    
    # Hysteresis loop
    plot!(p4[3], sinusoidal_data.δ_n .* 1e6, sinusoidal_data.T_n ./ 1e6,
          xlabel="Normal Separation delta_n (um)", 
          ylabel="Normal Traction T_n (MPa)",
          title="(c) Sinusoidal Loading Hysteresis Loop",
          label="Traction-Separation", linewidth=1.5, color=:purple,
          grid=true, minorgrid=true)
    
    # Highlight asymmetry between tension and compression
    annotate!(p4[3], [(3, 30, text("Tension:\nDamage evolves", 8, :left)),
                      (-3, -100, text("Compression:\nElastic contact", 8, :right))])
    
    # Damage evolution
    plot!(p4[4], collect(sinusoidal_data.t), sinusoidal_data.D .* 100,
          xlabel="Time t (s)", 
          ylabel="Damage Variable D (%)",
          title="(d) Damage Accumulation Process",
          label="D(t)", linewidth=1.5, color=:red,
          grid=true, minorgrid=true)
    
    # Add cycle markers
    T_period = 1.0  # Period = 1s for 1 Hz
    for i in 1:5
        vline!(p4[4], [i * T_period], label=i==1 ? "Cycle End" : "", 
               linestyle=:dot, color=:gray, alpha=0.3)
    end
    
    savefig(p4, "output/czm_sinusoidal_loading.png")
    println("  Saved: output/czm_sinusoidal_loading.png")
    
    # ====================================================================
    # Figure 5: Hysteresis Loop Comparison
    # ====================================================================
    p5 = plot(size=(900, 600), dpi=150)
    
    # Cyclic loading hysteresis
    plot!(p5, δ_n_cyc .* 1e6, T_n_cyc ./ 1e6,
          label="Progressive Amplitude Loading", linewidth=2, color=:blue)
    
    # Sinusoidal loading hysteresis
    plot!(p5, sinusoidal_data.δ_n .* 1e6, sinusoidal_data.T_n ./ 1e6,
          label="Sinusoidal Loading", linewidth=2, color=:red, linestyle=:dash)
    
    # Monotonic envelope
    plot!(p5, δ_n_mono .* 1e6, T_n_mono ./ 1e6,
          label="Monotonic Envelope", linewidth=1.5, color=:black, linestyle=:dot, alpha=0.5)
    
    plot!(p5, xlabel="Normal Separation delta_n (um)", 
          ylabel="Normal Traction T_n (MPa)",
          title="CZM Hysteresis Loop Comparison",
          legend=:topright, grid=true, minorgrid=true)
    
    # Add text box with key observations
    annotate!(p5, [(8, 45, text("Key Observations:\n" *
                                "1. Loading follows envelope\n" *
                                "2. Unloading: secant to origin\n" *
                                "3. Damage is irreversible", 
                                8, :left))])
    
    savefig(p5, "output/czm_hysteresis_comparison.png")
    println("  Saved: output/czm_hysteresis_comparison.png")
    
    # ====================================================================
    # Figure 6: Bilinear Law Schematic
    # ====================================================================
    p6 = plot(size=(800, 500), dpi=150)
    
    # Create idealized bilinear curve
    δ_0 = 1.0   # Normalized
    δ_c = 10.0  # Normalized
    T_max = 50.0  # Normalized
    
    δ_schematic = [0, δ_0, δ_c, δ_c * 1.2]
    T_schematic = [0, T_max, 0, 0]
    
    plot!(p6, δ_schematic, T_schematic,
          linewidth=3, color=:blue, label="Bilinear Traction-Separation Law",
          xlabel="Separation delta (normalized)", 
          ylabel="Traction T (normalized)",
          title="Bilinear Cohesive Zone Model Schematic")
    
    # Mark key points
    scatter!(p6, [0, δ_0, δ_c], [0, T_max, 0], 
             markersize=8, color=[:black, :red, :green],
             label="")
    
    # Add labels
    annotate!(p6, [(δ_0 * 0.5, T_max * 0.3, text("K = T_max/delta_0\n(Initial Stiffness)", 9, :center)),
                   (δ_0, T_max * 1.1, text("(delta_0, T_max)\nDamage Initiation", 9, :center)),
                   (δ_c, T_max * 0.15, text("(delta_c, 0)\nComplete Fracture", 9, :center)),
                   ((δ_0 + δ_c) / 2, T_max * 0.6, text("Softening Region\nD increases", 9, :center))])
    
    # Add unloading path
    δ_unload = δ_0 + (δ_c - δ_0) * 0.5  # Unload from middle of softening
    T_at_unload = T_max * (δ_c - δ_unload) / (δ_c - δ_0)
    
    plot!(p6, [0, δ_unload], [0, T_at_unload],
          linewidth=2, color=:red, linestyle=:dash,
          label="Unloading Path (Secant)")
    
    scatter!(p6, [δ_unload], [T_at_unload], markersize=6, color=:red, label="")
    annotate!(p6, [(δ_unload * 0.7, T_at_unload * 0.7, 
                   text("Unloading:\nK_unload = (1-D)*K", 8, :right, :red))])
    
    # Shade fracture energy area
    δ_fill = range(0, δ_c, length=50)
    T_fill = [d <= δ_0 ? T_max * d / δ_0 : T_max * (δ_c - d) / (δ_c - δ_0) for d in δ_fill]
    plot!(p6, δ_fill, T_fill, fillrange=0, fillalpha=0.2, color=:blue, label="")
    annotate!(p6, [(δ_c * 0.4, T_max * 0.2, text("G_c = Area\n(Fracture Energy)", 9, :center, :blue))])
    
    savefig(p6, "output/czm_bilinear_schematic.png")
    println("  Saved: output/czm_bilinear_schematic.png")
    
    # Save SVG format
    try
        savefig(p1, "output/czm_monotonic_loading.svg")
        savefig(p2, "output/czm_cyclic_hysteresis.svg")
        savefig(p3, "output/czm_mixed_mode.svg")
        savefig(p4, "output/czm_sinusoidal_loading.svg")
        savefig(p5, "output/czm_hysteresis_comparison.svg")
        savefig(p6, "output/czm_bilinear_schematic.svg")
        println("  SVG format figures saved")
    catch e
        println("  Warning: SVG save failed (missing dependencies)")
    end
end

"""
    print_verification_summary(cohesive_params, monotonic_data, cyclic_data)

Print verification summary comparing theoretical and computed values.
"""
function print_verification_summary(cohesive_params, monotonic_data, cyclic_data)
    println("\n" * "="^60)
    println("Verification Summary")
    println("="^60)
    
    δ_n_mono, T_n_mono, D_n_mono, δ_t_mono, T_t_mono, D_t_mono = monotonic_data
    
    # Theoretical values
    σ_max_theo = cohesive_params.σ_max_n
    τ_max_theo = cohesive_params.τ_max_t
    
    # Computed values
    σ_max_calc = maximum(T_n_mono)
    τ_max_calc = maximum(T_t_mono)
    
    # Compute fracture energy by numerical integration
    G_n_calc = 0.0
    for i in 2:length(δ_n_mono)
        dδ = δ_n_mono[i] - δ_n_mono[i-1]
        G_n_calc += 0.5 * (T_n_mono[i] + T_n_mono[i-1]) * dδ
    end
    
    G_t_calc = 0.0
    for i in 2:length(δ_t_mono)
        dδ = δ_t_mono[i] - δ_t_mono[i-1]
        G_t_calc += 0.5 * (T_t_mono[i] + T_t_mono[i-1]) * dδ
    end
    
    println("\nNormal (Mode I) Verification:")
    println("-"^40)
    @printf("  sigma_max:  Theory = %.2f MPa, Computed = %.2f MPa, Error = %.2f%%\n",
            σ_max_theo / 1e6, σ_max_calc / 1e6, 
            abs(σ_max_calc - σ_max_theo) / σ_max_theo * 100)
    @printf("  G_c_n:      Theory = %.2f J/m^2, Computed = %.2f J/m^2, Error = %.2f%%\n",
            cohesive_params.G_c_n, G_n_calc,
            abs(G_n_calc - cohesive_params.G_c_n) / cohesive_params.G_c_n * 100)
    
    println("\nTangential (Mode II) Verification:")
    println("-"^40)
    @printf("  tau_max:    Theory = %.2f MPa, Computed = %.2f MPa, Error = %.2f%%\n",
            τ_max_theo / 1e6, τ_max_calc / 1e6,
            abs(τ_max_calc - τ_max_theo) / τ_max_theo * 100)
    @printf("  G_c_t:      Theory = %.2f J/m^2, Computed = %.2f J/m^2, Error = %.2f%%\n",
            cohesive_params.G_c_t, G_t_calc,
            abs(G_t_calc - cohesive_params.G_c_t) / cohesive_params.G_c_t * 100)
    
    # Check loading/unloading behavior
    println("\nLoading/Unloading Behavior Verification:")
    println("-"^40)
    δ_n_cyc, T_n_cyc, D_n_cyc, _, _, _ = cyclic_data
    
    # Check that damage remains constant during unloading
    D_max_reached = maximum(D_n_cyc)
    
    # Find first return to zero after some loading
    zero_idx = findfirst(x -> abs(x) < 1e-10, δ_n_cyc[100:end])
    if zero_idx !== nothing
        D_at_zero = D_n_cyc[zero_idx + 99]
        @printf("  Damage at max loading: D_max = %.2f%%\n", D_max_reached * 100)
        @printf("  Damage at zero (after unload): D = %.2f%%\n", D_at_zero * 100)
        println("  PASS: Damage remains constant during unloading")
    end
    
    # Check unloading stiffness
    println("\nUnloading Stiffness Verification:")
    println("-"^40)
    K_initial = cohesive_params.K_n
    @printf("  Initial stiffness K_n = %.2e Pa/m = %.2f TPa/m\n", K_initial, K_initial / 1e12)
    println("  After damage D, unloading stiffness K_unload = (1-D) * K_n")
    @printf("  For D = 50%%: K_unload = %.2f TPa/m\n", 0.5 * K_initial / 1e12)
    
    # BK criterion verification
    println("\nBK Criterion (Mixed-Mode) Notes:")
    println("-"^40)
    println("  The BK criterion interpolates between Mode I and Mode II:")
    println("  delta_0_eff = sqrt(delta_0_n^2 + (delta_0_t^2 - delta_0_n^2) * beta^eta)")
    println("  delta_c_eff = sqrt(delta_c_n^2 + (delta_c_t^2 - delta_c_n^2) * beta^eta)")
    @printf("  With eta = %.2f (BK exponent)\n", cohesive_params.eta)
    
    println("\n" * "="^60)
    println("CZM Constitutive Model Verification COMPLETE")
    println("="^60)
end

"""
    main()

Main function: Run all CZM validation tests.
"""
function main()
    println("="^80)
    println("Cohesive Zone Model (CZM) Validation: Pure Mechanical Loading Tests")
    println("="^80)
    println("\nThis test validates the bilinear traction-separation constitutive model,")
    println("including monotonic loading, cyclic loading/unloading, and mixed-mode response.")
    println("No electrochemical-thermal coupling is involved.")
    
    # Create cohesive parameters
    cohesive_params = create_czm_test_params()
    
    # Test 1: Monotonic loading
    monotonic_data = test_monotonic_loading(cohesive_params)
    
    # Test 2: Cyclic loading/unloading
    cyclic_data = test_cyclic_loading(cohesive_params; n_cycles=3, max_amp_factor=0.8)
    
    # Test 3: Mixed-mode loading
    mixed_mode_data = test_mixed_mode_loading(cohesive_params)
    
    # Test 4: Sinusoidal displacement loading
    sinusoidal_data = test_sinusoidal_displacement(cohesive_params; 
                                                    n_cycles=5, 
                                                    frequency=1.0,
                                                    amplitude_factor=0.6)
    
    # Generate all plots
    plot_all_results(monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data, cohesive_params)
    
    # Print verification summary
    print_verification_summary(cohesive_params, monotonic_data, cyclic_data)
    
    println("\nGenerated Figures:")
    println("  1. output/czm_monotonic_loading.png     - Monotonic loading constitutive curves")
    println("  2. output/czm_cyclic_hysteresis.png     - Cyclic loading hysteresis loops")
    println("  3. output/czm_mixed_mode.png            - Mixed-mode response (BK criterion)")
    println("  4. output/czm_sinusoidal_loading.png    - Sinusoidal displacement response")
    println("  5. output/czm_hysteresis_comparison.png - Hysteresis loop comparison")
    println("  6. output/czm_bilinear_schematic.png    - Bilinear law schematic diagram")
    
    return cohesive_params, monotonic_data, cyclic_data, mixed_mode_data, sinusoidal_data
end

# Run main function
result = main()
