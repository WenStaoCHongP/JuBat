using LinearAlgebra, SparseArrays, Statistics, Plots
include("../src/JuBat.jl")

# 2D distributed thermal (pure thermal validation with custom volumetric heat source)
function main()
    # 1) Options and params (electrochem meshes are created but not used)
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    opt = JuBat.Option()
    Crates = 1
    i =  5*Crates
    opt.Current = x-> i 
    opt.model = "SPMe" 
    opt.Nn = 20; opt.Ns = 10; opt.Np = 20
    opt.Nrn = 15; opt.Nrp = 15
    opt.gsorder = 2
    opt.dimension = 1
    opt.time = [0.0 50]    # 总无量纲时间窗口（相对 t0）；仅用于热时间推进
    opt.dt = [0.5, 0.5]       # 固定无量纲时间步（相对 t0）
    opt.jacobi = "update"
    opt.solveType = "backward" 
    # Thermal
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    # Default true so the example runs with SPMe -> thermal coupling enabled
    electrochem_enabled = false
    case = JuBat.SetCase(param_dim, opt)
    y0 = JuBat.ModelInitialisation(case)
    yt = copy(y0)

    # 2) Thermal mesh (Q4 annulus)
    # Use collector-seeded mesh for thermal discretization
    # jellyroll_collector_seed_mesh expects nθ (number of angular bands) and gsorder
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    # 预计算元素中心与半径用于径向热源与后续可视化
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    # Precompute layer weights for thermal assembly (fallback to sampling if not collector-seeded)
    fks = try
        w = JuBat.jellyroll_get_layer_weights(mesh_th)
        w === nothing ? JuBat.jellyroll_element_layer_weights(mesh_th, param_dim; nsamples_per_dim=4, logic=:spiral) : w
    catch
        nothing
    end

    # 3) 初始变量容器（纯热）
    variables = Dict{String, Union{Array{Float64},Float64}}()
    T = fill(case.param.cell.T0, mesh_th.nlen)
    # 4) Time loop (Backward Euler for both electrochem and thermal)
    t = opt.time[1]
    Tend = opt.time[end]
    dt = opt.dt[1]
    times = Float64[]
    Tmean_hist = Float64[]
    # 热应力时间序列
    stress_mean_hist = Float64[]  # [Pa]
    # 自定义平均体积热源设置
    q_user = try parse(Float64, get(ENV, "SPME_Q", "1e-9")) catch; 1e-9 end
    q_units = get(ENV, "SPME_Q_UNITS", "SI")  # "SI" | "nd"
    istep = 0
    while t <= Tend + 1e-12
        istep += 1
    # 4.1（纯热）可选将均温写回 variables 以便热应力等后处理
    T_mean = mean(T)
    variables["temperature"] = T_mean

        # 4.3 写入热场与层权重
        variables["T_nodes"] = T
        if fks !== nothing
            variables["thermal2D layer_weights"] = fks
        end
            ne = size(mesh_th.element, 1)
            # If electrochem coupling is enabled, compute per-element heat sources via SPMe
            if electrochem_enabled
                # compute element areas and element-averaged temperatures (Te)
                ngs = length(mesh_th.gs.detJ)
                areas = zeros(Float64, ne)
                @inbounds for g in 1:ngs
                    e = mesh_th.gs.ele[g]
                    areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
                end
                A_global = sum(areas)
                # element-averaged temperature (nodal T -> element average)
                T_e = zeros(Float64, ne)
                @inbounds for e in 1:ne
                    nds = mesh_th.element[e, :]
                    T_e[e] = sum(T[nds]) / length(nds)
                end

                # total applied (nondim) current used by SPMe_variables: follow SPMe_variables convention
                param = case.param
                I_total = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ

                # Diagnostics before branch-current solve
                println("[diag] I_total (nd) = ", I_total)
                println("[diag] areas: ne=", ne, ", sum(areas)=", A_global, ", min=", minimum(areas), ", max=", maximum(areas))
                println("[diag] T_e stats (K): mean=", mean(T_e) * case.param_dim.scale.T_ref, ", min=", minimum(T_e) * case.param_dim.scale.T_ref, ", max=", maximum(T_e) * case.param_dim.scale.T_ref)

                # 先建立 SPMe 基线变量（按单元平均温度的全域均值），供热源公式所需的过电位/浓度等量
                I_seed = I_total./ne
                try
                    spme_baseline = JuBat.SPMe_variables(case, yt, t; I_app=I_seed, T_e=mean(T_e))
                    for (k,v) in pairs(spme_baseline)
                        variables[k] = v
                    end
                catch e
                    @warn "SPMe_variables failed in example" exception=(e, catch_backtrace())
                end

                # 诊断：计算并打印 SPMe 基线电压（nd 与 V）、允许区间
                try
                    phi = case.param.scale.phi
                    Vmin_nd = hasproperty(case.param_dim.cell, :v_l) ? case.param_dim.cell.v_l / phi : -Inf
                    Vmax_nd = hasproperty(case.param_dim.cell, :v_h) ? case.param_dim.cell.v_h / phi : Inf
                    if haskey(variables, "cell voltage")
                        V0_nd = variables["cell voltage"]
                        V0_V = V0_nd * phi
                        println("[diag] SPMe baseline V_cell: ", V0_nd, " (nd) = ", V0_V, " [V]; allowed [", Vmin_nd, ", ", Vmax_nd, "] nd -> [", Vmin_nd*phi, ", ", Vmax_nd*phi, "] V")
                    end
                    # 也打印零电流 OCV 供参考
                    spme_ocv = JuBat.SPMe_variables(case, yt, t; I_app=0.0, T_e=mean(T_e))
                    if haskey(spme_ocv, "cell voltage")
                        V_ocv_nd = spme_ocv["cell voltage"]
                        V_ocv_V = V_ocv_nd * phi
                        println("[diag] OCV @ I=0: ", V_ocv_nd, " (nd) = ", V_ocv_V, " [V]")
                    end
                catch e
                    @warn "Voltage diagnostics failed" exception=(e, catch_backtrace())
                end

                # call branch-current solver to get per-element currents I (nondim)
                variables, I_e, Vc = JuBat.solve_branch_currents_newton(case, variables, yt, t, I_total, areas, T_e)
                # ensure consistent naming: store per-element current under expected key
                variables["thermal2D element current"] = I_e
                # heatQ_Source 需要无量纲的分电流
                variables["cell current"] = I_total

                # Diagnostics after branch-current solve (use area fractions w)
                try
                    A_sum = sum(areas)
                    w = A_sum > 0 ? areas ./ A_sum : fill(1.0 / max(length(areas),1), length(areas))
                    sum_wIe = (I_e isa AbstractArray) ? sum(w .* I_e) : Float64(I_e) * sum(w)
                    sum_Ie_phys = case.param_dim.cell.I1C * sum_wIe
                    println("[diag] branch currents: sum(I_e)=", (I_e isa AbstractArray ? sum(I_e) : I_e), ", sum(w .* I_e) (nd)=", sum_wIe, ", sum(Ie_phys) [A]=", sum_Ie_phys, ", Vc=", Vc)
                catch
                    println("[diag] branch currents: failed to compute weighted sums; Vc=", Vc)
                end
                if haskey(variables, "thermal2D Vsolve status")
                    println("[diag] thermal2D Vsolve status: ", variables["thermal2D Vsolve status"])
                end

                # Compute heat sources using the centralized implementation in ThermalDistributed
                # It expects `variables["thermal2D element current"]` to be set (done above) and will
                # populate `variables["heat_source_fields"]` and the units code appropriately.
                variables = JuBat.heatQ_Source(case, variables, t, yt)

                # Diagnostics for heat sources
                if haskey(variables, "heat_source_fields")
                    hs = variables["heat_source_fields"]
                    println("[diag] heat_source_fields: len=", length(hs), ", mean=", mean(hs), ", min=", minimum(hs), ", max=", maximum(hs), ", units_code=", get(variables, "heat_source_units_code", "?") )
                else
                    println("[diag] heat_source_fields not set by heatQ_Source")
                end

            else
                # 组装径向线性递减热源：q(r) = q_user * (1 - (r - Rin)/(Rout - Rin))，外边界衰减至 0
                # 使用元素中心半径 r_centers 计算每个元素的热源强度
                ne = size(mesh_th.element, 1)
                den = max(Rout - Rin, 1e-12)
                profile = clamp.(1 .- (r_centers .- Rin) ./ den, 0.0, 1.0)
                if uppercase(q_units) == "SI"
                    variables["heat_source_fields"] = q_user .* profile
                    variables["heat_source_units_code"] = 1.0
                else
                    q_ref = case.param_dim.scale.q_th
                    variables["heat_source_fields"] = (q_user / max(q_ref, 1e-16)) .* profile
                    variables["heat_source_units_code"] = 0.0
                end
            end

    # 4.4 Assemble thermal system for this step and apply BC
    # 强制各向同性热导，便于得到环状分布
    # ensure Float64 to match variables Dict type Union
    variables["force_isotropic_k"] = 1.0
        MT, KT, FT = JuBat.ThermalDistributed2D(case, variables)
        JuBat.ThermalDistributed2D_BC(KT, FT, case, t)
        # Thermal backward Euler with proper scaling: dt_th = dt * (t0/t_th)
        scale = case.param_dim.scale
        dt_th = dt * (scale.t0 / max(scale.t_th, 1e-16))
        A = (1.0/dt_th) .* MT + KT
        rhsT = (1.0/dt_th) .* (MT * T) + FT
        # small diagonal regularization anchored to ambient (dimensionless)
        alpha = 1e-12
        if alpha > 0
            nT = size(A,1)
            @inbounds for i in 1:nT
                A[i,i] += alpha
            end
            T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref
            rhsT .+= alpha .* T_amb_nd
        end
        T = A \ rhsT

        # Record
        push!(times, t)
        push!(Tmean_hist, mean(T) * case.param_dim.scale.T_ref)  # convert to K

        # Update variables with latest T for thermal stress post-processing and record mean stress
        variables["T_nodes"] = T
        try
            variables = JuBat.thermal_stress(case, variables)
            if haskey(variables, "thermal2D element thermal stress")
                σe = variables["thermal2D element thermal stress"]
                push!(stress_mean_hist, mean(σe))
            else
                # fallback: push NaN to keep vector sizes aligned
                push!(stress_mean_hist, NaN)
            end
        catch e
            @warn "thermal_stress evaluation failed" exception=(e, catch_backtrace())
            push!(stress_mean_hist, NaN)
        end

        t += dt
    end

    # 5) Plots
    p1 = plot(times, Tmean_hist, xlabel="time (nd)", ylabel="mean T [K]", label="mean(T)")
    savefig(p1, "spme_thermal2d_Tmean.png")

    # Stress time history (plot in MPa for readability)
    p2 = plot(times, (stress_mean_hist ./ 1e6);
        xlabel="time (nd)", ylabel="mean thermal stress [MPa]", label="mean(σ_th)")
    savefig(p2, "spme_thermal2d_StressMean.png")

    # Plot element-averaged temperature at element centers to better visualise
    # collector-seeded band mesh (node layout is spokes/in-out curves)
    ne = size(mesh_th.element, 1)
    Te = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        Te[e] = mean(T[nds])
    end
    Te_K = Te .* case.param_dim.scale.T_ref
    # Interpolate nodal temperature to a regular grid using Gaussian kernel (RBF-like) for smooth visualization
    # We'll use node positions (not element centers) as the known field points
    xnod = mesh_th.node[:,1]
    ynod = mesh_th.node[:,2]
    Tnod_K = T .* case.param_dim.scale.T_ref

    # grid resolution (higher for finer image; tune if slow)
    nx = 800; ny = 800
    xs = range(minimum(xnod), stop=maximum(xnod), length=nx)
    ys = range(minimum(ynod), stop=maximum(ynod), length=ny)

    # choose Gaussian kernel width relative to grid spacing; smaller sigma -> preserve more detail
    dx = step(xs); dy = step(ys)
    sigma = 1.0 * max(dx, dy)   # around 1 grid-cell smoothing radius to enhance inner detail
    two_sigma2 = 2.0 * sigma^2

    Z = fill(NaN, ny, nx)

    # mask: only evaluate inside annulus [Rin, Rout]
    # Rin/Rout 已预取

    # pre-copy arrays for speed
    xn = xnod; yn = ynod; tn = Tnod_K
    nn = length(xn)

    @inbounds for j in 1:ny
        yv = ys[j]
        for i in 1:nx
            xv = xs[i]
            r = sqrt(xv^2 + yv^2)
            if r < Rin || r > Rout
                continue
            end
            # compute squared distances to all nodes
            dxv = xn .- xv
            dyv = yn .- yv
            d2 = dxv .* dxv .+ dyv .* dyv
            w = exp.(-d2 ./ two_sigma2)
            s = sum(w)
            if s == 0.0
                Z[j,i] = NaN
            else
                Z[j,i] = sum(w .* tn) / s
            end
        end
    end

    # compute color limits from valid grid points
    valid = .!isnan.(Z)
    Zvals = Z[valid]
    vmin = minimum(Zvals); vmax = maximum(Zvals)
    # apply a mild gamma stretch to improve contrast (gamma < 1 brightens lower values)
    gamma = 0.8
    Z_plot = copy(Z)
    if vmax > vmin
        @inbounds for j in 1:ny, i in 1:nx
            if !isnan(Z[j,i])
                Z_plot[j,i] = vmin + ((Z[j,i] - vmin)/(vmax - vmin))^gamma * (vmax - vmin)
            end
        end
    end

    # plot heatmap + contour and overlay node markers for reference
    plt = plot(size=(1200,1100), title="Final T field (Gaussian-smoothed nodal -> grid)")
    heatmap!(plt, xs, ys, Z_plot; aspect_ratio=1, color=:inferno, colorbar=true, xlabel="x", ylabel="y", clims=(vmin, vmax))
    contour!(plt, xs, ys, Z; levels=16, linewidth=0.5, linecolor=:black, alpha=0.6)
    # overlay node markers (small and semi-transparent)
    scatter!(plt, xnod, ynod; ms=1.5, color=:black, alpha=0.5, label=false)
    savefig(plt, "spme_thermal2d_Tfield.png")
    # export high-resolution files for print/publishing
    try
        # SVG (vector) - good for arbitrary scaling/printing
        savefig(plt, "spme_thermal2d_Tfield.svg")
        # For high-res PNG, recreate the plot at a large pixel size (approx 6in x 5.5in @300dpi -> 1800x1650)
        plt_hr = plot(size=(1800,1650), title="Final T field (Gaussian-smoothed nodal -> grid)")
        heatmap!(plt_hr, xs, ys, Z_plot; aspect_ratio=1, color=:inferno, colorbar=true, xlabel="x", ylabel="y", clims=(vmin, vmax))
        contour!(plt_hr, xs, ys, Z; levels=16, linewidth=0.5, linecolor=:black, alpha=0.6)
        scatter!(plt_hr, xnod, ynod; ms=1.5, color=:black, alpha=0.5, label=false)
        savefig(plt_hr, "spme_thermal2d_Tfield_300dpi.png")
        println("Saved high-res: spme_thermal2d_Tfield.svg and spme_thermal2d_Tfield_300dpi.png")
    catch e
        @warn "Failed to save high-res images" exception=(e, catch_backtrace())
    end

    println("Saved: spme_thermal2d_Tmean.png")
    println("Saved: spme_thermal2d_Tfield.png")

    # ======================
    # Final thermal stress field (element-based -> grid via RBF-like smoothing)
    # ======================
    σe = haskey(variables, "thermal2D element thermal stress") ? variables["thermal2D element thermal stress"] : zeros(Float64, ne)
    σe_MPa = σe ./ 1e6

    # Use element centers as sample points
    xc = centers[:,1]; yc = centers[:,2]
    # Build grid same as T plot
    Zs = fill(NaN, ny, nx)
    two_sigma2_s = two_sigma2  # reuse same smoothing radius
    @inbounds for j in 1:ny
        yv = ys[j]
        for i in 1:nx
            xv = xs[i]
            r = sqrt(xv^2 + yv^2)
            if r < Rin || r > Rout
                continue
            end
            dxv = xc .- xv
            dyv = yc .- yv
            d2 = dxv .* dxv .+ dyv .* dyv
            w = exp.(-d2 ./ two_sigma2_s)
            s = sum(w)
            if s == 0.0
                Zs[j,i] = NaN
            else
                Zs[j,i] = sum(w .* σe_MPa) / s
            end
        end
    end

    valid_s = .!isnan.(Zs)
    if any(valid_s)
        Zsvals = Zs[valid_s]
        vmin_s = minimum(Zsvals); vmax_s = maximum(Zsvals)
        gamma_s = 0.9
        Zs_plot = copy(Zs)
        if vmax_s > vmin_s
            @inbounds for j in 1:ny, i in 1:nx
                if !isnan(Zs[j,i])
                    Zs_plot[j,i] = vmin_s + ((Zs[j,i] - vmin_s)/(vmax_s - vmin_s))^gamma_s * (vmax_s - vmin_s)
                end
            end
        end

        plt_s = plot(size=(1200,1100), title="Final thermal stress field [MPa]")
        heatmap!(plt_s, xs, ys, Zs_plot; aspect_ratio=1, color=:viridis, colorbar=true, xlabel="x", ylabel="y", clims=(vmin_s, vmax_s))
        contour!(plt_s, xs, ys, Zs; levels=16, linewidth=0.5, linecolor=:black, alpha=0.6)
        # overlay element centers for reference
        scatter!(plt_s, xc, yc; ms=1.5, color=:black, alpha=0.4, label=false)
        savefig(plt_s, "spme_thermal2d_StressField.png")
        try
            savefig(plt_s, "spme_thermal2d_StressField.svg")
            plt_s_hr = plot(size=(1800,1650), title="Final thermal stress field [MPa]")
            heatmap!(plt_s_hr, xs, ys, Zs_plot; aspect_ratio=1, color=:viridis, colorbar=true, xlabel="x", ylabel="y", clims=(vmin_s, vmax_s))
            contour!(plt_s_hr, xs, ys, Zs; levels=16, linewidth=0.5, linecolor=:black, alpha=0.6)
            scatter!(plt_s_hr, xc, yc; ms=1.5, color=:black, alpha=0.4, label=false)
            savefig(plt_s_hr, "spme_thermal2d_StressField_300dpi.png")
            println("Saved high-res: spme_thermal2d_StressField.svg and spme_thermal2d_StressField_300dpi.png")
        catch e
            @warn "Failed to save high-res stress images" exception=(e, catch_backtrace())
        end
        println("Saved: spme_thermal2d_StressField.png")
    else
        @warn "No valid thermal stress values to plot (Zs all NaN)." 
    end
end

main()
