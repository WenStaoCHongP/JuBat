# tools/check_branch_currents_and_heat.jl
# 检查分电流求解、SPMe 基线变量与热源计算在两个时间步是否正确工作
# 用法: 在仓库根目录 (含 src/) 下运行：
# julia --project=. tools/check_branch_currents_and_heat.jl

using LinearAlgebra, SparseArrays, Statistics, Printf
include("../src/JuBat.jl")

function describe_value(name, value; max_entries=6)
    info = String[]
    push!(info, "[diag] $(name) type=$(typeof(value))")
    if isa(value, AbstractArray)
        push!(info, "size=$(size(value))")
        flat = vec(value)
        if !isempty(flat)
            sample = flat[1:min(end, max_entries)]
            push!(info, "first=$(join(string.(sample), ", "))")
        end
    elseif isa(value, Number)
        push!(info, "value=$(value)")
    end
    println(join(info, "; "))
end

function safe_vec(x)
    if isa(x, Number)
        return [Float64(x)]
    elseif isa(x, AbstractVector)
        return Vector{Float64}(x)
    elseif isa(x, AbstractArray)
        # try collapse first column
        try
            return Vector{Float64}(x[:,1])
        catch
            return vec(Float64.(x))
        end
    else
        return Float64[]
    end
end

function check_two_steps()
    println("=== JuBat coupling check: branch currents + SPMe + heat source (2 steps) ===")

    # 1) 简化参数/选项以加快执行
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    # 设置较小网格以便快速运行
    opt.Nn = 8; opt.Ns = 6; opt.Np = 8
    opt.Nrn = 8; opt.Nrp = 8
    opt.gsorder = 2
    opt.model = "SPMe"
    # constant current (dimensional) - choose modest value
    Crates = 1
    iA = 1.0 * Crates
    opt.Current = t -> iA

    # Thermal settings
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"

    # time steps: two steps (nd time)
    opt.time = [0.0, 1.0]
    opt.dt = [0.1]
    opt.czm_enabled = false

    case = JuBat.SetCase(param_dim, opt)

    # initial electrochem state
    y0 = JuBat.ModelInitialisation(case)
    yt = copy(y0)

    # thermal mesh (collector-seeded) - keep moderate nθ for speed
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # precompute fks (允许为 nothing)
    fks = try
        JuBat.get_element_layer_weights(mesh_th, param_dim; nsamples_per_dim=3, logic=:spiral)
    catch err
        @warn "get_element_layer_weights failed" err
        println("[diag] JuBat has get_element_layer_weights? ", isdefined(JuBat, :get_element_layer_weights))
        nothing
    end

    # initial variables and temperature field
    variables = Dict{String, Union{Array{Float64},Float64}}()
    if fks !== nothing
        variables["thermal2D layer_weights"] = fks
    end

    T = fill(case.param.cell.T0, mesh_th.nlen)

    # helper: compute element areas & element avg temperature
    function compute_areas_and_Te(mesh, Tnod)
        ne = size(mesh.element,1)
        ngs = length(mesh.gs.detJ)
        areas = zeros(Float64, ne)
        for g in 1:ngs
            e = mesh.gs.ele[g]
            areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
        end
        Te = zeros(Float64, ne)
        for e in 1:ne
            nds = mesh.element[e,:]
            Te[e] = sum(Tnod[nds]) / length(nds)
        end
        return areas, Te
    end

    areas, T_e = compute_areas_and_Te(mesh_th, T)
    ne = size(mesh_th.element,1)
    println("[diag] mesh_th element count=", ne, ", node count=", size(mesh_th.node,1))
    println("[diag] yt size=", size(yt))

    # run two timesteps: t in opt.time (we'll do first two entries)
    times = opt.time[1:min(2, length(opt.time))]
    for (step_idx, t) in enumerate(times)
        println("\n--- Step $step_idx: t=$t (nd) ---")
        variables["T_nodes"] = T
        variables["thermal2D element area"] = areas

        # total applied nondim current following earlier conventions
        param = case.param
        I_total = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
        println("I_total (nd) = ", I_total)

        # ensure baseline SPMe variables exist (with backtrace on failure)
        spme_ok = false
        try
            spme_baseline = JuBat.SPMe_variables(case, yt, t; I_app=I_total, T_e=mean(T_e))
            for (k,v) in pairs(spme_baseline)
                variables[k] = v
            end
            spme_ok = true
        catch err
            @error "SPMe_variables failed" exception=(err, catch_backtrace())
        end

        # call branch-current solver (with shape/type checks and backtrace on failure)
        bc_ok = false
        I_e = Float64[]; Vc = NaN
        try
            # pre-check: if layer weights exist, ensure shape roughly matches ne x 5
            fks = get(variables, "thermal2D layer_weights", nothing)
            if fks !== nothing
                # print diagnostics about fks for debugging
                try
                    println("[diag] fks typeof=", typeof(fks), ", eltype=", (isa(fks, AbstractArray) ? eltype(fks) : "N/A"), ", ndims=", (isa(fks, AbstractArray) ? ndims(fks) : 0), ", size=", (isa(fks, AbstractArray) ? size(fks) : ()))
                    if isa(fks, AbstractArray) && ndims(fks) >= 2
                        nprint = min(5, size(fks,1))
                        println("[diag] fks sample rows (first $nprint rows):")
                        for i in 1:nprint
                            println("  row $i: ", collect(fks[i, :]))
                        end
                        # print row-sum statistics to confirm each row sums to ~1.0
                        try
                            row_sums = vec(sum(fks, dims=2))
                            println("[diag] fks row sums: min=", minimum(row_sums), ", max=", maximum(row_sums), ", mean=", mean(row_sums), ", std=", std(row_sums))
                            # print first few row sums
                            println("[diag] fks sample row sums (first $nprint): ", row_sums[1:nprint])
                            # check deviation from 1
                            max_dev = maximum(abs.(row_sums .- 1.0))
                            println("[diag] fks max abs deviation from 1 = ", max_dev)
                        catch err
                            @warn "failed to compute fks row sums" err
                        end
                    else
                        println("[diag] fks value: ", fks)
                    end
                catch err
                    @warn "failed to print fks diagnostics" err
                end
                if !(fks === nothing) && isa(fks, AbstractArray)
                    if ndims(fks) == 2
                        if size(fks,1) != ne
                            @error "thermal2D layer_weights has wrong number of rows" expected=ne got=size(fks,1) exception=(ErrorException("layer_weights shape mismatch"), catch_backtrace())
                        end
                    end
                end
            end

            variables, I_e, Vc = JuBat.solve_branch_currents_newton(case, variables, yt, t, I_total, areas, T_e)
            variables["thermal2D element current"] = I_e
            variables["cell current"] = I_total
            bc_ok = true
        catch err
            @error "solve_branch_currents_newton failed" exception=(err, catch_backtrace())
        end

        # diagnostics for branch currents
        try
            sI = (I_e isa AbstractArray) ? sum(I_e) : Float64(I_e)
            # compute area fractions and physical branch currents for diagnostics
            A_sum = sum(areas)
            w = (A_sum > 0.0) ? areas ./ A_sum : fill(1.0 / max(length(areas),1), length(areas))
            sum_wIe = (I_e isa AbstractArray) ? sum(w .* I_e) : Float64(I_e) * sum(w)
            sum_Ie_phys = case.param_dim.cell.I1C * sum_wIe
            println("branch currents: length=", (I_e isa AbstractArray ? length(I_e) : 1), ", sum(I_e) = ", sI,
                    ", sum(w .* I_e) (nd)=", sum_wIe, ", sum(Ie_phys) [A] =", sum_Ie_phys)
            println("Vc = ", Vc)

            # 打印每个网格的电压（在 t=0 时打印详细信息）
            try
                if step_idx == 1
                    println("[diag] per-element voltages at t=0 (nd)")
                    V_elem_nd = fill(Vc, ne)
                    nprint = min(6, length(V_elem_nd))
                    println(" V_elem_nd (first $nprint): ", V_elem_nd[1:nprint])
                    # 若存在物理刻度，打印物理电压 [V]
                    if hasproperty(case, :param) && hasproperty(case.param, :scale)
                        phi_scale = case.param.scale.phi
                        V_elem = V_elem_nd .* phi_scale
                        println(" V_elem [V] (first $nprint): ", V_elem[1:nprint])
                    else
                        # 兜底：尝试使用 case.param.scale.phi（若可用）
                        try
                            phi_scale = case.param.scale.phi
                            V_elem = V_elem_nd .* phi_scale
                            println(" V_elem [V] (first $nprint): ", V_elem[1:nprint])
                        catch
                            println(" [diag] 无可用的 phi 缩放，跳过物理电压输出")
                        end
                    end

                    # 同时打印变量中的端电压（若存在）
                    cv = get(variables, "cell voltage", nothing)
                    if cv === nothing
                        println("[diag] variables['cell voltage'] not present")
                    else
                        try
                            if isa(cv, AbstractArray)
                                ncv = min(6, length(cv))
                                println("[diag] variables[\"cell voltage\"] (nd) first $ncv = ", cv[1:ncv])
                                println("[diag] variables[\"cell voltage\"] [V] first $ncv = ", (cv[1:ncv] .* case.param.scale.phi))
                            else
                                println("[diag] variables[\"cell voltage\"] (nd) = ", cv)
                                println("[diag] variables[\"cell voltage\"] [V] = ", cv * case.param.scale.phi)
                            end
                        catch err
                            @warn "打印 variables[\"cell voltage\"] 失败" err
                        end
                    end
                end
            catch err
                @warn "printing per-element voltages failed" err
            end
        catch err
            @warn "branch current diagnostics failed" err
        end

        # call heat source calculation (with detailed diagnostics and backtrace on failure)
        hs_ok = false
        try
            println("[diag] heatQ_Source inputs summary")
            println("[diag] variables keys=", sort(collect(keys(variables))))
            for key in ["thermal2D layer_weights", "thermal2D element current", "thermal2D element area", "T_nodes", "cell current"]
                val = get(variables, key, nothing)
                if val === nothing
                    println("[diag] missing key ", key)
                else
                    describe_value(key, val)
                end
            end
            describe_value("areas", areas)
            describe_value("T_e", T_e)
            variables = JuBat.heatQ_Source(case, variables, t, yt)
            hs_ok = get(variables, "heat_source_fields", nothing) !== nothing
        catch err
            @error "heatQ_Source failed" exception=(err, catch_backtrace())
            println("[diag] heatQ_Source failure diagnostics: mesh thermal2D element size=", size(case.mesh["thermal2D"].element))
            println("[diag] mesh thermal2D node size=", size(case.mesh["thermal2D"].node))
            println("[diag] yt size during failure=", size(yt))
            println("[diag] available variables after failure: ", sort(collect(keys(variables))))
        end

        # diagnostics for heat sources
        if hs_ok
            hs = variables["heat_source_fields"]
            try
                println("heat_source_fields: len=", length(hs), ", mean=", mean(hs), ", min=", minimum(hs), ", max=", maximum(hs))
                # simple checks
                n_ok = (length(hs) == ne)
                nan_ok = !any(x->(!isfinite(x)), hs)
                println("heat source length matches ne? ", n_ok, ", any NaN/Inf? ", !nan_ok ? "YES" : "NO")

                # If heat source components are present, verify they sum to total
                comp_keys = filter(k->startswith(k, "heat_source_") && k != "heat_source_fields" && k != "heat_source_units_code", collect(keys(variables)))
                if !isempty(comp_keys)
                    println("Found heat source component keys: ", comp_keys)
                    comps = Dict{String, Vector{Float64}}()
                    for k in comp_keys
                        v = variables[k]
                        try
                            comps[k] = safe_vec(v)
                        catch
                            comps[k] = Float64[]
                        end
                    end
                    # check shapes
                    all_ok = true
                    for (k,v) in pairs(comps)
                        if length(v) != ne
                            @warn "component length mismatch" key=k expected_ne=ne got_length=length(v)
                            all_ok = false
                        end
                    end
                    if all_ok
                        # elementwise sum
                        sum_comps = zeros(Float64, ne)
                        for v in values(comps)
                            sum_comps .+= v
                        end
                        diffs = abs.(sum_comps .- hs)
                        maxdiff = maximum(diffs)
                        println("component sum vs total: max abs diff = ", maxdiff)
                        if maxdiff > 1e-8
                            @warn "heat component sums do not match total heat_source_fields" maxdiff=maxdiff
                        end
                    end
                end

            catch err
                @error "heat source diagnostics error" exception=(err, catch_backtrace())
            end
        else
            println("heat_source_fields not set by heatQ_Source")
        end

        # print presence and types of key SPMe variables used by heatQ_Source
        # Check SPMe variables: use the canonical full names (do NOT rename variables in src)
        keys_of_interest = ["negative electrode overpotential", "positive electrode overpotential", "negative particle surface lithium concentration", "positive particle surface lithium concentration", "cell current"]
        println("Variables presence & types (canonical names):")
        for k in keys_of_interest
            v = get(variables, k, nothing)
            if v !== nothing
                tname = typeof(v)
                len = isa(v, AbstractArray) ? length(v) : 1
                println(" - $k : type=", tname, ", len=", len)
            else
                println(" - $k : MISSING")
            end
        end
        # (legacy short-name checks removed)  -- use canonical variable names only

        # simple PASS/FAIL for this step
        step_ok = spme_ok && bc_ok && hs_ok
        println("Step $step_idx result: ", step_ok ? "PASS" : "FAIL")
    end
    println("\n=== Check finished ===")
end

check_two_steps()
