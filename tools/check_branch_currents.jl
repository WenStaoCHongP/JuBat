# tools/check_branch_currents.jl
# 验证文件：对 `solve_branch_currents_newton` 的行为做两步（两个时间点）验证。
# 用途：打印 x_e, Vc, status, iters，并在失败时落回到合成模型以验证算法流。

using Printf
using Serialization

# 尝试加载库与函数（确保先加载定义 Case 的文件）
# 预先定义 have_real 以避免未定义错误
have_real = false
try
    include("../src/SetCase.jl")
    include("../src/Initialisation.jl")
    include("../src/SPMe.jl")
    @printf("Included project sources for real-case test.\n")
    have_real = true
catch e
    @printf("Failed to include project sources for real-case test (will use synthetic fallback): %s\n", sprint(showerror, e))
    have_real = false
end

# Wrapper to call solve_branch_currents_newton if available
function run_real_test()
    # Attempt to construct a minimal case using Initialisation / SetCase
    try
        # Try common entrypoints - adapt depending on project API
        # Many repo examples build a Case via SetCase or Initialisation; attempt a minimal path
        case = nothing
        if isdefined(Main, :SetCase) && hasmethod(Main.SetCase, :SetCase)
            # If SetCase provides a constructor, try default
            case = Main.SetCase()
        end
        if case === nothing && isdefined(Main, :Initialisation) && hasmethod(Main.Initialisation, :Initialisation)
            case = Main.Initialisation()
        end
        if case === nothing
            error("Cannot auto-create Case; falling back to synthetic model.")
        end

        # Create yt (state vector) from case.index size if available
        if hasproperty(case, :index) && isa(case.index, Dict)
            nvars = length(keys(case.index))
            yt = zeros(Float64, nvars)
        else
            error("Case lacks index information; fallback.")
        end

        # Prepare areas and temperatures (Te_prev)
        if haskey(case.mesh, "thermal2D")
            mesh = case.mesh["thermal2D"]
            ne = size(mesh.element, 1)
            areas = ones(Float64, ne) .* (1.0 / ne)
            Te_prev = ones(Float64, ne) .* case.param.cell.T0
        else
            # Fallback shapes
            ne = 10
            areas = ones(Float64, ne)
            Te_prev = fill(case.param.cell.T0, ne)
        end

        # choose two total currents to test
        I_tests = [0.1, -0.2]
        for (i, I_total) in enumerate(I_tests)
            t = 0.01 * i
            variables = Dict{String,Union{Array{Float64},Float64}}()
            # call solver - wrap in try to capture errors
            @printf("Running real-case test #%d: I_total=%.6f, t=%.6f\n", i, I_total, t)
            try
                vars_out, x, Vc = solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev)
                @printf("Result #%d: Vc=%.8g, sum(x)=%.8g, ne=%d\n", i, Vc, sum(x), length(x))
                @printf("  status=%.1f, iters=%.1f, converged=%.1f\n",
                        get(vars_out, "thermal2D Vsolve status", -1.0),
                        get(vars_out, "thermal2D Vsolve iters", -1.0),
                        get(vars_out, "thermal2D Vsolve converged", -1.0))
                # print sample of x
                for j in 1:min(8, length(x))
                    @printf("  x[%d]=%.6g\n", j, x[j])
                end
            catch e
                @printf("  Real-case solver threw: %s\n", sprint(showerror, e))
            end
        end
        return true
    catch e
        @printf("Real-case test failed: %s\n", sprint(showerror, e))
        return false
    end
end

# Synthetic model for algorithm verification
function synthetic_V_of_x(e, x_e, params)
    # simple non-linear per-branch model: V = Voc_e - R_e * tanh(k_e * x_e)
    Voc = params[:Voc][e]
    R = params[:R][e]
    k = params[:k][e]
    return Voc - R * tanh(k * x_e)
end

function run_synthetic_test()
    @printf("Running synthetic model test (algorithm verification)\n")
    ne = 12
    areas = ones(Float64, ne)
    A_global = sum(areas)
    # build synthetic parameters (per-branch different OCV/R)
    params = Dict(:Voc => rand(3.6:0.01:3.9, ne), :R => rand(0.01:0.01:0.05, ne), :k => fill(1.0, ne))

    # Implement local Ve_of_x using synthetic model
    function Ve_of_x_local(e::Int, x_e::Float64)
        return synthetic_V_of_x(e, x_e * (A_global / areas[e]), params)
    end

    # Implement solve_x_for_V local copy (same algorithm as in SPMe.jl)
    function solve_x_for_V_local(e::Int, V::Float64; x_hint::Union{Nothing,Float64}=nothing, I_total=1.0)
        bound = max(2.0*abs(I_total), 10.0)
        lo = -bound; hi = bound
        if x_hint !== nothing
            lo = min(lo, x_hint - 2.0*abs(I_total))
            hi = max(hi, x_hint + 2.0*abs(I_total))
        end
        r_lo = Ve_of_x_local(e, lo) - V
        r_hi = Ve_of_x_local(e, hi) - V
        expand = 0
        while r_lo * r_hi > 0.0 && expand < 8
            lo *= 2.0; hi *= 2.0
            r_lo = Ve_of_x_local(e, lo) - V
            r_hi = Ve_of_x_local(e, hi) - V
            expand += 1
        end
        if r_lo * r_hi > 0.0
            return I_total * (areas[e] / A_global)
        end
        xL, xU = lo, hi
        rL, rU = r_lo, r_hi
        for _ in 1:40
            xM = 0.5*(xL + xU)
            rM = Ve_of_x_local(e, xM) - V
            if abs(rM) <= 1e-8
                return xM
            end
            if rL * rM <= 0.0
                xU = xM; rU = rM
            else
                xL = xM; rL = rM
            end
            if abs(xU - xL) <= 1e-10 + 1e-8*max(abs(xU),abs(xL))
                return 0.5*(xU + xL)
            end
        end
        return 0.5*(xL + xU)
    end

    function sum_x_of_V_local(V::Float64, x_hint_vec::Union{Nothing,Vector{Float64}}, I_total)
        xvec = zeros(Float64, ne)
        for e in 1:ne
            x_hint = x_hint_vec === nothing ? nothing : x_hint_vec[e]
            xvec[e] = solve_x_for_V_local(e, V; x_hint=x_hint, I_total=I_total)
        end
        return sum(xvec), xvec
    end

    # Two time steps with different I_total
    I_tests = [0.5, 1.2]
    x_prev = nothing
    for (i, I_total) in enumerate(I_tests)
        # seed x with area-weighted
        x_seed = (areas ./ A_global) .* I_total
        V0_vals = [Ve_of_x_local(e, x_seed[e]) for e in 1:ne]
        V0 = sum(V0_vals) / ne
        dV = 0.05
        f0, x0v = sum_x_of_V_local(V0, x_prev, I_total)
        V1 = V0 + dV
        f1, x1v = sum_x_of_V_local(V1, x0v, I_total)
        tries = 0
        while f0 * f1 > 0.0 && tries < 10
            dV *= 2.0
            V1 = V0 + ((tries % 2 == 0) ? dV : -dV)
            f1, x1v = sum_x_of_V_local(V1, x0v, I_total)
            tries += 1
        end
        if f0 * f1 > 0.0
            # Could not bracket: use V0-based approximate x0v (normalize to I_total)
            x = x0v
            Vc = V0
            sx = sum(x)
            if sx != 0.0
                x .*= (I_total / sx)
            end
            status = "used-V0-approx"
            @printf("Synthetic step %d: could not bracket -> used V0-based approx (normalized)\n", i)
        else
            VL, VR = V0, V1
            fL, fR = f0, f1
            xL = x0v
            outer_iters = 0
            converged = false
            while outer_iters < 30
                VM = 0.5*(VL + VR)
                fM, xM = sum_x_of_V_local(VM, nothing, I_total)
                if abs(fM) <= max(1e-8, 1e-6*abs(I_total))
                    converged = true
                    VL = VM; fL = fM; xL = xM
                    break
                end
                if fL * fM <= 0.0
                    VR = VM; fR = fM
                else
                    VL = VM; fL = fM; xL = xM
                end
                outer_iters += 1
                if abs(VR - VL) <= 1e-10 + 1e-8*max(abs(VR),abs(VL))
                    break
                end
            end
            Vc = 0.5*(VL + VR)
            x = xL
            sx = sum(x)
            if sx != 0.0
                x .*= (I_total / sx)
            end
            status = converged ? "converged" : "bracket-limit"
        end
        @printf("Synthetic step %d: I_total=%.6g, Vc=%.8g, sum(x)=%.8g, status=%s\n", i, I_total, Vc, sum(x), status)
        @printf("  sample x[1..6]: %s\n", join(round.(x[1:min(end,6)], digits=6), ", "))
        # set x_prev for next step
        x_prev = x
    end
end

# Main driver
function main()
    ok = false
    if have_real
        try
            ok = run_real_test()
        catch e
            @printf("Real test raised: %s\n", sprint(showerror, e))
            ok = false
        end
    end
    if !ok
        run_synthetic_test()
    end
    @printf("check_branch_currents.jl completed.\n")
end

main()
