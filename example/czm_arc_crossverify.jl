"""
跨库验证：用 JuBat 非几何柱面弧长迭代算法复现 BifurcationKit.jl 的入门基准

案例来源：bifurcationkit/BifurcationKit.jl `test/simple_continuation.jl`
（PALC 伪弧长延拓的标准测试曲线）：
    F(x, p) = p·x + x³/3 + 0.01      (imperfect pitchfork, N=1)
平衡支 λ(x) = −x²/3 − 0.01/x（λ≡p），解析折叠点（极限点）由 ∂F/∂x=0 给出：
    x_f = (3·0.01/2)^(1/3) = 0.015^(1/3),  λ_f = −x_f²

本脚本从 BifurcationKit 同一起点（λ₀=−1.5 支上平衡点）出发，用 JuBat 的两条
校正路径追踪该曲线并跨越折叠点：
- 分裂路径：`cylindrical_arc_correction`（K⁻¹R / K⁻¹F 两解 + 二次根选择）
- 增广路径：`cylindrical_arc_bordered_correction`（Chan–Saad bordered，含约束行平衡）

验收：
1) 全程 |λ − f(x)| < 1e-10 且弧长约束 |Δu|²−r² < 1e-10
2) 追踪到的折叠点 λ_max 与解析值 −0.015^(2/3) 相对误差 < 1%（弧长分辨率内）
3) 折叠后支折返（dλ<0）——载荷控制不可达、弧长可达的区域
4) 两条校正路径给出同一条曲线

输出按 AGENTS.md §9.9 写入 `output/czm_arc_crossverify/`。
日期：2026-09-26
"""

using Printf
using LinearAlgebra
using SparseArrays
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

const OUTDIR = joinpath(@__DIR__, "..", "output", "czm_arc_crossverify")
mkpath(OUTDIR)

f_curve(x) = -x^2 / 3 - 0.01 / x          # 平衡支 λ = f(x)
df_curve(x) = -2x / 3 + 0.01 / x^2        # f'(x)（标量"切线"，同既有标量测试约定）
x_fold = 0.015^(1 / 3)
lambda_fold = -x_fold^2

function newton_start(lambda0::Float64)
    x = 0.01 / abs(lambda0)               # 一次近似（x³/3 项小）
    for _ in 1:60
        r = lambda0 - f_curve(x)
        abs(r) < 1e-14 && break
        x += r / df_curve(x)
    end
    return x
end

function trace_curve(::Val{:split}; radius=0.02, nsteps=60)
    x = newton_start(-1.5)
    lambda = -1.5
    alpha = 1.0
    previous_u = nothing
    previous_lambda = 0.0
    xs = [x]; lambdas = [lambda]
    max_res = 0.0; max_arc = 0.0
    for _ in 1:nsteps
        K = df_curve(x)
        tangent = [1.0 / K]
        dl = radius / norm(tangent)
        if previous_u !== nothing &&
           dot(tangent, previous_u) + alpha^2 * previous_lambda < 0.0
            dl = -dl
        end
        x0, lambda0 = x, lambda
        x += tangent[1] * dl
        lambda += dl
        direction_u = previous_u === nothing ? [x - x0] : previous_u
        direction_lambda = previous_u === nothing ? dl : previous_lambda
        for _ in 1:30
            R = lambda - f_curve(x)
            g = (x - x0)^2 - radius^2
            abs(R) < 1e-11 && abs(g) < 1e-11 && break
            K = df_curve(x)
            corr = JuBat.cylindrical_arc_correction([x - x0], [R / K], [1.0 / K],
                lambda, lambda0, radius^2, direction_u, direction_lambda, alpha)
            corr === nothing && error("split 校正失败")
            x += corr.delta_u[1]
            lambda += corr.delta_lambda
            max_res = max(max_res, abs(lambda - f_curve(x)))
            max_arc = max(max_arc, abs((x - x0)^2 - radius^2))
        end
        previous_u = [x - x0]
        previous_lambda = lambda - lambda0
        push!(xs, x); push!(lambdas, lambda)
    end
    return (xs=xs, lambdas=lambdas, max_res=max_res, max_arc=max_arc)
end

function trace_curve(::Val{:bordered}; radius=0.02, nsteps=60)
    x = newton_start(-1.5)
    lambda = -1.5
    previous_u = nothing
    previous_lambda = 0.0
    xs = [x]; lambdas = [lambda]
    max_res = 0.0; max_arc = 0.0
    for _ in 1:nsteps
        K = df_curve(x)
        tangent = [1.0 / K]
        dl = radius / norm(tangent)
        if previous_u !== nothing &&
           dot(tangent, previous_u) + previous_lambda < 0.0
            dl = -dl
        end
        x0, lambda0 = x, lambda
        x += tangent[1] * dl
        lambda += dl
        direction_u = previous_u === nothing ? [x - x0] : previous_u
        direction_lambda = previous_u === nothing ? dl : previous_lambda
        for _ in 1:30
            R = lambda - f_curve(x)
            g = (x - x0)^2 - radius^2
            abs(R) < 1e-11 && abs(g) < 1e-11 && break
            K_bc = spzeros(1, 1)
            K_bc[1, 1] = df_curve(x)
            corr = JuBat.cylindrical_arc_bordered_correction(K_bc, [R], [1.0],
                [x - x0], radius^2)
            corr === nothing && error("bordered 校正失败")
            cand_u = (x - x0) + corr.delta_u[1]
            cand_l = (lambda - lambda0) + corr.delta_lambda
            dot([cand_u], direction_u) + cand_l * direction_lambda >= 0.0 ||
                error("bordered 方向否决")
            x += corr.delta_u[1]
            lambda += corr.delta_lambda
            max_res = max(max_res, abs(lambda - f_curve(x)))
            max_arc = max(max_arc, abs((x - x0)^2 - radius^2))
        end
        previous_u = [x - x0]
        previous_lambda = lambda - lambda0
        push!(xs, x); push!(lambdas, lambda)
    end
    return (xs=xs, lambdas=lambdas, max_res=max_res, max_arc=max_arc)
end

function main()
    println("="^78)
    println("跨库验证：JuBat 柱面弧长 vs BifurcationKit simple_continuation 基准")
    println("F(x,p) = p·x + x³/3 + 0.01；解析折叠点 x_f=$(x_fold), λ_f=$(lambda_fold)")
    println("="^78)

    results = Dict{String, NamedTuple}()
    for (name, fn) in (("split", Val(:split)),
                       ("bordered", Val(:bordered)))
        r = trace_curve(fn)
        results[name] = r
        lam_max, i_max = findmax(r.lambdas)
        rel_err = abs(lam_max - lambda_fold) / abs(lambda_fold)
        beyond_fold = count(diff(r.lambdas) .< 0)
        @printf("\n[%s] 步数=%d  折叠后折返步数=%d\n", name, length(r.lambdas) - 1, beyond_fold)
        @printf("  追踪 λ_max=%.9f (x=%.9f) vs 解析 λ_f=%.9f → 相对误差 %.3e\n",
            lam_max, r.xs[i_max], lambda_fold, rel_err)
        @printf("  全程 |λ−f(x)|最大=%.3e  |Δu|²−r²最大=%.3e\n", r.max_res, r.max_arc)
        @printf("  支端点: (x,λ)=(%.6f, %.6f)\n", r.xs[end], r.lambdas[end])
    end

    n = min(length(results["split"].xs), length(results["bordered"].xs))
    dev = maximum(abs(results["split"].lambdas[i] - results["bordered"].lambdas[i]) for i in 1:n)
    @printf("\n两校正路径曲线最大偏差（前 %d 点）: %.3e\n", n, dev)

    open(joinpath(OUTDIR, "branch_curve.csv"), "w") do io
        println(io, "i,x_split,lambda_split,x_bordered,lambda_bordered,lambda_analytic")
        for i in 1:n
            xs_, ls_ = results["split"].xs[i], results["split"].lambdas[i]
            xb_, lb_ = results["bordered"].xs[i], results["bordered"].lambdas[i]
            println(io, "$i,$xs_,$ls_,$xb_,$lb_,$(f_curve(xs_))")
        end
    end
    println("\n曲线数据已保存: $(joinpath(OUTDIR, "branch_curve.csv"))")
    println("="^78)
end

main()
