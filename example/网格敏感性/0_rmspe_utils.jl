"""
    0_rmspe_utils.jl

网格敏感性分析共享工具函数。
被 Script 2-5 通过 `include("0_rmspe_utils.jl")` 引入。

Spec: docs/superpowers/specs/2026-04-29-grid-sensitivity-statistical-metrics-design.md §4
"""

using Statistics

"""
    rmspe(y, y_ref; rel_tol=1e-3) -> (rmspe_val, skip_rate)

计算相对均方根百分比误差。跳过 |y_ref| < rel_tol * max(|y_ref|) 的点。
返回 (RMSPE值, 跳过率)。若跳过率 > 50%，调用者应改用绝对误差。
"""
function rmspe(y, y_ref; rel_tol=1e-3)
    threshold = rel_tol * maximum(abs.(y_ref))
    mask = abs.(y_ref) .> threshold
    skip_rate = 1.0 - count(mask) / length(y_ref)
    count(mask) == 0 && return (NaN, 1.0)
    val = sqrt(mean(((y[mask] .- y_ref[mask]) ./ y_ref[mask]).^2)) * 100
    return (val, skip_rate)
end

"""
    spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)

计算空间场 RMSPE 的时间平均值。
T_hist: (nnode × nt) 矩阵
"""
function spatial_rmspe_over_time(T_hist, T_ref_hist; rel_tol=1e-3)
    nt = size(T_hist, 2)
    errs = Float64[]
    for k in 1:nt
        val, skip = rmspe(T_hist[:,k], T_ref_hist[:,k]; rel_tol)
        isnan(val) || push!(errs, val)
    end
    isempty(errs) && return 0.0
    return mean(errs)
end

"""
    area_error(x, y, x_ref, y_ref)

计算归一化曲线面积偏差（梯形积分）。
x 和 y 至少需要 2 个点。
"""
function area_error(x, y, x_ref, y_ref)
    length(x) < 2 && return NaN
    length(x_ref) < 2 && return NaN
    A  = abs(trapz(x, y))
    Ar = abs(trapz(x_ref, y_ref))
    Ar == 0 && return NaN
    return abs(A - Ar) / Ar * 100
end

"""
    trapz(x, y)

简单梯形积分。
"""
function trapz(x, y)
    return sum(0.5 .* diff(x) .* (y[2:end] .+ y[1:end-1]))
end

"""
    align_to_ref(t_cand, y_cand, t_ref)

将候选解插值到参考解的时间网格上，返回对齐后的 y_cand_aligned。
不依赖外部包，手写线性插值。超出 t_cand 范围的值用端点值填充。
"""
function align_to_ref(t_cand, y_cand, t_ref)
    return [let
        idx = searchsortedfirst(t_cand, t)
        if idx == 1
            y_cand[1]
        elseif idx > length(t_cand)
            y_cand[end]
        else
            frac = (t - t_cand[idx-1]) / (t_cand[idx] - t_cand[idx-1])
            y_cand[idx-1] + frac * (y_cand[idx] - y_cand[idx-1])
        end
    end for t in t_ref]
end
