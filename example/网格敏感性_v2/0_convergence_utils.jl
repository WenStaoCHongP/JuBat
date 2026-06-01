"""
    0_convergence_utils.jl

网格收敛性分析共享工具函数（GCI 框架）。
被 Script 2-5 通过 `include("0_convergence_utils.jl")` 引入。

Spec: docs/superpowers/specs/2026-06-01-grid-convergence-gci-design.md §5
"""

using Statistics, LinearAlgebra

# ── 离散范数 ──

"""
    l2_norm(f, f_ref) -> Float64

绝对 L2 范数：sqrt(mean((f - f_ref)^2))
"""
function l2_norm(f, f_ref)
    return sqrt(mean((f .- f_ref).^2))
end

"""
    l2_rel_norm(f, f_ref) -> Float64

归一化 L2 范数：||f - f_ref||_L2 / ||f_ref||_L2
分母 ||f_ref||_L2 = sqrt(mean(f_ref^2))，对非零参考场永远有定义。
"""
function l2_rel_norm(f, f_ref)
    norm_ref = sqrt(mean(f_ref.^2))
    norm_ref == 0 && return NaN
    return l2_norm(f, f_ref) / norm_ref
end

"""
    max_norm(f, f_ref) -> Float64

L∞ 最大范数：max|f - f_ref|
"""
function max_norm(f, f_ref)
    return maximum(abs.(f .- f_ref))
end

# ── GCI 框架 ──

"""
    observed_order(f1, f2, f3, r21, r32) -> Float64

计算观测收敛阶 p，使用 Celik et al. (2008) 的迭代法处理非恒定细化比。
f1: 最细网格标量值, f2: 中间, f3: 最粗
r21 = h2/h1, r32 = h3/h2
"""
function observed_order(f1, f2, f3, r21, r32)
    denom = f2 - f1
    numer = f3 - f2
    if abs(denom) < 1e-15 || abs(numer) < 1e-15
        return NaN
    end
    s = sign(denom / numer)
    if s < 0
        return NaN  # 非单调收敛
    end

    alpha = numer / denom

    # 初始值：恒定比近似
    p = abs(log(abs(alpha)) / log(r21))

    # 固定点迭代（Celik et al. 2008）
    for _ in 1:100
        r21p = r21^p
        r32p = r32^p
        if r21p - 1 < 1e-15 || r32p - 1 < 1e-15
            break
        end
        correction = log((r21p - 1) / (r32p - 1))
        p_new = abs(log(abs(alpha)) + correction) / log(r21)
        if abs(p_new - p) < 1e-6
            return p_new
        end
        p = p_new
        if p > 20 || isnan(p)
            return abs(log(abs(alpha)) / log(r21))
        end
    end
    return p
end

"""
    compute_gci(f_fine, f_coarse, r; p, Fs=1.25) -> Float64

计算 Grid Convergence Index [%]。
f_fine: 细网格标量值, f_coarse: 粗网格标量值
r: 细化比 = h_coarse / h_fine
p: 收敛阶
Fs: 安全系数（3+ 级网格用 1.25）
"""
function compute_gci(f_fine, f_coarse, r; p, Fs=1.25)
    abs(f_fine) < 1e-15 && return NaN
    epsilon = (f_coarse - f_fine) / f_fine
    denom = r^p - 1
    denom <= 0 && return NaN
    return Fs * abs(epsilon) / denom * 100
end

"""
    asymptotic_check(gci_12, gci_23, r21, p) -> Float64

渐近收敛检查。返回 GCI_23 / (r21^p * GCI_12)。
比值接近 1.0 表明解在渐近收敛区间内。
"""
function asymptotic_check(gci_12, gci_23, r21, p)
    gci_12 <= 0 && return NaN
    return gci_23 / (r21^p * gci_12)
end

"""
    effective_h(mesh) -> Float64

计算等效单元尺寸 h = sqrt(A_total / N_elem)。
接收 thermal2D mesh 对象（case.mesh["thermal2D"]）。
通过 Gauss 积分点计算单元面积（Mesh 没有 .area 字段）。
"""
function effective_h(mesh)
    ne = size(mesh.element, 1)
    A_total = 0.0
    for g in eachindex(mesh.gs.weight)
        A_total += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    return sqrt(A_total / ne)
end

# ── 场对齐与插值 ──

"""
    align_to_ref(t_cand, y_cand, t_ref)

将候选解插值到参考解的时间网格上，返回对齐后的 y_cand_aligned。
线性插值，超出范围用端点值填充。
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

"""
    interpolate_to_ref_field(coarse_vals, coarse_x, coarse_y,
                             ref_x, ref_y; k=4)

使用反距离加权 (IDW) 将粗网格场插值到参考节点上。
坐标为笛卡尔 (x, y)，直接从 mesh.node[:,1] 和 mesh.node[:,2] 获取。

注：JuBat Mesh 结构的 node 字段存储的已经是笛卡尔坐标（生成时从
极坐标 r*cos(θ), r*sin(θ) 转换），无需额外转换。

coarse_vals: 粗网格节点值 (n_coarse,)
coarse_x, coarse_y: 粗网格节点笛卡尔坐标
ref_x, ref_y: 参考网格节点笛卡尔坐标
k: 最近邻数量
"""
function interpolate_to_ref_field(coarse_vals, coarse_x, coarse_y,
                                  ref_x, ref_y; k=4)
    n_ref = length(ref_x)
    n_coarse = length(coarse_x)
    result = similar(coarse_vals, n_ref)
    k_eff = min(k, n_coarse)

    for i in 1:n_ref
        dists = sqrt.(
            (coarse_x .- ref_x[i]).^2 .+
            (coarse_y .- ref_y[i]).^2
        )
        idx = partialsortperm(dists, 1:k_eff)
        w = 1.0 ./ (dists[idx].^2 .+ 1e-30)
        result[i] = sum(w .* coarse_vals[idx]) / sum(w)
    end
    return result
end

# ── 辅助函数 ──

"""
    trapz(x, y)

梯形积分。
"""
function trapz(x, y)
    return sum(0.5 .* diff(x) .* (y[2:end] .+ y[1:end-1]))
end

"""
    area_error(x, y, x_ref, y_ref)

归一化曲线面积偏差 [%]。
"""
function area_error(x, y, x_ref, y_ref)
    length(x) < 2 && return NaN
    length(x_ref) < 2 && return NaN
    A  = abs(trapz(x, y))
    Ar = abs(trapz(x_ref, y_ref))
    Ar == 0 && return NaN
    return abs(A - Ar) / Ar * 100
end