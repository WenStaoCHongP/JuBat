# tools/czm_mesh_probe.jl
#
# D13 网格能力探针（spec 2026-08-20 §3.8，Batch 2''）：
#   nθ ∈ {80,360} × thin_subdiv ∈ {1,3}，四组扫描；
#   临界载荷因子 μ_crit（Cholesky 二分：K_mat − μ·K_G 失去正定处）
#   主模态（μ*−δ 逆幂迭代）与周向阶数 n（内边界模态位移 DFT 主峰）。
# 参考态（双跑，判定以压缩参考态为准）：
#   A 卷绕预应力 σ₀（u=0，Batch 2' 语义）；B 压缩本征应变 ε₀=−0.01（先解至平衡 u*，再取 K_G(u*)）。
# 只读诊断：不改求解路径。输出 → output/czm_mesh_probe/（AGENTS §9.9）。

using LinearAlgebra
using SparseArrays
using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

const OUT = joinpath(@__DIR__, "..", "output", "czm_mesh_probe")

"""BC 缩减（fix_inner 默认）：返回自由 dofs 与缩减矩阵。"""
function reduce_free(K::SparseMatrixCSC, bc_dofs::Vector{Int}, ndof::Int)
    isfree = trues(ndof)
    isfree[bc_dofs] .= false
    free = findall(isfree)
    return K[free, free], free
end

"""Cholesky 二分：K0 − μ·K1 失去正定的最小 μ（≥0）。μ_hi 倍增定位后 22 次二分。"""
function mu_crit_chol(K0, K1)
    pd(μ) = try
        cholesky(Symmetric(K0 - μ * K1); check = true); true
    catch
        false
    end
    pd(0.0) || return NaN          # 无预应力态即不正定（不应发生）
    hi = 1.0
    for _ in 1:60
        pd(hi) && (hi *= 2.0; continue)
        break
    end
    pd(hi) && return Inf           # 全稳定
    lo = 0.0
    for _ in 1:22
        mid = 0.5 * (lo + hi)
        pd(mid) ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

"""μ*−δ 处逆幂迭代主模态；返回模态向量与 Rayleigh 商 μ_mode。"""
function buckling_mode(K0, K1, μs; δ = 1e-6 * (1 + abs(μs)))
    A = Symmetric(K0 - (μs - δ) * K1)
    F = cholesky(A)
    v = ones(size(K0, 1))
    v ./= norm(v)
    for _ in 1:50
        w = F \ v
        v = w / norm(w)
    end
    μ_mode = (v' * K0 * v) / (v' * K1 * v)
    return v, μ_mode
end

"""内边界（最内螺旋）模态径向位移的周向阶数主峰 (n, 幅值比)。"""
function mode_circumferential_order(czm_mesh, v_full)
    sub = czm_mesh.czm_submesh.mesh
    # 螺旋 1（最内）= 前 nθn 个节点（层优先布局）；nθn 由首单元连通性反推
    # 单匝窗口：节点 k 的 θ=(k-1)·dθ、dθ=2π/nθ_per_turn —— 前 nθ_per_turn 个节点恰为第一匝
    nθ_per = sub.element[1, 2] - 1        # 每匝角节点数-? element[1,2]=N+1 ⟹ 每螺旋角节点数=N+1
    nθn = nθ_per - 1                      # 第一匝闭合窗口（去重复起点）
    inner = collect(1:nθn)
    nθn >= 16 || return (-1, NaN)   # nθ=8 时窗口 8 点不足，返回 -1（nθ=80/360 充足）
    order = sort(inner, by = k -> atan(sub.node[k, 2], sub.node[k, 1]))
    u_r = zeros(nθn)
    for (i, k) in enumerate(order)
        x, y = sub.node[k, 1], sub.node[k, 2]
        c, s_ = x / hypot(x, y), y / hypot(x, y)
        u_r[i] = c * v_full[2*k-1] + s_ * v_full[2*k]
    end
    u_r .-= sum(u_r) / nθn
    amp = Float64[]
    for n in 1:min(30, nθn ÷ 4)
        re = sum(u_r[i] * cos(2π * (i - 1) * n / nθn) for i in 1:nθn)
        im = sum(u_r[i] * sin(2π * (i - 1) * n / nθn) for i in 1:nθn)
        push!(amp, hypot(re, im))
    end
    p_ = argmax(amp)
    srt = sort(amp, rev = true)
    ratio = srt[2] > 0 ? srt[1] / srt[2] : Inf
    return (p_, ratio)
end

function run_config(nθ::Int, subdiv::Int, ref::String)
    t0 = time()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ = nθ, czm_enabled = true, gsorder = 2,
                                             thin_subdiv = subdiv)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    pc = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, pc)
    cm = case.czm_mesh
    ndof = 2 * cm.nnode
    ne = size(cm.bulk_element, 1)
    bc_dofs, _ = JuBat.extract_bc_dofs(cm, case.param; cache = cache)

    if ref == "prestress"
        σ0 = JuBat.winding_prestress_field(cm, case.param)
        u_star = zeros(ndof)
        _, K0, K1 = JuBat.assemble_bulk_residual_tangent(cm, u_star, pc;
            geo_nl = true, prestress = σ0, split_KG = true)
    else
        eig = (α_eff = 1.0, β_n = 0.0, β_p = 0.0, dT = fill(-0.01, ne),
               Δsn = zeros(ne), Δsp = zeros(ne))
        r, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
            α_eff = 1.0, β_n = 0.0, β_p = 0.0, dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn,
            Δsoc_p_elem = eig.Δsp, max_iter = 100, tol = 1e-8, iter_method = "basic",
            cache = cache, geo_nl = true, eigenstrain = eig)
        r.converged || (@printf("  [warn] 平衡未收敛 res=%.2e，结果参考性\n", r.residual_norm))
        u_star = r.displacement
        _, K0, K1 = JuBat.assemble_bulk_residual_tangent(cm, u_star, pc;
            geo_nl = true, eigenstrain = eig, split_KG = true)
    end

    # split 只含 bulk；补 cohesive 切线（缺它则层间解耦出现机构，K0 非正定）
    K_coh, _, _, _ = JuBat.assemble_czm_system(cm, u_star, pc;
        damage_states = cm.damage_states, geom_cache = cache.cohesive_geom, ws = cache.ws)
    K0 = K0 + K_coh
    K0r, free = reduce_free(K0, bc_dofs, ndof)
    K1r = K1[free, free]
    μ = mu_crit_chol(K0r, K1r)
    n_mode, ratio = (-1, NaN)
    μ_mode = NaN
    if isfinite(μ) && μ > 0
        v, μ_mode = buckling_mode(K0r, K1r, μ)
        v_full = zeros(ndof)
        v_full[free] .= v
        v_full ./= maximum(abs.(v_full))
        n_mode, ratio = mode_circumferential_order(cm, v_full)
    end
    dt = time() - t0
    @printf("  nθ=%3d subdiv=%d ref=%-10s ndof=%7d μ_crit=%12.4e mode_n=%3d ratio=%7.2f  [%.0fs]\n",
            nθ, subdiv, ref, ndof, μ, n_mode, ratio, dt)
    return (nθ = nθ, subdiv = subdiv, ref = ref, ndof = ndof, μ_crit = μ,
            mode_n = n_mode, amp_ratio = ratio, seconds = dt)
end

function main()
    mkpath(OUT)
    rows = []
    for ref in ("compression", "prestress"), subdiv in (1, 3), nθ in (80, 360)
        push!(rows, run_config(nθ, subdiv, ref))
    end
    open(joinpath(OUT, "probe_results.csv"), "w") do io
        println(io, "n_theta,thin_subdiv,reference,ndof,mu_crit,mode_n,amp_ratio,seconds")
        for r in rows
            @printf(io, "%d,%d,%s,%d,%.6e,%d,%.4f,%.1f\n",
                    r.nθ, r.subdiv, r.ref, r.ndof, r.μ_crit, r.mode_n, r.amp_ratio, r.seconds)
        end
    end
    csv_path = joinpath(OUT, "probe_results.csv")
    println("\nCSV → ", csv_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
