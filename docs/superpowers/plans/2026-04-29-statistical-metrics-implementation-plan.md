# 网格敏感性分析统计指标实施计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 `example/网格敏感性/` 下脚本 2-5 的点值对比指标替换为 RMSPE 统计指标

**Architecture:** 新增一个共享工具函数文件 `0_rmspe_utils.jl`，各脚本 include 该文件后替换后处理逻辑。仿真配置不变，只改汇总表、误差计算和绘图。

**Tech Stack:** Julia, Plots.jl, Statistics (stdlib)

**Spec:** `docs/superpowers/specs/2026-04-29-grid-sensitivity-statistical-metrics-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `example/网格敏感性/0_rmspe_utils.jl` | 共享工具函数：rmspe, spatial_rmspe_over_time, area_error, trapz, align_to_ref |
| Modify | `example/网格敏感性/2_electrochemical_mesh_sensitivity.jl` | 电化学 Track 后处理 |
| Modify | `example/网格敏感性/3_thermal_mesh_sensitivity.jl` | 热学 Track 后处理 |
| Modify | `example/网格敏感性/4_czm_mesh_sensitivity.jl` | CZM Track 后处理 |
| Modify | `example/网格敏感性/5_energy_conservation_check.jl` | 能量守恒后处理 |
| Modify | `docs/planning-with-files/网格敏感性分析/findings.md` | 更新发现记录 |
| Modify | `docs/planning-with-files/网格敏感性分析/progress.md` | 更新进度记录 |

---

## Chunk 1: 共享工具函数

### Task 1: 创建 `0_rmspe_utils.jl`

**Files:**
- Create: `example/网格敏感性/0_rmspe_utils.jl`

- [ ] **Step 1: 创建工具函数文件**

```julia
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
    return [begin
        idx = searchsortedfirst(t_cand, t)
        idx == 1 && (return y_cand[1])
        idx > length(t_cand) && (return y_cand[end])
        frac = (t - t_cand[idx-1]) / (t_cand[idx] - t_cand[idx-1])
        y_cand[idx-1] + frac * (y_cand[idx] - y_cand[idx-1])
    end for t in t_ref]
end
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性/0_rmspe_utils.jl"
git commit -m "feat: add shared RMSPE utility functions for mesh sensitivity scripts"
```

---

## Chunk 2: 电化学 Track (Script 2)

### Task 2: 修改 Script 2 后处理

**Files:**
- Modify: `example/网格敏感性/2_electrochemical_mesh_sensitivity.jl`

- [ ] **Step 1: 替换文件头部 — 添加 include 和 Statistics**

将文件头部 `using Printf, Plots` 修改为：

```julia
using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_rmspe_utils.jl"))
```

- [ ] **Step 2: 替换脚本 docstring**

将 docstring 中的指标描述改为：

```julia
"""
Script 2: 电化学网格敏感性分析

4 组 (Nn, Ns, Np) 下运行 SPMe + lumped thermal：
  - (40, 20, 40)  ← 参考解
  - (20, 10, 20)
  - (20,  5, 20)
  - (10,  5, 10)

指标：V(t) RMSPE、T(t) RMSPE、dT/dt(t) RMSPE
输出：收敛图
"""
```

- [ ] **Step 3: 替换 main() 中的参考解提取和汇总逻辑**

将第 76-121 行（从 `# ── 参考解 ──` 到汇总循环结束）替换为：

```julia
    # ── 参考解 ──
    ref = results[MESH_CONFIGS[1]]
    t_ref = ref["time [s]"]
    V_ref = ref["cell voltage [V]"]
    T_ref = ref["temperature [K]"]

    # ── 汇总表 (RMSPE) ──
    println("\n" * "=" ^ 70)
    println("电化学网格收敛性汇总 (RMSPE)")
    println("=" ^ 70)
    @printf("%-16s  %12s  %12s  %12s  %10s  %10s  %10s  %10s\n",
            "Mesh", "V(t) RMSPE%", "T(t) RMSPE%", "dT/dt RMSPE%",
            "V_skip%", "T_skip%", "dT_skip%", "Wall [s]")
    println("-" ^ 100)

    V_errors = Float64[]
    T_errors = Float64[]
    dTdt_errors = Float64[]

    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        t_c = r["time [s]"]
        V_c = r["cell voltage [V]"]
        T_c = r["temperature [K]"]

        # 对齐到参考时间网格
        V_aligned = align_to_ref(t_c, V_c, t_ref)
        T_aligned = align_to_ref(t_c, T_c, t_ref)

        # 电压曲线 RMSPE
        err_V, skip_V = i == 1 ? (0.0, 0.0) : rmspe(V_aligned, V_ref)

        # 温度曲线 RMSPE
        err_T, skip_T = i == 1 ? (0.0, 0.0) : rmspe(T_aligned, T_ref)

        # dT/dt 曲线 RMSPE
        if i == 1
            err_dTdt, skip_dTdt = (0.0, 0.0)
        else
            dTdt_ref_vals = diff(T_ref) ./ diff(t_ref)
            dTdt_cand_vals = diff(T_aligned) ./ diff(t_ref)
            err_dTdt, skip_dTdt = rmspe(dTdt_cand_vals, dTdt_ref_vals)
        end

        push!(V_errors, err_V)
        push!(T_errors, err_T)
        push!(dTdt_errors, err_dTdt)

        @printf("(%2d,%2d,%2d)       %12.4f  %12.4f  %12.4f  %10.1f  %10.1f  %10.1f  %10.2f\n",
                Nn, Ns, Np, err_V, err_T, err_dTdt,
                skip_V*100, skip_T*100, skip_dTdt*100, wall_times[i])
    end
```

- [ ] **Step 4: 替换图 3 的标签**

将图 3 的绘图部分（原第 153-168 行）替换为：

```julia
    # 图3: 收敛误差折线图 (RMSPE)
    xp = collect(1:4)
    p3 = plot(xlabel="Mesh No.", ylabel="RMSPE [%]",
              title="Electrochemical Mesh: Convergence Error (RMSPE)",
              xticks=xp, xlims=(0.5, 4.5),
              legend=:topright)
    plot!(p3, xp, V_errors, marker=:o, lw=2, label="V(t) RMSPE", color=:blue)
    plot!(p3, xp, T_errors, marker=:s, lw=2, label="T(t) RMSPE", color=:orange)
    plot!(p3, xp, dTdt_errors, marker=:d, lw=2, label="dT/dt RMSPE", color=:green)
    hline!([5.0], label="5% threshold", color=:black, ls=:dash, lw=2)
    hline!([1.0], label="1%", color=:gray, ls=:dot, lw=1)
    mesh_txt = "1→(40,20,40) ref\n2→(20,10,20)\n3→(20,5,20)\n4→(10,5,10)"
    annotate!(p3, 4.4, maximum([V_errors; T_errors; dTdt_errors]) * 0.95,
              text(mesh_txt, 8, :right, :top, :Courier))
    savefig(p3, joinpath(out_dir, "echem_convergence_error.png"))
```

- [ ] **Step 5: Commit**

```bash
git add "example/网格敏感性/2_electrochemical_mesh_sensitivity.jl"
git commit -m "refactor: replace point-value metrics with RMSPE in electrochemical script"
```

---

## Chunk 3: 热学 Track (Script 3)

### Task 3: 修改 Script 3 后处理

**Files:**
- Modify: `example/网格敏感性/3_thermal_mesh_sensitivity.jl`

- [ ] **Step 1: 替换文件头部**

在 `using Printf, Plots, Statistics` 之后添加：

```julia
include(joinpath(@__DIR__, "0_rmspe_utils.jl"))
```

- [ ] **Step 2: 替换 docstring**

```julia
"""
Script 3: 热网格敏感性分析

4 组 nθ = {20, 40, 80, 160} 下运行纯热模型（均匀体积热源 + 表面冷却），
基于 Jellyroll spiral mesh，使用 opt.model = "thermal" 直接求解热方程。

指标：T_max(t) RMSPE、T_range(t) RMSPE、空间场 RMSPE
输出：收敛图
"""
```

- [ ] **Step 3: 替换 main() 中的汇总逻辑**

将第 105-121 行（从 `# ── 参考解` 到汇总循环结束）替换为：

```julia
    # ── 参考解 (nθ=160) ──
    ref_idx = length(THERMAL_Nθ)
    ref_result = all_results[THERMAL_Nθ[ref_idx]]
    T_hist_ref = ref_result.T_hist
    t_ref = ref_result.times
    nt_ref = length(t_ref)

    # 参考解 T_max(t) 和 T_range(t) 曲线
    Tmax_ref_curve = [maximum(T_hist_ref[:, k]) for k in 1:nt_ref]
    Trange_ref_curve = [maximum(T_hist_ref[:, k]) - minimum(T_hist_ref[:, k]) for k in 1:nt_ref]

    println("\n" * "=" ^ 70)
    println("热网格收敛性汇总 (RMSPE)")
    println("=" ^ 70)
    @printf("%-10s  %12s  %12s  %12s  %10s\n",
            "nθ", "Tmax RMSPE%", "Trange RMSPE%", "Spatial RMSPE%", "Solve [s]")
    println("-" ^ 70)

    Tmax_errors = Float64[]
    Trange_errors = Float64[]
    Spatial_errors = Float64[]

    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        t_c = r.times
        nt_c = length(t_c)

        if i == ref_idx
            push!(Tmax_errors, 0.0)
            push!(Trange_errors, 0.0)
            push!(Spatial_errors, 0.0)
            @printf("%-10d  %12s  %12s  %12s  %10.2f\n",
                    nθ, "ref", "ref", "ref", all_wall[i])
            continue
        end

        # T_max(t) 曲线 RMSPE
        Tmax_cand = [maximum(T_hist[:, k]) for k in 1:nt_c]
        Tmax_cand_aligned = align_to_ref(t_c, Tmax_cand, t_ref)
        err_Tmax, _ = rmspe(Tmax_cand_aligned, Tmax_ref_curve)

        # T_range(t) 曲线 RMSPE
        Trange_cand = [maximum(T_hist[:, k]) - minimum(T_hist[:, k]) for k in 1:nt_c]
        Trange_cand_aligned = align_to_ref(t_c, Trange_cand, t_ref)
        err_Trange, _ = rmspe(Trange_cand_aligned, Trange_ref_curve)

        # 空间场 RMSPE（对齐时间步）
        nt_overlap = min(nt_c, nt_ref)
        err_spatial = spatial_rmspe_over_time(
            T_hist[:, 1:nt_overlap], T_hist_ref[:, 1:nt_overlap])

        push!(Tmax_errors, err_Tmax)
        push!(Trange_errors, err_Trange)
        push!(Spatial_errors, err_spatial)

        @printf("%-10d  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, err_Tmax, err_Trange, err_spatial, all_wall[i])
    end
```

- [ ] **Step 4: 替换图 3 和图 4**

将图 3（收敛误差）和图 4（空间均匀性）替换为：

```julia
    # 图3: RMSPE 收敛折线图
    xp = collect(1:length(THERMAL_Nθ))
    nθ_labels = ["$n" for n in THERMAL_Nθ]
    p3 = plot(xlabel="Mesh No.", ylabel="RMSPE [%]",
              title="Thermal Mesh: Convergence Error (RMSPE)",
              xticks=(xp, nθ_labels), legend=:topright)
    plot!(p3, xp, Tmax_errors, marker=:o, lw=2, label="T_max(t) RMSPE", color=:blue)
    plot!(p3, xp, Trange_errors, marker=:s, lw=2, label="T_range(t) RMSPE", color=:orange)
    plot!(p3, xp, Spatial_errors, marker=:d, lw=2, label="Spatial RMSPE", color=:green)
    hline!([5.0], label="5%", color=:red, ls=:dash)
    hline!([1.0], label="1%", color=:green, ls=:dash)
    mesh_txt = "1→nθ=20\n2→nθ=40\n3→nθ=80\n4→nθ=160 ref"
    annotate!(p3, maximum(xp) - 0.1, maximum(Tmax_errors) * 0.9,
              text(mesh_txt, 8, :right, :top, :Courier))
    savefig(p3, joinpath(out_dir, "thermal_convergence_error.png"))

    # 图4: 空间场 RMSPE 随 nθ 变化
    p4 = plot(THERMAL_Nθ, Spatial_errors,
              xlabel="nθ (per revolution)", ylabel="Spatial RMSPE [%]",
              title="Thermal Mesh: Spatial Field RMSPE",
              marker=:diamond, lw=2, color=:purple, label="Spatial RMSPE",
              yscale=:log10)
    hline!([5.0], label="5%", color=:red, ls=:dash)
    savefig(p4, joinpath(out_dir, "thermal_spatial_rmspe.png"))
```

- [ ] **Step 5: Commit**

```bash
git add "example/网格敏感性/3_thermal_mesh_sensitivity.jl"
git commit -m "refactor: replace point-value metrics with RMSPE in thermal script"
```

---

## Chunk 4: CZM Track (Script 4)

### Task 4: 修改 Script 4 后处理

**Files:**
- Modify: `example/网格敏感性/4_czm_mesh_sensitivity.jl`

- [ ] **Step 1: 替换文件头部**

将 `using Printf, Plots` 修改为：

```julia
using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_rmspe_utils.jl"))
```

- [ ] **Step 2: 替换 docstring**

```julia
"""
Script 4: CZM 网格敏感性分析

先运行 Script 1 确定 nθ 区间，再用 4 组 nθ 运行全耦合模型
（热-化学载荷驱动 CZM 损伤），对比 D_max(t)、n_fractured(t)、
δ_max_n(t) 的 RMSPE 和牵引-分离面积偏差。

指标：D_max(t) RMSPE、n_frac(t) RMSPE、δ_max_n(t) RMSPE、Traction-Sep 面积偏差
输出：收敛图
"""
```

- [ ] **Step 3: 替换 main() 中的汇总逻辑**

将第 121-140 行（从 `# ── 参考解` 到汇总循环结束）替换为：

```julia
    # ── 参考解 (最细) ──
    ref_idx = length(nθ_list)
    ref_result = all_results[nθ_list[ref_idx]]
    t_ref = ref_result["time [s]"]

    println("\n" * "=" ^ 70)
    println("CZM 网格收敛性汇总 (RMSPE)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %12s  %12s  %12s  %12s  %10s\n",
            "nθ", "n_coh",
            "D_max RMSPE%", "n_frac RMSPE%", "δ_max_n RMSPE%", "Area err%",
            "Wall [s]")
    println("-" ^ 90)

    D_rmspe = Float64[]
    nfrac_rmspe = Float64[]
    delta_rmspe = Float64[]
    area_errors = Float64[]

    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        t_c = r["time [s]"]
        n_coh = all_czm[nθ].n_cohesive

        if i == ref_idx
            push!(D_rmspe, 0.0)
            push!(nfrac_rmspe, 0.0)
            push!(delta_rmspe, 0.0)
            push!(area_errors, 0.0)
            @printf("%-10d  %10d  %12s  %12s  %12s  %12s  %10.2f\n",
                    nθ, n_coh, "ref", "ref", "ref", "ref", all_wall[i])
            continue
        end

        # D_max(t) RMSPE
        if haskey(r, "czm D_max") && haskey(ref_result, "czm D_max")
            D_aligned = align_to_ref(t_c, r["czm D_max"], t_ref)
            err_D, _ = rmspe(D_aligned, ref_result["czm D_max"])
        else
            err_D = NaN
        end

        # n_fractured(t) RMSPE
        if haskey(r, "czm n_fractured") && haskey(ref_result, "czm n_fractured")
            nf_aligned = align_to_ref(t_c, r["czm n_fractured"], t_ref)
            err_nf, _ = rmspe(nf_aligned, ref_result["czm n_fractured"])
        else
            err_nf = NaN
        end

        # δ_max_n(t) RMSPE
        if haskey(r, "czm δ_max_n [m]") && haskey(ref_result, "czm δ_max_n [m]")
            d_aligned = align_to_ref(t_c, r["czm δ_max_n [m]"], t_ref)
            err_d, _ = rmspe(d_aligned, ref_result["czm δ_max_n [m]"])
        else
            err_d = NaN
        end

        # 牵引-分离面积偏差（峰值损伤单元）
        err_area = NaN
        if haskey(r, "czm traction normal [Pa]") && haskey(ref_result, "czm traction normal [Pa]")
            traction_c = r["czm traction normal [Pa]"]
            sep_c = r["czm separation normal [m]"]
            traction_ref = ref_result["czm traction normal [Pa]"]
            sep_ref = ref_result["czm separation normal [m]"]

            if ndims(traction_c) == 2 && ndims(traction_ref) == 2
                # 选取最终时刻 D 值最大的单元
                D_c = haskey(r, "czm damage [0-1]") ? r["czm damage [0-1]"] : nothing
                D_r = haskey(ref_result, "czm damage [0-1]") ? ref_result["czm damage [0-1]"] : nothing

                peak_c = D_c !== nothing && ndims(D_c) == 2 ? argmax(D_c[:, end]) : 1
                peak_r = D_r !== nothing && ndims(D_r) == 2 ? argmax(D_r[:, end]) : 1

                err_area = area_error(
                    sep_c[peak_c, :], traction_c[peak_c, :],
                    sep_ref[peak_r, :], traction_ref[peak_r, :])
            end
        end

        push!(D_rmspe, err_D)
        push!(nfrac_rmspe, err_nf)
        push!(delta_rmspe, err_d)
        push!(area_errors, err_area)

        @printf("%-10d  %10d  %12.4f  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, n_coh, err_D, err_nf, err_d, err_area, all_wall[i])
    end
```

- [ ] **Step 4: 替换图 4（收敛误差）**

将图 4 部分替换为：

```julia
    # 图4: RMSPE 收敛误差
    if ref_idx > 1
        p4 = plot(xlabel="nθ (per revolution)", ylabel="RMSPE [%]",
                  title="CZM Mesh: Convergence Error (RMSPE)",
                  legend=:topright, yscale=:log10)
        plot!(p4, nθ_list[1:end-1], D_rmspe[1:end-1],
              marker=:o, lw=2, label="D_max(t)", color=:blue)
        plot!(p4, nθ_list[1:end-1], nfrac_rmspe[1:end-1],
              marker=:s, lw=2, label="n_frac(t)", color=:orange)
        plot!(p4, nθ_list[1:end-1], delta_rmspe[1:end-1],
              marker=:d, lw=2, label="δ_max_n(t)", color=:green)
        plot!(p4, nθ_list[1:end-1], area_errors[1:end-1],
              marker=:p, lw=2, label="Traction-Sep area", color=:purple)
        hline!(p4, [5.0], label="5%", color=:red, ls=:dash)
        savefig(p4, joinpath(out_dir, "czm_convergence_error.png"))
    end
```

- [ ] **Step 5: Commit**

```bash
git add "example/网格敏感性/4_czm_mesh_sensitivity.jl"
git commit -m "refactor: replace point-value metrics with RMSPE in CZM script"
```

---

## Chunk 5: 能量守恒 (Script 5)

### Task 5: 修改 Script 5 后处理

**Files:**
- Modify: `example/网格敏感性/5_energy_conservation_check.jl`

- [ ] **Step 1: 替换文件头部**

在 `using Printf, Plots, LinearAlgebra` 之后添加：

```julia
using Statistics
include(joinpath(@__DIR__, "0_rmspe_utils.jl"))
```

- [ ] **Step 2: 在第 179 行（`max ε_R` 输出之后）插入归一化 RMS 残余**

在 `@printf("  max ε_R    = %.4f %%\n", maximum(ε_R[2:end]))` 之后插入：

```julia
    # ── 归一化 RMS 残余 (spec §3.4) ──
    ε_R_rms = sqrt(mean(R[2:end].^2)) / abs(W_elec[end]) * 100
    @printf("  ε_R,rms   = %.4f %%\n", ε_R_rms)
```

- [ ] **Step 3: 在汇总表中添加 ε_R,rms 行**

在 `@printf("  %-25s %15.4f %%\n", "ε_R    (相对误差)", ε_R[end])` 之后插入：

```julia
    @printf("  %-25s %15.4f %%\n", "ε_R,rms (归一化RMS)", ε_R_rms)
```

- [ ] **Step 4: Commit**

```bash
git add "example/网格敏感性/5_energy_conservation_check.jl"
git commit -m "feat: add normalized RMS residual to energy conservation check"
```

---

## Chunk 6: 更新规划文档

### Task 6: 更新 findings.md 和 progress.md

**Files:**
- Modify: `docs/planning-with-files/网格敏感性分析/findings.md`
- Modify: `docs/planning-with-files/网格敏感性分析/progress.md`

- [ ] **Step 1: 在 findings.md 末尾追加统计指标设计记录**

```markdown
## 8. 统计指标体系升级（2026-04-29）

### 8.1 问题

原有三个 Track 的所有对比指标都基于单点值（snapshot），对局部波动敏感，无法反映整条曲线或空间场的整体收敛质量。

### 8.2 解决方案

将所有指标统一替换为基于 RMSPE 的统计量：

- **电化学 Track**: V(t) RMSPE、T(t) RMSPE、dT/dt(t) RMSPE
- **热学 Track**: T_max(t) RMSPE、T_range(t) RMSPE、空间场 RMSPE（时间平均）
- **CZM Track**: D_max(t) RMSPE、n_frac(t) RMSPE、δ_max_n(t) RMSPE、牵引-分离面积偏差
- **能量守恒**: 保留 ε_R(t) 瞬时值，新增归一化 RMS 残余 ε_R,rms

### 8.3 关键设计决策

- 误差公式：RMSPE（相对均方根百分比误差），带零点保护（threshold = 1e-3 * max|y_ref|）
- 验收阈值：统一 5%
- 时间对齐：手写线性插值 `align_to_ref`，不依赖外部包
- 牵引-分离面积偏差：选取最终时刻 D 值最大的单元进行对比
- 能量残余：归一化 RMS（非 RMSPE，避免除零）

### 8.4 旧指标放弃理由

- 角变化收敛：Bi_t ≈ 0.004 导致角变化极小，RMSPE 零点保护大量跳过
- 应力峰值：空间分布不均匀，单点意义有限
- 损伤起始时间：事件时间，RMSPE 不适用
- 载荷-位移曲线：需纯机械位移 BC，与电池实际驱动不符

### 8.5 参考文件

- Spec: `docs/superpowers/specs/2026-04-29-grid-sensitivity-statistical-metrics-design.md`
- Plan: `docs/superpowers/plans/2026-04-29-statistical-metrics-implementation-plan.md`
```

- [ ] **Step 2: 在 progress.md 末尾追加进度记录**

```markdown
### 2026-04-29 统计指标升级

**已完成：**

- 创建规格文档 `specs/2026-04-29-grid-sensitivity-statistical-metrics-design.md`
- 规格评审通过（2 轮，修复 3 严重 + 5 重要问题）
- 创建共享工具函数 `example/网格敏感性/0_rmspe_utils.jl`
- 修改 Script 2-5 后处理逻辑

**当前状态：**

- 所有脚本后处理已从点值对比替换为 RMSPE 统计指标
- 尚未运行修改后的脚本
- 下一步：运行脚本验证输出格式

**下一步：**

1. 运行 Script 2 验证电化学 RMSPE 输出
2. 运行 Script 3 验证热学 RMSPE 输出
3. 运行 Script 4 验证 CZM RMSPE 输出
4. 运行 Script 5 验证能量守恒 RMS 输出
5. 汇总结果回写到 findings.md
```

- [ ] **Step 3: Commit**

```bash
git add "docs/planning-with-files/网格敏感性分析/findings.md" "docs/planning-with-files/网格敏感性分析/progress.md"
git commit -m "docs: update findings and progress for statistical metrics upgrade"
```
