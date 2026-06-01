# GCI 网格收敛性分析实施计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 `example/网格敏感性_v2/` 目录下实现基于 Roache GCI + L2/H1 范数的网格收敛性分析脚本集

**Architecture:** 从原版 `example/网格敏感性/` 独立，新建 `example/网格敏感性_v2/` 目录。核心工具函数 `0_convergence_utils.jl` 提供 GCI 框架（Richardson 外推、观测收敛阶、渐近检查）和离散范数（L2、L∞）。Script 2-4 分别处理电化学、热学、CZM 三个 Track，每个脚本遵循统一的收敛分析流程。Script 5 处理能量守恒检查。

**Tech Stack:** Julia, JuBat 框架, Plots.jl, Printf, Statistics, LinearAlgebra

**Spec:** `docs/superpowers/specs/2026-06-01-grid-convergence-gci-design.md`

---

## File Structure

```
example/网格敏感性_v2/
├── 0_convergence_utils.jl              # 核心：GCI + 范数 + 插值工具函数
├── 1_cohesive_characteristic_length.jl  # 从原版复制，不变
├── 2_electrochemical_convergence.jl     # Script 2: 电化学 Track
├── 3_thermal_convergence.jl             # Script 3: 热学 Track
├── 4_czm_convergence.jl                 # Script 4: CZM Track
└── 5_energy_conservation.jl             # Script 5: 能量守恒检查
```

```
output/mesh_convergence/                 # 所有输出图
```

---

## Chunk 1: 工具函数

### Task 1: 创建目录和复制不变文件

**Files:**
- Create: `example/网格敏感性_v2/` (directory)
- Copy: `example/网格敏感性/1_cohesive_characteristic_length.jl` → `example/网格敏感性_v2/1_cohesive_characteristic_length.jl`

- [ ] **Step 1: 创建目录**

```bash
mkdir -p "example/网格敏感性_v2"
```

- [ ] **Step 2: 复制 Script 1**

```bash
cp "example/网格敏感性/1_cohesive_characteristic_length.jl" "example/网格敏感性_v2/1_cohesive_characteristic_length.jl"
```

- [ ] **Step 3: 创建输出目录**

```bash
mkdir -p output/mesh_convergence
```

- [ ] **Step 4: Commit**

```bash
git add "example/网格敏感性_v2/"
git commit -m "chore: create mesh_convergence_v2 directory and copy Script 1"
```

---

### Task 2: 实现 `0_convergence_utils.jl` — 离散范数函数

**Files:**
- Create: `example/网格敏感性_v2/0_convergence_utils.jl`

这是整个框架的基础。先实现离散范数部分。

- [ ] **Step 1: 创建文件头部和离散范数函数**

```julia
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
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/0_convergence_utils.jl"
git commit -m "feat: add discrete norm functions to convergence utils"
```

---

### Task 3: 实现 `0_convergence_utils.jl` — GCI 框架函数

**Files:**
- Modify: `example/网格敏感性_v2/0_convergence_utils.jl`

追加 GCI 框架核心函数。

- [ ] **Step 1: 追加 GCI 函数到文件末尾**

```julia
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
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/0_convergence_utils.jl"
git commit -m "feat: add GCI framework functions (observed_order, compute_gci, asymptotic_check, effective_h)"
```

---

### Task 4: 实现 `0_convergence_utils.jl` — 场对齐与插值函数

**Files:**
- Modify: `example/网格敏感性_v2/0_convergence_utils.jl`

- [ ] **Step 1: 追加对齐、插值和辅助函数**

```julia
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
所有坐标为笛卡尔 (x, y)（调用前已完成极坐标→笛卡尔转换）。

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
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/0_convergence_utils.jl"
git commit -m "feat: add field alignment, IDW interpolation, and utility functions"
```

---

## Chunk 2: 电化学 Track (Script 2)

### Task 5: 实现 Script 2 — 电化学网格收敛性分析

**Files:**
- Create: `example/网格敏感性_v2/2_electrochemical_convergence.jl`

**说明**：电化学 Track 使用 4 级 `{Nn, Ns, Np}` 网格配置，lumped thermal 模型。GCI 标量为 $V_{\text{end}}$ 和 $T_{\text{peak}}$。L2 范数用于 V(t)、T(t)、dT/dt(t) 曲线误差。无热网格，所以 `h` 近似为 $\sqrt{1/(Nn + Ns + Np)}$。

- [ ] **Step 1: 创建 Script 2**

```julia
"""
Script 2: 电化学网格收敛性分析 (GCI 框架)

4 组 (Nn, Ns, Np) 下运行 SPMe + lumped thermal：
  - (40, 20, 40)  ← 参考解
  - (20, 10, 20)
  - (20,  5, 20)
  - (10,  5, 10)

指标：V(t) L2_rel, T(t) L2_rel, dT/dt L2
GCI 标量：V_end, T_peak, max|dT/dt|
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const MESH_CONFIGS = [
    (40, 20, 40),  # Level 1 (参考解)
    (20, 10, 20),  # Level 2
    (20,  5, 20),  # Level 3
    (10,  5, 10),  # Level 4 (最粗)
]

function run_echem_case(param_dim, Nn, Ns, Np)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "lumped"
    opt.dimension = 1
    opt.Nn = Nn; opt.Ns = Ns; opt.Np = Np
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.per_element_spme = false
    opt.mechanicalmodel = "none"
    opt.czm_enabled = false

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    return result
end

"""电化学 Track 的近似 h（无热网格时）"""
function echem_h(Nn, Ns, Np)
    return 1.0 / sqrt(Nn + Ns + Np)
end

function main()
    println("=" ^ 70)
    println("Script 2: 电化学网格收敛性分析 (GCI)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    results = Dict{Tuple{Int,Int,Int}, Dict}()
    wall_times = Float64[]
    h_vals = Float64[]

    for (Nn, Ns, Np) in MESH_CONFIGS
        label = "($Nn, $Ns, $Np)"
        @printf("\n--- 运行 %s ---\n", label)
        t0 = time_ns()
        r = run_echem_case(param_dim, Nn, Ns, Np)
        dt_wall = (time_ns() - t0) * 1e-9
        push!(wall_times, dt_wall)
        results[(Nn, Ns, Np)] = r
        push!(h_vals, echem_h(Nn, Ns, Np))
        @printf("  耗时 %.2f s,  %d 步\n", dt_wall, length(r["time [s]"]))
    end

    # ── 参考解 ──
    ref = results[MESH_CONFIGS[1]]
    t_ref = ref["time [s]"]
    V_ref = ref["cell voltage [V]"]
    T_ref = ref["temperature [K]"]

    # ── 收敛误差表 ──
    V_l2 = Float64[]; T_l2 = Float64[]; dTdt_l2 = Float64[]
    V_linf = Float64[]; T_linf = Float64[]

    println("\n" * "=" ^ 70)
    println("电化学网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-16s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n",
            "Mesh", "h",
            "V L2_rel%", "T L2_rel%", "dT/dt L2",
            "V Linf", "T Linf", "Wall [s]")
    println("-" ^ 110)

    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        t_c = r["time [s]"]
        V_c = r["cell voltage [V]"]
        T_c = r["temperature [K]"]

        V_aligned = align_to_ref(t_c, V_c, t_ref)
        T_aligned = align_to_ref(t_c, T_c, t_ref)

        if i == 1
            push!(V_l2, 0.0); push!(T_l2, 0.0); push!(dTdt_l2, 0.0)
            push!(V_linf, 0.0); push!(T_linf, 0.0)
            @printf("(%2d,%2d,%2d)       %10.4f  %12s  %12s  %12s  %12s  %12s  %10.2f\n",
                    Nn, Ns, Np, h_vals[i], "ref", "ref", "ref", "ref", "ref", wall_times[i])
        else
            err_V_l2 = l2_rel_norm(V_aligned, V_ref) * 100
            err_T_l2 = l2_rel_norm(T_aligned, T_ref) * 100
            err_V_inf = max_norm(V_aligned, V_ref)
            err_T_inf = max_norm(T_aligned, T_ref)

            dTdt_ref_vals = diff(T_ref) ./ diff(t_ref)
            dTdt_cand_vals = diff(T_aligned) ./ diff(t_ref)
            err_dTdt = l2_norm(dTdt_cand_vals, dTdt_ref_vals)

            push!(V_l2, err_V_l2); push!(T_l2, err_T_l2); push!(dTdt_l2, err_dTdt)
            push!(V_linf, err_V_inf); push!(T_linf, err_T_inf)

            @printf("(%2d,%2d,%2d)       %10.4f  %12.4f  %12.4f  %12.4f  %12.6f  %12.6f  %10.2f\n",
                    Nn, Ns, Np, h_vals[i], err_V_l2, err_T_l2, err_dTdt,
                    err_V_inf, err_T_inf, wall_times[i])
        end
    end

    # ── GCI 计算 ──
    # GCI 标量：V_end, T_peak
    V_end_vals = [results[c]["cell voltage [V]"][end] for c in MESH_CONFIGS]
    T_peak_vals = [maximum(results[c]["temperature [K]"]) for c in MESH_CONFIGS]
    dTdt_max_vals = [let
        r = results[c]; t = r["time [s]"]; T = r["temperature [K]"]
        maximum(abs.(diff(T) ./ diff(t)))
    end for c in MESH_CONFIGS]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %12s  %8s  %12s  %12s  %12s\n",
            "Pair", "r", "epsilon", "p_obs", "GCI(V_end)", "GCI(T_peak)", "GCI(dTdt)")
    println("-" ^ 90)

    gci_results = []  # (r21, p, gci_V, gci_T, gci_dTdt)
    for i in 1:length(MESH_CONFIGS)-1
        r21 = h_vals[i+1] / h_vals[i]
        r32 = i+1 < length(MESH_CONFIGS) ? h_vals[i+2] / h_vals[i+1] : r21

        # 三组值计算 p（如果可用）
        if i + 2 <= length(MESH_CONFIGS)
            p_V = observed_order(V_end_vals[i], V_end_vals[i+1], V_end_vals[i+2], r21, r32)
            p_T = observed_order(T_peak_vals[i], T_peak_vals[i+1], T_peak_vals[i+2], r21, r32)
            p_dT = observed_order(dTdt_max_vals[i], dTdt_max_vals[i+1], dTdt_max_vals[i+2], r21, r32)
        else
            p_V = 2.0; p_T = 2.0; p_dTdt = 2.0  # 回退到理论值
        end

        gci_V = compute_gci(V_end_vals[i], V_end_vals[i+1], r21; p=isnan(p_V) ? 2.0 : p_V)
        gci_T = compute_gci(T_peak_vals[i], T_peak_vals[i+1], r21; p=isnan(p_T) ? 2.0 : p_T)
        gci_dT = compute_gci(dTdt_max_vals[i], dTdt_max_vals[i+1], r21; p=isnan(p_dTdt) ? 2.0 : p_dTdt)

        push!(gci_results, (r21, p_V, gci_V, gci_T, gci_dT))

        eps_V = abs(V_end_vals[i+1] - V_end_vals[i]) / abs(V_end_vals[i])
        @printf("L%d→L%d       %8.3f  %12.6f  %8.2f  %12.4f  %12.4f  %12.4f\n",
                i, i+1, r21, eps_V,
                isnan(p_V) ? -1.0 : p_V,
                gci_V, gci_T, gci_dT)
    end

    # ── 渐近收敛检查 ──
    if length(gci_results) >= 2
        println("\n渐近收敛检查:")
        for (name, gcis) in [("V_end", [r[3] for r in gci_results]),
                             ("T_peak", [r[4] for r in gci_results]),
                             ("dTdt", [r[5] for r in gci_results])]
            ratio = asymptotic_check(gcis[1], gcis[2],
                                     gci_results[1][1], gci_results[1][2])
            @printf("  %-8s: GCI_12/GCI_23 ratio = %.4f %s\n",
                    name, isnan(ratio) ? 0.0 : ratio,
                    isnan(ratio) ? "(insufficient data)" :
                    0.8 <= ratio <= 1.2 ? "✓ asymptotic" : "✗ not asymptotic")
        end
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    colors = [:blue, :orange, :green, :red]

    # 图1: 电压曲线对比
    p1 = plot(xlabel="Time [s]", ylabel="Cell Voltage [V]",
              title="Electrochemical Mesh: Voltage", legend=:topright)
    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        label = i == 1 ? "ref ($Nn,$Ns,$Np)" : "($Nn,$Ns,$Np)"
        ls = i == 1 ? :solid : :dash
        plot!(p1, r["time [s]"], r["cell voltage [V]"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p1, joinpath(out_dir, "echem_voltage_convergence.png"))

    # 图2: 温度曲线对比
    p2 = plot(xlabel="Time [s]", ylabel="Temperature [K]",
              title="Electrochemical Mesh: Temperature", legend=:bottomright)
    for (i, (Nn, Ns, Np)) in enumerate(MESH_CONFIGS)
        r = results[(Nn, Ns, Np)]
        label = i == 1 ? "ref ($Nn,$Ns,$Np)" : "($Nn,$Ns,$Np)"
        ls = i == 1 ? :solid : :dash
        plot!(p2, r["time [s]"], r["temperature [K]"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p2, joinpath(out_dir, "echem_temperature_convergence.png"))

    # 图3: log-log 收敛误差图（含理论斜率线）
    h_plot = h_vals[2:end]
    V_plot = V_l2[2:end]; T_plot = T_l2[2:end]

    p3 = plot(xlabel="h (element size)", ylabel="L2_rel error [%]",
              title="Electrochemical Mesh: Convergence",
              xscale=:log10, yscale=:log10, legend=:topright)

    plot!(p3, h_plot, V_plot, marker=:o, lw=2, label="V(t) L2_rel", color=:blue)
    plot!(p3, h_plot, T_plot, marker=:s, lw=2, label="T(t) L2_rel", color=:orange)

    # 理论斜率线 O(h^2)
    if length(h_plot) >= 2
        h_ref_line = [minimum(h_plot), maximum(h_plot)]
        # 基于最细网格点的误差值，画理论斜率
        e_ref_V = V_plot[1]
        e_ref_T = T_plot[1]
        h_ref_val = h_plot[1]
        theory_V = e_ref_V .* (h_ref_line ./ h_ref_val).^2
        theory_T = e_ref_T .* (h_ref_line ./ h_ref_val).^2
        plot!(p3, h_ref_line, theory_V, lw=1, ls=:dash, color=:blue, label="O(h²) V ref")
        plot!(p3, h_ref_line, theory_T, lw=1, ls=:dash, color=:orange, label="O(h²) T ref")
    end

    # 标注观测收敛阶
    if !isnan(p_V) || !isnan(p_T)
        annotate_text = ""
        !isnan(p_V) && (annotate_text *= "p_obs(V)=$(round(p_V, digits=2))\n")
        !isnan(p_T) && (annotate_text *= "p_obs(T)=$(round(p_T, digits=2))\n")
        annotate_text *= "p_theory=2.0"
        annotate!(p3, maximum(h_plot), minimum([V_plot; T_plot]) * 1.5,
                  text(annotate_text, 8, :left, :bottom, :Courier))
    end

    savefig(p3, joinpath(out_dir, "echem_error_convergence.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/2_electrochemical_convergence.jl"
git commit -m "feat: implement electrochemical Track convergence analysis with GCI"
```

---

## Chunk 3: 热学 Track (Script 3)

### Task 6: 实现 Script 3 — 热学网格收敛性分析

**Files:**
- Create: `example/网格敏感性_v2/3_thermal_convergence.jl`

**说明**：热学 Track 使用 4 级 `nθ`，纯热模型（均匀体积热源 + 表面冷却）。关键新特性是跨网格 IDW 插值（使用笛卡尔坐标）和空间场 L2 范数。

- [ ] **Step 1: 创建 Script 3**

```julia
"""
Script 3: 热学网格收敛性分析 (GCI 框架)

4 组 nθ = {20, 40, 80, 160} 下运行纯热模型（均匀体积热源 + 表面冷却），
基于 Jellyroll spiral mesh，使用 opt.model = "thermal" 直接求解热方程。

指标：空间场 L2_rel, T_max(t) L2_rel, T_range(t) L2
GCI 标量：T_max(t_end)
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

const THERMAL_Nθ = [20, 40, 80, 160]
const Q0 = 2.0e5  # 均匀体积热源 [W/m³]

function run_thermal_case(param_dim, nθ)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "distributed2D"
    opt.solveType = "Crank-Nicolson"
    opt.time = [0.0, 3600]
    opt.dt = [5, 5]
    opt.gsorder = 2
    opt.Nn = 5; opt.Ns = 3; opt.Np = 5
    opt.Nrn = 5; opt.Nrp = 5

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    # 切换为纯热模式
    case.opt.model = "thermal"

    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    nnode = mesh.nlen
    scale = param_dim.scale

    q0_nd = Q0 / scale.q

    variables = Dict{String,Any}()
    variables["thermal2D temperature at nodes"] = fill(param_dim.cell.T0 / scale.T_ref, nnode)
    variables["heat_source_fields"] = fill(q0_nd, ne)

    update_fn = (t, vars) -> begin
        vars["heat_source_fields"] = fill(q0_nd, ne)
    end

    solve_t0 = time_ns()
    result = JuBat.Solve(case;
        thermal_variables=variables,
        thermal_update_fn=update_fn,
        thermal_record=true)
    solve_wall = (time_ns() - solve_t0) * 1e-9

    T_nodes = result.T_nodes .* scale.T_ref
    times   = result.time .* scale.t0
    T_hist  = result.T_hist .* scale.T_ref

    return (T_nodes=T_nodes, times=times, T_hist=T_hist,
            mesh=mesh, mesh_data=mesh_data, solve_wall=solve_wall)
end

function main()
    println("=" ^ 70)
    println("Script 3: 热学网格收敛性分析 (GCI, q0=$(Q0) W/m³)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")

    all_results = Dict{Int, Any}()
    h_vals = Float64[]
    wall_times = Float64[]

    for nθ in THERMAL_Nθ
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        r = run_thermal_case(param_dim, nθ)
        all_results[nθ] = r

        h = effective_h(r.mesh)
        push!(h_vals, h)
        push!(wall_times, r.solve_wall)

        T_final = r.T_nodes
        @printf("  耗时 %.2f s,  T_max = %.3f K,  T_min = %.3f K,  h = %.6f\n",
            r.solve_wall, maximum(T_final), minimum(T_final), h)
    end

    # ── 参考解 (nθ=160) ──
    ref_idx = length(THERMAL_Nθ)
    ref_result = all_results[THERMAL_Nθ[ref_idx]]
    T_hist_ref = ref_result.T_hist
    t_ref = ref_result.times
    nt_ref = length(t_ref)
    mesh_ref = ref_result.mesh

    # 参考节点笛卡尔坐标
    ref_x = mesh_ref.node[:, 1]
    ref_y = mesh_ref.node[:, 2]

    Tmax_ref_curve = [maximum(T_hist_ref[:, k]) for k in 1:nt_ref]
    Trange_ref_curve = [maximum(T_hist_ref[:, k]) - minimum(T_hist_ref[:, k]) for k in 1:nt_ref]

    # ── 收敛误差表 ──
    Tmax_l2 = Float64[]; Trange_l2 = Float64[]; Spatial_l2 = Float64[]
    Tmax_linf = Float64[]

    println("\n" * "=" ^ 70)
    println("热学网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %12s  %12s  %12s  %12s  %10s\n",
            "nθ", "h", "Tmax L2_rel%", "Trange L2", "Spatial L2%", "Tmax Linf", "Wall [s]")
    println("-" ^ 100)

    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        t_c = r.times
        nt_c = length(t_c)
        mesh_c = r.mesh

        if i == ref_idx
            push!(Tmax_l2, 0.0); push!(Trange_l2, 0.0)
            push!(Spatial_l2, 0.0); push!(Tmax_linf, 0.0)
            @printf("%-10d  %10.6f  %12s  %12s  %12s  %12s  %10.2f\n",
                    nθ, h_vals[i], "ref", "ref", "ref", "ref", wall_times[i])
            continue
        end

        # T_max(t) 曲线 L2_rel
        Tmax_cand = [maximum(T_hist[:, k]) for k in 1:nt_c]
        Tmax_cand_aligned = align_to_ref(t_c, Tmax_cand, t_ref)
        err_Tmax = l2_rel_norm(Tmax_cand_aligned, Tmax_ref_curve) * 100
        err_Tmax_inf = max_norm(Tmax_cand_aligned, Tmax_ref_curve)

        # T_range(t) L2（绝对值）
        Trange_cand = [maximum(T_hist[:, k]) - minimum(T_hist[:, k]) for k in 1:nt_c]
        Trange_cand_aligned = align_to_ref(t_c, Trange_cand, t_ref)
        err_Trange = l2_norm(Trange_cand_aligned, Trange_ref_curve)

        # 空间场 L2_rel（需要 IDW 插值到参考网格）
        err_spatial = NaN
        nt_overlap = min(nt_c, nt_ref)
        if size(T_hist, 1) != size(T_hist_ref, 1)
            # 节点数不同，使用 IDW 插值
            coarse_x = mesh_c.node[:, 1]
            coarse_y = mesh_c.node[:, 2]
            spatial_errs = Float64[]
            for k in 1:nt_overlap
                T_cand_interp = interpolate_to_ref_field(
                    T_hist[:, k], coarse_x, coarse_y, ref_x, ref_y)
                val = l2_rel_norm(T_cand_interp, T_hist_ref[:, k])
                isnan(val) || push!(spatial_errs, val)
            end
            err_spatial = isempty(spatial_errs) ? NaN : mean(spatial_errs) * 100
        else
            # 节点数相同，直接计算
            spatial_errs = Float64[]
            for k in 1:nt_overlap
                val = l2_rel_norm(T_hist[:, k], T_hist_ref[:, k])
                isnan(val) || push!(spatial_errs, val)
            end
            err_spatial = isempty(spatial_errs) ? NaN : mean(spatial_errs) * 100
        end

        push!(Tmax_l2, err_Tmax); push!(Trange_l2, err_Trange)
        push!(Spatial_l2, err_spatial); push!(Tmax_linf, err_Tmax_inf)

        @printf("%-10d  %10.6f  %12.4f  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, h_vals[i], err_Tmax, err_Trange, err_spatial,
                err_Tmax_inf, wall_times[i])
    end

    # ── GCI 计算 ──
    Tmax_end_vals = [maximum(all_results[nθ].T_nodes) for nθ in THERMAL_Nθ]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %12s  %8s  %12s\n",
            "Pair", "r", "epsilon", "p_obs", "GCI(T_max)")
    println("-" ^ 60)

    gci_results = []
    for i in 1:length(THERMAL_Nθ)-1
        r21 = h_vals[i+1] / h_vals[i]
        r32 = i+1 < length(THERMAL_Nθ) ? h_vals[i+2] / h_vals[i+1] : r21

        if i + 2 <= length(THERMAL_Nθ)
            p_T = observed_order(Tmax_end_vals[i], Tmax_end_vals[i+1],
                                 Tmax_end_vals[i+2], r21, r32)
        else
            p_T = 2.0
        end

        gci_T = compute_gci(Tmax_end_vals[i], Tmax_end_vals[i+1], r21;
                            p=isnan(p_T) ? 2.0 : p_T)
        push!(gci_results, (r21, p_T, gci_T))

        eps_T = abs(Tmax_end_vals[i+1] - Tmax_end_vals[i]) / abs(Tmax_end_vals[i])
        @printf("L%d→L%d       %8.3f  %12.6f  %8.2f  %12.4f\n",
                i, i+1, r21, eps_T, isnan(p_T) ? -1.0 : p_T, gci_T)
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    colors = [:red, :orange, :green, :blue]

    # 图1: 温度场时间演化
    p1 = plot(xlabel="Time [s]", ylabel="Mean Temperature [K]",
              title="Thermal Mesh: Temperature (q0=$(Q0) W/m³)",
              legend=:bottomright)
    for (i, nθ) in enumerate(THERMAL_Nθ)
        r = all_results[nθ]
        T_hist = r.T_hist
        if !isempty(T_hist)
            nt = size(T_hist, 2)
            T_mean = [mean(T_hist[:, k]) for k in 1:nt]
            label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
            ls = i == ref_idx ? :solid : :dash
            plot!(p1, r.times[1:nt], T_mean, label=label, lw=1.5, color=colors[i], ls=ls)
        end
    end
    savefig(p1, joinpath(out_dir, "thermal_temperature_convergence.png"))

    # 图2: log-log 收敛误差图
    h_plot = h_vals[2:end]
    Tmax_plot = Tmax_l2[2:end]
    Spatial_plot = Spatial_l2[2:end]

    p3 = plot(xlabel="h (element size)", ylabel="L2_rel error [%]",
              title="Thermal Mesh: Convergence",
              xscale=:log10, yscale=:log10, legend=:topright)

    plot!(p3, h_plot, Tmax_plot, marker=:o, lw=2, label="T_max(t) L2_rel", color=:blue)
    plot!(p3, h_plot, Spatial_plot, marker=:d, lw=2, label="Spatial L2_rel", color=:green)

    # 理论斜率线 O(h^2) (Q4, k=1)
    if length(h_plot) >= 2 && Tmax_plot[1] > 0
        h_ref_line = [minimum(h_plot), maximum(h_plot)]
        e_ref = Tmax_plot[1]
        h_ref_val = h_plot[1]
        theory_line = e_ref .* (h_ref_line ./ h_ref_val).^2
        plot!(p3, h_ref_line, theory_line, lw=1, ls=:dash, color=:blue, label="O(h²) theory")
    end

    savefig(p3, joinpath(out_dir, "thermal_spatial_error.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/3_thermal_convergence.jl"
git commit -m "feat: implement thermal Track convergence analysis with GCI and IDW interpolation"
```

---

## Chunk 4: CZM Track (Script 4)

### Task 7: 实现 Script 4 — CZM 网格收敛性分析

**Files:**
- Create: `example/网格敏感性_v2/4_czm_convergence.jl`

**说明**：CZM Track 使用全耦合模型。主指标为断裂能耗散 $E_{\text{frac}}(t)$ 的 GCI。辅指标为 D_max(t) L2_rel、n_fractured(t) L2、牵引-分离面积偏差。

- [ ] **Step 1: 创建 Script 4**

```julia
"""
Script 4: CZM 网格收敛性分析 (GCI 框架)

先运行 Script 1 确定 nθ 区间，再用 4 组 nθ 运行全耦合模型。

主指标：断裂能耗散 E_frac(t) 的 GCI
辅指标：D_max(t) L2_rel, n_frac(t) L2, δ_max_n(t) L2_rel, 牵引-分离面积偏差
输出：收敛误差表、GCI 汇总表、log-log 收敛图
"""

using Printf, Plots, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

"""从参数计算 CZM nθ 候选值"""
function get_czm_nθ(param_dim)
    E_NE  = param_dim.NE.E;  t_NE = param_dim.NE.thickness
    E_PE  = param_dim.PE.E;  t_PE = param_dim.PE.thickness
    E_eff = (E_NE * t_NE + E_PE * t_PE) / (t_NE + t_PE)
    G_c   = param_dim.cohesive.G_c_n
    σ_max = param_dim.cohesive.σ_max_n
    l_c   = G_c * E_eff / σ_max^2

    R_in  = param_dim.cell.Rin
    R_out = param_dim.cell.Rout

    nθ_inner = ceil(Int, 2π * R_in  / l_c)
    nθ_outer = ceil(Int, 2π * R_out / l_c)

    step = max(1, round(Int, (nθ_outer - nθ_inner) / 3))
    nθ_list = [nθ_inner + i * step for i in 0:3]
    nθ_list[4] = nθ_outer

    return nθ_list, l_c, E_eff
end

function run_czm_case(param_dim, nθ)
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
    case.czm_mesh = czm_mesh

    result = JuBat.Solve(case)
    return result, czm_mesh, mesh_data
end

"""计算断裂能耗散 E_frac(t)"""
function compute_fracture_energy(result, czm_mesh, param_dim)
    G_c = param_dim.cohesive.G_c_n
    scale = param_dim.scale
    nt = length(result["time [s]"])

    coh_lengths = [elem.length * scale.L for elem in czm_mesh.cohesive_elements]

    E_frac = zeros(Float64, nt)
    if haskey(result, "czm damage [0-1]")
        D_all = result["czm damage [0-1]"]
        if ndims(D_all) == 2
            n_coh = min(length(coh_lengths), size(D_all, 1))
            for k in 1:min(nt, size(D_all, 2))
                for e in 1:n_coh
                    E_frac[k] += G_c * coh_lengths[e] * D_all[e, k]
                end
            end
        end
    end
    return E_frac
end

function main()
    println("=" ^ 70)
    println("Script 4: CZM 网格收敛性分析 (GCI)")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    nθ_list, l_c, E_eff = get_czm_nθ(param_dim)
    @printf("\nl_c = %.1f μm,  E_eff = %.2f GPa\n", l_c * 1e6, E_eff * 1e-9)
    @printf("CZM nθ 候选值: %s\n", string(nθ_list))

    all_results  = Dict{Int, Dict}()
    all_czm      = Dict{Int, Any}()
    all_wall     = Float64[]
    h_vals       = Float64[]
    E_frac_curves = Dict{Int, Vector{Float64}}()

    for nθ in nθ_list
        @printf("\n--- 运行 nθ = %d ---\n", nθ)
        t0 = time_ns()
        r, czm, mesh_data = run_czm_case(param_dim, nθ)
        dt_wall = (time_ns() - t0) * 1e-9

        all_results[nθ] = r
        all_czm[nθ] = czm

        h = effective_h(mesh_data.thermal2D)
        push!(h_vals, h)
        push!(all_wall, dt_wall)

        E_frac = compute_fracture_energy(r, czm, param_dim)
        E_frac_curves[nθ] = E_frac

        Dmax_end = haskey(r, "czm D_max") ? r["czm D_max"][end] : NaN
        nfrac_end = haskey(r, "czm n_fractured") ? r["czm n_fractured"][end] : NaN

        @printf("  耗时 %.2f s,  D_max = %.4f,  n_frac = %d,  E_frac = %.4e,  h = %.6f\n",
                dt_wall, Dmax_end, Int(nfrac_end), E_frac[end], h)
    end

    # ── 参考解 ──
    ref_idx = length(nθ_list)
    ref_result = all_results[nθ_list[ref_idx]]
    t_ref = ref_result["time [s]"]

    # ── 收敛误差表 ──
    D_l2 = Float64[]; nfrac_l2 = Float64[]; delta_l2 = Float64[]
    area_errors = Float64[]; Efrac_l2 = Float64[]

    println("\n" * "=" ^ 70)
    println("CZM 网格收敛性误差 (L2 范数)")
    println("=" ^ 70)
    @printf("%-10s  %10s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n",
            "nθ", "n_coh", "h",
            "E_frac L2", "D_max L2%", "n_frac L2", "δ L2%", "Area err%", "Wall [s]")
    println("-" ^ 115)

    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        t_c = r["time [s]"]
        n_coh = all_czm[nθ].n_cohesive

        if i == ref_idx
            push!(D_l2, 0.0); push!(nfrac_l2, 0.0); push!(delta_l2, 0.0)
            push!(area_errors, 0.0); push!(Efrac_l2, 0.0)
            @printf("%-10d  %10d  %10.6f  %12s  %12s  %12s  %12s  %12s  %10.2f\n",
                    nθ, n_coh, h_vals[i], "ref", "ref", "ref", "ref", "ref", all_wall[i])
            continue
        end

        # E_frac(t) L2（主指标）
        E_aligned = align_to_ref(t_c, E_frac_curves[nθ], t_ref)
        err_Efrac = l2_norm(E_aligned, E_frac_curves[nθ_list[ref_idx]])
        push!(Efrac_l2, err_Efrac)

        # D_max(t) L2_rel
        err_D = NaN
        if haskey(r, "czm D_max") && haskey(ref_result, "czm D_max")
            D_aligned = align_to_ref(t_c, r["czm D_max"], t_ref)
            err_D = l2_rel_norm(D_aligned, ref_result["czm D_max"]) * 100
        end
        push!(D_l2, err_D)

        # n_fractured(t) L2
        err_nf = NaN
        if haskey(r, "czm n_fractured") && haskey(ref_result, "czm n_fractured")
            nf_aligned = align_to_ref(t_c, r["czm n_fractured"], t_ref)
            err_nf = l2_norm(nf_aligned, ref_result["czm n_fractured"])
        end
        push!(nfrac_l2, err_nf)

        # δ_max_n(t) L2_rel
        err_d = NaN
        if haskey(r, "czm δ_max_n [m]") && haskey(ref_result, "czm δ_max_n [m]")
            d_aligned = align_to_ref(t_c, r["czm δ_max_n [m]"], t_ref)
            err_d = l2_rel_norm(d_aligned, ref_result["czm δ_max_n [m]"]) * 100
        end
        push!(delta_l2, err_d)

        # 牵引-分离面积偏差
        err_area = NaN
        if haskey(r, "czm traction normal [Pa]") && haskey(ref_result, "czm traction normal [Pa]")
            traction_c = r["czm traction normal [Pa]"]
            sep_c = r["czm separation normal [m]"]
            traction_ref = ref_result["czm traction normal [Pa]"]
            sep_ref = ref_result["czm separation normal [m]"]

            if ndims(traction_c) == 2 && ndims(traction_ref) == 2
                D_c = haskey(r, "czm damage [0-1]") ? r["czm damage [0-1]"] : nothing
                D_r = haskey(ref_result, "czm damage [0-1]") ? ref_result["czm damage [0-1]"] : nothing

                peak_c = D_c !== nothing && ndims(D_c) == 2 ? argmax(D_c[:, end]) : 1
                peak_r = D_r !== nothing && ndims(D_r) == 2 ? argmax(D_r[:, end]) : 1

                err_area = area_error(
                    sep_c[peak_c, :], traction_c[peak_c, :],
                    sep_ref[peak_r, :], traction_ref[peak_r, :])
            end
        end
        push!(area_errors, err_area)

        @printf("%-10d  %10d  %10.6f  %12.4e  %12.4f  %12.4f  %12.4f  %12.4f  %10.2f\n",
                nθ, n_coh, h_vals[i], err_Efrac, err_D, err_nf, err_d, err_area, all_wall[i])
    end

    # ── GCI 计算 ──
    Efrac_end_vals = [E_frac_curves[nθ][end] for nθ in nθ_list]
    Dmax_end_vals = [haskey(all_results[nθ], "czm D_max") ? all_results[nθ]["czm D_max"][end] : NaN for nθ in nθ_list]

    println("\n" * "=" ^ 70)
    println("GCI 汇总表")
    println("=" ^ 70)
    @printf("%-12s  %8s  %8s  %12s  %12s\n",
            "Pair", "r", "p_obs", "GCI(E_frac)", "GCI(D_max)")
    println("-" ^ 65)

    gci_results = []
    for i in 1:length(nθ_list)-1
        r21 = h_vals[i+1] / h_vals[i]
        r32 = i+1 < length(nθ_list) ? h_vals[i+2] / h_vals[i+1] : r21

        if i + 2 <= length(nθ_list)
            p_E = observed_order(Efrac_end_vals[i], Efrac_end_vals[i+1],
                                 Efrac_end_vals[i+2], r21, r32)
            p_D = observed_order(Dmax_end_vals[i], Dmax_end_vals[i+1],
                                 Dmax_end_vals[i+2], r21, r32)
        else
            p_E = NaN; p_D = NaN
        end

        gci_E = compute_gci(Efrac_end_vals[i], Efrac_end_vals[i+1], r21;
                            p=isnan(p_E) ? 2.0 : p_E)
        gci_D = isnan(Dmax_end_vals[i]) ? NaN :
                compute_gci(Dmax_end_vals[i], Dmax_end_vals[i+1], r21;
                            p=isnan(p_D) ? 2.0 : p_D)

        push!(gci_results, (r21, p_E, gci_E, gci_D))

        @printf("L%d→L%d       %8.3f  %8.2f  %12.4f  %12.4f\n",
                i, i+1, r21,
                isnan(p_E) ? -1.0 : p_E,
                gci_E, isnan(gci_D) ? NaN : gci_D)
    end

    # ── 绘图 ──
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    colors = [:red, :orange, :green, :blue]

    # 图1: D_max 演化对比
    p1 = plot(xlabel="Time [s]", ylabel="D_max",
              title="CZM Mesh: Damage Evolution", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        haskey(r, "czm D_max") || continue
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p1, r["time [s]"], r["czm D_max"],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p1, joinpath(out_dir, "czm_damage_evolution.png"))

    # 图2: 断裂能耗散对比（新增主指标图）
    p2 = plot(xlabel="Time [s]", ylabel="E_frac [J/m]",
              title="CZM Mesh: Fracture Energy Dissipation", legend=:topleft)
    for (i, nθ) in enumerate(nθ_list)
        r = all_results[nθ]
        label = i == ref_idx ? "ref nθ=$nθ" : "nθ=$nθ"
        ls = i == ref_idx ? :solid : :dash
        plot!(p2, r["time [s]"], E_frac_curves[nθ],
              label=label, lw=1.5, color=colors[i], ls=ls)
    end
    savefig(p2, joinpath(out_dir, "czm_fracture_energy.png"))

    # 图3: log-log 收敛误差图
    if ref_idx > 1
        h_plot = h_vals[1:end-1]
        Efrac_plot = Efrac_l2[1:end-1]
        D_plot = D_l2[1:end-1]

        p4 = plot(xlabel="h (element size)", ylabel="L2 error",
                  title="CZM Mesh: Convergence",
                  xscale=:log10, yscale=:log10, legend=:topright)
        plot!(p4, h_plot, Efrac_plot, marker=:o, lw=2,
              label="E_frac(t) L2", color=:blue)
        plot!(p4, h_plot, D_plot, marker=:s, lw=2,
              label="D_max(t) L2_rel%", color=:orange)

        # 理论斜率线
        if length(h_plot) >= 2 && Efrac_plot[1] > 0
            h_ref_line = [minimum(h_plot), maximum(h_plot)]
            e_ref = Efrac_plot[1]
            h_ref_val = h_plot[1]
            theory_line = e_ref .* (h_ref_line ./ h_ref_val).^2
            plot!(p4, h_ref_line, theory_line, lw=1, ls=:dash,
                  color=:blue, label="O(h²) theory")
        end

        savefig(p4, joinpath(out_dir, "czm_error_convergence.png"))
    end

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
end

main()
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/4_czm_convergence.jl"
git commit -m "feat: implement CZM Track convergence analysis with fracture energy GCI"
```

---

## Chunk 5: 能量守恒 (Script 5)

### Task 8: 实现 Script 5 — 能量守恒检查

**Files:**
- Create: `example/网格敏感性_v2/5_energy_conservation.jl`

**说明**：从原版 Script 5 重写，保留相同的物理公式。新增归一化 RMS 残余报告。单网格 (nθ=80) 运行。

- [ ] **Step 1: 创建 Script 5**

```julia
"""
Script 5: 全场耦合能量守恒检查

在全耦合算例（SPMe + distributed2D thermal + CZM）上验证：
  Q_generated = ΔE_thermal + Q_loss + ΔE_elastic + E_fracture

输出：R(t)、ε_R(t) 图、归一化 RMS 残余
"""

using Printf, Plots, LinearAlgebra, Statistics

include(joinpath(@__DIR__, "0_convergence_utils.jl"))

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

function main()
    println("=" ^ 70)
    println("Script 5: 全场耦合能量守恒检查")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    nθ = 80

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.dimension = 1
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.solveType = "Crank-Nicolson"
    opt.dtType = "auto"
    opt.dt = [0.5, 10.0]
    opt.time = [0.0, 3600]

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.czm_iter_method = "basic"
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3

    I1C = param_dim.cell.I1C
    opt.Current = x -> I1C

    println("\n[1] 创建案例 (nθ=$nθ, 全耦合)...")
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
    case.czm_mesh = czm_mesh
    n_coh = czm_mesh.n_cohesive

    @printf("  热单元: %d,  节点: %d,  CZM单元: %d\n", ne, mesh_th.nlen, n_coh)

    println("\n[2] 运行全耦合求解...")
    t0 = time_ns()
    result = JuBat.Solve(case)
    dt_wall = (time_ns() - t0) * 1e-9
    @printf("  耗时 %.2f s\n", dt_wall)

    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I_arr = result["cell current [A]"]
    nt = length(t)

    println("\n[3] 计算能量平衡...")

    # 电功 W_elec = ∫ V·I dt
    W_elec = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        W_elec[k] = W_elec[k-1] + 0.5 * (V[k]*I_arr[k] + V[k-1]*I_arr[k-1]) * dt_k
    end
    @printf("  W_elec(终) = %.4e J\n", W_elec[end])

    # 热能变化 ΔE_th
    scale = param_dim.scale
    ρ  = param_dim.cell.rho
    cp = param_dim.cell.heat_Q
    Vol = param_dim.cell.volume
    T_mean = result["temperature [K]"]
    T0 = T_mean[1]
    C_total = ρ * cp * Vol

    ΔE_th = C_total .* (T_mean .- T0)
    @printf("  ΔE_th(终)  = %.4e J  (T: %.2f → %.2f K)\n", ΔE_th[end], T0, T_mean[end])

    # 断裂能 E_frac
    G_c = param_dim.cohesive.G_c_n
    coh_lengths = [elem.length * scale.L for elem in czm_mesh.cohesive_elements]

    E_frac = zeros(Float64, nt)
    if haskey(result, "czm damage [0-1]")
        D_all = result["czm damage [0-1]"]
        if ndims(D_all) == 2
            for k in 1:nt
                for e in 1:min(n_coh, size(D_all, 1))
                    E_frac[k] += G_c * coh_lengths[e] * D_all[e, k]
                end
            end
        end
    end
    @printf("  E_frac(终) = %.4e J\n", E_frac[end])

    # 边界热损失 Q_loss
    h_cool   = param_dim.cell.h
    T_amb    = param_dim.cell.T_amb
    A_surface = param_dim.cell.cooling_surface

    Q_loss_inst = h_cool .* A_surface .* (T_mean .- T_amb)
    Q_loss = zeros(Float64, nt)
    for k in 2:nt
        dt_k = t[k] - t[k-1]
        Q_loss[k] = Q_loss[k-1] + 0.5 * (Q_loss_inst[k] + Q_loss_inst[k-1]) * dt_k
    end
    @printf("  Q_loss(终) = %.4e J\n", Q_loss[end])

    # 热源积分 Q_gen
    Q_gen = zeros(Float64, nt)
    if haskey(result, "total heat source [W]")
        Qs = result["total heat source [W]"]
        for k in 2:nt
            dt_k = t[k] - t[k-1]
            Q_gen[k] = Q_gen[k-1] + 0.5 * (Qs[k] + Qs[k-1]) * dt_k
        end
    end
    @printf("  Q_gen(终)  = %.4e J\n", Q_gen[end])

    println("\n[4] 能量残余计算...")
    R = Q_gen .- Q_loss .- ΔE_th .- E_frac

    ε_R = zeros(Float64, nt)
    for k in 1:nt
        denom = max(abs(Q_gen[k]), 1e-12)
        ε_R[k] = abs(R[k]) / denom * 100
    end

    ε_R_rms = sqrt(mean(R[2:end].^2)) / abs(W_elec[end]) * 100

    @printf("  R(终)      = %.4e J\n", R[end])
    @printf("  ε_R(终)    = %.4f %%\n", ε_R[end])
    @printf("  max ε_R    = %.4f %%\n", maximum(ε_R[2:end]))
    @printf("  ε_R,rms   = %.4f %%\n", ε_R_rms)

    println("\n" * "=" ^ 70)
    println("能量平衡汇总")
    println("=" ^ 70)
    @printf("  %-25s %15.4e J\n", "W_elec (电功)", W_elec[end])
    @printf("  %-25s %15.4e J\n", "Q_gen  (热源积分)", Q_gen[end])
    @printf("  %-25s %15.4e J\n", "ΔE_th  (热能变化)", ΔE_th[end])
    @printf("  %-25s %15.4e J\n", "Q_loss (边界热损失)", Q_loss[end])
    @printf("  %-25s %15.4e J\n", "E_frac (断裂能)", E_frac[end])
    @printf("  %-25s %15.4e J\n", "R      (残余)", R[end])
    @printf("  %-25s %15.4f %%\n", "ε_R    (相对误差)", ε_R[end])
    @printf("  %-25s %15.4f %%\n", "ε_R,rms (归一化RMS)", ε_R_rms)

    # ── 绘图 ──
    println("\n[5] 绘图...")
    out_dir = joinpath(root_dir, "output", "mesh_convergence")
    mkpath(out_dir)

    p1 = plot(xlabel="Time [s]", ylabel="Energy [J]",
              title="Energy Components", legend=:topleft)
    plot!(p1, t, Q_gen,     label="Q_gen (heat source)",    lw=2, color=:blue)
    plot!(p1, t, ΔE_th,    label="ΔE_th (thermal)",        lw=2, color=:red)
    plot!(p1, t, Q_loss,   label="Q_loss (boundary)",      lw=2, color=:green)
    plot!(p1, t, E_frac,   label="E_frac (fracture)",      lw=2, color=:purple)
    plot!(p1, t, W_elec,   label="W_elec (electrical)",     lw=2, color=:orange, ls=:dash)
    savefig(p1, joinpath(out_dir, "energy_components.png"))

    p2 = plot(xlabel="Time [s]", ylabel="Residual R [J]",
              title="Energy Balance Residual", legend=:topright)
    plot!(p2, t, R, label="R = Q_gen - Q_loss - ΔE_th - E_frac", lw=2, color=:black)
    hline!(p2, [0.0], label="", color=:gray, ls=:dot)
    savefig(p2, joinpath(out_dir, "energy_residual.png"))

    p3 = plot(xlabel="Time [s]", ylabel="Relative Error ε_R [%]",
              title="Energy Conservation Error", legend=:topright,
              yscale=:log10, ylims=(1e-6, max(100.0, maximum(ε_R[2:end]) * 2)))
    plot!(p3, t[2:end], ε_R[2:end], label="ε_R", lw=2, color=:blue)
    hline!(p3, [1.0],  label="1%",  color=:green, ls=:dash)
    hline!(p3, [5.0],  label="5%",  color=:red,   ls=:dash)
    savefig(p3, joinpath(out_dir, "energy_relative_error.png"))

    p4 = plot(xlabel="Time [s]", ylabel="Energy [J]",
              title="Energy Balance: Q_gen Decomposition", legend=:topleft)
    plot!(p4, t, Q_gen, label="Q_gen", lw=3, color=:black)
    plot!(p4, t, ΔE_th .+ Q_loss .+ E_frac, label="ΔE_th + Q_loss + E_frac",
          lw=2, color=:blue, ls=:dash)
    plot!(p4, t, ΔE_th .+ Q_loss .+ E_frac .+ R, label="+ R (should ≈ Q_gen)",
          lw=2, color=:red, ls=:dot)
    savefig(p4, joinpath(out_dir, "energy_decomposition.png"))

    @printf("\n图已保存到 %s\n", out_dir)
    println("=" ^ 70)
    println("能量守恒检查完成")
    println("=" ^ 70)
end

main()
```

- [ ] **Step 2: Commit**

```bash
git add "example/网格敏感性_v2/5_energy_conservation.jl"
git commit -m "feat: implement energy conservation check script with RMS residual"
```

---

## Chunk 6: 最终验证

### Task 9: 验证工具函数加载

**Files:** None (verification only)

- [ ] **Step 1: 验证 0_convergence_utils.jl 可正确加载**

```bash
cd "D:\OneDrive\Desktop\Jubat For Cursor\JuBat"
julia -e 'include("example/网格敏感性_v2/0_convergence_utils.jl"); println("OK: $(l2_norm([1.0,2.0], [1.1,2.1]))")'
```

Expected: `OK: 0.1` (approximately)

- [ ] **Step 2: 验证 GCI 函数**

```bash
julia -e '
include("example/网格敏感性_v2/0_convergence_utils.jl")
# 测试 observed_order: 精确二次收敛
f1, f2, f3 = 1.0, 1.25, 2.0  # h=0.5, 1.0, 2.0
r21, r32 = 2.0, 2.0
p = observed_order(f1, f2, f3, r21, r32)
println("p = $p (expected ~2.0)")
gci = compute_gci(f1, f2, r21; p=p)
println("GCI = $gci %")
'
```

Expected: `p ≈ 2.0`, `GCI > 0`

- [ ] **Step 3: 验证目录结构**

```bash
ls -la "example/网格敏感性_v2/"
```

Expected: 6 files (0-5)

- [ ] **Step 4: Commit (if any fixes were needed)**

```bash
git add -A "example/网格敏感性_v2/"
git commit -m "fix: address verification issues in convergence scripts"
```

---

### Task 10: 最终提交和文档更新

- [ ] **Step 1: 验证所有文件已提交**

```bash
git status
```

Expected: no uncommitted files in `example/网格敏感性_v2/`

- [ ] **Step 2: Final commit**

```bash
git add -A
git commit -m "feat: complete GCI-based grid convergence analysis suite (v2)

Implements the full mesh convergence analysis framework based on
Roache (1998) GCI methodology per spec 2026-06-01-grid-convergence-gci-design.md:
- 0_convergence_utils.jl: L2/L∞ norms, GCI, observed order, IDW interpolation
- 2_electrochemical_convergence.jl: V(t)/T(t) L2_rel + GCI on V_end/T_peak
- 3_thermal_convergence.jl: Spatial L2_rel + IDW cross-mesh interpolation + GCI
- 4_czm_convergence.jl: Fracture energy dissipation GCI (primary) + D_max L2
- 5_energy_conservation.jl: Energy balance residual check"
```
