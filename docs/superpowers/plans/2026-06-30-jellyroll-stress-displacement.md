# Jellyroll 二维应力/位移场后处理脚本 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 新建 `example/jellyroll_stress_displacement.jl`，在 Jellyroll 电池单次放电（1C, 3600s, nθ=360, czm 关闭）后，于用户指定时间节点调用 `thermal_diffusion_stress_2D` 计算二维应力场与位移场，并绘制一张大图。

**Architecture:** 纯后处理脚本，零侵入、不修改 `src/`。复用 `testexample.jl` 的求解器配置（唯一差异 `czm_enabled=false`），通过从 `result` 重建单时间节点 `variables` 字典喂给已存在的 `JuBat.thermal_diffusion_stress_2D(case, variables)`；用 `Plots.jl` (GR 后端) + 对每个 Q4 单元单独绘制 `Plots.Shape` 多边形按场值上色。

**Tech Stack:** Julia, JuBat (本地 `.JuBat`), `Printf`, `Plots.jl` (GR 后端)。无新依赖。

**Spec:** `docs/superpowers/specs/2026-06-30-jellyroll-stress-displacement-design.md`

---

## ⚠️ 关键领域注意点（实现者必读）

1. **尺度区分（最重要的领域知识）**：本脚本输出的应力/位移是 **极片/电极尺度（coating-scale, 二维宏观平面应力）**，由 `compute_effective_coating_modulus(case)` 返回的全叠合厚度加权有效模量 `E_eff / ν_eff / α_eff` 计算（典型 `E_eff ~ 5e8 Pa`，源于 `PE.E_coat / NE.E_coat`）。这与 **颗粒尺度** `Calstressdisp` 用的 `PE.E / NE.E`（~1e10 Pa）**完全不同**，二者不可混用、不可对比。脚本中**不要**任何引用 `param.PE.E / param.NE.E` 作为宏观应力模量的代码。**验收时 σ_vm 峰值应在 1e7–1e9 Pa，若接近 1e10 Pa 说明误用了颗粒模量。**

2. **量纲已核实**（spec §6.1）：传给 `thermal_diffusion_stress_2D` 的 `T_nodes` 是有量纲 K（取自 `result["thermal2D temperature at nodes [K]"]`）；`thermal2D element soc_n/p` 与函数内 `param.NE.cs0` 同为 cs_max 归一化（[0,1] 区间），**不要**乘除 `scale.cn_max`。

3. **函数返回值**：`thermal_diffusion_stress_2D(case, variables)` 返回一个新的 `variables` 字典（原字典 `copy`），从中读：
   - 单元常数场（长度 `ne`）：`"diffusion stress xx"`, `... yy"`, `... xy"`, `... vonMises"`（Pa）
   - 节点场（长度 `nlen`）：`"displacement x"`, `"displacement y"`（m）

4. **本脚本属于 `example/`**：仓库内 `example/` 下脚本无单元测试传统。本计划用 **in-script `@assert` + 末尾 smoke-run** 替代正式 TDD；每个任务末尾 commit。

5. **Plots 后端**：脚本头部预声明 `gr()`（GR 是 Plots 默认推荐、无 Python 依赖）。**不要**使用 `zcolor + Vector{Shape}` 组合——该组合在 GR 后端颜色映射不可靠。本计划采用"每单元一个 `plot!(Shape(...); color=...)` 调用"的稳健方案。

---

## Chunk 1: 脚本主体

### Task 1: 脚本骨架（文件头 + Option 配置 + 空 main）

**Files:**
- Create: `example/jellyroll_stress_displacement.jl`

**说明：** 先建立文件、模块加载、Plots 后端、注释头、与 testexample 对齐的 Option 配置（唯一差异 `czm_enabled=false`）。本 Task 的 `main()` 只含参数设置段，后续 Task 用 **明确的 anchor 字符串** 在 `# ===...# 5. ...` 段之前插入新代码（**不是替换**）。

- [ ] **Step 1: 创建脚本，写入文件头、模块加载、helper 函数、空 main**

写入 `example/jellyroll_stress_displacement.jl`：

```julia
"""
测试案例：Jellyroll 电池单次放电二维应力场/位移场后处理

功能：
- 多 SPMe 并行 + 二维分布式热模型（与 testexample 同工况，但关闭 CZM）
- 在用户指定的时间节点调用 thermal_diffusion_stress_2D 计算宏观应力/位移场
- 输出一张大图：行 = 场分量（σ_xx/σ_yy/σ_xy/σ_vm/|U|），列 = 时间节点

⚠️ 尺度说明：本脚本输出的是 **极片/电极尺度（coating-scale, 二维平面应力）** 的场，
   由全叠合厚度加权有效模量 E_eff（~5e8 Pa 量级）计算；
   与颗粒尺度 Calstressdisp（颗粒 E ~1e10 Pa）不同，二者不可混用。

日期：2026-06-30
"""

using Printf
ENV["PLOTS_DEFAULT_BACKEND"] = get(ENV, "JUBAT_PLOTS_BACKEND", "gr")
using Plots
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 用户可配置参数 ============================================================
# 指定要绘图的物理时间节点（秒）；脚本会自动在 result 时间序列里找最接近的索引
const PLOT_TIMES_S = [600, 1800, 3600]
# ===========================================================================

# --------------------------------------------------------------------------
# Helper: 提取每个 Q4 单元的四角坐标（首尾闭合，5 个点）
# --------------------------------------------------------------------------
function q4_element_polygons(mesh)
    ne = size(mesh.element, 1)
    xs = Vector{Vector{Float64}}(undef, ne)
    ys = Vector{Vector{Float64}}(undef, ne)
    @inbounds for e in 1:ne
        n1, n2, n3, n4 = mesh.element[e, 1], mesh.element[e, 2], mesh.element[e, 3], mesh.element[e, 4]
        xs[e] = [mesh.node[n1,1], mesh.node[n2,1], mesh.node[n3,1], mesh.node[n4,1], mesh.node[n1,1]]
        ys[e] = [mesh.node[n1,2], mesh.node[n2,2], mesh.node[n3,2], mesh.node[n4,2], mesh.node[n1,2]]
    end
    return xs, ys
end

# --------------------------------------------------------------------------
# Helper: 在子图 plt 上绘制单元常数场 field_elem（长度 ne）的填色云图。
# 采用"每个单元单独 plot!(Shape; color=...)"——GR 后端最稳健的方案。
# deform_xy = (U_x, U_y)：若提供，按 def_scale 放大叠加变形网格。
# 注：xs[e]/ys[e] 含 5 个点（第 5 = 第 1 闭合），所以 deform 循环用 mod1(k,4)
#     把第 5 个点也映射到 node 1，保持闭合。
# --------------------------------------------------------------------------
function plot_q4_field!(plt, mesh, field_elem, xs, ys; title="", cmap=:RdYlBu_r,
                        clim=nothing, deform_xy=nothing, def_scale=0.0)
    fmin, fmax = (clim === nothing) ? (minimum(field_elem), maximum(field_elem)) : clim
    if fmax - fmin < 1e-30
        fmax = fmin + 1.0  # 避免除零
    end
    cmap_grad = cgrad(cmap)
    for e in 1:length(field_elem)
        xe = copy(xs[e]); ye = copy(ys[e])
        if deform_xy !== nothing
            U_x, U_y = deform_xy
            for k in 1:5
                n = mesh.element[e, mod1(k, 4)]
                xe[k] += def_scale * U_x[n]
                ye[k] += def_scale * U_y[n]
            end
        end
        idx = (field_elem[e] - fmin) / (fmax - fmin)
        col = cmap_grad[clamp(idx, 0.0, 1.0)]
        plot!(plt, Plots.Shape(xe, ye); linecolor=:match, linewidth=0.1,
              color=col, label=false, colorbar_entry=false)
    end
    title!(plt, title); xlabel!(plt, "x [m]"); ylabel!(plt, "y [m]")
end

function main()
    println("="^80)
    println("Jellyroll 单次放电 二维应力/位移场后处理")
    println("="^80)

    # ====================================================================
    # 1. 参数设置（同 testexample，唯一差异：czm_enabled = false）
    # ====================================================================
    println("\n[1/5] 参数设置...")

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = JuBat.Option()
    Crates = 1.0
    i = 5 * Crates
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"   # 不启用颗粒尺度应力耦合

    opt.time = [0.0, 3600.0]
    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = false        # ← 关键差异：关闭 CZM，仅传统固体力学

    println("OK: 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f s\n", opt.time[end])
    println("  CZM: 关闭（仅传统固体力学）")

    # === ANCHOR_END_PARAMS ===
end

main()
```

- [ ] **Step 2: smoke-run 检查脚本可加载、main 可执行（仅参数段）**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia example/jellyroll_stress_displacement.jl`
Expected: 打印 "OK: 参数设置完成" 后正常退出。若报 `Package Plots not found`，先 `julia -e 'using Pkg; Pkg.add("Plots")'`；若 `Plots` 因 `GR` 编译失败，可设 `$env:JUBAT_PLOTS_BACKEND="pyplot"` 改用 pyplot（需 `Pkg.add("PyPlot")`）。

- [ ] **Step 3: Commit**

```bash
git add example/jellyroll_stress_displacement.jl
git commit -m "feat(example): 新增 Jellyroll 2D 应力脚本骨架（含 Option/helper）"
```

---

### Task 2: 网格创建 + 求解 + 时间节点选择

**Files:**
- Modify: `example/jellyroll_stress_displacement.jl`

**说明：** 在 `main()` 中、`# === ANCHOR_END_PARAMS ===` 这一行**之前**插入"网格 + 求解 + 时间节点选择"段。本 Task 通过查找 anchor 字符串 `# === ANCHOR_END_PARAMS ===` 来定位插入点。

- [ ] **Step 1: 在 anchor `# === ANCHOR_END_PARAMS ===` 之前插入网格创建与求解代码**

将 `# === ANCHOR_END_PARAMS ===` 这一行（不含其后的 `end`）替换为以下完整代码块（块末尾仍保留原 anchor 行，供 Task 3 继续在其之前插入）：

```julia
    # ====================================================================
    # 2. 创建案例与 Jellyroll 网格
    # ====================================================================
    println("\n[2/5] 创建案例与 Jellyroll 网格...")

    case = JuBat.SetCase(param_dim, opt)

    n_theta = 360
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # czm_enabled=false：不创建 case.czm_mesh

    ne   = size(mesh_th.element, 1)
    nlen = mesh_th.nlen
    println("OK: 网格创建完成")
    @printf("  n_theta=%d, ne=%d, nlen=%d\n", n_theta, ne, nlen)

    # 前置字段核查（避免后续 helper 误用字段名）
    @assert hasproperty(mesh_th, :node)    "thermal2D mesh 缺少 .node 字段，实际字段: $(propertynames(mesh_th))"
    @assert hasproperty(mesh_th, :element) "thermal2D mesh 缺少 .element 字段"
    @assert size(mesh_th.node, 2) >= 2     "mesh_th.node 列数 < 2"
    @assert size(mesh_th.element, 2) == 4  "Q4 单元应有 4 列"

    # ====================================================================
    # 3. 求解
    # ====================================================================
    println("\n[3/5] 运行 SPMe+热 求解器（无 CZM）...")

    t_wall_start = time_ns()
    result = JuBat.Solve(case)
    t_wall_s = (time_ns() - t_wall_start) * 1e-9
    println("OK: 求解完成")
    @printf("  wall-clock: %.2f s\n", t_wall_s)

    # ====================================================================
    # 4. 时间节点选择
    # ====================================================================
    println("\n[4/5] 选择时间节点...")

    t_s = result["time [s]"]
    @assert haskey(result, "thermal2D temperature at nodes [K]") "result 缺少 'thermal2D temperature at nodes [K]'，请确认 opt.thermalmodel == \"distributed2D\""
    @assert haskey(result, "thermal2D element soc_n") "result 缺少 'thermal2D element soc_n'"
    @assert haskey(result, "thermal2D element soc_p") "result 缺少 'thermal2D element soc_p'"

    ti_list = Int[]
    for t_target in PLOT_TIMES_S
        if t_target < t_s[1] || t_target > t_s[end]
            @warn "PLOT_TIMES_S 中 $(t_target) s 超出 [$(t_s[1]), $(t_s[end])] 范围，钳到端点" maxlog=1
        end
        ti = argmin(abs.(t_s .- t_target))
        push!(ti_list, ti)
    end
    unique!(ti_list)
    sort!(ti_list)

    println("OK: 选定时间节点")
    for ti in ti_list
        @printf("  ti=%d  t=%.1f s\n", ti, t_s[ti])
    end

    # === ANCHOR_END_PARAMS ===
```

> 插入完成后，`# === ANCHOR_END_PARAMS ===` 仍紧贴 `main()` 的 `end` 之前。

- [ ] **Step 2: smoke-run 检查"求解 + 节点选择"跑通**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia example/jellyroll_stress_displacement.jl`
Expected: 完整跑完一次 1C 放电，打印 3 行 `ti=...  t=... s`（分别接近 600/1800/3600）。
**预计 wall-clock：1–5 分钟（nθ=360，per_element_spme=true）。若需快速验证可临时把 `n_theta=360` 改成 120。**

- [ ] **Step 3: Commit**

```bash
git add example/jellyroll_stress_displacement.jl
git commit -m "feat(example): 网格创建+求解+时间节点选择（含字段前置核查）"
```

---

### Task 3: 全部节点循环 + sanity check + 绘图

**Files:**
- Modify: `example/jellyroll_stress_displacement.jl`

**说明：** 在 anchor `# === ANCHOR_END_PARAMS ===` 之前插入"场重建循环 + 数量级 sanity check + 大图绘制"段。本 Task 一次性写完，不依赖中间态替换。

- [ ] **Step 1: 在 `# === ANCHOR_END_PARAMS ===` 之前插入场重建与绘图代码**

将 `# === ANCHOR_END_PARAMS ===` 这一行替换为以下完整代码块（块末尾保留 anchor 行；本 task 之后 anchor 不再使用，仅留作标记）：

```julia
    # ====================================================================
    # 5. 在所有选定时间节点重建场 + 数量级 sanity check
    # ====================================================================
    println("\n[5/5] 重建二维场（$(length(ti_list)) 个时间节点）...")

    L_ref = case.param.scale.L

    # NamedTuple 向量（保留类型信息）
    FieldPack = NamedTuple{(:t, :σ_xx, :σ_yy, :σ_xy, :σ_vm, :U_x, :U_y),
                           Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64},
                                 Vector{Float64}, Vector{Float64}, Vector{Float64}}}
    fields_per_ti = FieldPack[]

    σ_vm_global_max = 0.0
    Umag_global_max = 0.0
    for ti in ti_list
        variables_ti = Dict{String, Union{Array{Float64},Float64}}(
            "T_nodes"                 => result["thermal2D temperature at nodes [K]"][:, ti],
            "thermal2D element soc_n" => result["thermal2D element soc_n"][:, ti],
            "thermal2D element soc_p" => result["thermal2D element soc_p"][:, ti],
        )
        @assert length(variables_ti["T_nodes"]) == nlen
        @assert length(variables_ti["thermal2D element soc_n"]) == ne

        vars_out = JuBat.thermal_diffusion_stress_2D(case, variables_ti)
        σ_xx = vars_out["diffusion stress xx"]
        σ_yy = vars_out["diffusion stress yy"]
        σ_xy = vars_out["diffusion stress xy"]
        σ_vm = vars_out["diffusion stress vonMises"]
        U_x  = vars_out["displacement x"]
        U_y  = vars_out["displacement y"]
        @assert length(σ_vm) == ne
        @assert length(U_x) == nlen

        Umag = sqrt.(U_x.^2 .+ U_y.^2)
        σ_vm_global_max = max(σ_vm_global_max, maximum(σ_vm))
        Umag_global_max = max(Umag_global_max, maximum(Umag))

        push!(fields_per_ti, (t=t_s[ti], σ_xx=σ_xx, σ_yy=σ_yy, σ_xy=σ_xy,
                              σ_vm=σ_vm, U_x=U_x, U_y=U_y))

        @printf("  t=%7.1f s: σ_vm_max=%.3e Pa (%.2f MPa), |U|_max=%.3e m (%.3f µm)\n",
                t_s[ti], maximum(σ_vm), maximum(σ_vm)*1e-6,
                maximum(Umag), maximum(Umag)*1e6)
    end

    # === sanity check：量级 + 量纲 ===
    ti_last = ti_list[end]
    soc_n_last = result["thermal2D element soc_n"][:, ti_last]
    @printf("  sanity: param.NE.cs0(normalized)=%.4f  (期望 [0,1], 典型~0.9)\n", case.param.NE.cs0)
    @printf("  sanity: param_dim.NE.cs0=%.1f [mol/m3]  (物理浓度, 仅对比)\n", case.param_dim.NE.cs0)
    @printf("  sanity: soc_n range=[%.4f, %.4f]  (期望 [0,1])\n",
            minimum(soc_n_last), maximum(soc_n_last))
    @printf("  sanity: σ_vm_global_max=%.3e Pa（极片尺度，期望 1e7–1e9；若 ~1e10 说明误用颗粒 E）\n",
            σ_vm_global_max)

    # 若所有节点位移为零，疑似力学求解失败
    @assert Umag_global_max > 0 "所有节点位移均为零，疑似力学求解失败（K_mech\\F_mech 异常），请检查 src/Mechanical.jl:289 catch 分支是否触发"

    # 自适应变形放大系数：使最大变形 ≈ 5% L_ref
    DEF_SCALE = 0.05 * L_ref / Umag_global_max
    @printf("  DEF_SCALE=%.3e（自适应，使 |U|_max×DEF_SCALE ≈ 5%% L_ref）\n", DEF_SCALE)

    # ====================================================================
    # 6. 绘图：一张大图，行=场分量，列=时间节点
    # ====================================================================
    println("\n[绘图]")

    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)

    nT = length(fields_per_ti)
    row_labels = ["σ_xx [MPa]", "σ_yy [MPa]", "σ_xy [MPa]", "σ_vm [MPa]", "|U| [μm]"]
    nrow = length(row_labels)
    subplots = Plot[]

    xs0, ys0 = q4_element_polygons(mesh_th)

    for col in 1:nT
        pack = fields_per_ti[col]
        Umag = sqrt.(pack.U_x.^2 .+ pack.U_y.^2)
        # 单位换算到 MPa / μm
        col_fields = [pack.σ_xx*1e-6, pack.σ_yy*1e-6, pack.σ_xy*1e-6,
                      pack.σ_vm*1e-6, Umag*1e6]
        for row in 1:nrow
            p = plot(legend=false)
            if row < 5
                clim = (minimum(col_fields[row]), maximum(col_fields[row]))
                plot_q4_field!(p, mesh_th, col_fields[row], xs0, ys0;
                    title=(@sprintf("t=%.0f s, %s", pack.t, row_labels[row])),
                    cmap=:RdYlBu_r, clim=clim)
            else
                clim = (minimum(col_fields[row]), maximum(col_fields[row]))
                plot_q4_field!(p, mesh_th, col_fields[row], xs0, ys0;
                    title=(@sprintf("t=%.0f s, %s (def ×%.0f)", pack.t, row_labels[row], DEF_SCALE)),
                    cmap=:viridis, clim=clim,
                    deform_xy=(pack.U_x, pack.U_y), def_scale=DEF_SCALE)
            end
            push!(subplots, p)
        end
    end

    p_combined = plot(subplots..., layout=(nrow, nT),
                      size=(450*nT, 350*nrow), left_margin=5mm, bottom_margin=5mm)
    out_path = joinpath(output_dir, "jellyroll_stress_displacement.png")
    savefig(p_combined, out_path)
    @printf("  大图已保存: %s\n", out_path)

    println("\n" * "="^80)
    @printf("完成：σ_vm 全局峰值 = %.2f MPa（极片/电极尺度）\n", σ_vm_global_max*1e-6)
    println("="^80)

    # === ANCHOR_END_PARAMS ===
```

- [ ] **Step 2: smoke-run 端到端**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia example/jellyroll_stress_displacement.jl`
Expected（按顺序）：
1. 求解完成（1–5 分钟）。
2. 3 行 `t=... σ_vm_max=... |U|_max=...`，**σ_vm_max 应在 1e7–1e9 Pa**（1e7–1e3 MPa），**|U|_max 应在 1e-7–1e-4 m**（0.1–100 µm）。
3. 4 行 sanity：`param.NE.cs0(normalized)` 在 [0,1]（典型 ~0.9）；`param_dim.NE.cs0` ~29866；`soc_n range` 在 [0,1]；`σ_vm_global_max` 在 1e7–1e9 Pa。
4. `DEF_SCALE=...` 一行。
5. 生成 `output/jellyroll_stress_displacement.png`。

**异常诊断**：
- 若 `σ_vm_max` ≈ 1e10 Pa → 误用颗粒模量 `PE.E`，回头核查 `thermal_diffusion_stress_2D` 内部用的是否 `compute_effective_coating_modulus` 返回值。
- 若 `soc_n range` 不在 [0,1] → 排查归一化路径（NormaliseParam）。
- 若图全空白 → 已用每单元单独 `plot!(Shape; color=...)`，理论上不会发生；若仍空白，临时设 `$env:JUBAT_PLOTS_BACKEND="pyplot"` 试一次。
- 若 assert "所有节点位移均为零" 触发 → 排查 `K_mech\F_mech` 在 src/Mechanical.jl:289 是否被 catch。

- [ ] **Step 3: 手工可视化检查**

打开 `output/jellyroll_stress_displacement.png`：
- 5 行 × N 列子图布局正确，每列标题含 `t=`。
- 应力子图（前 4 行）非全白/非全单色，有色彩梯度。
- 末行（位移）应能看到**轻微变形的螺旋网格**（最大变形 ≈ 5% 特征长度）。

- [ ] **Step 4: Commit**

```bash
git add example/jellyroll_stress_displacement.jl
git commit -m "feat(example): 全节点循环+sanity check+5×N 大图（极片尺度）"
```

---

### Task 4: 文档收尾

**Files:**
- Modify: `CLAUDE.md` 第 10 节示例文件表

**说明：** 在仓库示例索引中登记新脚本，便于他人发现。

- [ ] **Step 1: 在 CLAUDE.md §10 示例文件表中追加一行**

在 `| example/testexample.jl | 全耦合仿真 |` 之后追加：

```markdown
| `example/jellyroll_stress_displacement.jl` | 二维应力/位移场后处理（无 CZM） |
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs(claude.md): 登记 jellyroll_stress_displacement 示例"
```

---

## 验收清单（实现完成后逐项打勾）

- [ ] 脚本 `julia example/jellyroll_stress_displacement.jl` 一键跑通，无报错、无未捕获 `@warn` 之外的异常。
- [ ] 生成 `output/jellyroll_stress_displacement.png`，含 5×N 子图。
- [ ] `σ_vm_max` 量级在 1e7–1e9 Pa（极片尺度，**非**颗粒尺度 1e10）。
- [ ] `|U|_max` 量级在 1e-7–1e-4 m。
- [ ] `soc_n range` 落在 [0, 1]。
- [ ] `src/` 下文件零修改（`git diff main -- src/` 为空）。
