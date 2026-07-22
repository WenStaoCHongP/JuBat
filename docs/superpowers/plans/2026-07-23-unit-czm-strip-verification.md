# 单元级 CZM 条带验证 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现平直 8 层条带网格生成器（8×Q4 + 4×COH2D4，复用生产 `create_czm_mesh`），并交付两个独立验证脚本：位移 BC 验证双线性律、合成本征应变验证弹性段分离位移。

**Architecture:** `create_unit_czm_strip` 组装 `CzmSubmesh` + 覆盖整条带的 1 单元哑热网格，调用生产 `create_czm_mesh` 并硬断言拓扑。两脚本各自用增量 Newton + `apply_bc_czm(bc_dofs, bc_vals)`，不走环几何 `solve_czm_step`。脚本 2 用冷却（ΔT≤0）使固定端产生拉伸张开，与 1D 串联杆+弹簧解析解比对。

**Tech Stack:** Julia、JuBat（`ChooseCell("Jellyroll")`、`SetCase`、`compute_czm_params_per_interface`、`assemble_coupled_system` / `_full`、`update_damage_per_interface`）、`Test` 标准库。

**Spec:** `docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md`

## Global Constraints

- 参数：`ChooseCell("Jellyroll")` → `SetCase`；界面参数仅经 `compute_czm_params_per_interface`
- 层序自下而上：`[:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]`（与 `build_czm_submesh` 一致）
- 网格路径：方案 C——复用 `create_czm_mesh`，禁止手工平行实现节点复制
- 不修改生产 `create_czm_mesh` / `bilinear_*` 核心逻辑（仅新增文件与 `JuBat.jl` 接线）
- 求解：脚本内 Newton；禁止依赖 `identify_bc_nodes_czm`
- 脚本 2 本征应变契约：`ε₀[e] = α_eff·dT[e] + β_n·Δsoc_n[e] + β_p·Δsoc_p[e]`（生产 `assemble_thermal_chemical_load`）
- 脚本 2 加载：冷却为主（`ΔT ≤ 0`），保证固定端产生**拉伸**张开；幅值保持 `δ < δ_0`（弹性段）
- 位移↔分离：`u = δ̃ / Λ`（`Λ = params.Λ`）

## File Structure

| 文件 | 职责 |
|------|------|
| `src/CzmUnitMesh.jl` | `create_unit_czm_strip` + 哑热网格 + 硬断言；返回 `(czm_mesh, meta)` |
| `src/JuBat.jl` | `include` + `export create_unit_czm_strip` |
| `test/unit_czm_newton.jl` | 共享：增量 Newton、BC 组装、解析双线性辅助（被两脚本 include） |
| `test/unit_czm_strip_mesh.jl` | 网格拓扑 TDD / 回归 |
| `test/unit_czm_bilinear.jl` | 脚本 1：Mode I / 卸载 / Mode II |
| `test/unit_czm_eigenstrain.jl` | 脚本 2：合成斜坡 + 1D 解析 |
| 计划摘要 md | 勾选实施状态 |

---

### Task 1: 条带网格生成器（TDD）

**Files:**
- Create: `src/CzmUnitMesh.jl`
- Create: `test/unit_czm_strip_mesh.jl`
- Modify: `src/JuBat.jl`（`include` 在 `czm.jl` 之后；`export create_unit_czm_strip`）

**Interfaces:**
- Consumes: `CzmSubmesh`, `Mesh`, `GaussPoint`/`GetGS`, `create_czm_mesh`
- Produces:
  ```julia
  create_unit_czm_strip(param; width=nothing, y0=1.0, gsorder=2)
      -> (czm_mesh::CohesiveMesh, meta::NamedTuple)
  # meta 字段：
  #   y_interfaces::Vector{Float64}   # 长度 9，底→顶层界
  #   bottom_nodes::Vector{Int}       # 长度 2
  #   top_nodes::Vector{Int}          # 长度 2（原始顶边；复制后外层可能用副本，BC 用最终 mesh 上的顶边节点）
  #   top_nodes_after_czm::Vector{Int} # create_czm_mesh 后顶层 Q4 的上边 2 节点（施 BC 用这个）
  #   cohesive_ids::Vector{Int}       # 长度 4
  #   interface_types::Vector{Symbol} # 长度 4
  #   layer_materials::Vector{Symbol} # 长度 8
  #   width::Float64
  #   heights::Vector{Float64}        # 长度 8 层厚
  ```

- [ ] **Step 1: 写失败的网格测试**

创建 `test/unit_czm_strip_mesh.jl`：

```julia
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test

@testset "create_unit_czm_strip topology" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    param = case.param

    czm_mesh, meta = JuBat.create_unit_czm_strip(param; y0=1.0, gsorder=2)

    @test size(czm_mesh.bulk_element, 1) == 8
    @test czm_mesh.n_cohesive == 4
    types = [e.interface_type for e in czm_mesh.cohesive_elements]
    @test count(==(:PE_PCC), types) == 2
    @test count(==(:NE_NCC), types) == 2
    @test length(meta.layer_materials) == 8
    @test meta.layer_materials == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    @test length(meta.bottom_nodes) == 2
    @test length(meta.top_nodes_after_czm) == 2
    # 底边 y 最小、顶边 y 最大
    @test all(czm_mesh.node[n, 2] ≈ meta.y_interfaces[1] for n in meta.bottom_nodes)
    @test maximum(czm_mesh.node[meta.top_nodes_after_czm, 2]) ≈ meta.y_interfaces[end] atol=1e-12
end
```

- [ ] **Step 2: 运行确认失败**

Run: `julia --project=. test/unit_czm_strip_mesh.jl`  
Expected: 失败（`create_unit_czm_strip` 未定义）

- [ ] **Step 3: 实现 `src/CzmUnitMesh.jl`**

```julia
"""
    create_unit_czm_strip(param; width=nothing, y0=1.0, gsorder=2)

平直 8 层单元条带：8×Q4 + 经 create_czm_mesh 得到 4×COH2D4。
层序自下而上 PE→PCC→PE→SP→NE→NCC→NE→SP。底边 y=y0>0，节点 id 沿 x 递增。
"""
function create_unit_czm_strip(param; width=nothing, y0::Float64=1.0, gsorder::Int=2)
    layer_materials = [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    heights = [
        param.PE.thickness, param.PCC.thickness,
        param.PE.thickness, param.SP.thickness,
        param.NE.thickness, param.NCC.thickness,
        param.NE.thickness, param.SP.thickness,
    ]
    H = sum(heights)
    W = width === nothing ? H : Float64(width)
    y_interfaces = Vector{Float64}(undef, 9)
    y_interfaces[1] = y0
    for i in 1:8
        y_interfaces[i + 1] = y_interfaces[i] + heights[i]
    end

    # 节点：行优先，每行左→右（id 沿 x 递增）
    # 行 r=0..8 对应 y_interfaces[r+1]；列 c=0..1 对应 x=0,W
    nnode = 18
    node = zeros(Float64, nnode, 2)
    for r in 0:8
        for c in 0:1
            idx = r * 2 + c + 1
            node[idx, 1] = c == 0 ? 0.0 : W
            node[idx, 2] = y_interfaces[r + 1]
        end
    end

    # Q4：下层左、上层左、上层右、下层右（逆时针，inner=bottom）
    element = zeros(Int64, 8, 4)
    for e in 1:8
        bl = (e - 1) * 2 + 1   # bottom-left
        br = bl + 1
        tl = e * 2 + 1         # top-left
        tr = tl + 1
        element[e, :] = [bl, tl, tr, br]
    end

    gs = GetGS(element, node, gsorder, "Q4")
    bulk = Mesh("Q4", 2, node, nnode, element, gs)
    material_type = copy(layer_materials)
    winding_turn = ones(Int, 8)
    thermal_elem_map = ones(Int, 8)
    submesh = CzmSubmesh(bulk, material_type, winding_turn, thermal_elem_map)

    # 哑热网格：单 Q4 覆盖条带 bbox（供 build_thermal_to_czm_interp）
    pad = 1e-6 * max(W, H)
    tn = zeros(Float64, 4, 2)
    tn[1, :] = [-pad, y0 - pad]
    tn[2, :] = [W + pad, y0 - pad]
    tn[3, :] = [W + pad, y_interfaces[end] + pad]
    tn[4, :] = [-pad, y_interfaces[end] + pad]
    te = reshape(Int64[1, 2, 3, 4], 1, 4)
    tgs = GetGS(te, tn, gsorder, "Q4")
    dummy_thermal = Mesh("Q4", 2, tn, 4, te, tgs)

    czm_mesh = create_czm_mesh(submesh, dummy_thermal, param)

    # ---- 硬断言（方案 C）----
    @assert czm_mesh.n_cohesive == 4 "expected 4 cohesive, got $(czm_mesh.n_cohesive)"
    types = [e.interface_type for e in czm_mesh.cohesive_elements]
    @assert count(==(:PE_PCC), types) == 2 "PE_PCC count"
    @assert count(==(:NE_NCC), types) == 2 "NE_NCC count"
    for coh in czm_mesh.cohesive_elements
        n_lo, n_hi, n_hi_c, n_lo_c = coh.nodes
        @assert czm_mesh.node[n_lo, :] ≈ czm_mesh.node[n_lo_c, :] atol=1e-12
        @assert czm_mesh.node[n_hi, :] ≈ czm_mesh.node[n_hi_c, :] atol=1e-12
        @assert length(unique(coh.nodes)) == 4
        # 外层 bulk 共边应为副本
        e_out = coh.host_outer_elem
        outer_nodes = Set(czm_mesh.bulk_element[e_out, :])
        @assert n_lo_c in outer_nodes && n_hi_c in outer_nodes
        @assert !(n_lo in outer_nodes) && !(n_hi in outer_nodes)
    end

    bottom_nodes = [1, 2]
    # 顶层（层 8）在 create_czm_mesh 后：外层若被重写则取 bulk_element[8,:] 上边
    # 层 8 上边原始节点为 17,18；复制后读最终连接
    top_row_orig = [17, 18]
    e8 = czm_mesh.bulk_element[8, :]
    # 上边 = y 最大的两个节点
    ys = [czm_mesh.node[n, 2] for n in e8]
    ytop = maximum(ys)
    top_nodes_after = sort([n for n in e8 if abs(czm_mesh.node[n, 2] - ytop) < 1e-14])

    meta = (
        y_interfaces = y_interfaces,
        bottom_nodes = bottom_nodes,
        top_nodes = top_row_orig,
        top_nodes_after_czm = top_nodes_after,
        cohesive_ids = collect(1:4),
        interface_types = types,
        layer_materials = layer_materials,
        width = W,
        heights = heights,
    )
    return czm_mesh, meta
end
```

在 `src/JuBat.jl`：

```julia
# 紧接 include("czm.jl") 之后：
include("CzmUnitMesh.jl")
```

```julia
# export 区增加：
export create_unit_czm_strip
```

- [ ] **Step 4: 运行网格测试通过**

Run: `julia --project=. test/unit_czm_strip_mesh.jl`  
Expected: 全部 `@testset` PASS。若 `n_cohesive != 4`，检查层序/共边材料对与质心径向排序（`y0` 必须 >0）。

- [ ] **Step 5: Commit**

```bash
git add src/CzmUnitMesh.jl src/JuBat.jl test/unit_czm_strip_mesh.jl
git commit -m "feat(czm): add unit strip mesh generator (8 Q4 + 4 COH2D4)"
```

---

### Task 2: 共享增量 Newton 辅助

**Files:**
- Create: `test/unit_czm_newton.jl`

**Interfaces:**
- Consumes: `assemble_coupled_system`, `assemble_coupled_system_full`, `apply_bc_czm`, `update_damage_per_interface`, `clone_damage_states`（均在 JuBat 模块内）
- Produces:
  ```julia
  function unit_czm_newton_step!(czm_mesh, u, param_cache;
      bc_dofs, bc_vals, F_ext=nothing, F_thermo_chem=nothing,
      max_iter=50, tol=1e-8) -> (u, separations, tractions, converged, R_norm)

  function analytic_bilinear_T(δ, K, δ_0, δ_c, D_hist=0.0) -> T
  # 卸载：若 δ < δ_max_hist，T = (1-D)*K*δ
  ```

- [ ] **Step 1: 实现 `test/unit_czm_newton.jl`**

```julia
using LinearAlgebra

function unit_czm_newton_step!(czm_mesh, u::Vector{Float64}, param_cache;
        bc_dofs::Vector{Int}, bc_vals::Vector{Float64},
        F_ext::Union{Nothing,Vector{Float64}}=nothing,
        F_thermo_chem::Union{Nothing,Vector{Float64}}=nothing,
        max_iter::Int=80, tol::Float64=1e-8)

    ndof = length(u)
    F_e = F_ext === nothing ? zeros(ndof) : F_ext
    F_tc = F_thermo_chem === nothing ? zeros(ndof) : F_thermo_chem
    damage_states = czm_mesh.damage_states
    # 将位移 BC 硬写入初值，避免首步残差主导
    for (dof, val) in zip(bc_dofs, bc_vals)
        u[dof] = val
    end

    separations = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    tractions = Vector{Tuple{Float64,Float64}}(undef, czm_mesh.n_cohesive)
    converged = false
    R_norm = Inf

    for iter in 1:max_iter
        K, f_int, separations, tractions = JuBat.assemble_coupled_system(
            czm_mesh, u, param_cache; damage_states=damage_states)
        R = F_e + F_tc - f_int
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        R_norm = norm(R)
        if R_norm < tol
            converged = true
            damage_states = JuBat.update_damage_per_interface(
                czm_mesh, damage_states, separations, param_cache)
            czm_mesh.damage_states = damage_states
            break
        end
        K_bc, R_bc = JuBat.apply_bc_czm(K, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
        Δu = K_bc \ R_bc
        any(!isfinite, Δu) && break
        u .+= Δu
        for (dof, val) in zip(bc_dofs, bc_vals)
            u[dof] = val
        end
    end
    return u, separations, tractions, converged, R_norm
end

"""解析双线性（单调加载，无历史损伤时）。δ,δ_0,δ_c 同空间。"""
function analytic_bilinear_T(δ::Float64, K::Float64, δ_0::Float64, δ_c::Float64)
    δ <= 0 && return K * δ          # 压缩罚（与本构一致的简化）
    δ <= δ_0 && return K * δ
    δ >= δ_c && return 0.0
    # 软化：T = σ_max * (δ_c - δ)/(δ_c - δ_0)，σ_max = K*δ_0
    return (K * δ_0) * (δ_c - δ) / (δ_c - δ_0)
end
```

- [ ] **Step 2: Commit**

```bash
git add test/unit_czm_newton.jl
git commit -m "test(czm): shared unit-strip Newton helper"
```

---

### Task 3: 脚本 1 — 双线性位移验证

**Files:**
- Create: `test/unit_czm_bilinear.jl`
- Test: 同文件 `@testset`

**Interfaces:**
- Consumes: `create_unit_czm_strip`, `unit_czm_newton_step!`, `analytic_bilinear_T`, `compute_czm_params_per_interface`

- [ ] **Step 1: 写 `test/unit_czm_bilinear.jl`（三工况）**

```julia
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "unit_czm_newton.jl"))

function _setup_strip()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    czm_mesh, meta = JuBat.create_unit_czm_strip(case.param; y0=1.0)
    cache = JuBat.compute_czm_params_per_interface(case)
    return case, czm_mesh, meta, cache
end

function _bc_bottom_top(meta, czm_mesh; u_top_x=0.0, u_top_y=0.0)
    bc_dofs = Int[]; bc_vals = Float64[]
    for n in meta.bottom_nodes
        push!(bc_dofs, 2n-1); push!(bc_vals, 0.0)
        push!(bc_dofs, 2n);   push!(bc_vals, 0.0)
    end
    for n in meta.top_nodes_after_czm
        push!(bc_dofs, 2n-1); push!(bc_vals, u_top_x)
        push!(bc_dofs, 2n);   push!(bc_vals, u_top_y)
    end
    return bc_dofs, bc_vals
end

@testset "Mode I monotonic bilinear" begin
    case, czm_mesh, meta, cache = _setup_strip()
    # 重置损伤
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    # 位移空间目标分离（取第一个界面参数；两 PE_PCC 同参）
    δ_targets = [0.5 * pe.δ_0_n, 0.5 * (pe.δ_0_n + pe.δ_c_n), 1.2 * pe.δ_c_n]
    u = zeros(2 * czm_mesh.nnode)
    for δ_tgt in δ_targets
        u_y = δ_tgt / Λ   # 顶边相对底边张开；多界面串联时总张开分配到各界面+体变形
        # 为使界面进入目标段：用足够大的总位移；弹性段用小位移
        # 策略：逐步增大 u_y，记录平均 δ_n，与解析 T(δ) 比
        nothing
    end
    # ---- 实际加载：从 0 扫到 1.5*δ_c/Λ 的总顶部位移，取若干采样点 ----
    n_steps = 40
    u_max = 1.5 * pe.δ_c_n / Λ
    # 注意：4 个界面串联 + 体弹性，总位移 ≠ 单界面 δ。验收改为：
    # 对每个 cohesive，用其自身 (δ_n, T_n) 与 analytic_bilinear_T(δ_n, ...) 自洽。
    δ_hist = Float64[]; T_hist = Float64[]; D_hist = Float64[]
    for s in 1:n_steps
        u_y = u_max * s / n_steps
        bc_dofs, bc_vals = _bc_bottom_top(meta, czm_mesh; u_top_y=u_y)
        u, seps, tracts, ok, Rn = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        # 取第一个 PE_PCC 界面
        i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
        δn, _ = seps[i]; Tn, _ = tracts[i]
        push!(δ_hist, δn); push!(T_hist, Tn)
        push!(D_hist, czm_mesh.damage_states[i].D)
        T_ana = analytic_bilinear_T(δn, pe.K_n, pe.δ_0_n, pe.δ_c_n)
        if δn < pe.δ_0_n * 0.98
            @test Tn ≈ T_ana rtol=1e-2 atol=1e-8
        elseif δn < pe.δ_c_n
            @test Tn ≈ T_ana rtol=5e-2 atol=1e-6
        else
            @test abs(Tn) < 1e-3 * max(pe.K_n * pe.δ_0_n, 1.0)
            @test czm_mesh.damage_states[i].D > 0.99
        end
    end
    @test maximum(δ_hist) > pe.δ_0_n
    @test maximum(D_hist) > 0.99
end

@testset "Mode I unload reload" begin
    case, czm_mesh, meta, cache = _setup_strip()
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    u = zeros(2 * czm_mesh.nnode)
    # 加载到软化中点附近的总位移（经验：数倍 δ_0/Λ）
    u_peak = 0.6 * pe.δ_c_n / Λ
    for s in 1:20
        u_y = u_peak * s / 20
        bc_dofs, bc_vals = _bc_bottom_top(meta, czm_mesh; u_top_y=u_y)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
    end
    i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
    D_peak = czm_mesh.damage_states[i].D
    δ_peak = seps[i][1]
    @test 0.05 < D_peak < 0.99
    # 卸载
    for s in 1:20
        u_y = u_peak * (1 - s / 20)
        bc_dofs, bc_vals = _bc_bottom_top(meta, czm_mesh; u_top_y=u_y)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        @test czm_mesh.damage_states[i].D >= D_peak - 1e-12  # 单调
    end
    D_unload = czm_mesh.damage_states[i].D
    @test D_unload ≈ D_peak rtol=1e-8
    # 再加载未超过历史 δ_max 时 D 不变
    for s in 1:10
        u_y = 0.5 * u_peak * s / 10
        bc_dofs, bc_vals = _bc_bottom_top(meta, czm_mesh; u_top_y=u_y)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        @test czm_mesh.damage_states[i].D ≈ D_unload rtol=1e-8
    end
end

@testset "Mode II tangential" begin
    case, czm_mesh, meta, cache = _setup_strip()
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    u = zeros(2 * czm_mesh.nnode)
    u_max = 1.2 * pe.δ_c_t / Λ
    for s in 1:30
        u_x = u_max * s / 30
        bc_dofs, bc_vals = _bc_bottom_top(meta, czm_mesh; u_top_x=u_x, u_top_y=0.0)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
        δn, δt = seps[i]; Tn, Tt = tracts[i]
        @test abs(δn) < 0.2 * max(abs(δt), pe.δ_0_t)  # 法向保持较小
        if abs(δt) < pe.δ_0_t * 0.98
            @test Tt ≈ pe.K_t * δt rtol=5e-2 atol=1e-8
        end
    end
end

println("unit_czm_bilinear: ALL PASS")
```

**实现注意（执行时若失败按此调试）：**
- 多界面串联时，单界面 `δ` 远小于顶部位移；**不要**用 `u_y = δ_target/Λ` 期望单界面达到 `δ_target`，应用扫掠 + 自洽 `T(δ)` 比对（上文已写）。
- 若 Mode I 无法进入软化：增大 `u_max`（例如 `n_cohesive * 1.5 * δ_c / Λ`）。
- 若 Mode II 法向污染大：给左右节点加 `u_x` 仅顶/底相对、或约束侧边中点 `u_y`。

- [ ] **Step 2: 运行**

Run: `julia --project=. test/unit_czm_bilinear.jl`  
Expected: 三套 `@testset` PASS，打印 `ALL PASS`。

- [ ] **Step 3: Commit**

```bash
git add test/unit_czm_bilinear.jl
git commit -m "test(czm): unit-strip bilinear verification (Mode I/II + unload)"
```

---

### Task 4: 脚本 2 — 本征应变 + 1D 解析解

**Files:**
- Create: `test/unit_czm_eigenstrain.jl`

**Interfaces:**
- Consumes: `create_unit_czm_strip`, `assemble_thermal_chemical_load`, `unit_czm_newton_step!`, `moduli_of`
- Produces: 本地函数 `stack_1d_elastic_openings(heights, E, ε0, K_n_iface) -> (T, δs)`

**物理约束（必读）：** 底+顶固定时，**正**膨胀 → 压缩 → Mode I 不张开。必须用 **冷却**（`ΔT ≤ 0`）或使 `Σ h_i ε₀,i < 0`，产生拉伸。

- [ ] **Step 1: 实现解析解与测试脚本**

```julia
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "unit_czm_newton.jl"))

"""
1D 串联：8 杆 + 4 弹簧，固定端。
T * (Σ h/E + Σ 1/K) = -Σ h ε0
δ_j = T / K_j
要求 T>0（拉伸）时与 CZM 张开一致。
"""
function stack_1d_elastic_openings(
        heights::Vector{Float64},
        E::Vector{Float64},
        ε0::Vector{Float64},
        K_iface::Vector{Float64})
    @assert length(heights) == length(E) == length(ε0) == 8
    @assert length(K_iface) == 4
    compliance = sum(heights[i] / E[i] for i in 1:8) + sum(1 / K for K in K_iface)
    rhs = -sum(heights[i] * ε0[i] for i in 1:8)
    T = rhs / compliance
    δs = [T / K for K in K_iface]
    return T, δs
end

@testset "eigenstrain openings vs 1D analytic (elastic)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    param = case.param
    czm_mesh, meta = JuBat.create_unit_czm_strip(param; y0=1.0)
    cache = JuBat.compute_czm_params_per_interface(case)
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:4]

    α_eff = cache.by_interface[:PE_PCC].α
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0

    # 冷却斜坡：终点 ΔT* < 0，保持弹性
    ΔT_end = -0.02   # 无量纲温度增量；若 δ 过大则减小幅值
    n_steps = 10
    u = zeros(2 * czm_mesh.nnode)

    # 固定底+顶
    function bc_fixed_ends(meta)
        bc_dofs = Int[]; bc_vals = Float64[]
        for n in vcat(meta.bottom_nodes, meta.top_nodes_after_czm)
            push!(bc_dofs, 2n-1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2n);   push!(bc_vals, 0.0)
        end
        return bc_dofs, bc_vals
    end

    for s in 1:n_steps
        ΔT = ΔT_end * s / n_steps
        dT = fill(ΔT, 8)
        Δsoc_n = zeros(8)
        Δsoc_p = zeros(8)
        # 本步可不加 SOC；若需更大张开可对 NE 赋负向 Δsoc（收缩）

        F_tc = JuBat.assemble_thermal_chemical_load(
            czm_mesh, cache, α_eff, β_n, β_p, dT, Δsoc_n, Δsoc_p)
        bc_dofs, bc_vals = bc_fixed_ends(meta)
        u, seps, tracts, ok, Rn = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals,
            F_thermo_chem=F_tc)
        @test ok

        # 构造与 FEM 相同的 ε0
        ε0 = [α_eff * dT[e] + β_n * Δsoc_n[e] + β_p * Δsoc_p[e] for e in 1:8]
        E = Float64[JuBat.moduli_of(param, meta.layer_materials[e])[1] for e in 1:8]
        K_iface = Float64[
            cache.by_interface[czm_mesh.cohesive_elements[i].interface_type].K_n
            for i in 1:4]
        # 弹簧刚度需与虚功一致：解析 1D 用「单位厚度」下 T=Kδ；
        # FEM 的 K_n 为牵引-分离刚度 [应力/分离]。条带宽度 W 上合力 = ∫T dx ≈ T*W，
        # 等效轴向刚度 K_ax = K_n * W。1D 模型取单位出平面厚度、宽度 W：
        W = meta.width
        K_ax = K_iface .* W
        T_ana, δ_ana = stack_1d_elastic_openings(meta.heights, E, ε0, K_ax)
        # 分离空间：解析 δ 为位移空间张开，换到分离空间比 FEM
        for i in 1:4
            Λ = cache.by_interface[czm_mesh.cohesive_elements[i].interface_type].Λ
            δn_fem, _ = seps[i]
            δn_ana = Λ * δ_ana[i]
            @test δn_fem < cache.by_interface[:PE_PCC].δ_0_n  # 弹性段
            @test czm_mesh.damage_states[i].D ≈ 0.0 atol=1e-12
            @test δn_fem ≈ δn_ana rtol=0.05 atol=1e-10
        end
        @test T_ana > 0  # 冷却 → 拉伸
    end
end

println("unit_czm_eigenstrain: ALL PASS")
```

**调参说明：** 若 `δ` 超过 `δ_0`，减小 `|ΔT_end|`。若相对误差 >5%，检查：`K_ax = K_n*W`、平面应力泊松效应（可放宽至 `rtol=0.1`）、界面顺序与 `δ_ana` 索引一致。

- [ ] **Step 2: 运行**

Run: `julia --project=. test/unit_czm_eigenstrain.jl`  
Expected: PASS；`T_ana > 0`；`D≈0`。

- [ ] **Step 3: 若解析宽度因子有误，修正后重跑直至 PASS，再 Commit**

```bash
git add test/unit_czm_eigenstrain.jl
git commit -m "test(czm): unit-strip eigenstrain vs 1D analytic openings"
```

- [ ] **Step 4: 回写规格澄清（冷却加载）**

在 `docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md` §6.1 增加一句：

> 固定端约束下必须使 `Σ h_i ε₀,i < 0`（推荐冷却 `ΔT≤0`），以产生拉伸张开；正膨胀仅导致压缩接触，无法验证 Mode I 分离。

并更新计划摘要勾选状态。

```bash
git add docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md \
  "docs/planning-with-files/SPMe-热-内聚力耦合模型单元级验证/SPMe-热-内聚力耦合模型单元级验证.md"
git commit -m "docs(czm): note cooling eigenstrain for fixed-end strip opening"
```

---

### Task 5: 端到端验收

**Files:** 无新文件

- [ ] **Step 1: 依次运行全部相关测试**

```bash
julia --project=. test/unit_czm_strip_mesh.jl
julia --project=. test/unit_czm_bilinear.jl
julia --project=. test/unit_czm_eigenstrain.jl
```

Expected: 三个脚本均退出码 0，无 FAIL。

- [ ] **Step 2: 对照规格验收清单**

确认规格 §9 五项均满足：8+4 拓扑、双线性三工况、本征应变解析比对、可独立运行、未改生产本构。

- [ ] **Step 3: 更新计划摘要为已实现**

将 `docs/planning-with-files/SPMe-热-内聚力耦合模型单元级验证/SPMe-热-内聚力耦合模型单元级验证.md` 中状态改为：

```markdown
- [x] Brainstorming / 设计批准
- [x] 实施计划（writing-plans）
- [x] 实现与验证
```

```bash
git add "docs/planning-with-files/SPMe-热-内聚力耦合模型单元级验证/SPMe-热-内聚力耦合模型单元级验证.md"
git commit -m "docs(czm): mark unit-strip verification implementation complete"
```

---

## Spec Coverage Self-Check

| 规格条目 | 任务 |
|----------|------|
| §4 `create_unit_czm_strip` + 硬断言 | Task 1 |
| §3/§5 脚本内 Newton + 自定义 BC | Task 2–3 |
| §5 Mode I / 卸载 / Mode II | Task 3 |
| §6 合成本征应变 + 生产载荷契约 | Task 4 |
| §6 1D 解析弹性段 | Task 4 |
| §9 验收清单 | Task 5 |
| 不改生产本构 | Global + 各任务 Files 列表 |

## Placeholder / Consistency Notes

- `top_nodes_after_czm`：因 `create_czm_mesh` 可能复制顶层节点，BC 必须用复制后节点。
- 解析弹簧：`K_ax = K_n * W`（与条带宽度积分一致）；若测试显示系统偏差，优先查此因子而非乱改容差。
- `update_damage_per_interface` 未 export 时用 `JuBat.update_damage_per_interface`（模块内可见）。
