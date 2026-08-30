# 堆芯塌陷力学建模 Batch 2（K_G / C1）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现 `czm_geo_nonlinear=true` 的完全 Green-Lagrange 全 Lagrangian bulk 残差/切线（逐层各向同性 SVK + 标准初应力 K_G，无手工曲率项），接入 `assemble_bulk_residual_tangent` 的 `geo_nl` 槽位与 basic/load_substep 求解器，默认关零漂移，通过 spec §7 Batch 2 全部测试与 C1 验收门。

**Architecture:** 三段式。Task 1 在 `src/czm.jl` 新增 GL 单元装配 `gl_element_residual_tangent` 并启用 `assemble_bulk_residual_tangent` 的 `geo_nl` 槽位（本征应变内嵌：S = C:(E_GL − ε₀I)）；Task 2 建立 C1 退化/等价测试并把 Batch 1 夹具修正为生产一致的归一化网格；Task 3 把 `geo_nl`/`eigenstrain` 贯穿 `assemble_coupled_system` → `solve_czm_step`（basic + load_substep）→ 生产调用点 `update_czm_damage!`，`arc_length + geo_nl` 显式 `error` 留待 Batch 5；Task 4 跑三道门禁并记录。

**Tech Stack:** Julia 1.11.2；`Test`/`LinearAlgebra`/`SparseArrays`；既有 `IntQ4`（`src/Tools.jl:29`）。

**Spec:** `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md`（v1.3）§3.2（几何非线性）、§3.5（载荷约定衔接）、§4.1/§4.2、§7 Batch 2 行、§8；Theory/02 §1.3.7 与 v2.4 实现注记（D9）、Theory/07 式 (6.6)/(6.8)。

## Global Constraints

以下约束逐字来自 spec 与 `AGENTS.md`，每个 Task 的验收隐含包含本节。

- **运行环境**：Julia 1.11.2，单线程（`JULIA_NUM_THREADS=1`），`GKSwstype=100`，`--startup-file=no`。
- **强制行为基线**（AGENTS 9.6）：`example/testexample.jl` exit 0，全部冻结指标与 `output/testexample/testexample_results.png` SHA-256 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932` 一致（`czm_geo_nonlinear` 默认关，testexample 不受本批影响）。
- **兼容性契约**（spec §5 v1.3）：`czm_geo_nonlinear=false` 时行为、结果键、`tools/verify_czm_standalone.jl` 三方法快照与现状逐指标一致；全套测试全绿（现为 24/24，任何失败均为新回归）。
- **错误处理**（spec §6 / AGENTS 9.7）：未实现组合显式 `error`，不静默降级；`arc_length + geo_nl=true` → `error`（弧长扩展是 Batch 5，spec §4.1）。
- **D9 冻结**：物理 (x,y) 坐标完全 GL 全 Lagrangian，`K_G = ∫GᵀŜG dA` 标准初应力形式，无手工曲率项；`czm_geo_nonlinear=false` 严格走既有线性 B 矩阵代码路径（不是同一公式取零极限）。
- **提交粒度**：每 Task 一次提交，中文提交风格。

### 本计划的两项设计决策（供评审，批准后为实现契约）

**D-B2-1（本征应变内嵌 GL 残差）**：`geo_nl=true` 时 `S = C:(E_GL − ε₀I)`，ε₀ 逐单元复用 `assemble_thermal_chemical_load` 同式（`src/czm.jl:419`：`ε₀ = α_eff·dT[e] + β_n·Δsoc_n[e] + β_p·Δsoc_p[e]`，PE/NE 靠 Δsoc 零模式分工）。因此 **F_tc 在 geo_nl 路径不再单独装配**（避免双计），求解器残差变为 `R = F_ext − f_int^GL`；load_substep 的残差空间斜坡（`F_applied = f_int_ref + t·F_delta`，`src/CzmSolve.jl:520`）原样适用（F_target = F_ext）。依据：Theory/07 式 (6.6)（σ = C:(ε−εᵖ−ε*I)）一致；载荷约定已是 TL 累计制（`Δsoc = soc − cs0` 相对参考构型，`src/CouplingState.jl:506`；u 经 `czm_layout.u_prev` 累计）；spec §7 Batch 2"自由膨胀零应力"要求精确可测；spec §3.5 弧长 λ 缩放 Δε* 由此铺路。
**D-B2-2（geo_nl 的 eigenstrain 显式传递）**：`assemble_bulk_residual_tangent` 新增 `eigenstrain=nothing` 关键字（NamedTuple `(α_eff, β_n, β_p, dT, Δsn, Δsp)`，各为 `Float64`×3 + `Vector{Float64}`×3）；`nothing` 表示本次调用无本征应变（ε₀=0，运动学测试与未来纯接触场景的合法状态，非静默回退）。

### 单位契约（立项研究核定，见 findings"Batch 2 立项研究"节）

生产路径 `czm_mesh.node` 已按 x/L 归一（传 `case.param` 建网格：Rin*=11.11、Rout*=58.74）；u 为 L 归一位移（Λ=L/δ_czm=280 印证）；`moduli_of` 返回 σ_czm 归一模量。**GL 列式 F = I + ∂u/∂x 在生产构型上直接无量纲自洽，节点坐标原样使用，无需任何换算。** 测试夹具必须传 `case.param` 建网格（Batch 1 夹具误用 `param_dim`，Task 2 修正）。

## 范围与后续计划

本计划只覆盖 spec §9 的 `2(C1)`。后置：`2'`（卷绕预应力）需先确定 σ₀(r) 参数来源与形态（spec §3.7），单独小计划；`2''`（D13 网格探针）依赖本批切线/特征值能力且需薄层径向细分网格功能，另行立项。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/czm.jl` | 修改 | 新增 `gl_element_residual_tangent`；`assemble_bulk_residual_tangent` 启用 `geo_nl` 槽位 + `eigenstrain`；`assemble_coupled_system` 加 `geo_nl`/`eigenstrain` 透传 |
| `src/CzmSolve.jl` | 修改 | `solve_czm_step`/`solve_czm_basic_step`/`newton_raphson_czm`/`backtrack_line_search!` 加 `geo_nl`/`eigenstrain`；geo_nl 下跳过 F_tc、不用 K_bulk 缓存；arc_length+geo_nl 显式 error |
| `src/CouplingState.jl` | 修改 | `update_czm_damage!` 按-opt 传 `geo_nl`/`eigenstrain` |
| `test/test_czm_geometric_stiffness.jl` | 新建 | spec §7 Batch 2 前四项测试（FD/刚体转动/线性退化/自由膨胀/K_G 方向） |
| `test/test_czm_geo_c1.jl` | 新建 | C1：K_G→0 退化回 Batch 1 冻结解 + 求解器 geo_nl 冒烟 |
| `test/test_czm_mech_core.jl` | 修改 | 夹具改传 `case.param`（归一化网格，对齐生产） |

---

## Task 1: GL 单元装配与 `geo_nl` 槽位

**Files:**
- Modify: `src/czm.jl`（`gl_element_residual_tangent` 插在 `assemble_bulk_residual_tangent` 之前；后者替换 `geo_nl && error(...)` 分支）
- Test: `test/test_czm_geometric_stiffness.jl`

**Interfaces:**
- Consumes: `IntQ4(f, x_e, y_e; order)`（`src/Tools.jl:29`，回调 `(ξ,η,w,dNdx,dNdy,detJ)`）；`moduli_of(param, mt)`（`src/czm.jl:56`）；`assemble_bulk_stiffness`（`src/czm.jl:268`）。
- Produces:

```julia
gl_element_residual_tangent(x_e, y_e, u_e::Vector{Float64}, D_mat::Matrix{Float64},
                            ε0::Float64, gsorder::Int) -> (f_e::Vector{Float64}, K_e::Matrix{Float64})

assemble_bulk_residual_tangent(czm_mesh, u, param_cache, mech_state=nothing;
    geo_nl=false, plasticity=false, K_bulk_cached=nothing,
    eigenstrain=nothing)   # NamedTuple=(α_eff,β_n,β_p,dT,Δsn,Δsp)；nothing=ε₀≡0
```

- [ ] **Step 1: 写失败测试**

创建 `test/test_czm_geometric_stiffness.jl`：

```julia
using Test
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §7 Batch 2：切线 FD、大转角刚体运动零内力（D9 精确性质）、
# 线性退化（u=0 切线 ≡ 线性刚度）、自由膨胀零应力（D-B2-1）、K_G 压缩方向性。
# 夹具传 case.param（归一化网格，对齐生产；Batch 1 夹具误用 param_dim 已在别批修正）。

function build_geo_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    return case, param_cache
end

# 刚体旋转位移场：u = (R(φ) − I)·X（相对参考构型的大转角，φ=30°）
function rigid_rotation_u(node::Matrix{Float64}, φ::Float64)
    c, s = cos(φ), sin(φ)
    u = zeros(Float64, 2 * size(node, 1))
    for n in 1:size(node, 1)
        x, y = node[n, 1], node[n, 2]
        u[2*n-1] = (c - 1) * x - s * y
        u[2*n]   = s * x + (c - 1) * y
    end
    return u
end

# 均匀膨胀/压缩位移场：u = κ·X（κ 可为 √(1+2ε₀)−1 的精确自由膨胀幅值）
function uniform_scale_u(node::Matrix{Float64}, κ::Float64)
    u = zeros(Float64, 2 * size(node, 1))
    for n in 1:size(node, 1)
        u[2*n-1] = κ * node[n, 1]
        u[2*n]   = κ * node[n, 2]
    end
    return u
end

@testset "大转角刚体运动零内力（完全 GL 精确性质，D9）" begin
    case, param_cache = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u = rigid_rotation_u(czm_mesh.node, 30.0 * pi / 180)
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; geo_nl=true)
    @test norm(f_int) ≤ 1e-10 * norm(K_tan, Inf) * norm(u, Inf)
    @test !any(isnan, f_int) && !any(isnan, K_tan)
end

@testset "自由膨胀零应力（D-B2-1：ε₀ 均匀 + 精确膨胀位移）" begin
    case, param_cache = build_geo_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ε₀ = 1e-3
    κ = sqrt(1.0 + 2 * ε₀) - 1.0
    u = uniform_scale_u(czm_mesh.node, κ)
    eigenstrain = (α_eff=1.0, β_n=0.0, β_p=0.0,
                   dT=fill(ε₀, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; geo_nl=true, eigenstrain=eigenstrain)
    @test norm(f_int) ≤ 1e-10 * norm(K_tan, Inf) * norm(u, Inf)
end

@testset "零位移 GL 切线退化为线性刚度（patch/线性极限）" begin
    case, param_cache = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u0 = zeros(Float64, 2 * czm_mesh.nnode)
    _, K_geo = JuBat.assemble_bulk_residual_tangent(czm_mesh, u0, param_cache; geo_nl=true)
    K_lin = JuBat.assemble_bulk_stiffness(czm_mesh, param_cache)
    @test isapprox(K_geo, K_lin; rtol=1e-12, norm=(A, p) -> norm(Array(A), p))
end

@testset "切线有限差分（中等位移，含 K_G 与大梯度项）" begin
    case, param_cache = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u = zeros(Float64, 2 * czm_mesh.nnode)
    for n in 1:czm_mesh.nnode
        u[2*n-1] = 5e-3 * sin(3.0 * czm_mesh.node[n, 1] + 1.0)
        u[2*n]   = 5e-3 * cos(2.0 * czm_mesh.node[n, 2])
    end
    f0, K = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache; geo_nl=true)
    h = 1e-7
    for dof in [1, 2, 101, 5000, 2 * czm_mesh.nnode - 1]
        up = copy(u); up[dof] += h
        um = copy(u); um[dof] -= h
        fp, _ = JuBat.assemble_bulk_residual_tangent(czm_mesh, up, param_cache; geo_nl=true)
        fm, _ = JuBat.assemble_bulk_residual_tangent(czm_mesh, um, param_cache; geo_nl=true)
        fd = (fp .- fm) ./ (2 * h)
        @test isapprox(collect(K[:, dof]), fd; rtol=1e-6, atol=1e-8 * max(1.0, norm(fd)))
    end
end

@testset "K_G 方向性：均匀压缩降低压缩向刚度" begin
    case, param_cache = build_geo_fixture()
    czm_mesh = case.czm_mesh
    κ = -2e-3   # 压缩预载
    u = uniform_scale_u(czm_mesh.node, κ)
    d = uniform_scale_u(czm_mesh.node, 1.0)  # 压缩方向试探场
    K_lin = JuBat.assemble_bulk_stiffness(czm_mesh, param_cache)
    _, K_geo = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache; geo_nl=true)
    @test dot(d, K_geo * d) < dot(d, Array(K_lin) * d)
end
```

- [ ] **Step 2: 运行测试，确认失败**

Run: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. test/test_czm_geometric_stiffness.jl`
Expected: FAIL，`geo_nl=true 尚未实现（Batch 2…）`（Batch 1 冻结的显式 error）。

- [ ] **Step 3: 实现 GL 单元与 geo_nl 槽位**

`src/czm.jl`，在 `assemble_bulk_residual_tangent` 的文档字符串之前插入：

```julia
"""
    gl_element_residual_tangent(x_e, y_e, u_e, D_mat, ε0, gsorder) -> (f_e, K_e)

单个 Q4 的完全 Green-Lagrange 残差/切线（Batch 2，spec §3.2，D9）。坐标与位移均在
L 归一系（生产 `czm_mesh.node` 已 x/L 归一，u 为 L 归一位移），F = I + ∇u 无量纲。

- 应变（工程剪切约定，E_vec = [E11, E22, 2E12]，直接由 F = ½(FᵀF−I) 计算）；
- 本构 S = C:(E_vec − ε₀[1,1,0])（逐层各向同性 SVK，D_mat 同线性路径；D-B2-1）；
- f = ∫ B_GLᵀ S dA；K = ∫ B_GLᵀ C B_GL dA + ∫ Gᵀ Ŝ G dA（标准初应力 K_G，
  无手工曲率项——曲率由物理坐标网格几何携带）。
"""
function gl_element_residual_tangent(x_e, y_e, u_e::Vector{Float64},
                                     D_mat::Matrix{Float64}, ε0::Float64, gsorder::Int)
    f_e = zeros(Float64, 8)
    K_e = zeros(Float64, 8, 8)
    IntQ4(x_e, y_e; order=gsorder) do ξ, η, w, dNdx, dNdy, detJ
        # ∇u
        uxx = 0.0; uxy = 0.0; uyx = 0.0; uyy = 0.0
        for i in 1:4
            uxx += dNdx[i] * u_e[2*i-1]
            uxy += dNdy[i] * u_e[2*i-1]
            uyx += dNdx[i] * u_e[2*i]
            uyy += dNdy[i] * u_e[2*i]
        end
        # 完全 GL 应变（E = ½(FᵀF − I)）
        E11 = uxx + 0.5 * (uxx * uxx + uyx * uyx)
        E22 = uyy + 0.5 * (uxy * uxy + uyy * uyy)
        g12 = (uxy + uyx) + (uxx * uxy + uyx * uyy)
        S1 = D_mat[1, 1] * (E11 - ε0) + D_mat[1, 2] * (E22 - ε0)
        S2 = D_mat[1, 2] * (E11 - ε0) + D_mat[2, 2] * (E22 - ε0)
        S3 = D_mat[3, 3] * g12
        # B_GL = ∂E_vec/∂u_e（u=0 时退化为线性 B）
        B = zeros(Float64, 3, 8)
        for i in 1:4
            dx, dy = dNdx[i], dNdy[i]
            B[1, 2*i-1] = (1.0 + uxx) * dx
            B[1, 2*i]   = uyx * dx
            B[2, 2*i-1] = uxy * dy
            B[2, 2*i]   = (1.0 + uyy) * dy
            B[3, 2*i-1] = (1.0 + uxx) * dy + uxy * dx
            B[3, 2*i]   = (1.0 + uyy) * dx + uyx * dy
        end
        # 残差 f = ∫ Bᵀ S dA
        BtS = B' * [S1, S2, S3]
        axpy!(w * detJ, BtS, f_e)
        # 材料切线 ∫ Bᵀ C B dA
        DB = D_mat * B
        K_e .+= (w * detJ) .* (B' * DB)
        # 标准初应力 K_G = ∫ Gᵀ Ŝ G dA
        G = zeros(Float64, 4, 8)
        for i in 1:4
            G[1, 2*i-1] = dNdx[i]
            G[2, 2*i-1] = dNdy[i]
            G[3, 2*i]   = dNdx[i]
            G[4, 2*i]   = dNdy[i]
        end
        Sh = [S1 S3 0.0 0.0; S3 S2 0.0 0.0; 0.0 0.0 S1 S3; 0.0 0.0 S3 S2]
        K_e .+= (w * detJ) .* (G' * Sh * G)
    end
    return f_e, K_e
end
```

然后把 `assemble_bulk_residual_tangent` 中 `geo_nl && error(...)` 行替换为（`plasticity`/`mech_state` 两个 error 保持不变）：

```julia
    local K_tangent, f_int_bulk
    if geo_nl
        K_bulk_cached !== nothing && error(
            "assemble_bulk_residual_tangent: geo_nl=true 时切线依赖 u，不得传 K_bulk_cached（须逐迭代重组）。")
        α_eff = eigenstrain === nothing ? 0.0 : eigenstrain.α_eff
        β_n   = eigenstrain === nothing ? 0.0 : eigenstrain.β_n
        β_p   = eigenstrain === nothing ? 0.0 : eigenstrain.β_p
        dT_el = eigenstrain === nothing ? nothing : eigenstrain.dT
        Δsn   = eigenstrain === nothing ? nothing : eigenstrain.Δsn
        Δsp   = eigenstrain === nothing ? nothing : eigenstrain.Δsp
        if dT_el !== nothing
            (length(dT_el) == ne0 && length(Δsn) == ne0 && length(Δsp) == ne0) || throw(DimensionMismatch(
                "assemble_bulk_residual_tangent: eigenstrain 向量长度应为 bulk 单元数 $ne0"))
        end
        param = param_cache.param_ref
        submesh = czm_mesh.czm_submesh
        element = czm_mesh.bulk_element
        node = czm_mesh.node
        ne0 = size(element, 1)
        I_idx = Int64[]; J_idx = Int64[]; K_vals = Float64[]
        sizehint!(I_idx, ne0 * 64); sizehint!(J_idx, ne0 * 64); sizehint!(K_vals, ne0 * 64)
        f_gl = zeros(Float64, ndof)
        for e in 1:ne0
            E_e, ν_e = moduli_of(param, submesh.material_type[e])
            D_mat = E_e / (1.0 - ν_e^2) * [1.0 ν_e 0.0; ν_e 1.0 0.0; 0.0 0.0 (1.0 - ν_e) / 2.0]
            ε0 = 0.0
            if dT_el !== nothing
                ε0 = α_eff * dT_el[e] + β_n * Δsn[e] + β_p * Δsp[e]
            end
            elem_nodes = element[e, :]
            x_e = node[elem_nodes, 1]; y_e = node[elem_nodes, 2]
            u_e = zeros(Float64, 8)
            for (k, n) in enumerate(elem_nodes)
                u_e[2*k-1] = u[2*n-1]; u_e[2*k] = u[2*n]
            end
            f_e, K_e = gl_element_residual_tangent(x_e, y_e, u_e, D_mat, ε0, 2)
            dofs = Int64[]
            for n in elem_nodes
                push!(dofs, 2*n - 1); push!(dofs, 2*n)
            end
            for a in 1:8
                f_gl[dofs[a]] += f_e[a]
                for b in 1:8
                    push!(I_idx, dofs[a]); push!(J_idx, dofs[b]); push!(K_vals, K_e[a, b])
                end
            end
        end
        K_tangent = sparse(I_idx, J_idx, K_vals, ndof, ndof)
        f_int_bulk = f_gl
    else
        K_tangent = K_bulk_cached !== nothing ? K_bulk_cached :
                    assemble_bulk_stiffness(czm_mesh, param_cache)
        f_int_bulk = K_tangent * u
    end
```

实现说明（执行者按此落笔）：把原函数体三段 error 之后的**线弹性快路径两行**（`K_tangent = K_bulk_cached !== nothing ? …` 与 `f_int_bulk = K_tangent * u`）原样移入 `else` 分支；`ne0 = size(czm_mesh.bulk_element, 1)` 在维度检查前取得（GL 分支内亦用）；`axpy!` 来自 `LinearAlgebra`（模块级已导入，`src/JuBat.jl:2`）。GL 分支的稀疏三元组组装与既有 `assemble_bulk_stiffness`（:288-333）同模式。

- [ ] **Step 4: 运行测试，确认通过**

Run: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. test/test_czm_geometric_stiffness.jl`
Expected: 5 个 testset 全 PASS。若 FD 失败：优先检查 B_GL 第 3 行（工程剪切的因子 2）与 Ŝ 的 S3 位置；若刚体转动不为机器零：检查 E 是否含线性项丢失。修复前不得放宽容差。

- [ ] **Step 5: Batch 1 回归 + 模块加载**

Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl` → 8/8 testset PASS（geo_nl=false 路径逐位不变）。
Run: `julia --startup-file=no --project=. -e 'include("src/JuBat.jl")'` → 无警告加载。

- [ ] **Step 6: 提交**

```bash
git add src/czm.jl test/test_czm_geometric_stiffness.jl
git commit -m "feat(czm): 完全 Green-Lagrange bulk 残差/切线与 geo_nl 槽位（SVK + 标准初应力 K_G，本征应变内嵌）"
```

---

## Task 2: 夹具归一化修正与 C1 退化测试

**Files:**
- Modify: `test/test_czm_mech_core.jl`（`build_mech_core_fixture` 传 `case.param`）
- Test: `test/test_czm_geo_c1.jl`（新建）

**Interfaces:**
- Consumes: Task 1 的 `assemble_bulk_residual_tangent(geo_nl, eigenstrain)`。
- Produces: C1 验收证据（K_G→0 退化 + Batch 1 冻结解对照）；归一化夹具（后续批次复用）。

- [ ] **Step 1: 修正 Batch 1 夹具**

`test/test_czm_mech_core.jl` 的 `build_mech_core_fixture` 中，把

```julia
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, czm_enabled=true, gsorder=2)
```

替换为

```julia
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
```

并在函数末尾注释行说明：`# 归一化网格（生产路径传 case.param；Batch 1 曾误用 param_dim，本批对齐）`。运行该文件确认 8/8 testset 仍 PASS（等价断言与网格量纲无关，必须不变）。

- [ ] **Step 2: 写 C1 测试**

创建 `test/test_czm_geo_c1.jl`：

```julia
using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# C1（Theory/08 §7.5.1）：塑性关 + K_G→0 + 接触关 ⟹ 严格退回工况 R 线性理论。
# 数值判据：ε*→0 极限下 geo_nl 解趋近线性解；geo_nl=false 与 Batch 1 冻结路径逐位一致。

function build_c1_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    return case, param_cache, cache
end

@testset "C1：K_G→0 极限下 geo_nl 解趋近线性解" begin
    case, param_cache, cache = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(ndof)
    u0 = zeros(ndof)
    εtiny = 1e-6
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0,
           dT=fill(εtiny, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    # 线性参照：F_tc 显式（Batch 1 冻结路径）
    F_tc = JuBat.assemble_thermal_chemical_load(
        czm_mesh, param_cache, eig.α_eff, eig.β_n, eig.β_p, eig.dT, eig.Δsn, eig.Δsp)
    r_lin, _ = JuBat.solve_czm_step(
        czm_mesh, F_ext, param_cache, case.param, u0;
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache)
    # GL 路径：eigenstrain 内嵌
    r_geo, _ = JuBat.solve_czm_step(
        czm_mesh, F_ext, param_cache, case.param, u0;
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache,
        geo_nl=true, eigenstrain=eig)
    @test r_geo.converged && r_lin.converged
    @test isapprox(r_geo.displacement, r_lin.displacement; rtol=1e-4, atol=1e-12)
end

@testset "geo_nl=false 与 Batch 1 冻结解逐位一致（回归锚）" begin
    case, param_cache, cache = build_c1_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    eig = (α_eff=1.0, β_n=0.0, β_p=0.0, dT=fill(1e-4, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    r1, _ = JuBat.solve_czm_step(
        czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache)
    r2, _ = JuBat.solve_czm_step(
        czm_mesh, zeros(ndof), param_cache, case.param, zeros(ndof);
        α_eff=eig.α_eff, β_n=eig.β_n, β_p=eig.β_p,
        dT_elem=eig.dT, Δsoc_n_elem=eig.Δsn, Δsoc_p_elem=eig.Δsp,
        max_iter=100, tol=1e-10, iter_method="basic", cache=cache,
        geo_nl=false)
    @test r1.displacement == r2.displacement
end
```

- [ ] **Step 3: 运行，确认按预期失败**

Run: `julia --startup-file=no --project=. test/test_czm_geo_c1.jl`
Expected: 第一个 testset FAIL（`solve_czm_step` 无 `geo_nl` 关键字 → MethodError/unsupported keyword），第二个 testset PASS。

- [ ] **Step 4: 提交（仅夹具修正；C1 测试随 Task 3 变绿后一并提交，避免批次中途套件带红）**

```bash
git add test/test_czm_mech_core.jl
git commit -m "test(czm): mech_core 夹具对齐生产归一化网格（case.param 建网格）"
```

---

## Task 3: 求解器与生产路径接线

**Files:**
- Modify: `src/czm.jl`（`assemble_coupled_system` 签名与 bulk 分支）
- Modify: `src/CzmSolve.jl:651`（`solve_czm_step`）、`:107`（`backtrack_line_search!`）、basic 求解器（`:150-235` 区域）、`newton_raphson_czm`（`:481-640`）
- Modify: `src/CouplingState.jl:596`（生产调用点）

**Interfaces:**
- Consumes: Task 1 的 bulk 入口 `(geo_nl, eigenstrain)`。
- Produces:

```julia
assemble_coupled_system(...; geo_nl::Bool=false, eigenstrain=nothing)   # 新增两 kw，默认路径逐位不变
solve_czm_step(...; geo_nl::Bool=false, eigenstrain=nothing)            # 分发层透传
solve_czm_basic_step / newton_raphson_czm(...; geo_nl::Bool=false, eigenstrain=nothing)
# arc_length 方法收到 geo_nl=true → error("...Batch 5 弧长扩展...")
```

- [ ] **Step 1: `assemble_coupled_system` 透传**

签名追加 `geo_nl::Bool=false, eigenstrain=nothing`；bulk 分支改为：

```julia
    f_int_bulk, K_bulk = assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; K_bulk_cached=K_bulk_cached,
        geo_nl=geo_nl, eigenstrain=eigenstrain)
```

（`K_bulk_cached` 与 `geo_nl=true` 同时给出的冲突由 Task 1 的显式 error 拦截。）cohesive 分支与其余代码不动。

- [ ] **Step 2: 求解器分发与三个方法**

`solve_czm_step` 签名追加 `geo_nl::Bool=false, eigenstrain=nothing`：
- `load_substep`/`basic` 分支透传给 `newton_raphson_czm`/`solve_czm_basic_step`；
- `arc_length` 分支开头加：`geo_nl && error("solve_czm_step: arc_length + geo_nl 在 Batch 5 弧长扩展（spec §4.1）后支持，当前显式拒绝。")`。

`newton_raphson_czm`（:481，即 load_substep 方法）与 `solve_czm_basic_step`（:168；注意 `solve_czm_arc_length_step` **嵌套**在其体内 :265，本批不接线）：
1. 两个函数签名追加 `geo_nl::Bool=false, eigenstrain=nothing`；`backtrack_line_search!`（:107，装配点 :113）同样追加并透传；
2. 载荷准备改为条件式：

```julia
    if geo_nl
        F_thermo_chem_total = zeros(Float64, ndof)   # ε* 已内嵌 f_int^GL（D-B2-1），不再外载
        kb_cache = nothing                            # 切线依赖 u，禁用缓存
        eig_kwargs = (geo_nl=true, eigenstrain=eigenstrain)
    else
        F_thermo_chem_total = assemble_thermal_chemical_load(...)
        kb_cache = cache !== nothing ? cache.K_bulk : nothing
        eig_kwargs = ()
    end
```

3. 逐处替换 `assemble_coupled_system` 调用（用 `kb_cache` 替换 `K_bulk_cached` 变量并拼入 `eig_kwargs...`）：`newton_raphson_czm` 的 `:502`（f_int_ref）、`:530`（迭代内）、`:567`（线搜索 trial）、`:612`（收尾）；`solve_czm_basic_step` 的 `:197`（迭代内）、`:242`（收尾）与 `:226` 的 `backtrack_line_search!` 调用（其内部 `:113`）。**不得触碰**嵌套 `solve_czm_arc_length_step`（:307/:317/:348/:369/:439）——它只能经 `solve_czm_step` 分发进入，分发层已显式拒绝。

- [ ] **Step 3: 生产调用点**

`src/CouplingState.jl:596` 的 `solve_czm_step(...)` 调用追加：

```julia
        geo_nl = case.opt.czm_geo_nonlinear,
        eigenstrain = case.opt.czm_geo_nonlinear ?
            (α_eff=α_eff, β_n=β_n, β_p=β_p,
             dT=dT_elem, Δsn=Δsoc_n_elem, Δsp=Δsoc_p_elem) : nothing
```

- [ ] **Step 4: 定向与套件测试**

Run: `julia --startup-file=no --project=. test/test_czm_geo_c1.jl` → 2/2 testset PASS。
Run: `julia --startup-file=no --project=. test/test_czm_geometric_stiffness.jl` → 5/5 PASS。
Run: `julia --startup-file=no --project=. test/test_czm_mech_core.jl` → 8/8 PASS。
Run: `julia --startup-file=no --project=. test/test_assemble_coupled_system.jl` → PASS。

- [ ] **Step 5: 提交（含 Task 2 写好的 C1 测试，此时已变绿）**

```bash
git add src/czm.jl src/CzmSolve.jl src/CouplingState.jl test/test_czm_geo_c1.jl
git commit -m "feat(czm): geo_nl/eigenstrain 贯穿装配与 basic/load_substep 求解器，生产路径按 opt 接入（arc_length 显式拒绝至 Batch 5）；C1 测试转绿"
```

---

## Task 4: 三道门禁与批次记录

- [ ] **Step 1: 全套测试**

Run: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. test/runtests.jl`
Expected: 26/26 通过（24 既有 + 2 新文件），任何失败即停止。

- [ ] **Step 2: 基线快照门禁**

Run: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. tools/verify_czm_standalone.jl`
Expected: 与 `baseline_czm_standalone.md` 冻结表逐位一致（geo_nl 默认关，三方法行为不变）。

- [ ] **Step 3: 强制行为基线**

Run: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no example/testexample.jl`
Expected: exit 0；全部冻结指标与 PNG SHA-256 `4ba6207c…e932` 一致（路径 `output/testexample/testexample_results.png`）。

- [ ] **Step 4: 记录与提交**

`Simplify/baseline.md` 批次表追加 Batch 2 行（行数变化、26/26、快照一致、指标一致、PASS）；`docs/planning-with-files/30_堆芯塌陷力学建模/progress.md` 追加 Batch 2 小节（含 D-B2-1/D-B2-2 决策执行情况、C1 证据、门禁实测）；`index.md` 更新任务行。

```bash
git add Simplify/baseline.md "docs/planning-with-files/30_堆芯塌陷力学建模/progress.md" docs/planning-with-files/index.md
git commit -m "docs(baseline): Batch 2 门禁记录（C1 通过，快照与基线一致）"
```

**Batch 2 完成门**：spec §7 Batch 2 五项测试全绿 + C1 + 三道门禁全绿。未全绿不得宣称 C1，也不得进入 Batch 2'。

---

## 自评审记录

1. **spec 覆盖**：§3.2 完全 GL/TL/SVK/标准 K_G → Task 1；§3.2 末"false 走严格线性路径" → Task 1 Step 5 + Task 2 回归锚；§7 Batch 2 五项测试 → Task 1（FD/刚体/patch-线性退化/膨胀）+ Task 2（K_G→0 + 冻结解）；§4.1 CzmSolve"每次迭代重组切线 + 弹性快路径" → Task 3（geo_nl 禁缓存、默认路径缓存不变）；§4.1 "弧长扩展为全机械残动"属 Batch 5 → D-B2-2 显式 error。§3.5 弧长 λ 缩放 Δε* 由 D-B2-1 铺路（本批不实现）。
2. **决策登记**：D-B2-1（ε* 内嵌残差、geo_nl 下 F_tc 不装配）与 D-B2-2（eigenstrain 关键字、nothing≡零）为本计划新增设计决策，依据 findings"Batch 2 立项研究"；`assemble_bulk_residual_tangent` 因 D-B2-2 再加一个关键字（对 Batch 1 签名偏差清单的第四项，additive 默认值，不影响既有测试）。
3. **占位符扫描**：Task 1 Step 3 的 GL 分支给出了完整数据流但允许执行者按"包 else 分支"的明确说明落笔（非占位——改动边界、变量名、稀疏模式全部指定）；其余步骤均有完整代码或逐处机械替换清单。
4. **类型一致性**：`eigenstrain` NamedTuple 字段名（α_eff/β_n/β_p/dT/Δsn/Δsp）在 Task 1/2/3 三处一致；`gl_element_residual_tangent` 签名在 Interfaces/实现/测试一致；`solve_czm_step` 新关键字与 CouplingState 调用一致。
5. **已知风险**：GL 逐迭代重组使 geo_nl 开启时 CZM 步耗时上升（testexample 默认关不受影响）；spec §8 性能行的冒烟预算在 Batch 7。
