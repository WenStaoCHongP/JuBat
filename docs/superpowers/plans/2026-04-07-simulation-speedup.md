# Simulation Speedup Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **审查修订**: 2026-04-08 — 基于源码逐文件对照审查，修正了与既有架构冲突的内容，重新排列优先级。

**Goal:** Reduce JuBat multi-SPMe simulation wall-clock time by 15-25% through targeted optimization of the SPMe hot path, memory allocation reduction, and thermal assembly caching.

**Architecture:** Four-batch approach — (1) deepcopy 清理（即时见效）, (2) 静态计算缓存（类型安全）, (3) 原位操作与 @views, (4) SPMe 工作区（高收益高风险）. Each task is independently testable by running `example/testexample.jl` and comparing timing output.

**Tech Stack:** Julia 1.x, SparseArrays, Threads.@threads, FEM assembly

**Baseline (testexample.jl, czm_enabled=false, nθ=80):**
| Module | Time Ratio | Target After |
|--------|-----------|-------------|
| SPMe solve | 91.59% | < 88% |
| Thermal distributed | 7.85% | < 6% |
| Branch solver | 0.56% | 0.56% (unchanged) |

---

## 设计约束

> **重要**: 以下约束来自之前完成的结构体重构（`00_优化分析进度.md`），新增优化不得违反。

1. **禁止恢复 `thermal_extras::Dict{String,Any}`** — 该字段已从 `Case` struct 中刻意移除，所有缓存需使用类型安全的 struct 字段。
2. **Fail-Fast 原则** — 禁止 `haskey`/`get(dict, default)`/`try-catch` 吞错误，直接字段访问。
3. **main 分支不动** — SPM/SPMe/P2D 路径代码保持不变，仅影响 multi-SPMe 路径。

---

## File Structure

| File | Change Type | Responsibility |
|------|------------|---------------|
| `src/Solve.jl` | Modify | 替换全部 deepcopy → copy/直接赋值 |
| `src/Assemble.jl` | Modify | deepcopy(KI) → zeros；新增 Assemble! 预分配变体 |
| `src/CouplingState.jl` | Modify | MultiSPMeLayout 新增 areas 字段；新增 BoundaryEdgeCache struct + compute_boundary_edge_cache 函数 |
| `src/Initialisation.jl` | Modify | 缓存单元面积到 layout |
| `src/CallModel.jl` | Modify | 使用缓存面积；@views 状态提取；线程本地工作区 + copy_element_results |
| `src/ThermalDistributed.jl` | Modify | 原位热 BC；使用 BoundaryEdgeCache；使用 Assemble! |
| `src/SPMe.jl` | Modify | 接受 AbstractVector（兼容 view）；可选工作区 |
| `src/Jellyrollmodel.jl` | Modify | setup_thermal2D_mesh 中缓存边界边 |
| `src/Variables.jl` | Modify | 可选: create_element_workspace（仅 Chunk 4 需要） |

---

## Chunk 1: deepcopy 清理（低风险，即时见效）

### Task 1: 替换 Solve.jl 中全部 deepcopy

`deepcopy` 遍历整个对象图（对稀疏矩阵极慢）。`copy` 对 `Array{Float64}` 和 `SparseMatrixCSC` 产生浅拷贝，足够安全。

**代码审查发现**: 除计划原覆盖的 L245-247 外，Solve.jl 中还有 4 处 Float64 的 `deepcopy` 完全多余（Float64 是不可变值类型，直接赋值等价）。

**Files:**
- Modify: `src/Solve.jl`

- [ ] **Step 1: 替换 L245-247 三处 deepcopy（矩阵/向量）**

```julia
# OLD (line 245-247):
y_old = deepcopy(y_new)
K_old = deepcopy(K_new)
F_old = deepcopy(F_new)

# NEW:
y_old = copy(y_new)
K_old = copy(K_new)
F_old = copy(F_new)
```

**安全性**: `K_new`/`F_new` 来自 `CallModel` 每步新建，不会被原位修改。`y_new = vcat(y_c, y_phi)` 也是每步新建。

- [ ] **Step 2: 替换 4 处 Float64 deepcopy（L123, L230, L235, L239）**

```julia
# OLD (line 123):
dt = deepcopy(dt_min)
# NEW:
dt = dt_min

# OLD (line 230):
dt = deepcopy(dt_min)
# NEW:
dt = dt_min

# OLD (line 235):
dt = deepcopy(dt_temp)
# NEW:
dt = dt_temp

# OLD (line 239):
dt_temp = deepcopy(dt)
# NEW:
dt_temp = dt
```

- [ ] **Step 3: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果。

- [ ] **Step 4: Commit**

```bash
git add src/Solve.jl
git commit -m "perf: replace all deepcopy with copy/direct assign in Solve.jl"
```

---

### Task 2: 修复 Assemble.jl deepcopy(KI) → zeros

`deepcopy(KI)` 其中 `KI = zeros(Int64, ...)` — deepcopy 遍历整个对象图完全多余。

**Files:**
- Modify: `src/Assemble.jl:13`

- [ ] **Step 1: 替换 deepcopy**

```julia
# OLD (line 13):
KJ = deepcopy(KI)

# NEW:
KJ = zeros(Int64, length(KI))
```

- [ ] **Step 2: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果。

- [ ] **Step 3: Commit**

```bash
git add src/Assemble.jl
git commit -m "perf: replace deepcopy(KI) with zeros in Assemble.jl"
```

---

## Chunk 2: 静态计算缓存（需改 struct，类型安全）

### Task 3: 缓存单元面积到 MultiSPMeLayout

单元面积在 `CallModel_MultiSPMe` L27-34 每步重新计算，但网格不变。缓存到 layout 初始化一次。

**Files:**
- Modify: `src/CouplingState.jl` (MultiSPMeLayout 新增 `areas` 字段)
- Modify: `src/Initialisation.jl:110` (传入 mesh 计算面积)
- Modify: `src/Solve.jl:104` (传入 mesh 计算面积)
- Modify: `src/CallModel.jl:27-34` (使用缓存面积)

- [ ] **Step 1: MultiSPMeLayout 新增 areas 字段，更新两个构造器**

在 `src/CouplingState.jl` 更新 struct 和两个构造器：

```julia
struct MultiSPMeLayout
    ne::Int
    n_chem::Int
    nT::Int
    n_total::Int
    chem_range::UnitRange{Int}
    thermal_range::UnitRange{Int}
    areas::Vector{Float64}          # NEW: 预计算的单元面积（网格不变量）
end

# 保留原 3 参数构造器（兼容 Initialisation.jl / Solve.jl 中无 mesh 的调用点）
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)
    MultiSPMeLayout(ne, n_chem, nT,
                    ne * n_chem + nT,
                    1:(ne * n_chem),
                    (ne * n_chem + 1):(ne * n_chem + nT),
                    zeros(Float64, ne))   # areas 延迟到首次使用时填充
end

# 新增 4 参数构造器：接收 mesh 计算面积
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int, mesh_th)
    thermal_range = (ne * n_chem + 1):(ne * n_chem + nT)
    areas = zeros(Float64, ne)
    ngs = length(mesh_th.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh_th.gs.ele[g]
        areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    return MultiSPMeLayout(ne, n_chem, nT, ne * n_chem + nT,
                           1:(ne * n_chem), thermal_range, areas)
end
```

> **注意**: 3 参数构造器的 `areas=zeros(ne)` 为占位值。实际使用缓存的调用点（Initialisation.jl、Solve.jl）需在 Step 2/3 改为 4 参数版本。

- [ ] **Step 2: 更新 Initialisation.jl 中的布局构造**

```julia
# OLD (line 110):
case.layout = MultiSPMeLayout(ne, n_chem, nT)

# NEW:
case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
```

- [ ] **Step 3: 更新 Solve.jl 中的布局构造**

```julia
# OLD (line 104):
case.layout = MultiSPMeLayout(ne, n_chem, nT)

# NEW:
case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
```

- [ ] **Step 4: CallModel.jl 使用缓存面积**

```julia
# OLD (L27-34, 7 行):
A = zeros(Float64, ne)
ngs = length(mesh_th.gs.detJ)
@inbounds for g in 1:ngs
    e = mesh_th.gs.ele[g]
    A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
end
variables["thermal2D element area"] = A
areas = A

# NEW (2 行):
areas = layout.areas
variables["thermal2D element area"] = areas
```

- [ ] **Step 5: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果。

- [ ] **Step 6: Commit**

```bash
git add src/CouplingState.jl src/Initialisation.jl src/Solve.jl src/CallModel.jl
git commit -m "perf: cache element areas in MultiSPMeLayout (eliminates per-step recomputation)"
```

---

### Task 4: 缓存边界边列表到 MeshGeometry

`apply_convection_bc`（`ThermalDistributed.jl:66-76`）每步重建 `seen` Set 去重边界边。预计算一次存入 `MeshGeometry`。

**设计**: 在 `CouplingState.jl` 新增 `BoundaryEdgeCache` struct，存入 `MeshGeometry` 新字段。在 `setup_thermal2D_mesh` 阶段一次性计算。不使用 Dict。

**Files:**
- Modify: `src/CouplingState.jl` (新增 BoundaryEdgeCache struct；MeshGeometry 新增 boundary_edges 字段)
- Modify: `src/Jellyrollmodel.jl` (setup_thermal2D_mesh 中计算缓存)
- Modify: `src/ThermalDistributed.jl` (apply_convection_bc 使用缓存)

- [ ] **Step 1: 新增 BoundaryEdgeCache struct + compute_boundary_edge_cache 函数**

在 `src/CouplingState.jl` 添加：

```julia
"""
    BoundaryEdgeCache

预计算的外边界边列表（网格不变量），用于对流边界条件装配。
"""
struct BoundaryEdgeCache
    edges::Vector{Tuple{Int,Int}}   # (node_a, node_b) 对，a < b
    L_edge::Vector{Float64}         # 边长（无量纲）
end

"""
    compute_boundary_edge_cache(mesh, is_outer)

从网格和外部节点标记中提取去重的外边界边列表。
返回 BoundaryEdgeCache（边对 + 边长），仅计算一次。
"""
function compute_boundary_edge_cache(mesh, is_outer)
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    edges = Tuple{Int,Int}[]
    L_edge = Float64[]

    for e in 1:ne
        nodes = mesh.element[e, :]
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]),
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            (is_outer[a] && is_outer[b]) || continue
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            push!(edges, key)
            push!(L_edge, hypot(x[b] - x[a], y[b] - y[a]))
        end
    end
    return BoundaryEdgeCache(edges, L_edge)
end
```

- [ ] **Step 2: MeshGeometry 新增字段**

```julia
struct MeshGeometry
    element_layer::Vector{Int}
    is_inner_layer::Vector{Bool}
    layer_weights::Matrix{Float64}
    interface_pairs::Vector{Tuple{Int,Int}}
    czm_element_map::Dict{Int,Vector{Int}}
    inner_nodes::Vector{Int}
    outer_nodes::Vector{Int}
    boundary_edges::Union{Nothing, BoundaryEdgeCache}  # NEW: 预计算的边界边
end
```

- [ ] **Step 3: 在 setup_thermal2D_mesh 中计算缓存**

在 `src/Jellyrollmodel.jl` 的 `setup_thermal2D_mesh` 中，`MeshGeometry` 构造前计算边界边：

```julia
# 在构造 MeshGeometry 之前：
is_inner, is_outer = identify_boundary_nodes(mesh_th, case.param)
boundary_cache = compute_boundary_edge_cache(mesh_th, is_outer)

case.geometry = MeshGeometry(
    element_layer, is_inner_layer, layer_weights,
    interface_pairs, czm_element_map,
    inner_nodes, outer_nodes,
    boundary_cache   # NEW
)
```

- [ ] **Step 4: 修改 apply_convection_bc 使用缓存**

```julia
function apply_convection_bc(KT, FT, mesh, is_outer, case; edge_cache=nothing)
    K = copy(KT)
    F = copy(FT)
    Bi = case.param_dim.scale.h * case.param.cell.lambda_r
    if Bi == 0
        return K, F
    end

    param = case.param
    T_amb = param.cell.T_amb
    s_vals = (-0.577350269189626, 0.577350269189626)
    w_vals = (1.0, 1.0)

    # 使用缓存（来自 case.geometry.boundary_edges）
    if edge_cache === nothing
        # 兼容回退：从 is_outer 实时计算（ThermalRing2D_BC 等调用方传入 is_outer）
        if is_outer === nothing
            is_inner, is_outer = identify_boundary_nodes(mesh, case.param)
        end
        edge_cache = compute_boundary_edge_cache(mesh, is_outer)
    end

    for (idx, (a, b)) in enumerate(edge_cache.edges)
        J = edge_cache.L_edge[idx] / 2
        ke11, ke12, ke22 = 0.0, 0.0, 0.0
        fe1, fe2 = 0.0, 0.0

        for (s, w) in zip(s_vals, w_vals)
            N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
            wt = Bi * w * J
            ke11 += -wt * N1 * N1
            ke12 += -wt * N1 * N2
            ke22 += -wt * N2 * N2
            fe1 += wt * T_amb * N1
            fe2 += wt * T_amb * N2
        end
        K[a, a] += ke11; K[a, b] += ke12
        K[b, a] += ke12; K[b, b] += ke22
        F[a] += fe1; F[b] += fe2
    end
    return K, F
end
```

- [ ] **Step 5: ThermalDistributed2D_BC 传入缓存**

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    K = copy(KT)
    F = copy(FT)

    if case.opt.czm_enabled
        # ... (CZM 块不变)
    end

    # 从 geometry 获取预计算的边界边缓存
    edge_cache = case.geometry !== nothing ? case.geometry.boundary_edges : nothing
    K, F = apply_convection_bc(K, F, mesh, nothing, case; edge_cache=edge_cache)
    K, F = apply_cool_method(K, F, mesh, case)
    return K, F
end
```

- [ ] **Step 6: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果；热模型计时微降。

- [ ] **Step 7: Commit**

```bash
git add src/CouplingState.jl src/Jellyrollmodel.jl src/ThermalDistributed.jl
git commit -m "perf: cache boundary edge list in MeshGeometry (eliminates per-step Set rebuild)"
```

---

## Chunk 3: 原位操作与 @views（需验证调用链）

### Task 5: 使用 @views 替代 extract_element_state 拷贝

`extract_element_state`（`Initialisation.jl:121-124`）返回 `y[(offset+1):(offset+n_chem)]` 切片拷贝，`CallModel.jl:23` 再 `vec()` 又拷贝一次。用 `@view` 直接切片消除 ne 次拷贝。

**安全性**: `yt` 即 `y_old = copy(y_new)`（Task 1 实施后），view 指向的底层数组在并行循环期间不会被修改。

**Files:**
- Modify: `src/CallModel.jl:21-24`
- Modify: `src/SPMe.jl` (SPMe_element 接受 AbstractVector，移除 vec())

- [ ] **Step 1: CallModel.jl 使用 @views**

```julia
# OLD (L21-24):
yt_chem = Vector{Vector{Float64}}(undef, ne)
for e in 1:ne
    yt_chem[e] = vec(extract_element_state(yt, e, case.layout))
end

# NEW (使用 views，消除拷贝):
yt_chem = Vector{SubArray{Float64,1}}(undef, ne)
for e in 1:ne
    offset = (e - 1) * layout.n_chem
    yt_chem[e] = @view yt[(offset + 1):(offset + layout.n_chem)]
end
```

- [ ] **Step 2: SPMe_element 接受 AbstractVector，移除 vec()**

在 `src/SPMe.jl` 中，将 `SPMe_element` 的参数类型改为 `AbstractVector{Float64}`（兼容 view），并移除内部的 `vec(yt_e)` 调用：

```julia
# 旧签名:
function SPMe_element(case::Case, yt_e::Array{Float64}, ...)

# 新签名:
function SPMe_element(case::Case, yt_e::AbstractVector{Float64}, ...)

# 移除内部的:
yt_e_vec = vec(yt_e)   # 删除此行，所有 yt_e_vec → yt_e
```

- [ ] **Step 3: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果。

- [ ] **Step 4: Commit**

```bash
git add src/CallModel.jl src/SPMe.jl
git commit -m "perf: use @views for element state extraction, eliminate vec() copy"
```

---

### Task 6: 消除热边界条件三重 copy

**代码审查发现**: 当前热 BC 路径存在**三重 copy**:
1. `ThermalDistributed2D_BC` L189-190: `K = copy(KT); F = copy(FT)`
2. `apply_convection_bc` L50-51: `K = copy(KT); F = copy(FT)` — 第二次
3. `apply_cool_method` L111/124/146/172: 每分支 `copy(KT), copy(FT)` — 第三次

调用方 `CallModel.jl:136`: `KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)` — 返回值覆盖原变量，KT/FT 原始引用不被复用。**原位修改安全**。

**Files:**
- Modify: `src/ThermalDistributed.jl:49-219`

- [ ] **Step 1: 新增原位变体 apply_convection_bc!**

```julia
function apply_convection_bc!(K, F, case; edge_cache=nothing)
    # 原位版本：直接修改 K, F，不做 copy
    Bi = case.param_dim.scale.h * case.param.cell.lambda_r
    if Bi == 0
        return K, F
    end

    if edge_cache === nothing
        # 无缓存时回退到 copy 版本（不应发生在 Task 4 实施后）
        return apply_convection_bc(K, F, case.mesh["thermal2D"], nothing, case)
    end

    param = case.param
    T_amb = param.cell.T_amb
    s_vals = (-0.577350269189626, 0.577350269189626)
    w_vals = (1.0, 1.0)

    for (idx, (a, b)) in enumerate(edge_cache.edges)
        J = edge_cache.L_edge[idx] / 2
        ke11, ke12, ke22 = 0.0, 0.0, 0.0
        fe1, fe2 = 0.0, 0.0

        for (s, w) in zip(s_vals, w_vals)
            N1, N2 = 0.5 * (1 - s), 0.5 * (1 + s)
            wt = Bi * w * J
            ke11 += -wt * N1 * N1
            ke12 += -wt * N1 * N2
            ke22 += -wt * N2 * N2
            fe1 += wt * T_amb * N1
            fe2 += wt * T_amb * N2
        end
        K[a, a] += ke11; K[a, b] += ke12
        K[b, a] += ke12; K[b, b] += ke22
        F[a] += fe1; F[b] += fe2
    end
    return K, F
end
```

- [ ] **Step 2: 新增原位变体 apply_cool_method!**

```julia
function apply_cool_method!(K, F, mesh, case)
    # 原位版本：直接修改 K, F
    cool_method = case.opt.cool_method
    if cool_method == "none"
        return K, F
    elseif cool_method == "surface"
        # ... (surface 逻辑不变，但不做 copy)
    elseif cool_method == "tab"
        # ... (tab 逻辑不变，但不做 copy)
    end
    return K, F
end
```

- [ ] **Step 3: ThermalDistributed2D_BC 使用原位变体**

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    # 直接使用 KT, FT，不 copy（调用方已确认不复用原始值）
    K = KT
    F = FT

    if case.opt.czm_enabled
        # CZM 块不变（已在原位操作 K）
    end

    edge_cache = case.geometry !== nothing ? case.geometry.boundary_edges : nothing
    K, F = apply_convection_bc!(K, F, case; edge_cache=edge_cache)
    K, F = apply_cool_method!(K, F, mesh, case)
    return K, F
end
```

- [ ] **Step 4: 保留原 copy 版本作为兼容接口**

保留 `apply_convection_bc`（带 copy 的版本）供 `ThermalRing2D_BC` 等调用方使用，避免破坏 main 分支接口。

- [ ] **Step 5: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 完全一致的结果；热模型计时下降。

- [ ] **Step 6: Commit**

```bash
git add src/ThermalDistributed.jl
git commit -m "perf: in-place thermal BC (eliminates triple copy of KT/FT)"
```

---

## Chunk 4: SPMe 精简型线程工作区（高收益，中风险）

> **收益估计**: 消除 ~2400 次无用数组分配/步（StandardVariables 中 ~30 个 thermal2D 键），+ 对象池复用消除 ~1600 次有用数组分配/步。总分配量从 ~4050 次/步降至 ~60 次（预分配一次）。
> **风险**: 需维护 `create_element_workspace` + `SPMe_variables!` 两套代码，但与 `StandardVariables`/`SPMe_variables` 共享计算逻辑，非独立复制。

### Task 7: 精简型线程本地工作区预分配

**问题**: 每次 `SPMe_element` 调用 → `SPMe_variables` → `StandardVariables(case, 1)` 创建新 Dict + ~50 个预分配数组。其中 **~30 个是 `distributed2D` 专用的热变量数组**（`Variables.jl` L94-135），`SPMe_element` 内部从不使用。ne=80 个单元/步 × ~50 数组 = ~4000 次堆分配，其中 **~60% 完全无用**。

**方案**: 将对象池复用（方向 1 方案 A）与精简无用分配（方向 1 方案 B）结合：
1. 新增 `create_element_workspace(case)` — 仅创建 SPMe_element 需要的 ~30 个键（排除 ~30 个 thermal2D 键）
2. 每线程预分配一份 workspace，并行循环中复用
3. 新增 `SPMe_variables!` 原位变体，直接写入预分配 workspace

**Files:**
- Modify: `src/Variables.jl` (新增 `create_element_workspace`)
- Modify: `src/SPMe.jl` (新增 `SPMe_variables!` 原位变体，修改 `SPMe_element` 接受可选 workspace)
- Modify: `src/CallModel.jl` (使用精简型线程本地工作区)

**SPMe_element 调用链键使用分析**（基于 `SPMe.jl`, `ElectrolyteDiffusion.jl`, `Mechanical.jl`, `CallModel.jl` 源码审查）:

| 调用者 | 读取键 | 写入键 |
|--------|--------|--------|
| `SPMe_variables` 状态提取 (L139-141) | — | 通过 `case.index` 从 yt 覆写 ~7 个键 |
| `SPMe_variables` 计算 (L157-196) | cn/cp_surf, ce_n/p/sp, T | ~17 个计算结果键 |
| `ElectrolyteDiffusion` (L14-16, L23) | ce_n/p/sp_gs, temperature | — |
| `SPMe_BC` (L97-98) | j_n, j_p (interfacial current) | — |
| `Mechanicaloutput` (L6-11, L17-29) | c_n/p, eta_n/p, V_cell, T | ~11 个力学结果键 |
| `SPMe_element` (L48-49, L56-57) | gauss point conc, coupling coeff | element index |
| `CallModel` 辅助循环 (L118-127) | eta_n/p, cn/cp_surf, cn/csp_data | — |

**排除的 thermal2D 键**（StandardVariables L94-135，共 ~30 个，SPMe_element 内部不使用）:
> thermal2D temperature/history, heat_source_fields, q_rxn_ne/pe, q_rev_ne/pe, q_ohm_s/e_ne/pe, q_sp, q_pcc, q_ncc, element current, eta_n/p_e, dUdT_n/p_e, soc_n/p, voltages, OCV, active_mask, stress/diffusion/total/strain, n_cutoff, nearest_cutoff, margin, total heat source, displacement x/y

这些键由 `CallModel_MultiSPMe` 的辅助循环（L116-129）和 `compute_heat_sources` **从 `variables_elems[e]` 提取后写入全局 variables Dict**，而非在 SPMe_element 内部创建。

- [ ] **Step 1: 新增 `create_element_workspace(case)`**

在 `src/Variables.jl` 末尾新增精简型 workspace 工厂函数：

```julia
"""
    create_element_workspace(case)

创建精简型单元工作区 Dict，仅包含 SPMe_element 调用链实际需要的 ~30 个键。
排除 distributed2D 专用的 ~30 个 thermal2D 键（由 CallModel_MultiSPMe
在单元循环外独立管理）。比 StandardVariables(case, 1) 减少约 60% 数组分配。
"""
function create_element_workspace(case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nn = 1; Np = 1  # SPMe 模式
    Ne_ngs = case.opt.Nn * case.opt.gsorder
    Ne_pgs = case.opt.Np * case.opt.gsorder

    ws = Dict{String, Union{Array{Float64}, Float64}}()

    # ── 状态提取键（case.index 对应，SPMe_variables L139-141 从 yt 写入）──
    ws["negative particle lithium concentration"] = zeros(Float64, Nrn, 1)
    ws["positive particle lithium concentration"] = zeros(Float64, Nrp, 1)
    ws["negative particle surface lithium concentration"] = zeros(Float64, Nn, 1)
    ws["positive particle surface lithium concentration"] = zeros(Float64, Np, 1)

    # ── SPMe 专用键 ──
    if case.opt.model == "SPMe"
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        ws["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n, 1)
        ws["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p, 1)
        ws["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp, 1)
        # Gauss 点计算结果（ElectrolyteDiffusion 读取）
        ws["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs, 1)
    end

    ws["temperature"] = 0.0

    # ── SPMe_variables 计算结果键 ──
    ws["cell voltage"] = 0.0
    ws["time"] = 0.0
    ws["cell current"] = 0.0
    ws["negative electrode exchange current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode exchange current density"] = zeros(Float64, Np, 1)
    ws["negative electrode interfacial current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode interfacial current density"] = zeros(Float64, Np, 1)
    ws["negative electrode overpotential"] = zeros(Float64, Nn, 1)
    ws["positive electrode overpotential"] = zeros(Float64, Np, 1)
    ws["negative electrode open circuit potential"] = zeros(Float64, Nn, 1)
    ws["positive electrode open circuit potential"] = zeros(Float64, Np, 1)

    # ── Mechanicaloutput 结果键（条件）──
    if case.opt.mechanicalmodel == "full"
        ws["negative particle center radial stress"] = zeros(Float64, Nn, 1)
        ws["positive particle center radial stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface tangential stress"] = zeros(Float64, Nn, 1)
        ws["positive particle surface tangential stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface displacement"] = zeros(Float64, Nn, 1)
        ws["positive particle surface displacement"] = zeros(Float64, Np, 1)
        ws["negative particle concentration at gauss point"] = zeros(Float64, Nn * Nrn * case.opt.gsorder, 1)
        ws["positive particle concentration at gauss point"] = zeros(Float64, Np * Nrp * case.opt.gsorder, 1)
        ws["negative particle surface tangential stress at gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["positive particle surface tangential stress at gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["negative particle stress coupling diffusion coefficient"] = zeros(Float64, Nn, 1)
        ws["positive particle stress coupling diffusion coefficient"] = zeros(Float64, Np, 1)
    end

    # ── CZM 键（条件）──
    if case.opt.czm_enabled
        ws["negative electrode cohesive zone damage"] = zeros(Float64, Nn, 1)
        ws["positive electrode cohesive zone damage"] = zeros(Float64, Np, 1)
    end

    return ws
end
```

**数组数量**: mechanicalmodel="full" 时 ~30 个，否则 ~20 个。对比 StandardVariables 的 ~50 个，减少 40-60%。

- [ ] **Step 2: 新增 `SPMe_variables!` 原位变体**

在 `src/SPMe.jl` 中新增接受预分配 workspace 的原位变体。计算逻辑与 `SPMe_variables` (L147-196) **完全一致**，仅入口（不调用 StandardVariables）和状态提取（跳过 ws 不含的键）不同。

```julia
"""
    SPMe_variables!(ws, case, yt, t; I_app, T_e)

SPMe_variables 的原位变体：直接写入预分配 workspace，不创建新 Dict。
计算逻辑与 SPMe_variables 完全一致，仅省去 StandardVariables 分配。
"""
function SPMe_variables!(ws::Dict{String, Union{Array{Float64},Float64}},
                          case::Case, yt, t::Float64;
                          I_app::Union{Nothing,Float64}=nothing,
                          T_e::Union{Nothing,Float64}=nothing)
    param = case.param

    if isnothing(I_app)
        I_app = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    else
        I_app = Float64(I_app)
    end

    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = -I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]

    # 状态提取：仅覆写 ws 中已有的键（跳过 thermal2D 键）
    var_list = collect(keys(case.index))
    if T_e !== nothing
        var_list = filter(k -> k != "temperature", var_list)
    end
    for i in var_list
        if haskey(ws, i)
            val = yt[case.index[i]]
            if isa(ws[i], Array{Float64})
                ws[i][:] = val   # 原位覆盖
            else
                ws[i] = val
            end
        end
    end

    T = T_e === nothing ? yt[case.index["temperature"]] : T_e

    # ── 以下计算逻辑与 SPMe_variables L147-196 逐行一致 ──
    cn_surf = ws["negative particle surface lithium concentration"]
    cp_surf = ws["positive particle surface lithium concentration"]
    ce_n = ws["electrolyte lithium concentration in negative electrode"]
    ce_p = ws["electrolyte lithium concentration in positive electrode"]
    ce_sp = ws["electrolyte lithium concentration in separator"]

    ce_n_gs = sum(mesh_ne.gs.Ni .* ce_n[mesh_ne.element[mesh_ne.gs.ele, :]], dims=2)
    ce_p_gs = sum(mesh_pe.gs.Ni .* ce_p[mesh_pe.element[mesh_pe.gs.ele, :]], dims=2)
    ce_sp_gs = sum(mesh_sp.gs.Ni .* ce_sp[mesh_sp.element[mesh_sp.gs.ele, :]], dims=2)

    j0_n_gs = param.NE.k * Arrhenius(param.NE.Eac_k, T) .* abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5
    j0_p_gs = param.PE.k * Arrhenius(param.PE.Eac_k, T) .* abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5
    j0_n_av = IntV(j0_n_gs, mesh_ne) / param.NE.thickness
    j0_p_av = IntV(j0_p_gs, mesh_pe) / param.PE.thickness
    eta_n = 2.0 * T * asinh.(j_n / 2.0 / j0_n_av)
    eta_p = 2.0 * T * asinh.(j_p / 2.0 / j0_p_av)

    dphi_S = I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)
    kappa_ne_gs = param.EL.kappa(ce_n_gs, T) * param.NE.eps ^ param.NE.brugg
    kappa_pe_gs = param.EL.kappa(ce_p_gs, T) * param.PE.eps ^ param.PE.brugg
    kappa_sp_gs = param.EL.kappa(ce_sp_gs, T) * param.SP.eps ^ param.SP.brugg
    kappa_ne_av = IntV(kappa_ne_gs, mesh_ne) / param.NE.thickness
    kappa_pe_av = IntV(kappa_pe_gs, mesh_pe) / param.PE.thickness
    kappa_sp_av = IntV(kappa_sp_gs, mesh_sp) / param.SP.thickness
    R_EL = param.NE.thickness / (3.0 * kappa_ne_av) + param.SP.thickness / kappa_sp_av + param.PE.thickness / (3.0 * kappa_pe_av)
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    dphi_e = 2.0 * T * (1 - param.EL.tplus) * (csp_av - csn_av) / param.EL.ce0 - I_app * R_EL - dphi_S

    u_n = param.NE.U(cn_surf) .+ (T .- case.param.cell.T0) .* param.NE.dUdT(cn_surf)
    u_p = param.PE.U(cp_surf) .+ (T .- case.param.cell.T0) .* param.PE.dUdT(cp_surf)
    V_cell = u_p - u_n + eta_p - eta_n + dphi_e

    # ── 写入 workspace（与 SPMe_variables L180-200 对应）──
    ws["cell voltage"] = V_cell[1]
    ws["negative electrode exchange current density"] = j0_n_av
    ws["positive electrode exchange current density"] = j0_p_av
    ws["negative electrode interfacial current density"] = j_n
    ws["positive electrode interfacial current density"] = j_p
    ws["negative electrode overpotential"] = eta_n
    ws["positive electrode overpotential"] = eta_p
    ws["negative electrode open circuit potential"] = u_n
    ws["positive electrode open circuit potential"] = u_p
    ws["electrolyte lithium concentration at negative electrode Gauss point"] = ce_n_gs
    ws["electrolyte lithium concentration at positive electrode Gauss point"] = ce_p_gs
    ws["electrolyte lithium concentration at separator Gauss point"] = ce_sp_gs
    ws["time"] = t
    ws["temperature"] = T
    ws["cell current"] = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C

    return ws
end
```

**设计要点**:
- 与 `SPMe_variables` 共享完全一致的计算逻辑（L147-196），仅入口/出口不同
- `haskey(ws, i)` 用于过滤 `case.index` 中 workspace 不含的键（如 thermal2D 键），**非错误处理**
- 保留原始 `SPMe_variables` 不动（main 分支 `SPMe()` 路径仍使用）

- [ ] **Step 3: 修改 `SPMe_element` 接受可选 workspace**

```julia
# 旧签名:
function SPMe_element(case::Case, yt_e::Array{Float64}, t::Float64, e::Int;
                       I_e::Float64, T_e::Float64, jacobi::String="update")

# 新签名（与 Task 5 的 AbstractVector 兼容）:
function SPMe_element(case::Case, yt_e, t::Float64, e::Int;
                       I_e::Float64, T_e::Float64, jacobi::String="update",
                       workspace::Union{Nothing, Dict{String, Union{Array{Float64},Float64}}}=nothing)
    # 1) 使用 workspace 或创建新 variables
    if workspace !== nothing
        variables_e = SPMe_variables!(workspace, case, yt_e, t; I_app=I_e, T_e=T_e)
    else
        yt_e_vec = vec(yt_e)
        variables_e = SPMe_variables(case, yt_e_vec, t; I_app=I_e, T_e=T_e)
    end

    # ... 其余逻辑完全不变 (SPMe.jl L46-92) ...
end
```

**向后兼容**: `workspace=nothing`（默认）时退回原始路径，main 分支无需任何修改。

- [ ] **Step 4: 移除 CallModel.jl 中冗余 `sparse()` 转换**

`blockdiag` 对 `SparseMatrixCSC` 输入已返回 `SparseMatrixCSC`，外层 `sparse()` 是无操作。

```julia
# OLD (L92-93):
M_elems[e] = sparse(M_e)
K_elems[e] = sparse(K_e)

# NEW:
M_elems[e] = M_e
K_elems[e] = K_e
```

- [ ] **Step 5: 修改 `CallModel_MultiSPMe` 使用线程本地精简工作区**

在并行循环前预分配，循环中复用：

```julia
# 在并行循环前（L88 前）新增:
nthreads = Threads.nthreads()
ws_pool = [create_element_workspace(case) for _ in 1:nthreads]

# 替换并行循环（L89-96）:
t_spme_ns = time_ns()
Threads.@threads for e in 1:ne
    tid = Threads.threadid()
    ws_e = ws_pool[tid]
    M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e;
                                           I_e=I_e[e], T_e=Te_prev[e],
                                           jacobi=jacobi, workspace=ws_e)
    M_elems[e] = M_e   # Step 4 已移除 sparse()
    K_elems[e] = K_e
    F_elems[e] = vec(F_e)
    # 拷贝关键结果到独立 Dict（下游消费者需要 per-element 数据）
    variables_elems[e] = copy_element_results(vars_e)
end
t_spme_s = (time_ns() - t_spme_ns) * 1e-9
```

**`copy_element_results` 说明**:

`vars_e` 指向 `ws_e`（workspace 对象），在同一线程内下一个单元会复用 `ws_e` 并覆盖数据。因此必须拷贝关键结果到独立的 per-element Dict，供辅助循环和 `compute_heat_sources` 使用。

> ⚠️ **数据安全性关键**: workspace 中使用原位覆写 (`ws[i][:] = val`) 的键（状态提取键如 "negative particle lithium concentration"）与使用替换赋值 (`ws[key] = computed_value`) 的键行为不同：
> - **原位覆写的键**: 多个单元共享同一数组，必须 `copy()` 提取独立副本
> - **替换赋值的键**: 每次计算产生新数组，`vars_e[key]` 已指向独立数据，引用拷贝即可

```julia
"""
    copy_element_results(vars_e)

从 workspace Dict 中提取下游代码（compute_heat_sources + 辅助循环 + element voltages）
需要的键，返回轻量级独立 Dict。

状态提取键（通过 case.index 原位写入 workspace）必须 copy()，否则同线程下一个单元
覆写 workspace 后，先前单元的 Dict 引用会指向错误数据。
计算结果键（通过 = 赋值，每次创建新数组/标量）可安全引用拷贝。
"""
function copy_element_results(vars_e)
    Dict{String, Any}(
        # ── 计算结果键：每次 = 赋值创建新对象，引用安全 ──
        "negative electrode overpotential"             => vars_e["negative electrode overpotential"],
        "positive electrode overpotential"             => vars_e["positive electrode overpotential"],
        "cell voltage"                                 => vars_e["cell voltage"],
        "negative electrode exchange current density"  => vars_e["negative electrode exchange current density"],
        "positive electrode exchange current density"  => vars_e["positive electrode exchange current density"],
        "negative electrode interfacial current density" => vars_e["negative electrode interfacial current density"],
        "positive electrode interfacial current density" => vars_e["positive electrode interfacial current density"],
        "negative electrode open circuit potential"    => vars_e["negative electrode open circuit potential"],
        "positive electrode open circuit potential"    => vars_e["positive electrode open circuit potential"],
        "temperature"                                  => vars_e["temperature"],
        # ── 状态提取键：原位写入 workspace，必须 copy() ──
        "negative particle surface lithium concentration" => copy(vars_e["negative particle surface lithium concentration"]),
        "positive particle surface lithium concentration" => copy(vars_e["positive particle surface lithium concentration"]),
        "negative particle lithium concentration"      => copy(vars_e["negative particle lithium concentration"]),
        "positive particle lithium concentration"      => copy(vars_e["positive particle lithium concentration"]),
    )
end
```

> **注意**: 4 个 `copy()` 每步产生 4 × 81 = 324 次小数组拷贝（每个数组 10-20 个 Float64），相比消除的 ~3000+ 数组分配是可接受的代价。

- [ ] **Step 6: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: 电压/温度结果一致；SPMe 计时下降 15-25%。

- [ ] **Step 7: Commit**

```bash
git add src/Variables.jl src/SPMe.jl src/CallModel.jl
git commit -m "perf: slim thread-local workspace for SPMe_element (eliminate ~2400 useless allocs/step)"
```

### 预期收益分析

> **修正说明**: `SPMe_variables!` 对计算结果键（Gauss 点浓度、OCV 等）使用 `=` 赋值，每次创建新数组而非覆写预分配数组。因此实际节省低于仅按 StandardVariables 数组数计算的理论值。4 个状态提取键的 `copy()` 也引入 324 次/步小数组分配。

| 分配来源 | 优化前（每步） | 优化后（每步） | 节省 |
|----------|-------------|-------------|------|
| StandardVariables 数组 | 81 × 50 = **4050** | nthreads × 30 ≈ **60**（预分配）+ 81 × ~8 计算结果新数组 ≈ **648** | **~3342** |
| copy_element_results | 0 | 81 × 4 copy + 81 Dict ≈ **405** | -405（新增） |
| 冗余 sparse() 转换 | 160 | **0**（Step 4 移除） | **160** |
| Assemble KI/KJ/KV | 480 × 3 = 1440 | 1440（不变，见 Task 8） | 0 |
| **合计堆分配** | **~5650** | **~2553** | **~3100 (55%)** |

核心收益: 消除 4050 次无用的 StandardVariables 数组分配（其中 2400 次为 thermal2D 键），扣除计算结果新数组后净减 ~3400。Assemble 分配需由 Task 8 单独解决。

---

## Chunk 5: 低优先级优化（可延后）

### Task 8: 预分配 Assemble 缓冲区（P3）

`Assemble` 每次调用分配 KI/KJ/KV。可新增 `Assemble!` 预分配变体在热模型中使用。

> **注意**: 热模型仅占 7.85%，此优化收益约 0.3%。依赖 Chunk 2 的类型安全缓存方案（不用 Dict）。

**方案**: 在 `ThermalDistributed2D` 入口分配一次 KI/KJ/KV 缓冲区，内部多次 `Assemble!` 调用复用。缓冲区可在 `MeshGeometry` 新增字段或作为函数局部变量在首次调用时缓存。

- [ ] 待 Chunk 1-3 完成后根据实际 profiling 决定是否实施。

### Task 9: Variable_update! 优化（P3）

`Variable_update!` 每步遍历所有键两次。early-exit 扩展检查可减少一次遍历。

> **注意**: `Variable_update!` 不在热点路径中（91.59% SPMe + 7.85% Thermal 之外），实际收益 <0.1%。

- [ ] 仅当 profiling 显示此函数耗时显著时才实施。

---

## ~~已否决的任务~~

### ~~Task 0: 添加 thermal_extras 字段~~ — 否决

原计划要求恢复 `case.thermal_extras::Dict{String,Any}` 字段。该字段在结构体重构中已刻意移除（见 `00_优化分析进度.md`: "**消除 thermal_extras Dict**"）。恢复会违反 Fail-Fast 设计原则。所有缓存需求已改为使用类型安全的 struct 字段。

### ~~Task 5 (原): 消除代表性状态计算~~ — 降级为可选

原计划提议用 element 1 的状态代替平均值作为代表性状态。理论依据薄弱 — 分流求解器的 `compute_prefactors` 依赖代表性状态计算初始电压估计，element 1 可能位于温度/电流分布的极端位置。需在温度梯度大的场景下验证。**不推荐在核心路径实施**。

---

## Expected Results Summary

| 优先级 | Task | Target | Expected Savings | Risk |
|--------|------|--------|-----------------|------|
| **P0** | Task 1: Solve.jl deepcopy → copy | 消除 7 处不必要 deepcopy | ~2-3% | 低 |
| **P0** | Task 2: Assemble deepcopy(KI) | 消除 1 处 deepcopy | ~0.1% | 低 |
| **P0** | Task 3: 缓存单元面积 | 消除每步面积重算 | ~0.5% | 低 |
| **P1** | Task 4: 缓存边界边 | 消除每步 Set 重建 | ~0.3% | 中 |
| **P1** | Task 5: @views 状态提取 | 消除 ne 次拷贝 | ~1% | 低 |
| **P1** | Task 6: 原位热 BC | 消除三重 copy | ~0.5% | 中 |
| **P2** | Task 7: 精简型线程工作区 | 消除 ~3100 allocs/step (排除 ~30 无用 thermal2D 键 + 对象池复用，含 324 copy + 648 计算新数组) | **~12-20%** | 中 |
| P3 | Task 8: 预分配 Assemble | 减少热组装分配 | ~0.3% | 中 |
| P3 | Task 9: Variable_update! | 优化记录路径 | <0.1% | 低 |

**Chunk 1-3 保守估计: 5-7% wall-clock reduction**
**含 Chunk 4 (Task 7): 17-27% wall-clock reduction**

**Verification:** 每完成一个 Task 后运行 `example/testexample.jl` 对比:
1. 电压曲线 (误差 < 1e-6 V)
2. 温度演化 (误差 < 1e-4 K)
3. 计时分解 (SPMe 比例应从 91.59% 逐步下降)
