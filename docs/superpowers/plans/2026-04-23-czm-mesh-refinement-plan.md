# CZM Mesh Refinement Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three-mode CZM mesh refinement (default/multiple/auto) while keeping bulk mesh unchanged, with correct N:1 coupling.

**Architecture:** New `CZMMesh.jl` absorbs mesh generation from `czm.jl`. `czm.jl` retains only mechanics assembly. Refinement inserts intermediate CZM nodes and sub-elements, maps them back to thermal nodes for gap conductance, and rebuilds `czm_element_map` with N:1 entries.

**Tech Stack:** Julia, existing JuBat module system, SparseArrays.

**Spec:** `docs/superpowers/specs/2026-04-23-czm-mesh-refinement-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `src/Option.jl` | Modify | Add `czm_mesh_refine`, `czm_refine_factor` fields |
| `src/SetMesh.jl` | Modify | Change `cohesive_elements` to `Vector{CohesiveElement}`, `damage_states` to `Vector{DamageState}` |
| `src/CZMMesh.jl` | Create | Structs (`CohesiveElement`, `DamageState`), mesh generation, refinement, map rebuild |
| `src/czm.jl` | Modify | Remove migrated structs and `create_czm_mesh`; keep assembly only |
| `src/JuBat.jl` | Modify | Add `include("CZMMesh.jl")`, update exports |
| `src/Materialmatrix.jl` | Modify | AND/OR logic dispatch in `get_active_elements` |
| `src/ThermalDistributed.jl` | Modify | Use `thermal_node_bottom/top` in BC gap conductance |
| `src/CouplingState.jl` | Modify | Pass `opt` to `get_active_elements` for mode dispatch |

---

## Chunk 1: Foundation — Option, Structs, File Restructuring

### Task 1: Add Option fields

**Files:**
- Modify: `src/Option.jl:72-80`

- [ ] **Step 1: Add two fields to Option struct**

In `src/Option.jl`, after the existing CZM fields (after line 80 `czm_arc_length_alpha`), add:

```julia
    czm_mesh_refine::String = "default"   # "default" | "multiple" | "auto"
    czm_refine_factor::Int64 = 2          # "multiple" mode: sub-elements per original CZM segment
```

- [ ] **Step 2: Verify syntax**

Run: `julia -e 'include("src/Option.jl"); using Parameters; opt = Option(); println(opt.czm_mesh_refine, " ", opt.czm_refine_factor)'`
Expected: `default 2`

- [ ] **Step 3: Commit**

```bash
git add src/Option.jl
git commit -m "feat: add czm_mesh_refine and czm_refine_factor to Option"
```

---

### Task 2: Create CZMMesh.jl with migrated structs

**Files:**
- Create: `src/CZMMesh.jl`
- Modify: `src/czm.jl` (remove migrated code)
- Modify: `src/JuBat.jl` (include + exports)

- [ ] **Step 1: Create CZMMesh.jl with CohesiveElement (extended fields) and DamageState**

Create `src/CZMMesh.jl` containing:

```julia
# CZMMesh.jl — CZM mesh generation and refinement
# Structs migrated from czm.jl; mesh refinement functions added.

mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    nodes::Vector{Int64}              # [n1, n2, n3, n4]
    nodes_bottom::Vector{Int64}       # [n1, n2] 底面节点
    nodes_top::Vector{Int64}          # [n4, n3] 顶面节点
    length::Float64                   # 单元长度
    layer_idx::Int64                  # 层间界面索引
    # Refinement fields
    thermal_node_bottom::Vector{Int64}  # thermal mesh nodes for gap conductance
    thermal_node_top::Vector{Int64}     # thermal mesh nodes for gap conductance
    parent_czm_idx::Int64               # parent CZM element index (self if unrefined)
end

"""
    DamageState - 内聚力单元的损伤状态
"""
mutable struct DamageState <: AbstractDamageState
    D::Float64
    δ_max_n::Float64
    δ_max_t::Float64
    δ_max_eff::Float64
    fractured::Bool
    accumulated_damage::Float64
    DamageState() = new(0.0, 0.0, 0.0, 0.0, false, 0.0)
end
```

- [ ] **Step 2: Migrate create_czm_mesh from czm.jl**

Copy the entire `create_czm_mesh` function (lines 54–134 of current `src/czm.jl`) into `CZMMesh.jl`.
Modify the `CohesiveElement` constructor call (around line 98) to include the three new fields:

```julia
        coh_elem = CohesiveElement(
            i,
            [n_in_1, n_in_2, n_out_2, n_out_1],
            [n_in_1, n_in_2],
            [n_out_1, n_out_2],
            elem_length,
            1,
            # Refinement fields (default = self)
            [n_in_1, n_in_2],      # thermal_node_bottom = nodes_bottom
            [n_out_1, n_out_2],    # thermal_node_top = nodes_top
            i                       # parent_czm_idx = self
        )
```

- [ ] **Step 3: Remove migrated code from czm.jl**

Delete from `src/czm.jl`:
- The `CohesiveElement` struct definition (lines 1–8)
- The `DamageState` struct definition (lines 23–33)
- The `create_czm_mesh` function (lines 54–134)

Keep everything else in `czm.jl` (assembly, cache, BC functions).

- [ ] **Step 4: Update JuBat.jl include order**

In `src/JuBat.jl`, add `include("CZMMesh.jl")` before `include("czm.jl")`:

```julia
include("CZMMesh.jl")  # CZM mesh generation and refinement
include("czm.jl")      # CZM mechanics assembly (structs removed)
```

- [ ] **Step 5: Update JuBat.jl exports**

The existing exports (`CohesiveElement`, `DamageState`, `create_czm_mesh`) remain valid since they now come from `CZMMesh.jl`. Add new exports:

```julia
export create_czm_mesh_refined, compute_cohesive_zone_length
export refine_czm_mesh_multiple!, refine_czm_mesh_auto!
export rebuild_czm_element_map
```

- [ ] **Step 6: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK` (no errors)

- [ ] **Step 7: Commit**

```bash
git add src/CZMMesh.jl src/czm.jl src/JuBat.jl
git commit -m "refactor: create CZMMesh.jl, migrate structs and create_czm_mesh from czm.jl"
```

---

### Task 3: Update CohesiveMesh container types in SetMesh.jl

**Files:**
- Modify: `src/SetMesh.jl:31,36`

- [ ] **Step 1: Change container types**

In `src/SetMesh.jl`, change `CohesiveMesh` struct (lines 26–46):

Line 31: `cohesive_elements::Vector{AbstractCohesiveElement}` → `cohesive_elements::Vector{CohesiveElement}`
Line 36: `damage_states::Vector{AbstractDamageState}` → `damage_states::Vector{DamageState}`

Also update the inner constructor (lines 39–45) to use the concrete types:

```julia
    function CohesiveMesh()
        new(Mesh("Q4", 2, zeros(0,2), 0, zeros(Int64,0,4),
            GaussPoint(zeros(0,2), zeros(0,2), zeros(0), zeros(0), zeros(Int64,0), zeros(0,4), zeros(0,8), 2)),
            zeros(0, 2), 0, zeros(Int64, 0, 4),
            CohesiveElement[], 0, 0, Dict{Int64, Vector{Int64}}(),
            Vector{Vector{Tuple{Int64,Int64}}}(), DamageState[])
    end
```

- [ ] **Step 2: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/SetMesh.jl
git commit -m "refactor: use concrete CohesiveElement/DamageState in CohesiveMesh"
```

---

## Chunk 2: Refinement Core — multiple and auto modes

### Task 4: Implement refine_czm_mesh_multiple!

**Files:**
- Modify: `src/CZMMesh.jl`

- [ ] **Step 1: Implement the function**

Add to `src/CZMMesh.jl`:

```julia
"""
    refine_czm_mesh_multiple!(czm_mesh, n)

将每个原始 CZM 单元等分为 n 个子单元。就地修改 czm_mesh。

# 参数
- `czm_mesh`: 原始 1:1 CohesiveMesh（由 create_czm_mesh 生成）
- `n::Int64`: 每个原始单元的等分数
"""
function refine_czm_mesh_multiple!(czm_mesh::CohesiveMesh, n::Int64)
    n == 1 && return czm_mesh  # 无需加密

    orig_elements = czm_mesh.cohesive_elements
    n_orig = length(orig_elements)
    nnode_thermal = size(czm_mesh.node, 1)

    new_elements = CohesiveElement[]
    new_damage = DamageState[]
    new_nodes = copy(czm_mesh.node)  # 复制现有节点

    next_node_id = nnode_thermal + 1

    for (orig_idx, orig) in enumerate(orig_elements)
        n_b1, n_b2 = orig.nodes_bottom
        n_t1, n_t2 = orig.nodes_top
        L = orig.length
        sub_L = L / n

        # 父单元的热节点（用于间隙热导映射）
        thermal_nb = orig.thermal_node_bottom
        thermal_nt = orig.thermal_node_top

        # 边缘节点序列（bottom: n_b1 → n_b2, top: n_t1 → n_t2）
        bottom_nodes = [n_b1]
        top_nodes = [n_t1]

        for k in 1:(n - 1)
            α = k / n
            # Bottom edge intermediate node
            xb = (1 - α) * czm_mesh.node[n_b1, 1] + α * czm_mesh.node[n_b2, 1]
            yb = (1 - α) * czm_mesh.node[n_b1, 2] + α * czm_mesh.node[n_b2, 2]
            new_nodes = [new_nodes; [xb yb]]
            push!(bottom_nodes, next_node_id)
            next_node_id += 1

            # Top edge intermediate node
            xt = (1 - α) * czm_mesh.node[n_t1, 1] + α * czm_mesh.node[n_t2, 1]
            yt = (1 - α) * czm_mesh.node[n_t1, 2] + α * czm_mesh.node[n_t2, 2]
            new_nodes = [new_nodes; [xt yt]]
            push!(top_nodes, next_node_id)
            next_node_id += 1
        end
        push!(bottom_nodes, n_b2)
        push!(top_nodes, n_t2)

        # 创建子单元
        for k in 1:n
            sub_id = (orig_idx - 1) * n + k
            sub_elem = CohesiveElement(
                sub_id,
                [bottom_nodes[k], bottom_nodes[k + 1], top_nodes[k + 1], top_nodes[k]],
                [bottom_nodes[k], bottom_nodes[k + 1]],
                [top_nodes[k], top_nodes[k + 1]],
                sub_L,
                orig.layer_idx,
                thermal_nb,   # 所有子单元映射到父单元的热节点
                thermal_nt,
                orig_idx      # 父单元索引
            )
            push!(new_elements, sub_elem)
            push!(new_damage, DamageState())
        end
    end

    # 更新 czm_mesh
    czm_mesh.node = new_nodes
    czm_mesh.nnode = size(new_nodes, 1)
    czm_mesh.cohesive_elements = new_elements
    czm_mesh.n_cohesive = length(new_elements)
    czm_mesh.damage_states = new_damage

    return czm_mesh
end
```

- [ ] **Step 2: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/CZMMesh.jl
git commit -m "feat: add refine_czm_mesh_multiple! for fixed multiplier refinement"
```

---

### Task 5: Implement compute_cohesive_zone_length and refine_czm_mesh_auto!

**Files:**
- Modify: `src/CZMMesh.jl`

- [ ] **Step 1: Implement compute_cohesive_zone_length**

Add to `src/CZMMesh.jl`:

```julia
"""
    compute_cohesive_zone_length(param)

计算内聚力区特征长度 lc = Gc · E_eff / σ_max²（归一化量）。

# 参数
- `param`: 归一化参数对象（case.param）
"""
function compute_cohesive_zone_length(param)
    Gc = param.cohesive.G_c_n
    σ_max = param.cohesive.σ_max_n
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) /
            (param.NE.thickness + param.PE.thickness)
    return Gc * E_eff / σ_max^2
end
```

- [ ] **Step 2: Implement refine_czm_mesh_auto!**

Add to `src/CZMMesh.jl`:

```julia
"""
    refine_czm_mesh_auto!(czm_mesh, param)

自适应加密：每个原始 CZM 单元等分为 n_i 个子单元，使子单元长度 ≤ lc/3。
n_i = max(1, ceil(L_i / (lc / 3)))
"""
function refine_czm_mesh_auto!(czm_mesh::CohesiveMesh, param)
    lc = compute_cohesive_zone_length(param)
    target_length = lc / 3.0

    orig_elements = czm_mesh.cohesive_elements
    nnode_thermal = size(czm_mesh.node, 1)

    new_elements = CohesiveElement[]
    new_damage = DamageState[]
    new_nodes = copy(czm_mesh.node)

    next_node_id = nnode_thermal + 1

    for (orig_idx, orig) in enumerate(orig_elements)
        L = orig.length
        n_i = max(1, ceil(Int, L / target_length))
        sub_L = L / n_i

        n_b1, n_b2 = orig.nodes_bottom
        n_t1, n_t2 = orig.nodes_top
        thermal_nb = orig.thermal_node_bottom
        thermal_nt = orig.thermal_node_top

        bottom_nodes = [n_b1]
        top_nodes = [n_t1]

        for k in 1:(n_i - 1)
            α = k / n_i
            xb = (1 - α) * czm_mesh.node[n_b1, 1] + α * czm_mesh.node[n_b2, 1]
            yb = (1 - α) * czm_mesh.node[n_b1, 2] + α * czm_mesh.node[n_b2, 2]
            new_nodes = [new_nodes; [xb yb]]
            push!(bottom_nodes, next_node_id)
            next_node_id += 1

            xt = (1 - α) * czm_mesh.node[n_t1, 1] + α * czm_mesh.node[n_t2, 1]
            yt = (1 - α) * czm_mesh.node[n_t1, 2] + α * czm_mesh.node[n_t2, 2]
            new_nodes = [new_nodes; [xt yt]]
            push!(top_nodes, next_node_id)
            next_node_id += 1
        end
        push!(bottom_nodes, n_b2)
        push!(top_nodes, n_t2)

        sub_id_offset = length(new_elements)
        for k in 1:n_i
            sub_elem = CohesiveElement(
                sub_id_offset + k,
                [bottom_nodes[k], bottom_nodes[k + 1], top_nodes[k + 1], top_nodes[k]],
                [bottom_nodes[k], bottom_nodes[k + 1]],
                [top_nodes[k], top_nodes[k + 1]],
                sub_L,
                orig.layer_idx,
                thermal_nb,
                thermal_nt,
                orig_idx
            )
            push!(new_elements, sub_elem)
            push!(new_damage, DamageState())
        end
    end

    czm_mesh.node = new_nodes
    czm_mesh.nnode = size(new_nodes, 1)
    czm_mesh.cohesive_elements = new_elements
    czm_mesh.n_cohesive = length(new_elements)
    czm_mesh.damage_states = new_damage

    return czm_mesh
end
```

- [ ] **Step 3: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add src/CZMMesh.jl
git commit -m "feat: add compute_cohesive_zone_length and refine_czm_mesh_auto!"
```

---

### Task 6: Implement create_czm_mesh_refined and rebuild_czm_element_map

**Files:**
- Modify: `src/CZMMesh.jl`

- [ ] **Step 1: Implement unified entry point**

Add to `src/CZMMesh.jl`:

```julia
"""
    create_czm_mesh_refined(thermal_mesh, param_dim, opt; tol=1e-8)

统一 CZM 网格生成入口。根据 opt.czm_mesh_refine 选择模式：
- "default": 不加密（当前 1:1 实现）
- "multiple": 固定倍数加密
- "auto": 自适应加密（基于 lc）
"""
function create_czm_mesh_refined(thermal_mesh::Mesh, param_dim, opt; tol::Float64=1e-8)
    # Step 1: 生成基础 1:1 CZM 网格
    czm_mesh = create_czm_mesh(thermal_mesh, param_dim; tol=tol)

    mode = opt.czm_mesh_refine

    if mode == "default"
        # 不加密，直接返回
        return czm_mesh
    elseif mode == "multiple"
        refine_czm_mesh_multiple!(czm_mesh, opt.czm_refine_factor)
    elseif mode == "auto"
        # 需要归一化参数（param_dim 尚未归一化时需要 case.param）
        # 此处先检查 param_dim 是否有归一化后的参数
        error("auto mode requires normalized params — use create_czm_mesh_refined(case) variant")
    else
        error("unknown czm_mesh_refine mode: $mode (expected: default, multiple, auto)")
    end

    return czm_mesh
end

"""
    create_czm_mesh_refined(case::Case; tol=1e-8)

从 Case 对象创建加密 CZM 网格（支持 auto 模式，需要归一化参数）。
"""
function create_czm_mesh_refined(case::Case; tol::Float64=1e-8)
    thermal_mesh = case.mesh["thermal2D"]
    czm_mesh = create_czm_mesh(thermal_mesh, case.param_dim; tol=tol)

    mode = case.opt.czm_mesh_refine

    if mode == "default"
        return czm_mesh
    elseif mode == "multiple"
        refine_czm_mesh_multiple!(czm_mesh, case.opt.czm_refine_factor)
    elseif mode == "auto"
        refine_czm_mesh_auto!(czm_mesh, case.param)
    else
        error("unknown czm_mesh_refine mode: $mode")
    end

    return czm_mesh
end
```

- [ ] **Step 2: Implement rebuild_czm_element_map**

Add to `src/CZMMesh.jl`:

```julia
"""
    rebuild_czm_element_map(czm_mesh, original_czm_element_map)

基于 parent_czm_idx 从原始 1:1 czm_element_map 重建 N:1 映射。
O(N) 复杂度。
"""
function rebuild_czm_element_map(czm_mesh::CohesiveMesh, original_czm_element_map::Dict{Int, Vector{Int}})
    # 反向索引: parent_czm_idx → [thermal_elements...]
    parent_to_thermal = Dict{Int, Vector{Int}}()
    for (e, czm_list) in original_czm_element_map
        for czm_idx in czm_list
            if !haskey(parent_to_thermal, czm_idx)
                parent_to_thermal[czm_idx] = Int[]
            end
            push!(parent_to_thermal[czm_idx], e)
        end
    end

    new_map = Dict{Int, Vector{Int}}()
    for e in keys(original_czm_element_map)
        new_map[e] = Int[]
    end

    for (idx, elem) in enumerate(czm_mesh.cohesive_elements)
        for e in get(parent_to_thermal, elem.parent_czm_idx, Int[])
            push!(new_map[e], idx)
        end
    end

    return new_map
end
```

- [ ] **Step 3: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add src/CZMMesh.jl
git commit -m "feat: add create_czm_mesh_refined and rebuild_czm_element_map"
```

---

## Chunk 3: Integration — Coupling Updates

### Task 7: Update ThermalDistributed2D_BC for gap conductance

**Files:**
- Modify: `src/ThermalDistributed.jl:292-309`

- [ ] **Step 1: Replace nodes_bottom/top with thermal_node_bottom/top**

In `ThermalDistributed2D_BC` (around line 295–308), change:

```julia
        n_bot = czm_elem.nodes_bottom
        n_top = czm_elem.nodes_top
```

to:

```julia
        n_bot = czm_elem.thermal_node_bottom
        n_top = czm_elem.thermal_node_top
```

The rest of the loop (`for (nb, nt) in zip(n_bot, n_top)`, matrix assembly) stays unchanged.

- [ ] **Step 2: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/ThermalDistributed.jl
git commit -m "refactor: use thermal_node_bottom/top for gap conductance in BC"
```

---

### Task 8: Update get_active_elements with AND/OR dispatch

**Files:**
- Modify: `src/Materialmatrix.jl:345-363`
- Modify: `src/CouplingState.jl` (if get_active_elements signature needs opt)

- [ ] **Step 1: Add opt parameter to get_active_elements**

Change signature in `src/Materialmatrix.jl`:

```julia
function get_active_elements(czm_mesh::CohesiveMesh, mesh_data::MeshGeometry; czm_refine_mode::String = "default")
```

- [ ] **Step 2: Implement mode dispatch**

Replace the body (lines 346–363) with:

```julia
function get_active_elements(czm_mesh::CohesiveMesh, mesh_data::MeshGeometry; czm_refine_mode::String = "default")
    ne = length(mesh_data.element_layer)
    active = ones(Bool, ne)
    fractured_czm = get_fractured_elements(czm_mesh)

    if czm_refine_mode == "default"
        # OR logic: 任一 CZM 断裂 → 热单元停止（当前行为）
        for e in 1:ne
            if !mesh_data.is_inner_layer[e]
                continue
            end
            for czm_idx in get(mesh_data.czm_element_map, e, Int64[])
                if czm_idx in fractured_czm
                    active[e] = false
                    break
                end
            end
        end
    else
        # AND logic: 所有 CZM 都断裂 → 热单元停止（加密模式）
        for e in 1:ne
            if !mesh_data.is_inner_layer[e]
                continue
            end
            czm_list = get(mesh_data.czm_element_map, e, Int64[])
            if !isempty(czm_list) && all(c -> c in fractured_czm, czm_list)
                active[e] = false
            end
        end
    end

    return findall(active)
end
```

- [ ] **Step 3: Update caller in compute_heat_sources_with_czm**

In `src/ThermalDistributed.jl`, around line 519, change:

```julia
        active_elements = get_active_elements(czm_mesh, geom)
```

to:

```julia
        active_elements = get_active_elements(czm_mesh, geom; czm_refine_mode=case.opt.czm_mesh_refine)
```

- [ ] **Step 4: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add src/Materialmatrix.jl src/ThermalDistributed.jl
git commit -m "feat: add AND/OR dispatch in get_active_elements based on czm_mesh_refine mode"
```

---

### Task 9: Wire create_czm_mesh_refined into setup_thermal2D_mesh

**Files:**
- Modify: `src/Jellyrollmodel.jl` (where czm_mesh is created)

- [ ] **Step 1: Find and update CZM mesh creation call site**

In `src/Jellyrollmodel.jl`, locate where `create_czm_mesh` is called within `setup_thermal2D_mesh`. Replace:

```julia
case.czm_mesh = create_czm_mesh(mesh_th, case.param_dim)
```

with:

```julia
case.czm_mesh = create_czm_mesh_refined(case)
```

This uses the Case-based overload which has access to both `case.opt` (for mode) and `case.param` (for auto mode).

- [ ] **Step 2: Rebuild czm_element_map after refinement**

After the `create_czm_mesh_refined` call, if the mode is not `"default"`, rebuild the map:

```julia
if case.opt.czm_mesh_refine != "default"
    new_map = rebuild_czm_element_map(case.czm_mesh, mesh_data.czm_element_map)
    # Update case.geometry with the new N:1 map
    if case.geometry !== nothing
        case.geometry = MeshGeometry(
            case.geometry.element_layer,
            case.geometry.is_inner_layer,
            case.geometry.layer_weights,
            case.geometry.interface_pairs,
            new_map,  # 替换为 N:1 map
            case.geometry.inner_nodes,
            case.geometry.outer_nodes,
            case.geometry.boundary_edges
        )
    end
end
```

- [ ] **Step 3: Verify module loads**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add src/Jellyrollmodel.jl
git commit -m "feat: wire create_czm_mesh_refined into setup_thermal2D_mesh with map rebuild"
```

---

## Chunk 4: Verification

### Task 10: Smoke test — default mode unchanged

**Files:**
- Test: manual verification

- [ ] **Step 1: Run existing CZM example with default mode**

Run an existing CZM test/case with `opt.czm_mesh_refine = "default"` (the default). Verify:
- Module loads without error
- CZM mesh has same number of elements as before
- Simulation results match previous output

- [ ] **Step 2: Commit (if any fixes needed)**

```bash
git add -A
git commit -m "fix: address issues from default mode smoke test"
```

---

### Task 11: Smoke test — multiple mode

**Files:**
- Test: manual verification

- [ ] **Step 1: Run CZM example with multiple mode**

Set `opt.czm_mesh_refine = "multiple"`, `opt.czm_refine_factor = 2`.
Verify:
- CZM mesh element count = 2 × original count
- CZM mesh node count = original + (2-1) × 2 × n_orig
- Each sub-element has `thermal_node_bottom/top` = parent's thermal nodes
- Each sub-element has `parent_czm_idx` pointing to the correct original element
- Simulation runs without error

- [ ] **Step 2: Commit (if any fixes needed)**

```bash
git add -A
git commit -m "fix: address issues from multiple mode smoke test"
```

---

### Task 12: Smoke test — auto mode

**Files:**
- Test: manual verification

- [ ] **Step 1: Run CZM example with auto mode**

Set `opt.czm_mesh_refine = "auto"`.
Verify:
- `compute_cohesive_zone_length` returns a reasonable value
- All sub-element lengths ≤ lc/3
- Minimum n_i = 1 for already-small elements
- Simulation runs without error

- [ ] **Step 2: Commit (if any fixes needed)**

```bash
git add -A
git commit -m "fix: address issues from auto mode smoke test"
```
