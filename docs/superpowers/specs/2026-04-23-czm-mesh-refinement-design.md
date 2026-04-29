# CZM Mesh Refinement Design

**Date**: 2026-04-23
**Status**: Draft

## Problem

Current CZM mesh generation (`create_czm_mesh`) creates one CZM element per thermal element (1:1 mapping). This limits resolution of the cohesive zone — a single CZM element spans the entire thermal element face, which may be much larger than the cohesive zone length `lc = Gc·E / σ_max²`. Without sufficient CZM elements within `lc`, damage gradients cannot be resolved accurately.

## Goal

Add CZM mesh refinement while keeping the bulk (thermal) mesh unchanged. Support three refinement modes via Option fields. Correctly handle the resulting N:1 mapping (multiple CZM elements per thermal element) in electro-chemical deactivation, gap conductance, and the `czm_element_map`.

## Constraints

- Bulk mesh nodes and topology must not change.
- Refined CZM edge nodes must coincide exactly with thermal mesh nodes.
- Refinement is confined within a single original CZM segment — no cross-element refinement.
- `lc/3` convergence criterion is hardcoded, not a user parameter.

## Design

### 1. Option Fields (`src/Option.jl`)

Add two fields to the `Option` struct:

```julia
czm_mesh_refine::String = "default"   # "default" | "multiple" | "auto"
czm_refine_factor::Int64 = 2          # fixed multiplier for "multiple" mode
```

### 2. Three Refinement Modes

#### `default` — no refinement (current behavior)

1:1 mapping. `create_czm_mesh` runs unchanged. One CZM element per thermal element face.

#### `multiple` — fixed multiplier

1. Run `create_czm_mesh` to generate the base 1:1 CZM mesh.
2. For each original CZM element, insert `czm_refine_factor - 1` equally-spaced intermediate nodes on both bottom and top edges.
3. Create `czm_refine_factor` sub-elements, each with `length = original_length / n`.
4. Example: `czm_refine_factor = 3`, 80 original CZM → 240 refined CZM.

#### `auto` — adaptive refinement based on cohesive zone length

1. Run `create_czm_mesh` to generate the base 1:1 CZM mesh.
2. Compute `lc = Gc · E_eff / σ_max²` using normalized parameters.
3. For each original CZM element of length `L_i`:
   ```
   n_i = max(1, ceil(L_i / (lc / 3)))
   ```
4. Refine each element into `n_i` sub-elements (same procedure as `multiple`).
5. Different elements may get different `n_i`.

**`lc` parameter sources (all normalized):**
- `Gc = param.cohesive.G_c_n`
- `E_eff = (E_NE·t_NE + E_PE·t_PE) / (t_NE + t_PE)` (already computed by `compute_czm_effective_params`)
- `σ_max = param.cohesive.σ_max_n`

### 3. Node Strategy

- **Edge nodes**: reuse thermal mesh node IDs (no new nodes).
- **Intermediate nodes**: appended to the CZM node table, numbered from `nnode_thermal + 1`. Coordinates are linearly interpolated along bottom/top edges.
- **Intermediate nodes belong to CZM only** — they do not participate in the thermal system.

**Node count formula**: `nnode_czm = nnode_thermal + Σ(n_i - 1) × 2` (each sub-segment inserts `n-1` nodes, one on each edge).

**C3 修复**：细化函数必须就地更新 `czm_mesh.nnode` 和 `czm_mesh.node`：
- `czm_mesh.nnode = nnode_czm`
- `czm_mesh.node = [czm_mesh.node; new_node_coords]`（纵向扩展矩阵，新行对应中间节点坐标）

这确保下游 `CzmLayout(czm_mesh)` 计算 `ndof = 2 * nnode` 时得到正确维度，`assemble_czm_system` 中稀疏矩阵尺寸也正确。`czm_mesh.bulk_mesh` 保持不变（仍指向原始热网格，其 `nlen` 不变）。

### 4. File Structure

Create `src/CZMMesh.jl` — contains all mesh generation and refinement logic:

| Content | Source |
|---------|--------|
| `CohesiveElement` struct | Migrated from `czm.jl` (add `thermal_node_bottom/top`, `parent_czm_idx` fields) |
| `DamageState` struct | Migrated from `czm.jl` |
| `create_czm_mesh` | Migrated from `czm.jl` (base 1:1 generation) |
| `compute_cohesive_zone_length` | New — computes `lc` from normalized parameters |
| `refine_czm_mesh_multiple!` | New — fixed multiplier refinement |
| `refine_czm_mesh_auto!` | New — adaptive refinement |
| `create_czm_mesh_refined` | New — unified entry point, dispatches by `opt.czm_mesh_refine` |
| `rebuild_czm_element_map` | New — rebuilds N:1 map from `parent_czm_idx` |

`src/czm.jl` retains only mechanics assembly functions:

- `assemble_czm_system`, `assemble_bulk_stiffness`, `assemble_thermal_chemical_load`
- `build_czm_cache`, `ensure_czm_cache`
- `assemble_coupled_system`, `assemble_coupled_system_full`
- `apply_bc_czm`, `identify_bc_nodes_czm`

`src/JuBat.jl` update:
- Add `include("CZMMesh.jl")` before `include("czm.jl")`
- Remove `CohesiveElement`, `DamageState`, `create_czm_mesh` exports from `czm.jl` section
- Add new exports for `CZMMesh.jl` functions

### 5. CohesiveElement Extended Fields

**C4 修复**：`CohesiveElement` 定义从 `czm.jl` 迁移到 `CZMMesh.jl`。`CohesiveMesh.cohesive_elements` 的容器类型从 `Vector{AbstractCohesiveElement}` 改为 `Vector{CohesiveElement}`（在 `src/SetMesh.jl` 中修改），以支持直接字段访问新字段，无需通过 abstract 类型动态派发。

```julia
mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    nodes::Vector{Int64}              # [n1, n2, n3, n4] — may include intermediate nodes
    nodes_bottom::Vector{Int64}       # [n1, n2]
    nodes_top::Vector{Int64}          # [n4, n3]
    length::Float64
    layer_idx::Int64
    # --- New fields for refinement ---
    thermal_node_bottom::Vector{Int64}  # corresponding thermal mesh nodes (for gap conductance)
    thermal_node_top::Vector{Int64}     # corresponding thermal mesh nodes (for gap conductance)
    parent_czm_idx::Int64               # parent CZM element index (self for unrefined)
end
```

For `default` mode: `thermal_node_* = nodes_*`, `parent_czm_idx = id`.
For refined modes: `thermal_node_*` = parent's edge nodes (which are thermal mesh nodes), `parent_czm_idx` = original CZM index.

### 6. Electro-Chemical Deactivation (AND Logic)

**C1 修复**：AND 逻辑仅在启用加密（`czm_mesh_refine != "default"`）时生效。`default` 模式保持当前 OR 逻辑不变，避免静默改变现有仿真结果。

Modify `get_active_elements` in `src/Materialmatrix.jl`:

**`default` mode** — keep current OR logic unchanged:
```julia
# any-one-fractured → deactivate (current behavior, unchanged)
for czm_idx in get(mesh_data.czm_element_map, e, Int64[])
    if czm_idx in fractured_czm
        active[e] = false
        break
    end
end
```

**`multiple` / `auto` mode** — AND logic:
```julia
# Initialize: active = ones(Bool, ne) (default active)
# All-must-be-fractured → deactivate
czm_list = get(mesh_data.czm_element_map, e, Int64[])
if !isempty(czm_list) && all(c -> c in fractured_czm, czm_list)
    active[e] = false
end
```

**C2 修复**：保持 `active = ones(Bool, ne)` 初始化不变（默认活跃），避免空映射时误杀热单元。仅当 `czm_list` 非空且全部断裂时才设为 `false`。

### 7. Gap Conductance Assembly

Modify `ThermalDistributed2D_BC` in `src/ThermalDistributed.jl` (lines 292-309):

```julia
for (elem_idx, czm_elem) in enumerate(czm_mesh.cohesive_elements)
    state = czm_mesh.damage_states[elem_idx]
    D = state.D
    δ_n = state.δ_max_n
    h_eff_nd = compute_gap_conductance(D, δ_n, param.cohesive)
    coeff = h_eff_nd * czm_elem.length
    # Use thermal_node_bottom/top instead of nodes_bottom/top
    for (nb, nt) in zip(czm_elem.thermal_node_bottom, czm_elem.thermal_node_top)
        K[nb, nb] -= coeff
        K[nb, nt] += coeff
        K[nt, nb] += coeff
        K[nt, nt] -= coeff
    end
end
```

Multiple refined CZM elements contribute to the same thermal node pair — their `coeff × length` terms sum naturally.

### 8. `czm_element_map` Rebuild

In `CZMMesh.jl`, after refinement:

```julia
function rebuild_czm_element_map(czm_mesh, original_czm_element_map)
    # Build reverse index: parent_czm_idx → thermal elements (O(N) construction)
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

This produces a true N:1 map in O(N) time without modifying `Jellyrollmodel.jl`.

### 9. Unchanged Modules

| Module | Reason |
|--------|--------|
| `compute_czm_strain_inputs` | Based on `bulk_element` (thermal elements), unchanged |
| `assemble_czm_system` / `build_czm_cache` | Per-CZM-element assembly, naturally compatible |
| `CzmLayout` / `CZMAssemblyWorkspace` | Auto-adapts to new `n_coh` and `ndof` |
| `CzmPostProcess.jl` | Iterates `damage_states` by index, compatible |
| `solve_branch_currents` | `deactivated_elements` logic unchanged (AND logic in `get_active_elements` handles it) |

## Code Impact Summary

| File | Change | Priority |
|------|--------|----------|
| `src/Option.jl` | Add `czm_mesh_refine`, `czm_refine_factor` | P0 |
| `src/CZMMesh.jl` | New file: structs + refinement functions | P0 |
| `src/czm.jl` | Remove migrated code, keep assembly only | P0 |
| `src/JuBat.jl` | Update include order and exports | P0 |
| `src/Materialmatrix.jl` | AND logic in `get_active_elements` | P1 |
| `src/ThermalDistributed.jl` | Use `thermal_node_bottom/top` in BC assembly | P1 |
| `src/SetMesh.jl` | Change `cohesive_elements` type to `Vector{CohesiveElement}` | P1 |

## Verification

1. **Consistency**: uniform damage → default/multiple/auto produce identical thermal and voltage results.
2. **Auto lc**: verify `compute_cohesive_zone_length` returns correct value, `n_i ≥ 1`, all sub-elements ≤ `lc/3`.
3. **Gap conductance**: sum of refined CZM contributions equals original value when damage is uniform.
4. **AND deactivation**: thermal element stays active until all associated CZM elements reach D=1.
5. **Mesh sensitivity**: extend `example/网格敏感性/4_czm_mesh_sensitivity.jl` to test all three modes.
