# Variable Dataflow Unification: Thermal Dead Allocation Cleanup + CZM Variable Optimization

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 统一温度模块和 CZM 模块的变量数据流，使其对齐"Variables.jl 定义 → 模块写入 → PostProcessing.jl 还原"的设计规范。具体包括：(1) 清理温度模块 9 个死内存分配，(2) 统一 CZM 变量键名并补齐预分配，(3) 收敛 CZM 状态传递签名。

**Architecture:** 变量生命周期遵循 `StandardVariables (定义+预分配) → 模块内部写入 → czm_output_to_variables / Solve.jl 输出 → PostProcessing.jl 物理单位还原`。不新增独立状态 struct，仅在 CouplingState.jl 中增加 `CzmLayout` 轻量 struct 收敛签名。

**Tech Stack:** Julia, 现有 JuBat 模块体系，无外部依赖变更。

**Design Specs:**
- `docs/planning-with-files/温度数据流管理/findings.md`
- `docs/planning-with-files/内聚力模块变量优化/findings.md`

**Protected Constraints:**
1. PostProcessing.jl 中已有的物理单位还原公式不可改变（scale.q, scale.T_ref, scale.L, scale.E_n, scale.E_p, scale.r0）
2. Variable_update! 的动态扩展机制不可改变
3. update_czm_damage! 的损伤提交门控（只在收敛时提交）不可改变
4. 现有热源分量还原公式（`thermal2D q_* → thermal2D Q_* [W/m3]`）的映射关系不可改变
5. 电化学变量键不受任何影响

---

## File Structure

### Before (现状)
```
src/Variables.jl           270 行
  - 9 个 thermal2D 死分配 (~lines 97, 100, 123-127, 133-134)
  - 1 个双空格 bug (line 100)
  - CZM 块缺少 6 个场变量预分配 (lines 136-144)
src/CzmPostProcess.jl      117 行
  - czm_output_to_variables 使用 3 个与 StandardVariables 不一致的键名
src/PostProcessing.jl      ~77 行
  - 温度还原路径已统一 (lines 49-64)
src/Solve.jl               ~430 行
  - 热变量直接输出绕过 PostProcessing (lines 389-420)
src/CouplingState.jl       375 行
  - update_czm_damage! 签名有 6 个松散参数
src/SetCase.jl             ~112 行
  - Case 缺少 czm_layout 字段
```

### After (目标)
```
src/Variables.jl           ~260 行 (删除 9 行死分配, 新增 7 行 CZM 预分配, 修复双空格)
src/CzmPostProcess.jl      117 行 (统一 3 个键名)
src/PostProcessing.jl      ~90 行 (从 Solve.jl 迁入热变量输出+还原)
src/Solve.jl               ~400 行 (删除热变量直接输出代码块)
src/CouplingState.jl       ~385 行 (新增 CzmLayout struct)
src/SetCase.jl             ~114 行 (新增 czm_layout 字段)
```

---

## Chunk 1: Thermal Dead Allocation Cleanup (最低风险)

### Task 1: Remove 9 dead variable allocations from Variables.jl

**Files:**
- Modify: `src/Variables.jl`

**Dead allocations audit** (verified: defined but never written anywhere in src/):

| Line | Variable | Why dead |
|------|----------|----------|
| 97 | `thermal2D temperature` | Never written; element temp computed in Solve.jl:422-428 directly |
| 100 | `thermal2D temperature  at nodes history` | Double-space bug + never read/written |
| 123 | `thermal2D element thermal stress` | Never written |
| 124 | `thermal2D element diffusion stress` | Never written |
| 125 | `thermal2D element total stress` | Never written |
| 126 | `thermal2D element diffusion strain` | Never written |
| 127 | `thermal2D element thermal strain` | Never written |
| 133 | `thermal2D displacement x` | Never written (CZM uses `czm displacement x`) |
| 134 | `thermal2D displacement y` | Never written (CZM uses `czm displacement y`) |

- [ ] **Step 1: Delete 9 lines from `src/Variables.jl`**

Remove lines 97, 100, 123-127, 133-134. The remaining `distributed2D` block should look like:

```julia
    if case.opt.thermalmodel == "distributed2D"
        ne = size(case.mesh["thermal2D"].element, 1)
        nT = case.mesh["thermal2D"].nlen
        variables["thermal2D temperature at nodes"] = zeros(Float64, nT, num)
        variables["thermal2D temperature history"] = zeros(Float64, ne, num)
        variables["heat_source_fields"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_sp"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_pcc"] = zeros(Float64, ne, num)
        variables["thermal2D q_ncc"] = zeros(Float64, ne, num)
        variables["thermal2D element current"] = zeros(Float64, ne, num)
        variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p"] = zeros(Float64, ne, num)
        variables["thermal2D element voltages"] = zeros(Float64, ne, num)
        variables["thermal2D element OCV"] = zeros(Float64, ne, num)
        variables["thermal2D active_mask"] = zeros(Float64, ne, num)
        variables["thermal2D n_cutoff_elements"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_element"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_ocv"] = zeros(Float64, 1, num)
        variables["thermal2D margin_to_cutoff"] = zeros(Float64, 1, num)
        variables["total heat source"] = zeros(Float64, 1, num)
    end
```

Note: removed `thermal2D temperature` (line 97), `thermal2D temperature  at nodes history` (line 100 with double space), 5 stress/strain variables (lines 123-127), and 2 displacement variables (lines 133-134).

- [ ] **Step 2: Verify no code references the deleted keys**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
grep -r "thermal2D temperature\"" src/ || echo "OK: no bare thermal2D temperature reference"
grep -r "thermal2D temperature  at nodes history" src/ || echo "OK: no double-space key"
grep -r "thermal2D element thermal stress\|thermal2D element diffusion stress\|thermal2D element total stress\|thermal2D element diffusion strain\|thermal2D element thermal strain" src/ || echo "OK: no stress/strain references"
grep -r "\"thermal2D displacement" src/ || echo "OK: no displacement references"
```

Expected: All `OK` — no code references these deleted keys.

- [ ] **Step 3: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add src/Variables.jl
git commit -m "refactor: remove 9 dead thermal variable allocations from Variables.jl

Remove variables never written anywhere in the codebase:
- thermal2D temperature (element-level, never used)
- thermal2D temperature  at nodes history (double-space bug)
- 5 stress/strain variables (never written)
- thermal2D displacement x/y (never written; CZM uses czm displacement x/y)

Saves ~5-15 MB per simulation run."
```

---

## Chunk 2: CZM Variable Key Unification

### Task 2: Add missing CZM pre-allocations to StandardVariables

**Files:**
- Modify: `src/Variables.jl`

**Current CZM block** (lines 136-144) only pre-allocates 7 keys. `czm_output_to_variables` writes 10 additional keys that are not pre-allocated.

- [ ] **Step 1: Replace CZM block in `src/Variables.jl`**

Replace the existing `if case.opt.czm_enabled == true` block (lines 136-144) with:

```julia
    if case.opt.czm_enabled == true
        n_coh = case.czm_mesh.n_cohesive
        n_czm_nodes = size(case.czm_mesh.nodes, 1)

        # ── 标量统计 ──
        variables["czm D_max"] = zeros(Float64, 1, num)
        variables["czm D_mean"] = zeros(Float64, 1, num)
        variables["czm n_fractured"] = zeros(Float64, 1, num)
        variables["czm δ_max_n"] = zeros(Float64, 1, num)
        variables["czm δ_mean_n"] = zeros(Float64, 1, num)

        # ── 场变量 (per cohesive element / per node) ──
        variables["czm damage"] = zeros(Float64, n_coh, num)
        variables["czm displacement x"] = zeros(Float64, n_czm_nodes, num)
        variables["czm displacement y"] = zeros(Float64, n_czm_nodes, num)
        variables["czm traction normal"] = zeros(Float64, n_coh, num)
        variables["czm traction tangent"] = zeros(Float64, n_coh, num)
        variables["czm separation normal"] = zeros(Float64, n_coh, num)
        variables["czm separation tangent"] = zeros(Float64, n_coh, num)

        # ── 电极级损伤 ──
        variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
        variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
    end
```

- [ ] **Step 2: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK` (may need `case.czm_mesh` to be initialized before StandardVariables is called — verify this is the case)

**Caution:** If `StandardVariables` is called before `czm_mesh` is set on `case`, accessing `case.czm_mesh.n_cohesive` will error on `nothing`. Check the call order in Solve.jl. If czm_mesh is not yet available at StandardVariables time, use fallback dimensions.

- [ ] **Step 3: Commit**

```bash
git add src/Variables.jl
git commit -m "refactor: add CZM field variable pre-allocations to StandardVariables

Pre-allocate czm damage, displacement, traction, separation arrays
that were previously dynamically added by czm_output_to_variables.
Ensures Variable_update! can properly track all CZM variables."
```

### Task 3: Unify CZM key names in czm_output_to_variables

**Files:**
- Modify: `src/CzmPostProcess.jl`

**Key naming inconsistency** (3 keys):

| czm_output_to_variables (写入) | StandardVariables / PostProcessing (期望) |
|-------------------------------|----------------------------------------|
| `czm max damage` | `czm D_max` |
| `czm mean damage` | `czm D_mean` |
| `czm fractured elements` | `czm n_fractured` |

- [ ] **Step 1: Fix 3 key names in `czm_output_to_variables`**

In `src/CzmPostProcess.jl`, replace lines 111-114:

```julia
    # Before:
    new_variables["czm max damage"] = stats.max_D
    new_variables["czm mean damage"] = stats.mean_D
    new_variables["czm fractured elements"] = Float64(stats.n_fractured)

    # After:
    new_variables["czm D_max"] = stats.max_D
    new_variables["czm D_mean"] = stats.mean_D
    new_variables["czm n_fractured"] = Float64(stats.n_fractured)
```

- [ ] **Step 2: Search for any other references to old key names**

```bash
grep -rn "czm max damage\|czm mean damage\|czm fractured elements" src/ example/ tools/
```

Expected: Only the old lines in CzmPostProcess.jl should appear. If other files reference these keys, update them too.

- [ ] **Step 3: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 4: Commit**

```bash
git add src/CzmPostProcess.jl
git commit -m "refactor: unify CZM key names in czm_output_to_variables

czm max damage → czm D_max
czm mean damage → czm D_mean
czm fractured elements → czm n_fractured

Aligns with StandardVariables pre-allocation and PostProcessing restoration."
```

### Task 4: Align create_element_workspace CZM keys

**Files:**
- Modify: `src/Variables.jl`

**Current** (lines 219-223): workspace only has 2 CZM keys (electrode damage). Missing field variables that per-element code may need.

- [ ] **Step 1: Extend CZM block in `create_element_workspace`**

Replace the CZM conditional block (lines 219-223) with:

```julia
    # ── CZM 键（条件）──
    if case.opt.czm_enabled
        n_coh = case.czm_mesh.n_cohesive
        n_czm_nodes = size(case.czm_mesh.nodes, 1)
        ws["negative electrode cohesive zone damage"] = zeros(Float64, Nn, 1)
        ws["positive electrode cohesive zone damage"] = zeros(Float64, Np, 1)
        ws["czm damage"] = zeros(Float64, n_coh, 1)
        ws["czm D_max"] = 0.0
        ws["czm D_mean"] = 0.0
        ws["czm n_fractured"] = 0.0
    end
```

- [ ] **Step 2: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 3: Commit**

```bash
git add src/Variables.jl
git commit -m "refactor: extend CZM keys in create_element_workspace

Add czm damage, D_max, D_mean, n_fractured to per-element workspace."
```

---

## Chunk 3: Thermal Output Path Unification

### Task 5: Move thermal output from Solve.jl to PostProcessing.jl

**Files:**
- Modify: `src/PostProcessing.jl`
- Modify: `src/Solve.jl`

**Problem:** Solve.jl lines 389-420 directly output 13+ thermal variables, bypassing PostProcessing.jl's restoration pattern. PostProcessing.jl lines 49-64 handle heat source and node temperature restoration separately. This creates two output paths with inconsistent style.

**Strategy:** Migrate the thermal variable output from Solve.jl into PostProcessing.jl's `distributed2D` block, using consistent key naming (lowercase internal key → display name with physical units where applicable).

- [ ] **Step 1: Extend PostProcessing.jl `distributed2D` block**

After the existing `thermal2D temperature at nodes [K]` line (line 63), add:

```julia
        # ── 从 Solve.jl 迁入的热变量输出 ──
        result["thermal2D element current [A]"] = variables["thermal2D element current"][:, 1:v] * case.param_dim.cell.I1C
        result["thermal2D element voltages [V]"] = variables["thermal2D element voltages"][:, 1:v] * case.param.scale.phi
        result["thermal2D element OCV [V]"] = variables["thermal2D element OCV"][:, 1:v] * case.param.scale.phi
        result["thermal2D element soc_n"] = variables["thermal2D element soc_n"][:, 1:v]
        result["thermal2D element soc_p"] = variables["thermal2D element soc_p"][:, 1:v]
        result["thermal2D eta_n_e [V]"] = variables["thermal2D eta_n_e"][:, 1:v] * case.param.scale.phi
        result["thermal2D eta_p_e [V]"] = variables["thermal2D eta_p_e"][:, 1:v] * case.param.scale.phi
        result["thermal2D dUdT_n_e [V/K]"] = variables["thermal2D dUdT_n_e"][:, 1:v] * case.param.scale.phi / case.param_dim.scale.T_ref
        result["thermal2D dUdT_p_e [V/K]"] = variables["thermal2D dUdT_p_e"][:, 1:v] * case.param.scale.phi / case.param_dim.scale.T_ref
        result["thermal2D active_mask"] = variables["thermal2D active_mask"][:, 1:v]
        result["thermal2D n_cutoff_elements"] = variables["thermal2D n_cutoff_elements"][1, 1:v]
        result["thermal2D nearest_cutoff_element"] = variables["thermal2D nearest_cutoff_element"][1, 1:v]
        result["thermal2D nearest_cutoff_ocv [V]"] = variables["thermal2D nearest_cutoff_ocv"][1, 1:v] * case.param.scale.phi
        result["thermal2D margin_to_cutoff [V]"] = variables["thermal2D margin_to_cutoff"][1, 1:v] * case.param.scale.phi
        result["heat_source_fields"] = variables["heat_source_fields"][:, 1:v]
        result["total heat source [W]"] = vec(variables["total heat source"][1, 1:v]) * case.param_dim.scale.L^3 / case.param_dim.cell.volume
```

**Caution:** The scale factors above are estimates. Verify each one:
- `element current`: I_typ normalization → × I1C (A)
- `voltages`, `OCV`: phi normalization → × scale.phi (V)
- `soc_n/soc_p`: dimensionless → direct pass
- `eta_n_e/eta_p_e`: phi normalization → × scale.phi (V)
- `dUdT_n_e/dUdT_p_e`: phi/T_ref → × scale.phi / T_ref (V/K)
- `active_mask`, `n_cutoff_elements`: integer-like → direct pass
- `total heat source`: see Solve.jl:401 for existing formula

- [ ] **Step 2: Remove thermal output block from Solve.jl**

Remove the entire `try ... end` block at Solve.jl lines 389-420 that directly outputs thermal variables. Keep any non-thermal code in that region.

- [ ] **Step 3: Verify module loads and test**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 4: Run a thermal example and check output keys**

```bash
julia example/SPMe_Thermal_example.jl
```

Expected: Output should contain all thermal2D keys with physical units. Compare with previous run.

- [ ] **Step 5: Commit**

```bash
git add src/PostProcessing.jl src/Solve.jl
git commit -m "refactor: unify thermal output path through PostProcessing.jl

Move thermal variable output from Solve.jl direct-write to PostProcessing.jl
restoration block. All thermal outputs now follow the same pattern as
electrochemical and CZM variables. Adds physical unit restoration for
current, voltage, OCV, overpotential, and entropy coefficients."
```

---

## Chunk 4: CZM State Convergence

### Task 6: Add CzmLayout struct to CouplingState.jl

**Files:**
- Modify: `src/CouplingState.jl`
- Modify: `src/SetCase.jl`

- [ ] **Step 1: Add CzmLayout struct after CZMAssemblyCache in `src/CouplingState.jl`**

```julia
"""
    CzmLayout

CZM 求解的布局信息和跨时间步状态。对标电化学的 MultiSPMeLayout。
"""
mutable struct CzmLayout
    n_coh::Int                    # cohesive 单元数
    ndof::Int                     # 总位移 DOF 数 (2 * nnode)
    u_prev::Vector{Float64}       # 上一步位移场（跨时间步持有）
end

"""便捷构造器：从 czm_mesh 初始化"""
function CzmLayout(czm_mesh::CohesiveMesh)
    n_coh = czm_mesh.n_cohesive
    ndof = 2 * czm_mesh.nnode
    CzmLayout(n_coh, ndof, zeros(Float64, ndof))
end
```

- [ ] **Step 2: Add czm_layout field to Case struct in `src/SetCase.jl`**

```julia
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
    layout::Union{Nothing, MultiSPMeLayout}
    geometry::Union{Nothing, MeshGeometry}
    czm_mesh::Union{Nothing, CohesiveMesh}
    czm_cache::Union{Nothing, CZMAssemblyCache}
    czm_layout::Union{Nothing, CzmLayout}        # ← 新增
end
```

Update the 5-parameter convenience constructor:

```julia
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, nothing, nothing)
end
```

- [ ] **Step 3: Initialize CzmLayout when czm_mesh is created**

Find where `case.czm_mesh` is first assigned (likely in `SetCase` or a setup function) and add:

```julia
    case.czm_layout = CzmLayout(case.czm_mesh)
```

- [ ] **Step 4: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl src/SetCase.jl
git commit -m "refactor: add CzmLayout struct for CZM state convergence

CzmLayout holds n_coh, ndof, and u_prev (previous step displacement).
Stored in Case.czm_layout alongside existing czm_mesh and czm_cache.
Prepares for converging update_czm_damage! signature."
```

### Task 7: Converge update_czm_damage! signature

**Files:**
- Modify: `src/CouplingState.jl` (update_czm_damage!)
- Modify: `src/Solve.jl` (caller)

**Current signature:**
```julia
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
```

**Target signature:**
```julia
function update_czm_damage!(case, variables, T_nodes_carry)
```

- [ ] **Step 1: Update update_czm_damage! in `src/CouplingState.jl`**

Replace the function signature and internal parameter extraction:

```julia
function update_czm_damage!(case, variables, T_nodes_carry)
    czm_mesh = case.czm_mesh
    czm_params = case.param.cohesive
    czm_layout = case.czm_layout
    u_czm_prev = czm_layout.u_prev
    param = case.param

    # ... rest of function body unchanged ...

    # After solve, update u_prev if converged
    if result.converged
        czm_mesh.damage_states = updated_czm_mesh.damage_states
        czm_layout.u_prev = result.displacement
    end

    return result.displacement, result.converged
end
```

Keep the old 6-parameter signature as a thin wrapper for backward compatibility:

```julia
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
    case.czm_layout = CzmLayout(czm_mesh)
    case.czm_layout.u_prev = u_czm_prev !== nothing ? u_czm_prev : zeros(Float64, 2 * czm_mesh.nnode)
    u_czm, converged = update_czm_damage!(case, variables, T_nodes_carry)
    return u_czm, converged
end
```

- [ ] **Step 2: Update caller in `src/Solve.jl`**

Find the call site (approximately Solve.jl:276) and simplify:

```julia
    # Before:
    u_czm, czm_converged = update_czm_damage!(case.czm_mesh, case.param.cohesive, case, variables, T_nodes_carry, u_czm_prev)

    # After:
    u_czm, czm_converged = update_czm_damage!(case, variables, T_nodes_carry)
```

- [ ] **Step 3: Check for other callers**

```bash
grep -rn "update_czm_damage!" src/ tools/ example/
```

If other callers use the 6-parameter signature, the backward-compatible wrapper handles them.

- [ ] **Step 4: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl src/Solve.jl
git commit -m "refactor: converge update_czm_damage! to (case, variables, T_nodes_carry)

Extract czm_mesh, czm_params, u_prev from case struct.
Old 6-parameter signature preserved as thin wrapper for backward compatibility."
```

---

## Chunk 5: Final Verification

### Task 8: Verify all changes together

- [ ] **Step 1: Run module load check**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

- [ ] **Step 2: Run CZM baseline probe**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: All output values match pre-change baseline.

- [ ] **Step 3: Run thermal example**

```bash
julia example/SPMe_Thermal_example.jl
```

Expected: No errors, output contains all thermal2D keys.

- [ ] **Step 4: Run full coupling example**

```bash
julia example/testexample.jl
```

Expected: No errors.

- [ ] **Step 5: Verify exports**

```bash
julia -e '
include("src/JuBat.jl")
using .JuBat
for sym in [:update_czm_damage!, :czm_output_to_variables, :get_damage_statistics]
    @assert getfield(JuBat, sym) !== nothing
end
println("All exports OK")
'
```

- [ ] **Step 6: Update progress files**

Update `docs/planning-with-files/温度数据流管理/progress.md` and `docs/planning-with-files/内聚力模块变量优化/progress.md`.

- [ ] **Step 7: Final commit**

```bash
git add docs/planning-with-files/
git commit -m "docs: update progress for variable dataflow unification"
```

---

## Rollback Instructions

```bash
# View recent commits
git log --oneline -10

# Revert specific task
git revert <commit-hash>

# Full rollback
git reset --hard HEAD~8   # 8 tasks total
```

Allowed numerical tolerance: relative error < 1e-8 (determined by solver tol).
