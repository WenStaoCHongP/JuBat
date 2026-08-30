# CZM Module Refactoring: Split & Simplify Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 CZM 模块从当前的单体文件拆分为职责清晰的多文件结构，消除 ~170 行重复代码，保持所有公开接口和数值行为不变。

**Architecture:** 保留 `CzmSolve.jl` 作为求解器主入口，拆出 `CzmPostProcess.jl`（统计/后处理），并把 3 个耦合 helper 收入现有 `CouplingState.jl`（Julia module 级 include 保证函数调用在运行时解析，无需调整 include 顺序）。然后在 `CzmSolve.jl` 内部提取 3 个公共 helper 消除重复。

**Tech Stack:** Julia, 现有 JuBat 模块体系（`include` + `export` 模式），无外部依赖变更。

**Design Spec:** `docs/planning-with-files/04_内聚力模块拆分/task_plan.md`, `docs/planning-with-files/04_内聚力模块拆分/findings.md`

**Protected Constraints (不可破坏):**
1. 线搜索：`for _ in 1:8; α *= 0.5; end`
2. 失败回滚：`if !converged; u = u_start; damage_states = damage_start; end`
3. 损伤提交门控：`if result.converged; czm_mesh.damage_states = updated_czm_mesh.damage_states; end`
4. 载荷子步自适应：`step_size *= 0.5` / `step_size *= 1.25`
5. BC 处理：`newton_raphson_czm` 使用惩罚式残差 `R[dof] = val - u[dof]`；`solve_czm_basic_step` 使用零化式 `R[dof] = 0.0`；二者不可混用

---

## File Structure

### Before (现状)
```
src/CzmSolve.jl       1001 行, 19 函数 (调度+求解+后处理+耦合入口混杂)
src/czm.jl             659 行, 11 函数 (组装+网格+本构)
src/CouplingState.jl   187 行, 类型定义 + 缓存/工作区
src/JuBat.jl            79 行, include 顺序
```

### After (目标)
```
src/CzmSolve.jl       ~560 行, 求解器 + helpers (缩减 ~440 行)
src/CzmPostProcess.jl  ~120 行, 统计/后处理/损伤管理 (新文件)
src/CouplingState.jl   ~280 行, 类型 + 缓存 + 3 个耦合 helper (新增 ~90 行)
src/czm.jl             659 行, 不变
src/JuBat.jl            80 行, 新增 1 个 include (CzmPostProcess.jl)
```

---

## Chunk 1: Scaffolding & Baseline

### Task 1: Create Feature Branch & Record Baseline

**Files:**
- Reference: `tools/czm_baseline_probe.jl`

- [ ] **Step 1: Create feature branch**

```bash
git checkout -b czm-refactor
```

Expected: Switched to branch `czm-refactor`

- [ ] **Step 2: Run baseline probe**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
julia tools/czm_baseline_probe.jl
```

Expected: 输出 `BASELINE_START` ... `BASELINE_END` 段，记录三种方法的 `converged`, `iterations`, `residual_norm`, `D_max`, `D_mean`, `n_fractured`。将这些值复制到 `docs/planning-with-files/04_内聚力模块拆分/findings.md` 的 `## Baseline Snapshot` 段。

- [ ] **Step 3: Commit baseline**

```bash
git add tools/czm_baseline_probe.jl docs/planning-with-files/04_内聚力模块拆分/findings.md
git commit -m "chore: add CZM baseline probe for refactoring regression check"
```

---

## Chunk 2: File Split — CzmPostProcess.jl

### Task 2: Create CzmPostProcess.jl with 5 functions

**Files:**
- Create: `src/CzmPostProcess.jl`
- Modify: `src/CzmSolve.jl` (删除 5 个函数: 行 676-816)
- Modify: `src/JuBat.jl` (新增 1 行 include)

- [ ] **Step 1: Verify module loads before change**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 2: Create `src/CzmPostProcess.jl`**

将 `src/CzmSolve.jl` 行 676-816 的 5 个函数原封不动搬到新文件:
`get_damage_statistics`, `check_fracture_criterion`, `reset_damage_states`, `accumulate_cycle_damage`, `czm_output_to_variables`

文件头部注释:
```julia
# ========================================================================
# CZM Post-processing, Statistics & Damage Management
# ========================================================================
# 从 CzmSolve.jl 拆分出来的后处理、统计和损伤管理函数。
# 这些函数无状态、纯计算，独立于求解器主循环。
# 注: Statistics 已在 JuBat.jl module 级别 using，此处无需重复。
# ========================================================================
```

函数实现与 `src/CzmSolve.jl:684-816` 完全相同，不做任何修改。

- [ ] **Step 3: Delete the 5 functions from `src/CzmSolve.jl`**

删除行 676-816（从 `# 8. Damage statistics` 注释到 `czm_output_to_variables` 的 `end`）。

- [ ] **Step 4: Add include to `src/JuBat.jl`**

在 `include("CzmSolve.jl")` 之后添加:
```julia
include("CzmPostProcess.jl")  # CZM 后处理/统计/损伤管理
```

- [ ] **Step 5: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 6: Run baseline probe and compare**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: 输出与 Task 1 Step 2 基线值完全一致。

- [ ] **Step 7: Commit**

```bash
git add src/CzmPostProcess.jl src/CzmSolve.jl src/JuBat.jl
git commit -m "refactor: extract CzmPostProcess.jl from CzmSolve.jl

Move 5 functions (get_damage_statistics, check_fracture_criterion,
reset_damage_states, accumulate_cycle_damage, czm_output_to_variables)
to dedicated post-processing file. No behavior change."
```

---

## Chunk 3: Move Coupling Helpers into CouplingState.jl

### Task 3: Move 3 coupling functions to existing CouplingState.jl

**Files:**
- Modify: `src/CouplingState.jl` (在文件末尾追加 3 个函数)
- Modify: `src/CzmSolve.jl` (删除 3 个函数)

**Include 顺序说明:** `CouplingState.jl` 目前在 `czm.jl` 和 `CzmSolve.jl` 之前被 include。`update_czm_damage!` 调用 `solve_czm_step`（定义在 `CzmSolve.jl`）和 `ensure_czm_cache`（定义在 `czm.jl`）。这在 Julia 中合法：module 级 `include` 是文本插入，函数调用在运行时按符号名解析，不需要定义在调用点之前。只有类型定义（struct）才必须在首次使用前定义。3 个耦合函数的签名不使用任何 CZM 专有类型（参数均为 duck-typed），所以无需调整 include 顺序。

- [ ] **Step 1: Append coupling functions to `src/CouplingState.jl`**

在文件末尾（`CZMAssemblyCache` 的 `end` 之后）追加以下代码。函数实现与 `src/CzmSolve.jl:835-1001` 完全相同:

```julia

# ========================================================================
# CZM Coupling Helpers
# ========================================================================
# 参数计算、应变输入、损伤更新入口。
# 从 CzmSolve.jl 迁移至此，保持 CouplingState.jl 作为"状态+耦合"的统一归属点。
# 函数调用 solve_czm_step / ensure_czm_cache 等在运行时解析，无需调整 include 顺序。
# ========================================================================

"""
    compute_czm_effective_params(case)
"""
function compute_czm_effective_params(case::Case)
    param = case.param
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)
    ν_eff = (param.NE.nu * param.NE.thickness + param.PE.nu * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)
    α_eff = (param.NE.alphaT * param.NE.thickness + param.PE.alphaT * param.PE.thickness) /
        (param.NE.thickness + param.PE.thickness)
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0
    return E_eff, ν_eff, α_eff, β_n, β_p
end

"""
    compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)
"""
function compute_czm_strain_inputs(case::Case, variables::Dict, czm_mesh, T_nodes_carry)
    ne = size(czm_mesh.bulk_element, 1)
    param = case.param
    soc_ref_n = param.NE.cs0
    soc_ref_p = param.PE.cs0
    dT_elem = zeros(Float64, ne)
    Δsoc_n_elem = zeros(Float64, ne)
    Δsoc_p_elem = zeros(Float64, ne)

    if length(T_nodes_carry) >= czm_mesh.nnode
        for e in 1:ne
            nodes = czm_mesh.bulk_element[e, :]
            T_elem_nd = 0.0
            valid_nodes = 0
            for n in nodes
                if n <= length(T_nodes_carry)
                    T_elem_nd += T_nodes_carry[n]
                    valid_nodes += 1
                end
            end
            if valid_nodes > 0
                T_elem_nd /= valid_nodes
                dT_elem[e] = T_elem_nd - param.cell.T0
            end
        end
    end

    soc_n_elem = variables["thermal2D element soc_n"]
    soc_p_elem = variables["thermal2D element soc_p"]
    if isa(soc_n_elem, AbstractMatrix)
        soc_n_elem = soc_n_elem[:, end]
        soc_p_elem = soc_p_elem[:, end]
    end
    for e in 1:min(ne, length(soc_n_elem))
        Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
        Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
    end
    return dT_elem, Δsoc_n_elem, Δsoc_p_elem
end

"""
    update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)

更新 CZM 网格的损伤状态。被 Solve.jl:276 调用。
"""
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
    param = case.param
    czm_params.czm_model = case.opt.czm_model

    E_eff, ν_eff, α_eff, β_n, β_p = compute_czm_effective_params(case)
    cache = ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)
    dT_elem, Δsoc_n_elem, Δsoc_p_elem = compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)

    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)

    has_nan_T = any(isnan, T_nodes_carry)
    has_nan_soc_n = any(isnan, variables["thermal2D element soc_n"])
    has_nan_soc_p = any(isnan, variables["thermal2D element soc_p"])
    if has_nan_T || has_nan_soc_n || has_nan_soc_p || any(isnan, dT_elem) || any(isnan, Δsoc_n_elem) || any(isnan, Δsoc_p_elem)
        @warn "CZM inputs contain NaN" has_nan_T=has_nan_T has_nan_soc_n=has_nan_soc_n has_nan_soc_p=has_nan_soc_p
    end

    if u_czm_prev === nothing || length(u_czm_prev) != ndof
        u_czm_prev = zeros(Float64, ndof)
    elseif any(isnan, u_czm_prev)
        @warn "CZM u_czm_prev contains NaN, resetting to zeros"
        u_czm_prev = zeros(Float64, ndof)
    end

    result, updated_czm_mesh = solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, czm_params, param, u_czm_prev;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=case.opt.czm_max_iter, tol=case.opt.czm_tol,
        n_load_steps=case.opt.czm_load_steps, arc_length_alpha=case.opt.czm_arc_length_alpha,
        iter_method=case.opt.czm_iter_method, cache=cache)

    has_nan_disp = any(isnan, result.displacement)
    has_nan_damage = any(ds -> isnan(ds.D), updated_czm_mesh.damage_states)
    if has_nan_disp || has_nan_damage || !result.converged
        @warn "CZM solve issue" converged=result.converged iterations=result.iterations residual=round(result.residual_norm; digits=4)
    end

    if result.converged
        czm_mesh.damage_states = updated_czm_mesh.damage_states
    end

    return result.displacement, result.converged
end
```

- [ ] **Step 2: Delete the 3 functions from `src/CzmSolve.jl`**

删除 `compute_czm_effective_params`、`compute_czm_strain_inputs`、`update_czm_damage!` 及其上方的 `# ==...CZM 损伤更新` 注释块。即从约行 818 到文件末尾全部删除。

- [ ] **Step 3: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`。如果报 `UndefVarError: solve_czm_step not defined`，说明 include 顺序有问题——但理论上不应该，因为 Julia module 编译期不要求函数定义顺序。

- [ ] **Step 4: Run baseline probe and compare**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: 输出与基线值完全一致。

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl src/CzmSolve.jl
git commit -m "refactor: move CZM coupling helpers into CouplingState.jl

Move compute_czm_effective_params, compute_czm_strain_inputs,
update_czm_damage! to existing CouplingState.jl.
CzmSolve.jl now only owns solver logic (~670 lines). No behavior change."
```

---

## Chunk 4: Helper Extraction — clone_czm_mesh_with_new_states

### Task 4: Extract shared mesh cloning helper

**Files:**
- Modify: `src/CzmSolve.jl` (替换 `newton_raphson_czm` 中的内联克隆)
- Modify: `src/CzmPostProcess.jl` (替换 `reset_damage_states` 和 `accumulate_cycle_damage` 中的内联克隆)

- [ ] **Step 1: Add helper function to `src/CzmSolve.jl`**

在 `clone_czm_mesh_with_damage` 函数之后添加:

```julia
"""
    clone_czm_mesh_with_new_states(czm_mesh, new_damage_states)

通用 CZM 网格克隆 helper：用新的损伤状态创建 mesh 浅拷贝。
"""
function clone_czm_mesh_with_new_states(czm_mesh::CohesiveMesh, new_damage_states::Vector{<:AbstractDamageState})
    new_czm_mesh = CohesiveMesh()
    new_czm_mesh.bulk_mesh = czm_mesh.bulk_mesh
    new_czm_mesh.node = czm_mesh.node
    new_czm_mesh.nnode = czm_mesh.nnode
    new_czm_mesh.bulk_element = czm_mesh.bulk_element
    new_czm_mesh.cohesive_elements = czm_mesh.cohesive_elements
    new_czm_mesh.n_cohesive = czm_mesh.n_cohesive
    new_czm_mesh.n_layers = czm_mesh.n_layers
    new_czm_mesh.node_map = czm_mesh.node_map
    new_czm_mesh.interface_nodes = czm_mesh.interface_nodes
    new_czm_mesh.damage_states = new_damage_states
    return new_czm_mesh
end
```

- [ ] **Step 2: Refactor `newton_raphson_czm` in `src/CzmSolve.jl`**

将函数末尾的内联克隆代码（约 12 行 `new_czm_mesh = CohesiveMesh()` ... `new_czm_mesh.damage_states = damage_states`）替换为:
```julia
    new_czm_mesh = clone_czm_mesh_with_new_states(czm_mesh, damage_states)
```

- [ ] **Step 3: Refactor `reset_damage_states` in `src/CzmPostProcess.jl`**

将 `new_czm_mesh = CohesiveMesh()` ... `new_czm_mesh.damage_states = new_damage_states`（约 11 行）替换为:
```julia
    new_czm_mesh = clone_czm_mesh_with_new_states(czm_mesh, new_damage_states)
```

- [ ] **Step 4: Refactor `accumulate_cycle_damage` in `src/CzmPostProcess.jl`**

同样替换为:
```julia
    new_czm_mesh = clone_czm_mesh_with_new_states(czm_mesh, new_damage_states)
```

- [ ] **Step 5: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 6: Run baseline probe and compare**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: 输出与基线值完全一致。

- [ ] **Step 7: Commit**

```bash
git add src/CzmSolve.jl src/CzmPostProcess.jl
git commit -m "refactor: extract clone_czm_mesh_with_new_states helper

Eliminates 3 copies of mesh cloning boilerplate (~33 lines saved).
No behavior change."
```

---

## Chunk 5: Helper Extraction — extract_bc_dofs

### Task 5: Extract shared BC extraction helper

**Files:**
- Modify: `src/CzmSolve.jl` (替换 3 处 BC 提取逻辑)

- [ ] **Step 1: Add helper function to `src/CzmSolve.jl`**

在 `clone_czm_mesh_with_new_states` 之后添加:

```julia
"""
    extract_bc_dofs(czm_mesh, param; cache=nothing)

从 czm_mesh 提取 Dirichlet BC 的自由度列表和对应值。
优先使用缓存中的 bc_dofs/bc_vals，否则从 identify_bc_nodes_czm 重新计算。
"""
function extract_bc_dofs(czm_mesh::CohesiveMesh, param; cache::Union{Nothing, CZMAssemblyCache}=nothing)
    if cache !== nothing
        return cache.bc_dofs, cache.bc_vals
    end
    bc_nodes, _, _ = identify_bc_nodes_czm(czm_mesh, param)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        end
    end
    return bc_dofs, bc_vals
end
```

- [ ] **Step 2: Refactor `solve_czm_basic_step`**

将 BC 提取代码块（约 20 行 `if cache !== nothing ... else ... end`）替换为:
```julia
        bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)
```

- [ ] **Step 3: Refactor `solve_czm_arc_length_step`**

同上替换。

- [ ] **Step 4: Refactor `newton_raphson_czm`**

同上替换。

- [ ] **Step 5: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 6: Run baseline probe and compare**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: 输出与基线值完全一致。

- [ ] **Step 7: Commit**

```bash
git add src/CzmSolve.jl
git commit -m "refactor: extract extract_bc_dofs helper

Eliminates 3 copies of BC extraction logic (~60 lines saved).
No behavior change."
```

---

## Chunk 6: Helper Extraction — backtrack_line_search!

### Task 6: Extract line search helper (basic solver only)

**Files:**
- Modify: `src/CzmSolve.jl`

**关键约束:** `solve_czm_basic_step` 使用零化式残差 `R[dof] = 0.0`（通过 `apply_czm_dirichlet!`），而 `newton_raphson_czm` 使用惩罚式残差 `R[dof] = val - u[dof]`。二者语义不同，不可混用。因此 `backtrack_line_search!` **只服务 `solve_czm_basic_step`**，`newton_raphson_czm` 和 `solve_czm_arc_length_step` 保持原有线搜索不变。

- [ ] **Step 1: Add helper function to `src/CzmSolve.jl`**

在 `extract_bc_dofs` 之后添加:

```julia
"""
    backtrack_line_search!(u, Δu, czm_mesh, E_eff, ν_eff, cohesive_params,
                           damage_states, F_ext, F_thermo_chem, R_norm_current,
                           bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws;
                           max_halvings=8)

回溯线搜索（零化式 BC 残差）。返回 (u_new, R_new_norm, accepted, α_used)。
仅用于 solve_czm_basic_step。

不变量：accepted 时返回的 u_new 已含 BC 赋值，外部无需再执行 u = u + α*Δu。
        未 accepted 时返回原始 u（未修改），外部应 break。
"""
function backtrack_line_search!(u::Vector{Float64}, Δu::Vector{Float64},
                                czm_mesh::CohesiveMesh, E_eff::Float64, ν_eff::Float64,
                                cohesive_params::Cohesive, damage_states,
                                F_ext::Vector{Float64}, F_thermo_chem::Vector{Float64},
                                R_norm_current::Float64,
                                bc_dofs::Vector{Int64}, bc_vals::Vector{Float64},
                                K_bulk_cached, geom_cache, ws;
                                max_halvings::Int=8)
    α = 1.0
    for _ in 1:max_halvings
        u_trial = u + α * Δu
        apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

        _, f_int_trial, _, _ = assemble_coupled_system(
            czm_mesh, u_trial, E_eff, ν_eff, cohesive_params;
            damage_states=damage_states, K_bulk_cached=K_bulk_cached,
            geom_cache=geom_cache, ws=ws)

        R_trial = F_ext + F_thermo_chem - f_int_trial
        apply_czm_dirichlet!(R_trial, bc_dofs, zeros(length(bc_dofs)))

        R_trial_norm = norm(R_trial)
        if !isnan(R_trial_norm) && R_trial_norm < R_norm_current
            return u_trial, R_trial_norm, true, α
        end

        α *= 0.5
    end
    return u, R_norm_current, false, 0.0
end
```

- [ ] **Step 2: Refactor `solve_czm_basic_step`**

将线搜索代码块（约 22 行 `α = 1.0; ls_accepted = false; for _ in 1:8 ... end`）及后续位移更新替换为:

```julia
            u, R_norm_new, ls_accepted, α_used = backtrack_line_search!(
                u, Δu, czm_mesh, E_eff, ν_eff, cohesive_params,
                damage_states, F_ext, F_thermo_chem, R_norm,
                bc_dofs, bc_vals, K_bulk_cached, geom_cache, ws_basic)

            if !ls_accepted
                break
            end
```

注意: 删除原来的 `u = u + α * Δu` 和 `apply_czm_dirichlet!(u, bc_dofs, bc_vals)` 这两行，因为 helper 已在 accepted 时返回包含 BC 赋值的 `u_trial`。

- [ ] **Step 3: Verify module loads**

```bash
julia -e 'include("src/JuBat.jl"); using .JuBat; println("OK")'
```

Expected: `OK`

- [ ] **Step 4: Run baseline probe — focus on method=basic**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: `method=basic` 行的 `converged`, `iterations`, `residual_norm`, `D_max`, `D_mean` 与基线值完全一致。`method=load_substep` 和 `method=arc_length` 不受影响也应一致。

- [ ] **Step 5: Commit**

```bash
git add src/CzmSolve.jl
git commit -m "refactor: extract backtrack_line_search! for solve_czm_basic_step

Only applied to basic solver (zero-based BC residual).
newton_raphson_czm and arc_length keep their own line search
(penalty-based and arc-constrained respectively). No behavior change."
```

---

## Chunk 7: Final Verification & Cleanup

### Task 7: Final verification and line count

**Files:**
- Read: all modified files
- Modify: `docs/planning-with-files/04_内聚力模块拆分/progress.md`

- [ ] **Step 1: Run full baseline probe**

```bash
julia tools/czm_baseline_probe.jl
```

Expected: 所有输出值与基线一致。

- [ ] **Step 2: Count lines**

```bash
wc -l src/CzmSolve.jl src/CzmPostProcess.jl src/CouplingState.jl src/czm.jl
```

Expected approximate targets:
```
src/CzmSolve.jl       ~560 行 (原 1001)
src/CzmPostProcess.jl  ~120 行 (新)
src/CouplingState.jl   ~280 行 (原 187 + ~90 新)
src/czm.jl             659 行 (不变)
```

- [ ] **Step 3: Verify exports still work**

```bash
julia -e '
include("src/JuBat.jl")
using .JuBat
for sym in [:solve_czm_step, :newton_raphson_czm, :get_damage_statistics,
            :check_fracture_criterion, :reset_damage_states,
            :accumulate_cycle_damage, :czm_output_to_variables,
            :compute_czm_effective_params, :compute_czm_strain_inputs,
            :update_czm_damage!]
    @assert getfield(JuBat, sym) !== nothing
end
println("All exports OK")
'
```

Expected: `All exports OK`

- [ ] **Step 4: Update progress.md with final stats**

- [ ] **Step 5: Final commit**

```bash
git add docs/planning-with-files/04_内聚力模块拆分/progress.md
git commit -m "docs: update CZM refactoring progress — all phases complete"
```

---

## Rollback Instructions

```bash
# 查看提交历史
git log --oneline -10

# 回退到特定 Task 之前
git revert <commit-hash>

# 或者完全重置到 refactoring 之前
git checkout main
git branch -D czm-refactor
```

允许的数值偏差: 残差相对误差 < 1e-8 (由求解器 tol 决定)。
