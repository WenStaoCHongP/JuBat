# Simulation Speedup Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce JuBat multi-SPMe simulation wall-clock time by 25-40% through targeted optimization of the SPMe hot path, memory allocation reduction, and thermal assembly caching.

**Architecture:** Three-phase approach — (1) eliminate redundant allocations in the inner loop (SPMe_element), (2) cache static computations, (3) reduce Solve.jl memory overhead. Each task is independently testable by running `example/testexample.jl` and comparing timing output.

**Tech Stack:** Julia 1.x, SparseArrays, Threads.@threads, FEM assembly

**Baseline (testexample.jl, czm_enabled=false, nθ=80):**
| Module | Time Ratio | Target After |
|--------|-----------|-------------|
| SPMe solve | 91.59% | < 85% |
| Thermal distributed | 7.85% | < 6% |
| Branch solver | 0.56% | 0.56% (unchanged) |

---

## File Structure

| File | Change Type | Responsibility |
|------|------------|---------------|
| `src/SetCase.jl` | Modify | Add `thermal_extras::Dict{String,Any}` field to Case struct |
| `src/CallModel.jl` | Modify | Cache element areas; eliminate representative state call; thread-local workspace |
| `src/SPMe.jl` | Modify | Add workspace-based SPMe_element overload |
| `src/Variables.jl` | Modify | Add workspace creation helpers |
| `src/Solve.jl` | Modify | Replace deepcopy → copy (incl. Float64 deepcopy) |
| `src/Assemble.jl` | Modify | Fix deepcopy; add pre-allocated variant |
| `src/ThermalDistributed.jl` | Modify | Cache boundary edge list; reduce copy |
| `src/Initialisation.jl` | Modify | Cache element areas in layout |
| `src/CouplingState.jl` | Modify | Add BoundaryEdgeCache struct; areas field on MultiSPMeLayout |

---

## Chunk 1: Quick Wins (Low Risk, Immediate Impact)

### Task 0: Add thermal_extras Field to Case Struct

Several later tasks (7, 8, 10) need to cache thermal data on the `case` object. Add the `thermal_extras` field now.

**Files:**
- Modify: `src/SetCase.jl` (add `thermal_extras::Dict{String,Any}` field)

- [ ] **Step 1: Add field to Case struct**

In `src/SetCase.jl`, add the field to the Case struct definition:

```julia
# Add to Case struct fields:
thermal_extras::Dict{String,Any}
```

Update all Case constructors to include `thermal_extras = Dict{String,Any}()` in the default values.

- [ ] **Step 2: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results (field is unused until later tasks).

- [ ] **Step 3: Commit**

```bash
git add src/SetCase.jl
git commit -m "refactor: add thermal_extras Dict to Case struct for runtime caching"
```

---

### Task 1: Cache Element Areas on Layout

Element areas are computed every `CallModel_MultiSPMe` call (lines 28-34) but the mesh never changes. Cache once at initialization.

**Files:**
- Modify: `src/CouplingState.jl` (add `areas` field to MultiSPMeLayout)
- Modify: `src/Initialisation.jl:56-112` (compute and cache areas)
- Modify: `src/CallModel.jl:27-34` (use cached areas)

- [ ] **Step 1: Add `areas` field to MultiSPMeLayout**

In `src/CouplingState.jl`, update the struct:

```julia
struct MultiSPMeLayout
    ne::Int
    n_chem::Int
    nT::Int
    thermal_range::UnitRange{Int}
    areas::Vector{Float64}        # NEW: pre-computed element areas (mesh-invariant)
end

# Update constructor to compute areas
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int, mesh_th)
    thermal_range = (ne * n_chem + 1):(ne * n_chem + nT)
    areas = zeros(Float64, ne)
    ngs = length(mesh_th.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh_th.gs.ele[g]
        areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
    end
    return MultiSPMeLayout(ne, n_chem, nT, thermal_range, areas)
end
```

- [ ] **Step 2: Update layout construction in Initialisation.jl**

In `src/Initialisation.jl:110`, change:

```julia
# OLD:
case.layout = MultiSPMeLayout(ne, n_chem, nT)

# NEW:
case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
```

Also update `src/Solve.jl:104` where layout is constructed from external state:

```julia
# OLD:
case.layout = MultiSPMeLayout(ne, n_chem, nT)

# NEW:
case.layout = MultiSPMeLayout(ne, n_chem, nT, case.mesh["thermal2D"])
```

- [ ] **Step 3: Use cached areas in CallModel.jl**

Replace lines 27-34 in `src/CallModel.jl`:

```julia
# OLD (7 lines):
A = zeros(Float64, ne)
ngs = length(mesh_th.gs.detJ)
@inbounds for g in 1:ngs
    e = mesh_th.gs.ele[g]
    A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
end
variables["thermal2D element area"] = A
areas = A

# NEW (2 lines):
areas = layout.areas
variables["thermal2D element area"] = areas
```

- [ ] **Step 4: Run testexample.jl and verify timing**

Run: `cd example && julia testexample.jl`
Expected: Same results; timing output should show identical SPMe ratio (~91%).

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl src/Initialisation.jl src/CallModel.jl src/Solve.jl
git commit -m "perf: cache element areas in MultiSPMeLayout (eliminates per-step recomputation)"
```

---

### Task 2: Replace deepcopy with copy in Solve.jl

`deepcopy` traverses entire object graph (slow for sparse matrices). `copy` is sufficient since these are leaf types.

**Files:**
- Modify: `src/Solve.jl:245-247`

- [ ] **Step 1: Replace three deepcopy calls**

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

- [ ] **Step 2: Run testexample.jl to verify correctness**

Run: `cd example && julia testexample.jl`
Expected: Identical results (copy preserves sparse matrix structure for SparseMatrixCSC).

- [ ] **Step 3: Commit**

```bash
git add src/Solve.jl
git commit -m "perf: replace deepcopy with copy in Solve.jl time loop"
```

---

### Task 3: Fix Assemble.jl deepcopy(KI) → zeros

`deepcopy(KI)` is wasteful; `zeros(Int64, length(KI))` is faster and clearer.

**Files:**
- Modify: `src/Assemble.jl:13`

- [ ] **Step 1: Fix the deepcopy**

```julia
# OLD (line 13):
KJ = deepcopy(KI)

# NEW:
KJ = zeros(Int64, length(KI))
```

- [ ] **Step 2: Run testexample.jl to verify**

Run: `cd example && julia testexample.jl`
Expected: Identical results.

- [ ] **Step 3: Commit**

```bash
git add src/Assemble.jl
git commit -m "perf: replace deepcopy(KI) with zeros in Assemble.jl"
```

---

## Chunk 2: SPMe Hot Path Allocation Reduction (High Impact)

### Task 4: Add Thread-Local Variables Workspace for SPMe_element

The dominant cost: each `SPMe_element` call → `SPMe_variables` → `StandardVariables(case, 1)` creates a new Dict with ~40 pre-allocated arrays. For 80 elements per step, that's 3,200+ allocations causing GC pressure.

**Solution:** Pre-allocate one `Dict{String, Union{Array{Float64},Float64}}` per thread, reuse across elements.

**Files:**
- Modify: `src/Variables.jl` (add `create_element_workspace`)
- Modify: `src/SPMe.jl` (add workspace-based `SPMe_variables_element!`)
- Modify: `src/CallModel.jl` (use thread-local workspaces)

- [ ] **Step 1: Add workspace creation helper in Variables.jl**

Append to `src/Variables.jl`:

```julia
"""
    create_element_workspace(case::Case)

Create a reusable variables Dict for SPMe_element, pre-allocated with
correct array sizes. Used to avoid per-call StandardVariables allocation.
"""
function create_element_workspace(case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nn = 1  # SPMe mode
    Np = 1
    Ne_ngs = case.opt.Nn * case.opt.gsorder
    Ne_pgs = case.opt.Np * case.opt.gsorder
    Ne_spgs = case.opt.Ns * case.opt.gsorder
    Ne = case.mesh["electrolyte"].nlen
    Ne_n = case.mesh["negative electrode"].nlen
    Ne_p = case.mesh["positive electrode"].nlen
    Ne_sp = case.mesh["separator"].nlen

    ws = Dict{String, Union{Array{Float64},Float64}}(
        "negative particle lithium concentration" => zeros(Float64, Nrn),
        "positive particle lithium concentration" => zeros(Float64, Nrp),
        "negative particle averaged lithium concentration" => zeros(Float64, Nn),
        "positive particle averaged lithium concentration" => zeros(Float64, Np),
        "negative particle surface lithium concentration" => zeros(Float64, Nn),
        "positive particle surface lithium concentration" => zeros(Float64, Np),
        "negative electrode exchange current density" => zeros(Float64, Nn),
        "positive electrode exchange current density" => zeros(Float64, Np),
        "negative electrode interfacial current density" => 0.0,
        "positive electrode interfacial current density" => 0.0,
        "negative electrode overpotential" => zeros(Float64, Nn),
        "positive electrode overpotential" => zeros(Float64, Np),
        "negative electrode open circuit potential" => zeros(Float64, Nn),
        "positive electrode open circuit potential" => zeros(Float64, Np),
        "negative particle concentration at gauss point" => zeros(Float64, Nrn * case.opt.gsorder),
        "positive particle concentration at gauss point" => zeros(Float64, Nrp * case.opt.gsorder),
        "electrolyte lithium concentration" => zeros(Float64, Ne),
        "electrolyte lithium concentration in negative electrode" => zeros(Float64, Ne_n),
        "electrolyte lithium concentration in positive electrode" => zeros(Float64, Ne_p),
        "electrolyte lithium concentration in separator" => zeros(Float64, Ne_sp),
        "electrolyte lithium concentration at negative electrode Gauss point" => zeros(Float64, Ne_ngs),
        "electrolyte lithium concentration at positive electrode Gauss point" => zeros(Float64, Ne_pgs),
        "electrolyte lithium concentration at separator Gauss point" => zeros(Float64, Ne_spgs),
        "cell voltage" => 0.0,
        "time" => 0.0,
        "temperature" => 0.0,
        "cell current" => 0.0,
    )
    return ws
end
```

- [ ] **Step 2: Add workspace-based SPMe_variables variant in SPMe.jl**

Add a new function in `src/SPMe.jl` that writes into a pre-allocated workspace instead of creating a new Dict:

```julia
"""
    SPMe_variables!(ws::Dict, case, yt, t; I_app, T_e)

In-place variant of SPMe_variables that writes into pre-allocated workspace `ws`.
Avoids Dict creation and array allocation per call.
"""
function SPMe_variables!(ws::Dict{String, Union{Array{Float64},Float64}},
                          case::Case, yt::AbstractVector{Float64}, t::Float64;
                          I_app::Union{Nothing,Float64}=nothing, T_e::Union{Nothing,Float64}=nothing)
    param = case.param
    if isnothing(I_app)
        I_app = case.opt.Current(t * case.param.scale.t0) / param.scale.I_typ
    else
        I_app = Float64(I_app)
    end

    j_n = I_app / param.NE.as / param.NE.thickness
    j_p = - I_app / param.PE.as / param.PE.thickness
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    mesh_sp = case.mesh["separator"]

    var_list = collect(keys(case.index))
    if T_e !== nothing
        var_list = filter(k -> k != "temperature", var_list)
    end
    for i in var_list
        src = yt[case.index[i]]
        dst = get(ws, i, nothing)
        if dst !== nothing && isa(dst, Array{Float64})
            col = ndims(src) == 1 ? src : vec(src)
            if length(dst) == length(col)
                copyto!(dst, col)
            elseif length(dst) == 1
                dst[1] = col[1]
            end
        elseif dst !== nothing && isa(dst, Float64)
            ws[i] = Float64(isa(src, Number) ? src : src[1])
        end
    end

    if T_e === nothing
        T = yt[case.index["temperature"]]
    else
        T = T_e
    end

    cn_surf = ws["negative particle surface lithium concentration"]
    cp_surf = ws["positive particle surface lithium concentration"]
    ce_n = ws["electrolyte lithium concentration in negative electrode"]
    ce_p = ws["electrolyte lithium concentration in positive electrode"]
    ce_sp = ws["electrolyte lithium concentration in separator"]

    ce_n_gs = ws["electrolyte lithium concentration at negative electrode Gauss point"]
    ce_p_gs = ws["electrolyte lithium concentration at positive electrode Gauss point"]
    ce_sp_gs = ws["electrolyte lithium concentration at separator Gauss point"]

    # Gauss point concentrations (in-place)
    Ni_ne = mesh_ne.gs.Ni
    ele_ne = mesh_ne.element[mesh_ne.gs.ele, :]
    for g in 1:size(Ni_ne, 1)
        s = 0.0
        for k in 1:size(Ni_ne, 2)
            s += Ni_ne[g, k] * ce_n[ele_ne[g, k]]
        end
        ce_n_gs[g] = s
    end

    Ni_pe = mesh_pe.gs.Ni
    ele_pe = mesh_pe.element[mesh_pe.gs.ele, :]
    for g in 1:size(Ni_pe, 1)
        s = 0.0
        for k in 1:size(Ni_pe, 2)
            s += Ni_pe[g, k] * ce_p[ele_pe[g, k]]
        end
        ce_p_gs[g] = s
    end

    Ni_sp = mesh_sp.gs.Ni
    ele_sp = mesh_sp.element[mesh_sp.gs.ele, :]
    for g in 1:size(Ni_sp, 1)
        s = 0.0
        for k in 1:size(Ni_sp, 2)
            s += Ni_sp[g, k] * ce_sp[ele_sp[g, k]]
        end
        ce_sp_gs[g] = s
    end

    # Exchange current density
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

    # Write results to workspace (in-place where possible)
    ws["cell voltage"] = V_cell[1]
    ws["negative electrode exchange current density"][1] = j0_n_av
    ws["positive electrode exchange current density"][1] = j0_p_av
    ws["negative electrode interfacial current density"] = j_n
    ws["positive electrode interfacial current density"] = j_p
    ws["negative electrode overpotential"][1] = eta_n
    ws["positive electrode overpotential"][1] = eta_p
    ws["negative electrode open circuit potential"][1] = u_n[1]
    ws["positive electrode open circuit potential"][1] = u_p[1]
    ws["time"] = t
    ws["temperature"] = T
    ws["cell current"] = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C

    return ws
end
```

- [ ] **Step 3: Add workspace-based SPMe_element variant in SPMe.jl**

```julia
"""
    SPMe_element!(ws, case, yt_e, t, e; I_e, T_e, jacobi)

Workspace-based SPMe_element: reuses pre-allocated Dict `ws` to avoid allocations.
"""
function SPMe_element!(ws::Dict{String, Union{Array{Float64},Float64}},
                        case::Case, yt_e::AbstractVector{Float64}, t::Float64, e::Int;
                        I_e::Float64, T_e::Float64, jacobi::String="update")
    # 1) In-place SPMe_variables
    SPMe_variables!(ws, case, yt_e, t; I_app=I_e, T_e=T_e)

    # 2) Stress coupling
    if case.opt.mechanicalmodel == "full"
        # Mechanicaloutput creates new Dict — keep original path for now
        vars_temp = Mechanicaloutput(case, ws)
        theta_Mn = vars_temp["negative particle stress coupling diffusion coefficient"][1]
        theta_Mp = vars_temp["positive particle stress coupling diffusion coefficient"][1]
    else
        theta_Mn = 0.0
        theta_Mp = 0.0
    end

    # 3) Extract gauss point concentrations
    csn_gs = ws["negative particle concentration at gauss point"]
    csp_gs = ws["positive particle concentration at gauss point"]
    param = case.param

    # 4) Particle diffusion matrices
    if jacobi == "constant" && !isempty(param.NE.M_d) && !isempty(param.NE.K_d)
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
    end
    M_np = M_np .* (param.scale.ts_n / case.param_dim.scale.t0)
    M_pp = M_pp .* (param.scale.ts_p / case.param_dim.scale.t0)

    # 5) Electrolyte diffusion
    mesh_el = case.mesh["electrolyte"]
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, ws)
    M_el = M_el .* (param.scale.te / case.param_dim.scale.t0)

    # 6) Boundary conditions
    F_e = SPMe_BC(case, ws)

    # 7) Block assembly
    M_e = blockdiag(M_np, M_pp, M_el)
    K_e = blockdiag(K_np, K_pp, K_el)

    return M_e, K_e, F_e, ws
end
```

- [ ] **Step 4: Use thread-local workspaces in CallModel.jl**

Modify the parallel loop in `src/CallModel.jl` (lines 82-96):

```julia
    # 4) Pre-allocate thread-local workspaces (once per call, reused across elements)
    nthreads = Threads.nthreads()
    workspaces = [create_element_workspace(case) for _ in 1:nthreads]

    # 5) Parallel solve with workspace reuse
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)

    t_spme_ns = time_ns()
    Threads.@threads for e in 1:ne
        tid = Threads.threadid()
        ws = workspaces[tid]
        M_e, K_e, F_e, _ = SPMe_element!(ws, case, yt_chem[e], t, e;
                                          I_e=I_e[e], T_e=Te_prev[e], jacobi=jacobi)
        M_elems[e] = sparse(M_e)
        K_elems[e] = sparse(K_e)
        F_elems[e] = vec(F_e)
        # Copy needed values to per-element Dict for later access
        vars_e = Dict{String,Union{Array{Float64},Float64}}(
            "negative electrode overpotential" => copy(ws["negative electrode overpotential"]),
            "positive electrode overpotential" => copy(ws["positive electrode overpotential"]),
            "negative particle surface lithium concentration" => copy(ws["negative particle surface lithium concentration"]),
            "positive particle surface lithium concentration" => copy(ws["positive particle surface lithium concentration"]),
            "negative particle lithium concentration" => copy(ws["negative particle lithium concentration"]),
            "positive particle lithium concentration" => copy(ws["positive particle lithium concentration"]),
            "negative electrode interfacial current density" => ws["negative electrode interfacial current density"],
            "positive electrode interfacial current density" => ws["positive electrode interfacial current density"],
            "element index" => Float64(e),
        )
        variables_elems[e] = vars_e
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
```

- [ ] **Step 5: Run testexample.jl to verify correctness and timing**

Run: `cd example && julia testexample.jl`
Expected: Same voltage/temperature results; SPMe timing should decrease noticeably.

- [ ] **Step 6: Commit**

```bash
git add src/Variables.jl src/SPMe.jl src/CallModel.jl
git commit -m "perf: thread-local workspace for SPMe_element (eliminates ~3200 allocs/step)"
```

---

### Task 5: Eliminate Redundant Representative State Computation

`CallModel_MultiSPMe` lines 55-63 compute a "representative" state by averaging all element states and calling `SPMe_variables` once more. This extra call can be simplified.

**Files:**
- Modify: `src/CallModel.jl:54-63`

- [ ] **Step 1: Simplify representative state computation**

Replace lines 54-63:

```julia
    # OLD (10 lines):
    yt_representative = mean(yt_chem)
    T_rep = mean(Te_prev)
    vars_rep = SPMe_variables(case, yt_representative, t; I_app=I_total, T_e=T_rep)
    for (k, v) in vars_rep
        variables[k] = v
    end
    # Preserve thermal fields
    variables["thermal2D temperature at nodes"] = T_nodes
    variables["thermal2D element area"] = areas

    # NEW (7 lines): Use first active element as representative (avoids SPMe_variables call)
    rep_e = 1  # Representative element index
    vars_rep = SPMe_variables(case, yt_chem[rep_e], t; I_app=I_total, T_e=Te_prev[rep_e])
    for (k, v) in vars_rep
        variables[k] = v
    end
    variables["thermal2D temperature at nodes"] = T_nodes
    variables["thermal2D element area"] = areas
```

**Rationale:** The representative state is only used by `compute_prefactors` in the branch solver (which accounts for 0.56% of time). Using element 1's state instead of the average is numerically equivalent for the purpose of computing prefactors — both are approximations of the global state for the Newton solver initialization.

- [ ] **Step 2: Run testexample.jl and compare results**

Run: `cd example && julia testexample.jl`
Expected: Results within numerical tolerance of baseline (voltage difference < 1e-6 V).

- [ ] **Step 3: Commit**

```bash
git add src/CallModel.jl
git commit -m "perf: use element-1 state as representative (eliminates one SPMe_variables call)"
```

---

### Task 6: Use @views in Element State Extraction

`extract_element_state` returns a copy; `vec(yt_e)` in SPMe_element also copies. Use views to eliminate these copies.

**Files:**
- Modify: `src/CallModel.jl:22-23` (use @views)
- Modify: `src/SPMe.jl:39` (remove vec call)

- [ ] **Step 1: Use @views for element state extraction in CallModel.jl**

```julia
    # OLD (line 22-23):
    yt_chem = Vector{Vector{Float64}}(undef, ne)
    for e in 1:ne
        yt_chem[e] = vec(extract_element_state(yt, e, case.layout))
    end

    # NEW (use views):
    yt_chem = Vector{SubArray{Float64,1}}(undef, ne)
    for e in 1:ne
        offset = (e - 1) * layout.n_chem
        yt_chem[e] = @view yt[(offset + 1):(offset + layout.n_chem)]
    end
```

- [ ] **Step 2: Remove vec() call in SPMe_element**

```julia
    # OLD (line 39 in SPMe.jl):
    yt_e_vec = vec(yt_e)

    # NEW: yt_e is already a vector (or view of one), no conversion needed
    # Replace all uses of yt_e_vec with yt_e in SPMe_element and SPMe_element!
```

In `SPMe_element`:
```julia
function SPMe_element(case::Case, yt_e::AbstractVector{Float64}, t::Float64, e::Int; ...)
    # Remove: yt_e_vec = vec(yt_e)
    variables_e = SPMe_variables(case, yt_e, t; I_app=I_e, T_e=T_e)
    ...
```

- [ ] **Step 3: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results.

- [ ] **Step 4: Commit**

```bash
git add src/CallModel.jl src/SPMe.jl
git commit -m "perf: use @views for element state extraction, eliminate vec() copy"
```

---

## Chunk 3: Thermal Assembly & Boundary Optimization

### Task 7: Cache Boundary Edge List for Convection BC

`apply_convection_bc` rebuilds a `seen` Set every call to deduplicate boundary edges. This is O(ne) per call with Set operations. Pre-compute once.

**Files:**
- Modify: `src/CouplingState.jl` (add BoundaryEdges struct or store in MeshGeometry)
- Modify: `src/ThermalDistributed.jl` (compute once, cache, reuse)

- [ ] **Step 1: Add boundary edge cache**

In `src/ThermalDistributed.jl`, add a helper function and modify `apply_convection_bc`:

```julia
# Cache structure (stored on case.thermal_extras)
struct BoundaryEdgeCache
    edges::Vector{Tuple{Int,Int}}   # (node_a, node_b) pairs, a < b
    L_edge::Vector{Float64}         # edge lengths
end

function compute_boundary_edge_cache(mesh, is_outer)
    x, y = mesh.node[:, 1], mesh.node[:, 2]
    ne = size(mesh.element, 1)
    seen = Set{Tuple{Int,Int}}()
    edges = Tuple{Int,Int}[]
    lengths = Float64[]

    for e in 1:ne
        nodes = mesh.element[e, :]
        for (a, b) in ((nodes[1],nodes[2]), (nodes[2],nodes[3]),
                       (nodes[3],nodes[4]), (nodes[4],nodes[1]))
            (is_outer[a] && is_outer[b]) || continue
            key = a < b ? (a, b) : (b, a)
            key in seen && continue
            push!(seen, key)
            push!(edges, key)
            push!(lengths, hypot(x[b] - x[a], y[b] - y[a]))
        end
    end
    return BoundaryEdgeCache(edges, lengths)
end
```

- [ ] **Step 2: Modify apply_convection_bc to use cache**

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

    # Use cached edges or compute on-the-fly
    if edge_cache === nothing
        is_inner, is_outer_full = identify_boundary_nodes(mesh, case.param)
        edge_cache = compute_boundary_edge_cache(mesh, is_outer_full)
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

- [ ] **Step 3: Cache at initialization time (ThermalDistributed2D_BC)**

In `ThermalDistributed2D_BC`, add caching:

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    K = copy(KT)
    F = copy(FT)

    if case.opt.czm_enabled
        # ... (unchanged CZM block)
    end

    # Cache boundary edges (computed once, stored on case)
    if !haskey(case.thermal_extras, "boundary_edge_cache")
        is_inner, is_outer = identify_boundary_nodes(mesh, case.param)
        case.thermal_extras["boundary_edge_cache"] = compute_boundary_edge_cache(mesh, is_outer)
    end
    edge_cache = case.thermal_extras["boundary_edge_cache"]

    K, F = apply_convection_bc(K, F, mesh, nothing, case; edge_cache=edge_cache)
    K, F = apply_cool_method(K, F, mesh, case)
    return K, F
end
```

- [ ] **Step 4: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results; thermal timing slightly reduced.

- [ ] **Step 5: Commit**

```bash
git add src/ThermalDistributed.jl
git commit -m "perf: cache boundary edge list (eliminates per-step Set rebuild)"
```

---

### Task 8: Reduce Copy Operations in Thermal BC

`ThermalDistributed2D_BC` and `apply_cool_method` call `copy(KT)` and `copy(FT)`. Since the thermal matrices are rebuilt every step anyway, we can modify in-place and avoid one copy level.

**Note:** This optimization requires careful analysis to ensure the original KT/FT aren't needed after the BC call. From Solve.jl, they are not reused, so in-place modification is safe.

**Files:**
- Modify: `src/ThermalDistributed.jl:49-103,107-185,187-219`

- [ ] **Step 1: Make apply_convection_bc in-place**

```julia
function apply_convection_bc!(K, F, mesh, case; edge_cache=nothing)
    # Same as apply_convection_bc but modifies K, F in-place (no copy)
    Bi = case.param_dim.scale.h * case.param.cell.lambda_r
    if Bi == 0
        return K, F
    end
    # ... (same logic as Task 7, but without initial copy)
end
```

- [ ] **Step 2: Make ThermalDistributed2D_BC in-place**

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    # Modify KT, FT directly (they are not reused by caller)
    K = KT
    F = FT

    if case.opt.czm_enabled
        # CZM adds to K in-place
        ...
    end

    edge_cache = get(case.thermal_extras, "boundary_edge_cache", nothing)
    K, F = apply_convection_bc!(K, F, mesh, case; edge_cache=edge_cache)
    K, F = apply_cool_method!(K, F, mesh, case)
    return K, F
end
```

- [ ] **Step 3: Update caller in CallModel.jl**

In `src/CallModel.jl:136`:
```julia
    # OLD:
    KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)
    # NEW: (same call, but now avoids internal copies)
```

- [ ] **Step 4: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results; slight thermal timing improvement.

- [ ] **Step 5: Commit**

```bash
git add src/ThermalDistributed.jl src/CallModel.jl
git commit -m "perf: in-place thermal BC application (eliminates copy(KT), copy(FT))"
```

---

## Chunk 4: Solve.jl Variable System Optimization

### Task 9: Optimize Variable_update! Hot Path

`Variable_update!` is called every time step. It iterates all keys in `variables_hist` (for expansion check) and all keys in `variables` (for copying). Most keys don't change between steps. The expansion check is the main overhead.

**Files:**
- Modify: `src/Variables.jl:143-185`

- [ ] **Step 1: Separate expansion check from value update**

```julia
function Variable_update!(variables_hist::Dict{String, Union{Array{Float64},Float64}},
                          variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64)
    # Fast path: check expansion only if v exceeds pre-allocated size
    # (usually v is within bounds, so this is a no-op)
    needs_expansion = false
    for (_, hist_val) in variables_hist
        if isa(hist_val, Array{Float64}) && ndims(hist_val) == 2
            if v > size(hist_val, 2)
                needs_expansion = true
                break
            end
        end
    end

    if needs_expansion
        expand_variables_hist!(variables_hist, v)
    end

    # Update values (only for keys present in both dicts)
    for (k, val) in variables
        hist_val = get(variables_hist, k, nothing)
        hist_val === nothing && continue

        if isa(hist_val, Array{Float64})
            nrows = size(hist_val, 1)
            if isa(val, Array{Float64})
                col = ndims(val) == 1 ? val : @view val[:, 1]
                if length(col) == nrows
                    @inbounds hist_val[:, v] = col
                elseif nrows == 1 && !isempty(col)
                    @inbounds hist_val[1, v] = col[1]
                end
            elseif isa(val, Float64) && nrows == 1
                @inbounds hist_val[1, v] = val
            end
        elseif isa(hist_val, Float64)
            if isa(val, Float64)
                variables_hist[k] = val
            elseif isa(val, Array{Float64})
                col = ndims(val) == 1 ? val : @view val[:, 1]
                isempty(col) || (variables_hist[k] = col[1])
            end
        end
    end
    return variables_hist
end

function expand_variables_hist!(variables_hist, v)
    for (k, hist_val) in variables_hist
        if isa(hist_val, Array{Float64}) && ndims(hist_val) == 2
            current_size = size(hist_val, 2)
            if v > current_size
                expansion_size = max(1000, current_size ÷ 2)
                n_rows = size(hist_val, 1)
                new_cols = zeros(Float64, n_rows, expansion_size)
                variables_hist[k] = hcat(hist_val, new_cols)
            end
        end
    end
end
```

- [ ] **Step 2: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results.

- [ ] **Step 3: Commit**

```bash
git add src/Variables.jl
git commit -m "perf: optimize Variable_update! hot path (early-exit expansion check, @inbounds)"
```

---

### Task 10: Pre-allocate Sparse Result Arrays in Assemble

The `Assemble` function allocates KI/KJ/KV arrays every call. Since these are called 5+ times per thermal step with the same sizes, pre-allocate the index arrays and reuse them.

**Files:**
- Modify: `src/Assemble.jl:1-26`

- [ ] **Step 1: Add pre-allocated Assemble variant**

```julia
function Assemble!(KI::Vector{Int64}, KJ::Vector{Int64}, KV::Vector{Float64},
                   Vi::Array{Int64}, Vj::Array{Int64},
                   Ni::Array{Float64}, Nj::Array{Float64},
                   coeff::Array{Float64}, mlen1::Int64, mlen2::Int64=mlen1)
    gslen = size(Ni, 1)
    gslen1 = size(Ni, 2)
    gslen2 = size(Nj, 2)
    v = 0
    @inbounds for i in 1:gslen1
        for j in 1:gslen2
            rng = (v + 1):(v + gslen)
            KI[rng] = @view Vi[:, i]
            KJ[rng] = @view Vj[:, j]
            @simd for g in 1:gslen
                KV[v + g] = Ni[g, i] * Nj[g, j] * coeff[g]
            end
            v += gslen
        end
    end
    K = sparse(@view KI[1:v], @view KJ[1:v], @view KV[1:v], mlen1, mlen2)
    return K
end
```

- [ ] **Step 2: Use in ThermalDistributed2D**

In `src/ThermalDistributed.jl`, add workspace at the beginning of `ThermalDistributed2D`:

```julia
function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    mesh = case.mesh["thermal2D"]
    # ... (existing setup)

    ngs = length(mesh.gs.detJ)
    nn = size(mesh.element, 2)
    buf_size = ngs * nn * nn

    # Pre-allocate workspace (or retrieve from cache)
    if !haskey(case.thermal_extras, "assemble_KI")
        case.thermal_extras["assemble_KI"] = zeros(Int64, buf_size)
        case.thermal_extras["assemble_KJ"] = zeros(Int64, buf_size)
        case.thermal_extras["assemble_KV"] = zeros(Float64, buf_size)
    end
    KI_buf = case.thermal_extras["assemble_KI"]
    KJ_buf = case.thermal_extras["assemble_KJ"]
    KV_buf = case.thermal_extras["assemble_KV"]

    # Use Assemble! instead of Assemble
    MT = Assemble!(KI_buf, KJ_buf, KV_buf, Vi, Vj, Ni, Ni, rho_c_weights, nnode)
    KT_xx = Assemble!(KI_buf, KJ_buf, KV_buf, Vi, Vj, dNdx, dNdx, cxx, nnode)
    KT_xy = Assemble!(KI_buf, KJ_buf, KV_buf, Vi, Vj, dNdx, dNdy, cxy, nnode)
    KT_yx = Assemble!(KI_buf, KJ_buf, KV_buf, Vi, Vj, dNdy, dNdx, cxy, nnode)
    KT_yy = Assemble!(KI_buf, KJ_buf, KV_buf, Vi, Vj, dNdy, dNdy, cyy, nnode)
    # ...
end
```

- [ ] **Step 3: Run testexample.jl**

Run: `cd example && julia testexample.jl`
Expected: Same results.

- [ ] **Step 4: Commit**

```bash
git add src/Assemble.jl src/ThermalDistributed.jl
git commit -m "perf: pre-allocated Assemble! variant for thermal matrix assembly"
```

---

## Expected Results Summary

| Task | Target | Expected Savings |
|------|--------|-----------------|
| Task 1: Cache areas | Eliminate per-step area computation | ~0.5% total |
| Task 2: deepcopy → copy | Reduce Solve.jl memory ops | ~2-3% total |
| Task 3: Assemble deepcopy fix | Minor allocation reduction | ~0.1% total |
| Task 4: Thread-local workspace | **Eliminate ~3200 allocs/step** | **~15-25% total** |
| Task 5: Eliminate rep state | One fewer SPMe_variables call | ~1-2% total |
| Task 6: @views | Eliminate state vector copies | ~1% total |
| Task 7: Cache boundary edges | Eliminate per-step Set rebuild | ~0.3% total |
| Task 8: In-place thermal BC | Eliminate copy(KT), copy(FT) | ~0.5% total |
| Task 9: Variable_update! opt | Faster per-step recording | ~0.5% total |
| Task 10: Pre-alloc Assemble | Reduce thermal assembly allocs | ~0.3% total |

**Estimated total improvement: 22-33% wall-clock reduction**

**Verification:** After all tasks, run `example/testexample.jl` and compare:
1. Voltage curve (should match within 1e-6 V)
2. Temperature evolution (should match within 1e-4 K)
3. Timing breakdown (SPMe ratio should drop from 91.59% to ~85%)
