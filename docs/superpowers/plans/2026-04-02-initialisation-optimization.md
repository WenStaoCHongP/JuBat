# Initialisation.jl 类型安全重构 实施计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 `multi_spme_layout::Dict{String,Any}` 替换为类型安全的 struct 字段，消除所有 Dict 访问。

**Architecture:** 新建 `CouplingState.jl` 定义 `MultiSPMeLayout` 和 `MeshGeometry` 两个不可变 struct。修改 `Case` 结构体新增 `layout`、`geometry`、`czm_mesh`、`thermal_extras` 四个字段。重写 Initialisation.jl 中 4 个辅助函数为 snake_case 签名。连锁更新 Solve.jl、CycleSolver.jl、CycleData.jl、PostProcessing.jl、Jellyrollmodel.jl、ThermalDistributed.jl。

**Tech Stack:** Julia 1.x, 无外部依赖

**Fail-Fast 约束:** 禁止 `haskey`、`get(dict, key, default)`、`try/catch` 吞错误。直接访问，未初始化就让 Julia 抛异常。

---

## File Structure

| 文件 | 操作 | 职责 |
|------|------|------|
| `src/CouplingState.jl` | **新建** | `MultiSPMeLayout` + `MeshGeometry` struct 定义 |
| `src/SetCase.jl` | 修改 | Case struct 字段变更 + 兼容构造器 |
| `src/Initialisation.jl` | 重写 | 4 个辅助函数（snake_case + 直接传 layout） |
| `src/JuBat.jl` | 修改 | include + export 更新 |
| `src/Solve.jl` | 修改 | ~20 处 Dict 访问替换 |
| `src/CycleSolver.jl` | 修改 | 删除 `_ensure_multi_spme_layout!`，替换访问 |
| `src/CycleData.jl` | 修改 | 4 处 Dict 访问替换 |
| `src/PostProcessing.jl` | 修改 | 3 处 Dict 访问替换（规格遗漏） |
| `src/Jellyrollmodel.jl` | 修改 | 几何信息写入 `case.geometry` + layer_weights |
| `src/ThermalDistributed.jl` | 修改 | haskey → 直接访问 `case.czm_mesh` |

---

## Chunk 1: Foundation（结构体定义 + Case 重构）

### Task 1: 创建 CouplingState.jl

**Files:**
- Create: `src/CouplingState.jl`

- [ ] **Step 1: 创建文件并定义结构体**

```julia
# src/CouplingState.jl
"""
    CouplingState.jl — 类型安全的状态布局与网格几何定义
    替代 Dict{String,Any} 的 multi_spme_layout
"""

"""
多SPMe状态向量的布局索引。初始化后不可变。
"""
struct MultiSPMeLayout
    ne::Int                        # 热单元数
    n_chem::Int                    # 每单元电化学 DOF 数 = Nrn + Nrp + Nel
    nT::Int                        # 热节点 DOF 数
    n_total::Int                   # 全局状态向量总长 = ne*n_chem + nT
    chem_range::UnitRange{Int}     # 1:(ne*n_chem)
    thermal_range::UnitRange{Int}  # (ne*n_chem+1):(ne*n_chem+nT)
end

"""便捷构造器：自动计算 range 和 n_total"""
function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)
    MultiSPMeLayout(
        ne, n_chem, nT,
        ne * n_chem + nT,
        1:(ne * n_chem),
        (ne * n_chem + 1):(ne * n_chem + nT)
    )
end

"""
Jellyroll 网格的几何拓扑信息。构建后不可变。
"""
struct MeshGeometry
    element_layer::Vector{Int}                      # 每个单元的层类型 (1=NE, 2=SP, 3=PE, 4=NCC, 5=PCC)
    is_inner_layer::Vector{Bool}                    # 是否为内层
    layer_weights::Matrix{Float64}                  # ne × 5 层面积权重 [NE, SP, PE, PCC, NCC]
    interface_pairs::Vector{Tuple{Int,Int}}         # CZM 界面配对 (top_elem, bot_elem)
    czm_element_map::Dict{Int,Vector{Int}}          # 热单元号 → CZM 单元索引向量（一对多映射）
    inner_nodes::Vector{Int}                        # 内边界节点索引
    outer_nodes::Vector{Int}                        # 外边界节点索引
end
```

- [ ] **Step 2: 提交**

```bash
git add src/CouplingState.jl
git commit -m "feat: add CouplingState.jl with MultiSPMeLayout and MeshGeometry structs"
```

### Task 2: 重构 Case 结构体

**Files:**
- Modify: `src/SetCase.jl:92-105`

- [ ] **Step 1: 替换 Case struct 定义**

将 SetCase.jl 行 98-105 替换为：

```julia
mutable struct Case
    param_dim::Params                      # dimensional parameters
    param::Params                          # dimensionless parameters
    opt::Option                            # solver options
    mesh::Dict{String, Mesh}               # discretisation meshes
    index::Dict{String, Union{Array{Int64}, Int64}} # indices of unknowns
    layout::Union{Nothing, MultiSPMeLayout}   # 布局索引（初始化后不变）
    geometry::Union{Nothing, MeshGeometry}    # 几何拓扑（构建后不变）
    czm_mesh::Union{Nothing, CohesiveMesh}    # CZM 网格（演化但类型明确）
    thermal_extras::Dict{String,Any}          # 运行时状态暂存
end

# 5 参数兼容构造器（main 分支标准路径）
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, Dict{String,Any}())
end
```

- [ ] **Step 2: 更新 SetCase 函数调用**

将 SetCase.jl 行 93 的 `Dict{String,Any}()` 构造改为使用兼容构造器：

```julia
# 旧: case = Case(param_dim, param, opt, mesh, index, Dict{String,Any}())
# 新:
case = Case(param_dim, param, opt, mesh, index)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetCase.jl
git commit -m "refactor: replace multi_spme_layout Dict with typed fields in Case struct"
```

### Task 3: 更新 JuBat.jl 模块入口

**Files:**
- Modify: `src/JuBat.jl:4-38`

- [ ] **Step 1: 添加 CouplingState.jl include**

在 JuBat.jl 行 6（`include("SetCase.jl")`）之后添加：

```julia
include("CouplingState.jl")  # 类型安全的状态布局与网格几何
```

- [ ] **Step 2: 更新 export 列表**

将 JuBat.jl 行 37-38 替换为：

```julia
export ModelInitialisation_MultiSPMe, extract_element_state, get_thermal_dofs
export update_state
```

- [ ] **Step 3: 提交**

```bash
git add src/JuBat.jl src/CouplingState.jl src/SetCase.jl
git commit -m "feat: wire up CouplingState.jl in module, update exports"
```

---

## Chunk 2: Core（Initialisation.jl 重写）

### Task 4: 重写 Initialisation.jl 辅助函数

**Files:**
- Modify: `src/Initialisation.jl:49-240`

- [ ] **Step 1: 替换 ModelInitialisation_MultiSPMe 中的布局缓存**

将 Initialisation.jl 行 110-116 替换为：

```julia
# 旧:
#   empty!(case.multi_spme_layout)
#   case.multi_spme_layout["ne"] = ne
#   case.multi_spme_layout["n_chem"] = n_chem
#   case.multi_spme_layout["nT"] = nT
#   case.multi_spme_layout["n_total"] = length(y0)
#   case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
#   case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)

# 新 (1 行):
case.layout = MultiSPMeLayout(ne, n_chem, nT)
```

- [ ] **Step 2: 替换 MultiSPMe_extract_element_state**

将 Initialisation.jl 行 122-161 的整个函数（含文档字符串）替换为：

```julia
"""
    extract_element_state(y, e, layout) -> Vector{Float64}

从多SPMe全局状态向量中提取单个单元的电化学状态。
"""
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end
```

- [ ] **Step 3: 替换 MultiSPMe_get_thermal_dofs**

将 Initialisation.jl 行 164-193 的整个函数（含文档字符串）替换为：

```julia
"""
    get_thermal_dofs(y, layout) -> Vector{Float64}

从多SPMe全局状态向量中提取热场节点温度。
"""
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end
```

- [ ] **Step 4: 替换 MultiSPMe_update_state**

将 Initialisation.jl 行 196-239 的整个函数（含文档字符串）替换为：

```julia
"""
    update_state(y, layout; element_index, element_state, thermal_nodes) -> Vector{Float64}

更新多SPMe全局状态向量（返回新向量）。
"""
function update_state(y::AbstractVector, layout::MultiSPMeLayout;
                      element_index::Union{Nothing,Int}=nothing,
                      element_state::Union{Nothing,Vector{Float64}}=nothing,
                      thermal_nodes::Union{Nothing,Vector{Float64}}=nothing)
    y_new = copy(y)
    if element_index !== nothing
        @assert 1 <= element_index <= layout.ne "element_index $element_index out of range [1, $(layout.ne)]"
        @assert length(element_state) == layout.n_chem "element_state length $(length(element_state)) != n_chem $(layout.n_chem)"
        offset = (element_index - 1) * layout.n_chem
        y_new[(offset + 1):(offset + layout.n_chem)] .= element_state
    end
    if thermal_nodes !== nothing
        @assert length(thermal_nodes) == layout.nT "thermal_nodes length $(length(thermal_nodes)) != nT $(layout.nT)"
        y_new[layout.thermal_range] .= thermal_nodes
    end
    return y_new
end
```

- [ ] **Step 5: 提交**

```bash
git add src/Initialisation.jl
git commit -m "refactor: rewrite Initialisation helpers to typed layout parameter"
```

---

## Chunk 3: Consumers（连锁修改调用方）

### Task 5: 更新 Solve.jl

**Files:**
- Modify: `src/Solve.jl` — 多处改动

Solve.jl 中有最多的 Dict 访问需要替换。逐个替换：

- [ ] **Step 1: 替换运行时键访问（行 6-8）**

```julia
# 旧:
vars = case.multi_spme_layout["thermal_variables"]
update_fn = case.multi_spme_layout["thermal_update_fn"]
record = case.multi_spme_layout["thermal_record"]

# 新:
vars = case.thermal_extras["thermal_variables"]
update_fn = case.thermal_extras["thermal_update_fn"]
record = case.thermal_extras["thermal_record"]
```

- [ ] **Step 2: 替换 polar_mesh_data（行 30, 45）**

```julia
# 旧: case.multi_spme_layout["polar_mesh_data"]
# 新: case.thermal_extras["polar_mesh_data"]
```

- [ ] **Step 3: 替换布局构建（行 103-109）**

```julia
# 旧:
if isempty(case.multi_spme_layout)
    case.multi_spme_layout["ne"] = ne
    case.multi_spme_layout["n_chem"] = n_chem
    case.multi_spme_layout["nT"] = nT
    case.multi_spme_layout["n_total"] = expected_multi_len
    case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
    case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)
end

# 新:
if case.layout === nothing
    case.layout = MultiSPMeLayout(ne, n_chem, nT)
end
```

- [ ] **Step 4: 替换 CallModel_MultiSPMe 中的访问（行 432-452）**

```julia
# 行 432-433 旧:
if isempty(case.multi_spme_layout)
    error("CallModel_MultiSPMe requires populated multi_spme_layout...")
end
# 新:
if case.layout === nothing
    error("CallModel_MultiSPMe requires case.layout. Did you call ModelInitialisation_MultiSPMe?")
end

# 行 436 旧: layout = case.multi_spme_layout
# 新: layout = case.layout

# 行 448 旧: T_nodes = MultiSPMe_get_thermal_dofs(yt, case)
# 新: T_nodes = get_thermal_dofs(yt, case.layout)

# 行 452 旧: yt_chem[e] = vec(MultiSPMe_extract_element_state(yt, e, case))
# 新: yt_chem[e] = vec(extract_element_state(yt, e, case.layout))
```

- [ ] **Step 5: 替换自动初始化（行 584-608）**

```julia
# 行 584 旧: if should_use_multi_spme && isempty(case.multi_spme_layout)
# 新: if should_use_multi_spme && case.layout === nothing

# 行 597-603 旧: 6 行 Dict 赋值
# 新:
@warn "CallModel: layout is nothing but state vector matches multi-SPMe format, auto-initializing layout"
case.layout = MultiSPMeLayout(ne, n_chem, nT)

# 行 608 旧: multi_spme_enabled = should_use_multi_spme && !isempty(case.multi_spme_layout)
# 新: multi_spme_enabled = should_use_multi_spme && case.layout !== nothing
```

- [ ] **Step 6: 删除 layer_weights 的 try/catch 存储（行 171-178）**

```julia
# 旧:
if case.opt.collector_seeded
    try
        fks = jellyroll_element_properties(case.mesh["thermal2D"], case.param)[2]
        variables["thermal2D layer_weights"] = fks
    catch err
        @warn "Failed to set layer_weights: $err"
    end
end

# 新: 删除整个 if 块（layer_weights 已在 MeshGeometry 中管理）
```

- [ ] **Step 7: 提交**

```bash
git add src/Solve.jl
git commit -m "refactor: replace all multi_spme_layout Dict access in Solve.jl"
```

### Task 6: 更新 CycleSolver.jl

**Files:**
- Modify: `src/CycleSolver.jl:688-713`

- [ ] **Step 1: 删除 `_ensure_multi_spme_layout!` 函数**

删除 CycleSolver.jl 行 686-713 的整个函数（包括文档字符串）。

- [ ] **Step 2: 提交**

```bash
git add src/CycleSolver.jl
git commit -m "refactor: remove _ensure_multi_spme_layout! from CycleSolver"
```

### Task 7: 更新 CycleData.jl

**Files:**
- Modify: `src/CycleData.jl:65,73,77,149`

- [ ] **Step 1: 替换 4 处访问**

```julia
# 行 65 旧: _ensure_multi_spme_layout!(case)
# 新: (删除此行，或替换为)
if case.layout === nothing
    ne_cd = size(case.mesh["thermal2D"].element, 1)
    nT_cd = case.mesh["thermal2D"].nlen
    n_chem_cd = length(y_new) - nT_cd  # 从状态向量反推
    # 注意：需要从上下文获取 n_chem，或直接用 layout
    error("CycleData: case.layout is nothing, cannot extract thermal DOFs")
end

# 行 73 旧: thermal_range = case.multi_spme_layout["thermal_range"]
# 新: thermal_range = case.layout.thermal_range

# 行 77 旧: thermal_range = case.multi_spme_layout["thermal_range"]
# 新: thermal_range = case.layout.thermal_range

# 行 149 旧: T_nodes_carry = MultiSPMe_get_thermal_dofs(y_new, case)
# 新: T_nodes_carry = get_thermal_dofs(y_new, case.layout)
```

注意：行 65 原来调用 `_ensure_multi_spme_layout!`，按 fail-fast 原则，如果 layout 为 nothing 应直接报错而非静默构建。但需查看 CycleData.jl 行 65 上下文确认是否有 `ne`/`n_chem`/`nT` 可用。如果没有，应直接报错。

- [ ] **Step 2: 提交**

```bash
git add src/CycleData.jl
git commit -m "refactor: replace multi_spme_layout access in CycleData.jl"
```

### Task 8: 更新 PostProcessing.jl

**Files:**
- Modify: `src/PostProcessing.jl:100-102`

- [ ] **Step 1: 替换 3 处访问**

```julia
# 行 100 旧: _ensure_multi_spme_layout!(case)
# 新: (删除此行，layout 应已初始化)

# 行 101 旧: ne = case.multi_spme_layout["ne"]
# 新: ne = case.layout.ne

# 行 102 旧: n_chem = case.multi_spme_layout["n_chem"]
# 新: n_chem = case.layout.n_chem
```

按 fail-fast 原则，如果 `case.layout` 为 nothing，Julia 会自动抛出 `UndefRefError`，这正是预期行为。

- [ ] **Step 2: 提交**

```bash
git add src/PostProcessing.jl
git commit -m "refactor: replace multi_spme_layout access in PostProcessing.jl"
```

---

## Chunk 4: Geometry（几何数据迁移）

### Task 9: 更新 Jellyrollmodel.jl

**Files:**
- Modify: `src/Jellyrollmodel.jl:546-551`

- [ ] **Step 1: 替换几何信息写入**

将 Jellyrollmodel.jl 行 546-551 的 6 行 Dict 赋值替换为：

```julia
# 旧:
case_new.multi_spme_layout["interface_pairs"] = interface_pairs
case_new.multi_spme_layout["element_layer"] = mesh_data.element_layer
case_new.multi_spme_layout["is_inner_layer"] = mesh_data.is_inner_layer
case_new.multi_spme_layout["czm_element_map"] = mesh_data.czm_element_map
case_new.multi_spme_layout["inner_nodes"] = mesh_data.inner_nodes
case_new.multi_spme_layout["outer_nodes"] = mesh_data.outer_nodes

# 新:
_, layer_weights = jellyroll_element_properties(case_new.mesh["thermal2D"], case_new.param)
case_new.geometry = MeshGeometry(
    mesh_data.element_layer,
    mesh_data.is_inner_layer,
    layer_weights,
    interface_pairs,
    mesh_data.czm_element_map,
    mesh_data.inner_nodes,
    mesh_data.outer_nodes
)
```

注意：如果此处 `case_new.mesh["thermal2D"]` 尚未设置（取决于 `setup_thermal2D_mesh` 的执行顺序），需要确认 `jellyroll_element_properties` 调用时机。查看上下文确认 `thermal2D` mesh 已在此时可用。

- [ ] **Step 2: 提交**

```bash
git add src/Jellyrollmodel.jl
git commit -m "refactor: write geometry data to MeshGeometry struct in Jellyrollmodel"
```

### Task 10: 更新 ThermalDistributed.jl

**Files:**
- Modify: `src/ThermalDistributed.jl:192-193`

- [ ] **Step 1: 替换 czm_mesh 访问**

```julia
# 旧:
if case.opt.czm_enabled && haskey(case.multi_spme_layout, "czm_mesh")
    czm_mesh = case.multi_spme_layout["czm_mesh"]

# 新:
if case.opt.czm_enabled
    czm_mesh = case.czm_mesh
```

- [ ] **Step 2: 提交**

```bash
git add src/ThermalDistributed.jl
git commit -m "refactor: replace haskey check with direct czm_mesh access in ThermalDistributed"
```

---

## Chunk 5: 验证

### Task 11: 全局验证

- [ ] **Step 1: 搜索确认无残留 Dict 访问**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat"
grep -rn "multi_spme_layout" src/
grep -rn "MultiSPMe_extract_element_state\|MultiSPMe_get_thermal_dofs\|MultiSPMe_update_state" src/
grep -rn "_ensure_multi_spme_layout" src/
```

预期：**零匹配**。如果有残留，回到对应 Task 修复。

- [ ] **Step 2: 搜索确认无 haskey 残留**

```bash
grep -rn "haskey" src/ --include="*.jl"
```

预期：无 `haskey.*multi_spme_layout` 或 `haskey.*czm_mesh` 匹配。其他 haskey 使用（如有）不在本次范围内。

- [ ] **Step 3: 搜索确认无 try/catch 吞错误**

```bash
grep -rn "try" src/Solve.jl
```

预期：Solve.jl 中 `layer_weights` 的 try/catch 已删除。

- [ ] **Step 4: 手动冒烟测试**

在 Julia REPL 中执行：

```julia
include("src/JuBat.jl")
using .JuBat

# 测试 Case 构造
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
case = JuBat.SetCase(param_dim, opt)

# 确认 layout/geometry/czm_mesh 初始为 nothing
@assert case.layout === nothing
@assert case.geometry === nothing
@assert case.czm_mesh === nothing

# 测试 MultiSPMeLayout 构造
layout = JuBat.MultiSPMeLayout(80, 31, 160)
@assert layout.ne == 80
@assert layout.n_chem == 31
@assert layout.nT == 160
@assert layout.n_total == 80*31 + 160
@assert layout.chem_range == 1:(80*31)
@assert layout.thermal_range == (80*31+1):(80*31+160)
```

- [ ] **Step 5: 最终提交**

```bash
git add -A
git commit -m "refactor: complete Initialisation.jl type-safety migration

Replace multi_spme_layout Dict{String,Any} with typed structs:
- MultiSPMeLayout for layout indices
- MeshGeometry for geometric topology (incl. layer_weights)
- case.czm_mesh for CZM mesh
- case.thermal_extras for runtime state keys (temporary)

Fail-fast: no haskey, no try/catch swallowing errors.
"
```
