# Initialisation.jl 优化设计规格

> 日期: 2026-04-02
> 状态: 待实施
> 基于文档: 03_结构体重构方案.md, 05_Initialisation_jl优化方案.md

---

## 1. 目标

对 Initialisation.jl 及其关联文件执行类型安全重构，消除 `Dict{String,Any}` 的类型不安全问题，统一命名规范，简化布局构建逻辑。

## 2. 前置依赖 — 新建 `src/CouplingState.jl`

定义两个不可变结构体，替代 `multi_spme_layout::Dict{String,Any}` 中的布局索引和网格几何数据。

### 2.1 `MultiSPMeLayout`

```julia
struct MultiSPMeLayout
    ne::Int                        # 热单元数
    n_chem::Int                    # 每单元电化学 DOF 数
    nT::Int                        # 热节点 DOF 数
    n_total::Int                   # 全局状态向量总长
    chem_range::UnitRange{Int}     # 1:(ne*n_chem)
    thermal_range::UnitRange{Int}  # (ne*n_chem+1):(ne*n_chem+nT)
end

function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int)
    MultiSPMeLayout(
        ne, n_chem, nT,
        ne * n_chem + nT,
        1:(ne * n_chem),
        (ne * n_chem + 1):(ne * n_chem + nT)
    )
end
```

### 2.2 `MeshGeometry`

```julia
struct MeshGeometry
    element_layer::Vector{Int}
    is_inner_layer::BitVector
    layer_weights::Matrix{Float64}
    interface_pairs::Vector{Tuple{Int,Int}}
    czm_element_map::Dict{Int,Int}
    inner_nodes::Vector{Int}
    outer_nodes::Vector{Int}
end
```

## 3. Case 结构体修改 (`src/SetCase.jl`)

### 3.1 字段变更

```julia
# 旧
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
    multi_spme_layout::Dict{String,Any}
end

# 新
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
    layout::Union{Nothing, MultiSPMeLayout}
    geometry::Union{Nothing, MeshGeometry}
    czm_mesh::Union{Nothing, CohesiveMesh}
end
```

### 3.2 兼容构造器

```julia
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing)
end
```

### 3.3 SetCase 函数内

```julia
# 旧: case = Case(param_dim, param, opt, mesh, index, Dict{String,Any}())
# 新: case = Case(param_dim, param, opt, mesh, index)
```

## 4. Initialisation.jl 核心改动

### 4.1 `ModelInitialisation` — 不动

48 行代码完全保持原样。

### 4.2 `ModelInitialisation_MultiSPMe` — 布局构建简化

```julia
# 旧 (6 行):
empty!(case.multi_spme_layout)
case.multi_spme_layout["ne"] = ne
case.multi_spme_layout["n_chem"] = n_chem
case.multi_spme_layout["nT"] = nT
case.multi_spme_layout["n_total"] = length(y0)
case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)

# 新 (1 行):
case.layout = MultiSPMeLayout(ne, n_chem, nT)
```

函数名保持不变（`ModelInitialisation_MultiSPMe`），内部 Dict 访问全部替换为 `case.layout`。

### 4.3 `MultiSPMe_extract_element_state` → `extract_element_state`

```julia
# 旧:
function MultiSPMe_extract_element_state(y::Array{Float64}, e::Int, case::Case)
    y_vec = vec(y)
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    offset = (e - 1) * n_chem
    yt_e = y_vec[(offset + 1):(offset + n_chem)]
    return yt_e
end

# 新:
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end
```

### 4.4 `MultiSPMe_get_thermal_dofs` → `get_thermal_dofs`

```julia
# 旧:
function MultiSPMe_get_thermal_dofs(y::Array{Float64}, case::Case)
    y_vec = vec(y)
    layout = case.multi_spme_layout
    thermal_range = layout["thermal_range"]
    return y_vec[thermal_range]
end

# 新:
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end
```

### 4.5 `MultiSPMe_update_state` → `update_state`

```julia
# 旧: 35 行，含手动 Dict 检查
# 新:
function update_state(y::AbstractVector, layout::MultiSPMeLayout;
                      element_index::Union{Nothing,Int}=nothing,
                      element_state::Union{Nothing,Vector{Float64}}=nothing,
                      thermal_nodes::Union{Nothing,Vector{Float64}}=nothing)
    y_new = copy(y)
    if element_index !== nothing
        offset = (element_index - 1) * layout.n_chem
        y_new[(offset + 1):(offset + layout.n_chem)] .= element_state
    end
    if thermal_nodes !== nothing
        y_new[layout.thermal_range] .= thermal_nodes
    end
    return y_new
end
```

## 5. 连锁修改文件

### 5.1 `Solve.jl`

| 改动 | 旧 | 新 |
|------|-----|-----|
| 布局构建 | 6 行 Dict 赋值 | `case.layout = MultiSPMeLayout(ne, n_chem, nT)` |
| 空检测 | `isempty(case.multi_spme_layout)` | `case.layout === nothing` |
| 函数调用 | `MultiSPMe_extract_element_state(...)` | `extract_element_state(y, e, case.layout)` |
| 函数调用 | `MultiSPMe_get_thermal_dofs(...)` | `get_thermal_dofs(y, case.layout)` |
| 函数调用 | `MultiSPMe_update_state(...)` | `update_state(y, case.layout; ...)` |

### 5.2 `CycleSolver.jl`

- 删除 `_ensure_multi_spme_layout!` 函数（28 行）
- 所有 `case.multi_spme_layout[...]` 替换为 `case.layout.xxx`
- 空检测改为 `case.layout === nothing` 时 `case.layout = MultiSPMeLayout(ne, n_chem, nT)`

### 5.3 `Jellyrollmodel.jl`

```julia
# 旧: 6 行 Dict 赋值
case_new.multi_spme_layout["interface_pairs"] = interface_pairs
# ...

# 新:
case_new.geometry = MeshGeometry(
    mesh_data.element_layer,
    mesh_data.is_inner_layer,
    mesh_data.layer_weights,
    interface_pairs,
    mesh_data.czm_element_map,
    mesh_data.inner_nodes,
    mesh_data.outer_nodes
)
```

### 5.4 `Variables.jl`

- `case.multi_spme_layout["ne"]` → `case.layout.ne`
- `case.multi_spme_layout["n_total"]` → `case.layout.n_total`
- `case.multi_spme_layout["thermal_range"]` → `case.layout.thermal_range`

### 5.5 `Parallelsolution.jl`

- `case.multi_spme_layout[...]` → `case.layout.xxx` 或 `case.geometry.xxx`
- 函数调用改名同步

### 5.6 `ThermalDistributed.jl`

- `case.multi_spme_layout[...]` → `case.layout.xxx` 或 `case.geometry.xxx`

### 5.7 `SPMe.jl`

- `MultiSPMe_extract_element_state(...)` → `extract_element_state(y, e, case.layout)`

### 5.8 `CzmSolve.jl` / `mechanical.jl`

- `case.multi_spme_layout["czm_mesh"]` → `case.czm_mesh`

### 5.9 `JuBat.jl` (模块入口)

- 添加 `include("CouplingState.jl")` 到 SetCase.jl 之前
- 导出新函数: `extract_element_state`, `get_thermal_dofs`, `update_state`

## 6. 不动的部分

- `ModelInitialisation` 函数（48 行）完全不动
- `ModelInitialisation_MultiSPMe` 中获取电化学初始状态的核心逻辑不变
- 所有物理/数学计算逻辑不变
- P2D 模型相关代码不动
- `ThermalLumped` 不动

## 7. 预期效果

| 指标 | 重构前 | 重构后 |
|------|--------|--------|
| Dict 访问点 | ~20 处 | 0 处 |
| 布局构建重复 | 4 处 x 6 行 | 4 处 x 1 行 |
| `_ensure_multi_spme_layout!` | 28 行 | 删除 |
| 类型安全 | 无（运行时才发现拼写错误） | 编译期检查 |
| 函数命名 | `MultiSPMe_*` 前缀 | snake_case 无前缀 |

## 8. 迁移映射表

| 旧访问 | 新访问 | 类型 |
|--------|--------|------|
| `case.multi_spme_layout["ne"]` | `case.layout.ne` | `Int` |
| `case.multi_spme_layout["n_chem"]` | `case.layout.n_chem` | `Int` |
| `case.multi_spme_layout["nT"]` | `case.layout.nT` | `Int` |
| `case.multi_spme_layout["n_total"]` | `case.layout.n_total` | `Int` |
| `case.multi_spme_layout["chem_range"]` | `case.layout.chem_range` | `UnitRange{Int}` |
| `case.multi_spme_layout["thermal_range"]` | `case.layout.thermal_range` | `UnitRange{Int}` |
| `case.multi_spme_layout["element_layer"]` | `case.geometry.element_layer` | `Vector{Int}` |
| `case.multi_spme_layout["interface_pairs"]` | `case.geometry.interface_pairs` | `Vector{Tuple{Int,Int}}` |
| `case.multi_spme_layout["czm_element_map"]` | `case.geometry.czm_element_map` | `Dict{Int,Int}` |
| `case.multi_spme_layout["czm_mesh"]` | `case.czm_mesh` | `CohesiveMesh` |
| `case.multi_spme_layout["layer_weights"]` | `case.geometry.layer_weights` | `Matrix{Float64}` |

## 9. 测试验证

- 使用 `example/minimal_example.jl` 验证基础 SPM 路径不受影响
- 使用 `example/SPMe_Thermal_example.jl` 验证 SPMe + 热耦合路径
- 使用 `example/testexample.jl` 验证全耦合路径
- 确认所有 `case.layout` 和 `case.geometry` 访问不会在运行时抛出 `nothing` 错误
