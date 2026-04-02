# Initialisation.jl 优化设计规格

> 日期: 2026-04-02
> 状态: 待实施（已通过审查修订）
> 基于文档: 03_结构体重构方案.md, 05_Initialisation_jl优化方案.md
> 审查修订: v2 — 修复 11 处审查问题

---

## 1. 目标

对 Initialisation.jl 及其关联文件执行类型安全重构，消除 `multi_spme_layout::Dict{String,Any}` 中**布局索引**和**网格几何**部分的类型不安全问题，统一命名规范，简化布局构建逻辑。

**范围边界**：运行时状态键（`thermal_variables`, `thermal_update_fn`, `thermal_record`, `polar_mesh_data`）归入 `thermal_extras::Dict{String,Any}` 暂存，待后续 CouplingState.jl 完整重构时迁移到 `SimulationState`。

**设计约束（Fail-Fast）**：
- 禁止使用 `haskey`、`get(dict, key, default)` 等回退兜底模式
- 禁止 `try/catch` 吞没错误
- 直接访问字段/索引，如果状态未初始化就让错误立即暴露（Julia 的 `nothing` 访问和越界会自动抛出明确异常）
- 功能开关（如 `case.opt.czm_enabled`）仅控制是否进入某路径，不负责掩盖缺失数据

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
    element_layer::Vector{Int}                      # 每个单元的层类型
    is_inner_layer::Vector{Bool}                    # 是否为内层（与 JellyrollMesh 类型一致）
    interface_pairs::Vector{Tuple{Int,Int}}         # CZM 界面配对
    czm_element_map::Dict{Int,Vector{Int}}          # 热单元号 → CZM 单元索引向量（一对多映射）
    inner_nodes::Vector{Int}                        # 内边界节点索引
    outer_nodes::Vector{Int}                        # 外边界节点索引
end
```

**注意**：`layer_weights` 不纳入 `MeshGeometry`。`layer_weights` 由 `jellyroll_element_properties()` 函数单独计算（Jellyrollmodel.jl:245-316），当前通过 `variables["thermal2D layer_weights"]` 存储和访问（Solve.jl:172-178），该机制保持不变。

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
    # ---- 原字段 ----
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}

    # ---- 替代 multi_spme_layout ----
    layout::Union{Nothing, MultiSPMeLayout}         # 布局索引（初始化后不变）
    geometry::Union{Nothing, MeshGeometry}           # 几何拓扑（构建后不变）
    czm_mesh::Union{Nothing, CohesiveMesh}           # CZM 网格（演化但类型明确）

    # ---- 运行时状态暂存（待后续 SimulationState 重构） ----
    thermal_extras::Dict{String,Any}                 # thermal_variables, thermal_update_fn, thermal_record, polar_mesh_data
end
```

### 3.2 兼容构造器

```julia
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, Dict{String,Any}())
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

函数名保持 `ModelInitialisation_MultiSPMe` 不变（与 05 文档的 snake_case 改名建议有偏差，理由：避免同时改名和改类型造成过大 diff，降低出错风险）。

### 4.3 `MultiSPMe_extract_element_state` → `extract_element_state`

```julia
# 新:
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end
```

### 4.4 `MultiSPMe_get_thermal_dofs` → `get_thermal_dofs`

```julia
# 新:
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end
```

### 4.5 `MultiSPMe_update_state` → `update_state`

```julia
# 新:
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

保留了关键的范围/长度验证（使用 `@assert`），防止逻辑错误导致越界写入。

## 5. 连锁修改文件

### 5.1 `Solve.jl`

| 改动 | 旧 | 新 |
|------|-----|-----|
| 布局构建 (行 104-109) | 6 行 Dict 赋值 | `case.layout = MultiSPMeLayout(ne, n_chem, nT)` |
| 空检测 | `isempty(case.multi_spme_layout)` | `case.layout === nothing` |
| 函数调用 (行 448, 452) | `MultiSPMe_extract_element_state(y, e, case)` | `extract_element_state(y, e, case.layout)` |
| 函数调用 | `MultiSPMe_get_thermal_dofs(y, case)` | `get_thermal_dofs(y, case.layout)` |
| 函数调用 | `MultiSPMe_update_state(y, case; ...)` | `update_state(y, case.layout; ...)` |
| 运行时键 (行 6-8) | `case.multi_spme_layout["thermal_variables"]` 等 | `case.thermal_extras["thermal_variables"]` 等 |
| 运行时键 (行 30, 45) | `case.multi_spme_layout["polar_mesh_data"]` | `case.thermal_extras["polar_mesh_data"]` |

### 5.2 `CycleSolver.jl`

- 删除 `_ensure_multi_spme_layout!` 函数（28 行）
- 所有 `case.multi_spme_layout[...]` 替换为 `case.layout.xxx`
- 空检测改为 `case.layout === nothing` 时 `case.layout = MultiSPMeLayout(ne, n_chem, nT)`

### 5.3 `Jellyrollmodel.jl`

```julia
# 旧 (行 546-551): 6 行 Dict 赋值
case_new.multi_spme_layout["interface_pairs"] = interface_pairs
case_new.multi_spme_layout["element_layer"] = mesh_data.element_layer
case_new.multi_spme_layout["is_inner_layer"] = mesh_data.is_inner_layer
case_new.multi_spme_layout["czm_element_map"] = mesh_data.czm_element_map
case_new.multi_spme_layout["inner_nodes"] = mesh_data.inner_nodes
case_new.multi_spme_layout["outer_nodes"] = mesh_data.outer_nodes

# 新:
case_new.geometry = MeshGeometry(
    mesh_data.element_layer,
    mesh_data.is_inner_layer,
    interface_pairs,
    mesh_data.czm_element_map,
    mesh_data.inner_nodes,
    mesh_data.outer_nodes
)
```

注意：`layer_weights` 不在此处设置，保持现有 `variables["thermal2D layer_weights"]` 机制。

### 5.4 `CycleData.jl`

| 行号 | 旧 | 新 |
|------|-----|-----|
| 65 | `_ensure_multi_spme_layout!(case)` | `if case.layout === nothing; case.layout = MultiSPMeLayout(ne, n_chem, nT); end` |
| 73 | `case.multi_spme_layout["thermal_range"]` | `case.layout.thermal_range` |
| 77 | `case.multi_spme_layout["thermal_range"]` | `case.layout.thermal_range` |
| 149 | `MultiSPMe_get_thermal_dofs(y_new, case)` | `get_thermal_dofs(y_new, case.layout)` |

### 5.5 `ThermalDistributed.jl`

```julia
# 旧 (行 192-193):
if case.opt.czm_enabled && haskey(case.multi_spme_layout, "czm_mesh")
    czm_mesh = case.multi_spme_layout["czm_mesh"]

# 新（直接访问，czm_enabled 保证 czm_mesh 已初始化）:
if case.opt.czm_enabled
    czm_mesh = case.czm_mesh
```

### 5.6 `CzmSolve.jl` / `mechanical.jl`

- `case.multi_spme_layout["czm_mesh"]` → `case.czm_mesh`
- 所有 `haskey(case.multi_spme_layout, "czm_mesh")` 检查删除，改为功能开关 `case.opt.czm_enabled` 直接控制路径

### 5.7 `JuBat.jl` (模块入口)

- 添加 `include("CouplingState.jl")` 到 SetCase.jl 之前
- 导出新函数: `extract_element_state`, `get_thermal_dofs`, `update_state`

### 5.8 无需改动的文件（已确认无 `multi_spme_layout` 访问）

- `Variables.jl` — 无访问
- `Parallelsolution.jl` — 无访问
- `SPMe.jl` — 无 `MultiSPMe_extract_element_state` 调用（实际调用在 Solve.jl 中）

## 6. 不动的部分

- `ModelInitialisation` 函数（48 行）完全不动
- `ModelInitialisation_MultiSPMe` 中获取电化学初始状态的核心逻辑不变
- 所有物理/数学计算逻辑不变
- P2D 模型相关代码不动
- `ThermalLumped` 不动
- `layer_weights` 的计算和存储机制不变（仍通过 `jellyroll_element_properties()` 和 `variables` Dict）

## 7. 预期效果

| 指标 | 重构前 | 重构后 |
|------|--------|--------|
| `multi_spme_layout` Dict 访问 | ~20 处 | 0 处 |
| 布局构建重复 | 4 处 x 6 行 | 4 处 x 1 行 |
| `_ensure_multi_spme_layout!` | 28 行 | 删除 |
| 类型安全 | 无（运行时才发现拼写错误） | 编译期检查（布局/几何字段） |
| 辅助函数命名 | `MultiSPMe_*` 前缀 | snake_case 无前缀 |
| `thermal_extras` Dict | N/A | 4 个运行时键暂存（待后续迁移） |

## 8. 迁移映射表

### 8.1 布局索引 → `case.layout`

| 旧访问 | 新访问 | 类型 |
|--------|--------|------|
| `case.multi_spme_layout["ne"]` | `case.layout.ne` | `Int` |
| `case.multi_spme_layout["n_chem"]` | `case.layout.n_chem` | `Int` |
| `case.multi_spme_layout["nT"]` | `case.layout.nT` | `Int` |
| `case.multi_spme_layout["n_total"]` | `case.layout.n_total` | `Int` |
| `case.multi_spme_layout["chem_range"]` | `case.layout.chem_range` | `UnitRange{Int}` |
| `case.multi_spme_layout["thermal_range"]` | `case.layout.thermal_range` | `UnitRange{Int}` |

### 8.2 网格几何 → `case.geometry`

| 旧访问 | 新访问 | 类型 |
|--------|--------|------|
| `case.multi_spme_layout["element_layer"]` | `case.geometry.element_layer` | `Vector{Int}` |
| `case.multi_spme_layout["is_inner_layer"]` | `case.geometry.is_inner_layer` | `Vector{Bool}` |
| `case.multi_spme_layout["interface_pairs"]` | `case.geometry.interface_pairs` | `Vector{Tuple{Int,Int}}` |
| `case.multi_spme_layout["czm_element_map"]` | `case.geometry.czm_element_map` | `Dict{Int,Vector{Int}}` |
| `case.multi_spme_layout["inner_nodes"]` | `case.geometry.inner_nodes` | `Vector{Int}` |
| `case.multi_spme_layout["outer_nodes"]` | `case.geometry.outer_nodes` | `Vector{Int}` |

### 8.3 CZM 网格 → `case.czm_mesh`

| 旧访问 | 新访问 | 类型 |
|--------|--------|------|
| `case.multi_spme_layout["czm_mesh"]` | `case.czm_mesh` | `CohesiveMesh` |
| `haskey(case.multi_spme_layout, "czm_mesh")` | 删除检查，由 `case.opt.czm_enabled` 控制路径 | N/A |

### 8.4 运行时状态 → `case.thermal_extras`（暂存）

| 旧访问 | 新访问 | 说明 |
|--------|--------|------|
| `case.multi_spme_layout["thermal_variables"]` | `case.thermal_extras["thermal_variables"]` | 待迁入 SimulationState |
| `case.multi_spme_layout["thermal_update_fn"]` | `case.thermal_extras["thermal_update_fn"]` | 待迁入 SimulationState |
| `case.multi_spme_layout["thermal_record"]` | `case.thermal_extras["thermal_record"]` | 待迁入 SimulationState |
| `case.multi_spme_layout["polar_mesh_data"]` | `case.thermal_extras["polar_mesh_data"]` | 待迁入 SimulationState |

### 8.5 不迁移

| 键 | 原因 |
|----|------|
| `layer_weights` | 仍通过 `variables["thermal2D layer_weights"]` 管理，由 `jellyroll_element_properties()` 计算 |

## 9. 测试验证

### 9.1 回归测试

| 测试脚本 | 验证路径 | 通过标准 |
|----------|----------|----------|
| `example/minimal_example.jl` | 基础 SPM（不触发 multi-SPMe） | 运行无报错，输出与重构前一致 |
| `example/SPMe_Thermal_example.jl` | SPMe + 热耦合 | 电压/温度结果相对误差 < 1e-10 |
| `example/testexample.jl` | 全耦合 (SPMe+热+CZM) | 运行无报错，输出一致 |

### 9.2 具体验证项

- 重构前保存 `testexample.jl` 的输出结果作为基准
- 重构后对比所有数值结果在 1e-10 相对误差内
- 确认 `case.layout` 在 multi-SPMe 路径中非 `nothing`
- 确认 `case.czm_mesh` 在 CZM 路径中非 `nothing`
- 确认 `case.thermal_extras` 在 thermal 路径中包含所需键
