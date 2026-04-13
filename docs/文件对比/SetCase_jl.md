# SetCase.jl

## 文件状态
修改

## main分支
- **行数**: 106
- **主要函数/结构体列表**:
  | 行号 | 名称 | 说明 |
  |------|------|------|
  | 1 | `SetCase(param_dim::Params, opt::Option, y0::Array=[])` | 创建 Case 对象，设置网格和状态索引 |
  | 101 | `mutable struct Case` | 核心数据结构：param_dim, param, opt, mesh, index |

## Parameters_Design分支
- **行数**: 110
- **主要函数/结构体列表**:
  | 行号 | 名称 | 说明 |
  |------|------|------|
  | 1 | `SetCase(param_dim::Params, opt::Option, y0::Array=[])` | 创建 Case 对象（已修改，新增 distributed2D 索引） |
  | 97 | `mutable struct Case` | 核心数据结构（已修改，新增 layout/geometry/czm_mesh 类型化字段） |
  | 109 | `Case(param_dim, param, opt, mesh, index)` | 5参数兼容构造器 |

## 变更详情

### 新增函数
无

### 修改函数

#### `SetCase(param_dim::Params, opt::Option, y0::Array=[])`
- **变更 1**: 移除 `sP2D` 模型分支（`P2D || sP2D` 变为仅 `P2D`）

- **变更 2**: 新增 `distributed2D` 热模型索引：
  ```julia
  elseif opt.thermalmodel == "distributed2D"
      index["temperature"] = [v0 + 1]
  ```
  将第一个热 DOF 作为代表性温度索引。完整的温度自由度范围由 `MultiSPMeLayout` 管理。

- **变更 3**: Case 构造调用保持原始 5 参数（第 92 行）：
  - `Case(param_dim, param, opt, mesh, index)` — 通过 5 参数兼容构造器初始化
  - 新增字段自动初始化为 `nothing`

#### `mutable struct Case`
- **变更**: 替换 `multi_spme_layout::Dict{String,Any}` 为 3 个类型化字段（第 97-106 行）
- **原定义**:
  ```julia
  mutable struct Case
      param_dim::Params
      param::Params
      opt::Option
      mesh::Dict{String, Mesh}
      index::Dict{String, Union{Array{Int64}, Int64}}
  end
  ```
- **新定义**:
  ```julia
  mutable struct Case
      param_dim::Params                      # dimensional parameters
      param::Params                          # dimensionless parameters
      opt::Option                            # solver options
      mesh::Dict{String, Mesh}               # discretisation meshes
      index::Dict{String, Union{Array{Int64}, Int64}} # indices of unknowns
      layout::Union{Nothing, MultiSPMeLayout}   # 布局索引（初始化后不变）
      geometry::Union{Nothing, MeshGeometry}   # 几何拓扑（构建后不变）
      czm_mesh::Union{Nothing, CohesiveMesh}   # CZM 网格（演化但类型明确）
  end
  ```
- **新增兼容构造器**:
  ```julia
  function Case(param_dim, param, opt, mesh, index)
      Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing)
  end
  ```
- **设计原则**: Fail-fast — 直接访问字段，未初始化就让 Julia 抛异常（禁止 haskey / get(dict, default)）

### 删除函数
无

## 依赖关系

### 该文件依赖哪些其他文件
- `Parameters.jl` / 参数模块 — 使用 `Params`, `NormaliseParam`
- `Mesh.jl` — 使用 `Mesh`, `SetMesh`
- `Option.jl` — 使用 `Option` 结构体

### 哪些文件依赖该文件
- 几乎所有核心文件都依赖 `Case` 结构体：
  - `Solve.jl` — `Solve(case::Case)`, `CallModel(case::Case, ...)`
  - `Initialisation.jl` — `ModelInitialisation(case::Case)`, `ModelInitialisation_MultiSPMe(case::Case)`
  - `SPMe.jl` — `SPMe(case::Case, ...)`, `SPMe_element(case::Case, ...)`
  - `Parallelsolution.jl` — `solve_branch_currents_newton(case::Case, ...)`
  - `ThermalDistributed.jl` — 所有热学函数
  - `CycleSolver.jl` / `CycleData.jl` — 循环求解
  - `JuBat.jl` — export `SetCase`, `Case`

### 新增的外部依赖
无

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系
SetCase.jl 是耦合架构的**数据结构基础**：
- **multi-SPMe**: `Case` 结构体新增 `layout::Union{Nothing, MultiSPMeLayout}` 字段，提供类型安全的状态向量布局管理
- **distributed2D**: 新增 `distributed2D` 热模型的索引设置，使热 DOF 能被纳入标准状态向量
- **CZM**: `Case` 结构体新增 `czm_mesh::Union{Nothing, CohesiveMesh}` 字段，类型化 CZM 网格引用
- **几何**: `Case` 结构体新增 `geometry::Union{Nothing, MeshGeometry}` 字段，存储网格几何拓扑信息

### 哪些变更是耦合相关的
- `Case` 结构体新增 `layout`, `geometry`, `czm_mesh` 三个类型化字段 — 替代原 Dict，提供 fail-fast 语义
- `SetCase` 新增 `distributed2D` 索引 — distributed2D 热模型的状态向量支持
- 5 参数兼容构造器 — 确保 main 分支代码无需修改

### 哪些变更是独立的
- 移除 `sP2D` 模型分支 — 模型简化，与耦合无关
- `Case` 结构体字段注释改善 — 文档改进

## 后续变更 (2026-04-01)

- 移除了 `SetCase` 函数开头的 `if opt.model == "thermal"` 早期返回块（原第 7-12 行，共约 7 行）
- 该块是一个纯热模式快捷路径，绕过了正常的网格构造和参数归一化流程
- 移除后，所有模型类型（包括 `"thermal"`）均走正常的参数归一化和 Case 构建路径
- 文件行数从约 111 行减少到约 104 行

## 后续变更 (2026-04-07)

- **结构体重构**: `multi_spme_layout::Dict{String,Any}` → 3 个类型化字段 (`layout`, `geometry`, `czm_mesh`)
- **新增 5 参数兼容构造器**: `Case(param_dim, param, opt, mesh, index)` 自动将新字段初始化为 `nothing`
- **Fail-fast 设计**: 移除 `haskey` / `get(dict, default)` 模式，直接字段访问
- 行数从约 104 行增加到 110 行
