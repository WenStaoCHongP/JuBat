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
- **行数**: 111
- **主要函数/结构体列表**:
  | 行号 | 名称 | 说明 |
  |------|------|------|
  | 1 | `SetCase(param_dim::Params, opt::Option, y0::Array=[])` | 创建 Case 对象（已修改，新增 thermal 模式和 distributed2D 索引） |
  | 105 | `mutable struct Case` | 核心数据结构（已修改，新增 multi_spme_layout 字段） |

## 变更详情

### 新增函数
无

### 修改函数

#### `SetCase(param_dim::Params, opt::Option, y0::Array=[])`
- **变更 1**: 新增纯热模型早期返回路径（第 7-12 行）：
  ```julia
  if opt.model == "thermal"
      param = NormaliseParam(param_dim)
      mesh = Dict{String, Mesh}()
      index = Dict{String, Union{Array{Int64}, Int64}}()
      return Case(param_dim, param, opt, mesh, index, Dict{String,Any}())
  end
  ```
  当模型为 `"thermal"` 时，跳过电化学网格和索引创建，直接返回最小 Case 对象。

- **变更 2**: 移除 `sP2D` 模型分支（第 45 行 `P2D || sP2D` 变为仅 `P2D`）

- **变更 3**: 新增 `distributed2D` 热模型索引（第 85-88 行）：
  ```julia
  elseif opt.thermalmodel == "distributed2D"
      index["temperature"] = [v0 + 1]
  ```
  将第一个热 DOF 作为代表性温度索引。完整的温度自由度范围由 `multi_spme_layout` 管理。

- **变更 4**: 移除 `sP2D` 引用管理逻辑（原第 96-101 行）：
  - 删除 `opt.cite = vcat(opt.cite, "ai2024b")`（sP2D 引用）
  - 删除 `opt.cite = vcat(opt.cite, "ai2023")`（L3 网格引用）

- **变更 5**: Case 构造调用新增 `Dict{String,Any}()` 参数（第 98 行）：
  - **原**: `Case(param_dim, param, opt, mesh, index)`
  - **新**: `Case(param_dim, param, opt, mesh, index, Dict{String,Any}())`
  - 初始化空的 `multi_spme_layout` 字典

#### `mutable struct Case`
- **变更**: 新增字段 `multi_spme_layout::Dict{String,Any}`（第 111 行）
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
      multi_spme_layout::Dict{String,Any}    # multi-SPMe layout (populated after init)
  end
  ```
- **字段用途**: 存储 multi-SPMe 模式的结构布局信息，包括：
  - `ne`: 热单元数
  - `n_chem`: 单个单元电化学 DOF 数
  - `nT`: 热节点数
  - `n_total`: 全局状态向量总长度
  - `chem_range`: 电化学 DOF 范围
  - `thermal_range`: 热 DOF 范围

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
- **multi-SPMe**: `Case` 结构体新增 `multi_spme_layout` 字段，这是 multi-SPMe 模式全局状态管理的核心数据结构，被 `Initialisation.jl`（填充）、`Solve.jl`（读取）、`Parallelsolution.jl`（间接通过 Case 访问参数）共同使用
- **distributed2D**: 新增 `distributed2D` 热模型的索引设置，使热 DOF 能被纳入标准状态向量
- **纯热模式**: `"thermal"` 模型的早期返回路径，允许仅进行热学仿真（不创建电化学网格）

### 哪些变更是耦合相关的
- `Case` 结构体新增 `multi_spme_layout` 字段 — 所有 multi-SPMe 功能的基础
- `SetCase` 新增 `distributed2D` 索引 — distributed2D 热模型的状态向量支持
- `SetCase` 新增 `"thermal"` 模型早期返回 — 纯热学仿真路径
- Case 构造时传入空 `Dict{String,Any}()` — 确保 `multi_spme_layout` 字段被初始化

### 哪些变更是独立的
- 移除 `sP2D` 模型分支 — 模型简化，与耦合无关
- 移除引用管理逻辑（`opt.cite` 相关代码）— 代码清理
- `Case` 结构体字段注释改善 — 文档改进

## 后续变更 (2026-04-01)

- 移除了 `SetCase` 函数开头的 `if opt.model == "thermal"` 早期返回块（原第 7-12 行，共约 7 行）
- 该块是一个纯热模式快捷路径，绕过了正常的网格构造和参数归一化流程
- 移除后，所有模型类型（包括 `"thermal"`）均走正常的参数归一化和 Case 构建路径
- 文件行数从约 111 行减少到约 104 行
