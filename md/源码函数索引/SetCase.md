# SetCase.jl

- **源文件**: `src/SetCase.jl`
- **行数**: 116 行
- **函数/struct 计数**: 1 个 struct（`Case`）；2 个函数（`SetCase`、5 参数兼容构造器 `Case`）
- **职责**: 装配仿真案例对象 `case` —— 调用 `NormaliseParam` 完成参数归一化、按 `opt.model` 构建网格、建立状态向量索引
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/02_几何与网格.md`

## 数据结构

### `Case`（L100-L112，`mutable struct`）

仿真案例顶层容器，贯穿整个仿真生命周期。

| 字段 | 类型 | 含义 |
|------|------|------|
| `param_dim` | Params | 有量纲参数 |
| `param` | Params | 无量纲参数（由 `NormaliseParam` 产出） |
| `opt` | Option | 求解器配置 |
| `mesh` | Dict{String, Mesh} | 离散化网格（键如 "negative particle"、"electrolyte" 等） |
| `index` | Dict{String, Union{Array{Int64}, Int64}} | 状态向量中未知量的索引映射 |
| `layout` | Union{Nothing, MultiSPMeLayout} | 布局索引（初始化后不变；`CouplingState.jl` 定义） |
| `geometry` | Union{Nothing, MeshGeometry} | 几何拓扑（构建后不变） |
| `czm_mesh` | Union{Nothing, CohesiveMesh} | CZM 网格（演化但类型明确） |
| `czm_cache` | Union{Nothing, CZMAssemblyCache} | CZM 装配缓存（E/ν 变化时重建） |
| `czm_layout` | Union{Nothing, CzmLayout} | CZM 布局 + u_prev（跨时间步） |
| `czm_param_cache` | Union{Nothing, CzmParamCache} | per-interface 参数缓存（v5 新增） |

后 6 个字段（layout 之后）默认 `nothing`，由后续初始化函数（如 `setup_thermal2D_mesh`、`create_czm_mesh`、`build_czm_cache`）填入。

## 函数清单

### `SetCase`（L1-L97）

```julia
function SetCase(param_dim::Params, opt::Option, y0::Array=[])
```

**职责**: 案例装配入口 —— 参数归一化 → 构建网格 → 建立状态向量索引 → 返回 `Case` 对象。

**关键逻辑**:
- L8：调用 `NormaliseParam(param_dim)` 得到无量纲 `param`
- L10-L11：粘性正则化参数归一化 `tau_visc* = tau_visc / t0`（写入 `param.cohesive`）
- L13-L43：`opt.model == "SPM" || "SPMe"` 分支 —— 构建颗粒网格（`mesh_np`、`mesh_pp`）+ electrolyte 网格（仅 SPMe），通过 `SetMesh` 与 `PickElement` 切分 NE/SP/PE 子段
- L44-L78：`opt.model == "P2D"` 分支 —— 颗粒网格 × 电极单元数（`MultipleMesh`）形成 P2D 多颗粒网格
- L79-L85：热模型索引（lumped 单 DOF；distributed2D 占位第一个 thermal DOF）
- L86-L94：P2D 电位索引（电极电位 + electrolyte 电位）
- L95：调用 5 参数兼容构造器 `Case(param_dim, param, opt, mesh, index)`

**状态向量索引键（SPM/SPMe）**:
- `"negative particle lithium concentration"` / `"positive particle lithium concentration"`
- `"...surface lithium concentration"`
- SPMe 额外：`"electrolyte lithium concentration"` 及 NE/SP/PE 子段
- 热模型：`"temperature"`
- P2D 额外：`"negative/positive electrode potential"`、`"electrolyte potential..."`

**跨文件依赖**: `NormaliseParam`、`Params`（`SetParams.jl`）、`SetMesh`、`PickElement`、`MultipleMesh`、`Mesh`（`SetMesh.jl`）、`Option`（`Option.jl`）、`MultiSPMeLayout`/`MeshGeometry`/`CohesiveMesh`/`CZMAssemblyCache`/`CzmLayout`/`CzmParamCache`（`CouplingState.jl`、`czm.jl` 等）

### `Case`（5 参数兼容构造器，L115-L116）

```julia
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, nothing, nothing, nothing)
end
```

**职责**: 提供 5 参数便捷构造 —— 将 layout/geometry/czm_* 6 个可选字段统一填 `nothing`，委托给 11 参数主构造器。

**跨文件依赖**: `Case` 主 struct（同文件 L100）

## 省略项

- L2-L6：函数内 docstring 字符串
- L96：`return case`

### [DEBUG]

无

### [PLACEHOLDER]

无

> 说明：L51 / L29 中电解质网格 `space` 数组为字面几何坐标（非兜底占位）；L80-L84 distributed2D 分支仅用首 DOF 作代表温度，注释 L83 「Use the first thermal DOF as a representative temperature in the state vector」属设计说明而非兜底（实际 thermal2D 模型 DOF 由后续 `setup_thermal2D_mesh` 扩展），不标 PLACEHOLDER。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L13-L94 | `if opt.model == "SPM" \|\| opt.model == "SPMe" ... elseif opt.model == "P2D" ... end` 嵌套 3 层 + 多个 `if opt.model == "SPMe"` / `if opt.thermalmodel == "lumped" elseif "distributed2D"` / `if opt.model == "P2D"` 子分支 | 嵌套层级较深，且 SPM/SPMe 与 P2D 分支共享大量重复代码（mesh_np/mesh_pp/mesh_el 构造）。建议抽取 `build_mesh_spm(opt, param)` 与 `build_mesh_p2d(opt, param)` 子函数；状态索引构建可参数化 |
| L82-L85 | `elseif opt.thermalmodel == "distributed2D" # Use the first thermal DOF ... index["temperature"] = [v0 + 1]` | distributed2D 分支不递增 `v0`（与 lumped 不同），后续 thermal2D 装配会另算 DOF 数。控制流不直观，建议加注释说明 thermal DOF 由 thermal2D 模块接管 |
