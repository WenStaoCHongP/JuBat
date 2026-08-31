# Tools.jl

- **源文件**: `src/Tools.jl`
- **行数**: 167 行
- **函数/struct 计数**: 8 个（0 struct + 8 函数）
- **职责**: 通用数值与几何工具——Arrhenius 温度补偿、网格积分驱动、Q4 单元高斯积分回调接口、内/外边界节点识别、cohesive 单元法/切分离量计算、节点场→单元均值、热网格体积加权平均温度、Q4 中心梯度
- **相关技术文档**: `md/02_几何与网格.md`（边界节点识别）、`md/05_热模型_二维分布式.md`（体积平均温度）、`md/06_内聚力模型_CZM.md`（separation 计算）

## 数据结构

无独立 struct。函数均操作外部定义的类型（`Mesh`、`CohesiveMesh`、`CohesiveElement`）。

## 函数清单

### `Arrhenius(Eac::Float64, T::Union{Float64, Array{Float64}})` — L1-L9

Arrhenius 温度补偿 `exp.(Eac * (1 .- 1 ./ T))`。`Eac` 为已归一化（`Eac/(R·T0)`）活化能，`T` 为归一化温度（`T/T0`）。支持标量与数组输入。

### `IntV(x::Array{Float64}, mesh::Mesh)` — L11-L19

对域 `mesh` 积分 `sum(x .* mesh.gs.weight .* mesh.gs.detJ)`。`x` 为单元/积分点系数向量。

### `IntQ4(f::Function, x_e, y_e; order=2)` — L29-L45

通用 Q4 单元高斯积分驱动。回调签名 `f(ξ, η, w, dNdx, dNdy, detJ)`。`detJ ≤ 1e-15` 时报错（L38），不静默跳过积分点。返回 `nothing`——所有结果通过回调副作用累积。

### `identify_boundary_nodes(mesh, param, opt=nothing)` — L47-L69

识别内/外边界节点。返回 `(is_inner, is_outer)::Tuple{Vector{Bool}, Vector{Bool}}`。

- 通过 `typeof(mesh)` 分发 `nnode`（Mesh → `mesh.nlen`，CohesiveMesh → `mesh.nnode`）（L48-L54）
- 计算 Jellyroll 螺旋几何的内/外 θ 区间（L60-L63）
- 委托给 `edge_boundary`（L66-L67），对每个节点单独判定
- `opt` 参数当前未使用（L47 签名保留）

### `compute_separation(czm_mesh::CohesiveMesh, elem::AbstractCohesiveElement, u::Vector{Float64})` — L77-L110

Jellyroll 拓扑有向入口，计算单个 cohesive 单元的平均法/切向分离 `(δ_n, δ_t)`。

- 标架 `R` 由 `cohesive_local_frame` 构造：法向按 host-inner → host-outer 质心方向拓扑定向（L82），免疫节点排序/绕向约定
- 2 点 Gauss 积分（NCweight(2)）累加 `R * B_global * u_e`（L85-L107）
- 无 Λ 换算：结果停留在位移空间（L 归一）；装配路径中 δ_czm 归一分离见 `assemble_czm_system`（`src/czm.jl`）
- 旧 `compute_separation(elem, node, u)` 已删除；调用方必须传入 `CohesiveMesh`，不得再用节点边序猜测法向

### `element_nodal_mean(mesh::Mesh, nodal_values)` — L117-L125

节点场→单元均值的快速版（算术平均，非加权）。`@inbounds` 优化。对每个 Q4 单元的 4 个节点值求和后除以节点数。

### `thermal2D_volume_average_temperature(mesh::Mesh, T_nodes)` — L127-L143

热网格体积加权平均温度。双重循环（积分点 × 单元局部节点）累加积分权重 × 节点温度（L132-L141）；`return num / den` 直接相除（L142），零分母交由原生除法语义暴露。

### `q4_center_gradients(node::AbstractMatrix{<:Real}, elem_nodes)` — L155-L167

返回 Q4 形函数在单元中心 `(ξ=0, η=0)` 的空间梯度 `(dNdx, dNdy, detJ)`。`detJ ≤ 1e-12` 时报错（L162），不返回 `nothing`。返回的 `dNdx`/`dNdy` 为 `@view`，不应跨函数持有。

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

无。文件无 `println`/`print`/`@show`/非结构化 `@info`。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L47 | `function identify_boundary_nodes(mesh, param, opt=nothing)` 中 `opt` 参数未被函数体引用 | 未使用的可选参数——签名暴露但行为忽略；调用方若依赖 `opt` 控制行为会得到无差别结果。建议移除或文档化 |
| L64 | `tol = 1e-4` | magic number 容差；硬编码，未与网格尺度关联 |
| L142 | `return num / den` | 零分母时原生产生 NaN（无静默回退）；退化网格下结果非有限值，由上层有限性检查拦截 |
| L162 | `detJ > 1e-12 \|\| error(...)` | 退化单元判定阈值，magic number 1e-12；正常网格下不应触发，触发时报错而非返回 `nothing` |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L48 | `if mesh isa Mesh; nnode = mesh.nlen; elseif mesh isa CohesiveMesh; nnode = mesh.nnode; else; error(...)` | 多类型分发 + error——可改为在 `Mesh`/`CohesiveMesh` 上定义 `nnodes(mesh)` 方法，调用点 `nnode = nnodes(mesh)` |
