# Tools.jl

- **源文件**: `src/Tools.jl`
- **行数**: 192 行
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

通用 Q4 单元高斯积分驱动。回调签名 `f(ξ, η, w, dNdx, dNdy, detJ)`。`abs(detJ) < 1e-15` 跳过该积分点（退化保护，L38）。返回 `nothing`——所有结果通过回调副作用累积。

### `identify_boundary_nodes(mesh, param, opt=nothing)` — L47-L70

识别内/外边界节点。返回 `(is_inner, is_outer)::Tuple{Vector{Bool}, Vector{Bool}}`。

- 通过 `typeof(mesh)` 分发 `nnode`（Mesh → `mesh.nlen`，CohesiveMesh → `mesh.nnode`）（L48-L54）
- 计算 Jellyroll 螺旋几何的内/外 θ 区间（L61-L64）；`bval = max(b, 1e-12)` 防止除零（L59）
- 委托给 `edge_boundary`（L67-L68），对每个节点单独判定
- `opt` 参数当前未使用（L47 签名保留）

### `compute_separation(elem, node::Matrix{Float64}, u::Vector{Float64})` — L81-L133

计算单个 cohesive 单元的法/切向平均分离量 `(δ_n, δ_t)`。

- 从 `nodes_bottom`/`nodes_top` 取 4 节点，构造局部坐标系 `R = [n_vec t_vec]`（L94-L96）
- `L < 1e-15` 时返回 `(0.0, 0.0)`（退化保护，L90-L92）
- 2 点 Gauss 积分（NCweight）累加 `B_local * u_e`（L98-L125）
- `total_w == 0` 时返回 0（防御，L127-L130）

### `element_nodal_mean(mesh::Mesh, nodal_values)` — L140-L148

节点场→单元均值的快速版（算术平均，非加权）。`@inbounds` 优化。对每个 Q4 单元的 4 个节点值求和后除以节点数。

### `thermal2D_volume_average_temperature(mesh::Mesh, T_nodes)` — L150-L166

热网格体积加权平均温度。双重循环（积分点 × 单元局部节点）累加 `wi * T_nodes[nodes[i]]`（L160-L162），分母为 `den`。`den == 0` 时回退到 `mean(T_nodes)`（L165）。

### `q4_center_gradients(node::AbstractMatrix{<:Real}, elem_nodes)` — L178-L192

返回 Q4 形函数在单元中心 `(ξ=0, η=0)` 的空间梯度 `(dNdx, dNdy, detJ)`。`abs(detJ) < 1e-12` 时返回 `nothing`（退化保护，L185-L187）。返回的 `dNdx`/`dNdy` 为 `@view`，不应跨函数持有。

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

无。文件无 `println`/`print`/`@show`/非结构化 `@info`。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L47 | `function identify_boundary_nodes(mesh, param, opt=nothing)` 中 `opt` 参数未被函数体引用 | 未使用的可选参数——签名暴露但行为忽略；调用方若依赖 `opt` 控制行为会得到无差别结果。建议移除或文档化 |
| L59 | `bval = max(b, 1e-12)` | magic number 1e-12 + 兜底防除零；正常几何下 `b = t_repeat/(2π)` 远大于 1e-12，仅在退化参数（`t_repeat=0`）下生效 |
| L65 | `tol = 1e-4` | magic number 容差；硬编码，未与网格尺度关联 |
| L90 | `if L < 1e-15; return 0.0, 0.0` | 退化单元长度阈值，magic number 1e-15；正常网格下不应触发 |
| L127 | `if total_w > 0.0` 后 else 隐式返回 0 | 防御性兜底——`total_w` 为 Gauss 权重和，正常情况下恒为 2.0；仅在 `NCweight` 返回空数组时触发 |
| L165 | `den > 0.0 ? num / den : mean(T_nodes)` | 分母为 0 时回退到算术平均；退化网格（无积分点）下的兜底 |
| L185 | `if abs(detJ) < 1e-12; return nothing` | 退化单元判定，magic number 1e-12；调用方需处理 `nothing` 返回 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L48 | `if mesh isa Mesh; nnode = mesh.nlen; elseif mesh isa CohesiveMesh; nnode = mesh.nnode; else; error(...)` | 多类型分发 + error——可改为在 `Mesh`/`CohesiveMesh` 上定义 `nnodes(mesh)` 方法，调用点 `nnode = nnodes(mesh)` |
| L116 | `u_e[1] = u[2 * n1 - 1]; u_e[2] = u[2 * n1]; u_e[3] = u[2 * n2 - 1]; u_e[4] = u[2 * n2]; ...`（连续 8 行重复模板，跨 L116-L119） | 抽出 `gather_disp!(u_e, u, ns) = for (k, n) in enumerate(ns); u_e[2k-1] = u[2n-1]; u_e[2k] = u[2n]; end` helper |
