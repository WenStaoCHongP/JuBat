# SetMesh.jl

- **源文件**: `src/SetMesh.jl`
- **行数**: 768 行
- **函数/struct 计数**: 4 个 struct + 2 个 abstract type + 10 个独立函数
- **职责**: 有限元网格核心基础设施——`GaussPoint` / `Mesh` 数据结构、`CohesiveMesh` / `CzmSubmesh` CZM 网格容器、`SetMesh` / `Mesh1D` / `Mesh2D` 网格生成入口、`PickElement` / `CombineMesh` / `MultipleMesh` 网格操作、`GetGS` 高斯积分、`LagrangeBasis` / `ShapeFunction1D` / `ShapeFunction2D` 形函数、`GSweight` / `NCweight` 积分权重
- **相关技术文档**: `md/02_几何与网格.md`

## 数据结构

### `mutable struct GaussPoint` — L3-L12

高斯积分点缓存，存储物理/局部坐标、权重、Jacobi 行列式、单元映射、形函数及其导数。

- 字段：`x`（物理坐标 ngs×dim）、`xi`（局部坐标 ngs×dim）、`weight`、`detJ`、`ele`（高斯点 → 单元 id）、`Ni`（形函数 ngs×nnode_ele）、`dNidx`（形函数导数，Q4 前 4 列 dN/dx 后 4 列 dN/dy）、`order`

### `mutable struct Mesh` — L14-L21

通用网格容器，封装节点/单元/高斯点。

- 字段：`type`（"Q4" / "L2" / "L3" / "COH2D4" / "B8"）、`dimension`、`node`（nnode×dim）、`nlen`（节点数）、`element`（ne×nnode_ele）、`gs::GaussPoint`

### `abstract type AbstractCohesiveElement` — L25

CZM 内聚力单元抽象基类，定义于此（而非 `czm.jl`）避免 include 顺序问题（`CohesiveMesh` 引用此类型，而 SetMesh.jl 在 czm.jl 之前 include）。

### `abstract type AbstractDamageState` — L26

CZM 损伤状态抽象基类，同上 include 顺序考虑。

### `struct CzmSubmesh` — L43-L48

独立细化的 CZM 机械子网格（径向 8 层/卷绕圈），与粗热网格解耦，通过 `thermal_elem_map` 与 `thermal_to_czm` 矩阵耦合。

- 字段：`mesh::Mesh`（细化 Q4 网格）、`material_type::Vector{Symbol}`（:PE / :PCC / :SP / :NE / :NCC）、`winding_turn`（卷绕圈号）、`thermal_elem_map`（CZM 单元 → 粗热单元 id）
- 定义在 SetMesh.jl 而非 czm.jl 的原因：`CohesiveMesh` 引用 `CzmSubmesh`，而 SetMesh.jl 在 czm.jl 之前被 include

### `mutable struct CohesiveMesh` — L50-L74

CZM 装配网格容器，聚合 bulk 网格 + 内聚力单元 + 节点分层映射 + 损伤状态 + 子网格耦合矩阵。

- 字段：`bulk_mesh`、`node`（扩展后节点）、`nnode`、`bulk_element`（更新后的固体连接）、`cohesive_elements::Vector{AbstractCohesiveElement}`、`n_cohesive`、`n_layers`（分离面类型数，PE-PCC + NE-NCC = 2）、`node_map`（原节点 → 分层后的节点们）、`interface_nodes`（每界面的节点对）、`damage_states::Vector{AbstractDamageState}`、`czm_submesh`（v5 细化子网格，可为 nothing）、`thermal_to_czm`（v5 粗热 → 细 CZM 插值矩阵）、`cohesive_to_thermal`（v5 反向映射）
- 内部构造函数 `CohesiveMesh()`（L66-L73）：空初始化，所有数组为 0 维，新字段（v5）默认 `nothing`

## 函数清单

### `SetMesh(domain, num, type, gsorder) -> Mesh` — L76-L98

网格生成入口（dispatcher），按 `type` 路由到 `Mesh1D` / `Mesh2D`。

- `type ∈ ["L2", "L3"]` → `Mesh1D`
- `type == "Q4"` 或 `"COH2D4"` → `Mesh2D`
- 其他抛 `error`
- docstring 备注（L85）"There should be more sophisticated methods to build a mesh, to be finished"

### `Mesh2D(domain, num, type, gsorder) -> Mesh` — L100-L150

构建规则矩形区域的二维 Q4 网格。

- L107：`@assert length(domain) == 4`
- L108-L113：`num` 支持 Tuple / Vector 两种形式
- L114-L128：节点编号 `idx(i,j) = j·nnx + i + 1`（i,j 从 0 开始）
- L129-L146：单元连接 `(i,j) → (i+1,j) → (i+1,j+1) → (i,j+1)`，逐元素赋值避免 Julia 行切片广播限制
- L147-L148：调 `GetGS` 构造高斯点
- 跨文件依赖：`GetGS`

### `Mesh1D(domain, num, type, gsorder) -> Mesh` — L152-L189

构建一维 L2 / L3 网格。

- L156-L160：按 `type` 决定每单元节点数（L2 → 2，L3 → 3）
- L161-L174：节点坐标按 `num` 分段均匀生成，首段包含起点，后续段跳过首点避免重复
- L176-L185：单元连接按 `num` 分段顺序编号
- L186-L187：调 `GetGS` 构造高斯点
- 跨文件依赖：`GetGS`

### `PickElement(mesh, v) -> Mesh` — L191-L227

从网格中按索引抽取子集，重新编号节点/单元/高斯点。

- L202-L208：抽取指定单元，建立 `node_pool`（去重节点）和 `node_pool_pair`（旧 → 新节点编号）
- L210-L223：高斯点重新映射——按 `gsorder^dim` 每单元的点数，重排 `gs.ele` / `gs.x` / `gs.xi` / `gs.weight` / `gs.detJ` / `gs.Ni` / `gs.dNidx`
- L224-L226：节点数重算（`length(unique(element))`）并封装 `Mesh`
- 注：L224 重复赋值 `nlen`（L205 已赋值），可能是冗余

### `CombineMesh(meshes) -> Mesh` — L229-L298

合并多个同类型同维度同 gsorder 的网格为一个大网格。

- L236-L249：类型/维度/gsorder 一致性校验，不匹配抛 `error`
- L251-L265：预计算总单元数 / 节点数 / 高斯点数
- L267-L293：节点偏移 `v_node`，单元连接 `.+ v_node`，高斯点 `ele` 偏移 `v_ele`，逐字段拼接
- L295-L297：封装 `GaussPoint` 和 `Mesh`
- 跨文件依赖：`GaussPoint`、`Mesh`

### `MultipleMesh(mesh, n) -> Mesh` — L300-L313

复制单个网格 n 份并合并（调 `CombineMesh`）。

- 跨文件依赖：`CombineMesh`

### `GetGS(element, node, order, type, v) -> GaussPoint` — L315-L409

高斯积分点构造（核心积分基础设施）。

- L316-L363：`type == "COH2D4"` 分支——沿中线（节点 1-4 中点 ↔ 节点 2-3 中点）一维积分，`detJ = L/2`；`L < 1e-15` 时 `continue` 跳过零长度单元；调 `NCweight` 取 Newton-Cotes 权重
- L365-L378：其他类型分支——按 `type` 决定维度和节点索引 `points`
- L386-L399：单元循环——调 `LagrangeBasis` 取形函数，`J0 = dNdxi · node[sctr]`，`detJ = det(J0)`
- L400-L406：按 `type` 调 `ShapeFunction1D` 或 `ShapeFunction2D` 构造导数
- 跨文件依赖：`NCweight`、`GSweight`、`LagrangeBasis`、`ShapeFunction1D`、`ShapeFunction2D`

### `LagrangeBasis(type, dimen, coord) -> (N, dNdxi)` — L411-L464

Lagrange 形函数及其局部导数（单点求值，非批量）。

- L414-L420：`L2` / `L3` —— 2 节点线性形函数
- L421-L432：`Q4` —— 4 节点双线性形函数
- L433-L458：`B8` —— 8 节点六面体形函数
- L461-L462：返回 `N'` 和 `dNdxi'`（转置为行向量/矩阵）
- 注：L441-L442 `I1 = 1/2 - coord/2; I2 = 1/2 + coord/2` 使用整个 `coord` 数组（3 元素），依赖 Julia 广播

### `GSweight(order, dimen) -> (W, Q)` — L467-L623

Gauss-Legendre 积分点和权重（order 1-10）的张量积生成。

- L468-L470：`order > 10 || order < 0` 时调 `disp(...)`（非 `error`，可能静默）
- L475-L592：1D 点/权重的硬编码表（order 1-8 精确值，order 9-10 用 else 分支的 10 位精度值）
- L593-L621：张量积——dimen=1 直序、dimen=2 双重循环、dimen=3 三重循环
- 跨文件依赖：`disp`（外部，非标准 Julia 函数，可能未定义）

### `NCweight(order) -> (w, q)` — L625-L644

Newton-Cotes 积分权重（order 2-5），用于 COH2D4 单元。

- L626-L628：`order < 2 || order > 5` 抛 `error`
- L630-L642：硬编码 4 种阶的点/权重

### `ShapeFunction1D(element, type, node, xi, v) -> (Ni, dNidx)` — L646-L682

1D 单元（L2/L3）批量形函数及其全局导数。

- L647-L666：`L3` 分支——3 节点二次形函数（注释保留了两组等价表达式），导数乘 `dXdx = 2/ele_length` 转全局
- L667-L677：`L2` 分支——2 节点线性形函数
- L678-L680：其他类型抛 `error`
- 注：L658 `df1 = x -> x .- 1/2` 的运算符优先级（`1/2` 先于 `.-`），等价于 `x .- 0.5`，非 `x .- 1` 后除 2；L660 同理

### `ShapeFunction2D(element, type, node, xi, ele_map) -> (Ni, dNidx)` — L692-L767

2D 单元（Q4/COH2D4）批量形函数及其全局导数。

- L693-L718：`Q4` 分支——每高斯点计算局部形函数 `Nloc`、局部导数 `dN_dxi`、Jacobi `J = dN_dxi' · x_e`、全局导数 `dN_dx = dN_dxi · inv(J)`；输出 `dNidx` 前 4 列 dN/dx、后 4 列 dN/dy
- L719-L763：`COH2D4` 分支——沿中线一维形函数 `N1 = 0.5(1-ξ)`、`N2 = 0.5(1+ξ)`；`L < 1e-15` 时 `continue` 跳过零长度单元；切向导数 `dN/ds = (2/L)·dN/dξ`，再投影到 x/y
- L764-L766：其他类型抛 `error`

## 省略项

无。所有函数与 struct 均独立列出。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L469 | `disp("Order of quadrature too high for Gaussian Quadrature")` | `disp` 非 Julia 标准函数（MATLAB 残留），order > 10 或 < 0 时调用；可能抛 `UndefVarError` 或静默无输出，疑似未完成的错误处理 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L85 | `# There should be more sophisticated methods to build a mesh, to be finished`（docstring 中的 TODO 注释） | 显式 "to be finished" 标记，提示 `SetMesh` dispatcher 未来需扩展更复杂的网格生成方法 |
| L342 | `if L < 1e-15; continue; end`（COH2D4 零长度单元静默跳过，`GetGS` COH2D4 分支，跨 L342-L344） | 零长度 cohesive 单元被静默跳过，不报错；合法网格不应出现，但若出现会导致高斯点数组尾部为 0 值，下游积分异常且无警告 |
| L745 | `if L < 1e-15; continue; end`（COH2D4 零长度单元静默跳过，`ShapeFunction2D` COH2D4 分支，跨 L745-L747） | 同 L342，形函数层面的静默跳过；与 `GetGS` 的跳过配对，但两者阈值和位置独立，可能出现不一致 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|------|
| L88 | `if type in ["L2", "L3"] ... elseif type == "Q4" ... elseif type == "COH2D4" ... else error(...)`（4 分支 if-elseif 链，每分支单条件，跨 L88-L96） | 仅 4 分支单条件，未达 ≥3 条件阈值；可用 Dict/dispatch 替代但当前清晰度可接受 |
| L365 | `if type == "L2" ... elseif type == "L3" ... elseif type == "Q4" ... elseif type == "B8"` 后接 `total_num = size(element,1) * order ^ dimen`（4 分支设置 dimen/points，跨 L365-L377） | 同上，多分支但非嵌套条件链；可抽出 `element_dim_points(type) -> (dimen, points)` helper 提高可读性 |
