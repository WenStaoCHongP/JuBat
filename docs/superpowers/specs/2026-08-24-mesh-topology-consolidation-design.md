# 热网格与力学网格拓扑收敛设计

- 日期：2026-08-24
- 状态：设计已确认，待生成实施计划
- 前置审计：`docs/planning-with-files/32_热力网格拓扑优化审计/findings.md`

---

## 1. 背景与问题

JuBat 的 Jellyroll 模型当前同时持有五个拓扑对象，多个字段名与实际语义脱节，映射层散落且方向不明，可变损伤状态与不可变拓扑混装在同一个 `mutable struct` 中。

### 1.1 五个拓扑对象

| 对象 | 内容 | 是否参与求解 |
|---|---|---|
| `JellyrollMesh.thermal2D` | 未合并粗热 Q4 网格 | 否（仅构造期中间产物） |
| `JellyrollMesh.thermal2D_merged` | Φ 合并后的连续热网格 | 是（活动热网格） |
| `CzmSubmesh.mesh` | 未合并力学 Q4 网格 | 否 |
| `CzmSubmesh.mesh_bonded` | Φ 粘结后的力学 bulk 网格 | 是（cohesive 拓扑基础） |
| `CohesiveMesh` | 界面复制节点后的最终 bulk + COH2D4 | 是（最终力学求解拓扑） |

在 nθ=80 下，两个未合并网格分别常驻约 1.13 MB 和 8.85 MB；网格构造过程调用 `GetGS` 四次，其中两次产出的高斯数据没有任何生产消费者。

### 1.2 失效映射（正确性缺陷）

`JellyrollMesh.czm_element_map`（`src/Jellyrollmodel.jl:132-158`）在粗热网格构造阶段由相邻 Φ 段生成，值域是 `1:(length(interface_pairs)-1)`，语义是"跨匝 Φ 段附近的热单元关系"。

但 `src/CallModel.jl:85` 与 `src/Materialmatrix.jl:391` 把它的值直接与 `get_fractured_elements(czm_mesh)` 返回的真实 cohesive 单元 id 比较：

```julia
for czm_idx in get(geom.czm_element_map, e, Int64[])
    if czm_idx in fractured_czm
```

两者的值域和计数都不同（cohesive 单元数为 `4 * (length(theta)-1)`），该比较依赖编号偶然重叠，不是有效拓扑映射。当前工况 `D = 0`，`get_fractured_elements` 恒返回空集，该路径休眠；一旦发生断裂，它会停用错误的热/电化学单元。

同族缺陷还有一处：`src/CsvExport.jl:467` 的 `write_cohesive_driving_force_csv` 把 `geo.interface_pairs`（节点对）当作 `(e_top, e_bot)` 热单元对按 cohesive id 索引，而默认合并路径下该字段恒为空，导致 `cohesive_driving_force.csv` 长期只输出表头。本次不处理该项（见 §8）。

### 1.3 字段名与语义脱节

| 字段 | 声称语义 | 实际内容 |
|---|---|---|
| `CohesiveMesh.bulk_mesh` | "原始固体网格" | 指向 `czm_submesh.mesh_bonded`；最终求解拓扑是 `node` + `bulk_element`，二者不同 |
| `CohesiveMesh.n_layers` | 层数 | cohesive 本构类型数（恒为 2） |
| `MeshGeometry.interface_pairs` | 注释写 `(top_elem, bot_elem)` | 实际来源是节点对，且默认路径下为空 |
| `MeshGeometry.element_layer` | 注释写 `(1=NE, 2=SP, 3=PE, 4=NCC, 5=PCC)` | 实际是卷绕匝号 `floor((r_c - Rin) / layer) + 1`，与材料类型无关 |
| `CzmSubmesh.thermal_elem_map` | 目标不明确 | bulk 力学单元 → 父热单元 |

### 1.4 拓扑与可变状态混装

`CohesiveMesh` 同时持有拓扑（`node`、`bulk_element`、`cohesive_elements`、映射）和可变历史（`damage_states`）。更新损伤必须调用 `clone_czm_mesh_with_damage` 手工浅拷贝全部拓扑字段，每新增字段都可能漏拷贝。

`CZMAssemblyCache` 以 `objectid(czm_mesh)` 判断拓扑失效（`src/CouplingState.jl:243-244, 257`），但 `CohesiveMesh` 是 `mutable struct` 且内部数组可原位修改；发生原位变化时 `objectid` 不变，缓存的 `K_bulk`、DOF 映射、cohesive frame 和边界 DOF 可能全部过期。

### 1.5 构造入口分散

用户脚本需手工完成三步接线（`example/testexample.jl:79-85`）：

```julia
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, czm_enabled=true, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
```

第三步把活动热网格回传给 `create_czm_mesh`，而该参数只用于一次 `n_thermal` 越界检查。

### 1.6 死字段

以下字段在 `src/` 中无生产消费者（已完成仓库级检索）：

- `JellyrollMesh.merge_map`
- `JellyrollMesh.pos_tab_nodes` / `neg_tab_nodes`（极耳冷却在 `src/ThermalDistributed.jl:119` 用活动网格运行时重算）
- `JellyrollMesh.inner_nodes` / `outer_nodes` 及其在 `MeshGeometry` 中的副本（`ThermalRing2D_BC` 的 `outer_nodes` 实参来自 `variables["thermal2D outer_nodes"]`，由圆环验证脚本从 `src/ring.jl` 的网格数据提供，与 Jellyroll 路径无关）
- `CzmSubmesh.phi_keep`
- `CohesiveMesh.node_map` / `interface_nodes`
- `MeshGeometry.element_layer` / `is_inner_layer`（唯一消费者是 `get_active_elements`，本设计将其删除）

---

## 2. 范围

### 2.1 在范围内

1. 拓扑对象、字段、映射层的重塑
2. 不可变拓扑与可变状态分离，缓存键改为内容签名
3. 失效映射 `czm_element_map` 及其消费者的删除
4. 构造入口收敛为单一函数
5. 技术文档与源码函数索引同步

### 2.2 不在范围内

1. **网格生成算法与离散策略不变**：collector-seeded 螺旋布种、径向 8 层分层、nθ 共享机制、`thin_subdiv` 全部保持原样
2. **不接入损伤对电化学/热的反馈**：删除的是用错映射的那条通道，不新建正确通道
3. **不处理 `src/CsvExport.jl`**：用户已决定缓期（见 §8）
4. **不改变任何数值算法**：所有批次以零漂移为目标

### 2.3 兼容策略

用户已明确：**先建新结构，再直接删旧兼容**，不保留 facade，不设弃用窗口。仓库内部是唯一消费者，所有调用点与测试同步更新。

四个公开类型名 `JellyrollMesh` / `CzmSubmesh` / `CohesiveMesh` / `MeshGeometry` 全部保留，只重塑字段。真正误导人的是字段名而非类型名，保留类型名可消除约 20 个测试文件、examples 和 tools 中与零漂移验证无关的 diff 噪声。

---

## 3. 目标结构

### 3.1 三个拓扑对象

| 对象 | 唯一职责 |
|---|---|
| `JellyrollMesh.thermal2D` | 活动热网格（Φ 已合并的连续导热拓扑） |
| `CzmSubmesh.mesh` | 力学 bulk 网格（Φ 已粘结，8 层分层 Q4） |
| `CohesiveMesh` | 最终力学求解拓扑（界面复制节点后的 bulk + COH2D4） |

### 3.2 `JellyrollMesh`（14 → 2）

```julia
struct JellyrollMesh
    thermal2D::Mesh                            # 活动热网格（原 thermal2D_merged）
    czm_submesh::Union{Nothing, CzmSubmesh}
end
```

删除：`thermal2D`(未合并)、`merge_map`、`interface_pairs`、`czm_element_map`、`element_layer`、`is_inner_layer`、`inner_nodes`、`outer_nodes`、`pos_tab_nodes`、`neg_tab_nodes`、`ne`、`nnode`。

未合并热网格降为 `jellyroll_collector_seed_mesh` 的构造期局部变量。

### 3.3 `CzmSubmesh`（7 → 5）

```julia
struct CzmSubmesh
    mesh::Mesh                                 # Φ 粘结后的力学 bulk 网格（原 mesh_bonded）
    material_type::Vector{Symbol}
    winding_turn::Vector{Int}
    bulk_to_thermal::Vector{Int}               # 原 thermal_elem_map
    phi_seam::Vector{PhiSeam}                  # 原 phi_pairs，改为逻辑列表示
end

struct PhiSeam
    outer_col::Int    # 外螺旋角列号（θ）
    inner_col::Int    # 内螺旋角列号（θ + 2π）
    node::Int         # 合并后网格中的共享节点号
end
```

删除：`mesh`(未合并)、`phi_keep`、5 参数兼容构造器。

`phi_pairs` 现用未合并网格节点号表达，这是保留第二份完整力学网格的唯一技术理由。改用逻辑列号加合并后节点号后信息无损，且更贴合 AGENTS.md §9.8 对未来 SP–涂层接触装配复用该配对的要求。

### 3.4 `CohesiveMesh`（12 → 8，改为不可变）

```julia
struct CohesiveMesh
    node::Matrix{Float64}
    nnode::Int64
    bulk_element::Matrix{Int64}
    cohesive_elements::Vector{AbstractCohesiveElement}
    n_cohesive::Int64
    czm_submesh::CzmSubmesh
    cohesive_to_thermal::Vector{Int}
    topology_id::UInt64                        # 构造时按 node/element/cohesive 内容哈希
end
```

删除：`bulk_mesh`（高斯积分阶次改从 `czm_submesh.mesh.gs.order` 取）、`n_layers`、`node_map`、`interface_nodes`；`damage_states` 迁出至 `CzmLayout`。

`Union{Nothing, ...}` 包装去除：构造链保证 `czm_submesh` 与 `cohesive_to_thermal` 必然存在，空构造器（`src/CzmSolve.jl:32`）随之改写。

### 3.5 `CzmLayout`（+1）

```julia
mutable struct CzmLayout
    n_coh::Int
    ndof::Int
    u_prev::Vector{Float64}
    plastic_states::Union{Nothing, Matrix{PlasticState}}
    winding_prestress::Union{Nothing, Vector{NTuple{3, Float64}}}
    node_ref::Union{Nothing, Matrix{Float64}}
    damage_states::Vector{AbstractDamageState}   # 新增，自 CohesiveMesh 迁入
end
```

不新建类型：`CzmLayout` 本就是力学可变状态的容器（已持有位移、塑性、预应力、参考构型），损伤状态归入此处后 `CohesiveMesh` 自然成为纯不可变拓扑，`clone_czm_mesh_with_damage` 整个函数消失。

### 3.6 `MeshGeometry`（8 → 7）

```julia
struct MeshGeometry
    layer_weights::Matrix{Float64}             # ne × 5 层面积权重 [NE, SP, PE, PCC, NCC]
    boundary_edges::BoundaryEdgeCache
    inner_nodes::Vector{Int}                   # 活动网格上的内螺旋边界节点
    outer_nodes::Vector{Int}                   # 活动网格上的外螺旋边界节点
    pos_tab_nodes::Vector{Int}                 # 活动网格上的正极耳节点（保持 jellyroll_tab_node_indices 原始顺序）
    neg_tab_nodes::Vector{Int}                 # 活动网格上的负极耳节点（保持原始顺序）
    interface_pairs::Vector{Tuple{Int,Int}}    # 休眠字段，合并路径下恒为空；仅为 CsvExport 缓期保留，见 §8
end
```

删除：`element_layer`、`is_inner_layer`、`czm_element_map`。

**节点集来源改为活动网格。** 原 `inner_nodes`/`outer_nodes` 来自未合并热网格（3366 节点），而活动网格是合并后的（1763 节点），编号不通用，直接索引会取到错误节点。`setup_thermal2D_mesh` 已经在活动网格上调用 `identify_boundary_nodes`（`src/Jellyrollmodel.jl:564`），但 `is_inner` 算完即丢弃；本设计将其保留为 `findall(is_inner)`，并新增一次 `jellyroll_tab_node_indices(mesh_th, param)` 缓存极耳节点。

这同时消除了极耳冷却的每步重算：`apply_cool_method` 目前在每次 `ThermalDistributed2D_BC` 调用时重跑 `jellyroll_tab_node_indices`（每个极耳 80 次二分搜索后扫描全部节点），Crank-Nicolson 下每个时间步两次。

**顺序契约**：`apply_cool_method` 按 `unique(vcat(pos_idx, neg_idx))` 的顺序计算弧长权重。缓存必须原样保存 `jellyroll_tab_node_indices` 返回的两个数组，消费端仍执行相同的 `unique(vcat(...))`，不得预先合并、排序或去重，否则权重分配改变。

### 3.7 字段总账

| 类型 | 现在 | 之后 |
|---|---:|---:|
| `JellyrollMesh` | 14 | 2 |
| `CzmSubmesh` | 7 | 5 |
| `CohesiveMesh` | 12 | 8 |
| `MeshGeometry` | 8 | 7 |
| 合计 | 41 | 22 |

拓扑对象 5 → 3，`GetGS` 调用 4 → 2，nθ=80 常驻内存减少约 8.8 MB。

---

## 4. 批次 B1：删死重量

目标：删除无生产消费者的数据与恒等变换路径。全部四组均可独立论证零漂移。

### 4.1 B1.1 删除用错映射的损伤通道

当前存在两条"损伤 → 电化学/热"通道：

| 通道 | 映射来源 | 正确性 | 处置 |
|---|---|---|---|
| 断裂单元失效屏蔽 | `czm_element_map`（Φ 段序号） | 错误 | 删除 |
| 渐进式有效面积损失 | `cohesive_to_thermal`（真实 cohesive → 热单元） | 正确 | 保留 |

删除清单：

| 项 | 位置 |
|---|---|
| `JellyrollMesh.czm_element_map` 字段 | `src/Jellyrollmodel.jl:6` |
| `czm_element_map` 构造 | `src/Jellyrollmodel.jl:132-158` |
| `JellyrollMesh` 构造实参 | `src/Jellyrollmodel.jl:200` |
| `MeshGeometry.czm_element_map` 字段 | `src/CouplingState.jl:168` |
| `setup_thermal2D_mesh` 传参 | `src/Jellyrollmodel.jl:572` |
| `get_active_elements` 函数与 export | `src/Materialmatrix.jl:382-400`、`src/JuBat.jl:79` |
| `compute_heat_sources_with_czm` 函数与 export | `src/ThermalDistributed.jl:390-418`、`src/JuBat.jl:62` |
| `deactivated_elements` 生产块 | `src/CallModel.jl:80-94` |
| `solve_branch_currents` 的 `deactivated_elements` 形参 | `src/Parallelsolution.jl` |

`src/CallModel.jl:145` 的分支调用统一归并为 `compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas; per_element_spme=true)`。

`get_fractured_elements` 保留（是对损伤状态的正确查询），但在 B2 中改为作用于 `CzmLayout.damage_states`。

**零漂移论证**：基线档案记录 `czm_fractured_elements = 0`、`czm_d_max_percent = 0.0`，`get_fractured_elements` 恒返回空集。空集下 `deactivated_elements` 为空、`get_active_elements` 返回全集、`compute_heat_sources_with_czm` 的掩码循环不写入任何值，故等价于直接调用 `compute_heat_sources`。删除是恒等变换。

**边界声明**：该等价性**只在 D = 0 工况下成立**，不是永久等价。将来恢复断裂失效屏蔽时必须按 `cohesive_to_thermal` 构造反向索引 `thermal_to_cohesive::Vector{Vector{Int}}`，不得照抄被删逻辑（见 §7 被删路径登记表）。

### 4.2 B1.2 删除未合并力学网格

| 项 | 处置 |
|---|---|
| `CzmSubmesh.mesh`（未合并） | 删除 |
| `CzmSubmesh.phi_keep` | 删除 |
| `CzmSubmesh.mesh_bonded` | 改名 `mesh` |
| `CzmSubmesh.phi_pairs` | 改为 `phi_seam::Vector{PhiSeam}` |
| `CzmSubmesh.thermal_elem_map` | 改名 `bulk_to_thermal` |
| 5 参数兼容构造器 | 删除（`src/SetMesh.jl:53-57`） |
| `merge_phi_pairs` | 改为直接产出粘结拓扑与 seam 记录，不再先建中间 `Mesh` |
| `create_czm_mesh` 中 `sub_mesh = czm_submesh.mesh_bonded` | 改为 `czm_submesh.mesh` |

受影响测试：`test/test_czm_submesh.jl`、`test/test_czm_phi_merge.jl`、`test/test_czm_thin_subdiv.jl`、`test/test_create_czm_mesh.jl`。受影响工具：`tools/czm_mesh_probe.jl:66`。

**零漂移论证**：粘结网格的 `node_b` 与 `element_b` 完全由未合并 node/element 数组加 merge_map 决定，跳过中间 `Mesh` 对象不改变这两个数组；`GetGS(element_b, node_b, gsorder, "Q4")` 的输入逐位不变，输出逐位相同。被删的 `GetGS(element, node, gsorder, "Q4")`（`src/Jellyrollmodel.jl:743`）产出的高斯数据无生产消费者。

### 4.3 B1.3 删除未合并热网格与死字段

| 项 | 处置 |
|---|---|
| `JellyrollMesh.thermal2D`（未合并） | 降为构造期局部变量 |
| `merge_map`、`interface_pairs`、`inner_nodes`、`outer_nodes`、`pos_tab_nodes`、`neg_tab_nodes`、`ne`、`nnode`、`element_layer`、`is_inner_layer` | 从 `JellyrollMesh` 删除 |
| `MeshGeometry.element_layer`、`is_inner_layer` | 删除 |
| `MeshGeometry.inner_nodes`、`outer_nodes` | 改为活动网格来源 |
| `MeshGeometry.pos_tab_nodes`、`neg_tab_nodes` | 新增，活动网格来源 |
| `setup_thermal2D_mesh` 的 `use_merged` 关键字 | 删除（v2 起强制 `true`，`else` 分支不可达） |
| `apply_cool_method` 的 tab 分支 | 改为消费 `case.geometry.pos_tab_nodes` / `neg_tab_nodes` |

受影响 example：`example/coupled_czm_thermal_example.jl`（引用 `interface_pairs`、`is_inner_layer`、`n_layers`、`element_layer`、`get_active_elements`）、`example/网格敏感性_v2/4_czm_convergence.jl:135`。受影响工具：`tools/verify_czm_standalone.jl:66,137`。

**零漂移论证**：未合并热网格在构造期只被用于三件事——Φ 配对重合性校验、生成 merge_map、向 `build_czm_submesh` 提供 `size(element, 1)`——三者均只读 node/element，不触及高斯数据。被删的 `GetGS`（`src/Jellyrollmodel.jl:92`）无生产消费者。极耳节点改为活动网格上预计算，保持 §3.6 的顺序契约后，`apply_cool_method` 的弧长权重与 K/F 装配逐位不变。

### 4.4 B1 完成后的收益

- 拓扑对象 5 → 3
- `GetGS` 调用 4 → 2
- nθ=80 常驻内存减少约 8.8 MB
- 消除失效映射引发的断裂后错误屏蔽风险
- 极耳冷却从每步重算改为构造期一次

---

## 5. 批次 B2：不可变核与状态分离

### 5.1 拓扑不可变化

`CohesiveMesh` 从 `mutable struct` 改为 `struct`，按 §3.4 重塑字段，新增 `topology_id::UInt64`（构造时按 `node`、`bulk_element`、`cohesive_elements` 的节点/类型内容哈希）。

`CZMAssemblyCache.czm_mesh_id` 改名 `topology_id`，失效判据从 `objectid(czm_mesh)` 改为内容签名。空构造器 `CohesiveMesh()`（`src/CzmSolve.jl:32`）改写为显式构造或去除。

### 5.2 状态迁移

`damage_states` 从 `CohesiveMesh` 迁入 `CzmLayout`。随之：

| 项 | 处置 |
|---|---|
| `clone_czm_mesh_with_damage` | 删除（调用点 `src/CzmSolve.jl:295,498,680`、`src/CzmPostProcess.jl:61,92`） |
| `solve_czm_step` 返回值 | 不再返回 `updated_czm_mesh`，改为写 `case.czm_layout.damage_states` |
| `CycleSolver` 跨相位传递 | 从传 `czm_mesh` 改为传 `czm_layout` |
| `get_fractured_elements`、`compute_element_gap_conductance`、`compute_all_gap_conductances`、`map_czm_damage_to_thermal` | 改为从 `CzmLayout.damage_states` 取状态 |

### 5.3 风险

这是三个批次中**唯一可能改变物理行为**的项。损伤跨周期累积与失败步回滚现在依赖"传递或不传递 `czm_mesh` 对象"实现，迁移后必须保持完全一致的生命周期语义。

守门测试：`test/test_czm_multicycle_c4lite.jl`（跨相位状态、受控损伤）、失败步提交/回滚专项、`test/test_ensure_czm_cache.jl`、`test/test_cache_invariants.jl`、`test/test_cohesive_struct.jl`。

---

## 6. 批次 B3：构造入口收敛与文档

### 6.1 单一构造入口

```julia
setup_jellyroll_mesh!(case; nθ, gsorder=2, czm_enabled=false, thin_subdiv=1,
                      phase=0.0, tol=1e-8) -> case
```

一次写好 `case.mesh["thermal2D"]`、`case.geometry`、`case.czm_mesh`、`case.czm_layout`，取代现有三步接线。

`create_czm_mesh` 去掉 `thermal_mesh` 参数——它仅用于 `thermal_elem_of_outer <= n_thermal` 越界检查，该不变量改由 `build_czm_submesh` 在构造 `bulk_to_thermal` 时保证，并由 §7.2 的新增测试守护。

`jellyroll_collector_seed_mesh`、`setup_thermal2D_mesh`、`create_czm_mesh` 三个旧入口从公开 API 删除；`create_czm_mesh` 因 `src/CzmUnitMesh.jl:63` 的单元条带构造仍需要，降为内部函数。

### 6.2 文档同步

- `md/02_几何与网格.md`
- `md/对照/02_几何与网格对照.md`
- `md/源码函数索引/`：`Jellyrollmodel.md`、`SetMesh.md`、`CzmMesh.md`、`CouplingState.md`、`ThermalDistributed.md`、`Materialmatrix.md`、`CallModel.md`、`JuBat.md`
- `CLAUDE.md` / `AGENTS.md` §9.3（`phi_pairs` 表述改为 `phi_seam`）
- `md/07_界面热阻模型.md`、`md/10_参数传递与模块架构.md`（引用被删的失效屏蔽路径）

### 6.3 调用点更新

所有 `example/` 与 `tools/` 中的三步接线改为单一入口。

---

## 7. 验证

### 7.1 每批统一门禁

1. 模块加载：`include("src/JuBat.jl")` 成功
2. 该批受影响的专项测试
3. 全量 `test/runtests.jl`（当前 33 文件）
4. `example/testexample.jl` 强制基线（AGENTS.md §9.6）：
   - `exit_code = 0`
   - 网格 `theta_elements = 80`、`thermal_elements = 1682`、`thermal_nodes = 1763`
   - `time_steps = 19`
   - `metrics.toml` 全部科学指标按记录精度一致，含 `czm_max_normal_separation_m = 1.5174e-12`、`czm_converged_updates = 19`
   - PNG SHA-256 = `0946646ac91ef9493ea09a2ca199a7495573767b1ad2c7f7375919e8f290a447`
5. `tools/verify_czm_standalone.jl` 快照对照 `docs/planning-with-files/30_堆芯塌陷力学建模/baseline_czm_standalone.md`
6. 内存与构造探针：`Base.summarysize` 前后对比、`GetGS` 调用计数

固定运行环境：Julia 1.11.2、`JULIA_NUM_THREADS=1`、`GKSwstype=100`、`--startup-file=no`、`--project=.`。

**source manifest 不作为漂移判据**：基线档案记录 `source_manifest_tsv_sha256` 与 `source_file_count = 46`，每批都会修改 `src/`，清单哈希必然变化。判定漂移的只有退出码、网格/步数、科学指标、PNG SHA-256 四项。清单按 `Simplify/baseline/testexample/source_manifest.tsv` 既定算法（去表头与 aggregate 尾行，`path<TAB>sha256` 数据行以 LF 连接、末尾不加换行，取 SHA-256）机械刷新。

### 7.2 新增测试

被删路径守护的不变量必须由新测试接管：

| 测试 | 内容 | 批次 |
|---|---|---|
| 极耳缓存等价 | `cool_method="tab"` 下缓存的 `pos/neg_tab_nodes` 与运行时 `jellyroll_tab_node_indices(活动网格, param)` 逐位一致，且装配出的 K/F 逐位一致 | B1.3 |
| 活动网格边界节点集 | `inner_nodes`/`outer_nodes` 非空，节点半径分别落在 `Rin`/`Rout` 容差内 | B1.3 |
| Φ seam 完整性 | `length(phi_seam)` 与每匝角节点数派生值一致；每个 `node` 在粘结网格编号范围内；`outer_col`/`inner_col` 在角列范围内 | B1.2 |
| `cohesive_to_thermal` 值域 | 全部落在 `1:n_thermal`（接管被删的越界断言） | B1.2 / B3 |
| `bulk_to_thermal` 值域 | 全部落在 `1:n_thermal`，且每个 bulk 单元有且仅有一个父热单元 | B1.2 |
| `topology_id` 敏感性 | 改动任一节点坐标或连接后 id 必变；相同输入构造两次 id 相同 | B2 |
| 损伤生命周期 | 跨周期累积与失败步回滚在状态迁移前后行为一致 | B2 |

### 7.3 受影响的既有测试

| 测试 | 原因 | 批次 |
|---|---|---|
| `test/test_cache_invariants.jl:30-31` | 冻结了 `czm_element_map` 字符串检查 | B1.1 |
| `test/test_map_czm_damage.jl` | 证明保留的 `cohesive_to_thermal` 通道未变 | B1.1 |
| `test/test_czm_submesh.jl` | 冻结双网格字段契约与 5 参数构造器 | B1.2 |
| `test/test_czm_phi_merge.jl` | 冻结 `phi_pairs`/`phi_keep`/`mesh_bonded` 计数 | B1.2 |
| `test/test_czm_thin_subdiv.jl` | 比较 `czm_submesh.mesh` 与 `phi_pairs` | B1.2 |
| `test/test_create_czm_mesh.jl` | 构造、映射、cohesive 法向 | B1.2 / B3 |
| `test/smoke_thermal_bc.jl` | 冷却路径 | B1.3 |
| `test/test_cohesive_struct.jl` | `CohesiveMesh` 字段契约 | B2 |
| `test/test_czm_multicycle_c4lite.jl` | 跨相位状态 | B2 |
| `test/test_ensure_czm_cache.jl` | 缓存失效判据 | B2 |

### 7.4 停止条件

任一强制指标不一致时立即停止该批次，定位或回退，不以"数值接近"放行（AGENTS.md §9.6）。科学结果、网格/步数、PNG SHA-256 任一变化都表示该批改动不是恒等变换。

---

## 8. 被删路径登记表

用户已决定源码不保留死注释，被删路径的原貌与重建方式在此登记。

### 8.1 断裂单元失效屏蔽

- **被删位置**：`src/Materialmatrix.jl:382-400`（`get_active_elements`）、`src/ThermalDistributed.jl:390-418`（`compute_heat_sources_with_czm`）、`src/CallModel.jl:80-94`（`deactivated_elements` 生产）
- **原逻辑**：对每个热单元 `e`，遍历 `czm_element_map[e]` 得到一组索引，若其中任一索引出现在 `get_fractured_elements(czm_mesh)` 中，则将该热单元标记为失活；失活单元的热源置零，并从分流求解中剔除
- **为什么错**：`czm_element_map` 的值域是 Φ 段序号 `1:(n_pairs-1)`，`get_fractured_elements` 返回的是 cohesive 单元 id `1:(4*n_segments)`，两个编号空间不同，比较依赖偶然重叠
- **为什么删除是安全的**：当前工况 `D = 0`，`get_fractured_elements` 恒返回空集，整条路径是恒等变换
- **正确的重建方式**：由 `CohesiveMesh.cohesive_to_thermal`（cohesive 单元 → 热单元，一对一）构造反向索引 `thermal_to_cohesive::Vector{Vector{Int}}`（热单元 → 其上的全部 cohesive 单元），再用它做失活判定
- **重建前提**：损伤对电化学/热的反馈理论闭合，且用户明确授权（该反馈当前在范围外）

### 8.2 `MeshGeometry.interface_pairs`（本次保留）

- **现状**：合并路径下 `setup_thermal2D_mesh` 将其设为空数组，`src/CsvExport.jl:467` 据此构造 `coh_thermal_map`，结果恒为空，`cohesive_driving_force.csv` 只输出表头
- **为什么错**：`interface_pairs` 是未合并热网格上的**节点**对，而消费端按 `(e_top, e_bot)` **单元**对使用，且按 cohesive id 索引
- **本次处置**：保留该字段（恒为空），保持现有行为不变。用户决定缓期处理，原因是 `src/CsvExport.jl` 存在未提交的工作区修改
- **将来的正确做法**：改用 `CohesiveMesh.cohesive_to_thermal` 取单个父热单元（而非 top/bot 对），或删除该导出函数

### 8.3 `use_merged=false` 与界面热阻分支

- **现状**：`setup_thermal2D_mesh` 的 `use_merged` 自动选择逻辑与 `ThermalDistributed2D_BC` 中的界面热阻块自 2026-07-21 起被注释禁用（`src/Jellyrollmodel.jl:537-549`、`src/ThermalDistributed.jl:166-190`）
- **本次处置**：删除 `use_merged` 关键字与不可达的 `else` 分支。`ThermalDistributed2D_BC` 中的注释块属于 CZM–热双向耦合能力，不在本次范围内，保持原样
- **恢复条件**：界面热阻需要未合并热网格提供界面双侧节点。若将来恢复，须重新引入未合并热视图，且不得复用被删的 `inner_nodes`/`outer_nodes`（它们的编号对活动网格无效）

---

## 9. 待决项

以下两项需在生成实施计划前确定：

1. **工作区状态**：`src/` 现有 10 个文件已修改未提交（`CouplingState.jl`、`CsvExport.jl`、`CzmMesh.jl`、`CzmSolve.jl`、`CzmUnitMesh.jl`、`Jellyrollmodel.jl`、`JuBat.jl`、`SetMesh.jl`、`Tools.jl`、`czm.jl`）。基线 v4 即在此脏状态下冻结。是否先提交这些修改，以便每批有干净的起点与可归因的回退边界。
2. **CsvExport 缓期的解除时机**：`MeshGeometry.interface_pairs` 作为休眠字段保留，需确定何时随 CsvExport 一并处理。
