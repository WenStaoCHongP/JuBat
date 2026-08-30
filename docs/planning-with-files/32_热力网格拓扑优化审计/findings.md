# 热力网格拓扑优化审计：发现与决策

## 用户需求

- 判断当前热网格和力学网格复杂、混乱的拓扑能否优化。
- 优化应降低维护与理解成本，同时不破坏已有函数、科学行为和多物理场映射。

## 已知基线

- 项目采用 collector-seeded Jellyroll 热网格，并由共享角段派生力学/CZM 子网格。
- `thermal_elem_map` 负责力学 bulk 单元到父热单元的直接映射。
- `phi_pairs` 保存跨匝节点配对，但其存在不代表接触或摩擦方程已装配。
- 当前历史上同时出现未合并拓扑、Φ 完美粘结合并求解拓扑和 cohesive 双侧节点，需要先区分物理职责再决定是否合并对象。

## 初步判断

可以优化，但最有价值的方向预计不是强行只保留一个网格，而是：

1. 统一几何/拓扑的唯一所有者；
2. 将“基础拓扑”“求解视图”“场映射”拆成明确层次；
3. 消除重复存储或由多个对象分别推导同一事实；
4. 用构造时不变量替代消费端的隐式假设；
5. 通过兼容适配层分批迁移调用者。

## Phase 1 源码盘点：对象与构造链

### 当前对象层次

1. `JellyrollMesh.thermal2D`：两条螺旋边界组成的未合并粗热 Q4 网格。
2. `JellyrollMesh.thermal2D_merged`：按粗热 `interface_pairs` 合并重合 Φ 节点后的连续热网格；当前 `setup_thermal2D_mesh` 默认强制使用它。
3. `CzmSubmesh.mesh`：由同一 `theta` 派生、径向展开为 8 层（或 thin subdivision）的未合并力学 Q4 网格。
4. `CzmSubmesh.mesh_bonded`：将力学 `phi_pairs` 合并后的完美粘结求解拓扑。
5. `CohesiveMesh`：以 `mesh_bonded` 为 bulk 基础，再复制真实箔–涂层界面节点并改写外层 bulk 连接，形成最终 bulk + COH2D4 求解拓扑。

### 第一批可疑复杂度

- `JellyrollMesh` 同时保存两个完整 `Mesh`、`merge_map`、Φ 配对、热到“CZM”字典和边界元数据，混合了拓扑所有权、求解视图和后处理元数据。
- `CzmSubmesh` 同时保存未合并/合并两个完整 `Mesh`。源码注释表明未合并 `.mesh` 目前主要为旧插值 `mod` 定位保留，而真正力学求解消费 `mesh_bonded`；这是明显的迁移债务候选。
- `CohesiveMesh` 既保存 `bulk_mesh=mesh_bonded`，又保存经界面复制后的 `node`、`bulk_element`，所以 `bulk_mesh.node/element` 与最终求解节点/连接不相同；字段名“原始固体网格”容易被误用为最终 bulk 拓扑。
- 粗热 `czm_element_map` 是根据粗热 Φ 配对附近单元构造，索引对象是 Φ 段序号，不是后来 `create_czm_mesh` 创建的四类真实箔–涂层 cohesive 单元；名称可能制造错误语义。

### 当前物理必要性

- 粗热网格和分层力学网格不能直接合并：二者径向分辨率及材料语义不同。
- cohesive 节点复制是产生界面位移跳量的必要拓扑操作，不能删除。
- 未合并 Φ 配对仍需保留为物理身份/未来接触信息；但不一定需要为此永久保留第二个完整 Gauss 网格。
- 热网格当前强制使用合并拓扑，`use_merged=false` 和界面热阻路径属于被禁用但仍残留的兼容分支，需要先核对调用者再决定收敛接口。

## Phase 1 源码盘点：消费者与迁移债务

### 已确认生产消费者

- `case.mesh["thermal2D"]` 是热、电化学分流和结果输出的活动粗热网格，消费者很多，不能改名或替换其索引顺序而不提供迁移层。
- `thermal_elem_map` 是有效且高价值的父子拓扑：温差与正/负极 SOC 均通过它从粗热单元直接映射到分层力学 bulk 单元。
- `cohesive_to_thermal` 是真实 cohesive 单元到粗热单元的映射，当前被损伤面积代理路径消费；即使该反馈默认关闭，字段有明确语义和专项测试。
- 最终力学装配消费 `CohesiveMesh.node` + `bulk_element`；`bulk_mesh` 主要只提供高斯积分阶次。这证明 `bulk_mesh` 并不是最终求解拓扑本体。

### 明显迁移债务

- `T_nodes` 本身有明确下游：其父热单元平均值 `Te_prev` 驱动分流/SPMe/热源，父热单元平均温差经 `thermal_elem_map` 形成 `dT_czm` 并驱动热应变。无生产残差消费者的是额外生成的细力学节点场 `T_czm_nodes = thermal_to_czm * T_nodes`，而不是热节点温度本身。
- 用户已决定采用父热单元平均温度，不保留细力学节点温度场。因此本批不再把 `thermal_to_czm` 作为兼容诊断能力保留，而是删除该矩阵、返回字段及专用插值路径；必须用测试证明 `dT_czm` 和电化学 `Te_prev` 路径未被改变。
- 仓库调用点已完整枚举：生产结构字段在 `SetMesh.jl`，构造在 `CzmMesh.jl`，克隆在 `CzmSolve.jl`，计算/返回/helper 在 `CouplingState.jl`，公开导出在 `JuBat.jl`；直接冻结该旧契约的测试为 `test_thermal_to_czm_interp.jl`、`test_czm_strain_inputs.jl`、`test_czm_phi_merge.jl`、`test_czm_multicycle_c4lite.jl`、`test_cohesive_struct.jl`、`test_create_czm_mesh.jl`。
- 实现已删除上述字段、构造、克隆、计算、helper 和 export；原插值专项测试随能力删除。`test_czm_strain_inputs.jl` 改为逐个力学 bulk 单元核对 `dT_czm[e] = mean(T_nodes[parent Q4]) - T0`，同时显式断言不再返回 `T_czm_nodes`。
- `build_thermal_to_czm_interp` 强依赖 `CzmSubmesh.mesh` 的结构化未合并节点编号（用 `mod(i-1,n_theta_nodes)` 定位父热单元），测试还冻结了“矩阵行数等于未合并节点数”。这正是保留第二套完整未合并力学 `Mesh` 的主要技术债务，而非当前物理需要。
- `create_czm_mesh` 注释称插值会按 `phi_keep` 裁剪，但实际 `build_thermal_to_czm_interp` 返回未合并 `.mesh.nlen` 行，`compute_czm_strain_inputs` 也显式要求该行数；注释、字段名和实现契约不一致。
- `node_map` 除构造、克隆外没有生产消费者；`interface_nodes` 被固定为 `[[]]` 后仅克隆；二者是兼容字段候选。
- `n_layers` 实际固定为 cohesive 本构类型数 2，却仍被示例输出为“层数”；这是高风险语义错误，不只是命名不美观。
- `clone_czm_mesh_with_damage` 为更新损伤手工浅拷贝全部拓扑字段，只深拷贝状态。拓扑与可变历史状态耦合在同一 `mutable struct`，每新增字段都可能漏拷贝，历史上已出现此类修复需求。

### 初步优化抓手

1. 将不可变 `CzmTopology` 与可变 `CzmState` 分离，先通过兼容属性/构造器维持旧接口。
2. 删除未使用的细力学节点温度插值、常驻矩阵和专用 helper；热—力温度耦合唯一采用父热单元平均值与 `thermal_elem_map`。
3. 让积分阶次成为最终 bulk 拓扑的显式元数据，消除仅为 `gs.order` 保留的 `bulk_mesh` 重复对象。
4. 为遗留字段建立“实际语义—兼容名称—退役条件”表，避免一次性删字段打断外部调用者。

## Phase 2 风险分类：编号空间与失效映射

### P0：`czm_element_map` 与真实 cohesive 编号不相容

- `JellyrollMesh.czm_element_map` 在粗热网格构造阶段由 `interface_pairs` 的相邻 Φ 段生成，其值域是 `1:(length(interface_pairs)-1)`，表示跨匝 Φ 段附近关系。
- `get_fractured_elements(czm_mesh)` 返回的是 `CohesiveMesh.cohesive_elements` 的真实单元 id；这些单元由每个周向段的四个箔–涂层面组成，值域与计数均不同。
- `CallModel_MultiSPMe` 和 `get_active_elements` 却把两者直接用 `czm_idx in fractured_czm` 比较。该比较依赖编号偶然重叠，不是有效拓扑映射。
- 当前有效工况探针 `D=0`，因此该路径通常休眠；一旦发生断裂，它可能停用错误的热/电化学单元。优化前必须把该字段从“可清理冗余”升级为“先冻结行为、再修复语义”的独立问题。
- 正确映射来源应是已存在的 `cohesive_to_thermal`，可一次构造反向 `thermal_to_cohesive::Vector{Vector{Int}}`；不能继续使用粗热 Φ 段映射冒充 cohesive 映射。

### P1：`MeshGeometry` 混入活动网格无效或休眠字段

- `setup_thermal2D_mesh` 默认选择 `thermal2D_merged`，同时把活动 `interface_pairs` 设为空；但仍把从未合并热网格生成的 `inner_nodes`、`outer_nodes` 和 `czm_element_map` 写入 `MeshGeometry`。
- `inner_nodes`/`outer_nodes` 当前没有 `case.geometry` 生产消费者；若未来直接用于合并活动网格，节点编号可能失效。
- `MeshGeometry.interface_pairs` 字段注释称“top_elem/bot_elem”，实际来源是节点对；默认合并路径中又为空。名称、注释、活动状态三者不一致。
- `element_layer`、`is_inner_layer`、`layer_weights`、`boundary_edges` 有真实消费者，应保留并加强尺寸不变量。

### 与用户范围的关系

- 用户此前要求暂不做损伤对电化学和热模型的影响。本审计不会顺手修改上述失效屏蔽业务行为；应在优化计划中先隔离/标记错误映射，再由用户单独决定是禁用既有消费者还是按 `cohesive_to_thermal` 修复。
- 拓扑重构不能把这一逻辑变化夹带进纯结构整理批次，否则行为基线变化无法归因。

## 兼容边界与既有测试

- `CzmSubmesh`、`CohesiveMesh`、`MeshGeometry`、`build_thermal_to_czm_interp`、`create_czm_mesh` 均已公开导出；其中插值 helper/字段原本属于公开契约。用户本次明确选择不保留细力学节点温度场，因此仅该插值契约获准退役，其余公开类型与构造器继续保持。
- `CzmSubmesh` 的历史 5 参数构造器有明确回归测试；旧调用者没有 Φ 合并信息时要求 `mesh_bonded === mesh` 且 `phi_keep=1:nnode`。
- `test_czm_submesh.jl` 和 `test_czm_phi_merge.jl` 冻结了当前双网格字段级契约，包括未合并真实 `phi_pairs`、合并节点计数和 `phi_keep`。
- `test_thermal_to_czm_interp.jl` 只冻结已被用户退役的插值 helper，应随生产能力删除；不再以测试反向保留无目标物理用途的字段。
- `CohesiveMesh()` 空构造器继续保留；只删除获准退役的 `thermal_to_czm` 属性，`cohesive_to_thermal` 及其他字段契约不变。
- `MeshGeometry` 是公开类型，并被热模型、CSV、示例和测试直接读取；结构字段变更需要兼容构造/访问层，而不能一次性重排。

## 可安全实施顺序的约束

1. **先加不变量和语义测试**，不改变字段与行为。
2. **再引入新内部拓扑核**，旧公开结构从新核构造，冻结输出逐位一致。
3. **迁移生产消费者**到新核/映射 API，但保留旧字段只读兼容。
4. **证明仓库与外部迁移完成后**，才在单独破坏性版本中考虑删字段。

任何“把两个 `Mesh` 字段合成一个”或“删除未使用矩阵”的直接修改，都会先打断现有公开契约；优化必须是迁移，不是原地删减。

## 对象规模快照（nθ=80，CZM 开启）

| 对象 | 节点/尺寸 | `Base.summarysize` |
|---|---:|---:|
| 粗热未合并 `thermal2D` | 3366 节点 / 1682 单元 | 1,130,866 B |
| 粗热合并 `thermal2D_merged` | 1763 节点 / 1682 单元 | 1,105,218 B |
| 力学未合并 `CzmSubmesh.mesh` | 15147 节点 | 8,854,722 B |
| 力学 Φ 粘结 `mesh_bonded` | 13544 节点 | 8,829,074 B |
| 最终 `CohesiveMesh` | 20276 节点 / 13456 bulk / 6728 cohesive | 22,742,498 B |
| `thermal_to_czm` | 15147 × 1763 | 983,680 B |

说明：`summarysize` 会递归计入对象引用，不能把各行直接相加作为峰值内存；但它足以说明两个完整力学 `Mesh` 各自约 8.8 MB，而当前不进残差的插值矩阵接近 1 MB。最优先收益来自避免重复 Gauss 几何/节点/连接，而不是压缩必要的最终 cohesive 拓扑。

实施后同一 nθ=80 构造探针得到 `CohesiveMesh` 21,758,810 B；相对审计快照 22,742,498 B 减少 983,688 B，与被删除插值矩阵的 983,680 B 量级一致。`thermal_elem_map` 仍有 13,456 项，证明 bulk 父热单元映射继续常驻。

冻结示例在当前脏工作区中成功完成，网格、19 步、电压、容量、温度、零损伤及 19 次 CZM 收敛均与基线一致；最大法向分离为 `1.5174e-12 m`（冻结值 `6.6820e-15 m`），PNG 为 92,736 B / SHA-256 `0946646a...`（冻结 92,775 B / `272402bb...`）。任务开始前已存在的 `src/czm.jl` 修改把 cohesive 局部法向改为 `host_inner_elem → host_outer_elem` 定向，并直接替换 basic 装配与缓存中的 `R`，会改变分离符号/微小量；该既有物理修改是严格基线漂移的直接生产路径。本批不修改或回退它。

用户随后明确接受当前结果并要求“不用重跑，直接更新基线”。旧基线档案的顶层 ID 仍为 2026-08-15，source manifest 仍只有 46 个 Julia 文件，且其脚本哈希与当前 `example/testexample.jl` 不一致；直接重冻结必须同步刷新身份、结果、PNG 和当前运行源码清单，而不能只改一项分离量。

进一步核对：当前仍是 45 个 `src/**/*.jl` 加入口脚本，共 46 文件；过期的是清单成员/内容哈希，而不是文件总数。标准 `source_manifest_tsv_sha256` 算法已用旧档案反证恢复：去掉表头和 aggregate 尾行后，将 `path<TAB>sha256` 数据行以 LF 连接、末尾不加换行，再取 SHA-256；可精确复现旧值 `77cc5df3...`。当前入口脚本 SHA 为 `344700a7...`，当前 HEAD 为 `a7fecc6b...`。

v4 基线已完成：ID `testexample-20260824T043411-0600`，分离量 `1.5174e-12 m`，PNG 92,736 B / SHA-256 `0946646a...`，46 文件清单聚合 SHA-256 `b5986a34...`。档案显式记录 `reuse_completed_run_no_rerun`，并用实际 PNG 写入时间标识来源运行；没有声称发生第二次仿真。

## 缓存一致性风险

- `CZMAssemblyCache` 以 `objectid(czm_mesh)` 判断拓扑失效；但 `CohesiveMesh`、内部 `Mesh` 及其 `node`/`element` 数组均可原位修改。
- 同一对象发生原位连接、节点、cohesive 元数据或局部方向变化时，`objectid` 不变，缓存的 `K_bulk`、DOF 映射、cohesive frame 和边界 DOF 可能全部过期。
- 目标架构应让拓扑不可变，或为拓扑构造内容签名/递增版本；缓存键应绑定 `topology_id`，而不是可变容器身份。

## 建议目标架构

建议保留多个**求解视图**，但只保留一个**逻辑拓扑真相源**：

```text
SpiralGridSpec（唯一几何/离散真相）
  ├─ theta / segment ids / layer offsets / material sequence
  ├─ Φ logical pairs（逻辑列号，不依赖某个 Mesh 节点号）
  └─ parent segment metadata
          │
          ├─ ThermalGridView
          │    ├─ active Mesh（当前为 merged）
          │    ├─ raw→active node map
          │    └─ layer weights / boundary cache
          │
          └─ MechanicalTopology（不可变）
               ├─ final node + bulk connectivity
               ├─ material / winding / n_theta metadata
               ├─ cohesive elements + oriented interfaces
               └─ CouplingMaps
                    ├─ bulk_to_thermal
                    ├─ cohesive_to_thermal
                    ├─ thermal_to_cohesive（反向索引）
                    └─ 无细力学节点温度场；温度只经 bulk_to_thermal 传递

MechanicalState（可变、与拓扑分离）
  ├─ u / plastic states / prestress reference
  └─ cohesive damage histories
```

### 核心原则

- `SpiralGridSpec` 保存逻辑列、层和配对，不保存多个带 Gauss 数据的完整 `Mesh` 副本。
- 最终力学装配只消费 `MechanicalTopology`，不再同时读取 `bulk_mesh`、`bulk_element` 和 `czm_submesh.mesh_bonded` 三个事实来源。
- 场映射集中在 `CouplingMaps`，字段名包含 source/target，例如 `bulk_to_thermal`；不得再用 `czm_element_map` 这种目标不明确的名称。
- 状态克隆只克隆 `MechanicalState`；拓扑通过不可变引用共享，避免手工复制十余字段。
- 公开旧类型在迁移期作为 facade，属性仍可读；生产代码禁止新增对旧歧义字段的依赖。

## 分阶段优化方案

### Batch G0：契约与风险固化（零行为变化）

- 新增拓扑签名工具和构造时不变量：节点/单元计数、正 Jacobian、父映射值域、cohesive 法向、Φ 配对完整性。
- 新增测试证明 `czm_element_map` 与真实 cohesive 编号空间不同；只标记/隔离风险，不在本批修改损伤消费者。
- 修正文档和注释中 `interface_pairs`、`n_layers`、插值裁剪等错误语义；旧字段名暂不改。
- 门禁：所有矩阵、残差、结果键和 `testexample` 必须逐位一致。

### Batch G1：内部拓扑核与显式访问器（零行为变化）

- 引入不可变 `MechanicalTopology` 和 `CouplingMaps`，由现有构造链一次生成。
- 为旧 `CohesiveMesh` 增加只读访问器；旧字段与构造器继续存在。
- 新旧路径对节点坐标、最终 bulk 连接、cohesive 列表、材料标签和全部映射逐位比较。
- 暂不移除任何旧网格，先建立迁移锚点。

### Batch G2：生产消费者迁移与状态分离

- 装配、边界、预应力、`core_ovalization` 和 strain input 改读新拓扑访问器。
- 将损伤/位移/塑性历史移入独立 `MechanicalState`；`clone_czm_mesh_with_damage` 演进为只克隆状态。
- 缓存键改为不可变 `topology_id + param_cache.id + BC signature`，覆盖原位修改风险。
- 门禁：随机位移下 K/R/分离逐位或严格容差一致；失败步提交/回滚测试全绿。

### Batch G3：映射语义修复（需单独用户授权）

- 用 `cohesive_to_thermal` 构造 `thermal_to_cohesive`，替代错误的 Φ 段 `czm_element_map` 消费。
- 由于会改变发生断裂后的电化学/热失效屏蔽行为，必须与用户“暂不做损伤反馈”范围分开决策。
- 可选决策为：暂时显式禁用既有失效消费者，或在理论闭合后按真实映射恢复；不得夹带在结构重构中。

### Batch G4：温度路径收敛与重复网格压缩

- 生产力学载荷继续只用 `bulk_to_thermal`；删除节点温度插值矩阵、返回字段和公开 helper。
- 后续可用逻辑 Φ 配对 + raw→bonded map 代替常驻未合并 `Mesh`；本窄批暂不实施该更大拓扑变更。
- 在外部迁移期保留 facade；确认弃用窗口后才删除 `node_map`、`interface_nodes`、歧义 `n_layers` 等字段。
- 目标收益：nθ=80 时优先消除约 8.8 MB 的重复力学 `Mesh` 常驻数据和约 0.98 MB 的未使用插值矩阵，同时降低 Gauss 缓存重复构造。

### Batch G5：热网格返回对象收敛

- `JellyrollMesh` 内部改为 `SpiralGridSpec + active thermal view`；未合并热视图仅在兼容访问或诊断时生成。
- `MeshGeometry` 收敛为活动网格有效元数据，未合并节点 id 不再混入活动 merged 网格。
- 这是公开 API 影响最大的最后批次，应晚于所有生产消费者迁移，并提供明确弃用周期。

## 每批验证门

- 公开构造器/属性兼容测试；
- 拓扑计数、坐标、连接、材料、父子映射和法向签名；
- 单元 patch、刚体运动、自由膨胀和有限差分切线；
- CZM basic/load-substep/arc 专项与状态回滚；
- `test/runtests.jl` 全量；
- `example/testexample.jl` 网格/步数/科学指标/PNG SHA-256；
- 构造时间、峰值内存和 `summarysize` 前后对比；
- 任何科学输出漂移立即停止，不以“数值接近”放行。

## 待核实项

- 各结构字段定义、构造位置和消费位置。
- 合并/未合并网格的真实必要性。
- 热节点插值与热单元父子映射是否存在重复来源。
- 缓存失效条件是否覆盖拓扑及映射变化。
- 测试是否覆盖计数、方向、映射守恒和旧构造器兼容。

## 技术决策

| 决策 | 原因 |
|---|---|
| 不把所有网格直接合并成一个对象 | 热连续场、cohesive 双侧自由度和 Φ 粘结求解视图具有不同拓扑需求。 |
| 优先减少“双重真相” | 重复事实比对象数量本身更容易产生索引漂移和缓存错误。 |

## 资源

- `src/Jellyrollmodel.jl`
- `src/CzmMesh.jl`
- `src/SetCase.jl`
- `src/CouplingState.jl`
- `src/czm.jl`
- `md/02_几何与网格.md`
- `md/10_参数传递与模块架构.md`
