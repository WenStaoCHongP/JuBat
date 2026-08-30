# CZM 源码三文件拆分发现

## Requirements

- 以 `md/源码函数索引/czm.md` 和实际源码为依据。
- 把 `src/czm.jl` 拆为：
  - `src/CzmMesh.jl`
  - `src/Czm.jl`
  - `src/CzmBC.jl`
- 保持行为、API 与全耦合基线不变。

## Research Findings

- `src/czm.jl` 当前为 873 行，实际包含 2 个 struct 和 11 个函数；源码索引与定义数量一致。
- 网格职责：`CohesiveElement`、`create_czm_mesh`。
- 核心状态/装配职责：`DamageState`、`moduli_of`、五类装配函数及两个缓存函数。
- 边界职责：`apply_bc_czm`、`identify_bc_nodes_czm`。
- `build_czm_cache` 调用 `identify_bc_nodes_czm`，而 `create_czm_mesh` 初始化 `DamageState`。为了让每个跨文件符号在定义前已经加载，模块入口采用 `CzmBC.jl → Czm.jl → CzmMesh.jl`。
- 当前索引中的职责和行号需要按三个新源文件分别重建，而不是保留指向已删除 `src/czm.jl` 的单页索引。
- 当前活动文档中的旧路径包括 `AGENTS.md`、`md/03_边界条件.md`、`md/06_内聚力模型_CZM.md`、`md/07_界面热阻模型.md`、`md/14_粘性正则化.md` 以及 `md/对照/` 下两份对照文档；拆分后应按函数归属更新。
- 外部调用者通过 JuBat 顶层函数名使用这些 API，未发现依赖源文件局部 include 的测试或示例，因此公开名称保持不变即可兼容。
- 拆分后的实际规模为：`CzmMesh.jl` 182 行（1 struct + 1 函数）、`Czm.jl` 629 行（1 struct + 8 函数）、`CzmBC.jl` 63 行（2 函数）。
- Windows 文件系统大小写不敏感，但目录项已用大小写精确比较确认：旧名 `czm.jl` 不存在，新名 `Czm.jl` 存在。
- Git 的 `core.ignorecase=true` 仍把大小写变更显示为原跟踪路径的修改；最终需用两步 `git mv` 经临时名记录 `czm.jl → Czm.jl`，同理处理函数索引 `czm.md → Czm.md`。
- `md/源码函数索引/_索引.md` 目前仍只列 `czm.md`，并在统计表中把它作为单文件统计；需要拆成三个条目并保持总定义数为 2 struct + 11 function。
- `md/源码函数索引/JuBat.md` 需要把入口行数从 91 更新为 93，并把原单一 `czm.jl` include 改为三个新文件。
- `md/06_内聚力模型_CZM.md` 的代码位置表还列出实际不存在的 `cohesive_element_matrices`；已在本次路径更新中改为真实入口 `assemble_czm_system`。
- 仍含 `src/czm.jl` 的大量 `docs/superpowers/specs/`、旧 plans、`docs/comparison/` 和 `Simplify/baseline/source_manifest.tsv` 是历史设计/基线证据，不应伪装成当前路径而机械改写；当前入口、技术文档、源码索引和活动简化计划应更新。

## Technical Decisions

| Decision | Rationale |
|---|---|
| 先按定义间依赖而非单纯行号切割 | Julia include 顺序要求被引用的类型和函数先定义 |
| `CohesiveElement` 与 `create_czm_mesh` 归 `CzmMesh.jl` | 二者共同定义 CZM 网格拓扑及其构造 |
| `DamageState` 归 `Czm.jl` | 它是本构历史状态，不是网格结构 |
| 两个边界函数归 `CzmBC.jl` | 形成单一边界识别与施加职责 |
| include 顺序使用 `CzmBC.jl`、`Czm.jl`、`CzmMesh.jl` | BC 无本组三文件依赖；核心缓存依赖 BC；网格构造依赖核心的 `DamageState` |

## Issues Encountered

| Issue | Resolution |
|---|---|
| 暂无 | — |

## Resources

- `src/czm.jl`
- `md/源码函数索引/czm.md`
- `src/JuBat.jl`
