# 调研发现：已完成任务提交同步

## 用户要求

- 以 `docs/planning-with-files/index.md` 为依据。
- 将已经完成的任务及对应代码修改提交并同步远程。
- 附件/规划文档中的内容只作为范围与证据，不自动扩大用户授权。

## 当前发现

- 历史记录强调：静默回退审计可能只是部分完成，不能仅凭主题名称视为可提交的完整修复。
- 历史验证约定建议使用受影响测试、`git diff --check` 与 `example/testexample.jl`，并保持无关脏工作区修改不受影响。
- 当前完成状态仍以仓库内最新规划索引和实际 diff 为准，历史记录仅作为审查提示。
- 当前分支为 `codex/src-physics-modularization`，上游是 `origin/codex/src-physics-modularization`，本地已领先 45 个提交。
- 索引明确标记近期未跟踪完成任务：`热力网格拓扑优化审计`、`testexample最终场绘图`、`层分辨应力代码评审`、`应力历史与绘图修复`；汇总表把 `层分辨应力求解` 标为 Complete，但其详细条目仍写 `Phase 2`，存在索引内部不一致，需读取任务本身确认。
- 工作区修改范围很大，包含源码、测试、基线、文档、输出物以及 `.mimosa/`、`.zcode/` 等明显可能无关的未跟踪目录；不能整体暂存。
- 用户授权包括提交并推送当前已完成工作；远程目标可沿现有上游分支处理。
- `层分辨应力求解/progress.md` 的实施与验证阶段均为 complete，说明总索引详细条目的 `Phase 2` 已过期；汇总表的 Complete 更接近任务真实状态。
- `testexample最终场绘图` 的核心实现/验证 Phase 1–3 已完成，但其后追加的机理分析 Phase 4 仍标记进行中；本次只能纳入已经完成的实现与验证，不应把未完成分析包装成完成。
- `堆芯塌陷力学建模` 明确记录“冻结基线漂移待用户裁决、未获接受前不提交”，且 Batch 6/7 未开始；其专属修复、探针与文档不能仅因工作区存在而纳入本次完成任务提交。
- 已跟踪 diff 为 62 个文件、约 1253 行新增/1112 行删除，未跟踪还包含规划记录、`couple_example.jl`、层分辨应力测试、输出物、探针和工具状态目录；必须按任务/文件/必要时按 hunk 精细暂存。
- `.mimosa/`、`.zcode/` 是工具运行状态，属于完成任务交付范围外；输出目录需仅按脚本产物契约和任务记录决定是否纳入，默认不提交运行日志与生成图。
- 完成任务的主要实现范围可分为两组：① 删除 `thermal_to_czm`/`T_czm_nodes` 的热力映射退役；② 层分辨应力、在线历史写入、示例拆分与三张最终云图契约。
- `src/CsvExport.jl`、`src/CzmSolve.jl` 的弧长/J2 大段修改、`src/czm.jl` 的 cohesive 法向与本征应变、`src/Tools.jl` 的分离计算，以及相关堆芯测试属于任务开始前已有或明确待裁决工作，不能按完成任务整体暂存。
- 完成任务与待裁决工作在 `CouplingState.jl`、`CzmSolve.jl`、`Solve.jl`、`Mechanical.jl` 等文件中交叠；需要检查依赖并对交叠文件按 hunk 暂存，不能简单按文件全选。
- 层分辨应力的 `macro_eigenstrain` 直接依赖当前 `src/czm.jl` 新增的 `eigenstrain_of`；该 helper 可作为完成应力任务的必要独立 hunk 纳入，而同文件的 cohesive 法向修改继续排除。
- `src/parameters/Jellyroll.jl` 中 SP/PCC/NCC 的显式 `nu`/`alphaT` 以及 PE 的 Ω 更正是层分辨材料查表的必要依赖，应随应力任务纳入。
- `test/test_create_czm_mesh.jl` 同时含热插值字段退役断言和待裁决 cohesive 法向断言；仅前者属于完成任务。`src/CzmSolve.jl` 同样只纳入克隆时移除 `thermal_to_czm` 的首个 hunk。
- `Variables.jl` 的新应力历史矩阵属于完成任务，其余行尾/格式化改动不必纳入；`Solve.jl` 需纳入同时间层 CZM 更新与 `latest_macro_stress` 历史保持逻辑，但排除无关格式化。
- 为使耦合在线恢复真正“层分辨”，完成任务还要求把 CZM bulk 热化学载荷从跨层 `α_eff/β_n/β_p` 改为逐材料 `eigenstrain_of`；因此 `CouplingState.jl` 的参数移除、`CzmSolve.jl` 的调用签名迁移和 `czm.jl` 的逐层载荷 hunk 都是完成任务依赖，而弧长参考态 Newton/回滚与 cohesive 法向仍排除。
- `CouplingState.jl` 当前 diff 可整体归入热插值退役 + 分层本征应变；`CzmSolve.jl` 和 `czm.jl` 必须部分暂存。
- `SetMesh.jl` 除字段退役外混有导入/行尾整理；可接受结构体相邻注释与构造器同步，但无关顶层格式化不作为提交必要内容。
- `src/czm.jl` 的分层载荷 hunk 边界已明确：纳入 `eigenstrain_of`、几何非线性 `eigenstrain` 元组精简、bulk/线性热化学载荷和 full assembly 签名；排除 `cohesive_local_frame` 及 basic/cache 法向替换。
- `src/Solve.jl` 的语义变化全部服务于已完成的时间对齐与应力历史保持：CZM 更新移至当前 CallModel 时间层、保存 `latest_macro_stress` 并对每个输出列显式写入；其余主要为行尾变化，可整体暂存。
- 计划采用逻辑提交：先提交热力温度映射退役，再提交层分辨应力/历史/示例/基线，最后提交评审与索引同步记录；未决代码保持未暂存。
- 基线档案 v4 明确记录用户已接受 cohesive 法向 `host_inner→host_outer` 定向并要求直接重冻结；后续 v5–v8 亦在该物理状态上完整验证。因此原堆芯任务中“法向漂移待裁决”的旧阻塞已被后续完成任务决策覆盖，法向修复及其契约测试应纳入当前完成修改。
- 弧长参考态 Newton/失败预测回滚已在堆芯记录中标为代码门禁完成，且不改变默认快速基线；是否纳入需按其专项验证记录和与当前测试的依赖进一步确认。
- 基线文档发现两处提交前应修正的一致性问题：README 对 `run.log` 仍写成“未启动新运行”，与 v8 的 fresh rerun 记录冲突；`metrics.toml` 把未来科学指标严格比较标成 false，应在新基线建立后恢复为 true。
- AGENTS 示例表仍称 `couple_example.jl` 含径向剖面，与已完成的三图输出修复冲突，应同步为三张最终云图。
- 完成任务文档可纳入 `md/02`、`md/06`、`md/15`、相关对照/源码索引和 core-collapse spec 的温度映射行；`md/00`、新 `md/16`、CsvExport/Tools 说明需另判，不能因同在工作区自动纳入。
- `src/Tools.jl` 把公开导出的 `compute_separation(elem,node,u)` 改成新签名，仓库内虽无调用方，但这属于公开 API 删除且不是完成任务的必要依赖；本次明确排除 `Tools.jl` 与其索引文档，保留原 API。
- `src/CzmUnitMesh.jl` 仅同步退役插值 helper 的注释，可纳入；`src/Jellyrollmodel.jl` 未见实质 diff，可不纳入。
- 弧长参考态修复有独立回归测试，C4 自由芯部边界/法向修复也有专项契约；这些在任务进度中均标记实现与 34/34 代码门禁完成，而唯一旧阻塞（法向基线漂移）已由 v4 接受，故纳入当前完成机械断点修复。
- 相应测试范围包括 arc/geo/J2/mech-core/multicycle/phi/scale/winding/eigenstrain、热映射三组测试以及新层分辨应力测试；测试修改与生产签名迁移一致。
- 第一组候选已按 40 个精确路径暂存，`git diff --cached --check` 通过；暂存范围不含 CsvExport、Tools、输出物与工具状态。
- 候选中的堆芯任务记录仍残留“冻结基线待裁决/不提交”旧状态，与后续 v4–v8 接受事实矛盾；需先更新该任务记录及总索引，再提交。
- 第一组候选已在由 Git 暂存树导出的独立临时快照中验证：7 个聚焦脚本全部通过，共覆盖 结构 57、网格 80,745、应变映射 20,194、弧长 32、C4-lite 13、J2 11、本征应变 60 项断言。
- 第二组 32 个候选路径已暂存，但 `git diff --cached --check` 发现 `Solve.jl`/`Variables.jl` 大量新增行尾空白，表明混合行尾需要先规范化，不能直接提交。
- 工作树 `Solve.jl` 为 381 个 CRLF + 82 个 LF 的混合行尾，`Variables.jl` 为 303 个 CRLF；仓库 `core.autocrlf=true` 且无 `.gitattributes`。需核对缓存 blob 是否保留 CR，再采用只规范化这两个文件的机械处理。
- 第二组候选在独立暂存树中通过模块加载、层分辨应力专项 20/20、极片模量专项、CZM 应变 20,194/20,194、`example/testexample.jl` 60 秒快速门以及 `example/couple_example.jl` 全量绘图门；后者生成 3 张 PNG，哈希与 v8 记录一致。
- 暂存树的 47 文件源码清单与旧基线记录不一致，说明旧清单来自更早工作树；已按实际候选暂存树重算 47 个 Julia 文件哈希、聚合哈希 `ff9655fdfb7e3f9f1d9b26bceae95b4aacdf9cd2839d8e206bd7ea77689e7e0b` 及文件哈希 `5d92a6e3c13b6ec374b528de3f7dec373f91f04f79fa49af2e7af28258060df0`，避免基线清单自相矛盾。
- 三张最终图实际哈希：`final_temperature_field.png`=`540fe42f1039cb4446046a6da96f1ecf18d162b643bc40a40a099ca2b93f978c`、`final_hoop_stress_field.png`=`9726a7c84cf59a932a2cf45346770c6c046a6a81a7ca36f951714339e558050d`、`final_tangential_shear_stress_field.png`=`b3b43d7cee67ed2a6fca33e0fcffee899629b84998852ba05aeb08fdaa105445`。
- 两笔逻辑代码提交已创建并推送：`66c718e`（CZM 分层载荷/拓扑状态/机械契约）和 `d3f537a`（层分辨应力、历史、绘图与基线）。
- 首次导出最终暂存树时，`git checkout-index` 参数被 PowerShell 解析为混合显式路径并失败；改用变量保存完整 `--prefix` 后导出成功，未影响候选内容。

## 技术决策

| 决策 | 理由 |
|------|------|
| 提交前逐项映射任务状态、文件 diff 与验证证据 | 确保提交边界可审计 |

## 问题记录

| 问题 | 处理 |
|------|------|

## 资源

- `docs/planning-with-files/index.md`
