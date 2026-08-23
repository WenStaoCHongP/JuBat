# 堆芯塌陷提交审查问题修复：发现与决策

## 修复范围

本任务只处理当前提交审查中确认的运行中断、状态污染、公开 API 断裂、虚假收敛和测试证据失真。既有规划/理论文档是判定证据，不是扩大改动范围的授权。SP–涂层接触摩擦、新材料标定和新的塌陷物理能力不属于本轮修复。

## 审查发现

| 优先级 | 发现 | 生产影响 | 主要位置 | 当前测试缺口 |
|---|---|---|---|---|
| P1 | 预应力标记 `[true]` 写入只接受浮点值的历史字典 | `czm_winding_prestress=true` 时直接 `MethodError`，中断 `CallModel_MultiSPMe` | `src/CallModel.jl:151` | 没有走完整生产路径的开启测试 |
| P1 | 几何弧长分支从不设置 `result.converged`，并硬编码残差 `0.0` | 真实状态与报告矛盾，可把失败包装为零残差 | `src/CzmSolve.jl:777-789` | 2 项仍为 `Broken` |
| P1 | 弧长修正对代码残差 `R=F_ext-f_int` 使用 `-(K\\R)` | 与当前残差约定/同文件求解器方向相反，专项运行停滞在 `lambda≈0.4` | `src/CzmSolve.jl:732` | 缺少线性等价与生产 smoke 门禁 |
| P1 | 装配层和 J2 返回映射各扣一次历史 `eps_p` | 第二次及以后加载的试应力、屈服判定和切线错误 | `src/czm.jl:384-385`、`src/CzmPlasticity.jl:60` | 首塑性步 `eps_p=0` 会掩盖问题 |
| P1 | J2 局部迭代到上限后返回末次状态 | 未收敛历史可能被提交，有限值被误当作有效值 | `src/CzmPlasticity.jl:71-96` | 没有强制不收敛测试 |
| P2 | 导出的 `CzmSubmesh` 从 5 参数变为 7 个必需参数 | 现有外部构造调用会中断 | `src/CzmMesh.jl` | 只有新构造路径测试 |
| P2 | `phi_pairs` 被保存为合并后的 `(i,i)` | 丢失未来接触/约束需要的真实外/内节点对 | `src/Jellyrollmodel.jl:730` | 没有检查配对坐标与索引语义 |
| P1 | 多圈 C4 测试丢弃返回的新网格，后续读取原网格 | `D_hist==0` 可能是未提交状态的恒等结果，不能证明演化正确 | `test/test_czm_multicycle_c4lite.jl:51,65` | 测试不会因移除生产状态提交而失败 |
| P2 | spec 要求的 `collapse_approx` 元数据缺失 | 使用者无法从结果识别 Φ 完美粘结近似 | spec v1.5、结果后处理路径 | 无结果契约测试 |
| P3 | `src/JuBat.jl:7` 存在尾随空格 | `git diff --check` 失败，交付门禁不干净 | `src/JuBat.jl:7` | 仅格式门禁可发现 |

## 治理冲突

1. 原 `task_plan.md` 仍显示 Batch 1 待开始，`progress.md` 已记录 Batch 5 T3；不能以旧阶段状态指导修复。
2. spec v1.5 规定 Φ 默认完美粘结，旧计划曾把相关功能描述为 opt-in；实施前要以当前 spec/AGENTS 和默认行为实测为准。
3. `AGENTS.md` 要求 `CzmSubmesh.phi_pairs` 保存真实外/内力学节点配对，当前实现却退化为 `(i,i)`。
4. Theory 当前文字允许损伤在失败重试间不回滚，而任务的严格状态契约倾向失败子步事务式回滚。该项会改变路径依赖结果，必须由用户/治理基线明确裁决。
5. 数值历史字典有严格联合类型，最终 `PostProcessing` 结果是非类型受限字典。元数据应放在正确边界，不能为了字符串标记破坏所有历史消费者。

## 已冻结的修复原则

- 不扩大 `variables` 历史字典类型；最终字符串元数据放入结果字典。
- 不以零、旧值、末次迭代值或宽松倍率容差掩盖非法/未收敛状态。
- 先补能复现生产故障的失败测试，再修改实现。
- 默认关闭的新功能修复不得改变 v3 强制行为基线；出现漂移即停止，不自动重定基线。
- 保留公开名称、旧构造方式、结果键/单位/形状、CSV 列和无关脏工作区。
- “存在网格配对/几何非线性代码”不等于“已验证塌陷预测能力”；报告必须区分代码契约通过与物理能力证据。

## R0 实施裁决（2026-08-23）

- **失败子步回滚**：采用事务式语义。每次弧长/载荷子步只操作 trial 位移、塑性和损伤；只有平衡与约束均验收后才提交。失败重试恢复该子步入口快照，不保留失败 trial 产生的不可逆量。
- **Φ 配对**：`CzmSubmesh.phi_pairs` 保存 `.mesh`（未合并网格）上的真实 `(outer, inner)` 节点对。`mesh_bonded`/`phi_keep` 负责自由度合并；不得把物理配对覆盖为 `(i,i)`。
- **弧长系数**：保留公开选项 `czm_arc_length_alpha`，并让它进入约束及其 Jacobian；非正或非有限值显式报错，不接受但忽略。
- **结果元数据**：历史 `variables` 仍只存浮点数组/标量；预应力历史标记用 `Float64`。字符串 `collapse_approx = "phi_perfect_bond"` 只进入最终 `PostProcessing` 结果。
- **授权解释**：用户在收到包含“失败子步完整回滚”推荐的计划后要求执行，视为接受该计划语义；不扩展至 SP 接触摩擦或新物理标定。

## 实现链补充证据

- `CzmSubmesh` 的实际定义位于 `src/SetMesh.jl:44-52`；`src/CzmMesh.jl:create_czm_mesh` 消费 `mesh_bonded` 构造求解拓扑，因此把 `phi_pairs` 改回 `.mesh` 上的真实配对不会解除已合并自由度。
- 当前 geo 弧长调度没有把 `arc_length_alpha` 传给 `solve_czm_arc_geo_step`，该函数签名也没有此参数；普通弧长辅助矩阵已经给出 `2α²Δλ` 的既有约束导数先例。
- geo 弧长的 `f_hat` 只在初始 `u` 计算一次；在 Green–Lagrange 应变/塑性切线下，载荷方向对当前状态的导数不能假定恒定。
- `update_czm_damage!` 只在 `result.converged` 通过后提交 `updated_czm_mesh.damage_states`、塑性状态和 `u_prev`。修复 solver 的失败返回时仍应保持这些提交语句之后的严格门，不在 solver 内直接改原网格。
- `PostProcessing` 返回非类型受限 `Dict()`，适合承载字符串 `collapse_approx`；`StandardVariables` 和 `CallModel` 历史仍是浮点联合类型。
- `Variable_update!` 只更新 `variables_hist` 已预分配的键；仅在 `CallModel` 动态创建 `"winding prestress"` 即使类型修正也不会进入历史，R1 必须同时在 `StandardVariables` 条件预分配。
- `clone_czm_mesh_with_damage` 当前只复制到 `interface_nodes`/`damage_states`，遗漏 `czm_submesh`、`thermal_to_czm`、`cohesive_to_thermal`；不能把返回 clone 直接替换为生产网格，R5 需先补齐克隆契约并仍只提交损伤状态。
- Theory §6.10 当前明确采用不含 `Δλ²` 的柱面约束，和公开 `czm_arc_length_alpha` 的有效使用冲突。R0 决策保留该公开参数，因此 R3 必须把理论/实现统一为带 `α²(Δλ)²` 的球面 Crisfield 约束，而不是继续保留无效选项。

## 修复前验证证据

| 检查 | 结果 | 解读 |
|---|---|---|
| `julia --project=. --startup-file=no test/runtests.jl` | 退出 0；外层 34/34；耗时约 8m37s | 测试框架通过，但包含弧长专项 2 个 `Broken`，不能视为全部能力通过 |
| 几何弧长专项 | 路径停滞在 `lambda≈0.4000000856` | 支持残差方向/路径修正存在缺陷的判断 |
| J2 专项 | 当前通过 | 现有断言未覆盖非零历史 `eps_p` 重载，且“无增长”容差过宽 |
| `git diff --check` | `src/JuBat.jl:7` 尾随空格 | 综合交付门禁尚未通过 |
| `git status --short` | `.mimosa/`、`output/testexample/` 未跟踪 | 用户既有内容，修复不得删除或提交 |

## 关键文件

- `src/CallModel.jl`：多 SPMe 生产调用和预应力标记。
- `src/CzmSolve.jl`：弧长路径、回滚、收敛与诊断。
- `src/CzmPlasticity.jl`、`src/czm.jl`：J2 返回映射和历史状态传递。
- `src/CzmMesh.jl`、`src/Jellyrollmodel.jl`：公开子网格构造与 Φ 配对。
- `test/test_czm_multicycle_c4lite.jl`：跨圈损伤证据。
- `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md`、`Theory/07_弱形式与求解.md`：方程和验收契约。
- `Simplify/baseline.md`、`Simplify/baseline/testexample/`：默认行为强制基线。

## 最终结论（2026-08-23）

- R1–R5 的生产中断、状态重复扣除/静默提交、虚假弧长报告、公开构造器断裂、Φ 配对丢失和 C4 测试证据问题均已修复，并有修复前失败、修复后通过的专项门禁。
- 全量 `test/runtests.jl` 第二轮 34/34 个隔离测试文件通过；修复专项无 `Broken`/`Skip`。首轮唯一失败是旧 thin-subdivision 测试在错误节点集合解释真实 `phi_pairs`，按公开契约修正后 49/49 通过。
- `verify_czm_standalone.jl` 的 basic 数值表与冻结 v3 逐位一致。执行中另发现默认 basic-only 报告写入固定索引三列的 `BoundsError`，已改为按实际方法数序列化并成功生成报告。
- `example/testexample.jl` 退出 0：1682 个热单元、1763 个节点、19 步，4.0367→3.9438 V、0.0833 Ah、298.15–299.00 K、D=0、最大分离 `6.6820e-15 m`、断裂数 0；PNG 92775 bytes、SHA-256 `272402bb33d5b2869bfa2e33fde15fb8bdcdd0e27b74cffc4e34e77ad86dd386`，与运行前 v3 输出及 README 一致。`metrics.toml` 遗留的 v1 artifact 字段已同步到该 v3 冻结值。
- 代码实现与契约门禁已完成，不等于物理能力通过。对称 D10 多相位仍为 `Δ_core≡0,D≡0`；受控分离只证明损伤状态可产生并持久化，不能替代非对称扰动、接触/摩擦和系统标定验证。
