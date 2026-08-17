# Task Plan: Review JuBat src-simplification specification and plans

## Goal
在完成 simplification spec/plans 评审后，建立可复现的 `example/testexample.jl` 代码简化前基线，并要求后续每个修改批次与该基线结果一致。

## Current Phase
Phase 10（安全可执行范围完成；D1/D2 待重写）

## Phases

### Phase 1: Scope and document discovery
- [x] Read the governing anti-bloat plan and simplification spec
- [x] Inventory all plans governed by the spec
- [x] Establish review criteria and traceability map
- **Status:** complete

### Phase 2: Cross-document review
- [x] Check scope coverage, omissions, and contradictions
- [x] Check sequencing, dependencies, rollback, and verification
- [x] Compare material claims against the current source tree
- **Status:** complete

### Phase 3: Findings synthesis
- [x] Rank findings by severity and cite exact file/line evidence
- [x] Separate blocking issues from improvements and open questions
- [x] Record recommended spec/plan corrections
- **Status:** complete

### Phase 4: Verification and delivery
- [x] Re-check every reviewed plan against the baseline criteria
- [x] Ensure review conclusions are internally consistent
- [x] Deliver concise Chinese review report
- **Status:** complete

### Phase 5: testexample 基线冻结
- [x] 确认 Julia 运行时与脚本依赖
- [x] 记录源码版本、工作树状态和运行参数
- [x] 运行 `example/testexample.jl`
- [x] 保存控制台日志和关键数值结果
- [x] 定义后续一致性比较规则
- **Status:** complete

### Phase 6: 执行前校正与首个安全简化批次
- [x] 恢复执行上下文并确认技能/工具可用性
- [x] 根据评审结论跳过当前不安全的 D1/D2 方案
- [x] 固化静态规模指标与首批变更范围
- [x] 为 ThermalDistributed D3 双胞胎建立先验行为测试
- [x] 将 non-`!` API 改为原位实现的薄包装
- [x] 运行局部测试和 `example/testexample.jl` 基线回归
- [x] 记录变更量、结果与后续批次决策
- **Status:** complete

### Phase 7: 下一安全批次甄选
- [x] 重新核验低风险候选的真实调用者与公共 API
- [x] 完成 SetCase、参数集、ring 与 install 的保留审计
- **Status:** complete

### Phase 8: CouplingState 旧兼容入口清理
- [x] 核验 `update_czm_damage!` 两个方法及调用者
- [x] 删除无内部调用、未导出的六参数兼容方法
- [x] 运行 CZM 定向测试与强制全基线
- [x] 更新计划状态和批次指标
- **Status:** complete

### Phase 9: Solve 兜底路径审计
- [x] 核验 try/catch 的职责、调用路径与可观察行为
- [x] 仅处理可证明的 silent swallow，不改 CZM 数值恢复路径
- [x] 运行定向测试与强制全基线
- **Status:** complete

### Phase 10: 剩余计划收敛
- [x] 汇总剩余 Pending 计划并区分 audit、safe action、blocked action
- [x] 完成无源码变更的审计计划
- [x] 执行仍可由测试充分覆盖的安全简化 concern
- [x] 形成已完成/保留/阻塞的最终状态表
- **Status:** complete

### Phase 11: CsvExport 重复容错收敛
- [x] 用定向测试锁定成功/失败记账与 warning 行为
- [x] 将七处同构 try/catch 收敛为单一内部 helper
- [x] 运行定向测试、模块测试与强制全基线
- **Status:** complete

## Key Questions
1. Do the spec and plans faithfully implement the baseline plan's principles and scope?
2. Are the plans independently executable, correctly ordered, measurable, and regression-safe?
3. Would following them actually reduce code volume and conceptual duplication without breaking JuBat behavior?

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Treat `Simplify/reducing-ai-code-bloat-plan.md` as the governing baseline | The user explicitly requested review according to that plan |
| Review documentation only unless source inspection is needed for claim verification | The request is a review, not authorization to implement changes |
| Store planning-with-files artifacts under `docs/planning-with-files/<中文任务名>/` | User-defined project convention overrides the skill's generic project-root default |
| Use `example/testexample.jl` as the mandatory behavioral baseline for simplification batches | Explicit user requirement; no later cleanup is accepted if this baseline changes |
| Start implementation with ThermalDistributed D3 instead of D1/D2 | The review found D1/D2 plans behavior-changing/incomplete; D3 is a bounded duplicate-removal concern with a direct characterization test and unchanged public API |
| Do not use the unavailable `simplify` skill | It is not installed/exposed in this session; follow the governing Reduce-mode plan directly and retain planning-with-files tracking |
| Retain documented ChooseCell parameter options even without repository callers | They are public string API; repository grep cannot prove absence of external callers |
| Remove only the six-argument `update_czm_damage!` compatibility method in CouplingState task 1 | It is unexported, has zero callers, and only adapts obsolete externally supplied state to the active three-argument method; the active path and numerical core stay untouched |
| Retain `assemble_coupled_system_full` | It is explicitly exported; repository grep alone cannot authorize breaking an exposed CZM assembly API |
| Retain the Mechanical zero-displacement fallback this round | The only prescribed characterization is a 3600 s / nθ=360 postprocess example and the frozen baseline does not assert this failure path; deleting it without focused coverage would violate the safety rule |
| Implement CsvExport guard consolidation | Seven blocks are structurally identical, the behavior can be pinned with a small focused success/failure test, and the helper remains private |
| Replace the two thermal public helpers with their in-place variant implementations and delete the `!` variants | Explicit user direction supersedes the wrapper design; internal callers already own copied K/F matrices, and the full baseline will guard numerical identity |

## Errors Encountered

- 2026-08-05：基线档案验证时，Julia `-e` 表达式中的 Windows 路径引号被 PowerShell 剥离，出现 `UndefVarError: rawSimplify not defined`；改用通过 `ARGS[1]` 传递路径或直接使用 PowerShell 校验，避免重复该命令。
- 2026-08-05：首次从 `source_manifest.tsv` 反算聚合哈希得到 `cb0597...`，与预检值不一致；初步判断为 TSV 行筛选/换行重建方式与原始算法不同，需按生成时的“排序后的源文件路径 + TAB + 小写 SHA256，以 LF 拼接”算法直接复核。
- 2026-08-05：复核发现当前 Windows PowerShell/.NET 不支持 `Path.GetRelativePath`，且初版 `source_manifest.tsv` 将制表符误写成字面量 `` `t``；改用工作区绝对路径前缀截断生成相对路径，并通过补丁修正 TSV 分隔符。
- 2026-08-05：批量更新 plan 状态时误判 `CallModel.md` 的状态行格式，导致整批补丁验证失败且未产生部分写入；改为先读取每个文件的精确状态行，再分组应用补丁。
- 2026-08-05：剩余计划核验命令引用了不存在的 `example/czm_cycle_example.jl`（计划中的陈旧路径）；实际 CSV 示例为 `example/coupled_czm_thermal_example.jl`。后续不再使用陈旧路径。
- 2026-08-05：D3 原位语义测试错误假设 Jellyroll 默认侧面对流 Bi 非零，导致 `nnz(K)>0` 与 `F` 非零两项失败；当前参数的该边界贡献合法为零。改为验证对象同一性、有限值、变体符号已删除，并由全基线验证数值。
- 2026-08-05：最终全工作树 `git diff --check` 命中本任务范围外的既有尾随空格（另一 planning 目录与 `src/PostProcessing.jl`）；不修改用户既有内容，改为仅检查本轮源文件与简化计划目录，并单独检查新增未跟踪文件。
| Error | Attempt | Resolution |
|-------|---------|------------|
| Initial recursive PowerShell table omitted paths because of display formatting | 1 | Use `rg --files` and explicit string output for reliable inventory |
| Default `Get-Content` decoding rendered the Chinese spec as mojibake | 1 | Re-read all relevant documents explicitly as UTF-8 |
| `apply_patch` rejected empty move hunks for two planning files | 1 | Retry with a minimal log/resource update in every moved file |
| `Pkg.status(["Plots"])` version probe lost quoting and evaluated `Plots` as an undefined variable | 1 | Use `using Plots; Base.pkgversion(Plots)` instead |
| Direct Plots version probe still lost embedded string quotes through PowerShell | 2 | Remove all string literals and use Julia `@show` expressions |

## Notes
- Preserve existing project implementation and unrelated user changes.
- Cite absolute file paths and exact line numbers in the final review.
