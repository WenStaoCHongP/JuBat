# 堆芯塌陷提交审查问题修复：进度记录

## 2026-08-23 — 修复计划建立

### 当前阶段

计划完成，等待用户批准实施；源码和测试尚未修改。

### 已完成

- [x] 复核当前提交审查中的生产中断、状态、API 和测试证据问题。
- [x] 对照当前 `AGENTS.md`、原任务 plan/progress、spec v1.5 和 Theory 弧长方程。
- [x] 把修复拆为 R0–R6，按生产中断风险和依赖关系排序。
- [x] 固定“不扩大历史字典类型、不静默回退、不自动重定基线、保留旧 API”的原则。
- [x] 新建 `task_plan.md`、`findings.md`、`progress.md` 并同步总索引。

### 待开始

- [x] R0：裁决失败子步回滚、Φ 配对、弧长参数和结果元数据语义。
- [x] R1：预应力标记与最终结果元数据。
- [x] R2：J2 历史与局部收敛。
- [x] R3：几何弧长方向、导数、回滚和真实收敛。
- [x] R4：`CzmSubmesh` 兼容构造器与真实 `phi_pairs`。
- [x] R5：C4 跨圈状态提交和验证证据。
- [x] R6：全门禁、基线、文档与交付。

### 修复前证据（来自本次审查，不是本计划阶段的新测试）

- 全套测试：外层 34/34，仍有 2 个弧长 `Broken`。
- 几何弧长专项：停滞在 `lambda≈0.4000000856`。
- `git diff --check`：`src/JuBat.jl:7` 尾随空格。
- 工作区：仅见既有未跟踪 `.mimosa/`、`output/testexample/`。

### 本计划阶段测试

未运行新的源码测试；本轮只建立修复计划和索引，避免把计划编写误报为实现验证。

### 错误记录

无。

### 下一次恢复提示

1. 当前处于 R0，先读取本目录三个文件。
2. 先取得用户对失败弧长子步历史回滚语义的决定。
3. 在每个实现批次前，先添加能在修复前失败的生产路径测试。
4. 每批后运行专项测试；R6 再运行全套与 `testexample` 强制基线。
5. 不要改动或暂存 `.mimosa/`、`output/testexample/`。

## 2026-08-23 — 开始执行

### 当前阶段

R0 进行中：契约已裁决，正在建立修复前失败测试并同步治理文档。

### 已完成

- [x] 用户授权执行 R0–R6。
- [x] 重新核对 `HEAD 28e5c459abe1`；未发现新增的用户源码修改。
- [x] 冻结失败子步事务回滚、真实 `phi_pairs`、弧长系数有效使用和结果元数据边界。
- [x] 确认 `CzmSubmesh` 定义实际位于 `src/SetMesh.jl`，不是原审查表中的 `src/CzmMesh.jl`。

### 工作区保护

- 保留：`.mimosa/`、`output/testexample/`。
- 本任务已有修改：`docs/planning-with-files/index.md` 与本修复规划目录。

### 下一步

补充修复前失败测试，分别锁定生产标记类型、J2 非零历史、弧长报告、旧构造器、真实 Φ 配对和 C4 状态提交。

### 修复前失败测试证据

| 测试 | 修复前结果 | 锁定缺陷 |
|---|---|---|
| `test/unit_czm_j2.jl` 新增非零历史 testset | 失败；`f_hist` 恰为单扣参考内力的 2 倍 | `eps_p` 被重复扣除 |
| `test/test_czm_submesh.jl` | 失败；全部配对为 `(i,i)`，5 参数构造器 `MethodError` | 真实 Φ 映射丢失、公开构造契约断裂 |
| `test/test_czm_phi_merge.jl` | 失败；真实配对断言不成立 | 合并拓扑覆盖了物理配对 |
| `test/test_czm_multicycle_c4lite.jl` | 对称 D10 三相位真实提交后仍为零；受控分离产生非零 D，但 clone 的三项拓扑/映射字段均为 `nothing` | 零损伤是物理工况结论；clone 契约不完整 |
| `test/test_czm_winding_prestress.jl` | 生产 `CallModel_MultiSPMe` 在 `src/CallModel.jl:151` 抛出 `Vector{Bool}` 转换 `MethodError` | R1 运行中断已由完整生产路径复现 |
| `test/test_czm_arc_geo.jl` | 停在 `λ=0.4000000855585188`，`converged=false`、`residual_norm=Inf` | R3 符号、报告和步进缺陷已复现 |

以上失败均发生在审查所指位置，未出现新的无关环境错误。

### R0 结论

R0 complete；转入 R1。环境为 Julia 1.11.2、单线程、`GKSwstype=100`、`--project=.`、`--startup-file=no`。

### R1/R2 实施与验证

- R1：`StandardVariables` 条件预分配浮点历史键；`CallModel_MultiSPMe` 写 `[1.0]`；`PostProcessing` 输出 `collapse_approx` 与启用时的预应力历史。
- R1 门禁：`test/test_czm_winding_prestress.jl` 全部通过，新增生产路径 4/4。
- R2：删除单元入口对历史 `eps_p` 的第一次扣除；返回映射改为相对应力残差收敛检查，达到上限显式错误。
- R2 门禁：`test/unit_czm_j2.jl` 全部通过；`test/test_czm_j2_integration.jl` 四组全部通过，同载荷 κ 幂等门收紧至 `rtol=1e-6`。
- 当前转入 R3；旧 geo 弧长仍在集成测试中警告停滞，尚未误报为已修复。

### R3–R5 实施与验证

- R3：geo 路径改为带 `arc_length_alpha` 的球面弧长；代码残差下修正为 `K⁻¹R`；每次迭代重算 `f_hat`；精确落到 `lambda=1` 后复算实际残差；失败诊断携带历史。
- R3 门禁：`test/test_czm_arc_geo.jl` 四组全过，包括线性 basic 等价、alpha 契约、geo=false 零漂移和 23 项一自由度软化极限点路径。
- Theory §6.10 已同步球面约束、代码残差符号、当前状态载荷方向和失败 trial 全回滚。
- R4：恢复 5 参数 `CzmSubmesh` 构造；`phi_pairs` 改回未合并网格真实节点对；`merge_phi_pairs` 独立产生 bonded 网格并提供负向守卫。
- R4 门禁：`test_czm_submesh` 6436/6436；`test_czm_phi_merge` 三组 12/12。
- R5：clone 补齐 `czm_submesh`、`thermal_to_czm`、`cohesive_to_thermal`；C4 测试提交返回损伤状态并区分对称零损伤与受控非零损伤。
- R5 门禁：C4 多相位 6/6，受控损伤 5/5。
- 当前转入 R6 综合门禁。

### R6 静态审查

- `git diff --check`：通过。
- 修复专项文件已无 `@test_broken`/`@test_skip`。
- 源码 diff 逐项复核：历史容器类型未扩大；默认关闭的预应力只新增条件键；最终字符串元数据仅在 CZM 开启时输出；clone 只补回原对象已有拓扑/映射引用。
- Git 提示工作区 LF 将按配置恢复 CRLF；这是行尾策略提示，不是内容错误。
- `rg` 发现 `src/Jellyrollmodel.jl:281` 的既有尾随空格不在本轮 diff；先不扩大清理范围。原审查登记的 `src/JuBat.jl:7` 将按计划单独核实。
- 以 `e117fd2f4d48` 为 39 个待提交提交的 review base 运行范围 `--check`，唯一命中 `src/JuBat.jl:7`；已做单行尾随空格修复。

### 源码链路核对

- `create_czm_mesh` 使用 `mesh_bonded`，真实 `phi_pairs` 可作为独立物理映射保留。
- geo 弧长当前未接收 `arc_length_alpha`，且冻结初始 `f_hat`。
- `update_czm_damage!` 的永久状态提交已位于严格收敛检查之后；修复应保持这一边界。

### 错误记录

| 错误 | 尝试 | 处理 |
|---|---:|---|
| 读取不存在的 `test/test_czm_j2_coupling.jl` | 1 | 用 `rg --files test` 定位实际文件为 `test/test_czm_j2_integration.jl`，不重复原路径 |
| R2 首次使用绝对 `tol=1e-12` 检查 Pa 量级 J2 残差，FD testset 在 `7.45e-9 Pa` 报未收敛 | 1 | 改为以当前应力/屈服强度尺度归一化的相对残差；仍保留迭代上限和强制失败测试，不接受末次值静默返回 |
| 合成极限点测试首个补丁因目标块上下文不完全匹配而未应用 | 1 | 确认源码未被部分修改，读取实际行后改为小 hunk |
| 第二个补丁把同一文件写成两个 Update 操作，`apply_patch` 拒绝 | 2 | 合并为单个文件操作内的两个 hunk 后应用，不绕过 `apply_patch` |
| 首次全量测试在沙箱内无法写入 Julia compiled cache 的 pidfile | 1 | 按环境权限错误申请并获准用相同命令在沙箱外重跑；不归因于源码 |
| 首轮全量 34 个隔离测试文件中 33 个通过；`test_czm_thin_subdiv.jl` 仍用合并后求解网格解释真实 `phi_pairs` 索引，产生 33 个坐标失败 | 1 | 按 `CzmSubmesh` 契约改为在 `czm_submesh.mesh` 上核对配对坐标；不改变生产拓扑或放宽容差，随后执行聚焦与全量复测 |
| `verify_czm_standalone.jl` 默认只跑 basic，但报告序列化固定读取 `parts[1:3]`，数值校核完成后抛 `BoundsError` | 1 | 表头和数据行按实际 `methods` 动态生成；保留默认 basic-only 与 `JUBAT_PROBE_ALL=1` 三方法语义，不改任何求解参数或结果 |
| v3 `testexample` README 与现有/重跑 PNG 均为 `272402…`、92775 bytes，但 `metrics.toml` artifact 字段仍为 v1 的 `4ba620…`、96438 bytes | 1 | 科学指标本身已是 v3 且重跑逐项一致；仅把机器可读 artifact 元数据同步到已冻结 README 和重复运行结果，不把本轮源码变化误作重新定基线 |
| standalone 动态表头首次仍用右侧填充，basic-only 生成报告留下尾随空格；沙箱内 `git restore` 又因不能创建 `.git/index.lock` 失败 | 1 | 表头末列不再填充并重生成报告；测试生成的 tracked PNG 用获批的沙箱外定向 restore 恢复，不处理任何用户原有未跟踪内容 |

### 治理一致性

- spec §4.3/§6 已明确“失败回滚全部 trial 状态”和原子提交，与本次 R0 裁决一致。
- 原任务 `task_plan.md` 的当前阶段和默认 Φ 描述仍过期，将在 R0 文档同步中修正。

### R6 综合门禁与交付

- 修复专项全部通过，且专项文件无 `@test_broken`/`@test_skip`。
- 首轮全量：33/34；唯一失败为 `test_czm_thin_subdiv.jl` 使用错误索引空间。修正后聚焦测试 49/49；第二轮全量 34/34，耗时 8m40.2s，零失败、零 broken。
- standalone v3：Nodes 10144 / Bulk 6728 / Cohesive 3364；basic 2/8、total_iter=16，八个载荷水平的 OK/FAIL、D 和打印残差逐位一致。报告序列化修复后退出 0。
- `testexample`：退出 0；1682/1763 网格、19 步与全部打印科学指标一致；PNG 92775 bytes，SHA-256 `272402bb…d386`。仅同步机器可读 artifact 元数据，没有重新定义科学基线。
- `git diff --check` 与 `git diff e117fd2f4d48 --check` 均通过；`src/JuBat.jl:7` 尾随空格已清除。测试生成并截断的 tracked eigenstrain PNG 已定向恢复到任务开始时的 HEAD 内容。
- 旧 5 参数构造器、真实 `phi_pairs`、结果键/形状和生产 `MultiSPMe` 标记均有通过测试；本轮未修改 CSV 序列化文件或列定义。
- 用户既有 `.mimosa/`、`output/testexample/` 仍保留且未纳入修复文件清单。
- **状态：R0–R6 complete。** 对称 D10 未产生 `Δ_core`/损伤，保持物理能力未验证结论，不继续扩展到接触摩擦或新标定。
