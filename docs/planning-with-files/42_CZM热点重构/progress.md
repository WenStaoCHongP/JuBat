# CZM 热点重构：进展记录

## 会话：2026-08-31

### 阶段 1：任务基线与工作区审计

- **状态：** 进行中
- **开始：** 2026-08-31
- 已完成：
  - 检查当前技能与可调用工具，确认没有 `superpowers:*`。
  - 在本地技能目录中检索，只有 `planning-with-files`，未找到 Superpowers。
  - 检索内存索引和仓库中的 CZM/Superpowers 线索。
  - 读取 `planning-with-files` 完整说明及三份模板。
  - 创建本任务三件套。
  - 安装并验证 14 个官方 Superpowers 技能；本轮恢复后读取 `using-superpowers`、Codex 适配说明与 `brainstorming`。
  - 将任务分类为架构型重构：设计获用户批准前只做只读审计，不改源码。
  - 盘点当前分支、近期提交和工作区：tracked 文件无改动，HEAD 含 `98adb9a` CZM 热点定位探针。
  - 对照旧向量化计划与 2026-08-19 综合性能任务：确认已回退的 Task 7 不是本轮默认实现方向，LU 分解是首要复核对象。
  - 检查 `98adb9a`：仅新增只读 profile 探针，未改求解器。
  - 搜索当前 CZM 求解分支中的线性求解点，确认 basic/load-substep/arc-length/几何路径均需在设计前划清范围。
  - 完整阅读最新热点探针、历史根因再定位记录、basic 求解循环和 `apply_bc_czm`。
  - 确认首轮设计应以 basic 非几何路径为证据范围，矩阵复用必须比较施加 BC 后的实际数值内容。
  - 复核四层结构后的状态/缓存所有权：`MechState` 只承载演化状态，性能缓存应随 `CohesiveMesh`/workspace。
  - 复核统一分派与装配入口，决定不在缺乏证据时同时改动 load-substep/arc-length 分支或复做旧 Task 7。
  - 确认 `assembly_workspace(czm_mesh)` 跨时间步复用，LU 缓存可归入 workspace 而不污染 `MechState`。
  - 运行 Julia 1.11.2 稀疏解法等价性小实验：显式 `lu` 复用与直接 `A\\b` 未达到位级一致，历史可行性论证需重做。
  - 运行当前 `tools/probe_czm_hotspot.jl`：退出码 0；8,771 个样本栈中 6,831 个落在 `CzmSolve.jl:234`，约 77.88%，确认 basic 线性求解仍是主热点。
  - 发现探针百分比当前使用全部栈帧总数作分母，打印的 9.83% 不能作为墙钟占比；后续设计需修正统计口径。
  - 验证 `factorize(A)\\b` 与 `A\\b` 在 Julia 1.11.2 的对称/非对称稀疏小矩阵上逐位一致，确定不使用强制 `lu`。
  - 检查未跟踪 `output/ab_probe.jl`，确认其 API 已过时，仅保留为用户历史文件，不执行不修改。
  - 在真实 `nθ=80` 生产 `K_bc` 上验证 `factorize`：直接解与因子回代逐位相等；分解约 2.39 s、回代约 0.024 s，缓存命中有明确收益。
  - 发现已有 `41_位移协调源码审查`，将本任务三件套从冲突编号迁移为 `42_CZM热点重构` 并同步总索引。
  - 应用户追问审计逐元素数组路径：定位公开 `compute_separation` 与生产 `assemble_czm_system` 循环，确认后者已有几何/workspace/mul! 优化且当前 profile 占比很低。
  - 逐行阅读两条分离位移实现：记录公开工具函数的临时分配、生产路径的原位 workspace，以及可选的直接跳量/直接 `B_local` 标量化方向。
  - 按 systematic-debugging 区分性能与物理根因：日志证明零损伤近线性状态仍由 CZM 占据 88%–96% 时间，当前主因不是参数触发的损伤非线性，而是大稀疏系统每步重复分解。
- 已创建/修改：
  - `docs/planning-with-files/41_CZM热点重构/task_plan.md`
  - `docs/planning-with-files/41_CZM热点重构/findings.md`
  - `docs/planning-with-files/41_CZM热点重构/progress.md`

## 测试结果

| 测试 | 预期 | 实际 | 状态 |
|---|---|---|---|
| 当前技能目录检索 | 找到 Superpowers 或确认缺失 | 未找到 Superpowers | 已确认 |
| 任务文件名精确检索 | 定位热点重构计划 | 无直接匹配 | 改用正文/残留项检索 |
| 当前 CZM 热点探针 | 复核四层重构后主热点 | `K_bc\\R_bc` 命中 6,831/8,771 样本栈，约 77.88% | 通过 |
| 真实 `K_bc` 因子 A/B | 保持生产求解逐位一致 | `factorize(K_bc)\\R_bc == K_bc\\R_bc`，最大差 0 | 通过 |

## 错误日志

| 时间 | 错误 | 次数 | 处理 |
|---|---|---:|---|
| 2026-08-31 | 文件名搜索未定位“CZM 热点重构” | 1 | 改为读取既有向量化与综合性能计划正文、勾选状态和当前源码 |
| 2026-08-31 | Git 无权读取用户级 `.config/git/ignore` | 1 | 命令本身成功；保留 warning，不修改用户全局 Git 配置 |
| 2026-08-31 | 真实矩阵实验 `include_string` 的 `@__DIR__` 错误 | 1 | 未运行求解；下一次用 `abspath(tools/probe_czm_hotspot.jl)` 作为 include filename |
| 2026-08-31 | 新任务编号 41 与既有目录冲突 | 1 | 三件套迁移到 42，保留既有 41 不变 |
| 2026-09-01 | production 缓存接线测试缺少热化学单元数组 | 1 | 传入三组维度正确的显式零载荷数组，不改生产代码 |
| 2026-09-01 | 记录 fixture 错误的补丁上下文不匹配 | 1 | 查看实际文件尾部后以准确上下文重试 |
| 2026-09-01 | 尝试读取不存在的 `output/testexample/metrics.toml` | 1 | 改读冻结档案 `Simplify/baseline/testexample/metrics.toml` 并逐项对照控制台指标 |
| 2026-09-01 | scoped diff 检出 spec 头部 3 行尾随空格 | 1 | 清理这 3 行；保留无关 LF→CRLF 与用户级 ignore 权限 warning |

## 5 问题恢复检查

| 问题 | 回答 |
|---|---|
| 当前在哪？ | 阶段 1：任务基线与工作区审计 |
| 下一步去哪？ | 定位治理计划、盘点工作区、复核实际热点 |
| 目标是什么？ | 在不改变科学/API 契约的前提下完成 CZM 热点重构 |
| 已学到什么？ | 见 `findings.md` |
| 已做什么？ | 见本文件阶段 1 记录 |

## 会话：2026-09-01

### 阶段 2：方案 A 规范与实施计划

- **状态：** 完成
- 用户批准实施方案 A，并确认当前内聚力参数为正确实验值。
- 重新读取 `using-superpowers`、`planning-with-files`、`writing-plans`、`executing-plans`、`test-driven-development`、测试准则、worktree 与完成验证说明。
- 检查隔离状态：当前是普通 checkout，分支 `codex/src-physics-modularization`；保留现有任务规划与用户输出，后续执行 scoped 修改和验证。
- 明确本批不修改任何内聚力参数、本构、收敛判据、物理状态提交或公开结果。
- 对照既有 2026-08-19 性能 spec/plan 与当前 basic 求解、BC、机械核心测试，锁定缓存接线点和测试风格。
- 定位 workspace 定义与 basic 精确行号，确定新增缓存仍归属 mesh/workspace，矩阵判据覆盖尺寸、CSC 结构与数值内容。
- 创建方案 A 设计规格 `docs/superpowers/specs/2026-09-01-czm-basic-factorization-cache-design.md`。
- 创建对应实施计划 `docs/superpowers/plans/2026-09-01-czm-basic-factorization-cache.md`，下一步进入 TDD RED。
- 对 spec/plan 执行自检：修正 Julia `\` 运算符示例，并补齐数值失效、结构失效和生产入口测试的可执行代码。
- 自检结果：计划占位模式 0 处；发现并修正最后 1 处双反斜杠示例。
- 运行源码修改前 `example/testexample.jl`：退出码 0、19 步、墙钟 55.976 s、CZM 44.804 s（89.94%，2358.126 ms/call），科学锚点见 findings。
- 新建 `test/test_czm_factorization_cache.jl`，覆盖缓存复用、数值/结构失效和 production basic 接线；尚未修改生产源码。
- TDD RED：聚焦测试退出码 1，3 个失败均为预期的 helper/缓存字段不存在；0 个 error，fixture 与语法正常。
- 实现 `CZMAssemblyWorkspace` 的矩阵快照/分解字段，以及严格 CSC 内容比较和成功后提交缓存的求解 helper；尚未接入 basic 求解行。
- helper 行为测试已通过 8/8 与 5/5；production fixture 首次因缺少 `dT_elem/Δsoc_n_elem/Δsoc_p_elem` 出现 MethodError，已补为与 bulk 单元数一致的显式测试输入后重跑。
- production 接线 RED：helper 13 个断言通过，真实 basic 入口仅两个缓存非空断言失败，位移有限且无 error。
- 将 `solve_czm_basic_step` 的非几何求解接到缓存 helper；`geo_nl=true` 继续执行原 `K_bc \ R_bc`。
- GREEN：`test/test_czm_factorization_cache.jl` 退出码 0，三个测试集 8/8、5/5、4/4，共 17/17 通过；进入聚焦回归与科学基线验证。
- scoped diff 审查未发现范围扩张；直接 basic 回归目标为 `unit_czm_newton.jl`，另运行 `test_czm_mech_core.jl` 验证 workspace/装配契约。
- 聚焦回归：6 个测试文件依次运行并全部退出码 0；新缓存 17/17，机械核心所有测试集、geo C1、J2 integration 和 arc-geo 全部通过，`unit_czm_newton.jl` 也正常退出。
- 后置 `testexample.jl` 退出码 0、19 步，控制台科学锚点与前置基线一致；CZM 44.804→6.943 s（约 -84.50%），墙钟 55.976→18.237 s（约 -67.42%）。
- 对照 `Simplify/baseline/testexample/metrics.toml`：mesh、19 步和全部 result 记录精度一致；确认快速脚本本身不写 output metrics 文件。
- 最终 scoped 门：whitespace issue 0；参数文件差异 0；tracked 源码仅 `CouplingState.jl`、`CzmSolve.jl`。
- 任务完成，未提交 Git，保留当前 `codex/src-physics-modularization` 分支及所有无关用户文件原状。
