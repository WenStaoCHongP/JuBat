# JuBat 源码物理分层重构任务计划

## Goal

将 JuBat 从扁平单模块重构为标准 Julia package，形成
`Foundation -> Interfaces/Geometry -> Electrochemistry/Thermal/Mechanics -> Coupling -> Simulation -> DataIO`
的单向依赖，同时保持冻结的科学结果、结果键和 CSV 格式。

## Current Phase

Phase 1：标准 Julia package 与架构测试（in progress）

## Phases

### Phase 0：检查点、分支与基线
- [x] 将当前 tracked diff 和 untracked 文件保存到仓库外检查点
- [x] 创建 `codex/src-physics-modularization` 分支
- [x] 运行 22 个现有测试文件作为迁移前健康度记录（22/22 通过）
- [x] 运行用户指定的 7 个示例并冻结结果
- [x] 记录源码、运行时和基线哈希
- **Status:** complete

### Phase 1：标准 Julia package 与架构测试
- [x] 修复无热模型的 `temperature` 索引缺陷和 mechanical 示例旧绝对路径
- [ ] 新增 `Project.toml`、`Manifest.toml`、`test/runtests.jl`
- [ ] 建立九个顶层子模块及允许依赖规则
- [ ] 增加 package load、公开命名空间和依赖边测试
- **Status:** pending

### Phase 2：Foundation、Interfaces 与 Geometry
- [ ] 迁移通用 FEM、参数、网格和共享工具
- [ ] 分离通用 Mesh 与 CZM 专属网格类型
- [ ] 迁移 Jellyroll/Ring 几何
- **Status:** pending

### Phase 3：物理模块
- [ ] 迁移 Electrochemistry（SPM/SPMe/P2D/分流）
- [ ] 迁移 Thermal（lumped/distributed2D/polar/BC）
- [ ] 迁移 Mechanics（颗粒/宏观/CZM）
- [ ] 每层完成后运行定向测试与冻结基线
- **Status:** pending

### Phase 4：Coupling、Simulation 与 DataIO
- [ ] 用 `CoupledProblem` 和分域 context 取代旧 `Case`
- [ ] 由 Coupling 独占中央 `variables` 并实现最小输入适配
- [ ] 拆分 `CallModel_MultiSPMe` 并迁移时间积分/循环/IO
- **Status:** pending

### Phase 5：API 切换、示例与文档
- [ ] 顶层只导出子模块并删除临时 façade
- [ ] 测试和示例统一改用 `using JuBat`
- [ ] 更新 AGENTS、技术架构文档和函数索引
- [ ] 同步 `docs/planning-with-files/index.md`
- **Status:** pending

### Phase 6：最终验证
- [ ] `Pkg.test()` 通过
- [ ] 全部现有测试文件零失败（既有 broken 单独记录）
- [ ] `testexample` 严格科学指标和 PNG SHA-256 一致
- [ ] 静态架构规则与旧 API 清理检查通过
- **Status:** pending

## Locked Decisions

| Decision | Value |
|---|---|
| 模块边界 | 真实 Julia 子模块 |
| API | 允许一次破坏，采用命名空间 API |
| 范围 | 全部现有电化学、热、力学/CZM 模型 |
| 运行时状态 | 中央 `variables` 保留，由 Coupling 独占完整访问 |
| 不变量 | 参数视图、配置、布局、几何和缓存使用 struct |
| 结果 | 字符串键、单位、CSV 格式保持不变 |
| 状态/结果语义 | `src/` 中不得调用 `haskey(case.index, ...)` 或键集合成员关系猜测状态属性或固定结果维度；固定结果直接赋值，可选温度状态由热模型配置通过统一入口解析 |
| 温度结果形状 | `variables["temperature"]` 必须直接赋值为 `zeros(Float64, 1, num)`；无热模型不得读取或伪造 `case.index["temperature"]`，变量提取回退到 `case.param.cell.T0` |
| 已建立不变量 | `setup_thermal2D_mesh` 返回后 `case.geometry`、`layer_weights`、`boundary_edges` 必须存在；下游不得用 `?: nothing` 或现场重算掩盖初始化缺陷 |
| 可选缓存 | 只有函数签名确实允许缓存缺失时才保留 `cache === nothing` 分支；不存在于结构体的伪缓存字段不得用 `hasproperty` 探测 |
| 文件长度 | 无硬上限，按单一职责拆分 |
| 环境 | 标准 Julia package，Julia 1.11.2 基线 |

## Errors Encountered

| Timestamp | Error | Attempt | Resolution |
|---|---|---:|---|
| 2026-08-05 | 规划阶段 PowerShell 嵌套引号导致只读检索解析失败 | 1 | 改为分开检索和单引号模式；未修改文件 |
| 2026-08-05 | 规划阶段检索了不存在的 `scripts/` 目录导致 `rg` 返回错误 | 1 | 只检索已存在目录，并将跨平台索引工具纳入本计划 |
| 2026-08-05 | 检查点压缩未加载 `System.IO.Compression.ZipArchiveMode`，untracked zip 未生成 | 1 | 保留已成功生成的 tracked patch；改用同时加载 Compression 与 FileSystem 程序集的压缩路径 |
| 2026-08-05 | 沙箱禁止写入 `.git/refs`，无法创建 `codex/src-physics-modularization` | 1 | 工作树未变化；使用受控提升权限重试 `git switch -c` |
| 2026-08-05 | PowerShell 将 `foreach {...} | Format-Table` 解析为空管道 | 1 | 改为先收集结果数组，再单独格式化；未修改文件 |
| 2026-08-05 | `change_current.jl` 在 `thermalmodel="none"` 下因 `StandardVariables` 无条件读取 `case.index["temperature"]` 失败 | 1 | 记录为迁移前既有缺陷；继续运行其余指定脚本以识别共因，后续先加特征测试再修复 |
| 2026-08-05 | 更新阶段状态与 baseline README 的组合补丁因 `Current Phase` 上下文位置不匹配而整体失败 | 1 | 重新读取三个规划文件后按精确上下文应用；无部分写入 |
| 2026-08-05 | 并行只读检索中的可选 `rg` 返回非零，导致聚合调用无输出 | 1 | 改为逐项读取和 `Select-String -SimpleMatch`；无文件修改 |
| 2026-08-05 | `change_model.jl` 切换到 standalone SPMe 时，`SPMe_variables` 再次无条件读取温度索引 | 1 | 扩展聚焦测试，按 SPM/P2D 已有规则回退到 `case.param.cell.T0` |
| 2026-08-05 | Julia `-e` 依赖探测字符串的引号被 PowerShell 剥离，产生 ParseError | 1 | 不重试内联字符串；直接读取当前环境 TOML 获取 UUID/版本 |
| 2026-08-05 | 聚合 `rg` 正则在 PowerShell 中因双引号被重解释为未闭合分组 | 1 | 改用逐模式 `Select-String -SimpleMatch`；确认 `src/` 中三类旧温度索引存在性判断均为零 |
| 2026-08-05 | 记录温度审计结果的组合补丁因错误表文字与当前文件不完全一致而未应用 | 1 | 重新读取精确上下文后追加记录；无部分写入 |
| 2026-08-05 | 新增缓存不变量测试错误要求小网格 `boundary_edges.edges` 非空 | 1 | 修正为真实不变量：`boundary_edges !== nothing` 且边与边长数组尺寸一致；空边集合不等同于缓存缺失 |

## Baseline Gate

- Julia 1.11.2，Plots 1.40.9，1 thread，`GKSwstype=100`。
- 1682 thermal elements、1763 nodes、19 steps。
- 4.0367 -> 3.9438 V，0.0833 Ah，298.15--299.00 K。
- CZM damage/separation/fracture 全零。
- PNG 88744 bytes，SHA-256 `3538fe6ab336f9852e90566b17edbd2cd6c2c14b93c8eeff5aed01f7037df9d5`。

## Mandatory Example Regression Set

- `example/change_current.jl`
- `example/change_model.jl`
- `example/change_time_step.jl`
- `example/mechanical_example.jl`
- `example/minimal_example.jl`
- `example/testexample.jl`
- `example/thermal_example.jl`

用户确认上述 7 个示例即可作为旧示例回归范围；新增 package/架构测试仍需执行。
