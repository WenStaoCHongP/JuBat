# CZM 热插值尺寸问题进度

## 会话：2026-08-06

### Phase 3：理论复评与批准

- **状态：** 等待用户批准修改
- 已建立专项目录与诊断问题清单。
- 已追踪 `setup_thermal2D_mesh → create_czm_mesh → build_thermal_to_czm_interp → compute_czm_strain_inputs` 调用链。
- 已发现失败算例把未合并热网格用于构造插值矩阵，而求解状态来自默认合并热网格。
- 已用 Julia 对象探针精确复现两组失败尺寸，并确认合并前后单元数相同、单元映射索引仍合法。
- 已完整分类当前耦合调用、历史诊断探针和合法 standalone 调用；正确测试均使用 `case.mesh["thermal2D"]`。
- 已确认耦合示例还存在 CZM 构造早于活动热网格设置的顺序问题。
- 已从 2026-07-21 迁移记录确认正确 API 范式，无需新增公共包装入口。
- 已形成 16 个错误调用文件、1 个生产守卫文件和 1 个负向测试文件的待批准范围。
- 已完成 `Theory/04_CZM.md` 前 900 行复核：理论明确要求从唯一活动粗热场向 CZM 子网格插值，支持当前根因结论。
- 尚未修改任何业务源码。

## 测试记录

| 测试 | 结果 |
|---|---|
| `test/smoke_czm_redesign.jl`（前一批次） | `2524 != 1322`，在 `compute_czm_strain_inputs` 矩阵乘法处失败 |
| `example/testexample.jl`（前一批次） | `3366 != 1763`，在同一位置失败 |

## 错误日志

| 时间 | 错误 | 处理 |
|---|---|---|
| 2026-08-06 | 一条组合文档检索因部分模式无匹配返回 1 | 改用精确查询，未重复失败命令 |
| 2026-08-06 | 规划补丁使用了不存在的空表格锚点 | 读取当前文件并改用准确锚点 |
| 2026-08-06 | 全目录文档组合检索返回内容但最终退出码为 1 | 保留有效证据，改用 Julia 文件精确枚举 |
| 2026-08-06 | `julia --startup-file=no test/test_czm_strain_inputs.jl` 无法启动：当前 PATH 无 `julia` | 不重复原命令；改查项目配置与常见 Julia 1.11.2 安装路径后使用绝对路径 |
| 2026-08-06 | 首次绝对路径测试在 `CouplingState.jl:438` 遇到中文标点后的 Julia 字符串插值解析错误 | 将 `$n_thermal_node` 改为显式 `$(n_thermal_node)` 后重跑 |
| 2026-08-06 | PowerShell 原生 `julia -e` 长参数吞掉内部双引号，临时负向探针未解析 | 不重复 `-e` 形式；改由标准输入传递 Julia 源码 |
| 2026-08-06 | PowerShell 标准输入将临时探针的 Unicode 关键字 `nθ` 转成 `n?` | 用 ASCII Unicode 转义构造关键字符号后传递关键字 |
| 2026-08-06 | SOC 历史矩阵探针的临时断言不能用 `\u0394` 直接拼写 Julia 字段标识符 | 计算本身已完成；改用 `getproperty(..., Symbol("\u0394..."))` 验证返回字段 |
| 2026-08-06 | 当前 PowerShell/.NET 不支持静态 `SHA256.HashData`，聚合哈希输出为空 | 改用兼容的 `SHA256.Create().ComputeHash()` 实例方法 |
| 2026-08-06 | `Pkg.dependencies()` 尝试写用户 depot 日志而被权限拒绝 | 不请求扩权；改用 `using Plots; Base.pkgversion(Plots)` 只读获取版本 |
| 2026-08-06 | 当前 .NET 缺少 `Path.GetRelativePath`，首次新清单的路径列为空、聚合哈希无效 | 丢弃 `1539fe4b...`/`46217423...`；改用已验证工作区绝对前缀截取相对路径 |

## 理论复核进度

- 已完成 `Theory/04_CZM.md` 全文复核，重点核对 §3.6.0a–§3.6.3 的 C-skip-thermal 正向插值、父单元映射与损伤反向归并。
- 已完成 `Theory/06_热源.md` 全文复核，核对唯一活动温度场、损伤热反馈和热应变接口。
- 已按理论重新审阅诊断记录与当前 `CouplingState.jl`/`Jellyrollmodel.jl` 实现；确认根因不变，并识别出活动网格措辞、父单元静默跳过和 SOC 缺键伪载荷三项需要补正。

## 实施进度

- 用户已批准修改 `src/CouplingState.jl`；本阶段仅在该文件加入严格接口校验并删除静默处理，其他业务文件保持不动。
- 已确定校验位置与错误类型；保留历史矩阵最后一列语义，不保留缺键、错长或越界回退。
- 已修改 `src/CouplingState.jl::compute_czm_strain_inputs`：全部接口验证前置，删除父单元越界跳过和 SOC 缺键默认值；未触碰该文件中用户已有的其他修改。
- `git diff --check -- src/CouplingState.jl` 通过；仅有既有 Windows 行尾转换提示，无补丁格式错误。
- 正向测试 `test/test_czm_strain_inputs.jl`：7054/7054 通过。
- 进程内严格失败探针：错误节点空间、缺失 SOC、非法父单元、SOC 错长共 4/4 通过。
- SOC 历史矩阵最后一列正向探针：2/2 通过。
- `test/smoke_czm_redesign.jl` 按预期在乘法前抛出明确 `DimensionMismatch`：`M=(3969,2524)`、活动热节点 `1322`、`T_nodes=1322`；调用方尚未修复，因此该冒烟脚本当前退出码仍为 1。
- 本文件已完成，等待用户审阅；未修改其他业务文件。
- 用户已通过 `src/CouplingState.jl` 并批准修改 `test/smoke_czm_redesign.jl`；确认目标行当前无既有未提交差异，且 setup 已在 CZM 构造之前，只需替换热网格参数。
- 已将 `test/smoke_czm_redesign.jl` 的 CZM 构造参数改为 `case.mesh["thermal2D"]`，其余测试配置与断言不变。
- Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no` 冒烟测试通过：43 步，电压 `4.0367 → 3.9033 V`，`D_max(end)=0`，`δ_max_n(end)=6.7283e-14 m`，输出 `SMOKE OK`。
- 该测试文件当前为工作树未跟踪文件，因此 `git diff` 无输出；已通过内容读取确认目标行，并通过完整执行验证。
- 用户已通过冒烟测试文件并批准 `example/testexample.jl`。该文件已有用户修改：`time=60 s`、`n_theta=80`；这两项与冻结基线一致，必须保留。本次只修改 CZM 热网格参数。
- 已读取强制基线：19 步、1682 热单元、1763 热节点、最终电压 3.9438 V、最终容量 0.0833 Ah、温度 298.15–299.00 K、PNG SHA-256 `3538fe6a...f9d5`。
- 运行前脚本 SHA-256 `95fa0b...315a`、结果 PNG SHA-256 `3538fe6a...f9d5` 均与冻结基线一致。
- 已仅将 CZM 构造热网格参数改为 `case.mesh["thermal2D"]`，保留用户的 `time=60 s` 与 `n_theta=80`。
- Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no` 完整运行成功：1682 热单元、1763 热节点、19 步、初/末电压 4.0367/3.9438 V、容量 0.0833 Ah、温度 298.15–299.00 K、D_max/D_mean=0、断裂数 0；19 次 CZM 更新均收敛。
- 修复后最大法向分离为 `1.3527e-14 m`，旧基线为 `0.0000e+00 m`；结果 PNG 由 88744 bytes / `3538fe6a...f9d5` 变为 96292 bytes / `7a67c151...d915`。差异来自旧路径未成功执行 CZM、修复后得到真实的极小分离场，不应伪装成基线完全一致。
- 未更新 `Simplify/baseline/testexample/`；是否重建基线应在用户接受本文件后另行决定。
- 用户已接受 `example/testexample.jl` 的真实 CZM 分离结果，并明确授权修改基线档案。现有档案包含 `README.md`、`metrics.toml`、`source_manifest.tsv`、`preflight.log`、`run.log`，总入口为 `Simplify/baseline.md`；本次需要整体重建而非只改 PNG 哈希。
- 当前递归 `src/**/*.jl` 加入口脚本共 46 个文件，旧基线为 42 个；新拆分的活动模块必须纳入源码清单。
- 已验证旧清单标准 TSV 聚合算法：各 `path<TAB>sha256` 行以 LF 无尾换行拼接得到 `2d1a3c3d...`；字面量 `` `t`` 兼容算法得到 `40e8abd7...`。Plots 版本只读查询为 1.40.9。
- 新清单已正确生成：46 个文件，标准 TSV 聚合 SHA-256 `608286d19b5290b9f440f849a8f50f0af3d1ec3a4af4e5f0b1b49c7ee105dc4a`，字面量 `` `t`` 兼容 SHA-256 `dee4b05d085a64557b105e3c268838fc59ebcfbd02321f63b8210dba27f64bd0`。
- 正式基线捕获成功：开始 `2026-08-06T03:12:17.7098394-06:00`，结束 `2026-08-06T03:14:20.4414375-06:00`，外层 122.730 s，求解器 105.970 s，exit code 0；19 次 CZM 更新全部收敛，最大法向分离 `1.3527e-14 m`，未生成调试日志。
### 2026-08-06 基线源码清单更新重试

- 首次整体替换 `source_manifest.tsv` 未应用：补丁中的旧 `src/Variables.jl` 哈希锚点存在录入错误，文件未发生变化。
- 处理：按当前档案中的精确旧哈希修正锚点后重试，不改变目标内容。
- 已完成：`source_manifest.tsv` 更新为当前 46 个入口/源码文件，聚合 SHA-256 为 `608286d19b5290b9f440f849a8f50f0af3d1ec3a4af4e5f0b1b49c7ee105dc4a`。
- 待完成：把正式复跑输出写入 `run.log`，随后交叉校验全部基线文件。
- 当前桌面任务未附带可回读终端，旧 `run.log` 仍是 2026-08-05 记录；正式复跑的时间戳、总耗时、分阶段耗时和科学结果已经写入 `metrics.toml` 与本任务记录，可据此更新归档日志。
- 已确认未留存另一份临时控制台日志；正式运行只生成了新的结果 PNG。`metrics.toml` 与 README 已包含新的基线 ID、46 文件哈希、19/19 CZM 收敛、`1.3527e-14 m` 分离量和新 PNG 哈希。
- `run.log` 首次多段更新未应用：PowerShell 显示编码与补丁中的中文锚点不一致，补丁整体被拒绝且日志保持原样。下一次改用可稳定匹配的 ASCII 行分段更新，并单独定位中文结果行。
- 已完成 `run.log` 的 ASCII 与中文数值分段更新。首次并行校验因 Julia `-e` 的 PowerShell 引号转义错误中止，属于校验命令问题，不影响档案内容；改用 PowerShell 单引号包裹 Julia 表达式后重试。
- 第一轮交叉校验：`git diff --check` 通过；清单 46 行全部与当前文件 SHA-256 一致；标准聚合、兼容聚合、PNG 字节数与 PNG 哈希均与 `metrics.toml` 一致。路径顺序比较显示 false 是因为清单约定把示例入口放在 `src/` 之后，而全局排序把 `example/` 放在前面，需改为集合比较。
- Julia `-e` 仍被外层 PowerShell 剥离内部引号，TOML 尚未完成语法验证；下一次通过 here-string 把只读校验脚本送入 Julia 标准输入，避免命令行转义。
- 基线档案最终校验通过：TOML 可由 Julia 1.11.2 解析；46 个清单路径集合完全匹配且全部文件哈希一致；两种聚合哈希、PNG 96292 字节及其哈希一致；基线索引、README、preflight、run log 与 metrics 的关键标记一致；`git diff --check` 通过。
- 同步 `task_plan.md` 状态时使用了不存在的待办原文，组合补丁整体未应用；进度事实已单独记录，后续按计划文件现有措辞更新。
- 已将计划推进至 Phase 5。下一业务文件确定为 `example/coupled_czm_thermal_example.jl`：其第 104 行在活动热网格安装前用 `mesh_data.thermal2D` 构造 CZM，而第 118 行才安装并取得 `case.mesh["thermal2D"]`；当前等待用户按逐文件流程批准后再调整顺序与参数。
