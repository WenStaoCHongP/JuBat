# 静默回退逐文件修复记录

## 2026-08-15 `src/CycleSolver.jl` 初查

- 开始实施前确认三个目标文件均已有用户工作区改动：`CycleSolver.jl` 和 `CsvExport.jl` 为已跟踪修改，`CyclePostProcess.jl` 为未跟踪新文件。现有改动属于此前后处理职责拆分/严格 CSV 契约，本轮必须在其上做最小增量，不能恢复基线版本或覆盖这些内容。
- 精确接口调整：`CycleResult.rest1/rest2` 改为 `Union{PhaseResult,Nothing}`，构造循环结果时直接初始化为 `nothing`；正时长静置才赋 `PhaseResult`。`postprocess_cycle_result!` 的最高温度从必有的充/放电阶段开始，再显式纳入存在的静置阶段，避免让“阶段不存在”变成数值零。
- `postprocess_phase_result` 仅需把 `initial_state` 类型扩为 `Union{Dict,Nothing}`；只有 REST 分支需要并直接使用非空初始状态。首个放电阶段由 `Solve(initial_state=nothing)` 按模型初始化，后续阶段仍直接传递严格的 `final_state` 字典。
- CSV 的三个阶段循环必须共享“实际阶段列表”语义，顺序保持 `discharge → rest1（若存在）→ charge → rest2（若存在）`；实际存在但缺少 `solve_result` 的阶段仍由现有严格 helper 报错，不能用可选静置适配重新引入 `continue` 回退。
- 调用审计只发现示例对 `rest1/rest2` 的直接读取，现有示例休止时间均为正，正常路径仍得到 `PhaseResult`；零时长的新 `nothing` 语义需要新增针对性测试，避免以后又被默认对象替代。
- `Solve` 在 `initial_state === nothing` 时已有正式的模型初始化路径；非空状态则严格要求 `"y"` 数组，因此删除 CycleSolver 的占位字典不会丢失任何必要输入。`Solve` 返回的 `final_state` 同时保存真实 `y`、`T_nodes`、`V` 与局部阶段时间，`postprocess_phase_result` 负责把 `t_global` 改写为累计时间。
- 多 SPMe 温度自由度的权威位置是 `case.layout.thermal_range`，且 `final_state["T_nodes"]` 只是其副本。温度重置必须同时写 `y[thermal_range] .= case.param.cell.T0` 并重新同步 `T_nodes`；只清空副本不会影响下一阶段求解。
- `ModelInitialisation` 也会为非多 SPMe 的 lumped/distributed2D 模型追加热自由度，但其索引结构与多 SPMe 不同；本项目当前循环主路径是 Jellyroll 多 SPMe。温度重置 helper 应先以模型现有权威索引为准，不能猜测“末尾若干项”并误写 P2D 势自由度；需继续核对 `case.index["temperature"]` 的构造后再定通用范围。
- `SetCase` 对 lumped 温度给出单一 `case.index["temperature"]`，对 distributed2D 只记录首个热自由度；结合热网格 `nlen` 可构造连续节点区间。拟实现的重置因此覆盖：无热模型为明确无操作；多 SPMe 用 `layout.thermal_range`；lumped 用现有温度索引；非多 SPMe distributed2D 从记录的首索引取连续 `nT` 项。所有路径都直接改写 `final_state["y"]`，并同步 `T_nodes`。
- 两处重置调用也要求 `thermalmodel != "none"`，因此无热模型不会执行 helper，也不会再打印“温度场已重置”的虚假状态提示。
- 现有边界测试不构造真实 `CycleResult`，也没有零静置/温度重置契约覆盖。本轮不额外修改第四个持久测试文件，以保持批准的三业务文件边界；验证阶段用独立 Julia 合约探针覆盖可选阶段列表、汇总温度、首阶段 `nothing` 和热自由度重置，再运行既有边界测试及强制 `testexample` 基线。
- 强制 `testexample` 本次除最大法向分离和 PNG 外均匹配档案：本次分离 `1.2557e-14 m`，基线 `1.3527e-14 m`；本次 PNG SHA-256 `4BA6207C3CCF92DA5E37349EE335CF21A10A50B46A14CDA13DE95EEFA6CAE932`，基线 `7A67C151272EB1CE7A1BC3F18BC9181CFF9C03B9214385DF5D6129445A2FD915`。按强制规则不能判为通过，须先检查当前源码相对基线清单的漂移来源。
- 基线源码清单已明显过期：46 项中有 16 项当前哈希不同，包括 `example/testexample.jl` 本身，以及 CouplingState、Czm/CzmBC/CzmUnitMesh、Jellyrollmodel、Mechanical、Parallelsolution、SetMesh、ThermalDistributed、Tools、Variables 和本轮三个循环/导出文件。因求解/CZM/网格文件与入口脚本都已在基线后改变，当前 PNG 与分离量差异不能归因于本轮 CycleSolver 接口修改，也不能用旧档案判定本批的前后等价性。
- 当前工程已有本批前直接行为基准：CycleData 上一批完成时 `testexample` 同样得到最大分离 `1.2557e-14 m`、PNG 96438 bytes、SHA-256 `4BA6207C3CCF92DA5E37349EE335CF21A10A50B46A14CDA13DE95EEFA6CAE932`；更早 ThermalDistributed 批次还用本批前函数内存对照确认过同一哈希。当前本轮结果与该本批前工作树基准精确一致，旧归档差异源自已接受的热/力共享周向拓扑变更（力学段数 1760→1682），不是本轮接口修改。
- 当前待冻结清单仍为 46 个文件；标准 TSV 聚合 SHA-256 为 `77CC5DF31B1D1CE21C6329052546DB8B12BB9F2CC3B6E74DF7E1D76E36B9A0DF`，兼容字面量 `` `t`` 聚合为 `1CF823C5C77E23F85A8353B0AF0E34A3ECDB2D6C21866F2334266A287BD867DE`，入口脚本哈希为 `8C19E11619444BAC6954F01845EED676D281EC22B7AA40E2A88B7E367F06B9E9`。
- 为避免把上一次工具调用推算的开始时间写成精确实测，正式基线将再运行一次同一入口，在命令内部打印开始/结束时间和退出码；控制台结果用于完整更新 run.log，随后重新计算 PNG 与清单哈希。

## 2026-08-15 `src/SetParams.jl` 初查

- `ChooseCell` 只对五个已知字符串执行参数文件 `include`，没有未知类型分支；参数文件把 `param_dim` 写入模块全局。未知类型随后仍直接使用 `param_dim` 计算孔隙率、尺度和归一化锚点，因此行为取决于模块中是否已有上一次成功调用留下的全局参数：可能复用错误电芯，也可能在首次调用时报未定义变量。
- 五个合法公开名称为 `"LG M50"`、`"Northrop"`、`"Enertech"`、`"Jellyroll"`、`"Ring"`；默认无参调用等价于 `"LG M50"`，属于明确 API 默认值，不是未知字符串的回退。
- 参数文件均以 `param_dim = Params(...)` 结束，印证 `ChooseCell` 当前依赖 include 的模块全局副作用；最小修复应让未知名称在任何历史调用状态下都立即失败，不能改成默认电芯。
- 运行时已复现两种历史相关结果：新模块首次 `ChooseCell("__unknown__")` 抛 `UndefVarError: param_dim not defined`；先调用 `ChooseCell("Jellyroll")` 后再传未知字符串，则返回与 Jellyroll 完全相同的对象（`reused === valid` 为 true，I1C=5.0）。后一种会把拼写错误静默解释为上一次电芯，是本文件的结果失真点。
- 现有测试大量覆盖合法 `ChooseCell("Jellyroll")`，但未发现未知名称契约测试。最小验证应覆盖：未知名称在首次调用和已有合法全局状态后均抛同一明确异常；五个合法名称及无参默认返回仍保持原行为；`testexample` 强制基线保持新档案一致。
- 源码索引还列出 `cell.rho/heat_Q` 的显式厚度加权派生，以及缺失 E_coat/cohesive 时伴随 `@warn` 的非 CZM 兼容路径；它们不是本队列记录的“未知电芯复用旧对象”问题，本批不改材料派生、CZM 锚点或归一化公式，避免扩大批准范围。
- `src/SetParams.jl` 在本轮开始前已有未提交修改；拟修位置仅为 `ChooseCell` 型号分支末尾（当前约 L268），应叠加一个未知名称的直接异常，不覆盖现有参数/归一化改动。
- 修复完成：`ChooseCell` 型号选择链末尾现在直接抛出列明五种支持型号的 `ArgumentError`；未知名称不再读取模块中的旧 `param_dim`。首次未知调用与合法调用后的未知调用得到相同错误，五种显式合法型号和无参默认型号均保持可用。
- `example/testexample.jl` 以 Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no` 运行成功：1682 单元、1763 节点、19 步；4.0367→3.9438 V，0.0833 Ah，298.15–299.00 K，CZM 19/19 收敛。输出 PNG 为 96438 bytes，SHA-256 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932`，与正式基线完全一致。
- `git diff --check -- src/SetParams.jl` 通过；全工作区检查仍报告其他既有文档第 43–44 行尾随空格，不属于本批次。

## 2026-08-15 `src/Initialisation.jl` 初查

- `ModelInitialisation` 对电化学模型已有最终 `else`，但热模型只处理 `"lumped"` 与 `"distributed2D"`，其他任意字符串均直接跳过温度自由度追加；因此未知热模型与明确的 `"none"` 当前产生同样的纯电化学初态。
- 已确认工程合法热状态规则目前是 `"none"`、`"lumped"`、`"distributed2D"`：`StateAccess.representative_temperature` 对前者返回环境初温、对后两者读取温度状态，并已对其他字符串直接抛 `ArgumentError`。初始化层与该严格契约不一致。
- `ModelInitialisation_MultiSPMe` 会暂时把 `thermalmodel` 改为 `"none"` 以构造单元化学初态，这个内部明确值属于合法路径，不应被当作未知模型修复。
- `Option.jl` 对耦合状态初始化公开的热模型注释只列 `"none"`、`"lumped"`、`"distributed2D"`；`"ring2D_polar"` 仅出现在 `Solve.jl` 的独立纯热 `model == "thermal"` 分支，不经该电化学状态初始化入口，因此不应加入 `ModelInitialisation` 的合法值集合。
- 多 SPMe 初始化当前通过“保存字符串→改成 `"none"`→调用→恢复”获取纯化学初态；如果内部初始化抛错，恢复语句不会执行，原案例会被遗留成 `thermalmodel="none"`。这不是正常结果路径的回退，但会污染失败后的可复用案例；可通过不改写选项的纯化学初始化拆分，或至少 `try/finally` 保证恢复。是否纳入本文件同批修复需结合最小边界决定。
- 首轮运行时探针已在后续软作用域断言前确认：未知热模型生成的空默认状态与显式 `none` 完全相等，自定义 `opt.y0=[7,8]` 也会原样绕过热模型分支；因此若只在“自动生成 y0”的分支末尾加 `else`，仍不能覆盖自定义初态，合法值检查应位于函数入口。
- 函数作用域复测确认多 SPMe 失败污染：内部化学初始化抛 `TypeError` 后，原案例的 `thermalmodel` 从 `distributed2D` 留在了 `none`。因此本文件除未知值检查外，还存在一个明确的状态恢复缺口。
- `SPM`、`SPMe`、`P2D` 的变量计算最终都通过 `representative_temperature`，所以外部 `initial_state` 绕过初始化时，未知热模型一般会在首次模型计算处报错，而不是完成正常仿真；但直接调用公开的 `ModelInitialisation` 仍会返回错误布局，自定义 `opt.y0` 也仍会被接受，初始化入口本身的契约必须修正。
- `src/Initialisation.jl` 当前相对 Git 基线无未提交差异，可安全叠加本批最小修改；现有 `test/test_example_preconditions.jl` 是未跟踪用户文件，只可作为只读契约参考，未经批准不修改。
- 旧的简化审计把本文件列为“仅审查、无修改”，但该结论没有覆盖未知热模型的运行时等价性、`opt.y0` 绕过或多 SPMe 异常后的选项污染；当前实测证据优先于旧审计结论。
- 拟修改边界收敛为同一业务文件内两点：① `ModelInitialisation` 入口先限定 `none/lumped/distributed2D`，确保自动初态和自定义 `opt.y0` 都不能绕过；② `ModelInitialisation_MultiSPMe` 明确要求 `distributed2D`，并用 `try/finally` 恢复临时改写的热模型。两点都只暴露无效状态，不增加任何替代初态或回退值。
- 拟验证：未知值在空/非空 `opt.y0` 下均为同一 `ArgumentError`；三种标准模型初态长度/温度值保持；多 SPMe 只接受 distributed2D、正常布局不变、内部失败后选项恢复；随后运行现有初始化测试、`git diff --check -- src/Initialisation.jl` 与 `example/testexample.jl` 强制基线。
- 针对用户质疑重新核对后，临时改写 `case.opt.thermalmodel` 并非模型需求，只是为了复用标准初始化函数并抑制热自由度追加。多 SPMe 真正需要的是一个单元的 SPMe 化学向量 `[csn0; csp0; ce0]`，可在本文件内通过私有纯化学初始化 helper 直接获得，因此不应保留“改写再恢复”的方案。
- `opt.y0` 仅在 `Initialisation.jl` 被读取，未发现多 SPMe 专门文档；当前临时改写也不能从自定义 `y0` 中剔除温度，只会原样返回并复制。新 helper 可保持当前非空 `opt.y0` 的原样语义，同时让默认空 `y0` 直接构造化学部分，不使用截断、猜测长度或临时配置。
- 修订后的拟修改：删除多 SPMe 对 `thermalmodel` 的保存/改写/恢复；抽取仅构造默认 SPMe 化学初态的私有 helper，供标准 SPMe 空 `y0` 与多 SPMe 空 `y0` 共用；未知热模型入口检查仍保留。因不再存在临时状态，`try/finally` 也不需要。
- 修订方案已实现并通过聚焦契约：默认 SPMe 化学初态由 `_initialise_spme_chemical_state` 唯一构造；标准三种热模型的布局保持，多 SPMe 正常状态长度为 9081；未知/不适用模型直接失败且案例选项在成功、失败路径都不发生改变。
- `example/testexample.jl` 强制基线通过：1682 单元、1763 节点、19 步，4.0367→3.9438 V、0.0833 Ah、298.15–299.00 K、D=0、最大分离 `1.2557e-14 m`、断裂数 0；PNG 96438 bytes，SHA-256 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932`，与正式基线精确一致。
- 用户随后要求回退本文件修改；`src/Initialisation.jl` 已恢复到本批开始前内容，未知热模型静默等价于 `none` 及多 SPMe 临时改写选项的原行为均保留。Git 定向状态与差异为空。
- `opt.y0` 全仓调用复核：字段只在 `Option.jl` 定义、只在 `Initialisation.jl` 读取；生产、测试和示例均无赋值者。续算路径使用 `Solve(initial_state=...)`，多 SPMe 非均匀起点使用 `initial_soc_distribution`，因此删除该遗留字段不会切断仓内状态输入路径。
- `src/Option.jl` 与 `src/Initialisation.jl` 在本原子批次开始前均无 Git 差异，可直接实施且无需合并既有用户改动。
- 用户最终决定保留 `opt.y0` 公共字段与标准初始化检查；因此 `Option.jl` 和 `ModelInitialisation` 不属于最终业务差异，最终仅修改 `ModelInitialisation_MultiSPMe`，使其不再通过临时热模型值间接构造化学初态。
- 最终多 SPMe 语义保持：`opt.y0` 为空时按参数构造单元 `[csn0;csp0;ce0]`，非空时把 `vec(opt.y0)` 作为单元化学初态逐单元复制；两者随后统一追加热节点温度。唯一移除的是临时改写热模型及对标准初始化函数的间接调用。
- 最终强制基线通过：1682 单元、1763 节点、19 步，4.0367→3.9438 V、0.0833 Ah、298.15–299.00 K、D=0、最大分离 `1.2557e-14 m`、断裂数 0；PNG 96438 bytes，SHA-256 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932` 精确匹配。最终正文 diff 仅 `src/Initialisation.jl`；`Option.jl` 工作树/索引对象哈希同为 `58ceabc267928f0e7456db8d170ca77d41d659d3`。

- `solve_phase` 仍用 `get(initial_state, "t_global", 0.0)`，缺失全局时间时静默从零开始；随后又用空时间历史回退 `t_max`、缺失 `final_state` 回退空字典。后两段结果随后由严格的 `postprocess_phase_result` 重新读取，属于无效的旧回退/死代码，应删除并让严格后处理直接失败。
- `solve_cycling` 的首阶段状态仍构造为 `Dict("y"=>nothing, ...)`，与当前 `Solve` 的严格契约冲突：非空状态字典必须含数组状态。正确的首阶段语义应是显式 `nothing`，由 `Solve` 正常初始化；这需要 `solve_phase` 接受 `Union{Nothing,Dict}`，并联动放宽 `CyclePostProcess.postprocess_phase_result` 的非 REST 初态类型。
- `reset_T_each_cycle` 与 `reset_T_before_charge` 仅把冗余键 `current_state["T_nodes"]` 置为 `nothing`；当前 `Solve` 实际从 `current_state["y"]` 读取热自由度，因此温度并未重置，界面却打印“温度场已重置”。必须改为更新状态向量中的真实热自由度，不能继续把清空旁路键当作重置。
- 零时长静置通过默认 `PhaseResult()` 伪造完整阶段，仅补 duration/V/final_state，导致 `t_start/t_end/T_max/T_mean/D_*` 保持默认零值。应把可选静置阶段表达为 `nothing`，而非伪造结果；这会直接联动 `CycleResult` 字段与 `CyclePostProcess.postprocess_cycle_result!`。
- `solve_cycling` 对无效 `SOC_init` 仍静默跳过，和刚审阅通过的 `CycleData` 同类问题相同；可直接调用已有严格 `apply_initial_soc!`，无需新增范围检查。
- 仓库没有现成的循环温度重置 helper；多 SPMe 状态的真实热自由度由 `case.layout.thermal_range` 明确定义，初始归一化温度就是 `case.param.cell.T0`。因此温度重置可在 `CycleSolver.jl` 内直接、严格地更新 `current_state["y"][thermal_range]`，并同步 `"T_nodes"`，无需改动模型方程或增加回退。
- 把可选静置改成 `nothing` 会影响 `CsvExport.jl` 和示例对 `rest1/rest2` 的直接访问；当前严格 CSV 导出也要求每个阶段有 `solve_result`。因此“缺失阶段”不能只改字段类型后结束，必须同步定义可选阶段的导出语义，否则会把静默失真转成下游崩溃。
- 原专项审查的明确修复目标是 `solve_phase` 对缺失 `time [s]`/`final_state` 的回退；当前 `CyclePostProcess.jl` 已严格要求这些字段，因此删除 CycleSolver L159–165 的重复旧回退即可保持单一失败语义。
- 实测把 `PHASE_REST` 以 `t_max = 0.0` 交给现有 `Solve` 并不能表达“未发生的阶段”：求解器仍进入初始化/时间推进路径并产生历史数组扩展警告，探测脚本最终也未获得可用的零时长 `PhaseResult`。因此不能用一次零时长求解替代可选阶段语义，也不应重复该探测；未发生的静置应显式为 `nothing`。
- 拟议原子修复边界：`CycleSolver.jl` 删除阶段结果旧回退、首阶段传入真实的 `nothing`、直接应用 `SOC_init`、真实改写 `layout.thermal_range` 温度自由度，并把零时长静置保留为 `nothing`；直接依赖的 `CyclePostProcess.jl` 仅放宽首阶段后处理的初始状态类型并在循环汇总中忽略明确不存在的静置；`CsvExport.jl` 仅导出实际存在的阶段。三者是同一可选阶段接口契约，需用户明确批准后作为跨文件原子批次实施。
- 1 秒、nθ=20 的实际 `solve_cycling` 探针确认当前首阶段在进入求解前即抛 `ArgumentError`；根因正是 `current_state["y"] = nothing` 与严格 `Solve` 状态契约冲突，已从静态风险升级为当前必现故障。

## 2026-08-07 `src/CycleData.jl` 初查

- 专项原记录中的两项高风险回退已经存在于当前未提交工作区差异中：非空 `initial_state` 缺少 `"y"`/`"t_global"` 时已显式报错，不再重新初始化或把全局时间置零；未知 `solveType` 已抛出 `ArgumentError`，不再当作 backward。
- 文件仍有三处候选静默回退：无效 `SOC_init` 被条件分支静默忽略（L315–321）；热网格为空或总面积非正时，面积加权温度退回节点算术均值（L211–216）；`soc_mean` 在 `ne==0` 时写入 `0.0`（L229）。后两项均会掩盖理论上不允许的空热网格/零面积状态。
- `solve_phase_with_export`/`solve_cycle_with_export` 是公开 API，且分别有示例调用，不能删除或把本轮扩大成 CycleSolver 架构重构。
- `CycleSolver.jl` 对 `SOC_init` 使用同样的静默范围分支，因此若本轮只改 CycleData，会产生两个循环入口契约不一致；是否跨文件统一须在报告中明确并等待用户决定。
- `apply_initial_soc!` 已通过 `compute_cs0_from_soc` 原生拒绝 `[0,1]` 之外的 SOC，因此 CycleData 无需新增校验；删除外层范围分支、直接调用即可避免“无效 SOC 静默沿用旧浓度”。CycleSolver 的同类分支留待其下一文件批次处理，保持逐文件审阅。
- `solve_phase_with_export` 接受 `czm_mesh`/`czm_params`，但函数体不使用这两个参数，结束时无条件把 `PhaseResult.D_max`、`D_mean`、`ΔD_max` 写成 0。直接调用者若传入 CZM 网格会得到伪造的零损伤结果；这是本文件另一处独立的静默失真。
- 标准 `solve_phase` 会在阶段前后读取传入 CZM 网格的真实损伤统计；CycleData 应复用同一统计语义，但不在本轮展开 callback/求解器去重架构重构。
- 进一步确认 `solve_phase_with_export` 自己复刻时间循环、只调用 `CallModel`；CZM 损伤更新实际位于 `Solve` 的时间循环中，因此本函数即使挂载网格也不会推进损伤。单纯把最终统计从 0 改为读网格仍会伪装成“支持 CZM 但损伤不演化”。
- 因此本轮不能把 `czm_mesh` 接上后宣称功能正确。最小诚实契约应是：该在线热/SOC 导出入口明确不支持 CZM 演化，启用 CZM 时立即拒绝；真正支持需留到 CycleData/CycleSolver callback 去重扩展，不能用固定零值继续。
- 公开文档已把 CycleData 定义为电化学—热在线快照入口，并注明固定零损伤不替代完整 CZM 循环；两个现有示例也都设置 `czm_enabled=false`。因此拟保留公开关键字以避免签名破坏，但当 `case.opt.czm_enabled`、`czm_mesh` 或 `czm_params` 表示 CZM 请求时立即抛错，正常无 CZM 路径的零损伤才具有明确物理含义。

### 拟修改范围（仅 `src/CycleData.jl`）

1. `solve_phase_with_export` 拒绝任何 CZM 演化请求，避免忽略参数并伪造零损伤；不在本批复制 `Solve` 的 CZM 时间循环。
2. `solve_cycle_with_export` 直接调用 `apply_initial_soc!`，复用其既有范围失败语义，删除非法 SOC 时沿用旧浓度的静默分支。
3. 删除空热网格/零面积时改用节点均温的分支，始终按热单元拓扑计算面积加权均温；不新增面积有限性或正值判断。
4. `soc_mean` 直接使用 `mean(soc_n)`，空场依赖原生失败，不写入伪造的 `0.0`。
5. 已严格化的 `initial_state` 与 `solveType` 逻辑保持不变；`CycleSolver.jl` 的同类 SOC 分支留到下一文件单独报告。

### 拟验证

- Julia 1.11.2 语法和 `git diff --check`。
- 定向验证非法 SOC、CZM 请求、空 SOC 场不再静默继续。
- 运行 `test/test_cycle_data_import_removed.jl` 和一个短时无 CZM `solve_phase_with_export` 正常路径。
- 按项目约定运行 `example/testexample.jl` 强制行为基线；不更新基线档案。

### 实施结果

- 修改严格限定在 `src/CycleData.jl`；`CycleSolver.jl` 未改。
- 无 CZM 正常导出路径保持可用；任何 CZM 请求不再进入不会推进损伤的伪支持路径。
- 初始 SOC、面积均温和 SOC 均值均不再使用旧值、节点均值或物理零替代异常状态。
- 当前工程 `testexample` 科学结果及 PNG SHA-256 与本批前一致；旧归档基线仍按此前 Jellyroll 拓扑问题单独处理，本批未更新。

## 2026-08-07 ThermalDistributed 基线对照

- 当前 `testexample` 与档案在退出码、网格（1682/1763）、时间步（19）、电压、容量、温度范围、损伤统计上相同。
- 当前最大法向分离为 `1.2557e-14 m`，档案为 `1.3527e-14 m`；当前 PNG 为 96438 bytes、SHA-256 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932`，档案为 96292 bytes、SHA-256 `7a67c151272eb1ce7a1bc3f18bc9181cff9c03b9214385df5d6129445a2fd915`。
- 当前 `example/testexample.jl` 自身 SHA-256 为 `65e20c5d4844d7ffed2a3c135edddb79683424d7412156d835fadffd1d3acd0c`，档案记录为 `a72d33f1e25aa52bc44a504e35b90b2e321ffd381a5ad4801de63ac9b8c0ae27`，说明基线入口已先于本批发生变化；需通过工作区差异定位来源，不能直接把差异归因于 ThermalDistributed 本批。
- Git 差异确认当前 `example/testexample.jl` 相对基线提交已把仿真时长从 3600 s 改为 60 s、周向网格从 360 改为 80，并改用共享热/力周向拓扑；`src/Jellyrollmodel.jl` 也包含此前已批准的拓扑重构，因此旧 PNG/最大分离档案不是本批前状态。
- ThermalDistributed 本批四项在正常路径上不改变求解状态：`variables_elems` 在逐单元 SPMe 路径存在时选支相同；移除 `@inbounds` 只恢复边界检查；活动索引合法时直接赋值结果相同；`dot(q_total, areas)` 只写入后处理键 `"total heat source"`，源码检索显示求解器不读取该键反馈到状态方程。

- 采用同一当前工程状态，内存恢复本批修改前四处写法后重跑 `testexample`：退出码、1682/1763 网格、19 步、全部科学输出均与修改后相同；修改前后 PNG 均为 96438 bytes，SHA-256 均为 `4ba6207c3ccf92da5e37349ee335cf21a10a50b46a14cda13de95eefa6cae932`。因此本批满足严格前后行为一致；旧档案应在后续获批后按此前 Jellyroll 拓扑变更单独重建，不能归入本文件修改。

## 用户约束

- 一次只修改一个业务文件。
- 每个文件修改前必须报告具体问题并等待批准。
- 每个文件修改、验证完成后必须停止并等待用户审阅。

## 审查依据

- `docs/planning-with-files/静默回退专项审查/findings.md`

## 2026-08-06 恢复执行校准

- 原计划记录早于后续 CZM—热插值及热力网格拓扑修复，不能直接以旧“待审阅”状态判断当前源码。
- 当前 COH2D4 损伤父热单元映射已经改为拓扑直连；不再使用几何近邻搜索。数值探针确认 6728 个 cohesive 中宿主父映射不一致数为 0，且每个热单元恰有 4 个 cohesive。
- `CallModel.jl:12` 仍对越界 `cohesive_to_thermal` 静默跳过；这不是近邻回退，但会把损坏映射解释为零损伤，后续应纳入修复队列。
- 剩余队列必须以当前源码重新核对；工作树差异显示多个候选文件已被后续批次修改，不能重复执行旧方案。
- 当前 `Parallelsolution.jl` 已在 Newton 未收敛时直接报错，并对越界失效单元索引抛出 `ArgumentError`；用户此前明确认定总面积/全部有效面积为零属于理论根本错误，不为该不可能状态另加代码回退。
- 当前 `Mechanical.jl` 已移除线性求解失败后的零位移替代，并将退化 Q4 中心 Jacobian 从静默 `continue` 改为报错；两项原队列问题均已落地。
- 当前 `Variables.jl::Variable_update!` 已拒绝非二维历史矩阵、多列当前矩阵、行数不匹配、标量/数组契约错配，不再截断、取首值或保留预分配零列；`AGENTS.md` 已记录这些严格判断在完整测试证明冗余后才可独立简化。
- `CsvExport.jl` 仍有下一批明确回退：阶段汇总把 `NaN` 改写为 `0.0`；缺失阶段/时间/字段时 `continue` 或使用空矩阵；缺失面积写 `0.0`；CZM 快照按最短数组少写单元并对缺失分离/牵引补零。
- `CsvExport.jl` 中 `should_export_step`、`should_export_snapshot` 的 `continue` 是用户配置的主动抽样，保留；`write_csv_guarded!` 会发出警告并把失败文件加入 `files_skipped`，属于显式批量导出策略，也保留。
- 当前需要严格化的活动路径包括：cycle summary 的 SOH/NaN 处理；element/node 历史的 `solve_result`、时间、必需矩阵及网格行数；cohesive 快照的 5 个 `n_coh` 数组；节点位移的 `2*nnode` 长度。
- `cohesive_driving_force.csv` 当前因热/扩散参数尚未接入而显式打印原因并返回，后续主体不可达；本批不把这个未实现功能伪装成已支持，也不扩大到其模型设计。
- 当前 `Tools.jl::IntQ4` 已对零、负或近零 Jacobian 报错，`compute_separation` 也已对零长度 cohesive 边报错，原队列的两个主要问题已经落地。
- `Tools.jl` 仍有三处同类残余：体积平均温度在积分总权重非正时回退节点算术平均；`q4_center_gradients` 对奇异中心 Jacobian 返回 `nothing`，且使用 `abs(detJ)` 会接受负 Jacobian；`compute_separation` 在 Newton-Cotes 总权重非正时保留初始化的零分离。
- `thermal2D_volume_average_temperature` 由 `CallModel_MultiSPMe` 每步直接调用，网格积分失效时回退节点平均会产生看似合理的温度；其风险高于单纯诊断辅助函数。
- `q4_center_gradients` 当前唯一调用方 `Mechanical.jl` 已检查 `nothing` 并报错，但把失败留给调用方会允许未来调用遗漏；适合在工具函数内部统一采用正 Jacobian 契约。
- `Tools.jl` 修改后的双线性 CZM 测试暴露 `CzmUnitMesh.jl` 的条带 Q4 节点序为 `[bottom-left, top-left, top-right, bottom-right]`，在项目 Q4 基函数约定 `[bottom-left, bottom-right, top-right, top-left]` 下属于顺时针翻转，全部 32 个高斯点 detJ 为负。
- 数值探针确认这不是全项目的“负号约定”：Jellyroll 力学网格 13440 个高斯点和热网格 1680 个高斯点全部 detJ 为正；只有单元条带夹具为负。
- 因此不能把 `Tools.jl` 改回接受负 detJ；正确修复是另行调整 `CzmUnitMesh.jl` 的 Q4 连接顺序。该文件不在本批授权范围，需用户批准跨文件处理。

## 已完成：`src/CyclePostProcess.jl`

- `T_mean_end` 读取不存在的 `"thermal2D T_nodes [K]"`，当前会回退为初始温度。
- 多个必需的 `Solve` 结果字段缺失时被空数组、空字典或计划时长代替，阶段仍被包装成有效结果。
- 静置阶段缺失浓度状态时，方差和松弛率被写成零。
- 缺失终止原因时默认解释为时间终止。
- 充电容量为零时，未定义的库仑效率被写成 `0.0`。

### 已完成修改

1. 将 `time [s]`、`cell voltage [V]`、`final_state`、`termination_reason` 改为必需键，缺失或空历史立即报错。
2. 分布式热路径改读 `Solve.jl` 实际提供的 `thermal2D final temperature at nodes [K]`；用节点末温计算 `T_mean_end`，用单元温度历史计算 `T_max`。
3. 非分布式热路径使用通用 `temperature [K]` 历史，不再用缺失二维温度的方式回退到 `T0`。
4. `phase_termination_symbol` 保留 REST 阶段无条件归类为 `:time` 的既有契约；充放电阶段对未知终止原因报错。
5. 静置阶段要求初末状态均含 `y`；缺失时不再输出零方差/零松弛。
6. 充电容量非正时不再输出伪造的零库仑效率，而是暴露无效循环结果。

### 拟验证方法

- Julia 1.11.2 源码语法检查。
- 运行现有 `test/test_postprocessing_boundaries.jl`。
- 使用不写文件的 Julia 命令构造最小结果字典，验证正确末温、缺键失败和未知终止原因失败。
- 本步骤不修改测试文件；若需要新增永久测试，将作为独立文件另行报告并等待批准。

## 批准记录

| 文件 | 报告 | 用户批准 | 修改 | 用户审阅 |
|---|---|---|---|---|
| `src/CyclePostProcess.jl` | 已报告 | 已批准 | 已完成并验证 | 通过 |

## 等待审阅：`src/Solve.jl` + `src/CouplingState.jl` 联动批次

### 具体问题

1. 传入的 `initial_state` 缺少 `"y"` 时重新初始化模型；多 SPMe 状态长度不匹配时也只警告并重新初始化，阶段连续性被静默切断。
2. 多 SPMe 初始变量缺少热节点温度时使用空向量；时间步解不足以包含温度自由度时沿用旧的 `T_nodes_carry`，可能造成电化学状态和温度状态错步。
3. 纯热入口把任何既非 `ring2D_polar` 也非 `ring2D` 的 `thermalmodel` 当成 `distributed2D`，无效配置会运行另一种模型。
4. CZM 更新的任意异常仅写 `@debug` 后继续；返回 `converged=false` 时同样继续推进时间，损伤失效不会终止主求解。
5. CZM 快照把牵引力固定为零、迭代次数固定为 0、残差固定为 0.0；当前 `update_czm_damage!` 只返回位移与收敛标志，`Solve.jl` 无法单独生成这些真实诊断量。

### 已实施范围

1. 传入状态字典时强制要求 `"y"`；多 SPMe 状态长度不匹配时抛出 `DimensionMismatch`，不再重新初始化。
2. 多 SPMe 模式强制要求初始热节点温度存在且尺寸正确；每个时间步若解向量不足以提取全部温度自由度则立即报错，不再沿用旧温度。
3. 纯热入口仅显式接受 `ring2D_polar`、`ring2D` 和 `distributed2D`，其余值抛出 `ArgumentError`。
4. 移除吞掉 CZM 异常的 `catch`；CZM 返回未收敛时立即终止 `Solve`，同时保留耗时统计的 `finally` 语义。
5. 经用户授权联动修改 `CouplingState.jl`：`update_czm_damage!` 返回完整 `CZMResult`；`Solve.jl` 使用其中真实的位移、损伤、分离、牵引力、迭代次数和残差构造快照。
6. `CouplingState.jl` 对 CZM 非有限输入、非有限结果和不收敛直接报错，只在结果有效且收敛后提交损伤与位移状态。

### 拟验证方法

- Julia 1.11.2 源码语法检查。
- 使用不写文件的最小测试验证缺失 `"y"`、错误多 SPMe 状态长度和未知纯热模型均立即失败。
- 运行现有 CZM/热耦合相关测试，确认正常路径不受影响。
- 运行 `example/testexample.jl` 的既有强制行为基线；本步骤不修改测试文件。

## 批准记录（续）

| 文件 | 报告 | 用户批准 | 修改 | 用户审阅 |
|---|---|---|---|---|
| `src/Solve.jl` | 已报告 | 已批准 | 已完成并验证 | 等待 |
| `src/CouplingState.jl`（直接接口联动范围） | 已补充报告 | 已授权跨文件调整 | 已完成并验证 | 等待 |
| `src/CsvExport.jl` | 已报告 | 已批准 | 已完成并验证 | 通过 |
| `src/Tools.jl` + `src/CzmUnitMesh.jl`（验证暴露的直接拓扑依赖） | 已报告 | 已批准跨文件修正 | 已完成并验证 | 通过 |
| `src/CzmBC.jl` | 已报告 | 已批准 | 已完成并验证 | 通过 |
| `src/Jellyrollmodel.jl` | 已报告 | 已批准 | 已完成并验证 | 通过 |
| `src/ThermalDistributed.jl` | 已报告并按用户意见收缩 | 已批准 | 已完成并验证 | 通过 |
| `src/CycleData.jl` | 已报告 | 已批准 | 已完成并验证 | 通过 |

## `src/Tools.jl` 验证中的跨文件问题

- 本批已删除体积平均温度和分离平均的替代值分支，并使中心梯度拒绝奇异/翻转 Jacobian。
- 独立严格契约测试通过：正向 Q4 正常；顺时针和退化 Q4 报错；正常 cohesive 零位移得到零分离；零长度边报错。
- 重新核查当前 `Czm.jl`：`assemble_czm_system` 的缓存与非缓存路径，以及 `build_czm_cache` 的几何缓存构建路径，均已用 `L >= 1e-15 || error(...)` 拒绝零长度/非有限 cohesive 切向长度；原记录中的“零贡献继续组装”已不存在，因此本文件无需业务修改。
- `CzmBC.jl::apply_bc_czm` 以矩阵最大对角量级计算罚系数，但矩阵为空、全零、全零对角或 `dmax=NaN` 时会回退固定 `1e12`；这会把缺失/非法系统刚度伪装成可继续求解的边界约束，应改为要求有限且严格为正的刚度尺度。
- 同一入口使用 `zip(bc_dofs, bc_vals)`，长度不等时会静默截断较长一侧；只提供 `bc_dofs` 或只提供 `bc_vals`、两种 BC 接口同时提供、或完全不提供 BC 时，也会静默忽略部分/全部边界条件。应建立互斥且完整的输入契约，并校验 DOF/值等长。
- `bc_nodes` 中未知的边界类型当前无 `else` 分支，会被静默跳过；应仅接受 `:fixed_x`、`:fixed_y`、`:fixed_xy`，否则明确报错。
- 当前生产调用点均从 `CzmSolve.jl` 成对传入 `bc_dofs` 与 `bc_vals`，因此建立“恰好选择一种 BC 接口、DOF/值成对等长且非空”的严格契约不会改变合法求解路径。
- `CzmBC.jl` 已按批准范围严格化：矩阵必须非空、方阵且维数与载荷匹配；BC 模式必须互斥且完整；节点/DOF/数值/类型均预校验；罚系数只允许由有限正对角尺度计算，固定 `1e12` 回退已删除。
- `unit_czm_eigenstrain.jl` 的实际 BC 数组完整、等长、非空、无重复且全部为零，因此走的正常路径与修改前仍是同一 `penalty=1e6*dmax` 公式。其 58/60 失败来自最终四个界面均未形成正向张开（3 个为压缩、1 个近零），与新增输入校验无关；这是独立的载荷/法向/测试预期问题，本批不跨文件修复。
- `Jellyrollmodel.jl` 初筛发现：`nθ < 3` 会被 `max(3,nθ)` 静默提升；相位对齐后不足两个区间时会改用忽略 `phase` 的端点均分网格；单元层号小于 1 时会被 `max(1, ...)` 改成第一层。这三处都会把非法/不满足离散契约的输入或几何改写成可继续计算的网格。
- 极耳弧长换算中的 `max(1e-6, ...)`、二分容差下限及末尾区间中点需要进一步区分“数值算法初值/收敛近似”与“静默替代物理几何”；暂不把它们直接定性为应删除的回退。
- 原专项记录对应的两项当前已消失：层权重计算不再把退化 `dtheta` 提升为 `1e-10`，而是由 `A_total > 0` 拒绝零面积；`Jellyrollmodel.jl` 当前也已无 `thermal_elem_map` 的 `clamp`，该映射已在此前网格拓扑批次迁为严格父单元拓扑。因此不能按旧行号直接修改。
- 正常 `nθ=80` Jellyroll 探针得到 1682 个热单元，未截断层号范围为 1–22，`max(1, ...)` 在正常路径未触发；Φ 重合节点对为 1603，恰等于 1683 个周向节点减去一圈 80 段。删除这些静默改写预计不会影响该基线参数。
- `build_czm_submesh` 当前已严格要求 `n_segments == n_thermal`，并直接令每层 `thermal_elem_map[elem]=seg`；旧父单元索引夹边界问题确已修复。
- 但 `thermal_phi_pairs` 的上游仍由 `jellyroll_collector_seed_mesh` 用坐标容差双重循环挑选，未按均匀角节点的一圈索引偏移直接生成，也不校验应有配对数量。正常参数恰好得到正确的 1603 对，但容差或采样异常会静默少配/错配；这与此前已确认的 `outer(θ) ↔ inner(θ+2π)` 拓扑契约不完全一致。
- 负向探针确认静默行为：`nθ=1`、`nθ=3`、`nθ=-7` 均生成相同的 63 单元网格；`tol=-1e-8` 仍因 `tol^2` 合并重合节点，却因坐标配对条件 `abs(...) < tol` 得到 0 个 Φ 对，函数正常返回。这会产生“热网格已合并、力学 Φ 拓扑为空”的自相矛盾结果。
- 拟修改范围限定在 `jellyroll_collector_seed_mesh`：要求 `nθ>=3`、`tol` 有限且为正、有效相位采样至少形成一个区间（两个角节点），且卷绕跨度足以形成 Φ 对；按一圈 `nθ` 的列索引偏移直接生成并验证 Φ 对；合并映射只使用这些拓扑对；层号计算小于 1 时直接报错，不再夹成第一层。极耳二分求解的正常数值初值、容差与区间中点保留。
- 现有相关测试均使用合法分辨率 `nθ=20/40/60/80`，不存在依赖 `nθ<3` 自动提升的测试路径；CZM 子网格和建网测试可直接验证拓扑改动。
- `Jellyrollmodel.jl` 已按批准范围修改：删除低分辨率与相位端点重采样，Φ 对改为一圈列索引偏移并逐对验证重合，合并映射仅采用这些拓扑对，层号下界不再夹取；源码语法检查通过。
- `czm_element_map` 仍是活动兼容字段；本批使其沿累计 θ 的拓扑配对顺序构造，不再用 `atan` 将多圈节点折回 `[-π,π]` 后排序。缓存不变量测试通过。
- 源码残余检索确认低 `nθ` 提升、端点重采样、几何近邻配对、`atan` 重排序、`merged_to` 容差合并和层号 `max(1,...)` 均已删除；当前只保留显式输入/跨度判断及拓扑配对。
- 重新核查当前 `ThermalDistributed.jl`：原专项记录中的 `sigma_pcc/sigma_ncc = max(..., 1e-12)` 已不存在，集流体电导率现在直接参与除法，不能按旧方案机械删除回退。
- 但当前热源入口未要求固相、集流体及逐单元电解液有效电导率有限且为正；零值会产生 `Inf`，负值会产生负焦耳热并继续写入科学结果。虽然不再是“提升到微小正值”，仍需用严格参数/逐单元检查使根本材料或状态错误立即失败。
- `compute_heat_sources_with_czm` 对 `active_elements` 中越界索引使用范围判断后静默跳过，随后把所有未标记单元热源清零；损坏的活动单元映射会被解释成失活，属于同文件内确定的结果失真回退。
- 两个热源入口均为生产活动路径：`CallModel.jl` 在 CZM 启用时调用 `compute_heat_sources_with_czm`，否则调用 `compute_heat_sources`。`get_active_elements` 正常返回当前热单元布尔掩码的有效索引，因此可直接索引，越界时依赖原生 `BoundsError`，不需要额外回退或替代值。
- `apply_convection_bc` 的 `Bi==0` 和 `apply_cool_method` 的无极耳返回属于零换热/用户已确认的 `tab -> none` 模型极限；未知 `cool_method` 也按用户“操作人员只设置合法值”的前提保留，本批不扩大这些既定决策。
- `get_active_elements` 的当前实现最终返回 `findall(active)`，正常情况下索引天然位于 `1:ne`；因此删除范围判断、直接执行 `is_active[e]=true` 即可让未来映射错误原生失败，不需新增回退。
- `SetParams.jl` 只归一化 PCC/NCC 电导率及电解液电导函数，没有正值检查；标准 Jellyroll 参数为正，但错误/缺失参数会抵达本热源入口。
- 同一热源函数还有一处更直接的模型回退：当 `per_element_spme=true` 但 `variables_elems===nothing` 时，当前条件分支会转而使用全局电化学变量，等于把缺失的逐单元状态解释为所有热单元共用一个状态，应明确报错。
- 主循环带 `@inbounds`，因此至少应在循环前严格要求 `I_e`、`T_e`、`variables_elems`（逐单元模式）与热单元数一致；`areas` 也应严格等长，否则长度 1 会通过广播把同一面积用于所有单元总功率。
- `StandardVariables(case,1)` 为该入口提供 `ne×1` 热源工作数组，生产调用也总是 `num=1`；修改时应同时要求 `layer_weights` 为 `ne×5`、各热源工作数组长度为 `ne`，避免 `@inbounds` 掩盖容器尺寸错误。
- 标准 Jellyroll 归一化参数探针均为有限正值：`sig_n_eff=492.1771`、`sig_p_eff=0.36536`、`sigma_pcc=1.08355e8`、`sigma_ncc=1.81915e8`、`ce0=1.0`，严格检查不会改变正常基线。
- 拟修改范围限定在 `ThermalDistributed.jl`：入口严格校验关键数组/权重尺寸和逐单元状态存在性；固相与集流体有效电导率、每单元三段电解液平均电导率必须有限且为正；CZM 活动索引直接写掩码、不再越界跳过。不改变无冷却、无极耳或合法冷却方式语义，也不增加替代值。
- 用户否决上述电导率有限/正值判断：电导率属于材料数值定义和操作人员输入责任，不在代码层增加防呆。修订方案删除全部电导率检查，也不新增其他材料参数有效性判断。
- 为保持实现简洁，剩余修复改用原生失败语义：逐单元分支只判断 `per_element_spme` 后直接索引 `variables_elems[e]`，缺失状态自然失败；移除主循环 `@inbounds` 使尺寸错误原生失败；总功率用严格长度的内积替代可广播乘法；CZM 活动索引直接写掩码，越界自然失败。不新增显式防呆分支。
- `ThermalDistributed.jl` 已按收缩方案完成四处机械修改；残余检索确认未添加电导率或其他材料参数判断，源码语法检查通过。
- `test/smoke_thermal_bc.jl` 18/18 通过，`unit_czm_strip_mesh.jl` 拓扑检查 14/14 通过。
- `unit_czm_bilinear.jl` 因 `CzmUnitMesh.jl` 的 8 个体单元全部顺时针而失败；需要把连接从 `[bl, tl, tr, br]` 改为项目 Q4 约定的 `[bl, br, tr, tl]` 后才能完成回归。
- 经用户批准已修正条带连接顺序；条带 32/32 高斯点现均为正，detJ 范围 `1.1208e-9` 至 `7.95768e-9`。
- 最终验证：条带拓扑 14/14；Mode I 单调 90/90、卸载重载 83/83、Mode II 61/61；热边界 18/18；耦合装配 8/8；`git diff --check` 通过。

## 等待审阅：`src/CsvExport.jl`

### 已实施范围

1. cycle summary 强制 SOH 数量与循环数量一致；阶段中的 `NaN` 原样写入，不再改成 `0.0`。
2. element/node 导出在打开文件前要求每个阶段有字典型 `solve_result`、非空时间向量及尺寸精确的必需历史矩阵。
3. 单元面积必须与热单元数一致，写入时直接索引，不再对缺失面积补零。
4. cohesive damage 导出要求每个选中快照的 damage、法/切向分离和法/切向牵引均严格为 `n_cohesive` 长度。
5. node displacement 导出要求位移严格为 `2*nnode` 长度；不再按较短数组少写节点。
6. 删除 `safe_csv_matrix_value` 的 `NaN` 越界回退；主动抽样 `continue` 和显式 `files_skipped` 策略保持不变。
7. 所有数据契约检查在创建目标 CSV 前完成，避免失败后留下部分科学数据文件。

### 验证结论

- `git diff --check -- src/CsvExport.jl` 通过。
- `test/test_csv_export_guard.jl`：10/10 通过。
- 临时目录负向契约测试通过：NaN 原样输出；SOH、矩阵、CZM 数组和位移长度错误均明确失败，且未创建目标文件。

## 本批次验证结论

- 两个源码文件均通过 Julia 1.11.2 语法解析；`git diff --check` 通过。
- 缺失 `initial_state["y"]`、错误多 SPMe 状态长度和未知纯热模型的无文件行为测试通过，均能明确失败。
- CZM 接口、应变输入、CSV 防护、热到 CZM 插值、损伤映射和缓存测试全部通过。
- `test/smoke_czm_redesign.jl` 与强制基线 `example/testexample.jl` 均在首次 CZM 更新时暴露原有尺寸错误：`thermal_to_czm_interp` 的列数分别为 2524/3366，而热状态长度分别为 1322/1763。旧 `Solve.jl` 会吞掉该异常继续求解；当前代码正确终止，因此基线不再伪通过。
- 该尺寸错误属于下一项工程问题，本批次未擅自修改网格映射算法。

## `src/CyclePostProcess.jl` 实施结果

- 必需的时间、电压、终止原因与最终状态改为直接读取，缺失时立即报错；时间或电压历史为空同样报错。
- 分布式热路径改用实际输出的末态节点温度键计算 `T_mean_end`，用单元温度历史计算 `T_max`；非分布式路径要求通用温度历史存在。
- 状态方差计算要求有效状态向量和足够长度，REST 阶段不再把缺失状态解释为零松弛。
- 保留 REST 无条件映射为 `:time` 的既有契约；活动阶段的未知终止原因改为报错。
- 充电容量非正时抛出 `DomainError`，不再把未定义的库仑效率写成零。
