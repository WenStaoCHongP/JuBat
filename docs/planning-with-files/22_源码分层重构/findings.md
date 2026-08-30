# JuBat 源码物理分层重构调研记录

## Current Architecture

- `src/JuBat.jl` 将 36 个顶层源码文件 include 到同一命名空间并直接导出大量底层符号。
- `CallModel_MultiSPMe` 同时承担分流、电化学单元求解、热源、热矩阵、CZM 映射和全局拼装，是首要耦合热点。
- 旧 `Case` 同时持有电化学网格、热布局、CZM 网格和缓存，物理模块普遍接收完整 `Case`。
- 中央 `variables` 通过字符串键传递跨物理场数据；新架构保留其结果兼容性，但完整访问仅允许在 Coupling。
- SPM/SPMe 在矩阵装配时直接调用粒子力学，因此需通过 Interfaces 中立协议注入力学响应。
- `SetMesh.jl` 同时定义通用 Mesh 和 CZM 容器；迁移时必须拆开以解除基础层对力学概念的依赖。
- 根目录没有 `Project.toml`；现有测试/示例大多直接 include `src/JuBat.jl`。
- `md/源码函数索引` 基于旧提交 `96cee65`，行号和 ThermalDistributed 结构已经过期，且没有现成可复用生成器。

## Dependency Policy

- Foundation 无项目内部依赖。
- Interfaces 仅依赖 Foundation。
- Geometry 仅依赖 Foundation。
- Electrochemistry、Thermal、Mechanics 仅可依赖 Foundation、Interfaces、Geometry，不得互相依赖。
- Coupling 可依赖三个物理模块。
- Simulation 依赖 Coupling 和独立物理求解入口。
- DataIO 依赖 Simulation 的结果类型/格式，不得被求解层反向依赖。

## Worktree Constraint

- 当前分支 `czm-refactor`、HEAD `df8aab3e278243fa9efe0759b666f4e386628276`。
- 工作树包含大量用户既有 tracked/untracked 改动；不得 reset、stash 或混入不相关提交。
- 正式编辑前将 tracked binary diff 与 untracked 文件另存到仓库外检查点，再创建重构分支。

## Pre-existing Example Failures

- `change_current.jl` 使用默认 `thermalmodel="none"`。`SetCase` 在此模式不创建 `index["temperature"]`，但 `StandardVariables` 当前无条件读取该键，因此在进入时间循环前抛出 `KeyError("temperature")`。
- 该问题可能共同影响 change_model、change_time_step、minimal 和 mechanical 示例；需逐个运行确认。
- 实测确认 `change_current`、`change_model`、`change_time_step`、`mechanical_example`、`minimal_example` 均在同一位置失败。
- `SPM_variables` 和 `P2D_variables` 已有无温度 DOF 时使用 `case.param.cell.T0` 的正确分支；缺陷局限于 `StandardVariables` 的分配尺寸。
- 首轮修复后发现 `SPMe_variables` 与 workspace 变体没有采用 SPM/P2D 的温度回退规则；standalone SPMe 仍会读取缺失索引。
- `mechanical_example.jl` 使用不存在的 `C:\Users\user\Desktop\JuBat\src\data\`，仓库内实际数据位于 `src/data/pybamm_SPM_1C_Stress_t.csv`。
- `testexample` 通过且冻结 PNG 哈希一致；`thermal_example` 通过并生成两份 PDF。

## Temperature Result Constraint

- 根因：`thermalmodel="none"` 时状态向量没有温度自由度，因此 `case.index` 合法地不含 `"temperature"`；结果数组初始化却把状态索引误当作结果形状来源，导致 `KeyError("temperature")`。
- 结果语义：顶层 `variables["temperature"]` 是每个结果步的代表温度历史，不是温度自由度向量；对所有模型都保持 `1 × num`。
- 重要实现约束：固定结果数组的维度/属性由其公开结果语义决定，不得调用 `haskey(case.index, ...)` 进行条件判断或推导；`case.index` 只描述实际状态自由度。
- 温度直接赋值：`StandardVariables` 必须使用 `variables["temperature"] = zeros(Float64, 1, num)`，不再声明 `n_temperature`，也不通过 `case.index` 推导行数。
- 回退约束：无温度自由度时，SPM、SPMe、P2D 的变量提取使用 `case.param.cell.T0`；不得为了满足结果输出而向状态索引添加虚假温度自由度。
- 测试约束：静态检查 `StandardVariables` 使用直接赋值且不恢复原 `haskey` 维度判断；运行时对无热模型的 SPM、SPMe、P2D 验证 `haskey(case.index, "temperature") == false` 且温度结果尺寸为 `(1, num)`。
- 全源码审计：`SPM.jl`、`P2D.jl` 各有一处通过 `"temperature" in var_list` 猜测温度状态，`SPMe.jl` 有两处通过 `haskey(case.index, "temperature")` 猜测；现统一调用 `representative_temperature`，由 `thermalmodel` 明确选择环境初温、状态温度或调用方提供温度。
- 保留边界：工作区键过滤、CZM 节点复制和 CSV 映射中的 `haskey` 检查的是普通映射成员关系，不涉及状态属性或结果维度，保留其原语义。

## Optional-cache Ternary Audit

- `setup_thermal2D_mesh` 总会构造 `MeshGeometry`，其中包含 `layer_weights` 和 `boundary_edges`；`ThermalDistributed2D`、`ThermalDistributed2D_BC`、`compute_heat_sources` 都只能在热网格 setup 后使用，因此三处 `case.geometry !== nothing ? ... : ...` 属于掩盖初始化错误的冗余回退。
- `Case` 没有 `I_e_cache` 字段，也没有任何写入点；`CallModel_MultiSPMe` 的 `hasproperty(case, :I_e_cache) ? ... : nothing` 永远返回 `nothing`，已删除临时变量并在调用处分流求解器时明确传 `nothing`。
- `MeshGeometry` 只有 `setup_thermal2D_mesh` 中一个构造点，且该构造点总会生成 `BoundaryEdgeCache` 和 `czm_element_map`；因此 `boundary_edges::Union{Nothing,...}` 以及下游对 `czm_element_map` 的 `hasfield/geom !== nothing` 防御均被删除。
- CZM 函数公开签名允许 `cache=nothing`、`u0=nothing`、`F_ext=nothing`，且测试会直接调用这些无缓存路径；相关三元式承担真实默认值语义，不能按表面形式删除。

## Sources

- `Simplify/baseline/testexample/metrics.toml`
- `md/源码函数索引/_索引.md`
- `md/10_参数传递与模块架构.md`
- `src/JuBat.jl`
- `src/CallModel.jl`
- `src/SetCase.jl`
- `src/CouplingState.jl`
