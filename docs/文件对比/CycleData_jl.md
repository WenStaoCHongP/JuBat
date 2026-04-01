# CycleData.jl

## 文件状态: 新增 (A)

## main分支
- 不存在（全新文件）

## Parameters_Design分支
- 行数: 623
- 主要内容: 循环仿真数据导出功能，支持将电化学-热耦合仿真的温度场和 SOC 数据导出为 CSV，供后续 CZM 分析使用

## 结构列表

### `TimeStepData`
- `time::Float64` - 物理时间 (s)
- `phase::PhaseType` - 阶段类型 (充电/放电/静置)
- `V::Float64` - 电压 (V)
- `I::Float64` - 电流 (A)
- `T_nodes::Vector{Float64}` - 节点温度 (K)
- `T_max::Float64` - 最高温度 (K)
- `T_mean::Float64` - 平均温度 (K)
- `soc_n::Vector{Float64}` - 负极各单元 SOC
- `soc_p::Vector{Float64}` - 正极各单元 SOC
- `soc_mean::Float64` - 平均 SOC

### `CycleExportData`
- `cycle_idx::Int` - 循环编号
- `timesteps::Vector{TimeStepData}` - 所有时间步数据
- `node_coords::Matrix{Float64}` - 节点坐标 (nT x 2)
- `element_connectivity::Matrix{Int}` - 单元连接 (ne x 4)
- `ne::Int` - 单元数
- `nT::Int` - 节点数

## 函数列表

### `solve_phase_with_export(case, phase_type, t_max, I_current, V_limit, initial_state; czm_mesh, czm_params, dt_range, export_interval)`
- 核心功能: 执行单阶段仿真并导出每个时间步的数据
- 支持放电/充电/静置阶段
- 时间步策略: 自适应时间步 (`dt_min`, `dt_max`)
- 离散格式: Crank-Nicolson (`theta=0.5`) 或向前/向后 Euler
- 温度平均: 使用面积加权平均（而非简单算术平均）
- 导出间隔控制: `export_interval` 参数
- 奇偶步伪振荡抑制: 使用步间中点温度导出
- 终止条件: 时间达到上限 或 电压达到截止值
- 返回: `(PhaseResult, Vector{TimeStepData})`

**关键实现细节**:
1. 初始化状态支持 multi-SPMe 和标准 SPMe 两种模式
2. 每步先计算电化学，再提取温度场，最后导出数据
3. 步后重新调用 `CallModel` 以获取步后变量（避免混合步前/步后状态）
4. 容量累计: `capacity += abs(I_current) * dt_dim / 3600.0`
5. 自适应时间步基于误差估计 (`ErrorEstimation`)

### `solve_cycle_with_export(case, cycle_opt; verbose, export_interval)`
- 执行完整的充放电循环: 放电 -> 静置1 -> 充电 -> 静置2
- 使用 `CycleOption` 控制循环参数
- 支持 verbose 输出（进度信息）
- 返回: `CycleExportData`

### `export_cycle_data_to_csv(export_data, output_dir; prefix)`
- 将循环数据导出为 6 个 CSV 文件:
  1. `{prefix}_timesteps.csv` - 时间步汇总
  2. `{prefix}_T_nodes.csv` - 节点温度场
  3. `{prefix}_soc_n.csv` - 负极 SOC 场
  4. `{prefix}_soc_p.csv` - 正极 SOC 场
  5. `{prefix}_mesh_nodes.csv` - 网格节点坐标
  6. `{prefix}_mesh_elements.csv` - 网格单元连接

### `load_cycle_data_from_csv(input_dir; prefix)`
- 从 CSV 文件加载循环数据
- 返回 `Dict` 包含:
  - 时间步数据 (times, voltages, currents, T_max, T_mean, soc_mean)
  - 温度场矩阵 (T_nodes: n_steps x nT)
  - SOC 场矩阵 (soc_n, soc_p: n_steps x ne)
  - 网格信息 (node_coords, element_connectivity)

## 变更详情

### 全部新增
此文件在 main 分支上不存在，所有 623 行都是新增代码。

## 依赖关系

### 依赖
- `CycleSolver.jl`: `PhaseType`, `PHASE_CHARGE`, `PHASE_REST`, `PHASE_DISCHARGE`, `PhaseResult`, `CycleOption`, `apply_initial_soc!`, `ErrorEstimation`
- `Initialisation.jl`: `ModelInitialisation`, `ModelInitialisation_MultiSPMe`, `MultiSPMe_get_thermal_dofs`
- `Solve.jl`: `CallModel`
- `SetCase.jl`: `Case`
- `Printf`: `@printf`, `@sprintf`
- `Statistics`: `mean`
- `LinearAlgebra`: `dot`

### 被依赖
- 被 `JuBat.jl` include 和 export
- 数据可被外部 CZM 后处理工具使用

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 是

此文件是多物理场耦合框架的**数据输出层**，主要目的是:
1. **预计算数据导出**: 将 multi-SPMe + distributed2D 仿真的温度场和 SOC 数据导出为 CSV
2. **CZM 数据准备**: 导出的温度/SOC 数据可作为 CZM 独立分析的输入载荷
3. **多 SPMe 支持**: 完全支持 per-element SPMe 模式的状态提取和导出
4. **热模型集成**: 导出完整的 2D 温度场和网格信息

关键耦合点:
- 使用 `MultiSPMe_get_thermal_dofs` 提取 multi-SPMe 热自由度
- 使用面积加权温度平均（与分布式热模型一致）
- 导出数据格式设计支持 CZM 载荷输入需求

此文件体现了"电化学-热耦合"到"力学-CZM"之间的**数据流接口**设计。
