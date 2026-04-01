# Option.jl

## 文件状态
修改 (modified)

## main分支
- 行数: 28
- 主要函数/类型列表:
  - `Option` (mutable struct, @with_kw) -- 仿真选项配置，包含基本参数（Np, Nn, 网格、时间步、求解类型等）
  - 仅含 `thermalmodel` (none/lumped/full) 和 `mechanicalmodel` (none/full) 两个物理场开关

## Parameters_Design分支
- 行数: 80 (+52 行, +186%)
- 主要函数/类型列表:
  - `CycleOption` (mutable struct, @with_kw) -- 新增：充放电循环参数
  - `PhaseType` (enum) -- 新增：循环阶段类型枚举 (CHARGE/REST/DISCHARGE)
  - `Option` (mutable struct, @with_kw) -- 大幅扩展，新增热模块、CZM模块、调试等配置

## 变更详情

### 新增函数

#### 1. `CycleOption` (struct, 行 5-25)
```julia
@with_kw mutable struct CycleOption
    n_cycles::Int64 = 50
    t_charge::Float64 = 3600.0
    t_rest1::Float64 = 600.0
    t_discharge::Float64 = 3600.0
    t_rest2::Float64 = 600.0
    I_charge::Float64 = 5.0
    I_discharge::Float64 = 5.0
    V_upper::Float64 = 4.2
    V_lower::Float64 = 2.5
    SOC_init::Float64 = 0.05
    reset_T_each_cycle::Bool = true
    reset_T_before_charge::Bool = true
    dt_cycle::Vector{Float64} = [1.0, 10.0]
end
```
充放电循环参数结构体，包含循环次数、各阶段时长、电流、截止条件、SOC初值、温度重置选项和时间步长。

#### 2. `PhaseType` (enum, 行 30-34)
```julia
@enum PhaseType begin
    PHASE_CHARGE = 1
    PHASE_REST = 2
    PHASE_DISCHARGE = 3
end
```
循环阶段类型枚举，用于标识充电、静置、放电阶段。

### 修改函数

#### `Option` (struct)

**变更内容：**

1. **thermalmodel 注释更新**：`"none"` 的可选值从 `"none, lumped or full (not implemented)"` 改为 `"none, lumped, distributed2D"`

2. **cite 类型修正**：`cite::Array{String} = []` 改为 `cite::Vector{String} = String[]`（更规范的 Julia 写法）

3. **新增热模块字段 (6 个)**：
   - `thermal_enabled::Bool = false` -- 热模块总开关
   - `thermal_dim::String = "1D"` -- 热模型维度 ("1D"/"2D")
   - `thermalmeshType::String = "L2"` -- 热网格类型 (1D: L2/L3; 2D: Q4)
   - `cool_method::String = "tab"` -- 冷却方式 ("none"/"tab"/"surface")
   - `collector_seeded::Bool = false` -- 是否使用 collector-seeded 带状网格语义
   - `per_element_spme::Bool = false` -- 逐单元 SPMe 模式开关

4. **新增调试字段 (2 个)**：
   - `debug_coupling::Bool = false` -- 电-热耦合详细日志开关
   - `debug_log_path::String = "output/debug.log"` -- 调试日志路径

5. **新增 CZM 模块字段 (8 个)**：
   - `czm_model::String = "model1"` -- CZM 模型类型 ("model1"/"mix")
   - `czm_enabled::Bool = true` -- CZM 损伤模型开关
   - `czm_update_interval::Int64 = 1` -- 损伤更新间隔（时间步数）
   - `czm_soh_threshold::Float64 = 0.8` -- SOH 终止阈值
   - `czm_inner_exit_only::Bool = true` -- 断裂时仅内圈退出电化学反应
   - `czm_iter_method::String = "load_substep"` -- CZM 迭代方法 ("basic"/"load_substep"/"arc_length")
   - `czm_max_iter::Int64 = 100` -- 牛顿迭代最大步数
   - `czm_tol::Float64 = 1e-4` -- 收敛容差
   - `czm_load_steps::Int64 = 20` -- 载荷子步数
   - `czm_arc_length_alpha::Float64 = 1.0` -- 弧长法系数

### 删除函数
无删除函数。

## 依赖关系

### 该文件依赖哪些其他文件
- `Parameters.jl` (外部包) -- 提供 `@with_kw` 宏

### 哪些文件依赖该文件
- 几乎所有 `src/` 下的文件都通过 `case.opt` 访问 Option 字段，具体包括：
  - `src/Solve.jl` -- 检查 `thermalmodel`, `per_element_spme`, `czm_enabled` 等
  - `src/Variables.jl` -- 根据 `opt.thermalmodel`, `opt.czm_enabled` 分配变量
  - `src/ThermalDistributed.jl` -- 使用 `opt.cool_method`, `opt.czm_enabled`
  - `src/CycleSolver.jl` -- 使用 `CycleOption`, `PhaseType`
  - `src/CycleData.jl` -- 使用 `CycleOption`, `PhaseType`
  - `src/PostProcessing.jl` -- 使用 `PhaseType`
  - `src/JuBat.jl` -- 导出 `CycleOption`, `PhaseType`

### 新增的外部依赖
无新增外部包依赖（`Parameters.jl` 已在 main 中使用）。

## 耦合分析

### 与 multi-SPMe + distributed2D + CZM 耦合的关系
- **配置中枢**：Option.jl 定义了所有耦合模块的开关和参数，是整个多物理场耦合架构的控制面板。
- `per_element_spme` 是 multi-SPMe 架构的核心开关，决定了是否启用逐单元电化学模型。
- `thermalmodel = "distributed2D"` 配合 `per_element_spme = true` 实现电-热耦合。
- `czm_enabled` 和相关 CZM 参数控制热-力耦合。

### 哪些变更是耦合相关的
- 热模块字段 (6 个)：`thermal_enabled`, `thermal_dim`, `thermalmeshType`, `cool_method`, `collector_seeded`, `per_element_spme` -- 全部为耦合服务
- CZM 模块字段 (8 个)：全部为热-力耦合服务
- 调试字段 (2 个)：用于调试耦合行为
- `thermalmodel` 注释更新：反映 distributed2D 新模型

### 哪些变更是独立的
- `CycleOption` 和 `PhaseType` -- 循环仿真功能，独立于耦合本身（但循环会使用耦合）
- `cite` 类型修正 -- 纯代码风格改进
