# Option.jl

- **源文件**: `src/Option.jl`
- **行数**: 87 行
- **函数/struct 计数**: 0 个独立函数；2 个 struct（`CycleOption`、`Option`）；1 个 enum（`PhaseType`）
- **职责**: 定义仿真求解器配置 `Option`、循环仿真配置 `CycleOption`、循环阶段枚举 `PhaseType`
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/09_分流求解器.md`

## 数据结构

### `CycleOption`（L5-L25，`@with_kw mutable struct`）

充放电循环仿真的顶层配置。

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `n_cycles` | Int64 | 50 | 循环次数 |
| `t_charge` | Float64 | 3600.0 | 充电时长 [s] |
| `t_rest1` | Float64 | 600.0 | 充电后静置 [s] |
| `t_discharge` | Float64 | 3600.0 | 放电时长 [s] |
| `t_rest2` | Float64 | 600.0 | 放电后静置 [s] |
| `I_charge` | Float64 | 5.0 | 充电电流 [A]（正值） |
| `I_discharge` | Float64 | 5.0 | 放电电流 [A]（正值） |
| `V_upper` | Float64 | 4.2 | 充电截止电压 [V] |
| `V_lower` | Float64 | 2.5 | 放电截止电压 [V] |
| `SOC_init` | Float64 | 0.05 | 初始 SOC |
| `reset_T_each_cycle` | Bool | true | 每循环重置温度 |
| `reset_T_before_charge` | Bool | true | 循环内充电前重置温度 |
| `dt_cycle` | Vector{Float64} | [1.0, 10.0] | [dt_min, dt_max] |

### `PhaseType`（L30-L34，`@enum`）

循环阶段类型枚举：

- `PHASE_CHARGE = 1`：充电（I < 0）
- `PHASE_REST = 2`：静置（I = 0）
- `PHASE_DISCHARGE = 3`：放电（I > 0）

### `Option`（L35-L87，`@with_kw mutable struct`）

电池模型求解器的主配置结构。所有字段按源文件出现顺序列出。

#### 网格与离散化

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `Np` | Int64 | 10 | 正极 electrolyte 网格点数 |
| `Ns` | Int64 | 10 | 隔膜 electrolyte 网格点数 |
| `Nn` | Int64 | 10 | 负极 electrolyte 网格点数 |
| `Nrp` | Int64 | 10 | 正极颗粒径向网格点数 |
| `Nrn` | Int64 | 10 | 负极颗粒径向网格点数 |
| `meshType` | String | "L2" | 网格单元类型 |
| `gsorder` | Int64 | 4 | 高斯积分阶数 |
| `dimension` | Int64 | 1 | 维度 |

#### 模型与时间

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `model` | String | "SPM" | 模型类型（"SPM" / "SPMe" / "P2D"） |
| `time` | Array{Float64} | [0 3600] | 仿真时间区间 [s] |
| `Current` | Function | x->0 | 电流函数 I(t) [A] |
| `coupleMethod` | String | "fully coupled" | 耦合方法 |
| `coupleOrder` | Int64 | 0 | 耦合阶数 |
| `y0` | Array{Float64} | [] | 初始状态向量 |
| `dt` | Array{Float64} | [1, 100] | 时间步长范围 |
| `dtType` | String | "constant" | 时间步类型（"auto"/"manual"） |
| `dtThreshold` | Float64 | 0.01 | 时间步阈值 |
| `solveType` | String | "Crank-Nicolson" | 时间离散格式（"forward"/"backward"/"Crank-Nicolson"） |
| `outputType` | String | "auto" | 输出类型 |
| `jacobi` | String | "constant" | 雅可比策略（"constant"/"update"） |

#### 多物理场开关

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `thermalmodel` | String | "none" | 热模型（"none"/"lumped"/"distributed2D"） |
| `mechanicalmodel` | String | "none" | 力学模型（"none"/"full"） |
| `cite` | Vector{String} | String[] | 引用文献 |

#### 热学模块

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `thermal_enabled` | Bool | false | 是否启用热模块 |
| `thermal_dim` | String | "1D" | "1D" 或 "2D" |
| `thermalmeshType` | String | "L2" | 热网格类型（1D: L2/L3；2D: Q4） |
| `cool_method` | String | "tab" | 冷却方式（"none"/"tab"/"surface"） |
| `collector_seeded` | Bool | false | 启用 collector-seeded 带状网格语义（layer_weights） |
| `per_element_spme` | Bool | false | 允许逐单元传递 I_app 与 T 给 SPMe |
| `czm_model` | String | "model1" | CZM 模型（"model1"/"mix"） |

#### 调试

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `debug_coupling` | Bool | false | 打印电-热耦合详细日志 |
| `debug_log_path` | String | "output/debug.log" | 调试日志路径 |

#### CZM 模块

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `czm_enabled` | Bool | false | 启用 CZM 损伤模型 |
| `czm_update_interval` | Int64 | 1 | 损伤更新间隔（时间步数） |
| `czm_soh_threshold` | Float64 | 0.8 | SOH 终止阈值 |
| `czm_inner_exit_only` | Bool | true | 断裂时仅内圈单元退出电化学 |
| `czm_fix_inner` | Bool | true | 边界条件：是否固定内圈节点 |
| `czm_iter_method` | String | "basic" | 迭代方法（"basic"/"load_substep"/"arc_length"） |
| `czm_max_iter` | Int64 | 100 | 牛顿迭代最大步数 |
| `czm_tol` | Float64 | 1e-4 | 收敛容差 |
| `czm_load_steps` | Int64 | 2 | 载荷子步数（load_substep 模式） |
| `czm_arc_length_alpha` | Float64 | 1.0 | 弧长法系数（arc_length 模式） |
| `czm_viscous_enabled` | Bool | false | 粘性正则化开关 |
| `czm_visc_tau` | Float64 | 0.0 | 物理松弛时间 [s]（推荐 10~100） |
| `czm_area_loss_enabled` | Bool | false | 启用渐进式面积损失 |
| `czm_area_loss_threshold` | Float64 | 0.83 | 面积开始缩减的损伤阈值 |

## 函数清单

本文件无独立函数定义（仅 struct 与 enum 声明）。

## 省略项

- `using Parameters`（L1）：导入 `@with_kw` 宏
- 各字段注释（如 `# 循环次数`、`# 充电时长`）不计入条目

### [DEBUG]

无

> 注：`debug_coupling` (L69) 和 `debug_log_path` (L70) 是用户配置字段（非临时调试输出），不纳入 [DEBUG] 标注。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L48 | `Current::Function = x-> 0` | 默认电流为 0；用户必须显式设置否则仿真无意义。无兜底注释，但默认值本质上是占位 |
| L51 | `y0::Array{Float64} = []` | 空数组作为初始状态占位；由 `ModelInitialisation` 填充 |

### [COMPLEX-CHECK]

无
