# SetParams.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 319
- 主要结构列表:
  - `@with_kw mutable struct Electrode` - 电极参数（电化学 + 基础力学: E, nu, Omega）
  - `@with_kw mutable struct Separator` - 隔膜参数
  - `@with_kw mutable struct CurrentCollector` - 集流体参数
  - `@with_kw mutable struct Electrolyte` - 电解液参数
  - `@with_kw mutable struct Cell` - 电池级参数
  - `@with_kw mutable struct Tab` - 极耳参数（length, width, area）
  - `@with_kw mutable struct Binder` - 粘结剂参数
  - `@with_kw mutable struct Scale` - 归一化尺度（9个字段）
  - `@with_kw mutable struct Params` - 参数集合（不含 cohesive 字段）
- 函数列表:
  - `ChooseCell()` - 选择电池参数（支持 LG M50, Northrop, Enertech）
  - `NormaliseParam()` - 参数归一化

## Parameters_Design分支
- 行数: 440 (+121, +38%)
- 新增结构: `Cohesive`, 扩展 `Electrode`, `Cell`, `Tab`, `Scale`, `Params`

## 变更详情

### 新增结构体

#### `@with_kw mutable struct Cohesive` (25个字段)
**法向参数 (Mode I)**:
- `sigma_max_n` - 最大法向牵引力 [Pa]
- `delta_0_n` - 损伤起始分离位移 [m]
- `delta_c_n` - 临界（完全断裂）分离位移 [m]
- `G_c_n` - 法向断裂能 [J/m^2]
- `K_n` - 法向初始刚度 [Pa/m]

**切向参数 (Mode II)**:
- `tau_max_t` - 最大切向牵引力 [Pa]
- `delta_0_t`, `delta_c_t`, `G_c_t`, `K_t` - 对应切向参数

**混合模式参数**:
- `eta` - BK准则指数 [-]
- `czm_model` - 模型选择 ("model1" 或 "mix")

**界面热阻参数**:
- `h_c0` - 完好接触换热系数 [W/(m^2*K)]
- `k_air` - 空气热导率 [W/(m*K)]
- `lambda_m` - 平均自由程 [m]
- `beta` - 参数 [-]
- `threshold` - 阈值厚度 [m]

### 修改结构体

#### `Electrode`
新增字段:
- `alphaT::Float64` - 热膨胀系数 [1/K]（**新增**，原无）
- 字段重新排序: `E`, `nu`, `alphaT`, `Omega` 归入"力学/热学属性"分组
- 原 `Omega`, `nu`, `E` 散落在其他位置，现统一到物理意义相近的位置

#### `Cell`
新增字段 (Jellyroll 几何):
- `Rin::Float64` - 内半径
- `Rout::Float64` - 外半径
- `layer::Float64` - 层厚
- `lambda_r::Float64` - 径向导热率
- `lambda_t::Float64` - 周向导热率
- `Nr_th::Int` - 热网格径向分数
- `Ntheta_th::Int` - 热网格周向分数
- `n_windings::Int` - 卷绕圈数

#### `Tab`
新增字段:
- `h::Float64` - 换热系数
- `theta_pos::Vector{Float64}` - 正极耳角度位置
- `theta_neg::Vector{Float64}` - 负极耳角度位置

#### `Scale`
新增字段 (统一能量尺度 + CZM):
- `rho::Float64` - 密度尺度 [kg/m^3]
- `P_ref::Float64` - 功率参考值 = phi * I_typ
- `lambda::Float64` - 导热率尺度 = P_ref/(L*T_ref)
- `q::Float64` - 热源尺度 = P_ref/L^3
- `h::Float64` - Biot 数 = h_cell*L/lambda_r
- `sigma_czm::Float64` - CZM 牵引力参考值
- `delta_czm::Float64` - CZM 分离位移参考值
- `G_czm::Float64` - CZM 断裂能参考值
- `K_czm::Float64` - CZM 刚度参考值

#### `Params`
新增字段:
- `cohesive::Cohesive` - 内聚力模型参数（默认空 Cohesive()）

### 修改函数

#### `ChooseCell()`
- 使用 `joinpath(@__DIR__, "parameters")` 替代硬编码路径（更健壮）
- 新增 `"Jellyroll"` 和 `"Ring"` 电池类型支持
- 修改 `cell.rho` 和 `cell.heat_Q` 计算: 添加 else 分支保留原值（对 Jellyroll 等自定义参数集）
- 新增 Scale 字段计算: `rho`, `P_ref`, `lambda`, `h`, `q`, CZM 尺度 (`sigma_czm`, `delta_czm`, `G_czm`, `K_czm`)

#### `NormaliseParam()`
**重要归一化方案变更**:

1. **力学参数归一化调整**:
   - `PE.E`, `PE.nu`, `PE.Omega` 字段赋值位置调整（不影响值）
   - 新增 `alphaT` 归一化: `param.PE.alphaT = param_dim.PE.alphaT * scale.T_ref`

2. **新增热学参数归一化**（每层）:
   - `lambda` (导热率): `param.X.lambda = param_dim.X.lambda / scale.lambda`
   - `rho` (密度): `param.X.rho = param_dim.X.rho / scale.rho`
   - `heat_Q` (体积热容): `param.X.heat_Q = ... * scale.rho * scale.L^3 * scale.T_ref / (scale.t0 * scale.phi * scale.I_typ)`
   - 对 PE, NE, SP, PCC, NCC 五层都添加了归一化

3. **Cell 归一化变更**:
   - `area`: 从 `area * phi * I_typ / capacity` 改为 `area / area`（面积比，无量纲化方式改变）
   - `volume`: 归一化方式改变
   - 新增: `layer`, `lambda_r`, `lambda_t`, `Rin`, `Rout`, `width` 归一化
   - 新增: `rho` 归一化

4. **Tab 归一化**: 新增 `length`, `width`, `area`, `h` 的归一化

5. **CZM 归一化**:
   - 法向: `sigma_max_n`, `delta_0_n`, `delta_c_n`, `G_c_n`, `K_n`
   - 切向: `tau_max_t`, `delta_0_t`, `delta_c_t`, `G_c_t`, `K_t`
   - 界面热阻: `h_c0`, `k_air`, `lambda_m`, `beta`, `threshold`

6. **Base.invokelatest 包装**:
   - 所有闭包参数（`U`, `dUdT`, `De`, `kappa`）改用 `Base.invokelatest` 包装
   - 目的: 避免 Julia world-age 问题（当参数文件通过 `include` 加载时闭包可能不在同一个 world age）

7. **Cell 归一化方式重构**:
   - `area` 归一化: 旧 `area * phi * I_typ / capacity` -> 新 `area / area`（纯归一化，不包含物理量转换）
   - `volume` 归一化: 旧 `volume * phi / L * I_typ / capacity` -> 新 `volume / L^3`

## 依赖关系

### 被依赖关系
- 所有参数结构被 `NormaliseParam` 和参数文件使用
- `Cohesive` 被 `czm.jl`, `CzmSolve.jl` 使用
- `Scale` 的新字段被 `ThermalDistributed.jl`, `ThermalPolar2D.jl` 使用
- `Cell` 的 Jellyroll 字段被 `Jellyrollmodel.jl` 使用
- `Params` 是整个框架的核心数据结构

### 依赖关系
- 依赖 `@with_kw` 宏 (Parameters.jl)
- `ChooseCell` 通过 `include` 加载参数文件

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 是

此文件是参数传递和归一化的核心，为多物理场耦合提供统一的参数框架：
- `Cohesive` 结构支持 CZM 模型
- `Scale` 的统一能量尺度（`P_ref`, `lambda`, `q`）支持分布式热模型
- `Cell` 的 Jellyroll 几何参数支持螺旋网格
- `Tab` 的角度位置支持极耳冷却边界条件
- 热学参数归一化支持分层热源计算
- `Base.invokelatest` 修复确保参数闭包在多 SPMe 场景下正确工作

关键设计变更: 从"以电化学为主"的归一化方案转向"统一能量尺度"方案，使电化学和热模型使用一致的时间尺度和功率参考。
