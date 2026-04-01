# SPMe.jl

## 文件状态
修改

## main分支
- **行数**: 129
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `SPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)` | 标准 SPMe 模型求解（单全局模型） |
  | 38 | `SPMe_BC(case::Case, variables::Dict)` | SPMe 边界条件/源项 |
  | 63 | `SPMe_variables(case::Case, yt::Array{Float64}, t::Float64)` | 从状态向量提取 SPMe 物理变量 |

## Parameters_Design分支
- **行数**: 201
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `SPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)` | 标准 SPMe（已修改） |
  | 37 | `SPMe_element(case::Case, yt_e, t, e::Int; I_e, T_e, jacobi)` | 逐单元 SPMe 求解器（新增） |
  | 95 | `SPMe_BC(case::Case, variables::Dict)` | 边界条件（未变更） |
  | 120 | `SPMe_variables(case::Case, yt, t; I_app, T_e)` | 变量提取（已修改，新增外部覆盖参数） |

## 变更详情

### 新增函数

#### `SPMe_element(case::Case, yt_e::Array{Float64}, t::Float64, e::Int; I_e::Float64, T_e::Float64, jacobi::String="update")`
- **位置**: 第 37 行
- **功能**: 逐单元 SPMe 求解器，用于 multi-SPMe 架构中每个热单元的独立电化学求解
- **与 `SPMe()` 的关键区别**:
  - 接受外部电流 `I_e` 和温度 `T_e`（不从全局 `case.opt.Current` 和状态向量获取）
  - 调用 `SPMe_variables(case, yt_e; I_app=I_e, T_e=T_e)` 传入覆盖参数
  - 支持 `jacobi="constant"` 模式复用预存矩阵（`param.NE.M_d`, `param.NE.K_d`）
  - 在 `variables_e` 中添加 `"element index"` 字段
- **工作流程**:
  1. 向量化输入
  2. 调用 `SPMe_variables` 提取物理量（传入 I_e, T_e）
  3. 力学耦合（如果启用，计算应力耦合扩散系数 theta_Mn/Mp）
  4. 粒子扩散矩阵（负极/正极）
  5. 时间尺度归一化（`M_np *= ts_n / t0`, `M_pp *= ts_p / t0`）
  6. 电解液扩散矩阵
  7. 边界条件（源项）
  8. 装配局部系统矩阵 `blockdiag(M_np, M_pp, M_el)`

### 修改函数

#### `SPMe(case::Case, yt, t; jacobi)`
- **变更**: 时间尺度归一化从 `param_dim.scale.t0` 改为 `case.param_dim.scale.t0`（第 23-24 行）
  - **原**: `M_np = M_np .* param.scale.ts_n / param_dim.scale.t0`
  - **新**: `M_np = M_np .* param.scale.ts_n / case.param_dim.scale.t0`
  - 同样修改了 `M_pp` 和 `M_el`
- **影响**: 修正了变量作用域，使用 `case.param_dim` 而非未限定的 `param_dim`

#### `SPMe_variables(case::Case, yt::Array{Float64}, t::Float64; I_app=nothing, T_e=nothing)`
- **签名变更**: 新增关键字参数 `I_app::Union{Nothing,Float64}=nothing` 和 `T_e::Union{Nothing,Float64}=nothing`
- **新增逻辑**:
  - `I_app` 外部覆盖：若提供则直接使用（用于逐单元分流），否则从 `case.opt.Current` 计算
  - `T_e` 外部覆盖：若提供则直接使用（用于逐单元温度），否则从状态向量 `case.index["temperature"]` 读取
  - 当 `T_e` 被覆盖时，从 `var_list` 中过滤掉 `"temperature"` 键（避免从状态向量读取不存在的温度 DOF）
  - 当 `T_e` 为 `nothing` 时，回退到从状态向量读取温度
- **电导率计算改进**（第 162-172 行）:
  - **原**: 使用固定浓度 `param.EL.ce0` 计算 `kappa`（`kappa_ne = param.EL.kappa(param.EL.ce0, T) * eps^brugg`）
  - **新**: 使用高斯点浓度 `ce_n_gs`, `ce_p_gs`, `ce_sp_gs` 计算 `kappa`，然后积分平均
  - 新增 `IntV` 积分运算：`kappa_ne_av = IntV(kappa_ne_gs, mesh_ne) / param.NE.thickness`
  - 电阻计算使用平均电导率：`R_EL = t_NE/(3*kappa_ne_av) + t_SP/kappa_sp_av + t_PE/(3*kappa_pe_av)`
- **电解液电势公式修正**（第 173 行）:
  - **原**: `dphi_e = 2*T*(1-tplus)*(csp_av-csn_av)/ce0 .- I*R_EL .- dphi_S`（广播运算）
  - **新**: `dphi_e = 2*T*(1-tplus)*(csp_av-csn_av)/ce0 - I*R_EL - dphi_S`（标量运算，移除多余的广播点）
- **新增输出**（第 196-199 行）:
  - 在非 `distributed2D` 模式下设置 `"thermal2D element current"` 变量
  - 用于兼容单 SPMe + distributed2D 热模型路径

### 删除函数
无

## 依赖关系

### 该文件依赖哪些其他文件
- `SetCase.jl` — 使用 `Case` 结构体、`case.param`, `case.mesh`, `case.opt`, `case.index`
- `Variables.jl` — 调用 `StandardVariables()`
- `FEM.jl` (或等效) — 调用 `ElectrodeDiffusion()`, `ElectrolyteDiffusion()`, `IntV()`
- `Mechanical.jl` — 调用 `Mechanicaloutput()`（力学耦合）

### 哪些文件依赖该文件
- `Solve.jl` — `CallModel` 和 `CallModel_MultiSPMe` 调用 `SPMe()` 和 `SPMe_element()`
- `JuBat.jl` — export `SPMe`, `SPMe_element`

### 新增的外部依赖
无

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系
SPMe.jl 是耦合架构的**电化学核心**：
- **multi-SPMe**: `SPMe_element` 是逐单元求解的核心函数，接受外部指定的电流 `I_e` 和温度 `T_e`
- **distributed2D**: `SPMe_variables` 通过 `T_e` 参数接受单元温度，影响反应速率和电导率
- **CZM**: 通过 `Mechanicaloutput` 计算应力耦合扩散系数（theta_Mn/Mp），影响颗粒扩散

### 哪些变更是耦合相关的
- `SPMe_element` 整个新增函数 — multi-SPMe 架构的单元级求解器
- `SPMe_variables` 的 `I_app`/`T_e` 外部覆盖参数 — 使函数可被 `SPMe_element` 以逐单元模式调用
- `SPMe_variables` 的温度过滤逻辑（`var_list` 中移除 `"temperature"`） — 适配分布式热模型的状态向量结构
- 电导率计算改进（高斯点浓度积分平均） — 提高分流求解器精度

### 哪些变更是独立的
- 时间尺度归一化变量作用域修正（`param_dim` → `case.param_dim`） — bug 修复
- `dphi_e` 广播运算符修正 — 代码清理
- `"thermal2D element current"` 输出 — 兼容性输出
