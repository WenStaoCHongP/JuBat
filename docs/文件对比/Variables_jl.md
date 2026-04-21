# Variables.jl

## 文件状态
修改 (modified)

## main分支
- 行数: 105
- 主要函数列表:
  - `StandardVariables(case::Case, num::Int64)` -- 创建标准变量字典，根据模型类型 (SPM/SPMe/P2D/sP2D) 分配不同变量
  - `Variable_update!(variables_hist, variables, v)` -- 将第 v 步变量值写入历史记录

## Parameters_Design分支
- 行数: 181 (+76 行, +72%)
- 主要函数列表:
  - `StandardVariables(case::Case, num::Int64)` -- 扩展版，新增 distributed2D 热变量和 CZM 变量
  - `Variable_update!(variables_hist, variables, v)` -- 完全重写，增加动态数组扩展和更健壮的赋值逻辑

## 变更详情

### 新增函数
无新增函数，两个原有函数均有修改。

### 修改函数

#### 1. `StandardVariables(case::Case, num::Int64)`

**变更内容：**

1. **移除 sP2D 模型支持**：将 `"P2D" || "sP2D"` 条件统一为 `"P2D"`，在三个位置移除了对 `"sP2D"` 模型的判断（行 9, 52, 69）。

2. **温度变量无条件创建**：
   - main: `if "temperature" in collect(keys(case.index))` 条件判断后创建，否则 fallback 为 `zeros(1, num)`
   - HEAD: 直接 `variables["temperature"] = zeros(Float64, length(case.index["temperature"]), num)`，不再做条件判断

3. **新增 lumped 热变量** (thermalmodel == "lumped"):
   - `"thermal lumped internal heat"` -- 集总模型内部产热

4. **新增 distributed2D 热变量块** (thermalmodel == "distributed2D")，共 31 个新变量：
   - 节点温度场：`"T_nodes"`
   - 总热源：`"heat_source_fields"`, `"total heat source"`
   - 分层热源 (11 个)：`"thermal2D q_rxn_ne"`, `"thermal2D q_rev_ne"`, `"thermal2D q_ohm_s_ne"`, `"thermal2D q_ohm_e_ne"`, `"thermal2D q_sp"`, `"thermal2D q_rxn_pe"`, `"thermal2D q_rev_pe"`, `"thermal2D q_ohm_s_pe"`, `"thermal2D q_ohm_e_pe"`, `"thermal2D q_pcc"`, `"thermal2D q_ncc"`
   - 单元电化学变量 (8 个)：`"thermal2D element current"`, `"thermal2D eta_n_e"`, `"thermal2D eta_p_e"`, `"thermal2D dUdT_n_e"`, `"thermal2D dUdT_p_e"`, `"thermal2D element soc_n"`, `"thermal2D element soc_p"`, `"thermal2D element voltages"`, `"thermal2D element OCV"`
   - 活跃掩码与截止检测：`"thermal2D active_mask"`, `"thermal2D n_cutoff_elements"`, `"thermal2D nearest_cutoff_element"`, `"thermal2D nearest_cutoff_ocv"`, `"thermal2D margin_to_cutoff"`
   - 应力/应变：`"thermal2D element thermal stress"`, `"thermal2D element diffusion stress"`, `"thermal2D element total stress"`, `"thermal2D element diffusion strain"`, `"thermal2D element thermal strain"`
   - 位移场：`"thermal2D displacement x"`, `"thermal2D displacement y"`

5. **新增 CZM 变量** (czm_enabled == true):
   - `"negative electrode cohesive zone damage"` -- 负极内聚力损伤
   - `"positive electrode cohesive zone damage"` -- 正极内聚力损伤

#### 2. `Variable_update!(variables_hist, variables, v)`

**变更内容：** 完全重写，从简单的 `for` 循环改为更健壮的实现：

1. **动态数组扩展**：当时间步 `v` 超过预分配列数时，自动扩展 `expansion_size = max(1000, current_size / 2)` 列，并发出 `@warn` 警告。
2. **类型安全赋值**：根据 `hist_val` 和 `val` 的具体类型 (`Array{Float64}` vs `Float64`)，`ndims` 维度等做精确匹配赋值。
3. **避免 KeyError**：使用 `hist_keys = Set(keys(variables_hist))` 和 `k in hist_keys || continue` 避免向历史字典中写入新键。

### 删除函数
无删除函数。

## 依赖关系

### 该文件依赖哪些其他文件
- `src/Option.jl` -- 通过 `case.opt.model`, `case.opt.thermalmodel`, `case.opt.czm_enabled` 等字段决定分配哪些变量
- `src/SetCase.jl` -- 通过 `case::Case` 类型获取 `case.mesh`, `case.index`, `case.param` 等
- `src/Jellyrollmodel.jl` -- 通过 `case.mesh["thermal2D"]` 访问二维热网格（distributed2D 变量块）

### 哪些文件依赖该文件
- `src/Solve.jl` -- 调用 `StandardVariables` 创建变量，`Variable_update!` 更新历史
- `src/SPM.jl` -- 调用 `StandardVariables`
- `src/SPMe.jl` -- 调用 `StandardVariables`
- `src/P2D.jl` -- 调用 `StandardVariables`

### 新增的外部依赖
无新增外部包依赖。

## 耦合分析

### 与 multi-SPMe + distributed2D + CZM 耦合的关系
- **核心枢纽文件**：Variables.jl 是多物理场耦合的数据注册中心，所有耦合数据（热源、温度、损伤、应力等）都通过此文件的变量字典传递。
- distributed2D 的 31 个新变量直接支持逐单元 SPMe 架构下的电-热-力耦合数据存储。
- CZM 损伤变量 `cohesive zone damage` 为循环仿真中的退化计算提供存储。

### 哪些变更是耦合相关的
- distributed2D 变量块（全部 31 个变量）-- 电-热耦合
- CZM 变量（2 个损伤变量）-- 热-力耦合
- `Variable_update!` 的动态扩展 -- 为长时间循环仿真（耦合场景）提供支持
- 移除 sP2D -- 简化模型分支，聚焦于 SPMe + distributed2D 路径

### 哪些变更是独立的
- `Variable_update!` 的类型安全改进属于通用健壮性增强，独立于耦合需求
- 温度变量无条件创建属于代码简化，独立于耦合需求

## 后续变更 (2026-04-01)

- `T_nodes` 变量初始化从单一变量拆分为 4 个独立变量：
  - **旧** (1 行): `variables["T_nodes"] = zeros(Float64, nT, num)`
  - **新** (4 行):
    - `variables["thermal2D temperature"]` — 单元温度
    - `variables["thermal2D temperature at nodes"]` — 节点温度
    - `variables["thermal2D temperature history"]` — 单元温度历史
    - `variables["thermal2D temperature at nodes history"]` — 节点温度历史
- 此变更提供单元温度和节点温度的独立追踪，支持更精细的温度场分析
- 文件行数从约 181 行增加到约 185 行

## 后续变更 (2026-04-07)

- 无逻辑变更，文件行数约 184 行

## 后续变更 (2026-04-20)

- **新增 CZM 摘要统计变量**: `StandardVariables` 中 CZM 块新增 5 个变量：
  - `"czm D_max"` — 最大损伤值 (Float64, 1×num)
  - `"czm D_mean"` — 平均损伤值 (Float64, 1×num)
  - `"czm δ_max_n"` — 最大法向分离位移 (Float64, 1×num)
  - `"czm δ_mean_n"` — 平均法向分离位移 (Float64, 1×num)
  - `"czm n_fractured"` — 完全断裂单元数 (Float64, 1×num)
- 这些变量在 `CallModel_MultiSPMe` 中每步填充，用于 `PostProcessing` 输出和循环仿真监控
- 行数从约 184 行增加到约 191 行
