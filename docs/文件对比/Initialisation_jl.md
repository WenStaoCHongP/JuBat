# Initialisation.jl

## 文件状态
修改

## main分支
- **行数**: 43
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `ModelInitialisation(case::Case)` | 标准 SPMe/SPM/P2D 模型初始化，构建初始状态向量 y0 |

## Parameters_Design分支
- **行数**: 158
- **主要函数列表**:
  | 行号 | 函数签名 | 说明 |
  |------|----------|------|
  | 1 | `ModelInitialisation(case::Case)` | 标准初始化（已修改，新增 distributed2D 支持） |
  | 54 | `ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution=nothing)` | 多SPMe架构初始化，构建扩展状态向量 |
  | 121 | `extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)` | 从全局状态向量中提取单个单元的电化学状态 |
  | 132 | `get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)` | 从全局状态向量中提取热场节点温度 |
  | 142 | `update_state(y::AbstractVector, layout::MultiSPMeLayout; element_index, element_state, thermal_nodes)` | 更新全局状态向量（返回新向量） |

## 变更详情

### 新增函数

#### `ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution=nothing)`
- **位置**: 第 54 行
- **功能**: 为多 SPMe 并行架构初始化扩展状态向量
- **工作流程**:
  1. 获取热网格单元数 `ne` 和节点数 `nT`
  2. 临时将 `thermalmodel` 设为 `"none"`，调用标准 `ModelInitialisation` 获取单个单元的纯电化学初始状态
  3. 为每个单元复制/定制电化学初始状态（支持非均匀 SOC 分布）
  4. 追加热场初始温度 `T0_nodes`
  5. 组装全局状态向量 `y0 = [y0_chem_all; T0_nodes]`
  6. 构造 `case.layout = MultiSPMeLayout(ne, n_chem, nT)`（类型化替代 Dict 填充）

#### `extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout) -> Vector{Float64}`
- **位置**: 第 121 行
- **功能**: 从全局状态向量中提取第 `e` 个单元的局部电化学状态（cn_surf; cp_surf; ce）
- **设计**: 参数从 `case::Case` 改为 `layout::MultiSPMeLayout`（类型安全）

#### `get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout) -> Vector{Float64}`
- **位置**: 第 132 行
- **功能**: 从全局状态向量中提取热场节点温度 DOF
- **使用**: 在 `CallModel_MultiSPMe` 中每步调用

#### `update_state(y::AbstractVector, layout::MultiSPMeLayout; element_index, element_state, thermal_nodes)`
- **位置**: 第 142 行
- **功能**: 创建更新后的全局状态向量副本（非原地修改）
- **验证**: 包含元素索引范围检查、状态长度检查（assert-based）

### 修改函数

#### `ModelInitialisation(case::Case)`
- **变更 1**: 移除了 `sP2D` 模型分支（第 14 行 `P2D || sP2D` 变为仅 `P2D`）
- **变更 2**: 修正 P2D 模式中 `phis_p` 和 `phis_n` 的尺寸（原先是 `Np`/`Nn` 混用，改为 `Nn`/`Np`）
- **变更 3**: 新增 `distributed2D` 热模型初始化分支（第 27-30 行）：
  - 获取 `thermal2D` 网格节点数 `nT`
  - 用初始温度 `T0` 填充节点温度
  - 追加到状态向量 `y0`
- **变更 4**: 移除最后的 `sP2D` 条件判断（第 41 行）

### 删除函数
无

## 依赖关系

### 该文件依赖哪些其他文件
- `SetCase.jl` — 使用 `Case` 结构体、`case.mesh`、`case.param`、`case.opt`、`case.index`
- `Variables.jl` — 可能依赖标准变量定义

### 哪些文件依赖该文件
- `Solve.jl` — 调用 `ModelInitialisation()` 和 `ModelInitialisation_MultiSPMe()`
- `CallModel.jl` — 调用 `extract_element_state()` 和 `get_thermal_dofs()`
- `CycleData.jl` — 调用 `ModelInitialisation_MultiSPMe()` 和 `get_thermal_dofs()`
- `JuBat.jl` — export `ModelInitialisation`, `ModelInitialisation_MultiSPMe`, `extract_element_state`, `get_thermal_dofs`, `update_state`

### 新增的外部依赖
无

## 耦合分析

### 该文件与 multi-SPMe + distributed2D + CZM 耦合的关系
Initialisation.jl 是耦合架构的**数据准备层**：
- **multi-SPMe**: `ModelInitialisation_MultiSPMe` 构建多单元并行的扩展状态向量，并构造 `case.layout`（`MultiSPMeLayout`）
- **distributed2D**: 在标准初始化中新增 `distributed2D` 分支，将热节点温度 DOF 追加到状态向量
- **CZM**: 非均匀 SOC 分布的初始化接口（`initial_soc_distribution` 参数）为后续 CZM 仿真提供了预设损伤区域的入口

### 哪些变更是耦合相关的
- `ModelInitialisation_MultiSPMe` 整个新增函数 — multi-SPMe + distributed2D 的核心初始化
- `ModelInitialisation` 中 `distributed2D` 分支（第 27-30 行）— 热耦合初始化
- `extract_element_state` / `get_thermal_dofs` / `update_state` — 多 SPMe 状态管理的工具函数（参数使用 `MultiSPMeLayout` 替代 `Case`）

### 哪些变更是独立的
- P2D 模型中 `phis_p`/`phis_n` 尺寸修正（第 24-25 行）— bug 修复，与耦合无关
- 移除 `sP2D` 模型分支 — 模型简化，与耦合无关

## 后续变更 (2026-04-07)

- **函数重命名**: `MultiSPMe_extract_element_state` → `extract_element_state`, `MultiSPMe_get_thermal_dofs` → `get_thermal_dofs`, `MultiSPMe_update_state` → `update_state`
- **参数类型变更**: 辅助函数的参数从 `case::Case` 改为 `layout::MultiSPMeLayout`（类型安全）
- **布局缓存**: `case.multi_spme_layout[...] = ...` Dict 填充 → `case.layout = MultiSPMeLayout(ne, n_chem, nT)` 构造器
- 行数从约 239 行减少到 158 行（简化了冗余的类型检查和 Dict 访问代码）
