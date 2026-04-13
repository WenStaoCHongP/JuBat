# ThermalDistributed.jl

## 文件状态
新增 (new file)

## main分支
- 文件不存在于 main 分支

## Parameters_Design分支
- 行数: 425

## 后续变更 (2026-04-07)

- **P3 耦合泄漏修复**: `jellyroll_element_properties` 直接调用 → `case.geometry.layer_weights`（带 fallback）
  - `ThermalDistributed2D` 和 `ThermalDistributed2D_Ring` 两处均改为 `case.geometry !== nothing ? case.geometry.layer_weights : jellyroll_element_properties(...)[2]`
  - 消除 ThermalDistributed→Jellyrollmodel 的硬依赖（可选耦合）
- 主要函数列表:
  1. `ThermalDistributed2D(case, variables)` -- 2D 分布式热模型 FEM 矩阵组装（直角坐标）
  2. `apply_convection_bc(KT, FT, mesh, is_outer, case)` -- 对流边界条件施加（Biot 数）
  3. `apply_cool_method(KT, FT, mesh, case)` -- 冷却方式处理（none/tab/surface）
  4. `ThermalDistributed2D_BC(KT, FT, case, t)` -- 热模型边界条件入口（含 CZM 界面热阻）
  5. `ThermalDistributed2D_Ring(case, variables)` -- 环形极坐标 FEM 热模型（极坐标各向异性导热）
  6. `ThermalRing2D_BC(KT, FT, case, outer_nodes, t)` -- 环形热模型边界条件
  7. `compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme)` -- 分层热源计算
  8. `compute_heat_sources_with_czm(case, variables, variables_elems, I_e, T_e, areas, czm_mesh, mesh_data)` -- 含 CZM 损伤屏蔽的热源计算

## 变更详情

### 新增函数

#### 1. `ThermalDistributed2D(case::Case, variables)` (行 1-47)
- **功能**：基于 FEM 组装二维分布式热模型的质量矩阵 MT、刚度矩阵 KT 和载荷向量 FT
- **实现要点**：
  - 使用 `jellyroll_element_properties` 获取各单元的层权重
  - 通过 `thermal_capacity_weights_2d` 和 `thermal_anisotropic_conductivity_2d` 获取各向异性的热物性
  - 刚度矩阵由 4 个子矩阵 KT_xx, KT_xy, KT_yx, KT_yy 组合（各向异性）
  - 载荷向量由 `heat_source_fields` 映射到节点
  - 全部使用无量纲参数

#### 2. `apply_convection_bc(KT, FT, mesh, is_outer, case)` (行 49-105)
- **功能**：对边界节点施加对流边界条件（Newton 冷却定律）
- **实现要点**：
  - 使用 Biot 数 `Bi = h * lambda_r` 作为无量纲对流系数
  - 遍历所有外边界边，使用 2 点高斯积分
  - 修改刚度矩阵 K 和载荷向量 F
  - 防止重复处理同一条边（使用 `seen` 集合）

#### 3. `apply_cool_method(KT, FT, mesh, case)` (行 107-185)
- **功能**：根据冷却方式修改热矩阵
- **实现要点**：
  - `"none"`：无冷却，直接返回
  - `"surface"`：表面冷却，使用 `conv_factor = 2 * Bi / width` 修改全局矩阵
  - `"tab"`：极耳冷却，通过 `jellyroll_tab_node_indices` 定位极耳节点，按弧长加权施加冷却条件

#### 4. `ThermalDistributed2D_BC(KT, FT, case, t)` (行 187-219)
- **功能**：热模型边界条件入口函数，整合 CZM 界面热阻和对流冷却
- **实现要点**：
  - 如果 `czm_enabled` 且存在 `czm_mesh`，遍历所有内聚力单元，使用 `compute_gap_conductance` 计算损伤相关的有效热导
  - 在 cohesive 上下节点之间添加热耦合项
  - 调用 `apply_convection_bc` 和 `apply_cool_method`

#### 5. `ThermalDistributed2D_Ring(case, variables)` (行 221-270)
- **功能**：极坐标下的环形热模型 FEM 组装
- **实现要点**：
  - 将直角坐标梯度转换为极坐标 (dN/dr, dN/dtheta)
  - 使用各向异性导热率 `lambda_r` 和 `lambda_t`
  - 径向和切向独立组装后相加

#### 6. `ThermalRing2D_BC(KT, FT, case, outer_nodes, t)` (行 272-279)
- **功能**：环形热模型边界条件，对外边界施加对流

#### 7. `compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme)` (行 281-391)
- **功能**：逐单元计算分层热源（反应热 + 可逆热 + 欧姆热）
- **实现要点**：
  - 支持两种模式：`per_element_spme=true` 从 `variables_elems[e]` 读取单元级变量；否则从全局 `variables` 读取
  - 计算每层热源：NE (反应/可逆/固体欧姆/电解液欧姆), SP (欧姆), PE (反应/可逆/固体欧姆/电解液欧姆), PCC/NCC (集流体欧姆)
  - 按层权重 `fks[e,:]` 分配到单元
  - 温度相关电导率计算
  - 最终通过 `scale.L^3 / volume` 将热源转换为适当量纲
- **参数**：
  - `case` -- Case 结构
  - `variables` -- 全局变量字典
  - `variables_elems` -- 逐单元变量向量（可 Nothing）
  - `I_e` -- 单元电流向量
  - `T_e` -- 单元温度向量
  - `areas` -- 单元面积向量
  - `per_element_spme` -- 是否使用逐单元 SPMe

#### 8. `compute_heat_sources_with_czm(case, variables, variables_elems, I_e, T_e, areas, czm_mesh, mesh_data)` (行 393-425)
- **功能**：在 `compute_heat_sources` 基础上叠加 CZM 损伤屏蔽效果
- **实现要点**：
  - 先调用 `compute_heat_sources` 计算所有单元热源
  - 通过 `get_active_elements` 获取活跃单元
  - 将断裂（非活跃）单元的热源设为 0
  - 更新总功率（仅活跃单元）

### 修改函数
不适用（新文件）。

### 删除函数
不适用（新文件）。

## 依赖关系

### 该文件依赖哪些其他文件
- `src/Option.jl` -- 通过 `case.opt.cool_method`, `case.opt.czm_enabled` 等
- `src/SetCase.jl` -- `Case` 类型定义
- `src/Assemble.jl` -- `Assemble`, `Assemble1D` 矩阵组装函数
- `src/Jellyrollmodel.jl` -- `jellyroll_element_properties`, `jellyroll_tab_node_indices`, `identify_boundary_nodes`
- `src/czm.jl` / `src/CzmSolve.jl` -- `compute_gap_conductance`, CZM mesh 数据结构
- `src/Thermal.jl` (间接) -- 共享热模型接口概念

### 哪些文件依赖该文件
- `src/JuBat.jl` -- 导出所有公开函数
- `src/Solve.jl` -- 主求解器，调用 `ThermalDistributed2D`, `ThermalDistributed2D_BC`, `ThermalDistributed2D_Ring`, `ThermalRing2D_BC`, `compute_heat_sources`, `compute_heat_sources_with_czm`

### 新增的外部依赖
无新增外部包依赖。使用标准库 `SparseArrays`（可能通过 Assemble.jl 间接使用）。

## 耦合分析

### 与 multi-SPMe + distributed2D + CZM 耦合的关系
- **耦合实现的核心文件**：ThermalDistributed.jl 是整个电-热-力三场耦合的物理实现层。
- `compute_heat_sources` 实现了 SPMe 电化学 -> 热模型的单向耦合（电化学热源 -> 温度场）。
- `ThermalDistributed2D_BC` 中的 CZM 界面热阻实现了力学 -> 热模型的反向耦合（损伤 -> 层间热阻）。
- `compute_heat_sources_with_czm` 实现了力学 -> 电化学耦合（断裂单元退出电化学反应，热源置零）。
- `per_element_spme` 参数使得每个热单元可以有独立的电化学状态，是 multi-SPMe 架构的关键。

### 哪些变更是耦合相关的
- 全部 8 个函数均为耦合相关：
  - `ThermalDistributed2D` + `ThermalDistributed2D_Ring` -- 热-电耦合（温度场求解）
  - `apply_convection_bc` + `apply_cool_method` -- 热边界条件
  - `ThermalDistributed2D_BC` -- 热-力耦合（CZM 界面热阻）
  - `compute_heat_sources` -- 电-热耦合（热源计算）
  - `compute_heat_sources_with_czm` -- 电-热-力三场耦合

### 哪些变更是独立的
无独立变更（全部为耦合服务）。
