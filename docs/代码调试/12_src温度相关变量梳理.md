# src/ 温度相关变量梳理

> 生成日期: 2026-04-07
> 范围: `src/` 下所有 `.jl` 文件中的温度相关命名

---

## 1. 参数层 (物理量/无量纲参数)

### 1.1 `scale` 结构 — 归一化参考尺度

| 变量 | 类型 | 定义位置 | 值/含义 |
|------|------|----------|---------|
| `scale.T_ref` | Float64 | `SetParams.jl:179` | 温度参考值，默认 298 K |
| `scale.phi` | Float64 | `SetParams.jl:282` | 电势尺度 = T_ref × R / F |
| `scale.lambda` | Float64 | `SetParams.jl:295` | 导热率尺度 = P_ref / (L × T_ref) |

### 1.2 `cell` 参数 — 电池级温度参数

| 变量 | 类型 | 定义位置 | 归一化 | 含义 |
|------|------|----------|--------|------|
| `cell.T0` | Float64 | `SetParams.jl:399` | T0_dim / T_ref | 初始温度 (无量纲) |
| `cell.T_amb` | Float64 | `SetParams.jl:119,398` | T_amb_dim / T_ref | 环境温度 (无量纲) |

**参数文件中赋值（物理量 K）：**

| 参数集 | cell.T0 | cell.T_amb | 位置 |
|--------|---------|------------|------|
| Jellyroll | 298.15 | = T0 | `parameters/Jellyroll.jl:137-138` |
| LGM50 | 298 | 298 | `parameters/LGM50.jl:101-102` |
| Northrop | 298 | = T0 | `parameters/Northrop.jl:104-105` |
| Ring | 298.0 | 298.0 | `parameters/Ring.jl:55-56` |
| Enertech | 298.15 | — | `parameters/Enertech.jl:143` |

---

## 2. 索引层 (case.index)

| Key | 含义 | 设置位置 |
|-----|------|----------|
| `"temperature"` | 单温度 DOF 在状态向量中的索引 | `SetCase.jl:77,81` |

用于 SPM/SPMe/P2D 的单温度自由度（lumped 模式下）。

---

## 3. variables Dict — 核心温度键

### 3.1 标量温度（单值/全局）

| Dict Key | 维度 | 单位域 | 写入位置 | 读取位置 | 含义 |
|----------|------|--------|----------|----------|------|
| `"temperature"` | (1, num) | 无量纲 | SPM.jl:83, SPMe.jl:195, P2D.jl:283, CallModel.jl:148 | Thermal.jl:4, SPMe.jl:137,143, SPM.jl:59-60, P2D.jl:108,142,216, ElectrolyteDiffusion.jl:23, ElectrolytePotential.jl:17, Mechanical.jl:11,39 | 全局温度标量（单模型为 DOF 值；MultiSPMe 为体积平均温度） |

### 3.2 二维热场（节点温度）

| Dict Key | 维度 | 单位域 | 写入位置 | 读取位置 | 含义 |
|----------|------|--------|----------|----------|------|
| `"thermal2D temperature at nodes"` | (nT, num) | 无量纲 (T/T_ref) | CallModel.jl:47,62,149, Solve.jl:184,237 | Solve.jl:412-413,421, PostProcessing.jl:66 | 当前时刻节点温度场 |
| `"thermal2D temperature  at nodes history"` | (nT, num) | 无量纲 | Variables.jl:100 | — | 节点温度历史（注意双空格拼写） |
| `"thermal2D temperature"` | (ne, num) | 无量纲 | Variables.jl:97 | PostProcessing.jl:78 | 单元均温历史 |
| `"thermal2D temperature history"` | (ne, num) | 无量纲 | Variables.jl:99 | — | 单元温度历史 |
| `"T_nodes"` | — | 无量纲 | Solve.jl:25,60, CallModel间接 | Solve.jl:25,187, Mechanical.jl:174, CycleData.jl:54 | 当前节点温度的即时引用（非历史记录） |

### 3.3 热学边界/统计量

| Dict Key | 维度 | 单位域 | 含义 |
|----------|------|--------|------|
| `"thermal lumped internal heat"` | (1, num) | 无量纲 | 集总热模型的内部热源 |

---

## 4. 结果输出 Dict (result)

### 4.1 物理单位输出

| Key | 维度 | 单位 | 写入位置 | 含义 |
|-----|------|------|----------|------|
| `"temperature [K]"` | (v,) | K | PostProcessing.jl:6 | 标量温度时间序列 |
| `"thermal2D temperature at nodes [K]"` | (nT, v) | K | PostProcessing.jl:63, Solve.jl:413 | 节点温度场时间序列 |
| `"thermal2D temperature [K]"` | (ne, v) | K | Solve.jl:424 | 单元均温时间序列 |
| `"thermal2D final temperature at nodes [K]"` | (nT,) | K | Solve.jl:429 | 最终节点温度 |
| `"T_max"` | Float64 | K | CyclePostProcess.jl:54-61, CycleData.jl:174-190 | 最高温度 |
| `"T_mean"` | Float64 | K | CycleData.jl:174-199 | 平均温度 |

### 4.2 循环状态传递

| Key | 类型 | 单位域 | 写入位置 | 含义 |
|-----|------|--------|----------|------|
| `"T_nodes"` | Vector{Float64} 或 nothing | 无量纲 | Solve.jl:443, CycleData.jl:252, CycleSolver.jl:243,274,349 | 跨周期/跨相位温度场传递 |

---

## 5. 局部变量（函数内部）

### 5.1 节点温度类

| 变量名 | 类型 | 所在函数 | 含义 |
|--------|------|----------|------|
| `T_nodes` | Vector{Float64} | CallModel_MultiSPMe, Solve, CycleData | 节点温度向量（无量纲 T/T_ref） |
| `T_nodes_carry` | Vector{Float64} | Solve, CycleData | 跨时间步传递的当前温度场 |
| `T_nodes_input` | Vector{Float64} | Solve | 输入的初始温度场 |
| `T_nodes_nd` | Vector{Float64} | CycleSolver | 无量纲节点温度 |
| `T_nodes_out` | Vector{Float64} | CycleData | 导出用节点温度 |
| `T_nodes_K` | Vector{Float64} | CycleData, PostProcessing | 转换为 Kelvin 的节点温度 |
| `T_nodes_hist` | Matrix{Float64} | Solve | 节点温度历史 |
| `T_seed` | Vector{Float64} | Solve | 初始化种子温度 |
| `T_hist` | Matrix{Float64} | Solve (polar) | 极坐标模式温度历史 |

### 5.2 单元温度类

| 变量名 | 类型 | 所在函数 | 含义 |
|--------|------|----------|------|
| `Te_prev` | Vector{Float64} | CallModel_MultiSPMe, Solve, Parallelsolution | 各单元均温（上一步/当前步） |
| `T_e` | Float64 (kwarg) | SPMe_element, SPMe_variables, compute_heat_sources | 单个单元温度 |
| `T_rep` | Float64 | CallModel_MultiSPMe | 所有单元平均温度（代表性温度） |
| `T_elem` | Vector{Float64} | Mechanical, CycleSolver | 各单元温度 |
| `dT_elem` | Vector{Float64} | Mechanical, CycleSolver, czm, CzmSolve | 各单元温度变化 ΔT = T_elem - T0 |

### 5.3 标量温度类

| 变量名 | 类型 | 所在函数 | 含义 |
|--------|------|----------|------|
| `T` | Float64 或 Vector{Float64} | SPM, SPMe, P2D, Thermal, Electrolyte* | 当前温度（来自状态向量或变量） |
| `T_ref` | Float64 | CycleSolver, CycleData | 等于 cell.T0（物理量 K）或 scale.T_ref |
| `T_ref_scale` | Float64 | CycleSolver | = param_dim.scale.T_ref |
| `T_amb` | Float64 | ThermalDistributed, ThermalPolar2D | 环境温度（无量纲） |
| `T_amb_nd` | Float64 | ThermalPolar2D | 无量纲环境温度 |
| `T0` | Float64 | Mechanical, Initialisation, 参数文件 | 初始温度（无量纲或物理量视上下文） |
| `Tm` | Float64 | Solve | 体积平均温度（无量纲） |
| `T_max_phase` | Float64 | CycleData | 当前相位最高温度 (K) |
| `T_max_current` | Float64 | CycleData | 当前步最高温度 (K) |
| `T_max_K` | Float64 | CycleData | 转换为 Kelvin 的最高温度 |
| `T_mean_K` | Float64 | CycleData | 转换为 Kelvin 的平均温度 |

### 5.4 CZM 相关温度

| 变量名 | 类型 | 所在函数 | 含义 |
|--------|------|----------|------|
| `T_elem_nd` | Float64 (局部) | CycleSolver:595 | 单元无量纲温度累加器 |
| `T_elem_dim` | Float64 (局部) | CycleSolver:602 | 单元有量纲温度 = T_elem_nd × T_ref_scale |

---

## 6. 结构体字段

| 结构体 | 字段 | 类型 | 含义 |
|--------|------|------|------|
| `TimeStepData` | `T_nodes` | Vector{Float64} | 节点温度 (K) |
| `TimeStepData` | `T_max` | Float64 | 最高温度 (K) |
| `TimeStepData` | `T_mean` | Float64 | 平均温度 (K) |

---

## 7. 函数参数/关键字参数

| 函数签名 | 参数 | 类型 | 含义 |
|----------|------|------|------|
| `SPMe_element(...; I_e, T_e, jacobi)` | `T_e` | Float64 | 单元温度 |
| `SPMe_variables(...; I_app, T_e)` | `T_e` | Float64 或 Nothing | 单元温度（Nothing 时从 index 取） |
| `compute_heat_sources(..., I_e, T_e, areas)` | `T_e` | Vector{Float64} | 各单元均温向量 |
| `compute_heat_sources_with_czm(..., I_e, T_e, areas, ...)` | `T_e` | Vector{Float64} | 各单元均温向量 |
| `solve_branch_currents(..., Te_prev, ...)` | `Te_prev` | Vector{Float64} | 各单元温度 |
| `compute_element_coefficients(e, T_e, param, ...)` | `T_e` | Float64 | 单元温度 |
| `compute_all_coefficients(ne, Te_prev, ...)` | `Te_prev` | Vector{Float64} | 各单元温度 |
| `thermal2D_volume_average_temperature(mesh, T_nodes)` | `T_nodes` | AbstractVector | 节点温度 |
| `element_nodal_mean(mesh, T_nodes)` | — (nodal_values) | AbstractVector | 节点值（温度） |

---

## 8. 命名规律总结

### 8.1 量纲后缀约定

| 后缀 | 含义 | 示例 |
|------|------|------|
| (无后缀) | 无量纲 (T/T_ref) | `T_nodes`, `Te_prev`, `T_e` |
| `_K` | 物理量 (Kelvin) | `T_nodes_K`, `T_max_K`, `T_mean_K` |
| `_nd` | 无量纲 (显式标注) | `T_nodes_nd`, `T_amb_nd`, `T_elem_nd` |
| `_dim` | 有量纲/物理量 | `T_elem_dim` |
| `_ref` / `_ref_scale` | 参考量 | `T_ref`, `T_ref_scale` |
| `[K]` (Dict key) | 结果输出中的物理量 | `"temperature [K]"`, `"thermal2D temperature at nodes [K]"` |

### 8.2 语义前缀

| 前缀/模式 | 含义 | 示例 |
|-----------|------|------|
| `T_nodes` | 节点温度向量 | `T_nodes_carry`, `T_nodes_input` |
| `Te_prev` / `T_e` | 单元均温 | `Te_prev`, `T_e` |
| `T_elem` | 单元温度 | `T_elem`, `T_elem_K` |
| `dT_elem` | 温度变化 | `dT_elem` |
| `T_max` | 最高温度 | `T_max_phase`, `T_max_K` |
| `T_mean` | 平均温度 | `T_mean_K`, `T_mean_end` |

### 8.3 易混淆命名

| 变量 | 注意事项 |
|------|----------|
| `T_ref` | 在 CycleSolver 中 = cell.T0 (K)；在 CycleData 中 = scale.T_ref (K)。两者值相同（298 K），但来源不同 |
| `T0` | 在参数文件中是物理量 (K)，在 Mechanical 中是无量纲 param.cell.T0 |
| `cell.T0` | 归一化后是无量纲 T/T_ref；初始赋值时是物理量 K |
| `"temperature"` vs `"T_nodes"` | 前者是 variables 历史键（维度 1×num）；后者是即时状态引用 |
| `"thermal2D temperature  at nodes history"` | 注意双空格拼写（Variables.jl:100） |
| `T` (局部) | 不同函数中含义不同：SPM/SPMe 中是标量温度；Thermal/SPMe 中可能是向量 |
