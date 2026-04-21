---
name: Initialisation.jl 和 Solve.jl 重构计划
overview: 根据设计文档重构 Initialisation.jl 和 Solve.jl，通过新增 MultiSPMeState 结构体、合并初始化函数、整合 CallModel_MultiSPMe 到 CallModel，实现代码精简凝练。
todos:
  - id: "1"
    content: "Initialisation.jl: 新增 MultiSPMeState 结构体"
    status: pending
  - id: "2"
    content: "Initialisation.jl: 合并 ModelInitialisation 和 ModelInitialisation_MultiSPMe（内联逻辑，不用辅助函数）"
    status: pending
  - id: "3"
    content: "Initialisation.jl: 删除旧接口（extract/update/get_thermal_dofs），不保留兼容性"
    status: pending
  - id: "4"
    content: "Solve.jl: theta 合并到函数顶部一处赋值（删除行161-169重复）"
    status: pending
  - id: "5"
    content: "Solve.jl: 将 CallModel_MultiSPMe 内联到 CallModel"
    status: pending
  - id: "6"
    content: "Solve.jl: 使用 jellyroll_element_properties 获取面积，删除重复计算"
    status: pending
  - id: "7"
    content: "Solve.jl: 使用 Variables.jl 中的正确变量名（单元均温、节点温度）"
    status: pending
  - id: "8"
    content: "Solve.jl: 纯热模型分支提取为 solve_thermal_only() 独立函数"
    status: pending
  - id: "8.5"
    content: "Solve.jl: 调试代码提取为辅助函数（check_state_validity, check_temperature_field）"
    status: pending
  - id: "9"
    content: 更新所有调用点：CycleSolver.jl、PostProcessing.jl 等
    status: pending
  - id: "10"
    content: 验证：运行现有测试确保功能正常
    status: pending
isProject: false
---

## 重构范围

本计划基于设计文档：

- `docs/superpowers/specs/2026-03-24-multi-spme-state-refactor-design.md` - 状态结构与求解器重构

---

## 重要修改说明

### 变量名规范（来自 Variables.jl）


| 概念   | 变量名                                | 尺寸        |
| ---- | ---------------------------------- | --------- |
| 单元均温 | `"thermal2D temperature"`          | (ne, num) |
| 节点温度 | `"thermal2D temperature at nodes"` | (nT, num) |


### 面积计算（已存在于 Jellyrollmodel.jl）

```julia
# Jellyrollmodel.jl 第218-227行
function jellyroll_element_properties(mesh, param)
    ne = size(mesh.element, 1)
    areas = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    # ... 后续计算 layer_weights ...
    return areas, layer_weights
end
```

**重要**：不要再定义 `compute_element_areas`，直接调用 `jellyroll_element_properties(mesh, param)[1]` 获取面积。

---

## Initialisation.jl 修改（已完成）

### 1. 新增 MultiSPMeState 结构体（已完成）

### 2. 合并 ModelInitialisation 函数（已完成）

### 3. 删除旧接口（已完成）

---

## Solve.jl 修改（每项独立执行并测试）

### 4. theta 合并到 Solve 函数顶部

**当前问题**：theta 在两处独立赋值：

- thermal 分支（行46-54）
- electrochemical 分支（行161-169）

**修改方案**：将 theta 赋值移到 Solve 函数开头（行38 `if case.opt.model == "thermal"` 之前），thermal-only 的 early return 和 electrochemical 循环共享同一个 theta。

```julia
function Solve(case::Case; ...)
    # ... 日志初始化 ...
    result = nothing
    try
        # theta 时间离散系数（thermal 和 electrochemical 共用）
        if case.opt.solveType == "Crank-Nicolson"
            theta = 0.5
        elseif case.opt.solveType == "forward"
            theta = 0.0
        elseif case.opt.solveType == "backward"
            theta = 1.0
        else
            error("Error: $(case.opt.solveType) difference scheme has not been implemented!")
        end

        if case.opt.model == "thermal"
            # ... thermal-only 分支（不再需要重复赋值 theta）...
```

**删除代码**：行161-169 的 theta 赋值块。

**不提取 `get_theta()` 函数**：theta 赋值逻辑足够简单（3 个分支），合并到函数顶部一处赋值即可。

### 5. CallModel 内联 CallModel_MultiSPMe

**当前问题**：`CallModel_MultiSPMe`（行530-710）是独立函数，只有一个调用点在 `CallModel` 中（行744）。中间还有布局自动推断逻辑（行717-738）。

**修改方案**：将 `CallModel_MultiSPMe` 的核心逻辑直接搬到 `CallModel` 的 `multi_spme_enabled` 分支中。删除独立的 `CallModel_MultiSPMe` 函数。布局自动推断逻辑简化。

### 6. 面积计算统一用 jellyroll_element_properties

**当前重复**：

- `CallModel_MultiSPMe` 行576-583：手动循环 `A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]`
- `CallModel` distributed2D 分支行774-779：完全相同的循环

**修改方案**：两处统一替换为：

```julia
areas, layer_weights = jellyroll_element_properties(mesh_th, case.param)
variables["thermal2D element area"] = areas
```

### 7. 使用 Variables.jl 中的正确变量名

使用 `"thermal2D temperature"` 表示单元均温，`"thermal2D temperature at nodes"` 表示节点温度。

### 8. 纯热模型分支提取为 solve_thermal_only()

**当前问题**：行38-106 的 thermal-only 逻辑（约70行）内嵌在 Solve 函数中，增加函数复杂度。

**修改方案**：提取为独立函数：

```julia
function solve_thermal_only(case::Case; theta::Float64)
    # 行38-106 的全部逻辑
    # ...
    return (time = times, T_nodes = T_nodes, T_hist = T_hist)
end
```

Solve 中调用：

```julia
if case.opt.model == "thermal"
    return solve_thermal_only(case; theta=theta)
end
```

### 8.5 调试代码提取为辅助函数

**当前分散的调试代码**：

- Solve 行192-195：初始状态 NaN 检查
- Solve 行244-248：初始求解步骤异常检查
- CallModel_MultiSPMe 行559-567：温度场异常检查
- CallModel_MultiSPMe 行590-595：温度场 NaN 检查

**提取为**：

```julia
function check_state_validity(yt; context="")
    nan_count = sum(.!isfinite.(yt))
    if nan_count > 0
        prefix = isempty(context) ? "" : "[$context] "
        @warn "$(prefix)状态向量包含 $nan_count 个 NaN/Inf，长度 $(length(yt))"
    end
    return nan_count
end

function check_temperature_field(T_nodes; context="")
    nan_count = sum(.!isfinite.(T_nodes))
    abnormal = sum(abs.(T_nodes) .> 10.0)
    large_dev = sum(abs.(T_nodes .- 1.0) .> 5.0)
    if nan_count > 0 || abnormal > 0 || large_dev > 0
        T_min, T_max = extrema(T_nodes)
        prefix = isempty(context) ? "" : "[$context] "
        @warn "$(prefix)温度场异常" range=(T_min,T_max) nan=nan_count abnormal=abnormal large_dev=large_dev
    end
end
```

### 删除内容

- ~~`CallModel_MultiSPMe`~~ 独立函数 - 内联到 CallModel
- ~~theta 系数重复赋值~~（行161-169） - 合并到函数顶部
- ~~纯热模型分支内嵌逻辑~~（行38-106） - 提取为 `solve_thermal_only()`
- ~~重复的面积计算代码~~ - 使用 `jellyroll_element_properties`
- ~~分散的调试检查代码~~ - 提取为辅助函数

**预计行数减少**: 851 行 -> ~650 行（减少 200+ 行）

---

## 数据流示意

```mermaid
flowchart TB
    subgraph InitPhase["初始化阶段"]
        A[ModelInitialisation] --> B{multi_spme_mode?}
        B -->|是| C[MultiSPMeState<br/>Y_chem矩阵 + T_nodes向量]
        B -->|否| D[Vector{Float64}<br/>传统向量形式]
    end

    subgraph SolvePhase["求解阶段"]
        E[Solve] --> F[theta 赋值<br/>函数顶部一处]
        F --> G{thermal only?}
        G -->|是| H[solve_thermal_only]
        G -->|否| I{multi_spme_enabled?}
        I -->|是| J[CallModel<br/>多SPMe分支内联]
        I -->|否| K[CallModel<br/>原有分支]
        J --> L[jellyroll_element_properties<br/>获取面积]
        K --> L
    end
```



---

## 需要更新的调用点

### CycleSolver.jl


| 旧接口                                           | 新接口                         |
| --------------------------------------------- | --------------------------- |
| `ModelInitialisation_MultiSPMe(case)`         | `ModelInitialisation(case)` |
| `MultiSPMe_extract_element_state(y, e, case)` | `state.Y_chem[:, e]`        |
| `MultiSPMe_get_thermal_dofs(y, case)`         | `state.T_nodes`             |
| `MultiSPMe_update_state(y, case; ...)`        | 构造新 `MultiSPMeState`        |
| `case.multi_spme_layout["ne"]`                | `size(state.Y_chem, 2)`     |
| `case.multi_spme_layout["n_chem"]`            | `size(state.Y_chem, 1)`     |
| `case.multi_spme_layout["nT"]`                | `length(state.T_nodes)`     |


### PostProcessing.jl


| 旧接口                                       | 新接口                   |
| ----------------------------------------- | --------------------- |
| `case.multi_spme_layout["thermal_range"]` | 从 `MultiSPMeState` 获取 |
| 访问 `yt[thermal_range]`                    | `state.T_nodes`       |


### ThermalDistributed.jl

检查是否有调用 `compute_element_areas` 或 `compute_element_temperatures` 的代码，更新为内联或使用 `jellyroll_element_properties`。

---

## 保留内容

1. `multi_spme_layout` 中的 CZM 相关字段保留：
  - `czm_mesh`
  - `czm_element_map`
  - `interface_pairs`
  - `element_layer`
  - `is_inner_layer`
  - `inner_nodes`, `outer_nodes`
  - `thermal_variables`, `thermal_update_fn`, `thermal_record`
  - `polar_mesh_data`
2. 单 SPMe 模式（SPM, SPMe without distributed2D）保持原有向量形式
3. API 签名不变（`Solve`, `CallModel` 等）

---

## 验收标准

1. `Initialisation.jl` 行数减少 >= 100 行（从 240 行减少到 ~100 行）
2. `Solve.jl` 行数减少 >= 150 行（从 851 行减少到 ~650 行）
3. 初始化函数从 2 个减少到 1 个
4. `MultiSPMeState` 结构体正确定义
5. 删除 `multi_spme_layout` 中的 `ne`, `n_chem`, `nT`, `chem_range`, `thermal_range`
6. 不再定义 `compute_element_areas`（使用 `jellyroll_element_properties`）
7. 不新增 `get_theta()` 函数（theta 内联合并到 Solve 顶部）
8. 使用正确的变量名：`"thermal2D temperature"` 和 `"thermal2D temperature at nodes"`
9. 纯热模型逻辑在独立函数 `solve_thermal_only()` 中
10. 调试检查通过 `check_state_validity()` 和 `check_temperature_field()` 执行
11. 所有现有测试通过

---

## 实施顺序

### 阶段 1：Initialisation.jl 重构（已完成）

### 阶段 2：Solve.jl 重构（每项独立执行并测试）


| 步骤  | 任务                                   | 风险  | 测试                         |
| --- | ------------------------------------ | --- | -------------------------- |
| 4   | theta 合并到 Solve 顶部                   | 低   | 运行 SPMe_Thermal_example.jl |
| 5   | CallModel 内联 CallModel_MultiSPMe     | 中   | 运行全耦合 testexample.jl       |
| 6   | 面积计算统一用 jellyroll_element_properties | 低   | 对比热源计算结果                   |
| 7   | 使用 Variables.jl 中的正确变量名              | 低   | 检查输出变量                     |
| 8   | 纯热模型分支提取为 solve_thermal_only()       | 低   | 运行 thermal_verify.jl       |
| 8.5 | 调试代码提取为辅助函数                          | 低   | 运行任意仿真确认日志输出               |


### 阶段 3：更新调用点

1. 更新 `CycleSolver.jl`
2. 更新 `PostProcessing.jl`
3. 更新 `ThermalDistributed.jl`（如有需要）

### 阶段 4：验证

1. 运行所有测试
2. 对比重构前后数值结果

