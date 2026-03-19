# 代码风格统一和简化计划

## 目标
按照 main 分支代码编写风格（参见 `md/00_代码风格规范.md`），逐一检查 Parameters_Design 分支的新增代码，制定统一和简化计划。

## 基准规范
- **规范文档**: `md/00_代码风格规范.md`
- **基准分支**: main

---

## 阶段 1: 问题识别 [完成]

### 1.1 以 `_` 开头的内部函数 (需修改) ⚠️ 优先级高

根据 main 分支规范，函数应使用 **snake_case** 或 **功能前缀+功能描述**，不应使用 `_` 前缀。

#### PostProcessing.jl (9 个函数)

| 行号 | 原函数名 | 建议新名称 | 说明 |
|------|----------|-----------|------|
| 56 | `_phase_termination_symbol` | `GetPhaseTerminationSymbol` | 获取阶段截止符号 |
| 66 | `_state_concentration_variance` | `ComputeConcentrationVariance` | 计算浓度方差 |
| 97 | `_postprocess_phase_result` | `PostprocessPhaseResult` | 后处理阶段结果 |
| 157 | `_postprocess_cycle_result!` | `PostprocessCycleResult!` | 后处理循环结果 |
| 174 | `_append_cycle_result!` | `AppendCycleResult!` | 追加循环结果 |
| 192 | `_update_soh_and_capacity!` | `UpdateSohAndCapacity!` | 更新SOH和容量 |
| 203 | `_print_cycle_summary` | `PrintCycleSummary` | 打印循环摘要 |
| 209 | `_check_cycle_termination` | `CheckCycleTermination` | 检查循环终止条件 |
| 230 | `_print_cycling_summary` | `PrintCyclingSummary` | 打印循环总结 |

#### Parallelsolution.jl (11 个函数)

| 行号 | 原函数名 | 建议新名称 | 说明 |
|------|----------|-----------|------|
| 7 | `_debug_check_prefactors` | `DebugCheckPrefactors` | 调试：检查预因子 |
| 33 | `_debug_check_coefficients` | `DebugCheckCoefficients` | 调试：检查系数 |
| 55 | `_debug_check_initial_voltage` | `DebugCheckInitialVoltage` | 调试：检查初始电压 |
| 83 | `_compute_electrochemical_prefactors` | `ComputeElectrochemicalPrefactors` | 计算电化学预因子 |
| 119 | `_compute_element_coefficients` | `ComputeElementCoefficients` | 计算单元系数 |
| 152 | `_compute_all_coefficients` | `ComputeAllCoefficients` | 计算所有系数 |
| 162 | `_branch_voltage` | `ComputeBranchVoltage` | 计算支路电压 |
| 169 | `_branch_dVdI` | `ComputeBranchDVdI` | 计算支路电压对电流导数 |
| 179 | `_initialize_currents` | `InitializeBranchCurrents` | 初始化支路电流 |
| 194 | `_check_voltage_bounds` | `CheckVoltageBounds` | 检查电压边界 |
| 226 | `_detect_cutoff_elements` | `DetectCutoffElements` | 检测截止单元 |
| 294 | `_newton_iteration!` | `NewtonIteration!` | Newton迭代 |
| 391 | `_line_search` | `LineSearchBranchCurrents` | 线搜索 |

#### CycleSolver.jl (4 个函数)

| 行号 | 原函数名 | 建议新名称 | 说明 |
|------|----------|-----------|------|
| 551 | `_compute_czm_effective_params` | `ComputeCzmEffectiveParams` | 计算CZM有效参数 |
| 581 | `_compute_czm_strain_inputs` | `ComputeCzmStrainInputs` | 计算CZM应变输入 |
| 658 | `_update_czm_damage!` | `UpdateCzmDamage!` | 更新CZM损伤 |
| 708 | `_ensure_multi_spme_layout!` | `EnsureMultiSpmeLayout!` | 确保多SPMe布局 |

### 1.2 结构体风格问题

| 文件 | 结构体 | 问题 | 建议 |
|------|--------|------|------|
| CycleSolver.jl | PhaseResult | 无 `@with_kw` | 添加 `@with_kw` |
| CycleSolver.jl | CycleResult | 无 `@with_kw` | 添加 `@with_kw` |
| CycleSolver.jl | CyclingResult | 无 `@with_kw` | 添加 `@with_kw` |
| czm.jl | CohesiveElement | 无 `@with_kw` | 添加 `@with_kw` |
| czm.jl | DamageState | 无 `@with_kw` | 添加 `@with_kw` |

### 1.3 缺少 docstring 的核心函数

| 文件 | 函数 | 状态 |
|------|------|------|
| ThermalDistributed.jl | `ThermalDistributed2D` | 缺少 |
| ThermalDistributed.jl | `apply_convection_bc` | 缺少 |
| Jellyrollmodel.jl | `jellyroll_collector_seed_mesh` | 缺少 |
| Parallelsolution.jl | `solve_branch_currents_newton` | 需补充 |
| PostProcessing.jl | `PostProcessing` | 缺少 |

---

## 阶段 2: 详细代码审查 [完成]

### 2.1 新增核心模块检查结果

| 文件 | 行数 | `_`函数 | 结构体问题 | docstring | 总问题 |
|------|------|---------|-----------|-----------|--------|
| CycleSolver.jl | +729 | 4 | 3 | 部分缺失 | 7+ |
| CycleData.jl | +621 | 0 | 0 | OK | 0 |
| Parallelsolution.jl | +619 | 11 | 0 | 部分缺失 | 12+ |
| czm.jl | +503 | 0 | 2 | 部分缺失 | 2+ |
| CzmSolve.jl | +514 | 0 | 0 | 部分缺失 | 0+ |
| Jellyrollmodel.jl | +547 | 0 | 0 | 缺失 | 1+ |
| Materialmatrix.jl | +382 | 0 | 0 | 需检查 | ? |
| ThermalDistributed.jl | +390 | 0 | 0 | 缺失 | 1+ |
| ThermalPolar2D.jl | +142 | 0 | 0 | 需检查 | ? |
| Tools.jl | +178 | 0 | 0 | 需检查 | ? |
| PostProcessing.jl | +241 | 9 | 0 | 部分缺失 | 9+ |

**总计**: 25 个 `_` 前缀函数需重命名，5 个结构体需添加 `@with_kw`

---

## 阶段 3: 修改执行计划

### 3.1 函数重命名计划

#### 第一批：Parallelsolution.jl（优先级最高）

此文件被 Solve.jl 核心求解器调用，影响最大。

```julia
# 需要修改的调用链
Parallelsolution.jl (定义) → Solve.jl (调用) → 其他模块
```

**重命名列表**:
1. `_compute_electrochemical_prefactors` → `ComputeElectrochemicalPrefactors`
2. `_compute_element_coefficients` → `ComputeElementCoefficients`
3. `_compute_all_coefficients` → `ComputeAllCoefficients`
4. `_branch_voltage` → `ComputeBranchVoltage`
5. `_branch_dVdI` → `ComputeBranchDVdI`
6. `_initialize_currents` → `InitializeBranchCurrents`
7. `_check_voltage_bounds` → `CheckVoltageBounds`
8. `_detect_cutoff_elements` → `DetectCutoffElements`
9. `_newton_iteration!` → `NewtonIteration!`
10. `_line_search` → `LineSearchBranchCurrents`
11. `_debug_check_prefactors` → `DebugCheckPrefactors`
12. `_debug_check_coefficients` → `DebugCheckCoefficients`
13. `_debug_check_initial_voltage` → `DebugCheckInitialVoltage`

#### 第二批：CycleSolver.jl

**重命名列表**:
1. `_compute_czm_effective_params` → `ComputeCzmEffectiveParams`
2. `_compute_czm_strain_inputs` → `ComputeCzmStrainInputs`
3. `_update_czm_damage!` → `UpdateCzmDamage!`
4. `_ensure_multi_spme_layout!` → `EnsureMultiSpmeLayout!`

#### 第三批：PostProcessing.jl

**重命名列表**:
1. `_phase_termination_symbol` → `GetPhaseTerminationSymbol`
2. `_state_concentration_variance` → `ComputeConcentrationVariance`
3. `_postprocess_phase_result` → `PostprocessPhaseResult`
4. `_postprocess_cycle_result!` → `PostprocessCycleResult!`
5. `_append_cycle_result!` → `AppendCycleResult!`
6. `_update_soh_and_capacity!` → `UpdateSohAndCapacity!`
7. `_print_cycle_summary` → `PrintCycleSummary`
8. `_check_cycle_termination` → `CheckCycleTermination`
9. `_print_cycling_summary` → `PrintCyclingSummary`

### 3.2 结构体修改计划

#### CycleSolver.jl

```julia
# 修改前
mutable struct PhaseResult
    phase_type::PhaseType
    # ...
end

# 修改后
@with_kw mutable struct PhaseResult
    phase_type::PhaseType = PHASE_REST
    t_start::Float64 = 0.0
    # ...
end
```

同样应用于 `CycleResult` 和 `CyclingResult`。

#### czm.jl

```julia
# 修改前
mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    # ...
end

# 修改后
@with_kw mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64 = 0
    nodes::Vector{Int64} = Int64[]
    # ...
end
```

### 3.3 docstring 补充计划

#### ThermalDistributed.jl

```julia
"""
    ThermalDistributed2D(case, variables)

计算二维分布式热模型的质量矩阵、刚度矩阵和载荷向量。

# 输入
- `case::Case`: 案例对象，包含参数和网格
- `variables::Dict`: 当前变量状态

# 输出
- `MT`: 质量矩阵 (稀疏矩阵)
- `KT`: 刚度矩阵 (稀疏矩阵)
- `FT`: 载荷向量
"""
function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
```

#### Jellyrollmodel.jl

```julia
"""
    jellyroll_collector_seed_mesh(param; nθ=360, gsorder=2, phase=0.0, tol=1e-8)

生成 Jellyroll 电池的 collector-seeded 网格。

基于阿基米德螺旋线 r(θ) = a + bθ 生成网格。

# 参数
- `param`: 参数对象，包含几何信息
- `nθ::Int`: 周向网格分辨率 (默认 360)
- `gsorder::Int`: 高斯积分阶数 (默认 2)
- `phase::Float64`: 相位偏移 (默认 0.0)
- `tol::Float64`: 节点重合判断容差 (默认 1e-8)

# 返回
- `Mesh`: 包含节点、单元和高斯积分数据的网格对象
"""
function jellyroll_collector_seed_mesh(param; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0, tol::Float64=1e-8)
```

---

## 阶段 4: 执行顺序

### Step 1: 函数重命名（按文件顺序）

1. **Parallelsolution.jl** - 11 个函数
2. **CycleSolver.jl** - 4 个函数
3. **PostProcessing.jl** - 9 个函数

### Step 2: 结构体添加 @with_kw

1. **CycleSolver.jl** - PhaseResult, CycleResult, CyclingResult
2. **czm.jl** - CohesiveElement, DamageState

### Step 3: 添加 docstring

1. **ThermalDistributed.jl** - ThermalDistributed2D, apply_convection_bc
2. **Jellyrollmodel.jl** - jellyroll_collector_seed_mesh
3. **Parallelsolution.jl** - solve_branch_currents_newton

### Step 4: 验证

- [ ] 运行 `example/testexample.jl`
- [ ] 运行 `example/热模块验证/thermal_verify.jl`
- [ ] 运行 `example/czm/czm_cycle_example.jl`

---

## 状态

- [x] 阶段 1: 问题识别
- [x] 阶段 2: 详细代码审查
- [x] 阶段 3: 修改执行计划
- [ ] 阶段 4: 执行修改
- [ ] 阶段 5: 验证与测试

---

## 风险评估

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| 重命名后调用点遗漏 | 编译错误 | 全局搜索旧函数名确认无遗漏 |
| @with_kw 默认值不当 | 运行时错误 | 保持现有构造函数逻辑 |
| docstring 不准确 | 文档误导 | 对照代码仔细核对 |

---

*创建日期: 2026-03-19*
*最后更新: 2026-03-19*
