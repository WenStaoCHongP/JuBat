# 标识符改名重构设计

**日期**: 2026-07-07
**类型**: 重构（行为不变）
**作者**: brainstorm session
**状态**: 待审查

---

## 1. 背景与动机

`src/` 中存在多处**语义缩写**——变量名被压缩到失去原意，造成：

1. **可读性差**：`fks`、`Te_prev`、`T_nodes_carry` 让新读者难以理解含义
2. **跨文件不一致**：同一概念在不同文件用不同缩写（`mesh_th` vs `mesh`）
3. **搜索困难**：单字母或两字母变量名难以 grep
4. **代码评审负担**：每行都要回溯变量来源

本次重构目标：**对 11 类语义缩写进行改名**，统一命名规范，**完全不改变运行行为**。

## 2. 范围

### 2.1 在范围内（15 项）

**A. 标识符改名（10 项）**

| # | 当前名 | 目标名 | 性质 |
|---|--------|--------|------|
| 1 | `p` (= `case.param`) | `param` | 局部别名 |
| 2 | `mesh_th` | `mesh_thermal` | 跨 5 文件 |
| 3 | `ws` (自由变量) | `workspace` | 局部 |
| 4 | `ws_e` | `elem_workspace` | 局部 |
| 5 | `ws_pool` | `workspace_pool` | 局部 |
| 6 | `vars_e` | `elem_vars` | 局部 |
| 7 | `vars_hist` | `history_vars` | 跨文件参数 |
| 8 | `Te_prev` | `T_elem_prev` | 局部 |
| 9 | `fks` | `layer_weights` | 含义不明→明确 |
| 10 | `T_nodes_carry` | `T_nodes_step` | 跨文件 |

**B. 纯属性别名内联（5 项）**

| # | 位置 | 当前模式 | 处理 |
|---|------|----------|------|
| 11 | `CouplingState.jl:230-231` | `t_pe, t_ne, t_sp, t_pcc, t_ncc = param.X.thickness ...` | 内联到公式（保留 `Σt`，因 3× 分母复用） |
| 12 | `ThermalDistributed.jl:59` | `T_amb = param.cell.T_amb` | 后续 1 次引用直接用 `param.cell.T_amb` |
| 13 | `ThermalDistributed.jl:185` | `T_amb = param.cell.T_amb` | 同上 |
| 14 | `CsvExport.jl:478` | `T_ref = scale.T_ref` | 多处引用直接用 `scale.T_ref` |
| 15 | `CycleData.jl:41` | `T_ref = case.param_dim.scale.T_ref` | 多处引用直接用完整路径 |

**说明**：第 11-15 项与第 1 项（`p` → `param`）在 `compute_effective_coating_modulus` 内部协同修改。

### 2.2 不在范围内

以下命名**保持不变**：

- **物理量单字母**：`T`（温度）、`I`（电流）、`V`（电压）、`t`（时间）
- **矩阵命名**：`MT`、`KT`、`FT`、`M_e`、`K_e`、`F_e`、`M_chem`、`K_chem`
- **紧凑循环计数器**：`i`、`j`、`k`、`e`（单元）、`n`（节点）、`g`（高斯点）、`v`（紧凑步数累计，行业惯例）
- **希腊字母变体**：`α`、`β`、`ν`、`σ`、`δ`（保留 Unicode）
- **计算表达式命名**：如 `sig_n_eff = param.NE.sig * param.NE.eps_s`、`j_n = I_app / param.NE.as / param.NE.thickness`（有实际计算，非纯别名）
- **Dict 值提取**：如 `eta_n = vars["..."]`（属另一 spec 范畴）

### 2.3 关联但不在本次范围

**变量字典注册表重构**（`StandardVariables` / `create_element_workspace` / `copy_element_results` 三处重复消除）属于另一独立 spec，本次不涉及。

## 3. 详细改名映射

### 3.1 第 1 项 + 第 11 项：`p` → `param` 且 `t_xx` 内联

**位置**：`src/CouplingState.jl:228-240`（`compute_effective_coating_modulus` 函数）

**Before**:
```julia
function compute_effective_coating_modulus(case)
    p = case.param
    t_pe, t_ne = p.PE.thickness, p.NE.thickness
    t_sp, t_pcc, t_ncc = p.SP.thickness, p.PCC.thickness, p.NCC.thickness
    Σt = t_pe * 2 + t_ne * 2 + t_sp * 2 + t_pcc + t_ncc
    E_eff = (p.PE.E_coat*t_pe *2 + p.NE.E_coat*t_ne *2 +
             p.SP.E*t_sp *2 + p.PCC.E*t_pcc + p.NCC.E*t_ncc) / Σt
    ν_eff = (p.PE.nu_coat*t_pe*2 + p.NE.nu_coat*t_ne*2 +
             p.SP.nu*t_sp*2 + p.PCC.nu*t_pcc + p.NCC.nu*t_ncc) / Σt
    α_eff = (p.PE.alphaT*t_pe*2 + p.NE.alphaT*t_ne*2 +
             p.SP.alphaT*t_sp*2 + p.PCC.alphaT*t_pcc + p.NCC.alphaT*t_ncc) / Σt
    return E_eff, ν_eff, α_eff
end
```

**After**:
```julia
function compute_effective_coating_modulus(case)
    param = case.param
    # Σt 保留：3 个分母复用，是真 DRY
    Σt = param.PE.thickness*2 + param.NE.thickness*2 +
         param.SP.thickness*2 + param.PCC.thickness + param.NCC.thickness
    E_eff = (param.PE.E_coat*param.PE.thickness*2 +
             param.NE.E_coat*param.NE.thickness*2 +
             param.SP.E*param.SP.thickness*2 +
             param.PCC.E*param.PCC.thickness +
             param.NCC.E*param.NCC.thickness) / Σt
    ν_eff = (param.PE.nu_coat*param.PE.thickness*2 +
             param.NE.nu_coat*param.NE.thickness*2 +
             param.SP.nu*param.SP.thickness*2 +
             param.PCC.nu*param.PCC.thickness +
             param.NCC.nu*param.NCC.thickness) / Σt
    α_eff = (param.PE.alphaT*param.PE.thickness*2 +
             param.NE.alphaT*param.NE.thickness*2 +
             param.SP.alphaT*param.SP.thickness*2 +
             param.PCC.alphaT*param.PCC.thickness +
             param.NCC.alphaT*param.NCC.thickness) / Σt
    return E_eff, ν_eff, α_eff
end
```

**改动说明**：
- `p` → `param`（第 1 项）
- `t_pe/t_ne/t_sp/t_pcc/t_ncc` 全部内联为 `param.X.thickness`（第 11 项）
- `Σt` 保留，因为它是 3 个分母的真复用，不是简单别名

### 3.2 第 2 项：`mesh_th` → `mesh_thermal`

**位置**（共 6 处赋值，引用更多）：
- `src/Jellyrollmodel.jl:537, 540`
- `src/CallModel.jl:31`
- `src/CsvExport.jl:256, 329`
- `src/CycleData.jl:112, 269`
- `src/PostProcessing.jl:69`

**Before**: `mesh_th = case.mesh["thermal2D"]`
**After**: `mesh_thermal = case.mesh["thermal2D"]`

函数体内所有 `mesh_th.xxx` 引用一并改。

### 3.3 第 3-5 项：`ws` 系列改名

| 当前 | 目标 | 位置 |
|------|------|------|
| `ws`（自由变量） | `workspace` | Variables.jl:164；czm.jl:164；CzmSolve.jl:275, 503 |
| `ws_e` | `elem_workspace` | CallModel.jl:132 |
| `ws_pool` | `workspace_pool` | CallModel.jl:126 |

注：`workspace` 仍是局部变量名，**不**改成 `case.workspace` 之类的属性，避免引入新的 struct 字段。

### 3.4 第 6-7 项：`vars_*` 改名

| 当前 | 目标 | 位置 |
|------|------|------|
| `vars_e` | `elem_vars` | CallModel.jl:137, 167；ThermalDistributed.jl:413 |
| `vars_hist` | `history_vars` | Variables.jl:236（函数形参）；Solve.jl 调用处 |

**函数签名同时改**：
```julia
# Before
function Variable_update!(vars_hist::Dict{...}, vars::Dict{...}, v::Int64)

# After
function Variable_update!(history_vars::Dict{...}, current_vars::Dict{...}, v::Int64)
```

注：第三个参数 `v`（步索引）保留原命名（紧凑循环计数惯例）。

### 3.5 第 8 项：`Te_prev` → `T_elem_prev`

**位置**：`src/CallModel.jl:51`（及函数体内多处引用）

**Before**:
```julia
Te_prev = zeros(Float64, ne)
@inbounds for e in 1:ne
    nds = mesh_th.element[e, :]
    Te_prev[e] = sum(T_nodes[nds]) / length(nds)
end
```

**After**:
```julia
T_elem_prev = zeros(Float64, ne)
@inbounds for e in 1:ne
    nds = mesh_thermal.element[e, :]
    T_elem_prev[e] = sum(T_nodes[nds]) / length(nds)
end
```

**含义**：上一步的逐单元温度（elemental temperature）。

### 3.6 第 9 项：`fks` → `layer_weights`

**位置**：`src/ThermalDistributed.jl:15, 385`

**Before**:
```julia
fks = case.geometry !== nothing ? case.geometry.layer_weights : jellyroll_element_properties(mesh, case.param)[2]
```

**After**:
```julia
layer_weights = case.geometry !== nothing ? case.geometry.layer_weights : jellyroll_element_properties(mesh, case.param)[2]
```

**含义澄清**：`fks` 原意不明，实际是"5 元层面积权重矩阵"（`ne × 5`，对应 NE/SP/PE/PCC/NCC），改名后含义自明。

### 3.7 第 10 项：`T_nodes_carry` → `T_nodes_step`

**位置**：`src/CycleData.jl:71, 78, 81, 149, 152`；`src/Solve.jl:181, 228`

**Before**: `T_nodes_carry = ...` （"carry" 含义模糊——是携带？保留？传递？）
**After**: `T_nodes_step = ...` （明确表示"当前时间步的节点温度"）

### 3.8 第 11-15 项：纯属性别名内联

#### 3.8.1 第 11 项

已在 §3.1 与 `p` → `param` 合并展示。

#### 3.8.2 第 12-13 项：`T_amb = param.cell.T_amb` 内联

**位置**：`src/ThermalDistributed.jl:59` 和 `:185`

**Before**:
```julia
function apply_convection_bc(KT, FT, mesh, is_outer, case; edge_cache=nothing)
    ...
    param = case.param
    T_amb = param.cell.T_amb  # 已无量纲
    ...
    fe1 += wt * T_amb * N1
    fe2 += wt * T_amb * N2
end
```

**After**:
```julia
function apply_convection_bc(KT, FT, mesh, is_outer, case; edge_cache=nothing)
    ...
    param = case.param
    ...
    fe1 += wt * param.cell.T_amb * N1
    fe2 += wt * param.cell.T_amb * N2
end
```

**影响**：每个 `T_amb` 引用替换为 `param.cell.T_amb`，共约 4 处（2 函数 × 2 引用）。

#### 3.8.3 第 14 项：`T_ref = scale.T_ref` 内联

**位置**：`src/CsvExport.jl:478`

**Before**: `T_ref = scale.T_ref` 后续多处 `* T_ref`
**After**: 删除局部别名，所有引用替换为 `scale.T_ref`

#### 3.8.4 第 15 项：`T_ref = case.param_dim.scale.T_ref` 内联

**位置**：`src/CycleData.jl:41`

**Before**: `T_ref = case.param_dim.scale.T_ref`
**After**: 删除局部别名，引用替换为 `case.param_dim.scale.T_ref`

**注意**：若该函数内引用次数 ≥ 5，可考虑保留 `T_ref` 别名（防止行过长）。实施时按实际引用次数判断。

## 4. 实施约束

### 4.1 必须满足

- ✅ **运行行为完全不变**（纯改名，无逻辑改动）
- ✅ **测试套件全部通过**（若有现成测试）
- ✅ **逐文件改动、每次只改一项**，便于 review 与回滚
- ✅ **每处修改前先向用户展示 before/after**

### 4.2 不允许

- ❌ 借机"顺手"做其他重构（如键注册表、struct 化）
- ❌ 修改函数签名（仅改局部变量名/参数名）
- ❌ 引入新依赖或新文件
- ❌ 改动字符串字面量（变量键 `"cell voltage"` 等保持原样）

## 5. 验证方案

### 5.1 静态验证

每改完一项，运行：
```bash
julia -e 'include("src/JuBat.jl"); using .JuBat'
```
确认模块加载无错误。

### 5.2 行为验证

完成所有改名后，运行现有示例脚本对比结果：
```bash
julia example/testexample.jl
```
预期：输出与改名前**完全一致**（同一随机种子、同一结果文件 hash）。

### 5.3 回归检查

用 `git diff` 逐项审查：
- 每个改动仅是 identifier rename
- 没有遗漏的引用
- 没有意外改动的逻辑

## 6. 风险评估

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| 改名遗漏导致 `UndefVarError` | 中 | 编译期立即暴露 | 模块加载即可发现 |
| `v = 0` 误改为循环计数器 | 中 | 行为可能改变 | 严格按识别规则筛选 |
| 字符串键误改 | 低 | 历史数据丢失 | 仅改标识符，键字符串保持 |
| 跨文件变量名冲突 | 低 | 编译期暴露 | 逐文件改、单文件验证 |

**总体风险**：🟢 低（纯改名 + 编译期可捕获错误）。

## 7. 后续工作

完成本次改名后，可独立推进：

1. **变量字典注册表重构**（`StandardVariables` + `create_element_workspace` 合并）
2. **路径 A 第二阶段**：478 处 `variables["..."]` 调用点替换为 const 引用
3. **热模型后端统一**（独立 spec）
4. **求解器编排重构**（独立 spec）

## 8. 不在范围内的事项明确排除

- ❌ 不引入新文件（如 `VarKeys.jl`）
- ❌ 不引入 typed struct（路径 B 已排除）
- ❌ 不改 `variables` 字典类型（保持 `Dict{String, Union{Array{Float64}, Float64}}`）
- ❌ 不改 `StandardVariables` / `create_element_workspace` 的内部逻辑（仅改它们的形参/局部变量名）
- ❌ 不调用 `copy_element_results` 的键列表（保持硬编码）

---

**附**：本 spec 的实施应严格遵循"每处修改须向用户询问"的要求，按 11 项顺序逐项推进，每项改动后向用户报告并确认再进入下一项。
