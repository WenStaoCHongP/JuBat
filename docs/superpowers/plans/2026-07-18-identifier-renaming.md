# 标识符改名重构 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 完成 spec `docs/superpowers/specs/2026-07-07-identifier-renaming-design.md` 中的 13 项标识符改名 + 属性别名内联，行为零变化。

**Architecture:** 纯机械改名重构。按"风险递增"顺序分 9 个独立 chunk：先单文件局部改名（chunk 1-3），再跨文件参数改名（chunk 4-8），最后是涉及 kwarg + struct 字段的复合改名（chunk 9）。每个 chunk 完成后做静态验证（模块加载 + grep 残留检查）并独立提交，便于回滚。

**Tech Stack:** Julia 1.x, JuBat 项目（电池建模框架）

**Spec 参考:** `docs/superpowers/specs/2026-07-07-identifier-renaming-design.md`

**用户约定:** 每处修改前先与用户确认 before/after，得到 OK 再动手。

**行号时效性声明:** 本 plan 行号以 2026-07-18 的 HEAD 为基准。执行前请用 Grep 重新定位每个标识符（`\bNAME\b`），按实际命中行号为准。chunk 内"约 N 处"的描述为估算，实际以 grep 为准。

---

## 全局约束

### 每步通用验证命令

**模块加载验证**（每个 chunk 结束时必跑）：
```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```
期望输出：`load ok`

**残留检查**（每个 chunk 结束时必跑，把 `<OLD_NAME>` 换成当前 chunk 改的旧名）：
```bash
# 用 Grep 工具搜 src/ 全目录
# 期望：仅在注释/字符串字面量中出现；无未改的标识符引用
```

### 提交规范

每个 chunk 单独提交，commit message 格式：
```
refactor(rename): chunk N — <旧名> → <新名>

按 spec §3.X 完成 <范围>。
行为零变化，模块加载通过。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
```

### 风险与回滚

- 每个 chunk 独立 commit，发现回归立即 `git revert HEAD`
- 若 julia 加载失败，先 `git diff HEAD~1` 定位遗漏的引用
- 复杂 chunk（7、8、9）若用户审查后退回，整体回滚到上一个通过的 chunk

---

## Chunk 1: CouplingState.jl `p` → `param` + `t_xx` 内联（Items 1, 11）

**Files:**
- Modify: `src/CouplingState.jl:228-240`

**目标**：消除 `p`/`t_pe`/`t_ne`/`t_sp`/`t_pcc`/`t_ncc` 6 个局部短别名，保留 `Σt`（3 次分母复用）。

- [ ] **Step 1: 与用户确认 before/after**

向用户展示 spec §3.1 的 Before/After 代码块，得到 OK。

- [ ] **Step 2: 执行改名**

用 Edit 工具替换 `src/CouplingState.jl:228-240` 的整个 `compute_effective_coating_modulus` 函数体为新版本（见 spec §3.1 After 块）。

- [ ] **Step 3: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```
期望：`load ok`

- [ ] **Step 4: 残留检查**

用 Grep 工具搜 `src/CouplingState.jl`：
- Pattern: `\bp\b\s*=\s*case\.param` → 期望 0 匹配
- Pattern: `\bt_pe\b|\bt_ne\b|\bt_sp\b|\bt_pcc\b|\bt_ncc\b` → 期望 0 匹配

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl
git commit -m "refactor(rename): chunk 1 — p→param, t_xx inline in compute_effective_coating_modulus

按 spec §3.1。行为零变化，模块加载通过。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 2: ThermalDistributed.jl `T_amb` 内联（Items 12, 13）

**Files:**
- Modify: `src/ThermalDistributed.jl:59` 和 `:185`（两个函数各一处）

**目标**：消除 `T_amb = param.cell.T_amb` 别名，引用直接用 `param.cell.T_amb`。

- [ ] **Step 1: 与用户确认**

展示 spec §3.8.2 的 Before/After，得到 OK。

- [ ] **Step 2: 修改第 59 行函数（apply_convection_bc）**

读取 `src/ThermalDistributed.jl:49-96`，识别：
- 删除第 59 行 `T_amb = param.cell.T_amb`
- 把函数体内 `fe1 += wt * T_amb * N1` 和 `fe2 += wt * T_amb * N2` 改成 `param.cell.T_amb`

- [ ] **Step 3: 修改第 185 行函数（apply_convection_bc!）**

读取 `src/ThermalDistributed.jl:178-...`，做同样改动。

- [ ] **Step 4: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 5: 残留检查**

Grep `src/ThermalDistributed.jl`：
- Pattern: `T_amb\s*=\s*param\.cell\.T_amb` → 期望 0 匹配
- Pattern: `\bT_amb\b` → 仅在注释中可能出现，无未改引用

- [ ] **Step 6: Commit**

```bash
git add src/ThermalDistributed.jl
git commit -m "refactor(rename): chunk 2 — T_amb alias inlined

按 spec §3.8.2。行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 3: CallModel.jl + Parallelsolution.jl 局部变量改名（Items 4, 5, 8）

**Files:**
- Modify: `src/CallModel.jl:51, 126, 132`（Te_prev, ws_pool, ws_e 局部）
- Modify: `src/Parallelsolution.jl:87, 90, 328, 343, 358, 386`（**Te_prev 跨文件函数形参**，spec §3.5 修订后纳入）

**目标**：`ws_e`/`ws_pool` 仅在 CallModel 局部改名；`Te_prev` 涉及 2 个跨文件函数形参（§4.3 例外）。

- [ ] **Step 1: 与用户确认**

展示三个改名点的 Before/After（参考 spec §3.5 + §3.3），**特别说明 `Te_prev` 现在跨 Parallelsolution.jl**。

- [ ] **Step 2: Te_prev → T_elem_prev（CallModel.jl）**

读取 `src/CallModel.jl:51-55` 及函数体内引用。把：
- 第 51 行：`Te_prev = zeros(Float64, ne)` → `T_elem_prev = zeros(Float64, ne)`
- 函数体内所有 `Te_prev[...]` / `Te_prev` 引用同步改（包括传给 `solve_branch_currents` 的实参 `Te_prev` 位置）

- [ ] **Step 3: Te_prev → T_elem_prev（Parallelsolution.jl，函数签名同步）**

读取 `src/Parallelsolution.jl`：
- 第 87 行（定义）：`function compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref)` → `..., T_elem_prev, ...)`
- 第 90 行：`coeffs[e] = compute_element_coefficients(e, Te_prev[e], ...)` → `T_elem_prev[e]`
- 第 328 行（调用）：`solve_branch_currents(..., Te_prev, x_prev; ...)` → `..., T_elem_prev, x_prev; ...)`
- 第 343 行（docstring）：`` - `Te_prev`: 各单元温度 `` → `` - `T_elem_prev`: 各单元温度 ``
- 第 358 行（定义）：`function solve_branch_currents(..., Te_prev::Vector{Float64}, x_prev...)` → `..., T_elem_prev::Vector{Float64}, ...)`
- 第 386 行（体内）：`coeffs = compute_all_coefficients(ne, Te_prev, ...)` → `T_elem_prev`

- [ ] **Step 4: ws_pool → workspace_pool（CallModel.jl）**

读取 `src/CallModel.jl:125-138`。把：
- `ws_pool = [...]` → `workspace_pool = [...]`
- `ws_pool[tid]` → `workspace_pool[tid]`

- [ ] **Step 5: ws_e → elem_workspace（CallModel.jl）**

同一函数体内：
- `ws_e = ws_pool[tid]` → `elem_workspace = workspace_pool[tid]`（结合 Step 4）
- `workspace=ws_e` → `workspace=elem_workspace`（SPMe_element 调用的 kwarg）

- [ ] **Step 6: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 7: 残留检查（全 src/，不只单文件）**

Grep `src/`：
- `\bTe_prev\b` → 0 匹配
- `\bws_pool\b` → 0
- `\bws_e\b` → 0

- [ ] **Step 8: Commit**

```bash
git add src/CallModel.jl src/Parallelsolution.jl
git commit -m "refactor(rename): chunk 3 — Te_prev跨文件 + ws_pool/ws_e 局部

按 spec §3.5（修订版含 Parallelsolution.jl）+ §3.3 第 4-5 项 + §4.3。
涉及 compute_all_coefficients 与 solve_branch_currents 函数签名。
行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

**注**：本 chunk 同步更新 spec §3.5，把 Te_prev 从"局部"改为"跨文件 + 2 个函数形参"。

---

## Chunk 4: `vars_e` → `elem_vars`（Item 6）

**Files:**
- Modify: `src/CallModel.jl:137, 167, 242, 248-270`
- Modify: `src/ThermalDistributed.jl:413`

**目标**：消除 `vars_e` 局部变量 + `copy_element_results(vars_e)` 函数定义形参（§4.3 例外条款）。

- [ ] **Step 1: 与用户确认**

展示 spec §3.4 第 6 项表格 + 函数签名变化。

- [ ] **Step 2: 改 copy_element_results 函数定义**

读取 `src/CallModel.jl:241-270`（实际函数定义在第 248 行；第 242 行是 docstring `copy_element_results(vars_e)`，可同步改但非代码）：
- 第 248 行定义：`function copy_element_results(vars_e)` → `function copy_element_results(elem_vars)`
- 函数体 18+ 处 `vars_e["..."]` → `elem_vars["..."]`

- [ ] **Step 3: 改 CallModel_MultiSPMe 内的局部 vars_e**

读取 `src/CallModel.jl:137` 和 `:167`：
- `variables_elems[e] = copy_element_results(vars_e)` → `... = copy_element_results(elem_vars)`
- 第 167 行附近：`vars_e = variables_elems[e]` → `elem_vars = variables_elems[e]`
- 后续 `vars_e["..."]` 引用同步改

- [ ] **Step 4: 改 ThermalDistributed.jl:413**

读取 `src/ThermalDistributed.jl:410-420`。把 `vars_e = variables_elems[e]` → `elem_vars = variables_elems[e]`，后续引用同步改。

- [ ] **Step 5: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 6: 残留检查**

Grep `src/`：
- `\bvars_e\b` → 0 匹配（含函数定义、调用、函数体）

- [ ] **Step 7: Commit**

```bash
git add src/CallModel.jl src/ThermalDistributed.jl
git commit -m "refactor(rename): chunk 4 — vars_e → elem_vars (含函数形参)

按 spec §3.4 第 6 项 + §4.3 例外条款。行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 5: `variables_hist` → `history_vars`（Item 7）

**Files:**
- Modify: `src/Variables.jl:236`（函数定义形参 + 函数体内 11 处引用）
- Modify: `src/Solve.jl:152, 182, 236, 366`（4 处调用/赋值）

**目标**：消除 `variables_hist` 跨文件别名，函数签名同步改（§4.3 例外）。

- [ ] **Step 1: 与用户确认**

展示 spec §3.4 第 7 项 + 实际 grep 出来的 15 处位置清单。

- [ ] **Step 2: 改 Variables.jl Variable_update! 函数**

读取 `src/Variables.jl:236-278`。
- 第 236 行签名：`function Variable_update!(variables_hist::Dict{...}, variables::Dict{...}, v::Int64)` → `function Variable_update!(history_vars::Dict{...}, variables::Dict{...}, v::Int64)`
- 函数体内所有 `variables_hist` → `history_vars`（约 11 处）

**注意**：第二个参数 `variables` 和第三个参数 `v` 保留原命名（spec §3.4 末段注释）。

- [ ] **Step 3: 改 Solve.jl 4 处调用**

读取 `src/Solve.jl:152, 182, 236, 366`：
- 第 152 行：`variables_hist = StandardVariables(case, num)` → `history_vars = StandardVariables(case, num)`
- 第 182 行：`Variable_update!(variables_hist, variables, v)` → `Variable_update!(history_vars, variables, v)`
- 第 236 行：同上
- 第 366 行：`result = PostProcessing(case, variables_hist, v)` → `result = PostProcessing(case, history_vars, v)`

**注**：Solve.jl 实际行号可能 ±5（spec 已知偏差），以 grep 定位为准。

- [ ] **Step 4: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 5: 残留检查**

Grep `src/`：
- `\bvariables_hist\b` → 0 匹配

- [ ] **Step 6: Commit**

```bash
git add src/Variables.jl src/Solve.jl
git commit -m "refactor(rename): chunk 5 — variables_hist → history_vars

按 spec §3.4 第 7 项 + §4.3。涉及 Variable_update! 函数签名。
行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 6: `fks` → `layer_weights`（Item 9）

**Files:**
- Modify: `src/ThermalDistributed.jl:15, 385, 477-488`（赋值 + 函数体内 12 处数组索引）
- Modify: `src/Materialmatrix.jl:10, 16, 20, 26, 30, 35`（2 个函数定义形参 + 调用方）

**目标**：跨文件同步消除 `fks` 别名 + 函数参数（§4.3 例外）。

- [ ] **Step 1: 与用户确认**

展示 spec §3.6 + 7 处具体行号清单。

- [ ] **Step 2: 改 Materialmatrix.jl 函数定义 + 调用**

读取 `src/Materialmatrix.jl` 全文（约 420 行）。
- 第 10 行：`function thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)` → `..., layer_weights, ...)`
- 第 16 行：函数体内 `fks` 引用 → `layer_weights`
- 第 17 行：`ne = size(fks, 1)` → `size(layer_weights, 1)`（**易漏**，单独列）
- 第 20 行：调用 `thermal_capacity_weights_2d(param, fks, ele_of_gp, wJ)` → 同步改
- 第 26 行：`function thermal_anisotropic_conductivity_2d(param, fks, ...)` → `..., layer_weights, ...)`
- 第 30 行：函数体内 `fks` 引用 → `layer_weights`
- 第 31 行：`ne = size(fks, 1)` → `size(layer_weights, 1)`（**易漏**，单独列）
- 第 35 行：调用同步改

- [ ] **Step 3: 改 ThermalDistributed.jl**

读取 `src/ThermalDistributed.jl:15, 385, 477-488`。
- 第 15 行：`fks = case.geometry !== nothing ? ...` → `layer_weights = ...`
- 第 385 行：同上
- 第 477-488 行：函数体内 12 处 `fks[e, 1]` / `fks[e, 2]` / ... → `layer_weights[e, 1]` 等

- [ ] **Step 4: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 5: 残留检查**

Grep `src/`：
- `\bfks\b` → 0 匹配

- [ ] **Step 6: Commit**

```bash
git add src/Materialmatrix.jl src/ThermalDistributed.jl
git commit -m "refactor(rename): chunk 6 — fks → layer_weights (跨文件)

按 spec §3.6 + §4.3。涉及 Materialmatrix.jl 两个函数签名。
行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 7: `mesh_th` → `mesh_thermal`（Item 2）

**Files:** 跨 6 文件，共约 36 处引用（行号已 grep 核实，2026-07-18）
- Modify: `src/CouplingState.jl:31, 34, 36, 37`（函数形参 + 3 处体内）
- Modify: `src/CallModel.jl:31, 53, 151`（赋值 + 2 处引用）
- Modify: `src/CsvExport.jl:256, 257, 260, 329, 330, 331, 332`（两个函数，7 处）
- Modify: `src/CycleData.jl:112, 113, 118, 119, 120, 151, 179, 269, 270, 271, 397, 398`（12 处）
- Modify: `src/Jellyrollmodel.jl:537, 540, 544, 549, 550, 563, 564`（7 处）
- Modify: `src/PostProcessing.jl:69, 70, 76`（3 处）

**目标**：跨 6 文件统一改名，含 `MultiSPMeLayout` 函数第 4 参数（§4.3 例外）。

- [ ] **Step 1: 与用户确认**

展示 spec §3.2 + 上述 36 处位置清单（按文件分组）。

- [ ] **Step 2: 改 CouplingState.jl 函数签名**

读取 `src/CouplingState.jl:30-41`。
- 第 31 行：`function MultiSPMeLayout(ne::Int, n_chem::Int, nT::Int, mesh_th)` → `..., mesh_thermal)`
- 第 34, 36, 37 行：函数体内 `mesh_th.gs.detJ`、`mesh_th.gs.ele`、`mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]` → 全部改 `mesh_thermal.xxx`

- [ ] **Step 3: 改 CallModel.jl（3 处）**

读取 `src/CallModel.jl:31, 53, 151`。
- 第 31 行：`mesh_th = case.mesh["thermal2D"]` → `mesh_thermal = ...`
- 第 53 行：`mesh_th.element[e, :]` → `mesh_thermal.element[e, :]`
- 第 151 行：`compute_heat_sources_with_czm(..., mesh_th)` 实参 → `mesh_thermal`

- [ ] **Step 4: 改 CsvExport.jl（7 处，两个函数）**

读取 `src/CsvExport.jl:256-260` 和 `:329-332`。每个 `mesh_th.xxx` 改为 `mesh_thermal.xxx`。

- [ ] **Step 5: 改 CycleData.jl（12 处）**

读取 `src/CycleData.jl:112-120, 151, 179, 269-271, 397-398`。逐处改。

- [ ] **Step 6: 改 Jellyrollmodel.jl（7 处）**

读取 `src/Jellyrollmodel.jl:537-564`。逐处改。

- [ ] **Step 7: 改 PostProcessing.jl（3 处）**

读取 `src/PostProcessing.jl:69-76`。逐处改。

- [ ] **Step 8: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 9: 残留检查**

Grep `src/`：
- `\bmesh_th\b` → 0 匹配

- [ ] **Step 10: Commit**

```bash
git add src/CouplingState.jl src/CallModel.jl src/CsvExport.jl src/CycleData.jl src/Jellyrollmodel.jl src/PostProcessing.jl
git commit -m "refactor(rename): chunk 7 — mesh_th → mesh_thermal (跨 6 文件, 36 处)

按 spec §3.2 + §4.3。涉及 MultiSPMeLayout 函数签名。
行为零变化。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 8: `T_nodes_carry` → `T_nodes_step`（Item 10）

**Files:** 跨 3 文件，共约 33 处引用（行号已 grep 核实）
- Modify: `src/CycleData.jl`：71, 74, 78, 81, 99, 109, 131, 149, 152, 154, 173, 201, 244, 252（**14 处**）
- Modify: `src/Solve.jl`：176, 223, 275, 419, 421, 434（**6 处**；421/434 是 CSV 输出路径）
- Modify: `src/CouplingState.jl`：278, 287, 301, 307, 308, 338, 348(注释), 354, 370, 377, 458, 462, 471（**13 处**，含 2 个函数定义形参）

**目标**：跨 3 文件统一改名，含 `compute_czm_strain_inputs` 和 `update_czm_damage!` 两个函数的形参（§4.3 例外）。

- [ ] **Step 1: 与用户确认**

展示 spec §3.7 + 上述 33 处位置清单（按文件分组）。

- [ ] **Step 2: 改 CouplingState.jl 两个函数签名 + 引用（13 处）**

读取 `src/CouplingState.jl:277-472`。

**compute_czm_strain_inputs**：
- 第 278 行（调用方）：`compute_czm_strain_inputs(case, variables, czm_mesh, T_nodes_carry)` → `..., T_nodes_step)`
- 第 287 行（定义）：`function compute_czm_strain_inputs(case, variables::Dict, czm_mesh, T_nodes_carry)` → `..., T_nodes_step)`
- 第 301, 307, 308 行（函数体内）：`length(T_nodes_carry)`、`T_nodes_carry[n]` 等 → `T_nodes_step`

**update_czm_damage!**：
- 第 338 行（调用）：`update_czm_damage!(case, variables, T_nodes_carry)` → `..., T_nodes_step)`
- 第 348 行（注释）：手动改
- 第 354 行（定义）：`function update_czm_damage!(case, variables, T_nodes_carry)` → `..., T_nodes_step)`
- 第 370 行（调用 compute_czm_strain_inputs）：传参
- 第 377 行：`has_nan_T = any(isnan, T_nodes_carry)` → `T_nodes_step`
- 第 458, 462, 471 行（6 参数兼容版）：`update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)` → `..., T_nodes_step, u_czm_prev)`，以及函数体内转发

- [ ] **Step 3: 改 CycleData.jl（14 处）**

逐行读取并改：
- 第 71 行：`T_nodes_carry = copy(T_nodes_init)` → `T_nodes_step = copy(T_nodes_init)`
- 第 74 行：`y0[thermal_range] .= T_nodes_carry` → 同步改
- 第 78, 81 行：两个 `T_nodes_carry = y0[...]` 赋值
- 第 99 行：`T_max_phase = maximum(T_nodes_carry) * T_ref`
- 第 109 行：`T_nodes_prev_export = copy(T_nodes_carry)`
- 第 131 行：`case.param.cell.T0 = mean(T_nodes_carry)`
- 第 149, 152 行：两个 `T_nodes_carry = ...` 赋值
- 第 154 行：`T_max_current = maximum(T_nodes_carry) * T_ref`
- 第 173 行：`T_nodes_out = ... 0.5 .* (T_nodes_prev_export .+ T_nodes_carry)`
- 第 201 行：`T_nodes_prev_export = copy(T_nodes_carry)`
- 第 244 行：`result.T_mean_end = mean(T_nodes_carry) * T_ref`
- 第 252 行：`"T_nodes" => copy(T_nodes_carry)`

- [ ] **Step 4: 改 Solve.jl（6 处）**

- 第 176 行：`T_nodes_carry = get(variables, "thermal2D temperature at nodes", Float64[])` → `T_nodes_step = ...`
- 第 223 行：`T_nodes_carry = T_nodes` → `T_nodes_step = T_nodes`
- 第 275 行：`case, variables, T_nodes_carry` → `case, variables, T_nodes_step`
- 第 419 行：`if case.opt.per_element_spme && ... && !isempty(T_nodes_carry)` → 同步改
- 第 421 行：`result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref` → 同步改（**输出键的值**，影响 CSV 结果）
- 第 434 行：`"T_nodes" => copy(T_nodes_carry),` → 同步改（**输出字典**）

- [ ] **Step 5: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

- [ ] **Step 6: 残留检查**

Grep `src/`：
- `\bT_nodes_carry\b` → 0 匹配

- [ ] **Step 7: Commit**

```bash
git add src/CouplingState.jl src/CycleData.jl src/Solve.jl
git commit -m "refactor(rename): chunk 8 — T_nodes_carry → T_nodes_step (跨 3 文件, 33 处)

按 spec §3.7 + §4.3。涉及 2 个 CZM 函数签名。
行为零变化（注意 Solve.jl:421/434 是输出键的值，键名字符串保持不变）。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 9: `ws` → `workspace`（Item 3，含 CZM kwarg + struct 字段）

**Files:** 跨 4 文件，**约 95 处**引用（行号时效性见 plan 头部声明，**执行前务必用 grep 重新定位**）
- Modify: `src/Variables.jl`：`create_element_workspace` 函数（赋值 + 40+ 处 `ws["..."]`）
- Modify: `src/CouplingState.jl:170-187`：`CZMAssemblyCache` struct 字段定义 `ws::CZMAssemblyWorkspace`
- Modify: `src/czm.jl`：`assemble_czm_system`、`assemble_coupled_system` 等（kwarg + 38 处引用，含 `cache.ws`、`ws.f_int_coh`、`ws.K_coh`、`ws.u_e` 等字段访问）
- Modify: `src/CzmSolve.jl`：多个 Newton 求解函数（kwarg + 16 处引用，含 `ws_basic`、`cache.ws` 等）

**目标**：消除 `ws` 在 CZM 模块的三重身份（自由变量 / kwarg / struct 字段），同步改名保证一致（§4.3 例外）。

**重要**：本 chunk 涉及面最广，建议拆为 4 个子提交（每文件一个），便于回滚定位。

- [ ] **Step 1: 与用户确认**

展示 spec §3.3 第 3a/3b 项 + 强调：
- `CZMAssemblyCache.ws` struct 字段要改（影响 `cache.ws` 访问）
- 局部变量 `ws_basic`、`empty_ws` 是否同步改（建议保留——它们是临时构造的 workspace 实例，名字含语义，可保留或同步改，由用户决定）

- [ ] **Step 2: 改 Variables.jl（子提交 1）**

读取 `src/Variables.jl:157-234`（`create_element_workspace` 函数）。
- 第 164 行：`ws = Dict{...}()` → `workspace = Dict{...}()`
- 第 233 行：`return ws` → `return workspace`（**易漏**：不含 `ws[`，`replace_all` 不会覆盖）
- 函数体内所有 `ws["..."] = ...` → `workspace["..."] = ...`（约 40 处）

**实施技巧**：用 Edit 工具的 `replace_all=true`，把 `ws[` 全部替换为 `workspace[`。先确认函数体内 `ws[` 都是指本变量。**替换完后再手工补 164 行和 233 行**（这两行不含 `ws[`），最后 grep `\bws\b` 确认 Variables.jl 内 0 残留。

子提交：
```bash
git add src/Variables.jl
git commit -m "refactor(rename): chunk 9a — ws → workspace in Variables.jl"
```

- [ ] **Step 3: 改 CouplingState.jl struct 字段（子提交 2）**

读取 `src/CouplingState.jl:170-187`。
- 第 176 行：`ws::CZMAssemblyWorkspace` → `workspace::CZMAssemblyWorkspace`
- 第 183 行：`empty_ws = CZMAssemblyWorkspace(0, 0)` → **保留**（局部变量，名字含语义）
- 第 184 行：构造 `new(..., empty_ws, ...)` → 字段位置参数仍传 `empty_ws`，无需改

子提交：
```bash
git add src/CouplingState.jl
git commit -m "refactor(rename): chunk 9b — CZMAssemblyCache.ws field → workspace"
```

**注意**：改完 struct 字段后，czm.jl 和 CzmSolve.jl 中的 `cache.ws` 访问会立即编译失败——下一步必须紧跟。

- [ ] **Step 4: 改 czm.jl（子提交 3）**

读取 `src/czm.jl` 多段。**关键锚点**（执行前 grep `\bws\b` 重新核对）：
- 函数定义中 `ws::Union{Nothing, CZMAssemblyWorkspace}=nothing` kwarg（`assemble_czm_system`、`assemble_coupled_system` 等）
- 函数体内 `ws = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(...)` 赋值
- `cache.ws = CZMAssemblyWorkspace(...)` 字段写入（如第 540 行）
- `ws.K_coh`、`ws.f_int_coh`、`ws.u_e`、`ws.B_global`、`ws.B_local`、`ws.δ_local`、`ws.BL_dT`、`ws.BL_dT_B`、`ws.T_vec`、`ws.BLtT`、`ws.K_e`、`ws.f_int_e`、`ws.separations`、`ws.tractions` 等字段访问
- 调用 `geom_cache=geom_cache, ws=ws, visc_beta=visc_beta` 形式

**实施建议**：
1. 先 grep `\bws\b` 取 czm.jl 全部行号
2. 按函数分段 Edit（assemble_czm_system / assemble_coupled_system / 等）
3. 每段改完跑模块加载验证

子提交：
```bash
git add src/czm.jl
git commit -m "refactor(rename): chunk 9c — ws → workspace in czm.jl"
```

- [ ] **Step 5: 改 CzmSolve.jl（子提交 4）**

读取 `src/CzmSolve.jl` 多段。**关键锚点**：
- 函数定义 kwarg：`backtrack_line_search!`、Newton 迭代函数等
- `ws_basic = cache !== nothing ? cache.ws : CZMAssemblyWorkspace(...)` 局部赋值（**`ws_basic` 是否同步改由 Step 1 用户决定**）
- `cache.ws` 字段读取（如第 157, 252, 471 行）
- 调用 `assemble_coupled_system(...; ws=ws_basic, ...)` 或 `ws=ws`

子提交：
```bash
git add src/CzmSolve.jl
git commit -m "refactor(rename): chunk 9d — ws → workspace in CzmSolve.jl"
```

- [ ] **Step 6: 静态验证**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("load ok")'
```

若报 `Setfield` 或 `type CZMAssemblyCache has no field ws` 类错误，说明有遗漏的 `cache.ws` 访问，定位补改。

- [ ] **Step 7: 残留检查（修正 grep 模式）**

Grep `src/`，**用 word boundary 避免误报**：
- `\bws\b`（word boundary）→ 期望 0 匹配（除注释）
- `\.ws\b`（字段访问旧名）→ 期望 0
- `ws::Union{Nothing, CZMAssemblyWorkspace}` → 0
- `geom_cache=geom_cache, ws=` → 0（旧 kwarg 调用形式）
- `ws_basic` → **不期望 0**（如果用户 Step 1 选保留），这是允许的
- `empty_ws` → **不期望 0**（局部变量，保留）

**正向验证**：
- `\bworkspace\b` → 大量匹配（新名生效）
- `\.workspace\b`（字段访问新名）→ 大量匹配
- `workspace::Union{Nothing, CZMAssemblyWorkspace}` → 多处

- [ ] **Step 8: 行为验证（实际 Solve 调用，不只是 SetCase）**

```bash
julia --project=. -e '
include("src/JuBat.jl")
using .JuBat
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.model = "SPMe"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true
opt.czm_enabled = true
opt.mechanicalmodel = "full"
opt.time = [0, 1.0]   # 1 秒仿真，触发 CZM 实际调用
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
result = JuBat.Solve(case)
println("solve ok, V=$(result["cell voltage"][end])")
'
```
期望：`solve ok, V=...`，**实际跑过 CZM 求解路径**，确认字段访问无遗漏。

若 `Solve` 因参数缺失而失败（如 mesh 未正确设置），可降级回 SetCase 验证，但要在 commit 中标注"未做端到端验证"。

- [ ] **Step 9: 合并 commit 信息（可选）**

若希望 git log 简洁，可用 `git rebase -i HEAD~4` 把 4 个子提交合并为一个 `chunk 9`。否则保留 4 个子提交便于追溯。

---

## 最终验证

---

## 最终验证

- [ ] **Step F1: 跑全部 chunk 完成后的 grep**

Grep `src/`，确认以下旧名全部 0 匹配（除字符串字面量与注释）：
- `\bp\b\s*=\s*case\.param`
- `\bmesh_th\b`
- `\bws_e\b`, `\bws_pool\b`
- `\bvars_e\b`, `\bvariables_hist\b`
- `\bTe_prev\b`
- `\bfks\b`
- `\bT_nodes_carry\b`
- `\.ws\b`（字段访问）
- `\bt_pe\b`, `\bt_ne\b`, `\bt_sp\b`, `\bt_pcc\b`, `\bt_ncc\b`
- `T_amb\s*=\s*param\.cell\.T_amb`

- [ ] **Step F2: 模块加载 + smoke**

```bash
julia --project=. -e 'include("src/JuBat.jl"); using .JuBat; println("final load ok")'
```

- [ ] **Step F3: git log 确认 9 个 chunk 提交完整**

```bash
git log --oneline -10
```

期望：9 个 `refactor(rename): chunk N — ...` commit，从 1 到 9。

---

## 已知风险与降级

| 风险 | 触发条件 | 降级方案 |
|------|----------|----------|
| Chunk 9 字段访问遗漏 | `type has no field ws` 错误 | 全文 grep `\.ws` 逐处补改 |
| Chunk 8 函数签名链断 | `UndefVarError` 在 CZM 调用链 | 检查 6 参数 `update_czm_damage!` 是否漏改 |
| Chunk 5 grep 残留 | `variables_hist` 还在 Solve.jl | 实际行号可能 ±5，按 grep 重定位 |
| Chunk 7 跨文件漏改 | 6 文件之一遗漏 | grep `\bmesh_th\b` 全仓库扫描 |
| 任意 chunk 用户审查退回 | 用户要求调整 | `git revert HEAD`，按反馈修订 spec 后重做 |

**整体回滚**：若全计划放弃，`git revert HEAD~9..HEAD`（仅在有 9 个 chunk commit 的前提下）。

---

## 用户确认节奏

按用户"每处修改须向我询问"的要求，每个 chunk 在 Step 1（与用户确认）处暂停，等待 OK 再进入 Step 2。复杂 chunk（7、8、9）可拆为多段子任务，每段单独确认。

实施完成后通知用户做最终审查。
