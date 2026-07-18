# CZM 渐进式有效面积损失 实现计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在分流求解器中引入损伤调制面积权重，使 D > D_threshold 的热单元承载更少电流，产生渐进式容量损失和正反馈。

**Architecture:** 唯一干预点是 `solve_branch_currents` 的权重计算。新增两个 Option 字段控制开关和阈值，新增两个辅助函数（面积因子、D 映射）。权重修改后 I_e 分配自动变化，传导到 SPMe 电流密度和热源，无需修改其他模块。

**Tech Stack:** Julia, 现有 JuBat 框架

**Spec:** `docs/superpowers/specs/2026-05-26-czm-progressive-area-loss-design.md`

---

## Chunk 1: 核心实现

### Task 1: Option 结构体新增字段

**Files:**
- Modify: `src/Option.jl` (czm_visc_tau 字段之后)

- [ ] **Step 1: 在 Option 结构体中添加两个新字段**

在 `src/Option.jl` 的 `czm_visc_tau` 字段之后（约第 82 行 `end` 之前）添加：

```julia
    # CZM progressive area loss (渐进式有效面积损失)
    czm_area_loss_enabled::Bool = false      # 启用渐进式面积损失（D > threshold 时缩减有效面积）
    czm_area_loss_threshold::Float64 = 0.83  # 面积开始缩减的损伤阈值
```

- [ ] **Step 2: 验证语法正确**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e 'include("src/JuBat.jl"); using .JuBat; opt = Option(); println(opt.czm_area_loss_enabled, " ", opt.czm_area_loss_threshold)'`

Expected: `false 0.83`

- [ ] **Step 3: Commit**

```bash
git add src/Option.jl
git commit -m "feat: add czm_area_loss_enabled and czm_area_loss_threshold to Option"
```

---

### Task 2: effective_area_factor 辅助函数

**Files:**
- Modify: `src/Materialmatrix.jl` (`compute_all_gap_conductances` 函数之后)

- [ ] **Step 1: 在 Materialmatrix.jl 末尾添加函数**

在 `compute_all_gap_conductances` 函数（约第 405 行）之后添加：

```julia
"""
	effective_area_factor(D::Float64, D_threshold::Float64) -> Float64

计算热单元的有效面积比例因子。

当 D ≤ D_threshold 时返回 1.0（无缩减）；
当 D > D_threshold 时线性缩减至 D=1.0 时为 0.0。

公式: factor = (1 - D) / (1 - D_threshold)
"""
function effective_area_factor(D::Float64, D_threshold::Float64)
	D ≤ D_threshold && return 1.0
	return (1.0 - D) / (1.0 - D_threshold)
end
```

- [ ] **Step 2: 验证函数行为**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e 'include("src/JuBat.jl"); using .JuBat; @assert JuBat.effective_area_factor(0.0, 0.83) == 1.0; @assert JuBat.effective_area_factor(0.83, 0.83) == 1.0; @assert abs(JuBat.effective_area_factor(0.96, 0.83) - 0.2353) < 0.001; @assert abs(JuBat.effective_area_factor(0.99, 0.83) - 0.0588) < 0.001; @assert JuBat.effective_area_factor(1.0, 0.83) == 0.0; println("All assertions passed")'`

Expected: `All assertions passed`

- [ ] **Step 3: Commit**

```bash
git add src/Materialmatrix.jl
git commit -m "feat: add effective_area_factor for progressive area loss"
```

---

### Task 3: map_czm_damage_to_thermal 映射函数

**Files:**
- Modify: `src/CallModel.jl` (CZM 失效处理代码块之前)

- [ ] **Step 1: 在 CallModel.jl 的 CZM 失效处理代码之前添加函数**

在 `CallModel.jl` 第 60 行（`# 获取CZM失效单元列表` 注释之前）添加：

```julia
"""
	map_czm_damage_to_thermal(czm_mesh, geometry, ne) -> Vector{Float64}

将 CZM damage_states 映射为热单元级别的 D 数组。

假设：一个内聚力单元只影响对应的内侧热单元（1-to-1）。
仅遍历 is_inner_layer[e] = true 的内侧热单元。
"""
function map_czm_damage_to_thermal(czm_mesh, geometry, ne)
	D_elem = zeros(ne)
	for e in 1:ne
		czm_indices = get(geometry.czm_element_map, Int64[])
		if !isempty(czm_indices) && geometry.is_inner_layer[e]
			# 1-to-1: 内侧热单元只对应一个 CZM 单元，直接取其 D
			D_elem[e] = czm_mesh.damage_states[czm_indices[1]].D
		end
	end
	return D_elem
end
```

- [ ] **Step 2: 验证语法正确（模块加载无报错）**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e 'include("src/JuBat.jl"); println("Loaded OK")'`

Expected: `Loaded OK`

- [ ] **Step 3: Commit**

```bash
git add src/CallModel.jl
git commit -m "feat: add map_czm_damage_to_thermal for CZM→thermal D mapping"
```

---

### Task 4: 修改分流求解器权重计算

**Files:**
- Modify: `src/Parallelsolution.jl` (函数签名和权重计算部分)
- Modify: `src/CallModel.jl` (调用 solve_branch_currents 的位置)

- [ ] **Step 1: 修改 solve_branch_currents 函数签名，新增关键字参数**

在 `src/Parallelsolution.jl` 第 358 行，修改函数签名：

**当前**：
```julia
function solve_branch_currents(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing; deactivated_elements::Union{Nothing,Vector{Int64}}=nothing)
```

**修改为**：
```julia
function solve_branch_currents(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, yt::Array{Float64}, t::Float64, I_total::Float64, areas::Vector{Float64}, Te_prev::Vector{Float64}, x_prev::Union{Nothing,Vector{Float64}}=nothing; deactivated_elements::Union{Nothing,Vector{Int64}}=nothing, D_elem::Union{Nothing,Vector{Float64}}=nothing)
```

- [ ] **Step 2: 修改权重计算逻辑**

在函数体内（约第 360 行），替换权重计算：

**当前**：
```julia
    ne = length(areas)
    w = areas ./ sum(areas)
```

**修改为**：
```julia
    ne = length(areas)
    # 渐进式有效面积损失：损伤调制权重
    if case.opt.czm_area_loss_enabled && D_elem !== nothing
        A_eff = areas .* effective_area_factor.(D_elem, case.opt.czm_area_loss_threshold)
        w = A_eff ./ sum(A_eff)
    else
        w = areas ./ sum(areas)
    end
```

- [ ] **Step 3: 修改 CallModel.jl 中的调用点**

在 `src/CallModel.jl` 第 81 行（调用 `solve_branch_currents` 的位置）：

**当前**：
```julia
    variables, I_e, Vc = solve_branch_currents(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev; deactivated_elements=deactivated_elements)
```

**修改为**：
```julia
    # 计算渐进式面积损失的 D 映射（仅在启用时）
    D_elem_area_loss = nothing
    if case.opt.czm_area_loss_enabled && case.czm_mesh !== nothing && geom !== nothing && hasfield(typeof(geom), :czm_element_map)
        D_elem_area_loss = map_czm_damage_to_thermal(case.czm_mesh, geom, ne)
    end

    variables, I_e, Vc = solve_branch_currents(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev; deactivated_elements=deactivated_elements, D_elem=D_elem_area_loss)
```

注意：这段代码应放在 `deactivated_elements` 计算之后、`solve_branch_currents` 调用之前。即在第 78 行 `end` 之后、第 80 行 `t_branch_ns = time_ns()` 之前插入。

- [ ] **Step 4: 验证语法正确**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e 'include("src/JuBat.jl"); println("Loaded OK")'`

Expected: `Loaded OK`

- [ ] **Step 5: Commit**

```bash
git add src/Parallelsolution.jl src/CallModel.jl
git commit -m "feat: integrate progressive area loss into branch current solver weights"
```

---

### Task 5: 回归验证

**Files:**
- No new files

- [ ] **Step 1: 验证默认行为（czm_area_loss_enabled=false）不变**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e '
include("src/JuBat.jl")
using .JuBat
opt = Option()
opt.model = "SPMe"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true
opt.czm_enabled = true
opt.czm_area_loss_enabled = false  # 默认关闭
param = JuBat.ChooseCell("Jellyroll")
case = JuBat.SetCase(param, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=20, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param)
case.czm_mesh = czm_mesh
opt.time = [0.0, 10.0]
opt.dt = [1.0, 1.0]
result = JuBat.Solve(case)
println("Voltage: ", round(result["cell voltage [V]"][end], digits=3))
println("Regression OK")
'`

Expected: 正常输出电压值，无报错

- [ ] **Step 2: 验证启用后行为变化**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e '
include("src/JuBat.jl")
using .JuBat
opt = Option()
opt.model = "SPMe"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true
opt.czm_enabled = true
opt.czm_area_loss_enabled = true  # 启用渐进面积损失
opt.czm_area_loss_threshold = 0.83
param = JuBat.ChooseCell("Jellyroll")
case = JuBat.SetCase(param, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=20, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param)
case.czm_mesh = czm_mesh
opt.time = [0.0, 10.0]
opt.dt = [1.0, 1.0]
result = JuBat.Solve(case)
println("Voltage: ", round(result["cell voltage [V]"][end], digits=3))
println("Enabled OK")
'`

Expected: 正常输出电压值（可能与回归测试不同），无报错

- [ ] **Step 3: Commit regression verification (if any fixes needed)**

```bash
git add -A
git commit -m "fix: address regression issues from progressive area loss integration"
```

（仅在需要修复时执行此步骤）

---

## Chunk 2: 文档与参数记录

### Task 6: 更新 CLAUDE.md 和导出

**Files:**
- Modify: `CLAUDE.md` (Option 表格)
- Modify: `src/JuBat.jl` (export 行，如果有的话)

- [ ] **Step 1: 在 CLAUDE.md 的机械/CZM 选项表格中新增两行**

在 `CLAUDE.md` 第 5.3 节的 CZM 选项表格中添加：

```markdown
| `czm_area_loss_enabled` | false | 启用渐进式有效面积损失（D > threshold 时缩减面积） |
| `czm_area_loss_threshold` | 0.83 | 面积开始缩减的 D 阈值 |
```

- [ ] **Step 2: 确认导出（检查 JuBat.jl 是否需要新增 export）**

检查 `src/JuBat.jl` 中是否已有 `effective_area_factor` 和 `map_czm_damage_to_thermal` 的 export。如果没有，在 CZM 相关 export 块中添加：

```julia
export effective_area_factor, map_czm_damage_to_thermal
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md src/JuBat.jl
git commit -m "docs: add progressive area loss to CLAUDE.md and exports"
```

---

### Task 7: 最终提交整理

- [ ] **Step 1: 确认所有修改文件状态**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && git status`

Expected: working tree clean（所有修改已提交）

- [ ] **Step 2: 确认所有提交的 log**

Run: `cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && git log --oneline -6`

Expected: 看到 5-6 个新提交（Task 1-6 各一个）

- [ ] **Step 3: 更新设计文档状态**

将 `docs/superpowers/specs/2026-05-26-czm-progressive-area-loss-design.md` 的状态从 `审阅修正后` 改为 `已实现`。

```bash
git add docs/superpowers/specs/2026-05-26-czm-progressive-area-loss-design.md
git commit -m "docs: mark progressive area loss spec as implemented"
```
