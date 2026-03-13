# heat_Q 归一化修正实现计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修正 `heat_Q` 归一化公式，添加缺失的密度因子，确保密度和比热容正确分离归一化。

**Architecture:** 在 Scale struct 添加 `rho` 字段 → 修正五层热参数归一化公式 → 修改 Materialmatrix.jl 使用 `rho * heat_Q` 计算体积热容。

**Tech Stack:** Julia, JuliaUnit testing

**Spec Document:** `docs/superpowers/specs/2026-03-13-heat-Q-normalization-fix-design.md`

---

## File Structure

| 文件 | 操作 | 职责 |
|-----|------|------|
| `src/SetParams.jl` | 修改 | Scale struct 添加 rho 字段；修正归一化公式 |
| `src/Materialmatrix.jl` | 修改 | 体积热容计算改为 rho * heat_Q |
| `md/01_参数定义与归一化.md` | 修改 | 文档清理，统一术语 |
| `docs/plans/2026-03-13-code-implementation-plan.md` | 修改 | 删除/修正错误描述 |

---

## Chunk 1: Scale Struct 添加 rho 字段

### Task 1.1: 添加 Scale.rho 字段

**Files:**
- Modify: `src/SetParams.jl:175-212` (Scale struct 定义)

- [ ] **Step 1: 定位 Scale struct 定义**

在 `src/SetParams.jl` 中找到 Scale struct 定义（约第 175 行）。

- [ ] **Step 2: 添加 rho 字段**

在 Scale struct 中添加 `rho` 字段，放在热尺度参数区域：

```julia
@with_kw mutable struct Scale
    # ... existing fields ...
    E_p::Float64 = 0
    # --- Thermal scaling (统一能量尺度) ---
    rho::Float64 = 0          # 密度尺度 = 电池平均密度 [kg/m³]
    P_ref::Float64 = 0
    lambda::Float64 = 0          # 导热率尺度参数 P_ref/(L*T_ref) [W/(m·K)]
    # ... rest of fields ...
end
```

- [ ] **Step 3: 验证字段已存在**

确认 `ChooseCell` 函数中已有 `param_dim.scale.rho = param_dim.cell.rho`（约第 292 行）。此行无需修改。

- [ ] **Step 4: 提交**

```bash
git add src/SetParams.jl
git commit -m "feat(scale): add rho field to Scale struct for density normalization

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 2: 修正归一化公式

### Task 2.1: 修正 PE 比热容归一化

**Files:**
- Modify: `src/SetParams.jl` (NormaliseParam 函数)

- [ ] **Step 1: 定位 PE 归一化代码**

找到 `NormaliseParam` 函数中 PE 的 `heat_Q` 归一化（约第 336 行）。

当前代码:
```julia
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

- [ ] **Step 2: 修正 PE 归一化公式**

替换为:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
# 其中 P_ref = phi * I_typ
# 注: param.PE.rho 存储无量纲密度 ρ* = ρ/ρ_ref
# param.PE.heat_Q 存储无量纲比热容 c*
# 体积热容 (ρc)* = ρ* · c* 在使用时计算
param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                  (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetParams.jl
git commit -m "fix(normalize): correct PE heat_Q normalization with rho factor

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

### Task 2.2: 修正 NE 比热容归一化

**Files:**
- Modify: `src/SetParams.jl` (NormaliseParam 函数)

- [ ] **Step 1: 定位 NE 归一化代码**

找到 NE 的 `heat_Q` 归一化（约第 361 行）。

当前代码:
```julia
param.NE.heat_Q = param_dim.NE.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

- [ ] **Step 2: 修正 NE 归一化公式**

替换为:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
param.NE.heat_Q = param_dim.NE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                  (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetParams.jl
git commit -m "fix(normalize): correct NE heat_Q normalization with rho factor

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

### Task 2.3: 修正 SP 比热容归一化

**Files:**
- Modify: `src/SetParams.jl` (NormaliseParam 函数)

- [ ] **Step 1: 定位 SP 归一化代码**

找到 SP 的 `heat_Q` 归一化（约第 370 行）。

当前代码:
```julia
param.SP.heat_Q = param_dim.SP.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

- [ ] **Step 2: 修正 SP 归一化公式**

替换为:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
param.SP.heat_Q = param_dim.SP.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                  (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetParams.jl
git commit -m "fix(normalize): correct SP heat_Q normalization with rho factor

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

### Task 2.4: 修正 PCC 比热容归一化

**Files:**
- Modify: `src/SetParams.jl` (NormaliseParam 函数)

- [ ] **Step 1: 定位 PCC 归一化代码**

找到 PCC 的 `heat_Q` 归一化（约第 377 行）。

当前代码:
```julia
param.PCC.heat_Q = param_dim.PCC.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

- [ ] **Step 2: 修正 PCC 归一化公式**

替换为:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
param.PCC.heat_Q = param_dim.PCC.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                   (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetParams.jl
git commit -m "fix(normalize): correct PCC heat_Q normalization with rho factor

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

### Task 2.5: 修正 NCC 比热容归一化

**Files:**
- Modify: `src/SetParams.jl` (NormaliseParam 函数)

- [ ] **Step 1: 定位 NCC 归一化代码**

找到 NCC 的 `heat_Q` 归一化（约第 383 行）。

当前代码:
```julia
param.NCC.heat_Q = param_dim.NCC.heat_Q * param.scale.L^3 / (param.scale.lambda * param.scale.t0)
```

- [ ] **Step 2: 修正 NCC 归一化公式**

替换为:
```julia
# 比热容归一化: c* = c * ρ_ref * L³ * T_ref / (t0 * P_ref)
param.NCC.heat_Q = param_dim.NCC.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref /
                   (param.scale.t0 * param.scale.phi * param.scale.I_typ)
```

- [ ] **Step 3: 提交**

```bash
git add src/SetParams.jl
git commit -m "fix(normalize): correct NCC heat_Q normalization with rho factor

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 3: 修改体积热容计算

### Task 3.1: 修正 Materialmatrix.jl 体积热容计算

**Files:**
- Modify: `src/Materialmatrix.jl:16-23` (thermal_capacity_weights_2d 函数)

- [ ] **Step 1: 定位体积热容计算代码**

在 `src/Materialmatrix.jl` 中找到 `thermal_capacity_weights_2d` 函数（约第 16-23 行）。

当前代码:
```julia
function thermal_capacity_weights_2d(param::Params, fks::Matrix{Float64}, ele_of_gp::Vector{Int64}, wJ::Vector{Float64})
	ne = size(fks, 1)
	rho_c_e = zeros(Float64, ne)
	@inbounds for e in 1:ne
		rho_c_e[e] = fks[e, 1] * param.NE.rho + fks[e, 2] * param.SP.rho + fks[e, 3] * param.PE.rho + fks[e, 4] * param.PCC.rho + fks[e, 5] * param.NCC.rho
	end
	return rho_c_e[ele_of_gp] .* wJ
end
```

- [ ] **Step 2: 修正体积热容计算**

替换为:
```julia
function thermal_capacity_weights_2d(param::Params, fks::Matrix{Float64}, ele_of_gp::Vector{Int64}, wJ::Vector{Float64})
	ne = size(fks, 1)
	rho_c_e = zeros(Float64, ne)
	@inbounds for e in 1:ne
		# 体积热容 (ρc)* = ρ* · c*
		# param.layer.rho = 无量纲密度 ρ*
		# param.layer.heat_Q = 无量纲比热容 c*
		rho_c_e[e] = fks[e, 1] * (param.NE.rho * param.NE.heat_Q) +
		             fks[e, 2] * (param.SP.rho * param.SP.heat_Q) +
		             fks[e, 3] * (param.PE.rho * param.PE.heat_Q) +
		             fks[e, 4] * (param.PCC.rho * param.PCC.heat_Q) +
		             fks[e, 5] * (param.NCC.rho * param.NCC.heat_Q)
	end
	return rho_c_e[ele_of_gp] .* wJ
end
```

- [ ] **Step 3: 提交**

```bash
git add src/Materialmatrix.jl
git commit -m "fix(thermal): correct volumetric heat capacity calculation

Use rho * heat_Q to compute (ρc)* = ρ* · c* instead of incorrectly
treating param.layer.rho as volumetric heat capacity.

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 4: 文档修改

### Task 4.1: 更新 01_参数定义与归一化.md

**Files:**
- Modify: `md/01_参数定义与归一化.md`

- [ ] **Step 1: 审查文档矛盾点**

阅读文档，识别以下问题:
1. 体积热容被当作独立字段的描述
2. 比热容归一化公式缺少密度因子
3. `param.PE.rho` 被描述为存储体积热容的错误描述

- [ ] **Step 2: 修正比热容归一化公式描述**

确保文档公式与代码一致:
$$c^* = c \cdot \rho_{ref} \cdot \frac{L^3 \cdot T_{ref}}{t_0 \cdot P_{ref}}$$

- [ ] **Step 3: 明确字段语义**

更新表格和描述，明确:
- `param.layer.rho` = 无量纲密度 $\rho^*$
- `param.layer.heat_Q` = 无量纲比热容 $c^*$
- $(\rho c)^* = \rho^* \cdot c^*$ 在使用时计算

- [ ] **Step 4: 提交**

```bash
git add md/01_参数定义与归一化.md
git commit -m "docs: clarify rho and heat_Q semantics in normalization doc

- Remove contradictory volumetric heat capacity field descriptions
- Clarify that (ρc)* = ρ* · c* is computed at runtime
- Ensure formula matches code implementation

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

### Task 4.2: 更新 2026-03-13-code-implementation-plan.md

**Files:**
- Modify: `docs/plans/2026-03-13-code-implementation-plan.md`

- [ ] **Step 1: 定位错误描述**

找到关于 `param.PE.rho` 存储"体积热容"的错误描述。

- [ ] **Step 2: 添加取代说明**

在文档开头添加:
```markdown
> **注意**: 关于热参数归一化的描述已被 `docs/superpowers/specs/2026-03-13-heat-Q-normalization-fix-design.md` 取代。
> 如有矛盾，以设计文档为准。
```

- [ ] **Step 3: 提交**

```bash
git add docs/plans/2026-03-13-code-implementation-plan.md
git commit -m "docs: add superseded notice to implementation plan

Reference the new heat_Q normalization fix design spec.

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 5: 验证

### Task 5.1: 运行回归测试

**Files:**
- Test: 现有热模型验证算例

- [ ] **Step 1: 运行 SPMe-热耦合算例**

```bash
cd "D:\OneDrive\Desktop\Jubat For Cursor\JuBat"
julia --project=. example/电化学-热耦合验证/SPMe_Thermal_example.jl
```

Expected: 运行无错误，温度结果合理

- [ ] **Step 2: 验证无量纲参数量级**

检查输出中无量纲参数是否为 O(1) 量级:
- `param.PE.rho` 应约为 1.0
- `param.PE.heat_Q` 应约为 1.0
- `param.PE.rho * param.PE.heat_Q` 应约为 1.0

- [ ] **Step 3: 对比温度结果**

如果可能，与修改前的温度结果对比，误差应 < 0.1%

- [ ] **Step 4: 最终提交**

```bash
git add -A
git commit -m "test: verify heat_Q normalization fix

Regression test passed. Temperature results consistent.

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Summary

| Chunk | Tasks | Commits |
|-------|-------|---------|
| 1 | Scale.rho 字段 | 1 |
| 2 | 归一化公式 (5层) | 5 |
| 3 | Materialmatrix.jl | 1 |
| 4 | 文档 (2个文件) | 2 |
| 5 | 验证 | 1 |
| **Total** | **10 tasks** | **10 commits** |

**Prerequisites:**
- Chunk 1 必须先完成
- Chunk 2 依赖 Chunk 1
- Chunk 3 依赖 Chunk 2
- Chunk 4-5 可并行
