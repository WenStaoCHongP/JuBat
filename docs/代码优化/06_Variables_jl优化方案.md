# Variables.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/Variables.jl`
> 状态: 修改 (105→184 行, +75%)
> main 分支行数: 105

---

## 1. main 分支现状

### 1.1 `StandardVariables(case::Case, num::Int64)` — 94 行

设计模式：`Dict{String, Union{Array{Float64}, Float64}}`，按模型类型条件追加变量。

- SPM 基础变量: 32 个（浓度、应力、电位等）
- SPMe 追加: 7 个（电解液浓度）
- P2D 追加: 16 个（电位、Gauss 点变量）
- 温度: 1 个（有条件）

### 1.2 `Variable_update!(variables_hist, variables, v)` — 11 行

简单的 for 循环：

```julia
for i in var_list
    if haskey(variables_hist, i)
        variables_hist[i][:,v] = collect(variables[i])
    end
end
```

---

## 2. 当前分支变更

### 2.1 `StandardVariables` 新增变量块

**distributed2D 块** (行 94-135, 31 个变量):

```julia
if case.opt.thermalmodel == "distributed2D"
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    # 温度场 (5)
    variables["thermal2D temperature"] = zeros(Float64, ne, num)
    # ... 5 个温度相关 ...
    # 分层热源 (11)
    variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
    # ... 11 个热源 ...
    # 单元电化学 (9)
    variables["thermal2D element current"] = zeros(Float64, ne, num)
    # ... 9 个电化学 ...
    # 截止检测 (5)
    variables["thermal2D active_mask"] = zeros(Float64, ne, num)
    # ... 5 个截止 ...
    # 应力/应变 (7)
    variables["thermal2D element thermal stress"] = zeros(Float64, ne, num)
    # ... 7 个应力 + 位移 ...
end
```

**CZM 块** (行 136-139, 2 个变量):

```julia
if case.opt.czm_enabled == true
    variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
    variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
end
```

### 2.2 `Variable_update!` 完全重写 (行 143-184)

从 11 行简单 for 循环变为 43 行，增加：
- 动态数组扩展（当 `v` 超过预分配）
- 类型安全赋值（`Float64` vs `Array{Float64}` 分支）
- `hist_keys` Set 避免写入新键

---

## 3. 优化方案

### 3.1 约束

- **`StandardVariables` 前 86 行完全不动**（main 分支 SPM/SPMe/P2D 变量创建）
- distributed2D 块和 CZM 块保持内联，不提取为子函数
- `Variable_update!` 保持内联，不提取子函数

### 3.2 distributed2D 块保持内联

所有 distributed2D 和 CZM 变量创建逻辑保留在 `StandardVariables` 主函数内，仅做类型替换：

```julia
# 仅修改：case.multi_spme_layout → case.layout
# 其余保持不变
```

### 3.3 `Variable_update!` 保持内联

43 行的 `Variable_update!` 保留原样，不拆分。动态数组扩展和类型安全赋值逻辑均保持内联。

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| `StandardVariables` 行数 | 141 | 141（不变） |
| distributed2D 变量 | 42 行内联 | 42 行内联（不变） |
| 新增子函数 | 0 | 0（不新增） |
| `Variable_update!` 行数 | 43 | 43（不变） |
| main 分支代码改动 | N/A | 零（前 86 行不动） |
