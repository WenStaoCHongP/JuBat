# Variables.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/Variables.jl`
> 状态: 修改 (105→185 行, +76%)
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
    variables["thermal2D temperature"] = zeros(Float64, ne, num)
    variables["thermal2D temperature at nodes"] = zeros(Float64, nT, num)
    variables["thermal2D temperature history"] = zeros(Float64, ne, num)
    variables["thermal2D temperature  at nodes history"] = zeros(Float64, nT, num)
    variables["heat_source_fields"] = zeros(Float64, ne, num)
    # ... 11 个分层热源 ...
    variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
    variables["thermal2D q_rev_ne"] = zeros(Float64, ne, num)
    variables["thermal2D q_ohm_s_ne"] = zeros(Float64, ne, num)
    variables["thermal2D q_ohm_e_ne"] = zeros(Float64, ne, num)
    variables["thermal2D q_sp"] = zeros(Float64, ne, num)
    variables["thermal2D q_rxn_pe"] = zeros(Float64, ne, num)
    variables["thermal2D q_rev_pe"] = zeros(Float64, ne, num)
    variables["thermal2D q_ohm_s_pe"] = zeros(Float64, ne, num)
    variables["thermal2D q_ohm_e_pe"] = zeros(Float64, ne, num)
    variables["thermal2D q_pcc"] = zeros(Float64, ne, num)
    variables["thermal2D q_ncc"] = zeros(Float64, ne, num)
    # ... 8 个单元电化学变量 ...
    variables["thermal2D element current"] = zeros(Float64, ne, num)
    variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
    variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
    variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
    variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
    variables["thermal2D element soc_n"] = zeros(Float64, ne, num)
    variables["thermal2D element soc_p"] = zeros(Float64, ne, num)
    variables["thermal2D element voltages"] = zeros(Float64, ne, num)
    variables["thermal2D element OCV"] = zeros(Float64, ne, num)
    # ... 4 个截止检测 ...
    variables["thermal2D active_mask"] = zeros(Float64, ne, num)
    variables["thermal2D n_cutoff_elements"] = zeros(Float64, 1, num)
    variables["thermal2D nearest_cutoff_element"] = zeros(Float64, 1, num)
    variables["thermal2D nearest_cutoff_ocv"] = zeros(Float64, 1, num)
    variables["thermal2D margin_to_cutoff"] = zeros(Float64, 1, num)
    # ... 5 个应力/应变 ...
    variables["thermal2D element thermal stress"] = zeros(Float64, ne, num)
    variables["thermal2D element diffusion stress"] = zeros(Float64, ne, num)
    variables["thermal2D element total stress"] = zeros(Float64, ne, num)
    variables["thermal2D element diffusion strain"] = zeros(Float64, ne, num)
    variables["thermal2D element thermal strain"] = zeros(Float64, ne, num)
    variables["total heat source"] = zeros(Float64, 1, num)
    variables["thermal2D displacement x"] = zeros(Float64, nT, num)
    variables["thermal2D displacement y"] = zeros(Float64, nT, num)
end
```

**CZM 块** (行 136-139, 2 个变量):

```julia
if case.opt.czm_enabled == true
    variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
    variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
end
```

### 2.2 `Variable_update!` 完全重写 (行 143-185)

从 11 行简单 for 循环变为 43 行，增加：
- 动态数组扩展（当 `v` 超过预分配）
- 类型安全赋值（`Float64` vs `Array{Float64}` 分支）
- `hist_keys` Set 避免写入新键

---

## 3. 优化方案

### 3.1 约束

- **`StandardVariables` 前 86 行完全不动**（main 分支 SPM/SPMe/P2D 变量创建）
- 仅重构 distributed2D 块 (行 88-135) 和 CZM 块 (行 136-139)
- `Variable_update!` 可重构，但保留动态扩展功能

### 3.2 distributed2D 块提取为子函数

```julia
# ---- 主函数中，行 88 之后改为 ----
variables["temperature"] = zeros(Float64, length(case.index["temperature"]), num)

if case.opt.thermalmodel == "lumped"
    variables["thermal lumped internal heat"] = zeros(Float64, 1, num)
elseif case.opt.thermalmodel == "distributed2D"
    _add_distributed2d_variables!(variables, case, num)
end

if case.opt.czm_enabled
    _add_czm_variables!(variables, case, num)
end
return variables
```

### 3.3 `_add_distributed2d_variables!` 子函数

```julia
"""
添加 distributed2D 热模型所需的全部变量 (ne 个热单元)。
"""
function _add_distributed2d_variables!(variables::Dict, case::Case, num::Int)
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen

    # 温度场
    for (k, dims) in [
        ("thermal2D temperature",              ne),
        ("thermal2D temperature at nodes",     nT),
        ("thermal2D temperature history",      ne),
        ("thermal2D temperature  at nodes history", nT),
        ("T_nodes",                            nT),
    ]
        variables[k] = zeros(Float64, dims, num)
    end

    # 分层热源 (11 层)
    for layer in _HEAT_SOURCE_LAYERS
        variables["thermal2D q_$(layer)"] = zeros(Float64, ne, num)
    end
    variables["heat_source_fields"] = zeros(Float64, ne, num)
    variables["total heat source"] = zeros(Float64, 1, num)

    # 单元电化学变量 (9 个)
    for k in _ELEMENT_CHEMISTRY_VARIABLES
        variables[k] = zeros(Float64, ne, num)
    end

    # 截止检测 (4 个)
    for k in _CUTOFF_VARIABLES
        variables[k] = zeros(Float64, 1, num)
    end
    variables["thermal2D active_mask"] = zeros(Float64, ne, num)

    # 应力/应变 (5 个) + 位移 (2 个)
    for k in _STRAIN_VARIABLES
        variables[k] = zeros(Float64, ne, num)
    end
    for k in ("thermal2D displacement x", "thermal2D displacement y")
        variables[k] = zeros(Float64, nT, num)
    end
end

# 常量列表
const _HEAT_SOURCE_LAYERS = (
    "rxn_ne", "rev_ne", "ohm_s_ne", "ohm_e_ne", "sp",
    "rxn_pe", "rev_pe", "ohm_s_pe", "ohm_e_pe", "pcc", "ncc"
)

const _ELEMENT_CHEMISTRY_VARIABLES = (
    "thermal2D element current",
    "thermal2D eta_n_e", "thermal2D eta_p_e",
    "thermal2D dUdT_n_e", "thermal2D dUdT_p_e",
    "thermal2D element soc_n", "thermal2D element soc_p",
    "thermal2D element voltages", "thermal2D element OCV",
)

const _CUTOFF_VARIABLES = (
    "thermal2D n_cutoff_elements",
    "thermal2D nearest_cutoff_element",
    "thermal2D nearest_cutoff_ocv",
    "thermal2D margin_to_cutoff",
)

const _STRAIN_VARIABLES = (
    "thermal2D element thermal stress",
    "thermal2D element diffusion stress",
    "thermal2D element total stress",
    "thermal2D element diffusion strain",
    "thermal2D element thermal strain",
)
```

### 3.4 `_add_czm_variables!`

```julia
function _add_czm_variables!(variables::Dict, case::Case, num::Int)
    Nn = 1  # SPM/SPMe 单节点
    Np = 1
    if case.opt.model == "P2D"
        Nn = case.mesh["negative electrode"].nlen
        Np = case.mesh["positive electrode"].nlen
    end
    variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
    variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
end
```

### 3.5 `Variable_update!` 精简

```julia
function Variable_update!(variables_hist::Dict, variables::Dict, v::Int64)
    # 动态扩展
    _ensure_capacity!(variables_hist, v)

    # 赋值
    hist_keys = Set(keys(variables_hist))
    for (k, val) in pairs(variables)
        k in hist_keys || continue
        _write_variable!(variables_hist[k], val, v)
    end
    return variables_hist
end

function _ensure_capacity!(variables_hist::Dict, v::Int)
    for (k, hist) in pairs(variables_hist)
        if isa(hist, Matrix{Float64}) && size(hist, 2) < v
            expansion = max(1000, size(hist, 2) ÷ 2)
            variables_hist[k] = hcat(hist, zeros(size(hist, 1), expansion))
        end
    end
end

function _write_variable!(hist, val, v::Int)
    if isa(hist, Matrix{Float64})
        col = isa(val, Vector{Float64}) ? val : (isa(val, Matrix{Float64}) ? val[:, 1] : [val])
        n = min(length(col), size(hist, 1))
        hist[1:n, v] = col[1:n]
    elseif isa(hist, Float64) && isa(val, Float64)
        # 标量直接赋值（无时间列）
    end
end
```

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| `StandardVariables` 行数 | 141 | ~95 (基础86 + distributed2D调用 + czm调用) |
| distributed2D 变量行数 | 42 行连续 if | ~30 行（子函数） |
| 新增子函数 | 0 | 4 (`_add_distributed2d_variables!`, `_add_czm_variables!`, `_ensure_capacity!`, `_write_variable!`) |
| `Variable_update!` 行数 | 43 | ~25 |
| main 分支代码改动 | N/A | 零（前 86 行不动） |
