# Parallelsolution.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/Parallelsolution.jl`
> 状态: 新增 (619 行)

---

## 1. 现状分析

### 1.1 函数清单

| 函数 | 行范围 | 行数 | 职责 |
|------|--------|------|------|
| `_debug_check_prefactors` | 7-30 | 24 | 调试：检查预因子 NaN |
| `_debug_check_coefficients` | 33-52 | 20 | 调试：检查系数异常 |
| `_debug_check_initial_voltage` | 55-76 | 22 | 调试：检查初始电压 |
| `_compute_electrochemical_prefactors` | 83-116 | 34 | 计算交换电流密度、OCV、熵系数 |
| `_compute_element_coefficients` | 119-149 | 31 | 单元分流系数 (C1, C2, alpha, C5) |
| `_compute_all_coefficients` | 152-158 | 7 | 批量系数计算 |
| `_branch_voltage` | 162-166 | 5 | V_e = C1 + C2*I + alpha*T + C5*ln(I) |
| `_branch_dVdI` | 169-175 | 7 | dV_e/dI 解析导数 |
| `_initialize_currents` | 179-191 | 13 | 初始化电流猜测 |
| `_check_voltage_bounds` | 194-208 | 15 | 电压边界检查 |
| `_detect_cutoff_elements` | 226-287 | 62 | 截止电压检测 |
| `_newton_iteration!` | 294-388 | 95 | Newton-Raphson 主循环 |
| `_line_search` | 391-415 | 25 | 回溯线搜索 |
| `solve_branch_currents_newton` | 452-619 | 168 | 主入口 |

### 1.2 架构评价

**优点**：
- 辅助函数分解得好（`_compute_prefactors`, `_branch_voltage`, `_branch_dVdI`）
- Newton-Raphson 和线搜索逻辑清晰
- 内部函数统一 `_` 前缀

**问题**：
- `solve_branch_currents_newton` 仍有 168 行，含过多 active/inactive 分支逻辑
- `_detect_cutoff_elements` 62 行返回 7 字段 NamedTuple——可改为 struct
- 3 个 `_debug_check_*` 函数仅在 `debug_coupling=true` 时使用，可条件编译

---

## 2. 优化方案

### 2.1 约束

此文件是 Parameters_Design **新增**文件，无 main 分支代码，全部可重构。

### 2.2 `solve_branch_currents_newton` 拆分

```julia
# ===== 旧 (168 行): =====
function solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev;
                                       deactivated_elements=nothing)
    # 1. 提取物理量 (~20 行)
    # 2. 计算预因子 (~5 行)
    # 3. 计算系数 (~5 行)
    # 4. 截止检测 (~15 行)
    # 5. 初始化电流 (~5 行)
    # 6. Newton 迭代 (~5 行)
    # 7. 结果写入 variables (~40 行)
    # 8. 异常处理 (~20 行)
    # 9. 返回
end

# ===== 新 (~80 行): =====
function solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev;
                                deactivated_elements=nothing)
    # 1. 计算分流系数
    prefactors = compute_electrochemical_prefactors(variables, case.param,
                                                      case.mesh["negative electrode"],
                                                      case.mesh["positive electrode"])
    coeffs = compute_all_coefficients(case.layout.ne, Te_prev, case.param, prefactors,
                                       case.param_dim.scale.T_ref;
                                       debug_mode=case.opt.debug_coupling)

    # 2. 截止检测
    cutoff = detect_cutoff(coeffs, case.layout.ne, I_total, case.param_dim.scale.phi)
    active_mask = cutoff.active_mask
    if deactivated_elements !== nothing
        active_mask .&= .!deactivated_elements
    end

    # 3. 初始化 + 求解
    I_e = initialize_currents(case.layout.ne, areas, I_total, x_prev)
    w = vec(areas) ./ sum(areas)
    V_common = newton_solve!(I_e, coeffs, w, I_total; active_mask)

    # 4. 写回结果
    _write_branch_results!(variables, I_e, V_common, cutoff, coeffs)

    return variables, I_e, V_common
end
```

### 2.3 截止检测结果改用 struct

```julia
# 旧: 返回 NamedTuple，7 字段
# (cutoff_elements, cutoff_ocv, active_mask, n_cutoff, nearest_cutoff_element,
#  nearest_cutoff_ocv, margin_to_cutoff)

# 新:
struct CutoffInfo
    elements::Vector{Int}           # 已截止单元列表
    ocv::Vector{Float64}            # 对应 OCV
    active_mask::BitVector          # 活跃掩码
    n_cutoff::Int                   # 截止数量
    nearest_element::Int            # 最近截止单元
    nearest_ocv::Float64            # 最近截止 OCV
    margin::Float64                 # 到截止的裕度
end
```

### 2.4 主函数拆小后，保留的函数改名

| 旧名 | 新名 | 变化 |
|------|------|------|
| `solve_branch_currents_newton` | `solve_branch_currents` | 缩短，去掉 `_newton` |
| `_compute_electrochemical_prefactors` | `compute_electrochemical_prefactors` | 去掉 `_`（公开） |
| `_compute_element_coefficients` | `compute_element_coefficients` | 去掉 `_` |
| `_compute_all_coefficients` | `compute_all_coefficients` | 去掉 `_` |
| `_newton_iteration!` | `newton_solve!` | 缩短 |
| `_detect_cutoff_elements` | `detect_cutoff` | 缩短 |
| `_debug_check_*` | 保持不变 | 调试函数不改 |

### 2.5 写回结果提取为子函数

```julia
function _write_branch_results!(variables, I_e, V_common, cutoff::CutoffInfo, coeffs)
    ne = length(I_e)
    variables["thermal2D element current"] = I_e
    variables["thermal2D element voltages"] = [_branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
    variables["cell voltage"] = V_common
    variables["thermal2D active_mask"] = Float64.(cutoff.active_mask)
    variables["thermal2D n_cutoff_elements"] = Float64(cutoff.n_cutoff)
    variables["thermal2D nearest_cutoff_element"] = Float64(cutoff.nearest_element)
    variables["thermal2D nearest_cutoff_ocv"] = cutoff.nearest_ocv
    variables["thermal2D margin_to_cutoff"] = cutoff.margin
    return nothing
end
```

---

## 3. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| 总行数 | 619 | ~500 |
| 主函数行数 | 168 | ~80 |
| 截止检测返回 | NamedTuple | `CutoffInfo` struct |
| 内部函数前缀 | `_` (不一致) | 公开函数无前缀，内部保留 `_` |
