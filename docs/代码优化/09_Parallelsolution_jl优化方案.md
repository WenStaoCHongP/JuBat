# Parallelsolution.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/Parallelsolution.jl`
> 状态: 新增 (441 行)

---

## 1. 现状分析

### 1.1 函数清单（当前代码实际状态）

| 函数 | 行号 | 行数 | 职责 |
|------|------|------|------|
| `compute_prefactors` | 10 | 34 | 计算交换电流密度、OCV、熵系数 |
| `compute_element_coefficients` | 46 | 31 | 单元分流系数 (C1, C2, alpha, C5) |
| `compute_all_coefficients` | 74 | 7 | 批量系数计算 |
| `branch_voltage` | 84 | 5 | V_e = C1 + C2*I + alpha*T + C5*ln(I) |
| `branch_dVdI` | 91 | 7 | dV_e/dI 解析导数 |
| `initialize_currents` | 101 | 13 | 初始化电流猜测 |
| `check_voltage_bounds` | 116 | 15 | 电压边界检查 |
| `detect_cutoff_elements` | 148 | 62 | 截止电压检测 |
| `newton_iteration` | 216 | 95 | Newton-Raphson 主循环 |
| `line_search` | 313 | 25 | 回溯线搜索 |
| `solve_branch_currents` | 358-441 | 84 | 主入口 |

### 1.2 架构评价

**优点**：
- 辅助函数分解得好（`compute_prefactors`, `branch_voltage`, `branch_dVdI`）
- Newton-Raphson 和线搜索逻辑清晰
- 无 `_` 前缀

**问题**：
- `solve_branch_currents` 仍承担“截止检测 + 牛顿迭代调用 + 状态写回”等多职责，循环热点集中
- `detect_cutoff_elements` 已采用 `CutoffInfo`，仍可继续精简循环与写回路径

---

## 2. 优化方案

### 2.1 约束

此文件是 Parameters_Design **新增**文件，无 main 分支代码，全部可重构。
但**不新增函数**，仅做重命名和类型替换。

### 2.2 `solve_branch_currents` 命名与入口统一

```julia
# 当前统一入口:
solve_branch_currents(case, variables, yt, t, I_total, areas, Te_prev, x_prev;
                      deactivated_elements=nothing)
```

### 2.3 截止检测结果改用 struct

```julia
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

`detect_cutoff_elements` 返回值从 NamedTuple 改为 `CutoffInfo`。

### 2.4 结果写回逻辑保持内联

原方案提议提取 `write_branch_results!`——现在保持内联在 `solve_branch_currents` 中，不新增函数。

### 2.5 数据结构说明

`Parallelsolution.jl` 当前不直接依赖 `case.multi_spme_layout`，入口通过参数和 `variables` 字典传递状态。

### 2.6 循环嵌套优化参考（来自 11 文档）

本节补充 `11_其余文件优化方案汇总.md` 中 L/M 节对应到 `Parallelsolution.jl` 的可落地方案，重点覆盖：

1. 显式嵌套：`newton_iteration` 的 `iter -> e` 双层循环。
2. 隐式嵌套：`newton_iteration` 内调用 `line_search` 的 `attempt -> e` 双层循环。
3. 推导式/生成器：`mean([ ... for e in ...])`、`sum(... for e in ...)` 的隐式循环。

#### 2.6.1 显式+隐式嵌套：`newton_iteration` 与 `line_search` 对照

**旧逻辑（循环分散，重复计算较多）**

```julia
# newton_iteration 内
for iter in 1:max_iters
    for e in 1:ne
        V_e = branch_voltage(coeffs[e], I_e[e])
        F[e] = V_e - V
        dFdI[e] = branch_dVdI(coeffs[e], I_e[e])
    end

    # ... 计算 ΔV/ΔI ...

    λ, V_trial = line_search(I_e, V, ΔI, ΔV, I_trial, ne)
end

# line_search 内
for attempt in 1:max_attempts
    for e in 1:ne
        val = I_e[e] + λ * ΔI[e]
        # ... 校验 ...
    end
end
```

**新逻辑（不新增函数版本，强调中间量复用与早停）**

```julia
for iter in 1:max_iters
    # 1) 单次扫描同时更新 F 与 dFdI，避免二次遍历
    @inbounds for e in 1:ne
        Ve = branch_voltage(coeffs[e], I_e[e])
        dV = branch_dVdI(coeffs[e], I_e[e])
        F[e] = Ve - V
        dFdI[e] = abs(dV) < 1e-12 ? (dV == 0.0 ? -coeffs[e].C5 : sign(dV) * 1e-12) : dV
    end

    # 2) 收敛早停：先判断 res_V 与 res_I，再决定是否进入 line_search
    # ... res_V/res_I 计算与判定 ...

    # 3) 仅在需要时调用 line_search；line_search 内仅更新活跃索引
    λ, V_trial = line_search(I_e, V, ΔI, ΔV, I_trial, ne)
    λ == 0.0 && break
end
```

#### 2.6.2 截止检测与写回：旧/新逻辑对照

**旧逻辑（字段分散、返回结构弱类型）**

```julia
# detect_cutoff_elements 返回 NamedTuple/Dict 风格
active_mask, n_cutoff, cutoff_info = detect_cutoff_elements(...)

# solve_branch_currents 内分散写回
variables["thermal2D cutoff_elements"] = Float64.(cutoff_info["cutoff_elements"])
variables["thermal2D cutoff_ocv"] = cutoff_info["cutoff_ocv"]
```

**新逻辑（`CutoffInfo` 强类型 + 单路径写回）**

```julia
ci = detect_cutoff_elements(coeffs, ne, V_MIN, V_MAX, I_total, phi_scale)
active_mask = copy(ci.active_mask)

# 合并 CZM 失效后统一写回
variables["thermal2D n_cutoff_elements"] = Float64(ci.n_cutoff)
variables["thermal2D cutoff_elements"] = Float64.(ci.cutoff_elements)
variables["thermal2D cutoff_ocv"] = ci.cutoff_ocv
variables["thermal2D margin_to_cutoff"] = ci.margin
```

#### 2.6.3 推导式/生成器微优化：旧/新逻辑对照

**旧逻辑（隐式循环，临时对象较多）**

```julia
V = isempty(active_idx) ? mean([coeffs[e].C1 for e in 1:ne]) :
    mean([branch_voltage(coeffs[e], I_e[e]) for e in active_idx])

sx = sum(w[e] * I_e[e] for e in active_idx)
```

**新逻辑（显式累加，减少临时分配）**

```julia
if isempty(active_idx)
    s = 0.0
    @inbounds for e in 1:ne
        s += coeffs[e].C1
    end
    V = s / ne
else
    s = 0.0
    @inbounds for e in active_idx
        s += branch_voltage(coeffs[e], I_e[e])
    end
    V = s / length(active_idx)
end

sx = 0.0
@inbounds for e in active_idx
    sx += w[e] * I_e[e]
end
```

#### 2.6.4 嵌套优化建议汇总表（09 专用）

| 目标位置 | 嵌套类型 | 旧逻辑问题 | 新逻辑方向（不新增函数） | 预期收益 |
|---|---|---|---|---|
| `newton_iteration` (`iter->e`) | 显式 | 每迭代重复分散计算中间量 | 单次扫描复用 `F/dFdI`，收敛早停 | 迭代耗时下降 |
| `newton_iteration -> line_search` | 隐式调用链 | 外迭代内再触发 `attempt->e` 全量扫描 | `λ==0` 快速中断，仅必要时进入线搜索 | 无效尝试减少 |
| `detect_cutoff_elements + 写回` | 显式+隐式 | 截止信息字段散、写回路径长 | `CutoffInfo` + 单路径写回 | 可维护性提升 |
| `mean/sum` 推导式与生成器 | 隐式 | 临时对象与分配抖动 | 显式 `for` 累加替代 | 小幅降分配 |

---

## 3. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| 总行数 | 525 | ~510 |
| 主函数行数 | 152 | 84（重命名后已精简） |
| 截止检测返回 | NamedTuple | `CutoffInfo` struct |
| 新增函数 | 0 | 0 |

### 3.1 循环嵌套优化后的补充预期

| 指标 | 旧 | 新（目标） |
|------|-----|-----|
| Newton 单步中间量重复计算 | 较多 | 减少（复用 `F/dFdI`） |
| line_search 无效尝试 | 偏多 | 减少（早停与快速中断） |
| 推导式/生成器临时分配 | 存在 | 降低（显式累加） |
| 数值行为 | 基线 | 保持一致（需回归验证） |
