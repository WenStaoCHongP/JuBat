# PostProcessing.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/PostProcessing.jl`
> 状态: 修改 (48→312 行, +548%)
> main 分支行数: 48

---

## 1. main 分支现状

### 1.1 `PostProcessing(case, variables, v)` — 48 行

单一函数，职责清晰：无量纲变量 × scale → 物理单位输出。

结构：`result = Dict()` → 按 SPM/SPMe/P2D 分支填充 → 返回。

**设计简洁，第 1-47 行完全不动。**

---

## 2. 当前分支变更分析

### 2.1 变更分区

| 行范围 | 内容 | 性质 |
|--------|------|------|
| 1-47 | `PostProcessing()` 主函数 SPM/SPMe/P2D 分支 | main 分支，不动 |
| 48-64 | distributed2D 热源输出的 `elseif` | 新增，可重构 |
| 65-77 | CZM 输出（被注释掉） | 新增但未启用 |
| 79-312 | 9 个 cycling 辅助函数 + 绘图 | 新增，可重构 |

### 2.2 distributed2D 热源输出 (行 49-64)

```julia
elseif case.opt.thermalmodel == "distributed2D"
    result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] * case.param.scale.q
    result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] * case.param.scale.q
    # ... 11 行热源 + 节点温度还原 ...
end
```

问题：与 main 分支的 SPM/SPMe/P2D 分支并列在同一个 `if/elseif/elseif` 链中，增加了认知负担。

### 2.3 cycling 辅助函数 (行 79-312)

| 函数 | 行数 | 被谁调用 |
|------|------|---------|
| `_phase_termination_symbol` | 8 | `_postprocess_phase_result` |
| `_state_concentration_variance` | 30 | `_postprocess_phase_result` |
| `_postprocess_phase_result` | 60 | `CycleSolver.jl` |
| `_postprocess_cycle_result!` | 15 | `CycleSolver.jl` |
| `_append_cycle_result!` | 15 | `CycleSolver.jl` |
| `_update_soh_and_capacity!` | 10 | `CycleSolver.jl` |
| `_print_cycle_summary` | 5 | `CycleSolver.jl` |
| `_check_cycle_termination` | 20 | `CycleSolver.jl` |
| `_print_cycling_summary` | 16 | `CycleSolver.jl` |
| `plot_cycling_results` | 40 | 外部 |

**问题**：这些函数全部是 cycling 基础设施，与基础 `PostProcessing` 功能无关，却挤在同一文件中。

---

## 3. 优化方案

### 3.1 约束

- `PostProcessing()` 函数第 1-47 行（SPM/SPMe/P2D 分支）完全不动
- 仅重构第 48 行之后的 distributed2D 和 cycling 部分

### 3.2 distributed2D 热源提取为子函数

```julia
# PostProcessing() 主函数中，行 47 之后：
# 旧:
#   elseif case.opt.thermalmodel == "distributed2D"
#       result["thermal2D Q_rxn_NE [W/m3]"] = ...
#       ... 16 行 ...

# 新:
    if case.opt.thermalmodel == "distributed2D"
        _dimensionalize_distributed2d!(result, case, variables, v)
    end
```

子函数：

```julia
"""
将 distributed2D 的无量纲变量还原为物理单位。
"""
function _dimensionalize_distributed2d!(result::Dict, case::Case,
                                         variables::Dict, v::Int)
    q_scale = case.param.scale.q
    T_ref = case.param_dim.scale.T_ref

    # 11 层热源 (内部键 → 输出键 的映射)
    heat_layer_map = [
        ("thermal2D q_rxn_ne",  "thermal2D Q_rxn_NE [W/m3]"),
        ("thermal2D q_rev_ne",  "thermal2D Q_rev_NE [W/m3]"),
        ("thermal2D q_ohm_s_ne","thermal2D Q_ohm_s_NE [W/m3]"),
        ("thermal2D q_ohm_e_ne","thermal2D Q_ohm_e_NE [W/m3]"),
        ("thermal2D q_sp",      "thermal2D Q_SP [W/m3]"),
        ("thermal2D q_rxn_pe",  "thermal2D Q_rxn_PE [W/m3]"),
        ("thermal2D q_rev_pe",  "thermal2D Q_rev_PE [W/m3]"),
        ("thermal2D q_ohm_s_pe","thermal2D Q_ohm_s_PE [W/m3]"),
        ("thermal2D q_ohm_e_pe","thermal2D Q_ohm_e_PE [W/m3]"),
        ("thermal2D q_pcc",     "thermal2D Q_PCC [W/m3]"),
        ("thermal2D q_ncc",     "thermal2D Q_NCC [W/m3]"),
    ]
    for (internal_key, output_key) in heat_layer_map
        result[output_key] = variables[internal_key][:, 1:v] * q_scale
    end

    # 节点温度
    result["thermal2D temperature at nodes [K]"] = variables["T_nodes"][:, 1:v] * T_ref
end
```

### 3.3 cycling 辅助函数方案：保留在同一文件但分组标注

不新增文件，但在文件内用注释分隔为两个逻辑段：

```julia
# ══════════════════════════════════════════════════════
# 第一段: 基础后处理 (PostProcessing 主函数 + distributed2D)
# ══════════════════════════════════════════════════════
function PostProcessing(case, variables, v)
    # ... main 分支不动 ...
    # ... distributed2D 调用子函数 ...
end

function _dimensionalize_distributed2d!(result, case, variables, v)
    # ...
end

# ══════════════════════════════════════════════════════
# 第二段: 循环仿真辅助 (仅供 CycleSolver.jl 使用)
# ══════════════════════════════════════════════════════
function _phase_termination_symbol(...) end
function _state_concentration_variance(...) end
function _postprocess_phase_result(...) end
function _postprocess_cycle_result!(...) end
function _append_cycle_result!(...) end
function _update_soh_and_capacity!(...) end
function _print_cycle_summary(...) end
function _check_cycle_termination(...) end
function _print_cycling_summary(...) end
function plot_cycling_results(...) end
```

### 3.4 `_state_concentration_variance` 简化

```julia
# 旧 (30 行, 含嵌套 if 和手动 for 循环):
function _state_concentration_variance(case::Case, y_state)
    if y_state === nothing
        return 0.0, 0.0
    end
    y = vec(y_state)
    multi_spme = case.opt.model == "SPMe" && case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D"
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    if multi_spme
        _ensure_multi_spme_layout!(case)
        ne = case.multi_spme_layout["ne"]
        n_chem = case.multi_spme_layout["n_chem"]
        cs_n_all = Float64[]
        cs_p_all = Float64[]
        for e in 1:ne
            offset = (e - 1) * n_chem
            cs_n_e = y[(offset + 1):(offset + Nrn)]
            cs_p_e = y[(offset + Nrn + 1):(offset + Nrn + Nrp)]
            push!(cs_n_all, mean(cs_n_e))
            push!(cs_p_all, mean(cs_p_e))
        end
        return var(cs_n_all), var(cs_p_all)
    end
    cs_n = y[1:Nrn]
    cs_p = y[(Nrn + 1):(Nrn + Nrp)]
    return var(cs_n), var(cs_p)
end

# 新 (~15 行):
function _state_concentration_variance(case::Case, y_state)
    y_state === nothing && return 0.0, 0.0
    y = vec(y_state)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen

    if case.layout !== nothing
        ne, n_chem = case.layout.ne, case.layout.n_chem
        cs_n = [mean(y[((e-1)*n_chem+1):((e-1)*n_chem+Nrn)]) for e in 1:ne]
        cs_p = [mean(y[((e-1)*n_chem+Nrn+1):((e-1)*n_chem+Nrn+Nrp)]) for e in 1:ne]
        return var(cs_n), var(cs_p)
    end
    return var(y[1:Nrn]), var(y[(Nrn+1):(Nrn+Nrp)])
end
```

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| `PostProcessing` 主函数行数 | 77 (含 distributed2D) | ~55 (distributed2D 提取) |
| distributed2D 热源输出 | 16 行内联 | ~15 行子函数 (可复用) |
| `_state_concentration_variance` | 30 行 | ~15 行 |
| cycling 辅助函数 | 不变 (~230 行) | 不变 |
| main 分支 `PostProcessing` 行 1-47 | 不动 | 不动 |
