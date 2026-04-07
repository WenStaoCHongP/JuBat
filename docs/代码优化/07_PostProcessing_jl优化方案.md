# PostProcessing.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/PostProcessing.jl`
> 状态: 修改 (48→311 行, +548%)
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
| 79-311 | 9 个 cycling 辅助函数 + 绘图 | 新增 |

### 2.2 distributed2D 热源输出 (行 49-64)

```julia
elseif case.opt.thermalmodel == "distributed2D"
    result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] * case.param.scale.q
    result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] * case.param.scale.q
    # ... 11 行热源 + 节点温度还原 ...
end
```

### 2.3 cycling 辅助函数 (行 79-311)

| 函数 | 行数 | 被谁调用 |
|------|------|---------|
| `phase_termination_symbol` | 8 | `postprocess_phase_result` |
| `state_concentration_variance` | 30 | `postprocess_phase_result` |
| `postprocess_phase_result` | 60 | `CycleSolver.jl` |
| `postprocess_cycle_result!` | 15 | `CycleSolver.jl` |
| `append_cycle_result!` | 15 | `CycleSolver.jl` |
| `update_soh_and_capacity!` | 10 | `CycleSolver.jl` |
| `print_cycle_summary` | 5 | `CycleSolver.jl` |
| `check_cycle_termination` | 20 | `CycleSolver.jl` |
| `print_cycling_summary` | 16 | `CycleSolver.jl` |
| `plot_cycling_results` | 40 | 外部 |

---

## 3. 优化方案

### 3.1 约束

- `PostProcessing()` 函数第 1-47 行（SPM/SPMe/P2D 分支）完全不动
- 仅修改第 48 行之后的 distributed2D 部分
- 不新增函数，distributed2D 热源输出保持内联

### 3.2 distributed2D 热源输出保持内联

distributed2D 的 16 行热源还原逻辑保留在 `PostProcessing` 主函数的 `elseif` 分支内，不提取子函数。

### 3.3 `state_concentration_variance` 简化

```julia
# 旧 (30 行, 含手动 layout 访问):
function state_concentration_variance(case::Case, y_state)
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
        # ... 手动 for 循环 ...
    end
    # ...
end

# 新 (~15 行):
function state_concentration_variance(case::Case, y_state)
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

### 3.4 cycling 辅助函数

所有 cycling 辅助函数保持原样，不做拆分。仅在函数中涉及 `multi_spme_layout` 的地方替换为 `case.layout` 字段访问。

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| `PostProcessing` 主函数行数 | 77 (含 distributed2D) | 77 (不变) |
| distributed2D 热源输出 | 16 行内联 | 16 行内联（不变） |
| `state_concentration_variance` | 30 行 | ~15 行 |
| cycling 辅助函数 | 不变 (~230 行) | 不变（仅替换 layout 访问） |
| main 分支 `PostProcessing` 行 1-47 | 不动 | 不动 |
| 新增函数 | 0 | 0 |
