# Solve.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/Solve.jl`
> 状态: 修改 (155→850 行, +448%)
> main 分支行数: 155

---

## 1. main 分支现状

### 1.1 函数清单

| 函数 | 行数 | 职责 |
|------|------|------|
| `Solve(case)` | ~100 | 纯时间步进器 |
| `CallModel(case, yt, t; jacobi)` | ~25 | 薄调度器 (SPM/SPMe/P2D/sP2D) |
| `RecordMatrix!(case, M, K)` | 10 | 稀疏矩阵记录 |
| `ErrorEstimation(case, y_old, y_new, coeff)` | 20 | dt 自适应 |

### 1.2 核心设计模式

```
Solve: 纯时间步进器
  ├── 初始化 y0, dt, theta
  ├── while t <= t_end
  │     CallModel → (M, K, F, variables, y_phi)
  │     时间离散: y_new = Mt \ (Kt*y_old + Ft)
  │     ErrorEstimation → dt自适应
  │     Variable_update!
  │     电压截止 → break
  └── PostProcessing
```

**main 分支 Solve 约 100 行，架构清晰。**

---

## 2. 当前分支膨胀根因

| 新增功能 | 塞入位置 | 行数 | 应在位置 |
|----------|----------|------|----------|
| 文件日志重定向 | Solve 开头 | ~35 | 辅助函数 `setup_logging!` |
| MultiSPMe 初始化分支 | Solve 中部 | ~40 | Initialisation.jl |
| distributed2D 热场初始化 | Solve 中部 | ~40 | CallModel_MultiSPMe |
| 单元截止电压检测 | Solve 循环内 | ~60 | CallModel 返回终止信号 |
| 物理单位后处理 | Solve 尾部 | ~60 | PostProcessing |
| CallModel_MultiSPMe | Solve 下半部 | ~180 | CallModel.jl |
| CallModel 扩展 | Solve 下半部 | ~110 | CallModel.jl |

---

## 3. 优化方案

### 3.1 拆分策略

```
Solve.jl (851 → ~200 行)           拆出:
──────────────────────────────     ──────────────────
Solve()           ~200行           (保留了纯步进器精神)
RecordMatrix!      ~10行
ErrorEstimation    ~25行
setup_logging!     ~15行  ← 提取为独立辅助函数

                                   CallModel.jl (新, ~200 行):
                                   ──────────────────────────
                                   call_model()           ~30行
                                   call_model_multi_spme()~120行
                                   detect_cutoff()        ~30行

                                   PostProcessing.jl (已有):
                                   ──────────────────────────
                                   物理单位后处理部分从 Solve 尾部迁入
```

### 3.2 Solve() 重构后结构

```julia
function Solve(case::Case; initial_state=nothing, return_final_state=false)
    # ═══ 1. 初始化 (~40 行) ═══
    if case.opt.model == "thermal"
        return _solve_pure_thermal(case; initial_state, return_final_state)
    end

    y0 = isempty(case.opt.y0) ?
        (is_multi_spme(case) ?
            ModelInitialisation_MultiSPMe(case) :
            ModelInitialisation(case)) :
        case.opt.y0

    if initial_state !== nothing
        y0 = _restore_state(y0, initial_state, case)
    end

    dt, theta, t, t_end = _init_time_params(case)
    v, v_max = 1, case.opt.max_steps
    variables_hist = StandardVariables(case, v_max)

    if case.layout !== nothing && case.opt.thermalmodel == "distributed2D"
        T0 = fill(case.param.cell.T0, case.layout.nT)
        variables_hist["T_nodes"][:, 1] = T0
    end

    # ═══ 2. 时间步进循环 (~80 行) ═══
    termination_reason = "time_limit"
    yt = y0
    while t[1] < t_end[1] && v < v_max
        M, K, F, variables, y_phi = call_model(case, yt, t[1]; jacobi="left")

        # 时间离散
        Mt = M - theta * dt[1] * K
        Ft = M * yt + dt[1] * F
        y_new = Mt \ Ft
        y_new[isnan.(y_new)] .= 0.0  # NaN 安全

        # 提取温度 (distributed2D)
        if case.layout !== nothing && case.opt.thermalmodel == "distributed2D"
            T_nodes = y_new[case.layout.thermal_range]
            variables["T_nodes"] = T_nodes
        end

        # dt 自适应
        dt_new = ErrorEstimation(case, yt, y_new, theta)
        dt[1] = clamp(dt_new, case.opt.dt[1], case.opt.dt[2])

        # 电压截止检测
        V = variables["cell voltage"]
        if _check_voltage_cutoff(V, case.opt)
            termination_reason = "voltage_cutoff"
            break
        end

        Variable_update!(variables_hist, variables, v)
        yt = y_new
        t[1] += dt[1]
        v += 1
    end

    # ═══ 3. 后处理 (~30 行) ═══
    result = PostProcessing(case, variables_hist, v)
    result["termination_reason"] = termination_reason
    result["steps"] = v

    if return_final_state
        result["final_state"] = _pack_final_state(yt, variables, case)
    end

    return result
end
```

### 3.3 辅助函数提取

```julia
"""纯热模式入口"""
function _solve_pure_thermal(case; initial_state=nothing, return_final_state=false)
    # 从当前 Solve.jl 行 1-95 提取的纯热路径
    # ~80 行
end

"""判断是否为 multi-SPMe 模式"""
function is_multi_spme(case::Case)
    return case.opt.model == "SPMe" &&
           case.opt.per_element_spme &&
           case.opt.thermalmodel == "distributed2D"
end

"""从 initial_state Dict 恢复状态向量"""
function _restore_state(y0, initial_state, case)
    if haskey(initial_state, "y")
        y = initial_state["y"]
        if is_multi_spme(case) && case.layout === nothing
            ne = size(case.mesh["thermal2D"].element, 1)
            nT = case.mesh["thermal2D"].nlen
            case.layout = MultiSPMeLayout(ne, length(y) - nT, nT)
        end
        return vec(y)
    end
    return y0
end

"""电压截止检测"""
function _check_voltage_cutoff(V, opt)
    if opt.V_max !== nothing && V > opt.V_max
        return true
    end
    if opt.V_min !== nothing && V < opt.V_min
        return true
    end
    return false
end

"""打包最终状态"""
function _pack_final_state(yt, variables, case)
    state = Dict{String,Any}("y" => copy(yt))
    for k in ("cell voltage", "temperature", "T_nodes")
        if haskey(variables, k)
            state[k] = copy(variables[k])
        end
    end
    return state
end

"""初始化时间参数"""
function _init_time_params(case)
    dt = [case.opt.dt[1]]
    theta = case.opt.solveType == "Crank-Nicolson" ? 0.5 : 1.0
    t = [case.opt.time[1]]
    t_end = [case.opt.time[2]]
    return dt, theta, t, t_end
end
```

---

## 4. CallModel.jl (新文件)

### 4.1 `call_model` 统一调度

```julia
# CallModel.jl

"""
统一模型调度入口。给定 (case, yt, t) 返回 (M, K, F, variables, y_phi)。
"""
function call_model(case::Case, yt, t::Float64; jacobi::String="left")
    if is_multi_spme(case)
        return call_model_multi_spme(case, yt, t; jacobi)
    end

    # ═══ 标准路径: SPM / SPMe / P2D (不动 main 分支逻辑) ═══
    if case.opt.model == "SPM"
        M, K, F, variables = SPM(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "SPMe"
        if case.opt.thermalmodel == "distributed2D"
            return call_model_spme_distributed2d(case, yt, t; jacobi)
        end
        M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
        y_phi = Float64[]
    elseif case.opt.model == "P2D"
        M, K, F, variables, y_phi = P2D(case, yt, t, jacobi=jacobi)
    else
        error("Model $(case.opt.model) not implemented")
    end

    # 可选 lumped 热拼接
    if case.opt.thermalmodel == "lumped"
        MT, FT = ThermalLumped(case, variables)
        M = blockdiag(M, sparse(MT))
        K = blockdiag(K, sparse(zeros(1,1)))
        F = [F; FT]
    end

    return M, K, F, variables, y_phi
end
```

### 4.2 `call_model_multi_spme` (~120 行)

```julia
"""
Multi-SPMe 路径: 逐单元 SPMe + 分流 + 热源 + distributed2D FEM。
"""
function call_model_multi_spme(case::Case, yt, t::Float64; jacobi::String="left")
    layout = case.layout
    layout === nothing && error("call_model_multi_spme: layout is nothing")

    ne = layout.ne
    n_chem = layout.n_chem
    y_vec = vec(yt)

    # 1) 提取元素面积和均温
    A_elem, Te_prev = _compute_element_geometry(y_vec, case)

    # 2) 分流求解
    variables = Dict{String,Any}()
    variables, I_e, V_common = solve_branch_currents_newton(
        case, variables, y_vec, t,
        case.param.cell.I1C * case.opt.Current(t), A_elem, Te_prev, ones(ne)
    )

    # 3) 逐单元 SPMe (并行)
    M_chem = spzeros(ne * n_chem, ne * n_chem)
    K_chem = spzeros(ne * n_chem, ne * n_chem)
    F_chem = zeros(ne * n_chem)

    Threads.@threads for e in 1:ne
        yt_e = extract_element_state(y_vec, e, layout)
        Me, Ke, Fe, vars_e = SPMe_element(case, yt_e, t, e; I_e=I_e[e], T_e=Te_prev[e])
        offset = (e - 1) * n_chem
        M_chem[(offset+1):(offset+n_chem), (offset+1):(offset+n_chem)] = Me
        K_chem[(offset+1):(offset+n_chem), (offset+1):(offset+n_chem)] = Ke
        F_chem[(offset+1):(offset+n_chem)] = Fe
        _merge_element_variables!(variables, vars_e, e)
    end

    # 4) 热源计算
    T_nodes = y_vec[layout.thermal_range]
    q_fields = compute_heat_sources(case, variables, A_elem, I_e, Te_prev, T_nodes)

    # 5) 热学 FEM 装配 + BC
    MT, KT, FT = ThermalDistributed2D(case, T_nodes)
    ThermalDistributed2D_BC(case, MT, KT, FT, T_nodes)

    # 6) 填充热源到 F 向量
    FT .+= q_fields  # 简化，实际需映射到节点

    # 7) 全局装配
    M = blockdiag(M_chem, MT)
    K = blockdiag(K_chem, KT)
    F = [F_chem; FT]

    return M, K, F, variables, Float64[]
end
```

### 4.3 辅助函数

```julia
"""计算单元面积和均温"""
function _compute_element_geometry(y_vec, case)
    ne = case.layout.ne
    nT = case.layout.nT

    # 面积
    if haskey(case.mesh, "thermal2D")
        elements = case.mesh["thermal2D"].element
        A_elem = [compute_element_area(elements[e, :], case.mesh["thermal2D"]) for e in 1:ne]
    else
        A_elem = ones(ne)
    end

    # 均温
    T_nodes = y_vec[case.layout.thermal_range]
    Te_prev = element_nodal_mean(T_nodes, case.mesh["thermal2D"])

    return A_elem, Te_prev
end

"""合并单元变量到全局 variables Dict"""
function _merge_element_variables!(variables, vars_e, e::Int)
    for k in keys(vars_e)
        key = "thermal2D element $(k)"
        if !haskey(variables, key)
            variables[key] = Float64[]  # 首次创建占位
        end
    end
end
```

---

## 5. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| Solve.jl 行数 | 851 | ~200 |
| CallModel.jl 行数 | 0 (在 Solve 内) | ~200 |
| `Solve()` 主函数行数 | 502 | ~200 |
| `call_model()` | 110 (在 Solve 内) | ~30 |
| `call_model_multi_spme()` | 180 (在 Solve 内) | ~120 |
| main 分支 SPM/SPMe 标准路径 | N/A | 零改动 |
| `RecordMatrix!` | 不动 | 不动 |
| `ErrorEstimation` | 不动 | 不动 |
