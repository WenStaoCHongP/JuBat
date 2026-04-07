# CallModel.jl 调试代码记录（性能计时）

> 文件: src/CallModel.jl
> 状态: M (修改文件)
> 日期: 2026-04-07
> 目的: 统计多 SPMe 耦合路径中四个关键模块的单次调用耗时

---

## 1. 埋点范围

位置: CallModel_MultiSPMe(case, yt, t; jacobi)

新增四类计时:

1. 分流求解器
2. SPMe 并行求解
3. 热分布式模型（热源+热矩阵+边界）
4. CZM 模型（仅启用 CZM 时）

---

## 2. 关键代码片段

### 2.1 分流求解器计时

```julia
    t_branch_ns = time_ns()
    variables, I_e, Vc = solve_branch_currents(case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev; deactivated_elements=deactivated_elements)
    t_branch_s = (time_ns() - t_branch_ns) * 1e-9
```

### 2.2 SPMe 并行段计时

```julia
    t_spme_ns = time_ns()
    Threads.@threads for e in 1:ne
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e;I_e = I_e[e],T_e = Te_prev[e],jacobi = jacobi)
        ...
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
```

### 2.3 热模型和 CZM 计时

```julia
    t_thermal_ns = time_ns()
    t_czm_model_s = 0.0
    if case.opt.czm_enabled == true
        t_czm_ns = time_ns()
        variables = compute_heat_sources_with_czm(...)
        t_czm_model_s = (time_ns() - t_czm_ns) * 1e-9
    else
        variables = compute_heat_sources(...)
    end

    MT, KT, FT = ThermalDistributed2D(case, variables)
    KT, FT = ThermalDistributed2D_BC(KT, FT, case, t)
    t_thermal_s = (time_ns() - t_thermal_ns) * 1e-9
```

### 2.4 写入 variables（供 Solve 聚合）

```julia
    variables["timing branch solver [s]"] = t_branch_s
    variables["timing spme solve [s]"] = t_spme_s
    variables["timing thermal distributed [s]"] = t_thermal_s
    variables["timing czm model [s]"] = t_czm_model_s
```

---

## 3. 使用说明

1. CallModel 只负责“单次调用”计时，不负责跨步汇总。
2. 跨步汇总在 src/Solve.jl 中完成，最终写入 result。
3. 若 case.opt.czm_enabled=false，则 timing czm model [s] 恒为 0。

---

## 4. 注意事项

1. 计时单位统一为秒 (s)，使用 time_ns() 转换。
2. SPMe 计时包裹了 Threads.@threads 区段，反映并行段墙钟时间。
3. Thermal 计时包含热源更新与热矩阵装配+边界条件，属于“热分布式总开销”。
4. 当前埋点开销很小，可用于优化阶段对比，不建议长期用于超高频细粒度 profiling。
