# Solve.jl 调试代码记录

> 文件: `src/Solve.jl`
> 状态: M (修改文件)
> 移除行数: ~90 行

---

## 1. 文件日志重定向系统 (原 lines 2-35, 489-500)

**功能**: 将 stdout/stderr 重定向到 `output/debug_<timestamp>.log` 文件，通过 `opt.debug_to_file` 控制。

### 1.1 日志初始化 (原 lines 2-35)

```julia
    # 将调试输出统一写入 output 目录下的日志文件（可通过 opt.debug_to_file 关闭）
    enable_file_log = true
    if hasproperty(case.opt, :debug_to_file)
        try
            enable_file_log = Bool(getfield(case.opt, :debug_to_file))
        catch
            enable_file_log = true
        end
    end
    log_io = nothing
    log_file = ""
    old_out = stdout
    old_err = stderr
    if enable_file_log
        log_dir = normpath(joinpath(@__DIR__, "..", "output"))
        try
            isdir(log_dir) || mkpath(log_dir)
            # 使用 epoch 秒作为时间戳，避免依赖 Dates 标准库的 import
            timestamp = string(round(Int, time()))
            log_file = joinpath(log_dir, "debug_$(timestamp).log")
            log_io = open(log_file, "w")
            println(log_io, "===== JuBat Debug Log $(timestamp) =====")
            flush(log_io)
            redirect_stdout(log_io)
            redirect_stderr(log_io)
        catch err
            # 如果日志初始化失败，则回退为标准控制台输出
            try
                log_io !== nothing && close(log_io)
            catch
            end
            log_io = nothing
        end
    end
```

### 1.2 try/finally 包装 (原 line 36-37, 489-500)

**原结构**:
```julia
    result = nothing
    try
        # ... 502行主逻辑 ...
    finally
        if log_io !== nothing
            # 恢复标准输出/错误到控制台
            redirect_stdout(old_out)
            redirect_stderr(old_err)
            try
                close(log_io)
            catch
            end
            println("Debug log saved to $(log_file)")
        end
    end
    return result
```

**恢复**: 在 Solve 函数开头恢复以上代码块，将主逻辑包裹在 `try/finally` 中。

---

## 2. 状态向量 NaN 检查 (原 lines 191-195, 245-248)

### 2.1 初始状态 NaN 检查 (原 lines 191-195)

**功能**: 检查初始状态向量是否包含 NaN/Inf，打印警告。

```julia
    # DEBUG: 检查初始状态向量（只在有问题时打印）
    nan_in_y0 = sum(.!isfinite.(y0))
    if nan_in_y0 > 0
        @warn "初始状态向量包含 $(nan_in_y0) 个 NaN/Inf，长度 $(length(y0))"
    end
```

### 2.2 初始求解 NaN 检查 (原 lines 245-248)

**功能**: 检查初始求解步骤结果是否异常。

```julia
    # 检查初始求解步骤
    nan_in_yold = sum(.!isfinite.(y_old))
    if nan_in_yold > 0 || (multi_spme_enabled && maximum(abs.(y_old)) > 100.0)
        @warn "初始求解步骤异常: NaN=$(nan_in_yold), 范围=[$(minimum(y_old)), $(maximum(y_old))]"
    end
```

**恢复**: 粘贴回对应位置。

---

## 3. 元素温度追踪 (原 lines 226-238, 309-317, 467-471)

**功能**: 追踪指定元素的温度变化历史，用于调试可视化。通过环境变量 `JUBAT_TRACK_ELEM` 控制追踪哪个元素。

### 3.1 追踪初始化 (原 lines 226-238)

```julia
    # 跟踪元素温度（用于调试/作图）
    track_elem_index = 0
    T_elem_hist, time_hist = Float64[], Float64[]
    if case.opt.thermal_enabled
        ne_track = size(case.mesh["thermal2D"].element, 1)
        if ne_track > 0
            idx_env = try parse(Int, get(ENV, "JUBAT_TRACK_ELEM", "")) catch; nothing end
            track_elem_index = Int(clamp(idx_env !== nothing ? idx_env : round(ne_track/2), 1, ne_track))
            nodes_e0 = case.mesh["thermal2D"].element[track_elem_index, :]
            Te0 = length(T_nodes_carry) == case.mesh["thermal2D"].nlen ? sum(T_nodes_carry[nodes_e0]) / length(nodes_e0) : case.param.cell.T0
            push!(T_elem_hist, Te0)
            push!(time_hist, t0 * case.param.scale.t0)
        end
    end
```

### 3.2 循环内温度追踪 (原 lines 308-317)

```julia
                # 同步跟踪元素温度（在提交此时间步后、时间推进前记录）
                if track_elem_index > 0 && isa(variables["T_nodes"], Array{Float64})
                    nodes_e = case.mesh["thermal2D"].element[track_elem_index, :]
                    Tn_now = variables["T_nodes"]
                    if length(Tn_now) == case.mesh["thermal2D"].nlen
                        Te_now = sum(Tn_now[nodes_e]) / length(nodes_e)
                        push!(T_elem_hist, Te_now)
                        push!(time_hist, t * case.param.scale.t0)
                    end
                end
```

### 3.3 结果输出 (原 lines 467-471)

```julia
            if !isempty(T_elem_hist) && length(T_elem_hist) == length(time_hist)
                result["thermal2D tracked element index"] = track_elem_index
                result["thermal2D tracked element time [s]"] = time_hist
                result["thermal2D tracked element T [K]"] = T_elem_hist .* case.param_dim.scale.T_ref
            end
```

**恢复**: 按原位置顺序粘贴。需要 `track_elem_index`, `T_elem_hist`, `time_hist` 三个变量。

---

## 4. 截止电压 verbose 输出 (原 lines 368-373, 383, 391, 395)

**功能**: 首个截止单元和整体电压截止时打印详细信息。保留追踪变量和 break 逻辑，仅删除 println。

### 4.1 首个截止单元详细输出 (原 lines 368-373)

```julia
                # 打印首个截止单元信息
                println("\n[Solve] ★ 首个单元达到截止电压:")
                println("  时间: $(round(first_cutoff_time, digits=2)) s")
                println("  单元: $first_cutoff_element")
                println("  OCV: $(round(first_cutoff_ocv, digits=4)) V")
                println("  截止电压: $v_l V")
                println("  当前整体电压: $(round(V_cell, digits=4)) V")
```

### 4.2 所有单元截止输出 (原 line 383)

```julia
                println("\n[Solve] ★★ 所有单元 ($ne_total) 都已达到截止电压，终止仿真")
```

### 4.3 整体电压截止输出 (原 lines 391, 395)

```julia
            println("\n[Solve] ★ 整体电压 $(round(V_cell, digits=4)) V 低于截止电压 $v_l V")
            println("\n[Solve] ★ 整体电压 $(round(V_cell, digits=4)) V 高于截止电压 $v_h V")
```

**恢复**: 粘贴回对应位置。注意保留 `first_cutoff_detected` 等追踪变量（这些是功能逻辑，非调试）。

---

## 5. CallModel_MultiSPMe NaN 检查 (原 lines 541-548, 559-568, 590-595, 683-688)

### 5.1 状态向量 NaN 检查 (原 lines 541-548)

```julia
    # 检查输入状态向量
    nan_count = sum(.!isfinite.(yt))
    if nan_count > 0
        chem_range = 1:(ne * n_chem)
        thermal_range = (ne * n_chem + 1):(ne * n_chem + nT)
        nan_chem = sum(.!isfinite.(yt[chem_range]))
        nan_thermal = sum(.!isfinite.(yt[thermal_range]))
        @warn "CallModel_MultiSPMe 收到 NaN 状态向量" t=t total_nan=nan_count chem_nan=nan_chem thermal_nan=nan_thermal
    end
```

### 5.2 温度场异常检查 (原 lines 559-568)

```julia
    # 温度场检查
    if length(T_nodes) > 0
        nan_count_here = sum(.!isfinite.(T_nodes))
        abnormal_count = sum(abs.(T_nodes) .> 10.0)
        large_deviation_count = sum(abs.(T_nodes .- 1.0) .> 5.0)

        if nan_count_here > 0 || abnormal_count > 0 || large_deviation_count > 0
            T_min, T_max = extrema(T_nodes)
            @warn "温度场异常" t=t range=(T_min,T_max) nan=nan_count_here abnormal=abnormal_count large_dev=large_deviation_count
        end
    end
```

### 5.3 元素温度 NaN 检查 (原 lines 590-595)

```julia
      # 检查温度场
    nan_count_nodes = sum(.!isfinite.(T_nodes))
    nan_count_elem = sum(.!isfinite.(Te_prev))
    if nan_count_nodes > 0 || nan_count_elem > 0
        @warn "温度场包含 NaN/Inf" T_nodes_nan=nan_count_nodes Te_prev_nan=nan_count_elem
    end
```

### 5.4 热源量级检查 (原 lines 683-688)

```julia
    # 检查参数（只在初始调用时）
    if t < 1e-6
        (t_ratio < 0.001 || t_ratio > 1000.0) && @warn "时间尺度比异常" t_ratio=t_ratio
        q_elem = variables["heat_source_fields"]
        q_max = maximum(abs.(q_elem))
        q_max > 100.0 && @warn "无量纲热源过大" q_max=q_max q_range=extrema(q_elem)
    end
```

**恢复**: 各块粘贴回原位置。

---

## 6. 恢复依赖

恢复调试代码时需同步恢复:
1. 日志系统需在 Solve 函数开头恢复 `enable_file_log` 等变量
2. 温度追踪需恢复 `track_elem_index`, `T_elem_hist`, `time_hist` 变量
3. 截止电压 println 需在对应位置恢复
4. CallModel_MultiSPMe NaN 检查需在各位置恢复

---

## 7. 2026-04-07 新增性能计时调试代码（当前有效）

> 说明: 本节不是“移除记录”，而是当前已加入源码的性能计时调试代码，用于定位四类核心流程耗时占比。

### 7.1 Solve 内累计器（约 lines 150-170）

在 `Solve` 开头新增累计字典与累加函数:

```julia
    timing_totals = Dict{String,Float64}(
        "spme" => 0.0,
        "branch" => 0.0,
        "thermal" => 0.0,
        "czm" => 0.0,
    )
    timing_call_count = 0

    function accumulate_callmodel_timing!(totals::Dict{String,Float64}, vars::Dict{String, Union{Array{Float64},Float64}})
        totals["spme"] += get(vars, "timing spme solve [s]", 0.0)
        totals["branch"] += get(vars, "timing branch solver [s]", 0.0)
        totals["thermal"] += get(vars, "timing thermal distributed [s]", 0.0)
        totals["czm"] += get(vars, "timing czm model [s]", 0.0)
        return nothing
    end
```

并在每次 `CallModel(...)` 后执行累加。

### 7.2 Solve 结果写回（约 lines 330-355）

新增以下 `result` 字段，供后处理和脚本打印:

- `timing SPMe solve total [s]`
- `timing branch solver total [s]`
- `timing thermal distributed total [s]`
- `timing CZM model total [s]`
- `timing SPMe solve avg [ms]`
- `timing branch solver avg [ms]`
- `timing thermal distributed avg [ms]`
- `timing CZM model avg [ms]`
- `timing SPMe solve ratio [%]`
- `timing branch solver ratio [%]`
- `timing thermal distributed ratio [%]`
- `timing CZM model ratio [%]`
- `timing CallModel calls`

### 7.3 控制台摘要输出（约 lines 358-370）

当 `case.opt.debug_coupling == true` 时，打印:

```text
[Solve-Timing] CallModel 阶段累计耗时（用于优化定位）
```

并输出四类模块的总耗时、占比和平均耗时（ms/call）。

### 7.4 示例脚本联动

`example/testexample.jl` 已在结果提取阶段打印同一组耗时分解，便于直接比较“SPMe / 分流 / 热分布式 / CZM”四类流程。
