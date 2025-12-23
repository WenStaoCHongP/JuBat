# 热应力计算修复 - 更改日志

**日期**: 2025-12-23  
**版本**: 修复热应力计算中的温度历史问题

## 问题描述

### 症状
- 热应力计算结果始终为 0
- `stress_thermal_max_hist` 全为 0 或接近 0

### 根本原因
`result["thermal2D T_nodes [K]"]` 只保存了最终时刻的温度场，没有时间历史。在 `testexample.jl` 计算应力时，所有时间步都使用相同的温度值，导致：
- 温度变化 `ΔT = T(t) - T₀ = 0`
- 热应变 `ε_thermal = α·ΔT = 0`
- 热应力 `σ_thermal = E·ε_thermal/(1-ν) = 0`

## 修改清单

### 1. src/Variables.jl (第 122 行)

**位置**: Phase B - 多SPMe模式变量初始化

**修改前**:
```julia
variables["thermal2D element thermal strain"] = zeros(Float64, ne, num)
# 未初始化温度历史
```

**修改后**:
```julia
variables["thermal2D element thermal strain"] = zeros(Float64, ne, num)
# 添加温度场历史（用于热应力计算）
variables["thermal2D temperature"] = zeros(Float64, nT, num)
```

**作用**: 为多SPMe模式添加温度场历史记录的预分配空间

---

### 2. src/Solve.jl (第 207-210 行)

**位置**: 主时间步循环，提取温度自由度后

**修改前**:
```julia
T_nodes = y_c[(end - nT + 1):end]
variables["T_nodes"] = T_nodes
T_nodes_carry = T_nodes
# 未保存到历史
```

**修改后**:
```julia
T_nodes = y_c[(end - nT + 1):end]
variables["T_nodes"] = T_nodes
T_nodes_carry = T_nodes
# 将温度场复制到 thermal2D temperature 以便记录历史（用于热应力计算）
if haskey(variables_hist, "thermal2D temperature")
    Tref = case.param_dim.scale.T_ref
    variables["thermal2D temperature"] = T_nodes .* Tref
end
```

**作用**: 每个时间步将当前温度场（转换为有量纲）保存到 `thermal2D temperature` 变量中，以便 `Variable_update!` 函数自动记录到历史矩阵

---

### 3. example/testexample.jl (第 206-230 行)

**位置**: 应力历史计算循环

**修改前**:
```julia
# 获取SOC和温度历史
if haskey(result, "thermal2D element soc_n") && haskey(result, "thermal2D T_nodes [K]")
    soc_n_hist = result["thermal2D element soc_n"]
    soc_p_hist = result["thermal2D element soc_p"]
    
    for step in 1:num_steps
        # 温度场（所有步使用相同的最终温度）
        T_nodes_K = result["thermal2D T_nodes [K]"]
        variables_step["T_nodes"] = T_nodes_K ./ T_ref
```

**修改后**:
```julia
# 获取SOC和温度历史
if haskey(result, "thermal2D element soc_n") && 
   (haskey(result, "thermal2D temperature [K]") || haskey(result, "thermal2D T_nodes [K]"))
    soc_n_hist = result["thermal2D element soc_n"]
    soc_p_hist = result["thermal2D element soc_p"]
    
    # 获取温度历史（优先使用完整历史）
    T_nodes_hist_K = if haskey(result, "thermal2D temperature [K]")
        result["thermal2D temperature [K]"]  # (nT × num_steps)
    else
        # 后备方案：使用最终温度
        T_final = result["thermal2D T_nodes [K]"]
        repeat(T_final, 1, num_steps)
    end
    
    for step in 1:num_steps
        # 温度场（使用当前时间步的温度）
        T_nodes_K = T_nodes_hist_K[:, step]
        variables_step["T_nodes"] = T_nodes_K ./ T_ref
```

**作用**: 
- 优先使用完整的温度历史 `thermal2D temperature [K]`
- 每个时间步使用对应的温度场 `T_nodes_hist_K[:, step]`
- 提供向后兼容的后备方案

---

### 4. example/spme_thermal2d_example.jl (第 215 行)

**修改前**:
```julia
variables = JuBat.thermal_stress(case, variables)
```

**修改后**:
```julia
variables = JuBat.thermal_diffusion_stress_2D(case, variables)
```

**作用**: 修正函数名，`thermal_stress` 不存在，应为 `thermal_diffusion_stress_2D`

---

### 5. tools/plot_thermal_stress.jl (第 32-34 行)

**修改前**:
```julia
variables = JuBat.thermal_stress(case, variables)
if !haskey(variables, "thermal2D element thermal stress")
    println("thermal2D element thermal stress not computed by thermal_stress")
```

**修改后**:
```julia
variables = JuBat.thermal_diffusion_stress_2D(case, variables)
if !haskey(variables, "thermal2D element thermal stress")
    println("thermal2D element thermal stress not computed by thermal_diffusion_stress_2D")
```

**作用**: 修正函数名和错误提示信息

---

## 技术细节

### 数据流
1. **初始化** (`Variables.jl`):
   ```
   variables_hist["thermal2D temperature"] = zeros(nT, num_steps)
   ```

2. **每个时间步** (`Solve.jl`):
   ```
   T_nodes = y_c[化学DOF之后的nT个]  // 无量纲
   variables["thermal2D temperature"] = T_nodes * Tref  // 转换为 K
   Variable_update!(variables_hist, variables, v)  // 自动保存到历史矩阵
   ```

3. **结果输出** (`Solve.jl`):
   ```
   result["thermal2D temperature [K]"] = variables_hist["thermal2D temperature"][:, 1:v]
   ```

4. **后处理** (`testexample.jl`):
   ```
   T_nodes_hist_K = result["thermal2D temperature [K]"]  // (nT × num_steps)
   for step in 1:num_steps
       T_step = T_nodes_hist_K[:, step]  // 使用对应时刻的温度
       计算热应力...
   ```

### 热应力计算原理
```
ΔT(x,y,t) = T(x,y,t) - T₀
ε_thermal(x,y,t) = α_eff · ΔT(x,y,t)
σ_thermal(x,y,t) = E_eff / (1-ν_eff) · ε_thermal(x,y,t)
```

现在每个时间步 `t` 都使用该时刻的实际温度场 `T(x,y,t)`，而不是固定的最终温度。

---

## 验证方法

### 1. 运行测试
```bash
cd /workspace/example
julia testexample.jl
```

### 2. 检查输出
查看控制台输出，应该看到：
```
热应力峰值范围: [非零值, 非零值] MPa
```

### 3. 查看图表
检查生成的应力历史图：
- `output/testexample_total_stress.png`
- 热应力曲线应该非零且随时间变化

### 4. 数据验证
```julia
result = JuBat.Solve(case)
T_hist = result["thermal2D temperature [K]"]  # 应该是 (nT × num_steps) 矩阵
println(size(T_hist))  # 应该是 (节点数, 时间步数)
println(T_hist[1,1], " vs ", T_hist[1,end])  # 不同时刻应该不同
```

---

## 向后兼容性

- ✅ 简化耦合模式（已有 `thermal2D temperature`）：无影响
- ✅ 标准SPMe模式（无热计算）：无影响
- ✅ 旧版结果文件：通过后备方案兼容
- ✅ 内存开销：仅在需要时分配

---

## 影响范围

### 直接影响
- 热应力计算现在能正确反映温度变化
- 应力时间历史曲线正确

### 间接影响
- 温度场数据可用于其他时间相关的后处理
- 为未来的瞬态热-力耦合分析奠定基础

---

## 已知限制

1. **内存使用**: 温度历史需要额外的内存 `nT × num_steps × 8 bytes`
   - 对于大规模网格和长时间模拟，可能需要几百MB
   - 建议：减少 `num_steps` 或使用稀疏采样

2. **简化耦合模式**: 已有温度历史记录，但其他模式可能需要类似的修改

3. **初始温度**: 假设 `T₀ = param.cell.T0`，如果实际初始温度不同，需要额外处理

---

## 后续工作建议

1. **性能优化**: 
   - 考虑稀疏时间采样（不是每步都保存）
   - 实现流式输出到文件

2. **功能增强**:
   - 保存初始温度场 `T₀(x,y)` 用于更准确的 `ΔT` 计算
   - 支持用户指定参考温度

3. **测试覆盖**:
   - 添加单元测试验证温度历史记录
   - 添加集成测试验证热应力计算

---

**修改人**: AI Assistant (Claude Sonnet 4.5)  
**审核状态**: 待审核  
**测试状态**: 待测试
