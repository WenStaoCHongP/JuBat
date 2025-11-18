# 阶段4实现说明：Solve主循环修改

**完成日期**: 2025-11-17  
**版本**: v1.0  
**状态**: ✅ 已完成并通过测试

---

## 📋 概述

本文档说明阶段4对 `Solve` 主循环和结果记录系统的修改，使其支持多SPMe并行架构的完整时间推进仿真。

---

## 🎯 设计目标

### 主要目标
- ✅ 支持多SPMe模式的初始化
- ✅ 保持时间推进逻辑兼容
- ✅ 增强结果记录（逐单元变量）
- ✅ 完全向后兼容单SPMe模式

### 非目标（保持不变）
- ❌ 时间推进算法（Crank-Nicolson等）
- ❌ 误差估计和步长控制
- ❌ 温度提取和热应力计算
- ❌ 边界条件处理

---

## 🔧 修改1：Solve函数初始化切换

### 修改位置
`src/Solve.jl` 第9-26行（原第9-10行）

### 原代码
```julia
# initialisation 
y0 = ModelInitialisation(case)
```

### 新代码
```julia
# initialisation（根据模式选择初始化函数）
multi_spme_enabled = (
    case.opt.model == "SPMe" &&
    hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme &&
    case.opt.thermalmodel == "distributed2D" &&
    haskey(case.mesh, "thermal2D")
)

if multi_spme_enabled
    y0 = ModelInitialisation_MultiSPMe(case)
    if hasproperty(case.opt, :debug_multi_spme) && case.opt.debug_multi_spme
        println("[Solve] 多SPMe模式：状态向量维度 = $(length(y0))")
        layout = case.multi_spme_layout
        println("  ne = $(layout["ne"]), n_chem = $(layout["n_chem"]), nT = $(layout["nT"])")
    end
else
    y0 = ModelInitialisation(case)
end
```

### 设计说明

**1. 模式判断（5个条件）**

| 条件 | 说明 | 原因 |
|-----|------|------|
| `model == "SPMe"` | 模型必须是SPMe | 多SPMe仅支持SPMe |
| `per_element_spme == true` | 标志启用 | 用户显式启用 |
| `thermalmodel == "distributed2D"` | 必须是分布式热 | 多SPMe需要热单元 |
| `haskey(mesh, "thermal2D")` | 网格已创建 | 确保ne定义 |

**注意**: 相比 `CallModel` 的判断，这里**不检查** `multi_spme_layout`，因为它在初始化函数中才创建。

**2. 调试输出（可选）**

启用条件：
```julia
case.opt.debug_multi_spme = true
```

输出示例：
```
[Solve] 多SPMe模式：状态向量维度 = 3060
  ne = 50, n_chem = 60, nT = 66
```

用途：
- 验证初始化成功
- 检查状态向量维度
- 调试维度不匹配问题

---

## 🔧 修改2：结果记录增强

### 修改位置
`src/Variables.jl` 第102-112行

### 新增代码
```julia
# Phase B: 多SPMe模式的逐单元变量历史记录（如果启用）
if hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme && 
   case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
    ne = size(case.mesh["thermal2D"].element, 1)
    variables["thermal2D element current"] = zeros(Float64, ne, num)
    variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
    variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
    variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
    variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
    variables["thermal2D element voltages"] = zeros(Float64, ne, num)
end
```

### 记录的变量

| 变量名 | 物理意义 | 单位 | 维度 |
|-------|---------|------|------|
| `thermal2D element current` | 逐单元电流 | A (SI) | ne × num |
| `thermal2D eta_n_e` | 逐单元负极过电位 | V (SI) | ne × num |
| `thermal2D eta_p_e` | 逐单元正极过电位 | V (SI) | ne × num |
| `thermal2D dUdT_n_e` | 逐单元负极dU/dT | V/K | ne × num |
| `thermal2D dUdT_p_e` | 逐单元正极dU/dT | V/K | ne × num |
| `thermal2D element voltages` | 逐单元电压 | V (SI) | ne × num |

### 自动填充机制

**无需修改** `Variable_update!` 函数（第117-138行），它会自动处理：

```julia
function Variable_update!(variables_hist, variables, v)
    for k in keys(variables_hist)
        if haskey(variables, k)
            val = variables[k]
            # 自动提取并填充
            variables_hist[k][:, v] = val
        end
    end
end
```

**工作原理**:
1. `CallModel_MultiSPMe` 计算并写入 `variables["thermal2D element current"]` 等
2. `Variable_update!` 检测到 `variables_hist` 中有这些键
3. 自动提取并填充到第 `v` 列

---

## ✅ 验证：时间推进逻辑兼容性

### 关键洞察

**状态向量结构对比**:

```julia
# 单SPMe
y = [yt_chem[1:60]; T_nodes[1:nT]]
# 总长度: 60 + nT

# 多SPMe (ne=50)
y = [yt_e[1][1:60]; yt_e[2][1:60]; ...; yt_e[50][1:60]; T_nodes[1:nT]]
# 总长度: 50×60 + nT = 3000 + nT
```

**共同点**: `T_nodes` 位于最后 `nT` 个位置！

### 温度提取代码（无需修改）

**位置**: `src/Solve.jl` 第122-130行

```julia
if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D" && haskey(case.mesh, "thermal2D")
    nT = case.mesh["thermal2D"].nlen
    n_tot = size(M_new, 1)
    if length(y_c) == n_tot
        T_nodes = y_c[(end - nT + 1):end]  # ← 提取最后nT个
        variables["T_nodes"] = T_nodes
        T_nodes_carry = T_nodes
        variables = thermal_stress(case, variables)
    end
end
```

**验证**:
- 单SPMe: `y_c[(60+nT - nT + 1):(60+nT)] = y_c[61:end]` ✓
- 多SPMe: `y_c[(3000+nT - nT + 1):(3000+nT)] = y_c[3001:end]` ✓

**结论**: 代码完全兼容，无需修改！

---

## 📊 完整工作流程

### 初始化阶段

```julia
# 1. 判断模式
multi_spme_enabled = (...)

# 2. 选择初始化
if multi_spme_enabled
    y0 = ModelInitialisation_MultiSPMe(case)  # → [yt_e[1]; ...; yt_e[ne]; T_nodes]
else
    y0 = ModelInitialisation(case)  # → [yt_chem; T_nodes]
end

# 3. 创建变量历史
variables_hist = StandardVariables(case, num)
# 如果多SPMe，自动添加逐单元变量空间
```

---

### 时间推进阶段

```julia
while t <= t_end
    # 1. 电化学步（CallModel自动切换）
    M_new, K_new, F_new, variables, y_phi = CallModel(case, y_old, t, jacobi="update")
    # 单SPMe → 调用 SPMe
    # 多SPMe → 调用 CallModel_MultiSPMe
    
    # 2. 时间积分
    Mt = M_new - theta * K_new * dt
    Kt = (1 - theta) * K_old * dt + M_new
    Ft = theta * F_new * dt + (1 - theta) * F_old * dt
    y_c = Mt \ (Kt * y_old[vc] + Ft)
    y_new = vcat(y_c, y_phi)
    
    # 3. 提取温度（通用逻辑，单/多SPMe均适用）
    if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
        nT = case.mesh["thermal2D"].nlen
        T_nodes = y_c[(end - nT + 1):end]
        variables["T_nodes"] = T_nodes
    end
    
    # 4. 记录结果
    if (记录条件)
        Variable_update!(variables_hist, variables, v)
        # 自动记录所有变量，包括逐单元变量（如果有）
        v = v + 1
    end
    
    # 5. 更新
    y_old = y_new
    K_old = K_new
    F_old = F_new
    t += dt
end
```

---

### 后处理阶段

```julia
result = PostProcessing(case, variables_hist, v)

# 结果包含（多SPMe模式）：
result["time"]                          # 时间历史
result["cell voltage"]                  # 全局电压
result["cell current"]                  # 全局电流
result["temperature"]                   # 平均温度
result["thermal2D element current"]     # 逐单元电流 (ne × v)
result["thermal2D eta_n_e"]             # 逐单元负极η (ne × v)
result["thermal2D eta_p_e"]             # 逐单元正极η (ne × v)
result["thermal2D dUdT_n_e"]            # 逐单元负极dU/dT (ne × v)
result["thermal2D dUdT_p_e"]            # 逐单元正极dU/dT (ne × v)
result["thermal2D element voltages"]    # 逐单元电压 (ne × v)
```

---

## 🧪 使用示例

### 基本使用

```julia
using JuBat

# 1. 创建多SPMe案例
param_dim = ChooseCell("LG M50")
opt = Option()
opt.model = "SPMe"
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true  # ← 启用多SPMe
opt.Current = x -> 10.0
opt.time = [0 100]

case = SetCase(param_dim, opt)

# 2. 创建thermal2D网格
case.mesh["thermal2D"] = Mesh2D(nodes, elements, ...)

# 3. 求解（自动使用多SPMe）
result = Solve(case)

# 4. 基本结果
t = result["time"]
V = result["cell voltage"]
I = result["cell current"]
```

---

### 逐单元分析

```julia
# 提取逐单元电流
I_e_hist = result["thermal2D element current"]  # ne × num

# 绘制所有单元的电流演化
using Plots
plot(t, I_e_hist', legend=false, 
     xlabel="Time (s)", ylabel="Element Current (A)",
     title="Per-Element Current Evolution")

# 分析电流分布（最终时刻）
I_e_final = I_e_hist[:, end]
histogram(I_e_final, 
          xlabel="Current (A)", ylabel="Count",
          title="Current Distribution at t=$(t[end])s")

# 找出电流最大的单元
e_max = argmax(I_e_final)
println("Maximum current in element $e_max: $(I_e_final[e_max]) A")

# 绘制该单元的详细演化
plot(t, I_e_hist[e_max, :], label="Current")
plot!(twinx(), t, result["thermal2D eta_n_e"][e_max, :], 
      label="η_n", color=:red)
```

---

### 异质性分析

```julia
# 计算逐单元变量的时间演化统计
ne = size(I_e_hist, 1)
nt = length(t)

mean_I = mean(I_e_hist, dims=1)[:]
std_I = std(I_e_hist, dims=1)[:]
cv_I = std_I ./ mean_I  # 变异系数

# 绘制异质性演化
plot(t, cv_I, xlabel="Time (s)", ylabel="CV of Current",
     title="Heterogeneity Evolution", label="CV(I_e)")

# 分析过电位异质性
eta_n_hist = result["thermal2D eta_n_e"]
eta_p_hist = result["thermal2D eta_p_e"]
eta_total = eta_p_hist .- eta_n_hist  # 总过电位

cv_eta = std(eta_total, dims=1)[:] ./ abs.(mean(eta_total, dims=1)[:])
plot!(t, cv_eta, label="CV(η_total)")
```

---

### 空间分布可视化

```julia
# 假设thermal2D网格是矩形
nx, ny = 10, 5  # 网格尺寸
I_e_final = result["thermal2D element current"][:, end]
I_e_matrix = reshape(I_e_final, nx, ny)

# 热图
heatmap(I_e_matrix', 
        xlabel="X", ylabel="Y", 
        title="Current Distribution",
        colorbar_title="I (A)")

# 3D表面图
surface(I_e_matrix',
        xlabel="X", ylabel="Y", zlabel="I (A)",
        title="Current Field")
```

---

## 🔍 调试技巧

### 1. 启用调试输出

```julia
case.opt.debug_multi_spme = true
result = Solve(case)
```

输出：
```
[Solve] 多SPMe模式：状态向量维度 = 3060
  ne = 50, n_chem = 60, nT = 66
```

---

### 2. 检查状态向量结构

```julia
# 在初始化后
layout = case.multi_spme_layout
ne = layout["ne"]
n_chem = layout["n_chem"]
nT = layout["nT"]

@assert length(y0) == ne * n_chem + nT "状态向量维度不匹配"

# 验证T_nodes位置
T_nodes_check = y0[(end - nT + 1):end]
@assert all(T_nodes_check .≈ case.param.cell.T0) "T_nodes初始化错误"
```

---

### 3. 验证逐单元变量记录

```julia
result = Solve(case)

# 检查变量存在性
required_vars = [
    "thermal2D element current",
    "thermal2D eta_n_e",
    "thermal2D eta_p_e"
]

for var in required_vars
    if haskey(result, var)
        println("✓ $var: $(size(result[var]))")
    else
        println("✗ 缺少 $var")
    end
end

# 检查维度
ne = size(case.mesh["thermal2D"].element, 1)
num_steps = length(result["time"])

for var in required_vars
    if haskey(result, var)
        expected_size = (ne, num_steps)
        actual_size = size(result[var])
        @assert expected_size == actual_size "$var 维度错误"
    end
end
```

---

### 4. 验证物理合理性

```julia
# 电压范围
V = result["cell voltage"]
@assert all(2.0 .< V .< 5.0) "电压超出合理范围"

# 温度范围
if haskey(result, "temperature")
    T = result["temperature"]
    @assert all(250 .< T .< 400) "温度超出合理范围"
end

# 电流守恒
I_total = result["cell current"][end]
I_e = result["thermal2D element current"][:, end]
areas = case.thermal2D_element_area_cache
w = areas ./ sum(areas)
I_sum = sum(w .* I_e)
@assert abs(I_total - I_sum) < 1e-6 "电流不守恒"

# 热源非负（如果放电）
if I_total > 0  # 放电
    q_e = result["heat_source_fields"][:, end]
    @assert all(q_e .>= 0) "放电时热源应为正"
end
```

---

## 📈 性能考虑

### 内存占用估算

```julia
# 状态向量
mem_y = (ne * n_chem + nT) * 8  # 字节

# 矩阵（稀疏，假设10%非零元）
nnz = (ne * n_chem + nT)^2 * 0.1
mem_M = nnz * 8 * 2  # M和K

# 变量历史（逐单元变量）
num_steps = 100
mem_hist = ne * num_steps * 8 * 6  # 6个逐单元变量

# 总计
mem_total = mem_y + mem_M + mem_hist

println("预估内存占用:")
println("  状态向量: $(mem_y / 1024) KB")
println("  矩阵: $(mem_M / 1024 / 1024) MB")
println("  历史: $(mem_hist / 1024 / 1024) MB")
println("  总计: $(mem_total / 1024 / 1024) MB")
```

**示例**（ne=50, nT=66, num_steps=100）:
- 状态向量: ~24 KB
- 矩阵: ~7.5 MB
- 历史: ~2.4 MB
- 总计: ~10 MB

---

### 计算时间估算

```julia
# 单步时间（经验公式）
t_per_step = 1e-6 * (ne * n_chem + nT)^2  # 秒

# 总时间
t_total = t_per_step * num_steps

println("预估计算时间:")
println("  单步: $(t_per_step * 1000) ms")
println("  总计: $(t_total) s")
```

**示例**（ne=50, nT=66, num_steps=100）:
- 单步: ~9 ms
- 总计: ~0.9 s

**实际**: 通常比估算高2-5倍（非线性求解、I/O等）

---

## 🔗 依赖关系

### 阶段依赖

```
阶段1: SPMe_element
  ↓
阶段2: ModelInitialisation_MultiSPMe
  ↓
阶段3: CallModel_MultiSPMe
  ↓
阶段4: Solve主循环 ← 当前
```

### 函数依赖

**Solve函数依赖**:
- `ModelInitialisation_MultiSPMe` (阶段2)
- `CallModel` → `CallModel_MultiSPMe` (阶段3)
- `StandardVariables` (已修改)
- `Variable_update!` (已有，无需修改)
- `PostProcessing` (已有，无需修改)

---

## 📚 相关文档

- [阶段1_SPMe_element_实现说明.md](./阶段1_SPMe_element_实现说明.md)
- [阶段2_多SPMe初始化_实现说明.md](./阶段2_多SPMe初始化_实现说明.md)
- [阶段3_CallModel_MultiSPMe_实现说明.md](./阶段3_CallModel_MultiSPMe_实现说明.md)
- [阶段4完成总结.md](./阶段4完成总结.md)
- [多SPMe并行架构修改计划.md](./多SPMe并行架构修改计划.md)

---

**完成日期**: 2025-11-17  
**版本**: v1.0  
**作者**: AI Coding Assistant
