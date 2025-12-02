# 简化耦合模型 - 快速开始指南

## 🚀 5分钟快速开始

### 1. 检查环境
```bash
cd /workspace
julia --version  # 确保Julia可用
```

### 2. 运行测试
```bash
julia example/testexample_simple_coupling.jl
```

### 3. 查看结果
生成的图像文件：
- `simple_coupling_voltage.png` - 电压曲线
- `simple_coupling_temperature.png` - 温度演化
- `simple_coupling_heat_source.png` - 热源演化
- `simple_coupling_temperature_field.png` - 温度场分布
- `simple_coupling_T_vs_radius.png` - 温度-半径关系

---

## 📝 核心概念

### 什么是简化耦合模式？
- **单个SPMe**：只用总电流求解一次电化学
- **2D热传导**：保留完整的温度空间分布
- **均匀热源**：假设热源均匀分布
- **平均温度反馈**：用平均温度影响电化学

### 与多SPMe模式对比

| 特性 | 多SPMe | 简化耦合 | 差异 |
|-----|--------|---------|------|
| **电流分布** | ✅ 计算 | ❌ 假设均匀 | 失去空间异质性 |
| **温度分布** | ✅ 保留 | ✅ 保留 | 相同 |
| **内存需求** | 30 GB | 3 GB | **节省90%** |
| **计算时间** | 60分钟 | 6分钟 | **快10倍** |
| **适用场景** | 热失控 | 常规放电 | - |

---

## ⚙️ 参数配置

### 启用简化耦合
```julia
opt.simple_thermal_coupling = true  # 关键：启用简化耦合
opt.per_element_spme = false        # 关闭多SPMe模式
opt.thermalmodel = "distributed2D"  # 2D热传导
opt.thermal_enabled = true
```

### 调试模式
```julia
opt.debug_simple_coupling = true  # 打印详细信息
```

### 网格设置
```julia
# 因为内存需求低，可以用更密的网格
nθ = 80   # 周向单元数（多SPMe通常只能用40-60）
```

### 时间设置
```julia
opt.time = [0.0, 3600]     # 仿真时间 (s)
opt.dt = [0.5, 10.0]       # 时间步长 [min, max] (s)
opt.dtType = "auto"        # 自动步长
```

---

## 📊 结果提取

### 电化学变量
```julia
result["cell voltage"]         # 电压 (V)
result["average temperature"]  # 平均温度 (K)
result["total heat source"]    # 总热源 (W)
```

### 温度场
```julia
# 所有节点的温度历史 (nT × num_steps)
T_nodes = result["thermal2D temperature"]

# 最后时刻的温度场
T_final = T_nodes[:, end]

# 温度统计
T_avg = result["average temperature"][end]
T_max = maximum(T_final)
T_min = minimum(T_final)
```

---

## ⚠️ 常见问题

### Q1: 运行时报错 "UndefVarError: ThermalDistributed2D"
**A**: 确保`src/JuBat.jl`已经加载`ThermalDistributed.jl`
```julia
# 检查src/JuBat.jl中是否有：
include("ThermalDistributed.jl")
```

### Q2: 温度异常（NaN或极端值）
**A**: 检查初始化和参数设置
```julia
# 启用调试
opt.debug_simple_coupling = true

# 检查初始温度
println("T0 = ", case.param.cell.T0, " K")
```

### Q3: 内存仍然不够
**A**: 进一步减少网格或增大时间步
```julia
nθ = 40               # 减少周向单元数
opt.dt = [1.0, 20.0]  # 增大时间步长
```

### Q4: 计算太慢
**A**: 增大时间步或减少仿真时间
```julia
opt.time = [0.0, 1800]  # 减少仿真时间
opt.dt = [1.0, 10.0]    # 增大时间步
```

---

## 🔬 验证检查

### 物理合理性
```julia
# 温度范围应在 280-350K
T_avg[end] > 280 && T_avg[end] < 350

# 温升应为正
T_avg[end] > T_avg[1]

# 热源应为正（放电时）
mean(Q_total) > 0
```

### 能量守恒（粗略）
```julia
# 总放电能量
E_discharge = sum(result["cell voltage"] .* I_app .* dt_array)

# 总产热能量
E_heat = sum(result["total heat source"] .* dt_array)

# 比率应接近（考虑充放电效率）
ratio = E_heat / E_discharge
println("Heat/Discharge ratio: ", ratio)  # 应在0.05-0.20之间
```

---

## 📚 代码结构

### 新增文件
```
src/
├── ModelInitialisation_SimpleCoupling.jl  # 初始化
│   └── ModelInitialisation_SimpleCoupling()
│
├── CallModel_SimpleCoupling.jl            # 核心求解器
│   ├── CallModel_SimpleCoupling()
│   ├── compute_total_heat_source_simple()
│   └── distribute_heat_source_uniform()
│
example/
└── testexample_simple_coupling.jl         # 测试脚本
```

### 修改的文件
```
src/
├── Solve.jl        # 添加简化耦合分支
├── Variables.jl    # 添加变量定义
└── JuBat.jl        # 导出新函数
```

---

## 🎯 最佳实践

### 适合简化耦合的场景
✅ 常规充放电仿真  
✅ 温度演化预测  
✅ 快速参数扫描  
✅ 概念设计验证  

### 不适合简化耦合的场景
❌ 热失控扩散（需要电流分布）  
❌ 局部过热分析（需要精确热源分布）  
❌ 极耳效应研究（需要电流异质性）  

### 性能优化建议
1. **网格**：在内存允许下尽量密
   ```julia
   nθ = 80  # 简化模式可以用80甚至更高
   ```

2. **时间步**：在精度允许下尽量大
   ```julia
   opt.dt = [0.5, 10.0]  # 最大步长10秒
   ```

3. **并行**：如果需要参数扫描
   ```julia
   using Distributed
   @everywhere include("JuBat.jl")
   pmap(run_case, parameter_list)
   ```

---

## 📞 获取帮助

### 文档
- 详细实施报告：`简化耦合模型实施完成报告.md`
- 开发计划：`docs/简化耦合模型开发计划.md`
- 实施步骤：`docs/简化耦合模型实施步骤.md`

### 调试
```julia
# 启用所有调试输出
opt.debug_simple_coupling = true

# 检查状态向量
println("State vector length: ", length(y0))
println("Layout: ", case.simple_coupling_layout)

# 检查中间变量
println("Available variables: ", keys(variables))
```

---

## ✅ 核对清单

运行测试前确认：
- [ ] Julia环境可用
- [ ] 所有文件已创建/修改
- [ ] 参数设置合理
- [ ] 内存充足（至少5GB）

运行成功标志：
- [ ] 无编译错误
- [ ] 无运行时错误
- [ ] 生成5张图像
- [ ] 温度范围合理
- [ ] 总结报告输出

---

**祝您使用顺利！** 🎉

如有问题，请参考详细文档或检查代码注释。
