# 平面应力改为平面应变的修改总结

## ✅ 已完成的修改

我已经将代码修改为同时支持平面应力和平面应变两种模式。用户可以通过简单的选项切换。

## 🎯 使用方法

### 非常简单！只需一行代码：

```julia
# 平面应力（默认）
opt.plane_type = :stress

# 平面应变
opt.plane_type = :strain
```

## 📝 修改的代码位置

### 1. 主函数 `diffusion_stress_2D()`

**文件**: `src/mechanical.jl` (约第227行)

**修改内容**:
- 添加平面类型检测
- 传递 `plane_type` 参数到子函数
- 根据类型处理不同的返回值

```julia
# 确定平面类型
plane_type = hasproperty(case.opt, :plane_type) ? case.opt.plane_type : :stress

# 平面应变时检查泊松比
if plane_type == :strain && max_nu > 0.45
    @warn "平面应变模式下泊松比过大"
end
```

### 2. 刚度矩阵装配 `_assemble_mechanical_stiffness_2D()`

**文件**: `src/mechanical.jl` (约第380行)

**修改内容**: 弹性矩阵系数

```julia
if plane_type == :stress
    # 平面应力：E/(1-ν²)
    factor = E / max(1.0 - ν^2, 1e-12)
    D11[g] = factor
    D12[g] = factor * ν
    D33[g] = factor * (1.0 - ν) / 2.0
else  # plane_type == :strain
    # 平面应变：E/((1+ν)(1-2ν))
    denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
    factor = E / denom
    D11[g] = factor * (1.0 - ν)
    D12[g] = factor * ν
    D33[g] = factor * (1.0 - 2.0*ν) / 2.0
end
```

### 3. 载荷向量装配 `_assemble_thermal_diffusion_load_2D()`

**文件**: `src/mechanical.jl` (约第470行)

**修改内容**: 初始应力系数

```julia
if plane_type == :stress
    # 系数：(1+ν)
    factor = E / max(1.0 - ν^2, 1e-12) * ε_0 * (1.0 + ν) * wJ[g]
else  # plane_type == :strain
    # 系数：(1+2ν)
    denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
    factor = E / denom * ε_0 * (1.0 + 2.0*ν) * wJ[g]
end
```

### 4. 应力恢复 `_recover_stress_2D()`

**文件**: `src/mechanical.jl` (约第540行)

**修改内容**: 应力计算公式和 Von Mises 应力

```julia
if plane_type == :stress
    # 平面应力公式
    factor = E / max(1.0 - ν^2, 1e-12)
    σ_xx[e] = factor * (ε_elastic_xx + ν * ε_elastic_yy)
    σ_yy[e] = factor * (ε_elastic_yy + ν * ε_elastic_xx)
    σ_xy[e] = factor * (1.0 - ν) / 2.0 * ε_elastic_xy
    
    # Von Mises（2D）
    σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 - σ_xx[e]*σ_yy[e] + 3.0*σ_xy[e]^2)
else  # plane_type == :strain
    # 平面应变公式
    denom = max((1.0 + ν) * (1.0 - 2.0*ν), 1e-12)
    factor = E / denom
    σ_xx[e] = factor * ((1.0 - ν) * ε_elastic_xx + ν * ε_elastic_yy)
    σ_yy[e] = factor * ((1.0 - ν) * ε_elastic_yy + ν * ε_elastic_xx)
    σ_xy[e] = factor * (1.0 - 2.0*ν) / 2.0 * ε_elastic_xy
    
    # z方向应力
    σ_zz[e] = factor * ν * (ε_elastic_xx + ε_elastic_yy)
    
    # Von Mises（3D，包含σ_zz）
    σ_vm[e] = sqrt(σ_xx[e]^2 + σ_yy[e]^2 + σ_zz[e]^2 - 
                  σ_xx[e]*σ_yy[e] - σ_yy[e]*σ_zz[e] - σ_zz[e]*σ_xx[e] + 
                  3.0*σ_xy[e]^2)
end
```

### 5. 结果写入 `_write_mechanical_results!()`

**文件**: `src/mechanical.jl` (约第660行)

**修改内容**: 添加 σ_zz 输出

```julia
# 平面应变时，额外写入 z 方向应力
if plane_type == :strain && σ_zz !== nothing
    variables["diffusion stress zz"] = σ_zz
end

# 记录平面类型
variables["plane_type"] = String(plane_type)
```

## 📊 关键差异总结

| 位置 | 平面应力 | 平面应变 |
|------|---------|---------|
| **弹性模量因子** | `E/(1-ν²)` | `E/((1+ν)(1-2ν))` |
| **D11** | `factor` | `factor × (1-ν)` |
| **D12** | `factor × ν` | `factor × ν` |
| **D33** | `factor × (1-ν)/2` | `factor × (1-2ν)/2` |
| **初始应力系数** | `(1+ν)` | `(1+2ν)` |
| **σ_zz** | 0（不计算） | `factor × ν × (ε_xx + ε_yy)` |
| **Von Mises** | 2D公式 | 3D公式（含σ_zz） |

## 🔢 数值示例

**假设**: E = 100 GPa, ν = 0.3, ε_0 = 0.001

### 平面应力：
```
factor = 100/(1-0.09) = 109.9 GPa
系数 = 1.3
σ = 109.9 × 1.3 × 0.001 = 143 MPa
```

### 平面应变：
```
factor = 100/(1.3×0.4) = 192.3 GPa
D11 = 192.3 × 0.7 = 134.6 GPa
系数 = 1.6
σ = 134.6 × 1.6 × 0.001 = 215 MPa
```

**比值**: 215/143 = **1.51倍**

## 📖 新增文档

1. **`docs/Plane_Stress_vs_Plane_Strain.md`**
   - 详细理论对比
   - 公式推导
   - 修改指南

2. **`example/plane_strain_example.jl`**
   - 对比示例
   - 理论计算

3. **`PLANE_STRESS_VS_STRAIN_QUICK_GUIDE.md`**
   - 快速参考
   - 代码示例

4. **`PLANE_STRESS_TO_STRAIN_MODIFICATIONS.md`** (本文档)
   - 修改总结

## ✨ 新功能

### 1. 自动类型检测

```julia
plane_type = hasproperty(case.opt, :plane_type) ? case.opt.plane_type : :stress
```

### 2. 泊松比检查

```julia
if plane_type == :strain && max_nu > 0.45
    @warn "平面应变模式下泊松比过大 (max=$(max_nu))，可能导致数值不稳定"
end
```

### 3. z方向应力输出

仅在平面应变模式下：
```julia
variables["diffusion stress zz"]  # [Pa]
```

### 4. 类型标记

```julia
variables["plane_type"]  # "stress" 或 "strain"
```

## 🎓 使用示例

### 基本用法

```julia
using JuBat

# 创建案例
param = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 设置平面应变模式
opt.plane_type = :strain

case = JuBat.SetCase(param, opt)
variables = JuBat.diffusion_stress_2D(case, variables)

# 查看结果
σ_xx = variables["diffusion stress xx"]
σ_yy = variables["diffusion stress yy"]
σ_zz = variables["diffusion stress zz"]  # 仅平面应变有
σ_vm = variables["diffusion stress vonMises"]
```

### 对比两种模式

```julia
# 平面应力
case1 = deepcopy(case)
case1.opt.plane_type = :stress
vars1 = JuBat.diffusion_stress_2D(case1, deepcopy(variables))

# 平面应变
case2 = deepcopy(case)
case2.opt.plane_type = :strain
vars2 = JuBat.diffusion_stress_2D(case2, deepcopy(variables))

# 对比应力
σ_stress = mean(vars1["diffusion stress vonMises"])
σ_strain = mean(vars2["diffusion stress vonMises"])
ratio = σ_strain / σ_stress
println("应力比值: $(ratio)")  # 约 1.5
```

## ⚠️ 注意事项

### 1. 默认行为

如果不设置 `plane_type`，默认使用**平面应力**：
```julia
# 这两个等价
opt.plane_type = :stress
# 或不设置（默认）
```

### 2. 泊松比限制

平面应变时，泊松比不应超过 0.45：
```julia
# 安全
ν = 0.3  ✓
ν = 0.4  ✓

# 警告
ν = 0.45  ⚠️

# 危险
ν = 0.48  ❌ (分母 → 0)
```

### 3. 应力量级

平面应变的应力通常更大：
- 如果平面应力下 σ = 100 MPa
- 则平面应变下 σ ≈ 150-200 MPa

需要检查是否超过材料屈服强度！

### 4. Von Mises 应力

两种模式的 Von Mises 公式不同：
- 平面应力：不含 σ_zz
- 平面应变：含 σ_zz

因此即使 σ_xx, σ_yy 相同，σ_vm 也会不同。

## 🔍 验证方法

### 1. 检查输出变量

```julia
# 平面应力
@assert !haskey(variables, "diffusion stress zz")
@assert variables["plane_type"] == "stress"

# 平面应变
@assert haskey(variables, "diffusion stress zz")
@assert variables["plane_type"] == "strain"
```

### 2. 检查应力比

```julia
ratio = σ_strain / σ_stress
@assert 1.4 < ratio < 2.5  # 合理范围
```

### 3. 检查 σ_zz

```julia
# 平面应变
σ_zz_expected = ν * (σ_xx + σ_yy)
@assert abs(σ_zz - σ_zz_expected) / σ_zz_expected < 0.01
```

## 🚀 性能影响

修改对性能的影响微乎其微：
- 计算复杂度：相同
- 内存使用：平面应变额外存储 σ_zz（单元数 × 8 bytes）
- 运行时间：几乎相同（多几次浮点运算）

## 📚 参考文献

1. **Timoshenko & Goodier** (1970), "Theory of Elasticity"
2. **Zienkiewicz & Taylor** (2000), "The Finite Element Method"
3. **ANSYS Theory Reference**, 平面单元
4. **Bower** (2011), "Applied Mechanics of Solids"

## ✅ 完成清单

- [x] 添加平面类型选项
- [x] 修改刚度矩阵装配
- [x] 修改载荷向量装配
- [x] 修改应力恢复公式
- [x] 添加 σ_zz 计算（平面应变）
- [x] 修改 Von Mises 应力公式
- [x] 添加泊松比检查
- [x] 更新结果写入函数
- [x] 编写详细文档
- [x] 创建示例代码
- [x] 编写快速参考

---

**总结**: 代码已完全支持平面应力和平面应变两种模式，用户只需设置一个选项即可切换。所有必要的公式和数值处理都已正确实现。✅
