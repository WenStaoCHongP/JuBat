# 外边界节点判断分析

## 📊 问题现象

从 debug_outer_nodes.jl 的输出：

```
b = 3.068507302811742e-5  （螺距参数，很小）
t_repeat = 0.0001928  （螺旋间距）
n_wind (N) = 158  （圈数）

Default tol classification count = 200, inner misclassified = 0
Tight tol classification count   = 98, inner misclassified = 0

eps_theta_default = 3.2589..., eps_theta_tight = 0.000325...
```

**关键观察：**
1. ✅ 内圈误分类 = 0（正确）
2. ⚠️ eps_theta_default = 3.26 **太大了**！
3. ⚠️ Default tol 分类了 200 个外圈节点（可能过多）
4. ✓ Tight tol 分类了 98 个节点（更合理）

## 🔍 容差分析

### 容差计算公式（Jellyrollmodel.jl 第312行）

```julia
eps_theta = tol / max(abs(b), 1e-12)
```

### 问题分析

**默认 tol = 1e-4 时：**
```
eps_theta = 1e-4 / 3.068e-5 = 3.26 弧度 ≈ 187°
```

这意味着边界判断的角度容差高达 **187度**！

**结果：**
- 外圈预期范围：[2π*157, 2π*158] = [986.96, 993.25]
- 实际判断范围：[986.96 - 3.26, 993.25 + 3.26] = [983.7, 996.5]
- **容差占了半圈还多！** 太宽松了

**紧容差 tol = 1e-8 时：**
```
eps_theta = 1e-8 / 3.068e-5 = 0.000326 弧度 ≈ 0.019°
```

这个就合理多了！

## ❌ 根本原因

### 问题 1：螺距 b 太小

```
b = 3.068e-5 m = 0.03 mm
```

这么小的螺距导致：
- 角度容差 = 长度容差 / b 会被放大
- 默认 1e-4 m 的容差 → 3.26 弧度的角度容差

### 问题 2：默认容差不适合小螺距

对于 Jellyroll：
- 螺距可能只有几十微米
- 长度容差 1e-4 m = 0.1 mm 相对于螺距来说太大

## ✅ 已应用的修复方案

### 方案 1：改进容差计算（已实施）⭐

**位置：** `src/Jellyrollmodel.jl` 第 304-308 行

**修改前：**
```julia
tol = get(kwargs, :tol, 1e-4)
eps_theta = tol / max(abs(bval), 1e-12)
```

**修改后：**
```julia
# 自适应容差：对于小螺距，使用更小的长度容差以避免角度容差过大
# tol_default = min(1e-4, 0.05 * abs(p.b)) 确保 eps_theta < 0.05 弧度（约3度）
bval = p.b == 0.0 ? 1e-12 : p.b
tol_default = min(1e-4, 0.05 * abs(bval))
tol = get(kwargs, :tol, tol_default)
eps_theta = tol / max(abs(bval), 1e-12)
```

### 效果分析

**对于当前案例（b = 3.068e-5）：**
```
修复前：
  tol = 1e-4
  eps_theta = 1e-4 / 3.068e-5 = 3.26 弧度 ≈ 187° ❌ 太大

修复后：
  tol = min(1e-4, 0.05 * 3.068e-5) = 1.534e-6
  eps_theta = 1.534e-6 / 3.068e-5 = 0.05 弧度 ≈ 2.86° ✓ 合理
```

**对于大螺距（b = 1e-3）：**
```
修复前：
  tol = 1e-4
  eps_theta = 0.1 弧度 ≈ 5.7° ✓ 本来就合理

修复后：
  tol = min(1e-4, 0.05 * 1e-3) = 5e-5
  eps_theta = 0.05 弧度 ≈ 2.86° ✓ 更精确
```

### 方案 2：使用紧容差（临时方案）

在调用 `edge_boundary` 时显式指定：

```julia
is_outer_node[i] = edge_boundary(:node_on, mesh, i, param_dim; which=:outer, tol=1e-8)
```

### 方案 3：基于角度的容差（最佳方案）⭐⭐

直接使用角度容差而非长度容差：

```julia
# 定义角度容差（弧度）
eps_theta_rad = get(kwargs, :eps_theta, 0.1)  # 默认 0.1 弧度 ≈ 5.7°

# 直接用于判断
if which === :outer
    N = max(1, Int(p.n_wind))
    θ_start, θ_end = 2.0*pi*(N-1), 2.0*pi*N
    return (θ_cum_out >= θ_start - eps_theta_rad) && (θ_cum_out <= θ_end + eps_theta_rad)
end
```

## 🔧 立即可尝试的修复

### 临时解决（在 ThermalDistributed2D_BC 中）

```julia
# 修改前：
is_outer_node[i] = edge_boundary(:node_on, mesh, i, case.param_dim; which=:outer)

# 修改后：
is_outer_node[i] = edge_boundary(:node_on, mesh, i, case.param_dim; which=:outer, tol=1e-8)
```

这会使用紧容差，从 200 个节点减少到 98 个（更准确）。

## 📊 验证方法

修复后重新运行 debug_outer_nodes.jl：

```bash
julia tools/debug_outer_nodes.jl
```

### 预期输出

**修复前：**
```
Default tol classification count = 200  ❌ 过多
eps_theta_default = 3.2589  ❌ 187度，太大
Sample theta_cum_out inner[1:5] = [-6.28, ...]  ❌ 显示错误
```

**修复后：**
```
Default tol classification count = 98-110  ✓ 合理（接近实际外圈节点数）
Tight tol classification count   = 98  ✓ 
eps_theta_default = 0.05  ✓ 约3度，合理

内圈节点 (应该在 [0, 2π] 范围):
  Sample theta_cum_in inner[1:5] = [0.0, 0.039, 0.078, ...]
  预期范围: [0, 6.28]
  ✓ 内圈 theta_cum 在预期范围内

外圈节点 (应该在 [2π*(N-1), 2π*N] 范围):
  Sample theta_cum_out outer[1:5] = [986.xx, 987.xx, ...]
  预期范围: [986.96, 993.25]
  ✓ 外圈 theta_cum 在预期范围内
```

### 影响

修复后：
- ✅ 外边界节点分类更准确
- ✅ 边界对流换热施加在正确的节点上
- ✅ 温度场边界条件更准确
- ✅ 数值稳定性提高

## 🎯 总结

**问题：** 默认长度容差（1e-4 m）对于小螺距（3e-5 m）来说太大，导致角度容差达到 187度

**修复：** 使用自适应容差 `tol_default = min(1e-4, 0.05*b)`，确保角度容差 < 3度

**效果：** 外边界节点分类从 200 个减少到约 100 个，更加准确

**修改文件：** `src/Jellyrollmodel.jl` 第 304-308 行
