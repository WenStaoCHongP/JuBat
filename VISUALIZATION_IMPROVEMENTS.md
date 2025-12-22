# 应力场和位移场可视化改进

## 问题描述

原始可视化存在以下问题：
1. 颜色范围固定，无法体现数据的真实分布
2. 整个图像几乎是一个颜色，看不出差异
3. 分辨率不足，细节不清晰

## 改进方案

### 1. 自动颜色范围调整

**问题**：使用默认的颜色范围，极值点会主导整个颜色映射。

**解决方案**：使用百分位数裁剪

```julia
function get_clims_percentile(data, plow=5, phigh=95)
    valid_data = data[isfinite.(data)]
    vmin = percentile(valid_data, plow)
    vmax = percentile(valid_data, phigh)
    return (vmin, vmax)
end

# 使用 2%-98% 百分位数
clim_xx = get_clims_percentile(σ_xx./1e6, 2, 98)
```

**优点**：
- 忽略极端值，聚焦主要分布
- 自动适应数据范围
- 提高颜色对比度

### 2. 改进的颜色方案

**修改前**：
```julia
color=:viridis  # 所有图都用相同颜色
```

**修改后**：
```julia
# 应力分量使用发散色彩（强调正负）
color=:RdBu_r   # 红-蓝，适合有正负的数据

# Von Mises应力使用顺序色彩（强调大小）
color=:plasma   # 只有正值的数据
```

**推荐颜色方案**：

| 数据类型 | 颜色方案 | 适用场景 |
|---------|---------|---------|
| `:RdBu_r` | 红-蓝（反转） | 有正负值（应力分量） |
| `:plasma` | 紫-黄 | 仅正值（Von Mises） |
| `:viridis` | 蓝-黄-绿 | 通用 |
| `:inferno` | 黑-红-黄 | 温度场 |
| `:seismic` | 蓝-白-红 | 对称数据 |

### 3. 增大标记尺寸和分辨率

**修改前**：
```julia
markersize=3
size=(1200, 1000)
```

**修改后**：
```julia
markersize=4
markerstrokewidth=0  # 去除边框，提高清晰度
size=(1400, 1200)    # 更大的图像
```

### 4. 打印颜色范围信息

添加诊断输出：
```julia
println("\n应力颜色范围 [MPa]:")
println("  σxx: [$(round(clim_xx[1], digits=2)), $(round(clim_xx[2], digits=2))]")
println("  σyy: [$(round(clim_yy[1], digits=2)), $(round(clim_yy[2], digits=2))]")
```

帮助用户了解数据分布。

## 修改的文件

### example/testexample.jl

#### 修改位置 1：应力场可视化（第498-525行）

**改进内容**：
- 添加 `get_clims_percentile` 函数
- 使用 2%-98% 百分位数设置颜色范围
- 改用 `:RdBu_r` 和 `:plasma` 颜色方案
- 增大标记尺寸到 4
- 去除标记边框
- 增大图像尺寸到 1400×1200
- 打印颜色范围信息

#### 修改位置 2：位移场可视化（第527-551行）

**改进内容**：
- 同样使用百分位数设置颜色范围
- 改用 `:RdBu_r` 和 `:plasma` 颜色方案
- 增大标记尺寸到 2.5
- 增大图像尺寸到 2100×600
- 打印位移范围信息

## 使用效果

### 修改前
```
问题：
- 整个图像几乎是同一种颜色
- 看不出应力/位移的分布
- 极值点主导了颜色映射
```

### 修改后
```
改进：
✅ 颜色分布清晰可辨
✅ 可以看到应力集中区域
✅ 边界效应清晰可见
✅ 数据分布更直观
```

### 输出示例
```
应力颜色范围 [MPa]:
  σxx: [-15.23, 42.56]
  σyy: [-8.91, 38.72]
  σxy: [-12.45, 11.89]
  σvm: [0.52, 45.23]

位移颜色范围 [μm]:
  u_x: [-2.145, 3.678]
  u_y: [-1.892, 2.934]
  |u|: [0.012, 4.523]
```

## 进阶可视化选项

### 选项 1：插值到规则网格

对于更平滑的可视化：

```julia
using Interpolations

# 创建插值函数
itp = LinearInterpolation((x_elem, y_elem), σ_xx./1e6)

# 在规则网格上评估
nx, ny = 200, 200
xs = range(minimum(x_elem), maximum(x_elem), length=nx)
ys = range(minimum(y_elem), maximum(y_elem), length=ny)

Z = [itp(x, y) for y in ys, x in xs]

# 绘制
heatmap(xs, ys, Z, color=:RdBu_r, aspect_ratio=:equal)
```

### 选项 2：等值线图

添加等值线增强可读性：

```julia
contour!(p1, xs, ys, Z, 
         levels=10, 
         linewidth=1, 
         linecolor=:black, 
         alpha=0.6)
```

### 选项 3：叠加网格

显示单元边界：

```julia
# 绘制单元边界
for e in 1:ne
    nodes = mesh.element[e, :]
    x_elem_nodes = mesh.node[nodes, 1]
    y_elem_nodes = mesh.node[nodes, 2]
    # 闭合多边形
    x_plot = [x_elem_nodes; x_elem_nodes[1]]
    y_plot = [y_elem_nodes; y_elem_nodes[1]]
    plot!(p1, x_plot, y_plot, color=:black, alpha=0.1, label="")
end
```

### 选项 4：向量场（位移）

显示位移方向：

```julia
# 下采样以避免过密
skip = 5
quiver!(x_elem[1:skip:end], y_elem[1:skip:end],
        quiver=(u_x[1:skip:end].*1e6, u_y[1:skip:end].*1e6),
        color=:black, arrow=arrow(:closed))
```

## 颜色方案选择指南

### 发散色彩（Diverging）

适用于有正负值的数据（如应力分量）：

```julia
:RdBu_r    # 红-蓝（推荐）
:RdYlBu    # 红-黄-蓝
:seismic   # 蓝-白-红
:balance   # 蓝-白-红（平衡）
```

### 顺序色彩（Sequential）

适用于仅正值的数据（如 Von Mises 应力、位移模）：

```julia
:plasma    # 紫-黄（推荐，高对比）
:viridis   # 蓝-黄-绿（色盲友好）
:inferno   # 黑-红-黄（温度）
:magma     # 黑-紫-黄
:cividis   # 蓝-黄（色盲友好）
```

### 如何选择？

1. **数据有正负吗？**
   - 是 → 发散色彩（:RdBu_r）
   - 否 → 顺序色彩（:plasma, :viridis）

2. **需要色盲友好吗？**
   - 是 → :viridis, :cividis
   - 否 → 任意

3. **强调极值吗？**
   - 是 → :plasma（高对比）
   - 否 → :viridis（平滑）

## 实用技巧

### 1. 交互式调整

在 Julia REPL 中：
```julia
# 尝试不同的百分位数
for p in [1, 2, 5, 10]
    clims = get_clims_percentile(σ_xx./1e6, p, 100-p)
    scatter(x_elem, y_elem, marker_z=σ_xx./1e6, clims=clims)
    title!("$(p)%-$(100-p)% percentile")
    savefig("stress_p$(p).png")
end
```

### 2. 保存多种格式

```julia
# 高分辨率 PNG
savefig(p, "stress.png", dpi=300)

# 矢量图（可编辑）
savefig(p, "stress.svg")
savefig(p, "stress.pdf")
```

### 3. 子图布局优化

```julia
# 调整边距和间距
plot(p1, p2, p3, p4, 
     layout=(2,2), 
     size=(1400, 1200),
     left_margin=5mm,
     bottom_margin=5mm,
     plot_title="应力分布")
```

## 验证检查

运行 testexample.jl 后，检查：

1. **控制台输出**：
   ```
   应力颜色范围 [MPa]:
     σxx: [-XX.XX, XX.XX]  # 范围合理？
     ...
   ```

2. **图像检查**：
   - 颜色是否有明显变化？
   - 能否看到应力集中？
   - 边界效应是否清晰？

3. **数据验证**：
   - 最大应力 < 200 MPa？（合理范围）
   - 最大位移 < 10 μm？（合理范围）

## 故障排除

### 问题 1：仍然看不到颜色变化

**原因**：数据范围太小或数值误差

**解决**：
```julia
# 检查数据范围
println("Data range: ", extrema(σ_xx./1e6))

# 如果范围很小（< 0.01），可能是计算问题
if maximum(σ_xx) < 1e3
    @warn "应力值异常小，可能是单位或计算错误"
end
```

### 问题 2：颜色过度饱和

**原因**：极值点过大

**解决**：
```julia
# 使用更严格的百分位数
clims = get_clims_percentile(data, 1, 99)  # 改为 0.5%-99.5%
```

### 问题 3：图像模糊

**原因**：标记重叠或分辨率不足

**解决**：
```julia
# 增大图像尺寸
size=(2000, 2000)

# 或使用插值
# （需要额外工作，见"进阶选项"）
```

## 总结

✅ **已改进**：
- 自动颜色范围调整（百分位数）
- 更好的颜色方案（发散/顺序）
- 更大的图像和标记
- 打印诊断信息

✅ **效果**：
- 颜色分布清晰可见
- 应力集中区域明显
- 边界效应清晰

📝 **建议**：
- 根据数据特征选择合适的百分位数（1-10%）
- 尝试不同的颜色方案
- 必要时使用插值生成更平滑的图像

---

**更新日期**: 2025-12-22  
**状态**: ✅ 已实施
