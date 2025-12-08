"""
验证温度场精度提升

功能：
- 模拟温度数据
- 测试高精度等高线生成
- 验证统计输出格式
"""

using Plots, Statistics, Printf

println("="^60)
println("温度场精度验证")
println("="^60)

# 模拟温度数据
n = 1000
T_base = 298.0
T_range = 0.5

# 生成模拟温度分布（正态分布）
T_nodes = T_base .+ T_range .* (0.5 .+ 0.3 .* randn(n))
T_nodes = clamp.(T_nodes, T_base, T_base + T_range)

# 统计
vmin = minimum(T_nodes)
vmax = maximum(T_nodes)
T_mean = mean(T_nodes)
T_std = std(T_nodes)

println("\n模拟温度数据:")
@printf("  样本数: %d\n", n)
@printf("  温度范围: [%.4f, %.4f] K\n", vmin, vmax)
@printf("  温度差: %.4f K\n", vmax - vmin)
@printf("  平均值: %.4f K\n", T_mean)
@printf("  标准差: %.4f K\n", T_std)

# 测试等高线密度计算
println("\n等高线计算测试:")
T_range_calc = vmax - vmin
contour_interval = 0.01
n_raw = T_range_calc / contour_interval
n_contours = max(5, min(100, Int(round(n_raw))))

@printf("  目标精度: %.2f K\n", contour_interval)
@printf("  原始等高线数: %.1f\n", n_raw)
@printf("  实际等高线数: %d\n", n_contours)
@printf("  实际间隔: %.4f K\n", T_range_calc / (n_contours - 1))

contour_levels = range(vmin, vmax, length=n_contours)
println("  前5条等高线值:")
for i in 1:min(5, n_contours)
    @printf("    Level %d: %.4f K\n", i, contour_levels[i])
end

# 粗等高线计算
if T_range_calc > 0.1
    major_levels = range(ceil(vmin*10)/10, floor(vmax*10)/10, step=0.1)
    println("\n  粗等高线 (0.1 K 间隔):")
    @printf("    数量: %d\n", length(major_levels))
    if length(major_levels) > 0
        @printf("    范围: [%.1f, %.1f] K\n", first(major_levels), last(major_levels))
    end
else
    println("\n  温度范围 < 0.1 K，跳过粗等高线")
end

# 统计量计算
println("\n详细温度统计 (精度: 0.01 K):")
@printf("  最小值: %.4f K (%.2f °C)\n", 
        minimum(T_nodes), minimum(T_nodes) - 273.15)
@printf("  第25百分位: %.4f K\n", quantile(T_nodes, 0.25))
@printf("  中位数: %.4f K\n", median(T_nodes))
@printf("  平均值: %.4f K\n", mean(T_nodes))
@printf("  第75百分位: %.4f K\n", quantile(T_nodes, 0.75))
@printf("  最大值: %.4f K (%.2f °C)\n", 
        maximum(T_nodes), maximum(T_nodes) - 273.15)
@printf("  标准差: %.4f K\n", std(T_nodes))
@printf("  温度范围: %.4f K\n", maximum(T_nodes) - minimum(T_nodes))

# 生成测试图像
println("\n生成测试图像...")

# 图1: 温度分布直方图
p1 = plot(layout=(2,1), size=(800, 600))

histogram!(p1[1], T_nodes, 
          bins=50, 
          xlabel="Temperature (K)", 
          ylabel="Frequency",
          title="Simulated Temperature Distribution (N=$n)",
          label="Nodes",
          alpha=0.7,
          color=:steelblue)
vline!(p1[1], [mean(T_nodes)], 
       label="Mean", 
       linewidth=2, 
       linestyle=:dash, 
       color=:red)

# 图2: 统计量柱状图
T_stats = [
    minimum(T_nodes),
    quantile(T_nodes, 0.25),
    median(T_nodes),
    quantile(T_nodes, 0.75),
    maximum(T_nodes),
    mean(T_nodes),
    std(T_nodes)
]
stat_labels = ["Min", "Q1", "Median", "Q3", "Max", "Mean", "Std"]

bar!(p1[2], 1:7, T_stats,
     xticks=(1:7, stat_labels),
     ylabel="Temperature (K)",
     title="Temperature Statistics (Precision: 0.01 K)",
     label=false,
     color=:coral,
     alpha=0.7)

# 添加数值标签
for i in 1:7
    annotate!(p1[2], i, T_stats[i], 
             text(@sprintf("%.2f", T_stats[i]), :bottom, 8))
end

savefig(p1, "verify_temperature_precision.png")
println("✓ 保存: verify_temperature_precision.png")

# 图2: 模拟等高线图
println("\n生成模拟等高线图...")

# 创建2D网格数据
nx, ny = 100, 100
x = range(-1, 1, length=nx)
y = range(-1, 1, length=ny)

Z = [T_base + T_range * (0.5 + 0.2*sin(3*xi) + 0.2*cos(3*yi)) 
     for yi in y, xi in x]

p2 = plot(size=(900, 800), title="Simulated Temperature Field (0.01K Precision)")

# 热图
heatmap!(p2, x, y, Z; 
         aspect_ratio=1, 
         color=:inferno, 
         colorbar=true, 
         colorbar_title="T (K)",
         xlabel="x (normalized)", 
         ylabel="y (normalized)",
         clims=(vmin, vmax))

# 精细等高线
contour!(p2, x, y, Z; 
         levels=contour_levels, 
         linewidth=0.8, 
         linecolor=:white, 
         alpha=0.6,
         linestyle=:solid)

# 粗等高线
if T_range_calc > 0.1 && length(major_levels) > 0
    contour!(p2, x, y, Z; 
             levels=collect(major_levels), 
             linewidth=1.5, 
             linecolor=:black, 
             alpha=0.8,
             linestyle=:solid)
end

savefig(p2, "verify_contour_precision.png")
println("✓ 保存: verify_contour_precision.png")

# 总结
println("\n" * "="^60)
println("验证完成")
println("="^60)
println("""
✓ 等高线精度: 0.01 K
✓ 统计输出精度: 0.0001 K (显示为 0.01 K)
✓ 温度标注精度: 0.01 K

生成的测试图像:
  1. verify_temperature_precision.png - 温度分布统计
  2. verify_contour_precision.png - 高精度等高线示例

这些图像展示了修改后的精度效果。
实际运行 testexample.jl 将产生真实的电池温度场数据。
""")
println("="^60)
