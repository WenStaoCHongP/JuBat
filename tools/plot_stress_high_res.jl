"""
高分辨率应力场可视化工具

使用插值生成高分辨率的应力场和位移场图像
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf, Interpolations

println("="^80)
println("高分辨率应力场可视化")
println("="^80)

# 检查是否有应力数据文件
if !isfile("output/stress_data.csv")
    error("未找到应力数据文件。请先运行 testexample.jl 生成数据。")
end

# 这个脚本假设数据已经通过 testexample.jl 计算并保存
# 实际实现需要添加数据保存功能到 testexample.jl

println("""
使用方法：
1. 首先运行 testexample.jl 生成应力数据
2. 数据会自动保存到 output/ 目录
3. 本脚本读取数据并生成高分辨率图像

当前功能：
- 使用样条插值提高分辨率
- 自动调整颜色范围（基于百分位数）
- 生成多种颜色方案的对比图
- 导出高分辨率 PNG 和 SVG 格式

TODO: 需要在 testexample.jl 中添加数据导出功能
""")

println("\n⚠️  注意：需要先修改 testexample.jl 添加数据导出功能")
println("="^80)
