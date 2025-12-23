#!/usr/bin/env julia
# 测试热应力修复：验证温度场历史是否正确保存和使用

println("=" ^ 60)
println("热应力修复验证测试")
println("=" ^ 60)

# 1. 检查 Variables.jl 是否正确添加了温度历史记录
println("\n[1/4] 检查 Variables.jl...")
variables_content = read("src/Variables.jl", String)
if occursin("variables[\"thermal2D temperature\"] = zeros(Float64, nT, num)", variables_content)
    println("✓ Variables.jl 已添加 thermal2D temperature 历史记录")
    # 检查是否在多SPMe模式下添加
    if occursin("per_element_spme", split(variables_content, "thermal2D temperature")[1])
        println("✓ 多SPMe模式下正确初始化")
    else
        println("⚠ 可能未在多SPMe模式下初始化")
    end
else
    println("✗ Variables.jl 未找到温度历史记录初始化")
end

# 2. 检查 Solve.jl 是否正确保存温度到历史
println("\n[2/4] 检查 Solve.jl...")
solve_content = read("src/Solve.jl", String)
if occursin("variables[\"thermal2D temperature\"] = T_nodes .* Tref", solve_content)
    println("✓ Solve.jl 已添加温度场历史保存")
    if occursin("haskey(variables_hist, \"thermal2D temperature\")", solve_content)
        println("✓ 使用了安全检查")
    end
else
    println("✗ Solve.jl 未找到温度历史保存代码")
end

# 3. 检查 testexample.jl 是否正确使用温度历史
println("\n[3/4] 检查 testexample.jl...")
testexample_content = read("example/testexample.jl", String)
if occursin("result[\"thermal2D temperature [K]\"]", testexample_content)
    println("✓ testexample.jl 已修改为使用温度历史")
    if occursin("T_nodes_hist_K[:, step]", testexample_content)
        println("✓ 正确使用对应时间步的温度")
    else
        println("⚠ 可能未正确索引时间步")
    end
else
    println("⚠ testexample.jl 可能仍使用最终温度")
end

# 4. 检查函数名修复
println("\n[4/4] 检查函数名修复...")
spme_example_content = read("example/spme_thermal2d_example.jl", String)
plot_stress_content = read("tools/plot_thermal_stress.jl", String)

fixed_count = 0
if occursin("thermal_diffusion_stress_2D", spme_example_content)
    println("✓ spme_thermal2d_example.jl 函数名已修复")
    fixed_count += 1
else
    println("✗ spme_thermal2d_example.jl 仍使用错误函数名")
end

if occursin("thermal_diffusion_stress_2D", plot_stress_content)
    println("✓ plot_thermal_stress.jl 函数名已修复")
    fixed_count += 1
else
    println("✗ plot_thermal_stress.jl 仍使用错误函数名")
end

# 总结
println("\n" * "=" ^ 60)
println("测试总结")
println("=" ^ 60)
println("修复内容：")
println("  1. Variables.jl - 添加温度场历史记录初始化")
println("  2. Solve.jl - 每步保存温度场到历史")
println("  3. testexample.jl - 使用温度历史计算应力")
println("  4. 函数名修复 - thermal_stress → thermal_diffusion_stress_2D")
println()
println("预期效果：")
println("  • 每个时间步使用对应的温度场")
println("  • 温度变化 ΔT = T(t) - T₀ 正确计算")
println("  • 热应力 σ_thermal ≠ 0")
println("  • 应力历史反映真实的热-力耦合")
println()
println("建议：")
println("  运行 example/testexample.jl 验证热应力是否非零")
println("=" ^ 60)
