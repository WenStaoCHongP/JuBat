include("../src/JuBat.jl")

"""
诊断边界节点遗漏问题
"""

param_dim = JuBat.ChooseCell("Jellyroll")
p = JuBat.jellyroll_spiral_params(param_dim)
mesh = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)

println("="^80)
println("边界节点遗漏问题诊断")
println("="^80)

# 计算网格实际覆盖的θ范围
s_in = 0.0
s_out = p.t_repeat
θ0_mesh = max(0.0, (p.Rin - p.a - s_in) / p.b)
θ1_mesh = min((p.Rout - p.a - s_out) / p.b, (p.Rout - p.a) / p.b)

println("\n1. 网格实际θ范围:")
println("   θ0 = $θ0_mesh rad")
println("   θ1 = $θ1_mesh rad ($(θ1_mesh * 180 / π)°)")
println("   覆盖圈数 = $(θ1_mesh / (2π))")

# ThermalDistributed使用的θ范围
N = Int(p.n_wind)
θ_in_thermal = (0.0, 2.0*π)
θ_out_thermal = (2.0*π*(N-1), 2.0*π*N)

println("\n2. ThermalDistributed中使用的θ范围:")
println("   n_wind = $N (向下取整)")
println("   内边界: $(θ_in_thermal)")
println("   外边界: $(θ_out_thermal)")
println("   外边界范围: [$(θ_out_thermal[1] / (2π)), $(θ_out_thermal[2] / (2π))] 圈")

# 关键问题：终点不一致
θ_diff = θ1_mesh - θ_out_thermal[2]
println("\n3. ⚠️  问题所在:")
println("   网格终点 θ1 = $θ1_mesh")
println("   外边界终点 2π*N = $(θ_out_thermal[2])")
println("   差值 = $θ_diff rad = $(θ_diff / (2π)) 圈")
println("   ==> 最外圈有 $(θ_diff / (2π)) 圈的节点会被遗漏！")

# 计算实际遗漏的节点
x = mesh.node[:, 1]
y = mesh.node[:, 2]
r = hypot.(x, y)

# 外螺旋的所有节点
outer_spiral_nodes = Int[]
for i in 1:mesh.nlen
    θ_cum = (r[i] - p.a - p.t_repeat) / p.b
    # 检查是否在螺旋线上
    r_theo = p.a + p.b * θ_cum + p.t_repeat
    dist = abs(r[i] - r_theo)
    if dist < 1e-4 && θ0_mesh <= θ_cum <= θ1_mesh
        push!(outer_spiral_nodes, i)
    end
end

# 使用ThermalDistributed范围识别的节点
outer_thermal_nodes = Int[]
for i in 1:mesh.nlen
    θ_cum = (r[i] - p.a - p.t_repeat) / p.b
    r_theo = p.a + p.b * θ_cum + p.t_repeat
    dist = abs(r[i] - r_theo)
    if dist < 1e-4 && θ_out_thermal[1] <= θ_cum <= θ_out_thermal[2]
        push!(outer_thermal_nodes, i)
    end
end

missing_nodes = setdiff(outer_spiral_nodes, outer_thermal_nodes)

println("\n4. 节点统计:")
println("   网格上所有外螺旋节点数: $(length(outer_spiral_nodes))")
println("   ThermalDistributed识别到的: $(length(outer_thermal_nodes))")
println("   遗漏节点数: $(length(missing_nodes))")

if !isempty(missing_nodes)
    println("\n5. 遗漏节点分析:")
    θ_cum_missing = [(r[i] - p.a - p.t_repeat) / p.b for i in missing_nodes]
    println("   遗漏节点的θ_cum范围: [$(minimum(θ_cum_missing)), $(maximum(θ_cum_missing))]")
    println("   对应圈数: [$(minimum(θ_cum_missing) / (2π)), $(maximum(θ_cum_missing) / (2π))]")
    println("   这些节点在 θ ∈ (2π*N, θ1_mesh] 范围内")
end

println("\n"*"="^80)
println("解决方案")
println("="^80)
println("\n方案1: 修改ThermalDistributed中的θ范围计算")
println("   外边界应使用: (2π*(N-1), θ1_mesh)")
println("   而不是:      (2π*(N-1), 2π*N)")
println("\n方案2: 使用网格的实际θ范围")
println("   提供函数计算网格实际覆盖的θ范围")
println("   让ThermalDistributed自动使用正确的范围")
println("\n"*"="^80)
