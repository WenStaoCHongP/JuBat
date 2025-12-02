include("../src/JuBat.jl")

"""
验证边界节点识别修复 - 正确版本
只识别第一圈和最后一圈
"""

param_dim = JuBat.ChooseCell("Jellyroll")
nθ = 160
mesh = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
p = JuBat.jellyroll_spiral_params(param_dim)

println("="^80)
println("边界节点识别验证 - 只识别第一圈/最后一圈")
println("="^80)

# 计算网格实际θ范围
s_in = 0.0
s_out = p.t_repeat
bval = max(p.b, 1e-12)
θ0_mesh = max(0.0, (p.Rin - p.a - s_in) / bval)
θ1_mesh = min((p.Rout - p.a - s_out) / bval, (p.Rout - p.a) / bval)

println("\n网格参数:")
println("  节点总数: $(mesh.nlen)")
println("  分段数 nθ: $nθ")
println("  网格θ范围: [$θ0_mesh, $θ1_mesh]")
println("  覆盖圈数: $(θ1_mesh / (2π))")

# 内边界应该识别的θ范围
θ_inner_expected = (θ0_mesh, min(θ0_mesh + 2.0*π, θ1_mesh))
println("\n内边界（第一圈）:")
println("  θ范围: $(θ_inner_expected)")
println("  圈数范围: [$(θ_inner_expected[1]/(2π)), $(θ_inner_expected[2]/(2π))]")

# 外边界应该识别的θ范围
θ_outer_expected = (max(θ1_mesh - 2.0*π, θ0_mesh), θ1_mesh)
println("\n外边界（最后一圈）:")
println("  θ范围: $(θ_outer_expected)")
println("  圈数范围: [$(θ_outer_expected[1]/(2π)), $(θ_outer_expected[2]/(2π))]")

println("\n"*"="^80)
println("测试1: 使用edge_boundary默认参数识别")
println("="^80)

# 使用默认参数识别（不传theta_range）
inner_nodes_default = Int[]
outer_nodes_default = Int[]

for i in 1:mesh.nlen
    if JuBat.edge_boundary(mesh, i, param_dim; which=:inner)
        push!(inner_nodes_default, i)
    end
    if JuBat.edge_boundary(mesh, i, param_dim; which=:outer)
        push!(outer_nodes_default, i)
    end
end

println("\n使用默认参数:")
println("  识别到的内边界节点数: $(length(inner_nodes_default))")
println("  识别到的外边界节点数: $(length(outer_nodes_default))")

# 检查网格拓扑结构中的节点
# 对于collector_seed_mesh，前nθ+1个节点应该在内螺旋线上
# 后nθ+1个节点应该在外螺旋线上
inner_topo = 1:(nθ+1)
outer_topo = (nθ+2):(2*(nθ+1))

inner_topo_matched = length(intersect(inner_nodes_default, inner_topo))
outer_topo_matched = length(intersect(outer_nodes_default, outer_topo))

println("\n与网格拓扑对比:")
println("  内圈节点（拓扑）: 1-$(nθ+1) (共$(nθ+1)个)")
println("  外圈节点（拓扑）: $(nθ+2)-$(2*(nθ+1)) (共$(nθ+1)个)")
println("  识别到的内圈节点数: $inner_topo_matched / $(nθ+1)")
println("  识别到的外圈节点数: $outer_topo_matched / $(nθ+1)")

# 预期：如果网格覆盖超过1圈，只有第一圈的内螺旋节点被识别
# 因为collector_seed_mesh可能采样多圈，但我们只要第一圈
x = mesh.node[:, 1]
y = mesh.node[:, 2]
r = hypot.(x, y)

# 检查识别到的节点的θ范围
if !isempty(inner_nodes_default)
    θ_cum_inner = [(r[i] - p.a - 0.0) / bval for i in inner_nodes_default]
    println("\n内边界节点θ范围:")
    println("  实际: [$(minimum(θ_cum_inner)), $(maximum(θ_cum_inner))]")
    println("  预期: $(θ_inner_expected)")
    println("  圈数: [$(minimum(θ_cum_inner)/(2π)), $(maximum(θ_cum_inner)/(2π))]")
end

if !isempty(outer_nodes_default)
    θ_cum_outer = [(r[i] - p.a - p.t_repeat) / bval for i in outer_nodes_default]
    println("\n外边界节点θ范围:")
    println("  实际: [$(minimum(θ_cum_outer)), $(maximum(θ_cum_outer))]")
    println("  预期: $(θ_outer_expected)")
    println("  圈数: [$(minimum(θ_cum_outer)/(2π)), $(maximum(θ_cum_outer)/(2π))]")
end

println("\n"*"="^80)
println("测试2: 检查节点是否只在第一圈/最后一圈")
println("="^80)

# 验证内边界节点都在第一圈
inner_in_first_turn = true
for i in inner_nodes_default
    θ_cum = (r[i] - p.a - 0.0) / bval
    if θ_cum < θ_inner_expected[1] || θ_cum > θ_inner_expected[2]
        inner_in_first_turn = false
        println("  ⚠️  内边界节点 $i 的θ_cum=$θ_cum 超出预期范围")
        break
    end
end

if inner_in_first_turn
    println("  ✓ 所有内边界节点都在第一圈范围内")
else
    println("  ⚠️  部分内边界节点超出第一圈范围")
end

# 验证外边界节点都在最后一圈
outer_in_last_turn = true
for i in outer_nodes_default
    θ_cum = (r[i] - p.a - p.t_repeat) / bval
    if θ_cum < θ_outer_expected[1] || θ_cum > θ_outer_expected[2]
        outer_in_last_turn = false
        println("  ⚠️  外边界节点 $i 的θ_cum=$θ_cum 超出预期范围")
        break
    end
end

if outer_in_last_turn
    println("  ✓ 所有外边界节点都在最后一圈范围内")
else
    println("  ⚠️  部分外边界节点超出最后一圈范围")
end

println("\n"*"="^80)
println("测试3: 检查是否有遗漏的第一圈/最后一圈节点")
println("="^80)

# 找出所有应该在第一圈的内螺旋节点
should_be_inner = Int[]
for i in 1:mesh.nlen
    θ_cum = (r[i] - p.a - 0.0) / bval
    r_theo = p.a + p.b * θ_cum + 0.0
    dist = abs(r[i] - r_theo)
    if dist < 1e-4 && θ_inner_expected[1] <= θ_cum <= θ_inner_expected[2]
        push!(should_be_inner, i)
    end
end

# 找出所有应该在最后一圈的外螺旋节点
should_be_outer = Int[]
for i in 1:mesh.nlen
    θ_cum = (r[i] - p.a - p.t_repeat) / bval
    r_theo = p.a + p.b * θ_cum + p.t_repeat
    dist = abs(r[i] - r_theo)
    if dist < 1e-4 && θ_outer_expected[1] <= θ_cum <= θ_outer_expected[2]
        push!(should_be_outer, i)
    end
end

missing_inner = setdiff(should_be_inner, inner_nodes_default)
missing_outer = setdiff(should_be_outer, outer_nodes_default)

println("\n应识别的节点数:")
println("  内边界（第一圈）: $(length(should_be_inner))")
println("  外边界（最后一圈）: $(length(should_be_outer))")

println("\n实际识别的节点数:")
println("  内边界: $(length(inner_nodes_default))")
println("  外边界: $(length(outer_nodes_default))")

if isempty(missing_inner)
    println("\n✓ 内边界：所有第一圈节点都被正确识别")
else
    println("\n⚠️  内边界遗漏节点: $(length(missing_inner)) 个")
    if length(missing_inner) <= 10
        println("  遗漏节点索引: $(collect(missing_inner))")
    end
end

if isempty(missing_outer)
    println("✓ 外边界：所有最后一圈节点都被正确识别")
else
    println("⚠️  外边界遗漏节点: $(length(missing_outer)) 个")
    if length(missing_outer) <= 10
        println("  遗漏节点索引: $(collect(missing_outer))")
    end
end

println("\n"*"="^80)
println("总结")
println("="^80)

all_pass = inner_in_first_turn && 
           outer_in_last_turn && 
           isempty(missing_inner) && 
           isempty(missing_outer) &&
           (length(inner_nodes_default) == length(should_be_inner)) &&
           (length(outer_nodes_default) == length(should_be_outer))

if all_pass
    println("\n✅ 所有测试通过！")
    println("\n边界节点识别正确:")
    println("  - 内边界只识别第一圈: $(length(inner_nodes_default)) 个节点")
    println("  - 外边界只识别最后一圈: $(length(outer_nodes_default)) 个节点")
    println("  - 无遗漏、无超出范围的节点")
else
    println("\n⚠️  部分测试未通过:")
    if !inner_in_first_turn
        println("  - 内边界有节点超出第一圈范围")
    end
    if !outer_in_last_turn
        println("  - 外边界有节点超出最后一圈范围")
    end
    if !isempty(missing_inner)
        println("  - 内边界遗漏 $(length(missing_inner)) 个节点")
    end
    if !isempty(missing_outer)
        println("  - 外边界遗漏 $(length(missing_outer)) 个节点")
    end
end

println("\n"*"="^80)
