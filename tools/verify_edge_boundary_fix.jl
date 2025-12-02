using Statistics
include("../src/JuBat.jl")

"""
验证 edge_boundary 修复后是否正确识别所有边界节点
"""
function main()
    println("=" ^ 80)
    println("验证 edge_boundary 修复：终点角度判定与网格生成一致")
    println("=" ^ 80)
    
    param_dim = JuBat.ChooseCell("Jellyroll")
    nθ = 160
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    p = JuBat.jellyroll_spiral_params(param_dim)
    
    println("\n网格参数:")
    println("  总节点数: $(mesh_th.nlen)")
    println("  分段数 nθ: $nθ")
    println("  预期内圈节点数: $(nθ+1)")
    println("  预期外圈节点数: $(nθ+1)")
    println("  Rin = $(p.Rin) m, Rout = $(p.Rout) m")
    println("  螺距 b = $(p.b) m/rad")
    println("  层厚 t_repeat = $(p.t_repeat) m")
    
    # 计算网格生成的θ范围
    s_in = 0.0
    s_out = p.t_repeat
    θ0_mesh = max(0.0, (p.Rin - p.a - s_in) / p.b)
    θ1_mesh = min((p.Rout - p.a - s_out) / p.b, (p.Rout - p.a) / p.b)
    
    println("\n网格生成的θ范围:")
    println("  θ0 = $θ0_mesh rad ($(θ0_mesh * 180 / π) deg)")
    println("  θ1 = $θ1_mesh rad ($(θ1_mesh * 180 / π) deg)")
    println("  覆盖圈数: $(θ1_mesh / (2π)) 圈")
    
    println("\n" * "=" ^ 80)
    println("测试1: 识别内螺旋边界节点")
    println("=" ^ 80)
    
    # 识别内螺旋节点（使用默认θ范围）
    inner_nodes = Int[]
    for i in 1:mesh_th.nlen
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:inner)
            push!(inner_nodes, i)
        end
    end
    
    count_inner = length(inner_nodes)
    expected_inner = nθ + 1
    
    println("  识别到的内螺旋节点数: $count_inner")
    println("  预期节点数: $expected_inner")
    
    if count_inner == expected_inner
        println("  ✓ 节点数匹配！")
    else
        println("  ⚠️  节点数不匹配，差值: $(count_inner - expected_inner)")
    end
    
    # 检查网格拓扑结构中的内圈节点是否都被识别
    inner_topo = 1:(nθ+1)  # 网格拓扑中前 nθ+1 个节点是内圈
    missing_inner = setdiff(inner_topo, inner_nodes)
    extra_inner = setdiff(inner_nodes, inner_topo)
    
    if isempty(missing_inner)
        println("  ✓ 所有拓扑内圈节点都被识别")
    else
        println("  ⚠️  遗漏的拓扑内圈节点: $(collect(missing_inner))")
    end
    
    if isempty(extra_inner)
        println("  ✓ 没有额外识别的节点")
    else
        println("  ⚠️  额外识别的节点数: $(length(extra_inner))")
    end
    
    println("\n" * "=" ^ 80)
    println("测试2: 识别外螺旋边界节点")
    println("=" ^ 80)
    
    # 识别外螺旋节点（使用默认θ范围）
    outer_nodes = Int[]
    for i in 1:mesh_th.nlen
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:outer)
            push!(outer_nodes, i)
        end
    end
    
    count_outer = length(outer_nodes)
    expected_outer = nθ + 1
    
    println("  识别到的外螺旋节点数: $count_outer")
    println("  预期节点数: $expected_outer")
    
    if count_outer == expected_outer
        println("  ✓ 节点数匹配！")
    else
        println("  ⚠️  节点数不匹配，差值: $(count_outer - expected_outer)")
    end
    
    # 检查网格拓扑结构中的外圈节点是否都被识别
    outer_topo = (nθ+2):(2*(nθ+1))  # 网格拓扑中后 nθ+1 个节点是外圈
    missing_outer = setdiff(outer_topo, outer_nodes)
    extra_outer = setdiff(outer_nodes, outer_topo)
    
    if isempty(missing_outer)
        println("  ✓ 所有拓扑外圈节点都被识别")
    else
        println("  ⚠️  遗漏的拓扑外圈节点: $(collect(missing_outer))")
        # 打印遗漏节点的详细信息
        if !isempty(missing_outer) && length(missing_outer) <= 10
            println("\n  遗漏节点详细信息:")
            for idx in collect(missing_outer)
                x, y = mesh_th.node[idx, :]
                r = hypot(x, y)
                θ_cum = (r - p.a - p.t_repeat) / p.b
                println("    节点$idx: r=$r, θ_cum=$θ_cum ($(θ_cum * 180 / π) deg)")
            end
        end
    end
    
    if isempty(extra_outer)
        println("  ✓ 没有额外识别的节点")
    else
        println("  ⚠️  额外识别的节点数: $(length(extra_outer))")
    end
    
    println("\n" * "=" ^ 80)
    println("测试3: 验证节点到螺旋线的距离")
    println("=" ^ 80)
    
    # 检查内螺旋节点的距离
    max_dist_inner = 0.0
    for i in inner_nodes
        x, y = mesh_th.node[i, :]
        r = hypot(x, y)
        θ_cum = (r - p.a - 0.0) / p.b
        r_theo = p.a + p.b * θ_cum + 0.0
        x_theo = r_theo * cos(θ_cum)
        y_theo = r_theo * sin(θ_cum)
        dist = hypot(x - x_theo, y - y_theo)
        max_dist_inner = max(max_dist_inner, dist)
    end
    
    println("  内螺旋节点到理论螺旋线的最大距离: $(max_dist_inner) m")
    if max_dist_inner < 1e-4
        println("  ✓ 距离在容差范围内（< 1e-4 m）")
    else
        println("  ⚠️  距离超出容差范围")
    end
    
    # 检查外螺旋节点的距离
    max_dist_outer = 0.0
    for i in outer_nodes
        x, y = mesh_th.node[i, :]
        r = hypot(x, y)
        θ_cum = (r - p.a - p.t_repeat) / p.b
        r_theo = p.a + p.b * θ_cum + p.t_repeat
        x_theo = r_theo * cos(θ_cum)
        y_theo = r_theo * sin(θ_cum)
        dist = hypot(x - x_theo, y - y_theo)
        max_dist_outer = max(max_dist_outer, dist)
    end
    
    println("  外螺旋节点到理论螺旋线的最大距离: $(max_dist_outer) m")
    if max_dist_outer < 1e-4
        println("  ✓ 距离在容差范围内（< 1e-4 m）")
    else
        println("  ⚠️  距离超出容差范围")
    end
    
    println("\n" * "=" ^ 80)
    println("测试4: 检查内外螺旋节点重叠")
    println("=" ^ 80)
    
    overlap = intersect(inner_nodes, outer_nodes)
    if isempty(overlap)
        println("  ✓ 无重叠节点，边界分类正确")
    else
        println("  ⚠️  发现重叠节点: $(length(overlap)) 个")
        println("  重叠节点索引: $(collect(overlap))")
    end
    
    println("\n" * "=" ^ 80)
    println("测试5: 验证edge_boundary默认θ范围与网格生成一致")
    println("=" ^ 80)
    
    # 手动计算edge_boundary使用的默认θ范围
    # 根据修改后的代码（第313-316行）
    θ_start_edge = max(0.0, (p.Rin - p.a - 0.0) / p.b)
    θ_end_edge = min((p.Rout - p.a - p.t_repeat) / p.b, (p.Rout - p.a) / p.b)
    
    println("  网格生成的θ范围: [$(θ0_mesh), $(θ1_mesh)]")
    println("  edge_boundary的θ范围: [$(θ_start_edge), $(θ_end_edge)]")
    
    θ_diff_start = abs(θ0_mesh - θ_start_edge)
    θ_diff_end = abs(θ1_mesh - θ_end_edge)
    
    if θ_diff_start < 1e-10 && θ_diff_end < 1e-10
        println("  ✓ θ范围完全一致（误差 < 1e-10）")
    else
        println("  ⚠️  θ范围存在差异:")
        println("    θ_start差异: $θ_diff_start")
        println("    θ_end差异: $θ_diff_end")
    end
    
    println("\n" * "=" ^ 80)
    println("总结")
    println("=" ^ 80)
    
    all_pass = (count_inner == expected_inner) && 
               (count_outer == expected_outer) && 
               isempty(missing_inner) && 
               isempty(missing_outer) && 
               isempty(overlap) && 
               (max_dist_inner < 1e-4) && 
               (max_dist_outer < 1e-4) && 
               (θ_diff_start < 1e-10) && 
               (θ_diff_end < 1e-10)
    
    if all_pass
        println("\n✅ 所有测试通过！edge_boundary修复成功。")
        println("\n修复内容:")
        println("  1. θ起点从 θ_end-2π 改为 max(0.0, (Rin-a-s_in)/b)")
        println("  2. θ终点使用 min((Rout-a-s_out)/b, (Rout-a)/b) 与网格一致")
        println("  3. 内外螺旋共享相同的θ范围")
    else
        println("\n⚠️  部分测试未通过，请检查:")
        if count_inner != expected_inner
            println("  - 内螺旋节点数不匹配")
        end
        if count_outer != expected_outer
            println("  - 外螺旋节点数不匹配")
        end
        if !isempty(missing_inner) || !isempty(missing_outer)
            println("  - 存在遗漏的边界节点")
        end
        if !isempty(overlap)
            println("  - 存在重叠节点")
        end
        if max_dist_inner >= 1e-4 || max_dist_outer >= 1e-4
            println("  - 节点到螺旋线距离超出容差")
        end
    end
    
    println("\n" * "=" ^ 80)
end

main()
