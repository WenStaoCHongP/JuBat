using Statistics
include("../src/JuBat.jl")

function main()
    println("=" ^ 80)
    println("精确边界识别测试：基于螺旋线方程")
    println("=" ^ 80)

    param_dim = JuBat.ChooseCell("Jellyroll")
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)
    p = JuBat.jellyroll_spiral_params(param_dim)

    println("\n网格信息:")
    println("  总节点数: $(mesh_th.nlen)")
    println("  Rin = $(p.Rin) m, Rout = $(p.Rout) m")
    println("  螺距 b = $(p.b) m/rad")
    println("  层厚 t_repeat = $(p.t_repeat) m")
    println("  总圈数 N = $(p.n_wind)")

    x_all = mesh_th.node[:,1]
    y_all = mesh_th.node[:,2]
    r_all = hypot.(x_all, y_all)

    println("\n节点半径分布:")
    println("  最小半径: $(minimum(r_all)) m")
    println("  最大半径: $(maximum(r_all)) m")

    println("\n" * "=" ^ 80)
    println("测试1: 识别第1圈的内螺旋（θ ∈ [0, 2π]）")
    println("=" ^ 80)

    θ_range_inner = (0.0, 2.0*pi)
    is_inner_precise = [JuBat.edge_boundary(mesh_th, i, param_dim;
                                            which=:inner, theta_range=θ_range_inner, tol=1e-4)
                        for i in 1:mesh_th.nlen]

    count_inner = count(identity, is_inner_precise)
    println("  识别到的内边界节点数: $count_inner")

    if count_inner > 0
        r_inner = r_all[is_inner_precise]
        println("  内边界节点半径范围: [$(minimum(r_inner)), $(maximum(r_inner))] m")

        θ_cum_inner = [(r_all[i] - p.a) / p.b for i in 1:mesh_th.nlen if is_inner_precise[i]]
        println("  θ_cum 范围: [$(minimum(θ_cum_inner)), $(maximum(θ_cum_inner))]")
        println("  预期范围: [0.0, $(2*pi)]")

        max_dist = 0.0
        for i in 1:mesh_th.nlen
            if is_inner_precise[i]
                x, y = mesh_th.node[i,1], mesh_th.node[i,2]
                r = hypot(x, y)
                θ_cum = (r - p.a) / p.b
                r_theo = p.a + p.b * θ_cum
                x_theo = r_theo * cos(θ_cum)
                y_theo = r_theo * sin(θ_cum)
                dist = hypot(x - x_theo, y - y_theo)
                max_dist = max(max_dist, dist)
            end
        end
        println("  节点到螺旋线的最大距离: $(max_dist) m")
    end

    println("\n" * "=" ^ 80)
    println("测试2: 识别第N圈的外螺旋（θ ∈ [2π(N-1), 2π*N]）")
    println("=" ^ 80)

    N = Int(p.n_wind)
    θ_range_outer = (2.0*pi*(N-1), 2.0*pi*N)
    println("  第N圈范围: θ ∈ [$(θ_range_outer[1]), $(θ_range_outer[2])]")

    is_outer_precise = [JuBat.edge_boundary(mesh_th, i, param_dim;
                                            which=:outer, theta_range=θ_range_outer, tol=1e-4)
                        for i in 1:mesh_th.nlen]

    count_outer = count(identity, is_outer_precise)
    println("  识别到的外边界节点数: $count_outer")

    if count_outer > 0
        r_outer = r_all[is_outer_precise]
        println("  外边界节点半径范围: [$(minimum(r_outer)), $(maximum(r_outer))] m")

        θ_cum_outer = [(r_all[i] - p.a - p.t_repeat) / p.b for i in 1:mesh_th.nlen if is_outer_precise[i]]
        println("  θ_cum 范围: [$(minimum(θ_cum_outer)), $(maximum(θ_cum_outer))]")
        println("  预期范围: [$(2*pi*(N-1)), $(2*pi*N)]")

        max_dist = 0.0
        for i in 1:mesh_th.nlen
            if is_outer_precise[i]
                x, y = mesh_th.node[i,1], mesh_th.node[i,2]
                r = hypot(x, y)
                θ_cum = (r - p.a - p.t_repeat) / p.b
                r_theo = p.a + p.b * θ_cum + p.t_repeat
                x_theo = r_theo * cos(θ_cum)
                y_theo = r_theo * sin(θ_cum)
                dist = hypot(x - x_theo, y - y_theo)
                max_dist = max(max_dist, dist)
            end
        end
        println("  节点到螺旋线的最大距离: $(max_dist) m")
    else
        println("  ⚠️  没有识别到外边界节点！")
        println("  原因：collector_seed_mesh 只覆盖第0圈，不覆盖第N圈")
    end

    println("\n" * "=" ^ 80)
    println("测试3: 识别第0圈的外螺旋（θ ∈ [0, 2π]）- collector_seed_mesh 实际覆盖")
    println("=" ^ 80)

    θ_range_outer_0 = (0.0, 2.0*pi)
    is_outer_0 = [JuBat.edge_boundary(mesh_th, i, param_dim;
                                      which=:outer, theta_range=θ_range_outer_0, tol=1e-4)
                  for i in 1:mesh_th.nlen]

    count_outer_0 = count(identity, is_outer_0)
    println("  识别到的外螺旋节点数（第0圈）: $count_outer_0")

    if count_outer_0 > 0
        r_outer_0 = r_all[is_outer_0]
        println("  外螺旋节点半径范围: [$(minimum(r_outer_0)), $(maximum(r_outer_0))] m")

        θ_cum_outer_0 = [(r_all[i] - p.a - p.t_repeat) / p.b for i in 1:mesh_th.nlen if is_outer_0[i]]
        println("  θ_cum 范围: [$(minimum(θ_cum_outer_0)), $(maximum(θ_cum_outer_0))]")
        println("  预期范围: [0.0, $(2*pi)]")

        max_dist = 0.0
        for i in 1:mesh_th.nlen
            if is_outer_0[i]
                x, y = mesh_th.node[i,1], mesh_th.node[i,2]
                r = hypot(x, y)
                θ_cum = (r - p.a - p.t_repeat) / p.b
                r_theo = p.a + p.b * θ_cum + p.t_repeat
                x_theo = r_theo * cos(θ_cum)
                y_theo = r_theo * sin(θ_cum)
                dist = hypot(x - x_theo, y - y_theo)
                max_dist = max(max_dist, dist)
            end
        end
        println("  节点到螺旋线的最大距离: $(max_dist) m")
    end

    println("\n" * "=" ^ 80)
    println("测试4: 检查重叠（同时被识别为内边界和外边界的节点）")
    println("=" ^ 80)

    overlap = count(i -> is_inner_precise[i] && is_outer_0[i], 1:mesh_th.nlen)
    println("  重叠节点数（第1圈内螺旋 ∩ 第0圈外螺旋）: $overlap")
    if overlap == 0
        println("  ✓ 无重叠，边界分类正确")
    end

    println("\n" * "=" ^ 80)
    println("测试5: 与 collector_seed_mesh 网格拓扑结构对比")
    println("=" ^ 80)

    inner_curve_range = 1:(160+1)
    outer_curve_range = (160+2):(2*(160+1))

    count_inner_topo = 161
    count_outer_topo = 161

    println("  网格拓扑结构:")
    println("    内圈节点（拓扑）: 1-161 (共161个)")
    println("    外圈节点（拓扑）: 162-322 (共161个)")

    println("\n  精确识别结果（第0圈）:")
    println("    内螺旋（θ ∈ [0, 2π]）: $count_inner 个节点")
    println("    外螺旋（θ ∈ [0, 2π]）: $count_outer_0 个节点")

    match_inner = count(i -> is_inner_precise[i], inner_curve_range)
    match_outer = count(i -> is_outer_0[i], outer_curve_range)

    println("\n  匹配度:")
    println("    内圈节点中被识别为内螺旋的: $match_inner / $count_inner_topo")
    println("    外圈节点中被识别为外螺旋的: $match_outer / $count_outer_topo")

    missing_inner = setdiff(inner_curve_range, findall(is_inner_precise))
    extra_inner = setdiff(findall(is_inner_precise), inner_curve_range)
    missing_outer = setdiff(outer_curve_range, findall(is_outer_0))
    extra_outer = setdiff(findall(is_outer_0), outer_curve_range)

    if !isempty(missing_inner)
        println("    ⚠ Missing inner indices: $(collect(missing_inner))")
    end
    if !isempty(extra_inner)
        println("    ⚠ Extra inner indices: $(collect(extra_inner)[1:min(end,10)]) ...")
    end
    if !isempty(missing_outer)
        println("    ⚠ Missing outer indices: $(collect(missing_outer))")
    end
    if !isempty(extra_outer)
        preview = collect(extra_outer)[1:min(end,10)]
        println("    ⚠ Extra outer indices (first 10): $preview")
    end

    if match_inner == count_inner_topo && match_outer == count_outer_topo
        println("    ✓ 完全匹配！精确识别正确。")
    elseif match_inner == count_inner && match_outer == count_outer_0
        println("    ✓ 高度匹配（可能有边界节点未精确在螺旋线上）")
    else
        println("    ⚠️ 部分匹配，请检查容差设置")
    end

    println("\n" * "=" ^ 80)
    println("总结")
    println("=" ^ 80)

    println("\n✓ 精确边界识别方法特点:")
    println("  1. 基于螺旋线方程 r(θ) = a + b*θ (+ t_repeat)")
    println("  2. 计算节点到理论螺旋线的距离，验证 dist < tol")
    println("  3. 可指定任意 θ 范围，如 [0, 2π]（第1圈）或 [2π(N-1), 2π*N]（第N圈）")

    println("\n⚠️ collector_seed_mesh 的限制:")
    println("  - 网格只覆盖第0圈（或某一圈），θ ∈ [0, 2π]")
    println("  - 内外圈在同一圈层，只是半径不同")
    println("  - 无法识别第1圈内螺旋和第N圈外螺旋（网格不覆盖）")

    println("\n💡 建议:")
    if count_outer == 0
        println("  - 对于 collector_seed_mesh，使用 θ_range=(0, 2π) 识别外螺旋")
        println("  - 或创建覆盖全螺旋的网格（θ ∈ [0, 2π*N]）")
    end

    println("\n" * "=" ^ 80)
end

main()