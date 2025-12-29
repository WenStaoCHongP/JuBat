"""
诊断工具：检查极耳角度归一化问题

问题焦点：
1. theta_neg = [20π] 经过归一化后是否正确映射到网格节点
2. 网格的累计角度范围 θ_cum_out 是否覆盖到 20π
3. 归一化逻辑 (while θ0 > θc_max; θ0 -= 2π; end) 是否导致错误映射
"""

using LinearAlgebra, SparseArrays, Statistics, Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function diagnose_tab_angle_normalization()
    println("="^80)
    println("极耳角度归一化诊断")
    println("="^80)
    
    # 设置参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.tab.theta_pos = [0.0]
    param_dim.tab.theta_neg = [20.0 * π]
    
    println("\n[1] 用户设置:")
    println("  theta_pos = $(param_dim.tab.theta_pos)")
    println("  theta_neg = $(param_dim.tab.theta_neg) = $(round(20π, digits=2)) 弧度")
    println("  tab.width = $(param_dim.tab.width) m")
    
    # 创建网格
    nθ = 80
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    
    # 获取螺旋参数
    p = JuBat.jellyroll_spiral_params(param_dim)
    println("\n[2] 螺旋参数:")
    println("  a = $(p.a) m")
    println("  b = $(p.b) m")
    println("  t_repeat = $(p.t_repeat) m")
    println("  Rin = $(p.Rin) m")
    println("  Rout = $(p.Rout) m")
    
    # 计算网格节点的累计角度
    nn = size(mesh_th.node, 1)
    θ_cum_in = [(hypot(mesh_th.node[i,1], mesh_th.node[i,2]) - p.a) / p.b for i in 1:nn]
    θ_cum_out = [(hypot(mesh_th.node[i,1], mesh_th.node[i,2]) - p.a - p.t_repeat) / p.b for i in 1:nn]
    
    println("\n[3] 网格节点累计角度范围:")
    println("  内螺旋 θ_cum_in:")
    println("    范围: [$(round(minimum(θ_cum_in), digits=3)), $(round(maximum(θ_cum_in), digits=3))] 弧度")
    println("    范围: [$(round(minimum(θ_cum_in)/(2π), digits=2)), $(round(maximum(θ_cum_in)/(2π), digits=2))] 圈")
    println("  外螺旋 θ_cum_out:")
    println("    范围: [$(round(minimum(θ_cum_out), digits=3)), $(round(maximum(θ_cum_out), digits=3))] 弧度")
    println("    范围: [$(round(minimum(θ_cum_out)/(2π), digits=2)), $(round(maximum(θ_cum_out)/(2π), digits=2))] 圈")
    
    # 模拟归一化过程
    θ0_orig = 20π
    θc_min = minimum(θ_cum_out)
    θc_max = maximum(θ_cum_out)
    
    println("\n[4] 归一化过程模拟:")
    println("  原始角度 θ0_orig = $(round(θ0_orig, digits=2)) 弧度 = $(round(θ0_orig/(2π), digits=1)) 圈")
    println("  网格角度范围 [θc_min, θc_max] = [$(round(θc_min, digits=3)), $(round(θc_max, digits=3))]")
    
    θ0 = Float64(θ0_orig)
    iterations = 0
    while θ0 > θc_max
        θ0 -= 2π
        iterations += 1
        if iterations <= 3
            println("    归一化步骤 $iterations: θ0 = $(round(θ0, digits=2)) (减去 2π)")
        end
    end
    if iterations > 3
        println("    ... (共 $iterations 次迭代)")
    end
    
    println("  归一化后 θ0 = $(round(θ0, digits=3)) 弧度 = $(round(θ0/(2π), digits=2)) 圈")
    
    # 检查归一化后的角度是否在范围内
    if θ0 >= θc_min && θ0 <= θc_max
        println("  ✓ 归一化后的角度在网格范围内")
    else
        @warn "⚠️  归一化后的角度仍然超出网格范围！"
        println("    这可能导致无法识别极耳节点")
    end
    
    # 计算角度增量 Δθ
    tab_width = param_dim.tab.width
    delta_theta_fn = (θ, w) -> JuBat._delta_theta_from_width(p.a, p.b, θ, w)
    Δθ = delta_theta_fn(θ0, tab_width)
    
    println("\n[5] 极耳识别范围:")
    println("  tab.width = $(tab_width) m")
    println("  角度增量 Δθ = $(round(Δθ, digits=4)) 弧度")
    println("  负极耳识别范围（reverse_range=true）:")
    println("    [θ0 - Δθ, θ0] = [$(round(θ0 - Δθ, digits=3)), $(round(θ0, digits=3))]")
    
    # 计算该角度范围对应的半径
    r_start = p.a + p.b * (θ0 - Δθ) + p.t_repeat
    r_end = p.a + p.b * θ0 + p.t_repeat
    println("  对应半径范围: [$(round(r_start, digits=6)), $(round(r_end, digits=6))] m")
    
    # 检查是否与 Rout 接近
    if abs(r_end - p.Rout) < 0.001
        println("  ✓ 终点半径接近 Rout，符合外螺旋终点预期")
    else
        @warn "⚠️  终点半径 r_end = $(round(r_end, digits=6)) m 不接近 Rout = $(p.Rout) m"
        println("    差值: $(abs(r_end - p.Rout)) m")
    end
    
    # 统计符合条件的节点数
    count = 0
    for i in 1:nn
        r = hypot(mesh_th.node[i,1], mesh_th.node[i,2])
        θ_cum = θ_cum_out[i]
        if (p.Rin - 1e-8 <= r <= p.Rout + 1e-8) && 
           ((θ0 - Δθ) <= θ_cum <= θ0)
            count += 1
        end
    end
    
    println("\n[6] 识别结果:")
    println("  符合条件的节点数: $count")
    
    if count == 0
        @warn "❌ 未识别到任何极耳节点！"
        println("\n可能原因分析:")
        
        # 原因1：归一化导致角度映射到错误位置
        if iterations > 5
            println("  1. 归一化迭代次数过多（$iterations 次）")
            println("     原始角度 20π 被映射到了网格内的某个位置")
            println("     但该位置的半径 r = $(round(r_end, digits=6)) m 可能不在外螺旋终点")
        end
        
        # 原因2：网格角度范围不覆盖归一化后的位置
        if θ0 - Δθ < θc_min
            println("  2. 识别范围 [$(round(θ0 - Δθ, digits=3)), $(round(θ0, digits=3))] 的起点小于网格最小角度 $(round(θc_min, digits=3))")
            println("     部分范围超出了网格覆盖")
        end
        
        # 原因3：半径不匹配
        if r_end < p.Rout - 0.001 || r_start > p.Rout + 0.001
            println("  3. 半径范围 [$(round(r_start, digits=6)), $(round(r_end, digits=6))] 不在 Rout = $(p.Rout) m 附近")
            println("     可能是归一化导致识别到了内层螺旋的节点")
        end
        
        println("\n🔧 解决方案:")
        println("  问题根源：归一化逻辑将 20π 映射到网格覆盖范围内，但这会改变物理位置")
        println("  对于阿基米德螺旋，累计角度 20π 和 $(round(θ0/(2π), digits=1)) 圈的半径不同！")
        println("\n建议修改:")
        println("  1. 使用实际网格覆盖的角度范围内的值")
        println("     例如：theta_neg = [$(round(θc_max, digits=2))]")
        println("  2. 或者扩展网格的角度范围以覆盖 20π")
        println("     （需要修改 jellyroll_collector_seed_mesh 函数）")
        
    elseif count < 10
        @warn "⚠️  识别到的节点数过少（$count < 10）"
        println("  这可能导致边界条件施加不充分")
        println("\n建议:")
        println("  1. 增大 tab.width")
        println("  2. 或增大网格密度 nθ")
    else
        println("  ✓ 识别的节点数合理")
    end
    
    println("\n" * "="^80)
    println("诊断完成")
    println("="^80)
end

# 运行诊断
diagnose_tab_angle_normalization()
