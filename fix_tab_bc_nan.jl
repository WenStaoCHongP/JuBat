"""
极耳热边界条件NaN问题快速修复

用法：
在主程序 testexample.jl 中，在创建 opt 后添加：

    include("fix_tab_bc_nan.jl")
    apply_tab_bc_fix!(opt, param_dim, case)

作者：AI Assistant
日期：2025-12-29
"""

"""
应用极耳BC修复方案

参数：
- opt: Option对象
- param_dim: Params对象
- case: Case对象（可选，用于网格诊断）
"""
function apply_tab_bc_fix!(opt, param_dim, case=nothing)
    println("\n" * "="^70)
    println("应用极耳热边界条件修复")
    println("="^70)
    
    # ========================================================================
    # 修复1：降低惩罚系数（核心修复）
    # ========================================================================
    
    original_penalty = hasproperty(opt, :tab_penalty) ? opt.tab_penalty : 1e12
    new_penalty = 1e10  # 降低到1e10
    
    opt.tab_penalty = new_penalty
    
    println("\n✓ 修复1：降低惩罚系数")
    println("  原值: $(original_penalty)")
    println("  新值: $(new_penalty)")
    println("  说明: 降低刚度对比，提高数值稳定性")
    
    # ========================================================================
    # 修复2：设置合理的升温速率
    # ========================================================================
    
    if !hasproperty(opt, :tab_heating_rate)
        opt.tab_heating_rate = 0.0  # 默认为0（散热边界）
        println("\n✓ 修复2：设置升温速率")
        println("  值: 0.0 K/s（极耳保持环境温度）")
    else
        println("\n✓ 修复2：保留用户设置的升温速率")
        println("  值: $(opt.tab_heating_rate) K/s")
    end
    
    # ========================================================================
    # 修复3：启用调试模式
    # ========================================================================
    
    opt.debug_coupling = true
    println("\n✓ 修复3：启用调试模式")
    println("  将输出极耳节点识别信息")
    
    # ========================================================================
    # 诊断：检查极耳参数
    # ========================================================================
    
    println("\n" * "-"^70)
    println("极耳参数诊断")
    println("-"^70)
    
    # 检查极耳角度
    n_pos_angles = length(param_dim.tab.theta_pos)
    n_neg_angles = length(param_dim.tab.theta_neg)
    
    println("  极耳宽度: $(param_dim.tab.width*1e3) mm")
    println("  正极耳角度数: $n_pos_angles")
    if n_pos_angles > 0
        println("    角度值: $(param_dim.tab.theta_pos)")
        println("    角度值(度): $([rad2deg(θ) for θ in param_dim.tab.theta_pos])")
    end
    println("  负极耳角度数: $n_neg_angles")
    if n_neg_angles > 0
        println("    角度值: $(param_dim.tab.theta_neg)")
        println("    角度值(度): $([rad2deg(θ) for θ in param_dim.tab.theta_neg])")
    end
    
    # 警告：无极耳
    if n_pos_angles == 0 && n_neg_angles == 0
        @warn "未设置极耳角度，极耳BC不会生效"
        println("\n  ⚠️  极耳角度为空，极耳边界条件将被跳过")
        println("  如需启用极耳BC，请设置:")
        println("    param_dim.tab.theta_pos = [0.0, π]")
        println("    param_dim.tab.theta_neg = [π/2, 3π/2]")
    end
    
    # 警告：极耳过多
    if n_pos_angles + n_neg_angles > 8
        @warn "极耳数量过多 ($(n_pos_angles + n_neg_angles) 个)"
        println("  ⚠️  极耳数量较多，可能增加计算负担")
        println("  推荐：每侧1-4个极耳")
    end
    
    # ========================================================================
    # 诊断：检查网格（如果提供了case）
    # ========================================================================
    
    if case !== nothing && haskey(case.mesh, "thermal2D")
        println("\n" * "-"^70)
        println("网格诊断")
        println("-"^70)
        
        mesh_th = case.mesh["thermal2D"]
        nn = size(mesh_th.node, 1)
        ne = size(mesh_th.element, 1)
        
        println("  节点数: $nn")
        println("  单元数: $ne")
        
        # 尝试识别极耳节点
        if n_pos_angles > 0 || n_neg_angles > 0
            try
                # 需要 JuBat 模块
                if @isdefined JuBat
                    pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
                    n_pos = length(pos_idx)
                    n_neg = length(neg_idx)
                    n_tab = n_pos + n_neg
                    
                    println("\n  极耳节点识别:")
                    println("    正极耳: $n_pos 节点 ($(round(100*n_pos/nn, digits=1))%)")
                    println("    负极耳: $n_neg 节点 ($(round(100*n_neg/nn, digits=1))%)")
                    println("    合计: $n_tab 节点 ($(round(100*n_tab/nn, digits=1))%)")
                    
                    # 判断合理性
                    tab_ratio = n_tab / nn
                    if tab_ratio > 0.3
                        @warn "极耳节点占比过高 ($(round(tab_ratio*100, digits=1))%)"
                        println("\n  ⚠️  极耳节点占比异常高 (>30%)")
                        println("  可能原因:")
                        println("    1. tab.width 过大")
                        println("    2. 极耳角度设置错误")
                        println("    3. 网格角度分辨率不足")
                        println("  建议:")
                        println("    - 减小 tab.width (当前: $(param_dim.tab.width*1e3) mm)")
                        println("    - 检查极耳角度是否合理")
                    elseif tab_ratio < 0.01 && n_tab > 0
                        @warn "极耳节点占比过低 ($(round(tab_ratio*100, digits=1))%)"
                        println("\n  ⚠️  极耳节点占比过低 (<1%)")
                        println("  可能原因: 网格角度分辨率不足")
                        println("  建议: 增加 n_angular 到 60 以上")
                    elseif n_tab == 0
                        @warn "未识别到极耳节点"
                        println("\n  ⚠️  极耳节点数为0")
                        println("  可能原因:")
                        println("    1. 极耳角度超出网格覆盖范围")
                        println("    2. tab.width 过小")
                        println("    3. 网格未正确生成")
                    else
                        println("\n  ✅ 极耳节点占比正常 (5-15%)")
                    end
                    
                    # 可视化建议
                    if n_tab > 0
                        println("\n  提示: 可使用以下代码可视化极耳节点:")
                        println("    using Plots")
                        println("    x, y = mesh_th.node[:, 1], mesh_th.node[:, 2]")
                        println("    scatter(x, y, ms=2, alpha=0.3, label=\"网格\")")
                        println("    scatter!(x[pos_idx], y[pos_idx], ms=4, color=:red, label=\"正极耳\")")
                        println("    scatter!(x[neg_idx], y[neg_idx], ms=4, color=:blue, label=\"负极耳\")")
                        println("    savefig(\"tab_nodes.png\")")
                    end
                else
                    println("  (跳过极耳节点识别，JuBat模块未加载)")
                end
            catch e
                println("  (跳过极耳节点识别，出错: $e)")
            end
        end
    end
    
    # ========================================================================
    # 总结
    # ========================================================================
    
    println("\n" * "="^70)
    println("修复应用完成")
    println("="^70)
    println("\n下一步:")
    println("  1. 运行仿真，观察是否仍出现NaN")
    println("  2. 检查调试输出中的极耳节点信息")
    println("  3. 如问题仍存在，尝试进一步降低 opt.tab_penalty")
    println("\n当前配置:")
    println("  opt.tab_penalty = $(opt.tab_penalty)")
    println("  opt.tab_heating_rate = $(opt.tab_heating_rate)")
    println("  opt.debug_coupling = $(opt.debug_coupling)")
    println("="^70)
    println()
end

# 导出函数
export apply_tab_bc_fix!

println("✓ 极耳BC修复模块已加载")
println("  使用: apply_tab_bc_fix!(opt, param_dim, case)")
