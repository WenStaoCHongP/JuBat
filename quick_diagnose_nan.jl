"""
快速诊断：第一步NaN问题

用途：定位极耳BC导致NaN的根本原因

使用方法：
在主程序中，创建网格后立即调用：
    include("quick_diagnose_nan.jl")
    quick_diagnose_nan(mesh_th, param_dim, case)

作者：AI Assistant
日期：2025-12-29
"""

function quick_diagnose_nan(mesh, param_dim, case=nothing)
    println("\n" * "="^70)
    println("🚨 快速诊断：极耳BC第一步NaN问题")
    println("="^70)
    
    # ========================================================================
    # 1. 网格基本信息
    # ========================================================================
    println("\n1. 网格信息:")
    nn = size(mesh.node, 1)
    ne = size(mesh.element, 1)
    println("  节点数: $nn")
    println("  单元数: $ne")
    
    # ========================================================================
    # 2. 极耳参数检查
    # ========================================================================
    println("\n2. 极耳参数:")
    println("  width: $(param_dim.tab.width*1e3) mm")
    println("  theta_pos: $(param_dim.tab.theta_pos)")
    println("  theta_neg: $(param_dim.tab.theta_neg)")
    
    if isempty(param_dim.tab.theta_pos) && isempty(param_dim.tab.theta_neg)
        println("\n  ℹ️  极耳角度为空，BC不会应用")
        return
    end
    
    # ========================================================================
    # 3. 极耳节点识别
    # ========================================================================
    println("\n3. 极耳节点识别:")
    
    try
        if !@isdefined(JuBat)
            @error "JuBat模块未加载"
            return
        end
        
        pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh, param_dim)
        n_pos = length(pos_idx)
        n_neg = length(neg_idx)
        n_tab = n_pos + n_neg
        
        println("  正极耳: $n_pos 节点 ($(round(100*n_pos/nn, digits=1))%)")
        println("  负极耳: $n_neg 节点 ($(round(100*n_neg/nn, digits=1))%)")
        println("  合计: $n_tab 节点 ($(round(100*n_tab/nn, digits=1))%)")
        
        # ====================================================================
        # 4. 问题诊断
        # ====================================================================
        println("\n4. 问题诊断:")
        
        issues = []
        critical = false
        
        # 检查1：极耳节点为0
        if n_tab == 0
            push!(issues, "⚠️  极耳节点数为0，BC将被跳过")
        end
        
        # 检查2：极耳节点占比过高（关键！）
        tab_ratio = n_tab / nn
        if tab_ratio > 0.5
            push!(issues, "🔴 致命：极耳节点占比>50% → 识别错误！")
            critical = true
        elseif tab_ratio > 0.3
            push!(issues, "🔴 严重：极耳节点占比>30% → 可能识别错误")
            critical = true
        elseif tab_ratio > 0.15
            push!(issues, "⚠️  警告：极耳节点占比>15% → 偏高")
        end
        
        # 检查3：极耳节点等于全部节点
        if n_tab == nn
            push!(issues, "🔴 致命：所有节点都被识别为极耳！")
            push!(issues, "     → 这会导致整个矩阵被消元法破坏")
            push!(issues, "     → 必然产生NaN")
            critical = true
        end
        
        # 输出诊断结果
        if isempty(issues)
            println("  ✅ 极耳节点识别正常")
        else
            for issue in issues
                println("  $issue")
            end
        end
        
        # ====================================================================
        # 5. 根本原因分析
        # ====================================================================
        if critical
            println("\n" * "="^70)
            println("🔍 根本原因分析")
            println("="^70)
            
            println("\n极耳节点识别异常的可能原因：")
            
            # 原因1：width过大
            width_mm = param_dim.tab.width * 1e3
            if width_mm > 100
                println("  🔴 原因1：tab.width过大 ($(width_mm) mm)")
                println("     建议：将width减小到20-50 mm")
                println("     修复：param_dim.tab.width = 40e-3")
            end
            
            # 原因2：theta范围过大
            if length(param_dim.tab.theta_pos) > 10
                println("  🔴 原因2：theta_pos角度过多 ($(length(param_dim.tab.theta_pos))个)")
                println("     建议：减少极耳数量到1-4个")
            end
            
            if length(param_dim.tab.theta_neg) > 10
                println("  🔴 原因3：theta_neg角度过多 ($(length(param_dim.tab.theta_neg))个)")
                println("     建议：减少极耳数量到1-4个")
            end
            
            # 原因3：theta设置错误
            all_theta = vcat(param_dim.tab.theta_pos, param_dim.tab.theta_neg)
            if length(all_theta) > 0
                theta_range = maximum(all_theta) - minimum(all_theta)
                if theta_range > 2π
                    println("  🔴 原因4：theta角度范围过大 ($(round(theta_range, digits=2)) rad)")
                    println("     可能包含重复或错误的角度")
                end
            end
            
            # ================================================================
            # 6. 立即修复建议
            # ================================================================
            println("\n" * "="^70)
            println("💡 立即修复建议")
            println("="^70)
            
            println("\n选项1：临时禁用极耳BC（验证是否为极耳BC导致）")
            println("  param_dim.tab.theta_pos = []")
            println("  param_dim.tab.theta_neg = []")
            
            println("\n选项2：修正极耳参数")
            println("  # 减小宽度")
            println("  param_dim.tab.width = 40e-3  # 40 mm")
            println("  ")
            println("  # 设置合理的极耳角度（示例）")
            println("  param_dim.tab.theta_pos = [0.0, π]       # 2个正极耳")
            println("  param_dim.tab.theta_neg = [π/2, 3π/2]    # 2个负极耳")
            
            println("\n选项3：检查网格生成")
            println("  # 可能网格参数与极耳参数不匹配")
            println("  # 检查 n_angular 是否足够（推荐≥60）")
            
        else
            # 极耳识别正常，但仍有NaN → 其他问题
            println("\n" * "="^70)
            println("🔍 其他可能原因")
            println("="^70)
            
            println("\n极耳节点识别正常，但仍出现NaN，可能原因：")
            println("  1. 消元法实现有bug")
            println("  2. 矩阵初始状态异常")
            println("  3. 初始温度场未正确设置")
            println("  4. 稀疏矩阵操作问题")
            
            println("\n建议添加详细诊断：")
            println("  参见：诊断极耳BC第一步NaN问题.md")
        end
        
    catch e
        @error "诊断失败" exception=e
    end
    
    println("\n" * "="^70)
end

# 导出函数
export quick_diagnose_nan

println("✓ 快速诊断模块已加载")
println("  使用: quick_diagnose_nan(mesh_th, param_dim, case)")
