"""
诊断工具：检查极耳边界条件的数值问题

该工具用于诊断：
1. 极耳角度是否在网格覆盖范围内
2. 识别的极耳节点数量
3. 惩罚法导致的矩阵条件数问题
4. T_ref 和边界条件计算的数值稳定性
"""

using LinearAlgebra, SparseArrays, Statistics, Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function diagnose_tab_bc_issue()
    println("="^80)
    println("极耳边界条件数值问题诊断")
    println("="^80)
    
    # ========================================================================
    # 步骤 1：设置参数（与 testexample.jl 相同）
    # ========================================================================
    println("\n[步骤 1] 参数设置...")
    
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    
    # ⚠️ 用户的极耳设置
    param_dim.tab.theta_pos = [0.0]
    param_dim.tab.theta_neg = [20.0 * π]  # 约 62.83 弧度
    
    println("  极耳角度设置:")
    println("    theta_pos = $(param_dim.tab.theta_pos)")
    println("    theta_neg = $(param_dim.tab.theta_neg)")
    println("    20π = $(20π) ≈ $(round(20π, digits=2)) 弧度 ≈ $(round(20π/(2π), digits=1)) 圈")
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.Current = x -> 5.0
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.time = [0.0, 60.0]
    opt.dt = [0.5, 10]
    
    case = JuBat.SetCase(param_dim, opt)
    
    # ========================================================================
    # 步骤 2：创建网格并检查角度范围
    # ========================================================================
    println("\n[步骤 2] 创建 Jellyroll 网格并检查角度覆盖范围...")
    
    nθ = 80
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("  网格信息:")
    println("    总单元数 ne = $ne")
    println("    总节点数 nT = $nT")
    println("    周向单元数 nθ = $nθ")
    
    # 计算节点的累计角度范围
    p = JuBat.jellyroll_spiral_params(param_dim)
    nn = size(mesh_th.node, 1)
    θ_cum_in = [(hypot(mesh_th.node[i,1], mesh_th.node[i,2]) - p.a) / p.b for i in 1:nn]
    θ_cum_out = [(hypot(mesh_th.node[i,1], mesh_th.node[i,2]) - p.a - p.t_repeat) / p.b for i in 1:nn]
    
    println("\n  螺旋参数:")
    println("    a = $(p.a) m")
    println("    b = $(p.b) m")
    println("    t_repeat = $(p.t_repeat) m")
    println("    Rin = $(p.Rin) m")
    println("    Rout = $(p.Rout) m")
    
    println("\n  节点累计角度范围:")
    println("    内螺旋（正极耳）: [$(round(minimum(θ_cum_in), digits=2)), $(round(maximum(θ_cum_in), digits=2))] 弧度")
    println("                      [$(round(minimum(θ_cum_in)/(2π), digits=2)), $(round(maximum(θ_cum_in)/(2π), digits=2))] 圈")
    println("    外螺旋（负极耳）: [$(round(minimum(θ_cum_out), digits=2)), $(round(maximum(θ_cum_out), digits=2))] 弧度")
    println("                      [$(round(minimum(θ_cum_out)/(2π), digits=2)), $(round(maximum(θ_cum_out)/(2π), digits=2))] 圈")
    
    # ========================================================================
    # 步骤 3：检查极耳节点识别
    # ========================================================================
    println("\n[步骤 3] 识别极耳节点...")
    
    try
        pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
        
        println("  识别结果:")
        println("    正极耳节点数: $(length(pos_idx))")
        println("    负极耳节点数: $(length(neg_idx))")
        println("    总极耳节点数: $(length(unique(vcat(pos_idx, neg_idx))))")
        
        if isempty(pos_idx) && isempty(neg_idx)
            @warn "❌ 未识别到任何极耳节点！"
            println("\n可能的原因:")
            println("  1. theta_neg = 20π 超出了网格的角度覆盖范围")
            println("  2. 极耳宽度 tab.width = $(param_dim.tab.width) m 过小")
            println("\n建议修改:")
            println("  1. 将 theta_neg 改为合理值，例如 theta_neg = [0.0]（与 theta_pos 相同位置）")
            println("  2. 或者检查网格是否覆盖到该角度范围")
            return
        elseif length(unique(vcat(pos_idx, neg_idx))) > nT * 0.5
            @warn "⚠️  极耳节点数量过多（> 50% 总节点数）！"
            println("  这可能导致过度约束，建议减小 tab.width 或调整 theta 值")
        end
        
        # 打印部分极耳节点的坐标（用于验证）
        if !isempty(pos_idx)
            println("\n  正极耳节点示例（前5个）:")
            for (i, idx) in enumerate(pos_idx[1:min(5, length(pos_idx))])
                x, y = mesh_th.node[idx, 1], mesh_th.node[idx, 2]
                r = hypot(x, y)
                θ = atan(y, x)
                println("    节点 $idx: (x=$(round(x, digits=6)), y=$(round(y, digits=6))), r=$(round(r, digits=6)), θ=$(round(θ, digits=3))")
            end
        end
        
        if !isempty(neg_idx)
            println("\n  负极耳节点示例（前5个）:")
            for (i, idx) in enumerate(neg_idx[1:min(5, length(neg_idx))])
                x, y = mesh_th.node[idx, 1], mesh_th.node[idx, 2]
                r = hypot(x, y)
                θ = atan(y, x)
                println("    节点 $idx: (x=$(round(x, digits=6)), y=$(round(y, digits=6))), r=$(round(r, digits=6)), θ=$(round(θ, digits=3))")
            end
        end
        
    catch err
        @error "❌ 极耳节点识别失败！" exception=(err, catch_backtrace())
        return
    end
    
    # ========================================================================
    # 步骤 4：检查 T_ref 和边界条件计算
    # ========================================================================
    println("\n[步骤 4] 检查温度尺度和边界条件计算...")
    
    scale = case.param_dim.scale
    T_amb = case.param_dim.cell.T_amb
    T_ref = scale.T_ref
    
    println("  温度参数:")
    println("    T_ref = $T_ref K")
    println("    T_amb = $T_amb K")
    println("    T0 = $(case.param_dim.cell.T0) K")
    
    if T_ref <= 0 || isnan(T_ref) || isinf(T_ref)
        @error "❌ T_ref 值异常！" value=T_ref
        return
    end
    
    t = 0.0
    rate_Ks = 0.1
    penalty = 1e12
    
    T_amb_nd = T_amb / T_ref
    T_tab_nd = T_amb_nd + (rate_Ks * t) / T_ref
    
    println("\n  边界条件计算（t = 0）:")
    println("    T_amb_nd = T_amb / T_ref = $T_amb / $T_ref = $T_amb_nd")
    println("    T_tab_nd = T_amb_nd + (rate_Ks * t) / T_ref = $T_tab_nd")
    println("    penalty = $penalty")
    
    if isnan(T_amb_nd) || isinf(T_amb_nd) || isnan(T_tab_nd) || isinf(T_tab_nd)
        @error "❌ 边界条件计算结果异常！"
        return
    else
        println("    ✓ 边界条件计算正常")
    end
    
    # ========================================================================
    # 步骤 5：模拟矩阵装配并检查条件数
    # ========================================================================
    println("\n[步骤 5] 模拟热刚度矩阵装配并检查条件数...")
    
    # 创建简化的测试矩阵
    KT_test = spzeros(Float64, nT, nT)
    FT_test = zeros(Float64, nT)
    
    # 添加一些虚拟的扩散项（模拟正常的热传导矩阵）
    for i in 1:min(100, nT)
        KT_test[i, i] = 1.0
        if i > 1
            KT_test[i, i-1] = -0.5
            KT_test[i-1, i] = -0.5
        end
    end
    
    println("  装配前:")
    println("    矩阵维度: $(size(KT_test))")
    println("    非零元素数: $(nnz(KT_test))")
    
    if nnz(KT_test) > 0
        vals = nonzeros(KT_test)
        println("    矩阵元素范围: [$(minimum(vals)), $(maximum(vals))]")
        
        # 尝试计算条件数（仅对小矩阵）
        if nT <= 1000
            K_dense = Matrix(KT_test[1:min(100, nT), 1:min(100, nT)])
            if rank(K_dense) == size(K_dense, 1)
                cond_before = cond(K_dense)
                println("    条件数（装配前，局部）: $(round(cond_before, sigdigits=3))")
            end
        end
    end
    
    # 应用惩罚法（模拟极耳边界条件）
    pos_idx, neg_idx = JuBat.jellyroll_tab_node_indices(mesh_th, param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))
    
    if !isempty(tab_nodes)
        println("\n  应用惩罚法到 $(length(tab_nodes)) 个极耳节点...")
        for n in tab_nodes
            KT_test[n, n] += penalty
            FT_test[n] += penalty * T_tab_nd
        end
        
        println("  装配后:")
        vals_after = nonzeros(KT_test)
        println("    矩阵元素范围: [$(minimum(vals_after)), $(maximum(vals_after))]")
        println("    最大/最小比值: $(maximum(vals_after) / max(minimum(vals_after), 1e-300))")
        
        # 检查矩阵条件数
        if nT <= 1000
            K_dense_after = Matrix(KT_test[1:min(100, nT), 1:min(100, nT)])
            if rank(K_dense_after) == size(K_dense_after, 1)
                cond_after = cond(K_dense_after)
                println("    条件数（装配后，局部）: $(round(cond_after, sigdigits=3))")
                
                if cond_after > 1e10
                    @warn "⚠️  矩阵条件数非常大（> 1e10），可能导致求解不稳定"
                    println("\n建议:")
                    println("  1. 减小惩罚值 penalty（例如从 1e12 改为 1e8 或 1e6）")
                    println("  2. 使用更稳定的求解器（例如迭代求解器）")
                    println("  3. 检查极耳节点是否过多或分布不合理")
                end
            end
        end
        
        # 检查 FT 向量
        println("\n  载荷向量 FT:")
        println("    范围: [$(minimum(FT_test)), $(maximum(FT_test))]")
        println("    极耳节点的 FT 值（示例）: $(FT_test[tab_nodes[1]])")
        
        if any(isnan.(FT_test)) || any(isinf.(FT_test))
            @error "❌ FT 向量包含 NaN 或 Inf！"
        else
            println("    ✓ FT 向量正常")
        end
    end
    
    # ========================================================================
    # 步骤 6：总结和建议
    # ========================================================================
    println("\n" * "="^80)
    println("诊断总结")
    println("="^80)
    
    println("\n✓ T_ref 初始化正常: $T_ref K")
    println("✓ 边界条件计算正常: T_amb_nd=$T_amb_nd, T_tab_nd=$T_tab_nd")
    
    if isempty(tab_nodes)
        println("\n❌ 关键问题：未识别到极耳节点")
        println("   原因：theta_neg = 20π 超出网格角度范围")
        println("\n解决方案：")
        println("  将极耳角度改为网格覆盖范围内的值，例如：")
        println("    param_dim.tab.theta_pos = [0.0]")
        println("    param_dim.tab.theta_neg = [π]  # 或其他合理值")
    else
        println("\n⚠️  潜在问题：矩阵条件数过大")
        println("   原因：惩罚值 penalty = 1e12 过大")
        println("\n解决方案：")
        println("  1. 在 testexample.jl 中添加：")
        println("     opt.tab_penalty = 1e8  # 或 1e6")
        println("  2. 或者减少极耳节点数量（减小 tab.width）")
    end
    
    println("\n" * "="^80)
end

# 运行诊断
diagnose_tab_bc_issue()
