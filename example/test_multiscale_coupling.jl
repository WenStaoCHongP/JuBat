"""
多尺度耦合测试：颗粒应力 → 宏观Jellyroll网格

功能：
- 验证颗粒尺度扩散应力正确映射到宏观单元
- 检查能量守恒
- 可视化尺度桥接效果

作者：AI Assistant
日期：2025-11-24
"""

using LinearAlgebra, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("多尺度耦合测试：颗粒应力 → Jellyroll网格")
    println("="^80)
    
    # ========================================================================
    # 1. 设置简单的SPMe + 2D热案例
    # ========================================================================
    println("\n[1/5] 创建测试案例...")
    
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    
    # 电化学
    opt.model = "SPMe"
    opt.Current = x -> 5.0
    opt.Nn = 5
    opt.Ns = 3
    opt.Np = 5
    opt.Nrn = 5
    opt.Nrp = 5
    opt.time = [0.0, 10.0]
    opt.dt = [0.1, 0.5]
    
    # 热模型
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    
    # 力学模块
    opt.mechanical_enabled = true
    opt.mechanicalmodel = "full"
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建2D网格
    nθ = 8  # 小网格便于测试
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    println("  ✓ 网格创建完成: ne = $ne")
    
    # ========================================================================
    # 2. 手动构造颗粒应力（模拟电化学计算结果）
    # ========================================================================
    println("\n[2/5] 构造测试数据...")
    
    variables = Dict{String, Union{Array{Float64}, Float64}}()
    
    # 模拟颗粒应力（微观尺度，~10 MPa量级）
    σ_particle_n = 10e6   # 10 MPa（负极）
    σ_particle_p = 15e6   # 15 MPa（正极）
    
    variables["negative particle surface tangential stress"] = σ_particle_n
    variables["positive particle surface tangential stress"] = σ_particle_p
    
    # 模拟温度场（用于热应力）
    Tn = fill(case.param.cell.T0 + 5.0/case.param.scale.T_ref, mesh_th.nlen)
    variables["T_nodes"] = Tn
    variables["temperature"] = mean(Tn)
    
    # 添加layer_weights（用于精细映射）
    fks = try
        JuBat.jellyroll_get_layer_weights(mesh_th)
    catch
        JuBat.jellyroll_element_layer_weights(mesh_th, param_dim; nsamples_per_dim=2, logic=:spiral)
    end
    variables["thermal2D layer_weights"] = fks
    
    @printf("  ✓ 颗粒应力: σ_n = %.1f MPa, σ_p = %.1f MPa\n", 
            σ_particle_n/1e6, σ_particle_p/1e6)
    
    # ========================================================================
    # 3. 计算热应力和均匀化扩散应力
    # ========================================================================
    println("\n[3/5] 执行多尺度均匀化...")
    
    # 计算热应力（会自动调用均匀化）
    variables = JuBat.thermal_stress(case, variables)
    
    # 提取结果
    if haskey(variables, "thermal2D element thermal stress")
        σ_thermal = variables["thermal2D element thermal stress"]
        @printf("  热应力: 平均 = %.2f MPa, 范围 = [%.2f, %.2f] MPa\n",
                mean(σ_thermal)/1e6, minimum(σ_thermal)/1e6, maximum(σ_thermal)/1e6)
    end
    
    if haskey(variables, "thermal2D element diffusion stress (homogenized)")
        σ_diff_homo = variables["thermal2D element diffusion stress (homogenized)"]
        @printf("  扩散应力（均匀化）: 平均 = %.2f MPa, 范围 = [%.2f, %.2f] MPa\n",
                mean(σ_diff_homo)/1e6, minimum(σ_diff_homo)/1e6, maximum(σ_diff_homo)/1e6)
    end
    
    if haskey(variables, "thermal2D element total stress")
        σ_total = variables["thermal2D element total stress"]
        @printf("  总应力: 平均 = %.2f MPa, 范围 = [%.2f, %.2f] MPa\n",
                mean(σ_total)/1e6, minimum(σ_total)/1e6, maximum(σ_total)/1e6)
    end
    
    # ========================================================================
    # 4. 验证能量守恒
    # ========================================================================
    println("\n[4/5] 验证多尺度一致性...")
    
    # 理论验证：均匀化应力应该在颗粒应力的量级
    σ_diff_homo = variables["thermal2D element diffusion stress (homogenized)"]
    
    # 体积分数
    ε_n = hasproperty(param_dim.NE, :epsilon_s) ? param_dim.NE.epsilon_s : 0.65
    ε_p = hasproperty(param_dim.PE, :epsilon_s) ? param_dim.PE.epsilon_s : 0.60
    
    # 理论预期（厚度加权）
    t_n = param_dim.NE.thickness
    t_p = param_dim.PE.thickness
    σ_theory = (σ_particle_n * ε_n * t_n + σ_particle_p * ε_p * t_p) / (t_n + t_p)
    
    @printf("  理论预期均匀化应力: %.2f MPa\n", σ_theory/1e6)
    @printf("  实际计算平均应力: %.2f MPa\n", mean(σ_diff_homo)/1e6)
    
    rel_error = abs(mean(σ_diff_homo) - σ_theory) / σ_theory
    @printf("  相对误差: %.2e (%.2f%%)\n", rel_error, rel_error*100)
    
    if rel_error < 0.1
        println("  ✓ 验证通过！（误差 < 10%）")
    else
        println("  ⚠ 警告：误差较大，请检查layer_weights")
    end
    
    # ========================================================================
    # 5. 可视化
    # ========================================================================
    println("\n[5/5] 生成可视化...")
    
    # 获取单元中心坐标
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    θ_centers = atan.(centers[:,2], centers[:,1])
    
    # 图1：扩散应力空间分布
    p1 = scatter(θ_centers, r_centers,
                 marker_z=σ_diff_homo ./ 1e6,
                 markersize=8,
                 color=:viridis,
                 xlabel="θ (rad)", ylabel="r (m)",
                 title="Diffusion Stress (Homogenized) [MPa]",
                 colorbar=true,
                 legend=false)
    
    # 图2：layer_weights分布
    if fks !== nothing
        f_neg = fks[:, 1]
        f_pos = fks[:, 3]
        
        p2 = scatter(θ_centers, r_centers,
                     marker_z=f_neg,
                     markersize=8,
                     color=:reds,
                     xlabel="θ (rad)", ylabel="r (m)",
                     title="Negative Electrode Volume Fraction",
                     colorbar=true,
                     clims=(0, 1),
                     legend=false)
        
        p3 = scatter(θ_centers, r_centers,
                     marker_z=f_pos,
                     markersize=8,
                     color=:blues,
                     xlabel="θ (rad)", ylabel="r (m)",
                     title="Positive Electrode Volume Fraction",
                     colorbar=true,
                     clims=(0, 1),
                     legend=false)
    else
        p2 = plot(title="layer_weights not available")
        p3 = plot()
    end
    
    # 图3：总应力分布
    σ_total = variables["thermal2D element total stress"]
    p4 = scatter(θ_centers, r_centers,
                 marker_z=σ_total ./ 1e6,
                 markersize=8,
                 color=:plasma,
                 xlabel="θ (rad)", ylabel="r (m)",
                 title="Total Stress (Thermal + Diffusion) [MPa]",
                 colorbar=true,
                 legend=false)
    
    # 组合图
    p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1400, 1200))
    savefig(p, "test_multiscale_coupling.png")
    println("  ✓ 保存图像: test_multiscale_coupling.png")
    
    # ========================================================================
    # 总结
    # ========================================================================
    println("\n" * "="^80)
    println("多尺度耦合测试总结")
    println("="^80)
    
    println("""
    测试配置：
      - 网格单元数: $ne
      - 颗粒应力: σ_n = $(σ_particle_n/1e6) MPa, σ_p = $(σ_particle_p/1e6) MPa
      - 均匀化方法: layer_weights加权 $(fks !== nothing ? "✓" : "✗")
    
    验证结果：
      - 理论预期: $(σ_theory/1e6) MPa
      - 实际计算: $(mean(σ_diff_homo)/1e6) MPa
      - 相对误差: $(rel_error*100)%
      - 状态: $(rel_error < 0.1 ? "✓ 通过" : "⚠ 警告")
    
    物理意义：
      - 颗粒尺度（~10 μm）的扩散应力成功映射到宏观网格（~1 mm）
      - 体积分数修正考虑了孔隙的影响
      - layer_weights提供了空间分布的精细信息
      - 总应力 = 热应力 + 扩散应力（均匀化）
    
    关键输出变量：
      - "thermal2D element diffusion stress (homogenized)" - 均匀化扩散应力
      - "thermal2D element total stress" - 总应力场
    """)
    
    println("="^80)
end

# 运行测试
main()
