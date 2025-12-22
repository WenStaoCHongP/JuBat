"""
测试力学边界条件识别

此脚本验证内外螺旋边界节点的识别是否正确
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf

# 强制重新加载模块
if isdefined(Main, :JuBat)
    println("检测到已加载的 JuBat 模块，强制重新加载...")
    @eval Main JuBat = nothing
    GC.gc()
end

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function test_mechanical_bc()
    println("="^80)
    println("力学边界条件识别测试")
    println("="^80)
    
    # 创建简单的 Jellyroll 案例
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建网格
    nθ = 40
    println("\n创建 Jellyroll 网格 (nθ=$nθ)...")
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("  总单元数: $ne")
    println("  总节点数: $nT")
    
    # 测试边界节点识别
    println("\n识别边界节点...")
    
    # 直接调用边界识别逻辑（复制自 mechanical.jl）
    nnode = mesh_th.nlen
    bc_nodes = Dict{Int, Symbol}()
    
    # 获取螺旋参数
    pgeo = JuBat.jellyroll_spiral_params(param_dim)
    s_in = 0.0
    s_out = pgeo.t_repeat
    bval = max(pgeo.b, 1e-12)
    θ0_mesh = max(0.0, (pgeo.Rin - pgeo.a - s_in) / bval)
    θ1_mesh = min((pgeo.Rout - pgeo.a - s_out) / bval, (pgeo.Rout - pgeo.a) / bval)
    
    # 角度范围
    θ_in_range = (θ0_mesh, min(θ0_mesh + 2.0*π, θ1_mesh))
    θ_out_range = (max(θ1_mesh - 2.0*π, θ0_mesh), θ1_mesh)
    tol = 1e-4
    
    println("  内边界角度范围: [$(θ_in_range[1]), $(θ_in_range[2])] rad")
    println("  外边界角度范围: [$(θ_out_range[1]), $(θ_out_range[2])] rad")
    
    # 识别内边界节点（第一圈）
    inner_count = 0
    for i in 1:nnode
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:inner, theta_range=θ_in_range, tol=tol)
            bc_nodes[i] = :fixed_xy
            inner_count += 1
        end
    end
    
    # 识别外边界节点（最后一圈）
    outer_count = 0
    for i in 1:nnode
        if JuBat.edge_boundary(mesh_th, i, param_dim; which=:outer, theta_range=θ_out_range, tol=tol)
            bc_nodes[i] = :fixed_xy
            outer_count += 1
        end
    end
    
    println("  [力学边界条件] 内边界固定节点: $inner_count, 外边界固定节点: $outer_count")
    
    n_bc = length(bc_nodes)
    println("\n边界条件统计:")
    println("  固定节点总数: $n_bc")
    
    # 统计不同类型的边界条件
    n_fixed_xy = count(v -> v == :fixed_xy, values(bc_nodes))
    n_fixed_x = count(v -> v == :fixed_x, values(bc_nodes))
    n_fixed_y = count(v -> v == :fixed_y, values(bc_nodes))
    
    println("  固定 x,y: $n_fixed_xy")
    println("  仅固定 x: $n_fixed_x")
    println("  仅固定 y: $n_fixed_y")
    
    # 可视化边界节点
    println("\n生成边界节点可视化...")
    
    # 获取所有节点坐标
    x_all = mesh_th.node[:, 1]
    y_all = mesh_th.node[:, 2]
    
    # 获取边界节点坐标
    bc_indices = collect(keys(bc_nodes))
    x_bc = mesh_th.node[bc_indices, 1]
    y_bc = mesh_th.node[bc_indices, 2]
    
    # 绘图
    p = scatter(x_all, y_all, 
                markersize=2, 
                alpha=0.3,
                color=:lightgray,
                label="所有节点",
                aspect_ratio=:equal,
                xlabel="x [m]",
                ylabel="y [m]",
                title="力学边界条件：内外螺旋边界固定")
    
    scatter!(p, x_bc, y_bc,
             markersize=4,
             color=:red,
             label="固定节点 (内外边界)",
             markershape=:circle)
    
    # 添加参考圆
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    θ_ref = range(0, 2π, length=100)
    
    plot!(p, Rin .* cos.(θ_ref), Rin .* sin.(θ_ref),
          color=:blue, linestyle=:dash, label="内半径", linewidth=2)
    plot!(p, Rout .* cos.(θ_ref), Rout .* sin.(θ_ref),
          color=:green, linestyle=:dash, label="外半径", linewidth=2)
    
    output_dir = joinpath(@__DIR__, "..", "output")
    isdir(output_dir) || mkpath(output_dir)
    output_file = joinpath(output_dir, "mechanical_bc_nodes.png")
    savefig(p, output_file)
    println("  ✓ 保存: $output_file")
    
    # 统计边界节点的径向分布
    println("\n边界节点径向分布:")
    r_bc = hypot.(x_bc, y_bc)
    println("  半径范围: [$(minimum(r_bc)), $(maximum(r_bc))] m")
    println("  平均半径: $(mean(r_bc)) m")
    println("  Rin: $Rin m")
    println("  Rout: $Rout m")
    
    # 检查是否有节点接近内外边界
    tol = 1e-3
    n_near_inner = count(r -> abs(r - Rin) < tol, r_bc)
    n_near_outer = count(r -> abs(r - Rout) < tol, r_bc)
    println("\n节点分布验证 (容差 $tol m):")
    println("  接近内半径: $n_near_inner")
    println("  接近外半径: $n_near_outer")
    
    if n_near_inner > 0 && n_near_outer > 0
        println("\n✅ 测试通过：成功识别内外边界节点")
    else
        println("\n⚠️  警告：边界节点识别可能不正确")
    end
    
    println("\n" * "="^80)
    println("测试完成")
    println("="^80)
    
    return bc_nodes
end

# 运行测试
test_mechanical_bc()
