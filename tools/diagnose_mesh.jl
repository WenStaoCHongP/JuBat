using Statistics
include("../src/JuBat.jl")

println("=" ^ 80)
println("网格结构诊断")
println("=" ^ 80)

param_dim = JuBat.ChooseCell("Jellyroll")

# 测试不同的 gsorder
for gsorder in [1, 2]
    println("\n" * "-" ^ 80)
    println("gsorder = $gsorder")
    println("-" ^ 80)
    
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=gsorder)
    
    println("网格基本信息:")
    println("  总节点数: $(mesh_th.nlen)")
    println("  总单元数: $(size(mesh_th.element, 1))")
    println("  单元类型: $(mesh_th.type)")
    println("  维度: $(mesh_th.dimension)")
    
    # 检查单元连接性
    ne = size(mesh_th.element, 1)
    nodes_per_elem = size(mesh_th.element, 2)
    println("  每个单元的节点数: $nodes_per_elem")
    
    # 检查前几个单元的连接
    println("\n前3个单元的节点连接:")
    for i in 1:min(3, ne)
        println("    单元 $i: $(mesh_th.element[i,:])")
    end
    
    # 检查节点分布
    x_all = mesh_th.node[:,1]
    y_all = mesh_th.node[:,2]
    r_all = hypot.(x_all, y_all)
    
    println("\n节点分布:")
    println("  半径范围: [$(minimum(r_all)), $(maximum(r_all))] m")
    
    # 尝试识别边界节点（使用默认 θ 范围）
    p = JuBat.jellyroll_spiral_params(param_dim)
    
    # 识别第0圈内螺旋
    is_inner = [JuBat.edge_boundary(mesh_th, i, param_dim; 
                                     which=:inner, theta_range=(0.0, 2π), tol=1e-4) 
                for i in 1:mesh_th.nlen]
    count_inner = count(identity, is_inner)
    
    # 识别第0圈外螺旋
    is_outer = [JuBat.edge_boundary(mesh_th, i, param_dim; 
                                     which=:outer, theta_range=(0.0, 2π), tol=1e-4) 
                for i in 1:mesh_th.nlen]
    count_outer = count(identity, is_outer)
    
    println("\n边界识别结果（tol=1e-4）:")
    println("  内螺旋节点数: $count_inner")
    println("  外螺旋节点数: $count_outer")
    
    # 尝试不同的容差
    if gsorder == 2
        println("\n尝试不同的容差:")
        for tol in [1e-3, 1e-2, 1e-1]
            is_outer_tol = [JuBat.edge_boundary(mesh_th, i, param_dim; 
                                                which=:outer, theta_range=(0.0, 2π), tol=tol) 
                           for i in 1:mesh_th.nlen]
            count_outer_tol = count(identity, is_outer_tol)
            println("    tol=$tol: $count_outer_tol 个外螺旋节点")
        end
    end
end

println("\n" * "=" ^ 80)
println("结论")
println("=" ^ 80)

println("\n如果 gsorder=2 导致节点数增加，可能的原因:")
println("  1. 高阶单元（Q8/Q9）：每条边上有中间节点")
println("  2. 网格细化：自动插入额外节点")
println("  3. 高斯积分点被存储为节点")

println("\n建议:")
println("  1. 使用 gsorder=1（线性单元）")
println("  2. 或者调整容差以适应高阶单元的节点分布")
println("  3. 或者修改边界识别逻辑以适配高阶单元")
