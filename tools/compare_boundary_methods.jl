using Statistics
include("../src/JuBat.jl")

param_dim = JuBat.ChooseCell("Jellyroll")
mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)
p = JuBat.jellyroll_spiral_params(param_dim)

println("=" ^ 80)
println("边界识别方法对比：collector_seed_mesh 网格（nθ=160）")
println("=" ^ 80)

println("\n网格信息:")
println("  总节点数: $(mesh_th.nlen)")
println("  Rin = $(p.Rin) m, Rout = $(p.Rout) m")
println("  螺距 b = $(p.b) m/rad")
println("  层厚 t_repeat = $(p.t_repeat) m")
println("  总圈数 N = $(p.n_wind)")

# 计算节点半径
x_all = mesh_th.node[:,1]
y_all = mesh_th.node[:,2]
r_all = hypot.(x_all, y_all)

println("\n节点半径分布:")
println("  最小半径: $(minimum(r_all)) m")
println("  最大半径: $(maximum(r_all)) m")

# 方法 1：基于螺旋线的边界识别（原方法，刚修复）
println("\n" * "=" ^ 80)
println("方法 1：基于螺旋线 θ 范围的边界识别 (:node_on)")
println("=" ^ 80)

is_inner_spiral = [JuBat.edge_boundary(:node_on, mesh_th, i, param_dim; which=:inner) for i in 1:mesh_th.nlen]
is_outer_spiral = [JuBat.edge_boundary(:node_on, mesh_th, i, param_dim; which=:outer) for i in 1:mesh_th.nlen]

count_inner_spiral = count(identity, is_inner_spiral)
count_outer_spiral = count(identity, is_outer_spiral)

println("  内螺旋线节点数: $count_inner_spiral")
println("  外螺旋线节点数: $count_outer_spiral")

# 检查重叠
overlap_spiral = count(i -> is_inner_spiral[i] && is_outer_spiral[i], 1:mesh_th.nlen)
println("  重叠节点数: $overlap_spiral")

# 半径范围
if count_inner_spiral > 0
    r_inner_spiral = r_all[is_inner_spiral]
    println("  内螺旋线节点半径范围: [$(minimum(r_inner_spiral)), $(maximum(r_inner_spiral))] m")
end
if count_outer_spiral > 0
    r_outer_spiral = r_all[is_outer_spiral]
    println("  外螺旋线节点半径范围: [$(minimum(r_outer_spiral)), $(maximum(r_outer_spiral))] m")
end

println("\n  说明：")
println("    - 内螺旋线：θ_cum_in ∈ [0, 2π]（第1圈）")
println("    - 外螺旋线：mod(θ_cum_out, 2π) ∈ [0, 2π]（任意一圈）")
println("    - 对于 collector_seed_mesh：内外螺旋在同一 θ 段，只是半径不同")

# 方法 2：基于半径的边界识别（新方法）
println("\n" * "=" ^ 80)
println("方法 2：基于半径的边界识别 (:radial)")
println("=" ^ 80)

# 使用默认容差 1e-4
is_inner_radial = [JuBat.edge_boundary(:radial, mesh_th, i, param_dim; which=:inner, tol=1e-4) for i in 1:mesh_th.nlen]
is_outer_radial = [JuBat.edge_boundary(:radial, mesh_th, i, param_dim; which=:outer, tol=1e-4) for i in 1:mesh_th.nlen]

count_inner_radial = count(identity, is_inner_radial)
count_outer_radial = count(identity, is_outer_radial)

println("  内边界节点数 (|r - Rin| < 1e-4): $count_inner_radial")
println("  外边界节点数 (|r - Rout| < 1e-4): $count_outer_radial")

# 检查重叠
overlap_radial = count(i -> is_inner_radial[i] && is_outer_radial[i], 1:mesh_th.nlen)
println("  重叠节点数: $overlap_radial")

# 半径范围
if count_inner_radial > 0
    r_inner_radial = r_all[is_inner_radial]
    println("  内边界节点半径范围: [$(minimum(r_inner_radial)), $(maximum(r_inner_radial))] m")
end
if count_outer_radial > 0
    r_outer_radial = r_all[is_outer_radial]
    println("  外边界节点半径范围: [$(minimum(r_outer_radial)), $(maximum(r_outer_radial))] m")
end

println("\n  说明：")
println("    - 内边界：r ≈ Rin（靠近内半径）")
println("    - 外边界：r ≈ Rout（靠近外半径）")
println("    - 与 θ 范围无关，纯粹基于半径距离")

# 尝试更宽松的容差
println("\n  尝试更宽松的容差 (tol = 1e-3):")
is_inner_radial_loose = [JuBat.edge_boundary(:radial, mesh_th, i, param_dim; which=:inner, tol=1e-3) for i in 1:mesh_th.nlen]
is_outer_radial_loose = [JuBat.edge_boundary(:radial, mesh_th, i, param_dim; which=:outer, tol=1e-3) for i in 1:mesh_th.nlen]
println("    内边界节点数: $(count(identity, is_inner_radial_loose))")
println("    外边界节点数: $(count(identity, is_outer_radial_loose))")

# 对比总结
println("\n" * "=" ^ 80)
println("总结与建议")
println("=" ^ 80)

println("\n1. collector_seed_mesh 网格的结构特点:")
println("   - 内圈节点（1-161）：r_in = a + b*θ，θ ∈ [0, 2π]")
println("   - 外圈节点（162-322）：r_out = a + b*θ + t_repeat，θ ∈ [0, 2π]（相同的 θ！）")
println("   - 网格只覆盖一个 θ 段，无法区分"第1圈"和"第N圈"")

println("\n2. 您的需求理解:")
println("   如果您希望：")
println("   - 内边界 = 第1圈的内螺旋（θ ∈ [0, 2π]）")
println("   - 外边界 = 第N圈的外螺旋（θ ∈ [2π(N-1), 2π*N]）")
println("   → collector_seed_mesh 无法实现（网格本身不覆盖第N圈）")

println("\n3. 推荐方案:")
println("   使用方法2（基于半径）：")
println("   - 内边界：r ≈ Rin（物理上的最内层）")
println("   - 外边界：r ≈ Rout（物理上的最外层）")
println("   - 适用于任何网格类型（collector_seed、inscribed等）")
println("   - 容差可调（建议 tol = 1e-3 ~ 1e-4）")

println("\n4. 使用示例:")
println("   # 施加内边界条件")
println("   for i in 1:mesh.nlen")
println("       if edge_boundary(:radial, mesh, i, param_dim; which=:inner, tol=1e-3)")
println("           # 内边界节点，施加边界条件")
println("       end")
println("   end")
println("\n   # 施加外边界条件（对流换热等）")
println("   for i in 1:mesh.nlen")
println("       if edge_boundary(:radial, mesh, i, param_dim; which=:outer, tol=1e-3)")
println("           # 外边界节点，施加对流换热")
println("       end")
println("   end")

println("\n" * "=" ^ 80)
