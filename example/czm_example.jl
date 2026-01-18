# czm_example.jl - 内聚力模型示例
# 
# 演示功能：
# 1. 基于热网格创建内聚力网格
# 2. 施加载荷进行非线性分析
# 3. 观察损伤演化过程
# 4. 多步加载模拟循环损伤累积

using Pkg
Pkg.activate(".")

using LinearAlgebra
using SparseArrays
using Plots
using Statistics

# 加载JuBat模块
include("../src/JuBat.jl")
using .JuBat

println("="^60)
println("内聚力模型 (CZM) 示例")
println("="^60)

# ========================================================================
# 1. 参数设置
# ========================================================================

println("\n[1] 加载参数...")

# 使用 Jellyroll 参数（包含内聚力参数）
param_dim = ChooseCell("Jellyroll")
param = NormaliseParam(param_dim)

# 打印内聚力参数
println("\n内聚力参数（有量纲）：")
println("  法向最大牵引力 σ_max_n = $(param_dim.cohesive.σ_max_n / 1e6) MPa")
println("  法向临界位移 δ_c_n = $(param_dim.cohesive.δ_c_n * 1e6) μm")
println("  法向断裂能 G_c_n = $(param_dim.cohesive.G_c_n) J/m²")
println("  切向最大牵引力 τ_max_t = $(param_dim.cohesive.τ_max_t / 1e6) MPa")
println("  切向临界位移 δ_c_t = $(param_dim.cohesive.δ_c_t * 1e6) μm")
println("  BK指数 η = $(param_dim.cohesive.eta)")

# ========================================================================
# 2. 创建热网格
# ========================================================================

println("\n[2] 创建热网格...")

# 使用较少的单元数便于测试
nθ = 60  # 每圈60个单元

thermal_mesh = jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)

println("  热网格信息：")
println("    节点数：$(thermal_mesh.nlen)")
println("    单元数：$(size(thermal_mesh.element, 1))")
println("    单元类型：$(thermal_mesh.type)")

# ========================================================================
# 3. 创建内聚力网格
# ========================================================================

println("\n[3] 创建内聚力网格...")

czm_mesh = create_czm_mesh(thermal_mesh, param_dim; nθ_per_turn=nθ)

println("  内聚力网格信息：")
println("    扩展后节点数：$(czm_mesh.nnode)")
println("    固体单元数：$(size(czm_mesh.bulk_element, 1))")
println("    内聚力单元数：$(czm_mesh.n_cohesive)")
println("    层数：$(czm_mesh.n_layers)")
println("    层间界面数：$(czm_mesh.n_layers - 1)")

if czm_mesh.n_cohesive == 0
    println("\n警告：没有创建内聚力单元。可能需要调整网格参数。")
    println("当前热网格可能只有一层。")
end

# ========================================================================
# 4. 设置边界条件和载荷
# ========================================================================

println("\n[4] 设置边界条件和载荷...")

nnode = czm_mesh.nnode
ndof = 2 * nnode

# 外力向量（初始为零）
F_ext = zeros(Float64, ndof)

# 材料参数（使用平均值）
E_eff = 0.5 * (param_dim.NE.E + param_dim.PE.E)
ν_eff = 0.5 * (param_dim.NE.nu + param_dim.PE.nu)

println("  有效杨氏模量 E = $(E_eff / 1e9) GPa")
println("  有效泊松比 ν = $(ν_eff)")

# ========================================================================
# 5. 单步分析测试
# ========================================================================

println("\n[5] 单步分析测试...")

if czm_mesh.n_cohesive > 0
    # 施加小载荷进行测试
    # 在外边界节点施加径向位移载荷
    
    # 识别外边界节点（简化：使用半径最大的节点）
    r_nodes = [hypot(czm_mesh.node[i, 1], czm_mesh.node[i, 2]) for i in 1:nnode]
    r_max = maximum(r_nodes)
    outer_nodes = findall(r -> r > 0.95 * r_max, r_nodes)
    
    # 施加径向拉力
    load_magnitude = 1e6  # 1 MPa 等效
    for n in outer_nodes
        x, y = czm_mesh.node[n, 1], czm_mesh.node[n, 2]
        r = hypot(x, y)
        if r > 1e-10
            # 径向单位向量
            nx, ny = x / r, y / r
            F_ext[2*n - 1] += load_magnitude * nx
            F_ext[2*n] += load_magnitude * ny
        end
    end
    
    println("  施加载荷：在外边界 $(length(outer_nodes)) 个节点上施加径向拉力")
    
    # 执行牛顿-拉弗森求解
    result = newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, 
                                param_dim.cohesive, param_dim;
                                max_iter=30, tol=1e-6)
    
    println("\n  求解结果：")
    println("    收敛：$(result.converged)")
    println("    迭代次数：$(result.iterations)")
    println("    残差范数：$(result.residual_norm)")
    
    # 损伤统计
    stats = get_damage_statistics(czm_mesh)
    println("\n  损伤统计：")
    println("    最大损伤：$(round(stats.max_D * 100, digits=2))%")
    println("    平均损伤：$(round(stats.mean_D * 100, digits=2))%")
    println("    断裂单元数：$(stats.n_fractured)")
    
    # 位移统计
    u_mag = [hypot(result.displacement[2*i-1], result.displacement[2*i]) for i in 1:nnode]
    println("\n  位移统计：")
    println("    最大位移：$(maximum(u_mag) * 1e6) μm")
    println("    平均位移：$(mean(u_mag) * 1e6) μm")
else
    println("  跳过（无内聚力单元）")
end

# ========================================================================
# 6. 多步加载模拟
# ========================================================================

println("\n[6] 多步加载模拟...")

if czm_mesh.n_cohesive > 0
    # 重置损伤状态
    reset_damage_states!(czm_mesh)
    
    # 加载步数
    n_steps = 10
    load_factors = range(0.1, 2.0, length=n_steps)
    
    # 存储历史
    damage_history = zeros(Float64, n_steps)
    displacement_history = zeros(Float64, n_steps)
    
    u_prev = zeros(Float64, ndof)
    
    println("  执行 $(n_steps) 个加载步...")
    
    for (i, λ) in enumerate(load_factors)
        F_step = λ * F_ext
        
        result = solve_czm_step(czm_mesh, F_step, E_eff, ν_eff,
                               param_dim.cohesive, param_dim, u_prev;
                               max_iter=30, tol=1e-6)
        
        if result.converged
            u_prev = result.displacement
        end
        
        stats = get_damage_statistics(czm_mesh)
        damage_history[i] = stats.max_D
        
        u_mag = [hypot(result.displacement[2*j-1], result.displacement[2*j]) for j in 1:nnode]
        displacement_history[i] = maximum(u_mag)
        
        print("    步 $i: λ=$(round(λ, digits=2)), D_max=$(round(stats.max_D*100, digits=1))%, ")
        println("u_max=$(round(maximum(u_mag)*1e6, digits=2)) μm")
        
        # 检查断裂
        is_fractured, info = check_fracture_criterion(czm_mesh; threshold=0.99)
        if is_fractured
            println("    *** 检测到宏观断裂！***")
            break
        end
    end
    
    # 绘制损伤演化曲线
    println("\n  绘制损伤演化曲线...")
    
    p1 = plot(load_factors[1:length(damage_history)], damage_history * 100,
              xlabel="载荷因子 λ", ylabel="最大损伤 (%)",
              title="损伤演化", legend=false, lw=2, marker=:circle)
    
    p2 = plot(load_factors[1:length(displacement_history)], displacement_history * 1e6,
              xlabel="载荷因子 λ", ylabel="最大位移 (μm)",
              title="位移响应", legend=false, lw=2, marker=:circle)
    
    p = plot(p1, p2, layout=(1, 2), size=(800, 350))
    
    savefig(p, "output/czm_damage_evolution.png")
    println("  图像已保存到 output/czm_damage_evolution.png")
else
    println("  跳过（无内聚力单元）")
end

# ========================================================================
# 7. 双线性本构测试
# ========================================================================

println("\n[7] 双线性本构曲线测试...")

# 创建一个测试用的损伤状态
test_state = DamageState()
cohesive_params = param_dim.cohesive

# 生成分离位移范围
δ_range = range(0, 2 * cohesive_params.δ_c_n, length=100)
T_n_curve = Float64[]
D_curve = Float64[]

for δ in δ_range
    T_n, T_t, D = bilinear_traction(δ, 0.0, test_state, cohesive_params; update=true)
    push!(T_n_curve, T_n)
    push!(D_curve, D)
end

# 绘制牵引力-分离曲线
p3 = plot(δ_range * 1e6, T_n_curve / 1e6,
          xlabel="法向分离位移 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="双线性牵引力-分离曲线", legend=false, lw=2)

# 标记关键点
vline!([cohesive_params.δ_0_n * 1e6], linestyle=:dash, color=:gray, label="δ_0")
vline!([cohesive_params.δ_c_n * 1e6], linestyle=:dash, color=:red, label="δ_c")

p4 = plot(δ_range * 1e6, D_curve * 100,
          xlabel="法向分离位移 δ_n (μm)", ylabel="损伤变量 D (%)",
          title="损伤演化曲线", legend=false, lw=2, color=:red)

p_constitutive = plot(p3, p4, layout=(1, 2), size=(800, 350))
savefig(p_constitutive, "output/czm_constitutive.png")
println("  图像已保存到 output/czm_constitutive.png")

# ========================================================================
# 8. 加卸载测试
# ========================================================================

println("\n[8] 加卸载行为测试...")

# 重置损伤状态
test_state_cycle = DamageState()

# 定义加卸载历程：加载到某点，卸载，再加载
δ_history = vcat(
    range(0, 0.6 * cohesive_params.δ_c_n, length=20),  # 加载
    range(0.6 * cohesive_params.δ_c_n, 0.2 * cohesive_params.δ_c_n, length=10),  # 卸载
    range(0.2 * cohesive_params.δ_c_n, 0.8 * cohesive_params.δ_c_n, length=15),  # 再加载
    range(0.8 * cohesive_params.δ_c_n, 0, length=10)   # 完全卸载
)

T_cycle = Float64[]
D_cycle = Float64[]

for δ in δ_history
    T_n, _, D = bilinear_traction(δ, 0.0, test_state_cycle, cohesive_params; update=true)
    push!(T_cycle, T_n)
    push!(D_cycle, D)
end

p5 = plot(collect(δ_history) * 1e6, T_cycle / 1e6,
          xlabel="法向分离位移 δ_n (μm)", ylabel="法向牵引力 T_n (MPa)",
          title="加卸载行为", legend=false, lw=2)

# 添加箭头表示加载方向
scatter!([δ_history[20] * 1e6], [T_cycle[20] / 1e6], marker=:circle, color=:red, label="卸载起点")
scatter!([δ_history[30] * 1e6], [T_cycle[30] / 1e6], marker=:diamond, color=:blue, label="再加载起点")

savefig(p5, "output/czm_loading_unloading.png")
println("  图像已保存到 output/czm_loading_unloading.png")

# ========================================================================
# 总结
# ========================================================================

println("\n" * "="^60)
println("CZM 模块功能验证完成！")
println("="^60)

println("\n生成的文件：")
println("  - output/czm_damage_evolution.png：损伤演化曲线")
println("  - output/czm_constitutive.png：本构关系曲线")
println("  - output/czm_loading_unloading.png：加卸载行为曲线")

println("\n模块功能清单：")
println("  ✓ 数据结构：CohesiveElement, CohesiveMesh, DamageState, CZMResult")
println("  ✓ 网格划分：基于热网格的分层内聚力网格生成")
println("  ✓ 双线性本构：牵引力-分离关系、混合模式（BK准则）")
println("  ✓ 损伤演化：加卸载准则、历史依赖性")
println("  ✓ 单元计算：形函数、分离位移、切线刚度")
println("  ✓ 系统组装：固体-内聚力耦合")
println("  ✓ 牛顿-拉弗森求解：非线性迭代")
println("  ✓ 损伤统计与断裂判据")

println("\n下一步开发建议：")
println("  - 与电化学-热耦合模型集成")
println("  - 实现充放电循环下的损伤累积")
println("  - 添加温度依赖的内聚力参数")
println("  - 优化大规模网格的求解效率")
