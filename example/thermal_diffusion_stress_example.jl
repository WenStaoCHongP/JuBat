# ========================================================================
# 宏观热-扩散应力计算示例
# ========================================================================
#
# 本示例展示如何使用 JuBat 计算宏观层面的热应力和扩散应力
# 
# 理论基础：
# 1. 热应力：由温度变化引起，ε_th = α * ΔT
# 2. 扩散应力：由锂浓度（SOC）变化引起，ε_diff = β_c * c_s_max * ΔSOC
# 3. 总初始应变：ε_0 = ε_th + ε_diff
# 4. 应力平衡方程：∇·σ = 0
# 5. 本构关系：σ = D * (ε - ε_0)
#
# 输出：
# - 2D应力场：σ_xx, σ_yy, σ_xy
# - Von Mises等效应力
# - 位移场：u_x, u_y
# ========================================================================

using Plots
using Printf

# 包含 JuBat 主模块
include("../src/JuBat.jl")

# ========================================================================
# 1. 参数设置
# ========================================================================

println("="^70)
println("宏观热-扩散应力计算示例")
println("="^70)

# 选择电池参数
param_dim = JuBat.ChooseCell("LG M50")

# 创建选项
opt = JuBat.Option()
opt.model = "SPMe"  # 使用 SPMe 模型

# 电流设置（1C 放电）
Crate = 1.0
I_app = Crate * param_dim.cell.I1C
opt.Current = t -> I_app

# 时间设置
opt.time = [0.0, 600.0]  # 10分钟模拟
opt.dt = [1.0, 5.0]  # 自适应时间步长

# 启用热模型（可选，用于生成温度场）
opt.thermalmodel = "lumped"  # 简化起见先用集总热模型

# 启用力学模型
opt.mechanicalmodel = "full"

println("\n模型设置:")
println("  电池型号: LG M50")
println("  电化学模型: $(opt.model)")
println("  热模型: $(opt.thermalmodel)")
println("  力学模型: $(opt.mechanicalmodel)")
println("  放电倍率: $(Crate)C")
println("  模拟时间: $(opt.time[2]) s")

# ========================================================================
# 2. 创建案例并求解电化学模型
# ========================================================================

println("\n" * "="^70)
println("步骤 1: 求解电化学模型")
println("="^70)

case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)

println("✓ 电化学模型求解完成")
println("  时间步数: $(length(result["time [s]"]))")
println("  最终电压: $(result["cell voltage[V]"][end]) V")
println("  最终温度: $(result["temperature[K]"][end]) K")

# ========================================================================
# 3. 创建2D热网格（果冻卷结构）
# ========================================================================

println("\n" * "="^70)
println("步骤 2: 创建2D热-力网格")
println("="^70)

# 检查是否有果冻卷参数
if !hasproperty(param_dim.cell, :Rin) || !hasproperty(param_dim.cell, :Rout)
    println("\n警告: LG M50 不是果冻卷电池，创建简化的2D网格...")
    
    # 对于非果冻卷电池，创建简化的矩形网格用于演示
    # 这里我们需要手动创建一个简单的Q4网格
    
    # 使用简化的层状结构
    Lx = param_dim.NE.thickness + param_dim.SP.thickness + param_dim.PE.thickness
    Ly = 0.1  # 高度 (m)
    nx = 20
    ny = 10
    
    println("  创建矩形网格:")
    println("  - 宽度: $(Lx*1e6) μm (负极+隔膜+正极)")
    println("  - 高度: $(Ly*1e3) mm")
    println("  - 单元数: $(nx) × $(ny)")
    
    # 创建简单的Q4网格
    mesh_th = JuBat.create_rectangular_Q4_mesh(Lx, Ly, nx, ny)
    
else
    # 果冻卷电池，使用专用网格生成函数
    println("  创建果冻卷2D网格...")
    
    # 网格参数
    nx_theta = 60  # 周向单元数
    gsorder = 2    # 高斯积分阶数
    
    mesh_th = JuBat.jellyroll_Q4_mesh(param_dim; 
                                       nx=nx_theta, 
                                       gsorder=gsorder,
                                       crop_mode=:collector_seeded)
    
    println("✓ 果冻卷网格创建完成")
    println("  单元数: $(size(mesh_th.element, 1))")
    println("  节点数: $(mesh_th.nlen)")
end

# 添加到案例
case.mesh["thermal2D"] = mesh_th

# ========================================================================
# 4. 准备变量字典
# ========================================================================

println("\n" * "="^70)
println("步骤 3: 准备温度场和SOC分布")
println("="^70)

# 提取最后时刻的状态
variables = Dict{String, Union{Array{Float64},Float64}}()

# 温度场（假设均匀温度 + 小扰动）
T_final = result["temperature[K]"][end]
T_ref = param_dim.scale.T_ref
T_nd = T_final / T_ref
T0 = 298.0 / T_ref

# 创建温度场（单元中心温度略高）
ne = size(mesh_th.element, 1)
T_nodes = fill(T_nd, mesh_th.nlen)

# 给内部节点添加温度梯度（模拟加热）
if hasproperty(param_dim.cell, :Rin)
    # 果冻卷：径向温度梯度
    x = mesh_th.node[:, 1]
    y = mesh_th.node[:, 2]
    r = hypot.(x, y)
    r_min, r_max = minimum(r), maximum(r)
    T_gradient = 5.0 / T_ref  # 5K 温差
    T_nodes = T_nd .+ T_gradient .* (r .- r_min) ./ (r_max - r_min)
else
    # 矩形：轴向温度梯度
    x = mesh_th.node[:, 1]
    x_min, x_max = minimum(x), maximum(x)
    T_gradient = 5.0 / T_ref
    T_nodes = T_nd .+ T_gradient .* (x .- x_min) ./ (x_max - x_min)
end

variables["T_nodes"] = T_nodes

println("✓ 温度场创建完成")
println("  平均温度: $(mean(T_nodes) * T_ref) K")
println("  温度范围: $(minimum(T_nodes) * T_ref) - $(maximum(T_nodes) * T_ref) K")
println("  温差: $((maximum(T_nodes) - minimum(T_nodes)) * T_ref) K")

# SOC分布（从电化学模型获取）
if haskey(result, "negative electrode lithium concentration")
    cs_n = result["negative particle lithium concentration"][:, end]
    cs_p = result["positive particle lithium concentration"][:, end]
    variables["negative electrode lithium concentration"] = cs_n
    variables["positive electrode lithium concentration"] = cs_p
    
    SOC_global = mean(cs_n) / param_dim.NE.cs_max
    println("✓ SOC分布创建完成")
    println("  全局SOC: $(SOC_global)")
else
    SOC_global = 0.7  # 假设70% SOC
    println("  使用假设的SOC: $(SOC_global)")
end

# ========================================================================
# 5. 计算宏观热-扩散应力
# ========================================================================

println("\n" * "="^70)
println("步骤 4: 计算宏观热-扩散应力")
println("="^70)

# 调用宏观应力计算函数
try
    variables = JuBat.diffusion_stress_2D(case, variables)
    
    println("✓ 应力场计算完成")
    
    # 提取结果
    σ_xx = variables["diffusion stress xx"]
    σ_yy = variables["diffusion stress yy"]
    σ_xy = variables["diffusion stress xy"]
    σ_vm = variables["diffusion stress vonMises"]
    u_x = variables["displacement x"]
    u_y = variables["displacement y"]
    
    # 统计信息
    println("\n应力统计 [MPa]:")
    println("  σ_xx: min=$(minimum(σ_xx)/1e6), max=$(maximum(σ_xx)/1e6), mean=$(mean(σ_xx)/1e6)")
    println("  σ_yy: min=$(minimum(σ_yy)/1e6), max=$(maximum(σ_yy)/1e6), mean=$(mean(σ_yy)/1e6)")
    println("  σ_xy: min=$(minimum(σ_xy)/1e6), max=$(maximum(σ_xy)/1e6), mean=$(mean(σ_xy)/1e6)")
    println("  σ_vm: min=$(minimum(σ_vm)/1e6), max=$(maximum(σ_vm)/1e6), mean=$(mean(σ_vm)/1e6)")
    
    println("\n位移统计 [μm]:")
    println("  u_x: min=$(minimum(u_x)*1e6), max=$(maximum(u_x)*1e6), mean=$(mean(u_x)*1e6)")
    println("  u_y: min=$(minimum(u_y)*1e6), max=$(maximum(u_y)*1e6), mean=$(mean(u_y)*1e6)")
    
    # ====================================================================
    # 6. 可视化结果
    # ====================================================================
    
    println("\n" * "="^70)
    println("步骤 5: 可视化应力场")
    println("="^70)
    
    # 计算单元中心坐标
    x_elem = zeros(Float64, ne)
    y_elem = zeros(Float64, ne)
    for e in 1:ne
        nodes = mesh_th.element[e, :]
        x_elem[e] = mean(mesh_th.node[nodes, 1])
        y_elem[e] = mean(mesh_th.node[nodes, 2])
    end
    
    # 创建图形
    p1 = scatter(x_elem, y_elem, marker_z=σ_xx./1e6, 
                 color=:viridis, markersize=3,
                 xlabel="x [m]", ylabel="y [m]",
                 title="σxx [MPa]", colorbar=true,
                 aspect_ratio=:equal)
    
    p2 = scatter(x_elem, y_elem, marker_z=σ_yy./1e6,
                 color=:viridis, markersize=3,
                 xlabel="x [m]", ylabel="y [m]",
                 title="σyy [MPa]", colorbar=true,
                 aspect_ratio=:equal)
    
    p3 = scatter(x_elem, y_elem, marker_z=σ_xy./1e6,
                 color=:viridis, markersize=3,
                 xlabel="x [m]", ylabel="y [m]",
                 title="σxy [MPa]", colorbar=true,
                 aspect_ratio=:equal)
    
    p4 = scatter(x_elem, y_elem, marker_z=σ_vm./1e6,
                 color=:plasma, markersize=3,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Von Mises Stress [MPa]", colorbar=true,
                 aspect_ratio=:equal)
    
    plot_stress = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 1000))
    savefig(plot_stress, "output/thermal_diffusion_stress_field.png")
    println("✓ 应力场图保存至: output/thermal_diffusion_stress_field.png")
    
    # 位移场
    p5 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2], 
                 marker_z=u_x.*1e6, 
                 color=:viridis, markersize=2,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement u_x [μm]", colorbar=true,
                 aspect_ratio=:equal)
    
    p6 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2],
                 marker_z=u_y.*1e6,
                 color=:viridis, markersize=2,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement u_y [μm]", colorbar=true,
                 aspect_ratio=:equal)
    
    u_mag = hypot.(u_x, u_y)
    p7 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2],
                 marker_z=u_mag.*1e6,
                 color=:plasma, markersize=2,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement Magnitude [μm]", colorbar=true,
                 aspect_ratio=:equal)
    
    plot_disp = plot(p5, p6, p7, layout=(1,3), size=(1800, 500))
    savefig(plot_disp, "output/thermal_diffusion_displacement_field.png")
    println("✓ 位移场图保存至: output/thermal_diffusion_displacement_field.png")
    
    # 温度场
    T_elem_plot = zeros(Float64, ne)
    for e in 1:ne
        nodes = mesh_th.element[e, :]
        T_elem_plot[e] = mean(T_nodes[nodes]) * T_ref
    end
    
    p8 = scatter(x_elem, y_elem, marker_z=T_elem_plot,
                 color=:hot, markersize=3,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Temperature [K]", colorbar=true,
                 aspect_ratio=:equal)
    
    savefig(p8, "output/thermal_field.png")
    println("✓ 温度场图保存至: output/thermal_field.png")
    
catch e
    println("❌ 应力计算失败:")
    println(e)
    rethrow(e)
end

# ========================================================================
# 7. 总结
# ========================================================================

println("\n" * "="^70)
println("计算完成总结")
println("="^70)
println("""
本示例展示了如何计算宏观层面的热-扩散应力：

1. 理论基础：
   - 热应力：ε_th = α * ΔT
   - 扩散应力：ε_diff = β_c * c_s_max * ΔSOC
   - 总初始应变：ε_0 = ε_th + ε_diff

2. 有限元方法：
   - 使用Q4单元进行2D空间离散
   - 平面应力假设
   - 应力平衡方程：∇·σ = 0
   - 本构关系：σ = E/(1-ν²) * [(ε - ε_0) + ν(ε - ε_0)]

3. 输出结果：
   - 应力场：σ_xx, σ_yy, σ_xy
   - Von Mises等效应力：σ_vm
   - 位移场：u_x, u_y

4. 应用场景：
   - 电池结构完整性分析
   - 失效模式预测（裂纹、分层）
   - 优化设计（材料选择、几何参数）

参考文档：
- docs/Diffusion_Stress_Macroscale_Theory.md
""")

println("\n完成! 🎉")
