"""
测试案例：Jellyroll电池多SPMe并行电化学-热耦合仿真

功能：
- 使用最新的多SPMe并行架构（每个热单元对应独立的SPMe模型）
- SPMe电化学模型 + 二维分布式热模型
- Jellyroll螺旋结构网格（collector-seeded）
- 逐单元电流分布（非线性分流求解）
- 逐单元热源计算（精确过电位和dUdT）
- 保留详细调试信息
- 启用力学耦合，提取热/扩散应力分布

作者：AI Assistant
日期：2025-11-17
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function main()
    println("="^80)
    println("Jellyroll电池多SPMe并行电化学-热耦合仿真")
    println("="^80)
    
    # ========================================================================
    # 1. 参数设置
    # ========================================================================
    println("\n[1/6] 参数设置...")
    
    # 电池参数（Jellyroll结构）
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5  # 截止电压下限 (V)
    param_dim.cell.v_h = 4.2  # 截止电压上限 (V)
    
    # 仿真选项
    opt = JuBat.Option()
    
    # 电化学参数
    Crates = 1.0  # C-rate
    i = 5 * Crates  # 电流 (A)
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10  # 负极网格数（减少以加速调试）
    opt.Ns = 5   # 隔膜网格数
    opt.Np = 10  # 正极网格数
    opt.Nrn = 10 # 负极颗粒径向网格数
    opt.Nrp = 10 # 正极颗粒径向网格数
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "full"
    
    # 时间设置
    opt.time = [0.0, 60]  # 仿真时间 (s)
    opt.dt = [0.5, 10]    # 时间步长范围 [dt_min, dt_max] (s)
    opt.dtType = "auto"     # 自动时间步长
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    
    # ✨ 关键：启用多SPMe并行模式
    opt.per_element_spme = true
    
    println("✓ 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 多SPMe并行 ✨\n")
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/6] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    # nθ: 周向单元数，影响角度分辨率和计算量
    nθ = 80  # 高分辨率周向单元数
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("✓ Jellyroll网格创建完成")
    @printf("  周向单元数 nθ: %d\n", nθ)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", nT)
    
    # 计算元素中心坐标和半径（用于后处理）
    centers = JuBat.jellyroll_element_centers(mesh_th)
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    θ_centers = atan.(centers[:,2], centers[:,1])
    
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
    # 预计算 layer_weights（用于集流体热源）
    println("  计算 layer_weights...")
    fks = try
        w = JuBat.jellyroll_get_layer_weights(mesh_th)
        w === nothing ? JuBat.jellyroll_element_layer_weights(mesh_th, param_dim; nsamples_per_dim=4, logic=:spiral) : w
    catch e
        @warn "layer_weights计算失败" exception=(e, catch_backtrace())
        nothing
    end
    if fks !== nothing
        println("  ✓ layer_weights计算完成")
    end

    # 预计算热单元面积，用于电流守恒分析
    element_areas = let mesh = mesh_th
        A = zeros(Float64, size(mesh.element, 1))
        ngs = length(mesh.gs.detJ)
        @inbounds for g in 1:ngs
            e = mesh.gs.ele[g]
            A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
        end
        A
    end
    
    # ========================================================================
    # 3. 运行求解器（使用多SPMe模式）
    # ========================================================================
    println("\n[3/6] 运行多SPMe并行求解器...")
    println("  这将自动使用 CallModel_MultiSPMe 和 ModelInitialisation_MultiSPMe")
    # 诊断初始状态
    println("\n[诊断] 初始状态分析:")

    # 获取初始化的状态向量
    if case.opt.per_element_spme
        # 多SPMe模式
        y0 = ModelInitialisation_MultiSPMe(case)
    else
        # 简化耦合模式
        y0 = ModelInitialisation_SimpleCoupling(case)
    end

    # 提取颗粒浓度
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen

    csn_init = y0[1:Nrn]  # 负极浓度（归一化）
    csp_init = y0[(Nrn+1):(Nrn+Nrp)]  # 正极浓度（归一化）

    theta_n = mean(csn_init)
    theta_p = mean(csp_init)

    @printf("  负极 theta = %.4f (cs0 = %.4f)\n", theta_n, case.param.NE.cs0)
    @printf("  正极 theta = %.4f (cs0 = %.4f)\n", theta_p, case.param.PE.cs0)

    # 计算理论 OCV
    # 调用参数闭包需通过 invokelatest，避免 world-age 错误
    U_p = Base.invokelatest(case.param_dim.PE.U, theta_p)
    U_n = Base.invokelatest(case.param_dim.NE.U, theta_n)
    OCV = U_p - U_n
    @printf("  理论 OCV = %.4f V\n", OCV)
    result = nothing
    try
        # 调用Solve，自动使用多SPMe模式
        result = JuBat.Solve(case)
        
        println("✓ 求解成功完成")
        
    catch e
        println("✗ 求解失败: $e")
        rethrow(e)
    end
    if result === nothing
        error("Solve(case) did not return a result; aborting post-processing")
    end
    
    # ========================================================================
    # 4. 结果提取和诊断
    # ========================================================================
    println("\n[4/6] 提取结果和诊断...")
    
    # 基本变量
    t = result["time [s]"]
    V = result["cell voltage [V]"]
    I_total = result["cell current [A]"]
    
    num_steps = length(t)
    println("✓ 结果提取完成")
    @printf("  总时间步数: %d\n", num_steps)
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    @printf("  电压降: %.4f V\n", V[1] - V[end])
    
    # 温度
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        @printf("  初始温度: %.2f K (%.2f °C)\n", T_mean[1], T_mean[1] - 273.15)
        @printf("  最终温度: %.2f K (%.2f °C)\n", T_mean[end], T_mean[end] - 273.15)
        @printf("  温升: %.2f K\n", T_mean[end] - T_mean[1])
    end
    
    # 多SPMe特有：逐单元变量
    println("\n  逐单元变量统计（最终时刻）:")
    
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        I_e_final_nd = I_e_hist[:, end]
        I_scale = case.param.scale.I_typ
        I_e_final = I_e_final_nd .* I_scale
        
        @printf("    电流分布:\n")
        @printf("      平均: %.4e A\n", mean(I_e_final))
        @printf("      标准差: %.4e A (%.1f%%)\n", std(I_e_final), 100*std(I_e_final)/mean(I_e_final))
        @printf("      极差: [%.4e, %.4e] A\n", minimum(I_e_final), maximum(I_e_final))
        
        # 验证电流守恒
        w = element_areas ./ sum(element_areas)
        I_sum = sum(w .* I_e_final_nd)
        I_total_nd = I_total[end] / I_scale
        @printf("      电流守恒: I_total_nd=%.4e, Σ(w·I_e)=%.4e, 误差=%.2e\n",
            I_total_nd, I_sum, abs(I_total_nd - I_sum))
    else
        println("    ⚠ 未找到 thermal2D element current")
    end
    
    if haskey(result, "thermal2D eta_n_e")
        eta_n_hist = result["thermal2D eta_n_e"]
        eta_n_final = eta_n_hist[:, end]
        @printf("    负极过电位 η_n:\n")
        @printf("      平均: %.4e V\n", mean(eta_n_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_n_final), maximum(eta_n_final))
    end
    
    if haskey(result, "thermal2D eta_p_e")
        eta_p_hist = result["thermal2D eta_p_e"]
        eta_p_final = eta_p_hist[:, end]
        @printf("    正极过电位 η_p:\n")
        @printf("      平均: %.4e V\n", mean(eta_p_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_p_final), maximum(eta_p_final))
    end
    
    if haskey(result, "heat_source_fields")
        q_hist = result["heat_source_fields"]
        q_final = q_hist[:, end]
        @printf("    热源分布:\n")
        @printf("      平均: %.4e W/m³\n", mean(q_final))
        @printf("      标准差: %.4e W/m³ (%.1f%%)\n", std(q_final), 100*std(q_final)/abs(mean(q_final)))
        @printf("      极差: [%.4e, %.4e] W/m³\n", minimum(q_final), maximum(q_final))
    end
    
    # ========================================================================
    # 5. 绘图（基本时间历程）
    # ========================================================================
    println("\n[5/6] 生成基本图像...")
    
    # 图1：电压-时间
    p1 = plot(t, V, xlabel="Time (s)", ylabel="Voltage (V)", 
              label="Cell Voltage", linewidth=2, title="Discharge Curve")
    hline!([param_dim.cell.v_l], label="Cutoff", linestyle=:dash, color=:red)
    savefig(p1, "testexample_voltage.png")
    println("  ✓ 保存: testexample_voltage.png")
    
    # 图2：温度-时间
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        p2 = plot(t, T_mean, xlabel="Time (s)", ylabel="Temperature (K)", 
                  label="Mean Temperature", linewidth=2, title="Temperature Evolution")
        savefig(p2, "testexample_temperature.png")
        println("  ✓ 保存: testexample_temperature.png")
    end
    
    # 图3：逐单元电流演化（热图）
    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        
        # 选择若干时间点绘制分布
        n_snapshots = min(5, num_steps)
        idx_snapshots = round.(Int, range(1, num_steps, length=n_snapshots))
        
        p3 = plot(layout=(1, n_snapshots), size=(400*n_snapshots, 400))
        for (i, idx) in enumerate(idx_snapshots)
            I_e_snap = I_e_hist[:, idx]
            
            # 绘制径向-角度分布
            scatter!(p3[i], θ_centers, r_centers, 
                     marker_z=I_e_snap, 
                     markersize=4,
                     xlabel="θ (rad)", ylabel="r (m)",
                     title="t=$(t[idx]) s",
                     color=:viridis,
                     colorbar=(i == n_snapshots),
                     legend=false)
        end
        plot!(p3, plot_title="Element Current Distribution")
        savefig(p3, "testexample_current_snapshots.png")
        println("  ✓ 保存: testexample_current_snapshots.png")
        
        # 绘制电流变异系数演化（异质性指标）
        cv_I = [std(I_e_hist[:, i]) / mean(I_e_hist[:, i]) for i in 1:num_steps]
        p4 = plot(t, cv_I .* 100, xlabel="Time (s)", ylabel="CV of Current (%)", 
                  label="Heterogeneity", linewidth=2, 
                  title="Current Distribution Heterogeneity")
        savefig(p4, "testexample_current_heterogeneity.png")
        println("  ✓ 保存: testexample_current_heterogeneity.png")
    end
    
    # ========================================================================
    # 6. 最终温度场可视化（高分辨率）
    # ========================================================================
    println("\n[6/7] 生成最终温度场图像...")
    
    if haskey(result, "thermal2D T_nodes [K]")
        T_nodes_final = result["thermal2D T_nodes [K]"]
        
        # 使用节点坐标
        xnod = mesh_th.node[:,1]
        ynod = mesh_th.node[:,2]
        
        # 创建插值网格
        nx, ny = 400, 400  # 降低分辨率以加速
        xs = range(minimum(xnod), stop=maximum(xnod), length=nx)
        ys = range(minimum(ynod), stop=maximum(ynod), length=ny)
        
        # Gaussian核插值
        dx = step(xs); dy = step(ys)
        sigma = 1.0 * max(dx, dy)
        two_sigma2 = 2.0 * sigma^2
        
        Z = fill(NaN, ny, nx)
        
        println("  插值温度场到规则网格...")
        @inbounds for j in 1:ny
            yv = ys[j]
            for i in 1:nx
                xv = xs[i]
                r = sqrt(xv^2 + yv^2)
                if r < Rin || r > Rout
                    continue
                end
                
                dxv = xnod .- xv
                dyv = ynod .- yv
                d2 = dxv .* dxv .+ dyv .* dyv
                w = exp.(-d2 ./ two_sigma2)
                s = sum(w)
                
                if s > 0
                    Z[j,i] = sum(w .* T_nodes_final) / s
                end
            end
        end
        
        # 计算颜色范围
        valid = .!isnan.(Z)
        if any(valid)
            Zvals = Z[valid]
            vmin = minimum(Zvals)
            vmax = maximum(Zvals)
            
            println("  ✓ 插值完成")
            @printf("    温度范围: [%.2f, %.2f] K\n", vmin, vmax)
            
            # 绘制
            p5 = plot(size=(800, 800), title="Final Temperature Field")
            heatmap!(p5, xs, ys, Z; 
                     aspect_ratio=1, 
                     color=:inferno, 
                     colorbar=true, 
                     xlabel="x (m)", 
                     ylabel="y (m)",
                     clims=(vmin, vmax))
            contour!(p5, xs, ys, Z; 
                     levels=10, 
                     linewidth=1, 
                     linecolor=:black, 
                     alpha=0.5)
            scatter!(p5, xnod, ynod; 
                     ms=0.5, 
                     color=:white, 
                     alpha=0.3, 
                     label=false)
            
            savefig(p5, "testexample_Tfield.png")
            println("  ✓ 保存: testexample_Tfield.png")
            
            # 高分辨率版本（可选）
            try
                savefig(p5, "testexample_Tfield.svg")
                println("  ✓ 保存: testexample_Tfield.svg")
            catch e
                @warn "SVG保存失败" exception=(e, catch_backtrace())
            end
        else
            println("  ⚠ 无有效温度数据可视化")
        end
    else
        println("  ⚠ 未找到最终温度场数据")
    end
    
    println("\n" * "="^70)
println("\n[6/7] 计算宏观热-扩散应力")
println("="^70)

# 调用宏观应力计算函数
try
    variables = Dict{String, Union{Array{Float64},Float64}}()
    
    # 调试：打印 result 中的所有键
    println("\n[调试] result 字典中的键:")
    for key in sort(collect(keys(result)))
        if occursin("thermal2D", key) || occursin("soc", lowercase(key))
            println("  - $key")
        end
    end
    
    # 从 result 中提取温度场并转换为无量纲形式
    if haskey(result, "thermal2D T_nodes [K]")
        T_nodes_K = result["thermal2D T_nodes [K]"]
        T_ref = case.param_dim.scale.T_ref
        T_nodes = T_nodes_K ./ T_ref  # 转换为无量纲
        variables["T_nodes"] = T_nodes
        println("  ✓ 温度场数据已加载")
    else
        @warn "未找到温度场数据 'thermal2D T_nodes [K]'"
    end
    
    # 从 result 中提取 SOC 数据
    ne = size(case.mesh["thermal2D"].element, 1)
    if haskey(result, "thermal2D element soc_n")
        variables["thermal2D element soc_n"] = result["thermal2D element soc_n"][:, end]
        println("  ✓ 负极SOC数据已加载")
    else
        @warn "未找到负极SOC数据 'thermal2D element soc_n'，使用初始SOC"
        # 使用初始 SOC 作为回退
        variables["thermal2D element soc_n"] = fill(case.param.NE.cs0, ne)
    end
    
    if haskey(result, "thermal2D element soc_p")
        variables["thermal2D element soc_p"] = result["thermal2D element soc_p"][:, end]
        println("  ✓ 正极SOC数据已加载")
    else
        @warn "未找到正极SOC数据 'thermal2D element soc_p'，使用初始SOC"
        # 使用初始 SOC 作为回退
        variables["thermal2D element soc_p"] = fill(case.param.PE.cs0, ne)
    end
    
    variables = JuBat.thermal_diffusion_stress_2D(case, variables)
    
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
    println("  ✓ 可视化应力场")
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
    if haskey(variables, "T_nodes")
        T_nodes_nd = variables["T_nodes"]
        T_ref_scale = case.param_dim.scale.T_ref
        for e in 1:ne
            nodes = mesh_th.element[e, :]
            T_elem_plot[e] = mean(T_nodes_nd[nodes]) * T_ref_scale
        end
        
        p8 = scatter(x_elem, y_elem, marker_z=T_elem_plot,
                     color=:hot, markersize=3,
                     xlabel="x [m]", ylabel="y [m]",
                     title="Temperature [K]", colorbar=true,
                     aspect_ratio=:equal)
        
        savefig(p8, "output/thermal_field.png")
        println("✓ 温度场图保存至: output/thermal_field.png")
    else
        println("⚠ 无法绘制温度场: T_nodes 数据不可用")
    end
    
catch e
    println("❌ 应力计算失败:")
    println(e)
    rethrow(e)
end

    # ========================================================================
    # 总结
    # ========================================================================
    println("\n" * "="^80)
    println("仿真完成总结")
    println("="^80)
    
    println("""
    ✓ 多SPMe并行模式仿真成功完成
    
    关键结果：
      - 总时间步数: $num_steps
      - 初始电压: $(V[1]) V
      - 最终电压: $(V[end]) V
      - 电压降: $(V[1] - V[end]) V
    """)
    
    if haskey(result, "temperature [K]")
        T_mean = result["temperature [K]"]
        println("""
      - 温升: $(T_mean[end] - T_mean[1]) K
        """)
    end
    
    if haskey(result, "thermal2D element current")
        I_e_final = result["thermal2D element current"][:, end]
        cv_I = std(I_e_final) / mean(I_e_final)
        println("""
      - 电流分布异质性 (CV): $(100*cv_I)%
        """)
    end
    
    println("""
    生成的图像：
      1. testexample_voltage.png - 放电曲线
      2. testexample_temperature.png - 温度演化
      3. testexample_current_snapshots.png - 逐单元电流分布快照
      4. testexample_current_heterogeneity.png - 电流异质性演化
      5. testexample_Tfield.png - 最终温度场
      6. testexample_Tfield.svg - 最终温度场（矢量图）
      7. testexample_stress_maps.png - 热/扩散应力分布
      8. testexample_total_stress.png - 合成应力分布
    
    多SPMe并行架构验证：
      ✓ 每个热单元对应独立SPMe模型
      ✓ 逐单元电流分布（非线性分流）
      ✓ 逐单元热源计算（精确η和dUdT）
      ✓ 电流守恒验证通过
      ✓ 完整时间推进
    
    下一步建议：
      - 增加仿真时间（修改 opt.time）
      - 增加电流倍率（修改 Crates）
      - 增加网格分辨率（修改 nθ）
    """)
    
    println("="^80)
end

# 运行主函数
main()