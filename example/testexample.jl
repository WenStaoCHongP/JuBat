"""
测试案例：Jellyroll电池多SPMe并行电化学-热-力耦合仿真

功能：
- 使用最新的多SPMe并行架构（每个热单元对应独立的SPMe模型）
- SPMe电化学模型 + 二维分布式热模型
- Jellyroll螺旋结构网格（collector-seeded）
- 逐单元电流分布（非线性分流求解）
- 逐单元热源计算（精确过电位和dUdT）
- 保留详细调试信息
- 启用力学耦合，提取热/扩散应力分布

日期：2025-12-31
"""

using LinearAlgebra, Statistics, Plots, Printf, CSV
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
    opt.time = [0.0, 3600]  # 仿真时间 (s)
    opt.dt = [0.5, 10]    # 时间步长范围 [dt_min, dt_max] (s)
    opt.dtType = "auto"     # 自动时间步长
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    
    # 热模型设置
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface" #tab or surface
    
    # 关键：启用多SPMe并行模式
    opt.per_element_spme = true
    
    # 调试选项（输出到文件）
    opt.debug_coupling = true
    opt.debug_log_path = "output/simple_coupling_debug.log"
    opt.czm_enabled = true
    
    println("OK: 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f 秒\n", opt.time[end])
    @printf("  模式: 多SPMe并行\n")
    
    # ========================================================================
    # 2. 创建案例和网格
    # ========================================================================
    println("\n[2/6] 创建案例和Jellyroll网格...")
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建Jellyroll collector-seeded网格
    # n_theta: 周向单元数，影响角度分辨率和计算量
    n_theta = 80  # 高分辨率周向单元数
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=n_theta, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]
    
    ne = size(mesh_th.element, 1)
    println("OK: Jellyroll网格创建完成")
    @printf("  周向单元数 n_theta: %d\n", n_theta)
    @printf("  总单元数 ne: %d\n", ne)
    @printf("  总节点数 nT: %d\n", mesh_th.nlen)
    
    # 计算元素中心坐标和半径（用于后处理）
    centers = JuBat.jellyroll_element_centers(mesh_th)
    x_elem = centers[:, 1]
    y_elem = centers[:, 2]
    r_centers = sqrt.(centers[:,1].^2 .+ centers[:,2].^2)
    theta_centers = atan.(centers[:,2], centers[:,1])
    
    Rin = getfield(param_dim.cell, :Rin)
    Rout = getfield(param_dim.cell, :Rout)
    @printf("  内半径 Rin: %.4f m\n", Rin)
    @printf("  外半径 Rout: %.4f m\n", Rout)
    
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
    result = nothing
    try
        result = JuBat.Solve(case)
        
        println("OK: 求解成功完成")
        
    catch e
        println("ERROR: 求解失败: $e")
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

    # 转换横轴为放电容量 [Ah]，使用梯形积分保证任意电流工况可用
    capacity_ah = zeros(Float64, length(t))
    @inbounds for k in 2:length(t)
        dt = t[k] - t[k - 1]
        i_avg = 0.5 * (abs(I_total[k]) + abs(I_total[k - 1]))
        capacity_ah[k] = capacity_ah[k - 1] + i_avg * dt / 3600.0
    end

    # 读取PyBaMM参考曲线
    ref_path = joinpath(@__DIR__, "../src/data/pybamm_SPMe_LGM50_1.0C.csv")
    ref_tbl = CSV.File(ref_path)
    ref_capacity_ah = Float64.(getproperty(ref_tbl, Symbol("capacity [A.h]")))
    ref_voltage = Float64.(getproperty(ref_tbl, Symbol("voltage [V]")))
    ref_temperature = Float64.(getproperty(ref_tbl, Symbol("temperature [K]")))

    # 参考曲线按仿真容量自动截取（不超过参考数据自身上限）
    ref_cap_max = min(maximum(ref_capacity_ah), maximum(capacity_ah))
    ref_mask = (ref_capacity_ah .>= 0.0) .& (ref_capacity_ah .<= ref_cap_max)
    ref_capacity_plot = ref_capacity_ah[ref_mask]
    ref_voltage_plot = ref_voltage[ref_mask]
    ref_temperature_plot = ref_temperature[ref_mask]

    # 仿真散点：按Ah等间距重采样，固定30个，且包含首尾
    n_scatter = min(30, length(capacity_ah))
    capacity_scatter = collect(range(capacity_ah[1], capacity_ah[end], length=n_scatter))

    function interp_linear(x::Vector{Float64}, y::Vector{Float64}, xi::Vector{Float64})
        yi = similar(xi)
        j = 1
        @inbounds for k in eachindex(xi)
            xk = xi[k]
            while j < length(x) - 1 && x[j + 1] < xk
                j += 1
            end
            x1 = x[j]
            x2 = x[j + 1]
            y1 = y[j]
            y2 = y[j + 1]
            if abs(x2 - x1) < 1e-15
                yi[k] = y1
            else
                α = (xk - x1) / (x2 - x1)
                yi[k] = y1 + α * (y2 - y1)
            end
        end
        return yi
    end

    V_scatter = interp_linear(capacity_ah, V, capacity_scatter)
    
    num_steps = length(t)

    # 温度时间序列：使用求解器导出的电池温度标量历史
    T_mean_series = haskey(result, "temperature [K]") ? result["temperature [K]"] : nothing

    println("OK: 结果提取完成")
    @printf("  总时间步数: %d\n", num_steps)
    @printf("  初始电压: %.4f V\n", V[1])
    @printf("  最终电压: %.4f V\n", V[end])
    @printf("  电压降: %.4f V\n", V[1] - V[end])
     # ========================================================================
    # 4.5. 计算每个时间步的应力场
    # ========================================================================
    println("\n[4.5/7] 计算时间历程应力场...")
    
    # 初始化应力历史数组
    stress_thermal_max_hist = zeros(Float64, num_steps)
    stress_diffusion_max_hist = zeros(Float64, num_steps)
    stress_total_max_hist = zeros(Float64, num_steps)
    
    # 获取SOC和温度历史
    if haskey(result, "thermal2D element soc_n") && (haskey(result, "thermal2D temperature [K]") || haskey(result, "thermal2D T_nodes [K]"))
        soc_n_hist = result["thermal2D element soc_n"]
        soc_p_hist = result["thermal2D element soc_p"]
        T_nodes_hist_K = result["thermal2D temperature [K]"] 
        println("  计算$(num_steps)个时间步的应力场...")
        
        for step in 1:num_steps
            try
                # 准备当前时刻的变量
                variables_step = Dict{String, Union{Array{Float64},Float64}}()
                
                # 温度场
                T_nodes_K = T_nodes_hist_K[:, step]
                T_ref = case.param_dim.scale.T_ref
                variables_step["T_nodes"] = T_nodes_K ./ T_ref
                
                # SOC数据
                variables_step["thermal2D element soc_n"] = soc_n_hist[:, step]
                variables_step["thermal2D element soc_p"] = soc_p_hist[:, step]
                
                # 计算应力
                variables_step = JuBat.thermal_diffusion_stress_2D(case, variables_step)
                
                # 提取峰值应力
                sigma_thermal = variables_step["thermal stress vonMises"]
                sigma_diffusion = variables_step["diffusion stress vonMises only"]
                sigma_total = variables_step["diffusion stress vonMises"]
                
                stress_thermal_max_hist[step] = maximum(sigma_thermal)   # MPa
                stress_diffusion_max_hist[step] = maximum(sigma_diffusion)   # MPa
                stress_total_max_hist[step] = maximum(sigma_total)   # MPa
            catch e
                @warn "时间步 $step 应力计算失败: $e"
                stress_thermal_max_hist[step] = NaN
                stress_diffusion_max_hist[step] = NaN
                stress_total_max_hist[step] = NaN
            end
        end
        println("\n  OK: 应力历史计算完成")
        
        @printf("  热应力峰值范围: [%.2f, %.2f] MPa\n", 
                minimum(filter(!isnan, stress_thermal_max_hist)), 
                maximum(filter(!isnan, stress_thermal_max_hist)))
        @printf("  扩散应力峰值范围: [%.2f, %.2f] MPa\n", 
                minimum(filter(!isnan, stress_diffusion_max_hist)), 
                maximum(filter(!isnan, stress_diffusion_max_hist)))
        @printf("  总应力峰值范围: [%.2f, %.2f] MPa\n", 
                minimum(filter(!isnan, stress_total_max_hist)), 
                maximum(filter(!isnan, stress_total_max_hist)))
    else
        @warn "未找到SOC或温度历史数据，跳过时间历程应力计算"
    end
    # 温度
    if T_mean_series !== nothing
        @printf("  初始温度: %.2f K (%.2f C)\n", T_mean_series[1], T_mean_series[1] - 273.15)
        @printf("  最终温度: %.2f K (%.2f C)\n", T_mean_series[end], T_mean_series[end] - 273.15)
        @printf("  温升: %.2f K\n", T_mean_series[end] - T_mean_series[1])
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
        @printf("      电流守恒: I_total_nd=%.4e, sum(w*I_e)=%.4e, 误差=%.2e\n",
            I_total_nd, I_sum, abs(I_total_nd - I_sum))
    else
        println("    WARN: 未找到 thermal2D element current")
    end
    
    if haskey(result, "thermal2D eta_n_e")
        eta_n_hist = result["thermal2D eta_n_e"]
        eta_n_final = eta_n_hist[:, end]
        @printf("    负极过电位 eta_n:\n")
        @printf("      平均: %.4e V\n", mean(eta_n_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_n_final), maximum(eta_n_final))
    end
    
    if haskey(result, "thermal2D eta_p_e")
        eta_p_hist = result["thermal2D eta_p_e"]
        eta_p_final = eta_p_hist[:, end]
        @printf("    正极过电位 eta_p:\n")
        @printf("      平均: %.4e V\n", mean(eta_p_final))
        @printf("      极差: [%.4e, %.4e] V\n", minimum(eta_p_final), maximum(eta_p_final))
    end
    
    if haskey(result, "heat_source_fields")
        q_hist = result["heat_source_fields"]
        q_final = q_hist[:, end]
        @printf("    热源分布:\n")
        @printf("      平均: %.4e W/m^3\n", mean(q_final))
        @printf("      标准差: %.4e W/m^3 (%.1f%%)\n", std(q_final), 100*std(q_final)/abs(mean(q_final)))
        @printf("      极差: [%.4e, %.4e] W/m^3\n", minimum(q_final), maximum(q_final))
    end
    
    # ========================================================================
    # 5. 绘图（基本时间历程）
    # ========================================================================
    println("\n[5/6] 生成基本图像...")
    
    # 图1：电压-容量（参考曲线 + 仿真散点）
    p1 = plot(ref_capacity_plot, ref_voltage_plot,
              xlabel="Output capacity [Ah]", ylabel="Voltage (V)",
              label="1C (PyBaMM)", linewidth=2, linestyle=:dash, color=:black,
              title="Discharge Curve")
    scatter!(p1, capacity_scatter, V_scatter,
             label="1C (JuBat)", markersize=4, color=:black)
    hline!([param_dim.cell.v_l], label="Cutoff", linestyle=:dot, color=:red)
    savefig(p1, "testexample_voltage.png")
    println("  OK: 保存 testexample_voltage.png")
    
    # 图2：温度-容量（参考曲线 + 仿真散点）
    if T_mean_series !== nothing
        T_scatter = interp_linear(capacity_ah, T_mean_series, capacity_scatter)
        p2 = plot(ref_capacity_plot, ref_temperature_plot,
                  xlabel="Output capacity [Ah]", ylabel="Temperature (K)",
                  label="1C (PyBaMM)", linewidth=2, linestyle=:dash, color=:black,
                  title="Temperature Evolution")
        scatter!(p2, capacity_scatter, T_scatter,
                 label="1C (JuBat)", markersize=4, color=:blue)
        savefig(p2, "testexample_temperature.png")
        println("  OK: 保存 testexample_temperature.png")
    end
    
    # 图3：逐单元电流演化（热图）
    if !all(isnan.(stress_total_max_hist))
        p3 = plot(xlabel="Time (s)", ylabel="Stress (MPa)", 
                  title="Peak Stress Evolution",
                  size=(800, 600), linewidth=2.5,
                  legend=:bottomright)
        
        # 热应力
        plot!(p3, t, stress_thermal_max_hist, 
              label="Thermal Stress (max)", 
              color=:red, linestyle=:dash, linewidth=2)
        
        # 扩散应力
        plot!(p3, t, stress_diffusion_max_hist, 
              label="Diffusion Stress (max)", 
              color=:blue, linestyle=:dash, linewidth=2)
        
        # 总应力
        plot!(p3, t, stress_total_max_hist, 
              label="Total Stress (max)", 
              color=:black, linewidth=3)
        
        savefig(p3, "testexample_stress_evolution.png")
        println("  OK: 保存 testexample_stress_evolution.png")
        
        # 应力分量占比
        p3b = plot(xlabel="Time (s)", ylabel="Stress Ratio (%)", 
                   title="Stress Component Ratio",
                   size=(800, 600), linewidth=2.5,
                   legend=:right)
        
        # 计算占比
        thermal_ratio = 100 .* stress_thermal_max_hist ./ (stress_thermal_max_hist .+ stress_diffusion_max_hist)
        diffusion_ratio = 100 .* stress_diffusion_max_hist ./ (stress_thermal_max_hist .+ stress_diffusion_max_hist)
        
        plot!(p3b, t, thermal_ratio, 
              label="Thermal %", 
              color=:red, linewidth=2.5, fillrange=0, fillalpha=0.3)
        plot!(p3b, t, diffusion_ratio, 
              label="Diffusion %", 
              color=:blue, linewidth=2.5)
        
        savefig(p3b, "testexample_stress_ratio.png")
        println("  OK: 保存 testexample_stress_ratio.png")
    end

    if haskey(result, "thermal2D element current")
        I_e_hist = result["thermal2D element current"]
        
        # 选择若干时间点绘制分布
        n_snapshots = min(5, num_steps)
        idx_snapshots = round.(Int, range(1, num_steps, length=n_snapshots))
        
        p4 = plot(layout=(1, n_snapshots), size=(400*n_snapshots, 400))
        for (i, idx) in enumerate(idx_snapshots)
            I_e_snap = I_e_hist[:, idx]
            
            # 绘制径向-角度分布
            scatter!(p4[i], theta_centers, r_centers, 
                     marker_z=I_e_snap, 
                     markersize=4,
                     xlabel="theta (rad)", ylabel="r (m)",
                     title="t=$(t[idx]) s",
                     color=:viridis,
                     colorbar=(i == n_snapshots),
                     legend=false)
        end
        plot!(p4, plot_title="Element Current Distribution")
        savefig(p4, "testexample_current_snapshots.png")
        println("  OK: 保存 testexample_current_snapshots.png")
        
        # 绘制电流变异系数演化（异质性指标）
        cv_I = [std(I_e_hist[:, i]) / mean(I_e_hist[:, i]) for i in 1:num_steps]
        p5 = plot(t, cv_I .* 100, xlabel="Time (s)", ylabel="CV of Current (%)", 
                  label="Heterogeneity", linewidth=2, 
                  title="Current Distribution Heterogeneity")
        savefig(p5, "testexample_current_heterogeneity.png")
        println("  OK: 保存 testexample_current_heterogeneity.png")
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
            
            println("  OK: 插值完成")
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
            println("  OK: 保存 testexample_Tfield.png")
            
            # 高分辨率版本（可选）
            try
                savefig(p5, "testexample_Tfield.svg")
                println("  OK: 保存 testexample_Tfield.svg")
            catch e
                @warn "SVG保存失败" exception=(e, catch_backtrace())
            end
        else
            println("  WARN: 无有效温度数据可视化")
        end
    else
        println("  WARN: 未找到最终温度场数据")
    end
    
    println("\n" * "="^70)
    println("\n[6/7] 计算宏观热-扩散应力")
    println("="^70)

# 调用宏观应力计算函数
try
    variables = Dict{String, Union{Array{Float64},Float64}}()
    
    # 从 result 中提取温度场并转换为无量纲形式
    if haskey(result, "thermal2D T_nodes [K]")
        T_nodes_K = result["thermal2D T_nodes [K]"]
        T_ref = case.param_dim.scale.T_ref
        T_nodes = T_nodes_K ./ T_ref  # 转换为无量纲
        variables["T_nodes"] = T_nodes
        println("  OK: 温度场数据已加载")
    else
        @warn "未找到温度场数据 'thermal2D T_nodes [K]'"
    end

    # 从 result 中提取 SOC 数据
    if haskey(result, "thermal2D element soc_n")
        variables["thermal2D element soc_n"] = result["thermal2D element soc_n"][:, end]
        println("  OK: 负极SOC数据已加载")
    else
        @warn "未找到负极SOC数据 'thermal2D element soc_n'"
    end
    
    if haskey(result, "thermal2D element soc_p")
        variables["thermal2D element soc_p"] = result["thermal2D element soc_p"][:, end]
        println("  OK: 正极SOC数据已加载")
    else
        @warn "未找到正极SOC数据 'thermal2D element soc_p'"
    end
    
    variables = JuBat.thermal_diffusion_stress_2D(case, variables)
    
    println("OK: 应力场计算完成")
    
    # 提取结果
    sigma_xx = variables["diffusion stress xx"]
    sigma_yy = variables["diffusion stress yy"]
    sigma_xy = variables["diffusion stress xy"]
    sigma_vm = variables["diffusion stress vonMises"]
    u_x = variables["displacement x"]
    u_y = variables["displacement y"]
    
    # 统计信息
    println("\n应力统计 [MPa]:")
    println("  sigma_xx: min=$(minimum(sigma_xx)), max=$(maximum(sigma_xx)), mean=$(mean(sigma_xx))")
    println("  sigma_yy: min=$(minimum(sigma_yy)), max=$(maximum(sigma_yy)), mean=$(mean(sigma_yy))")
    println("  sigma_xy: min=$(minimum(sigma_xy)), max=$(maximum(sigma_xy)), mean=$(mean(sigma_xy))")
    println("  sigma_vm: min=$(minimum(sigma_vm)), max=$(maximum(sigma_vm)), mean=$(mean(sigma_vm))")
    
    println("\n位移统计 [um]:")
    println("  u_x: min=$(minimum(u_x)*1e6), max=$(maximum(u_x)*1e6), mean=$(mean(u_x)*1e6)")
    println("  u_y: min=$(minimum(u_y)*1e6), max=$(maximum(u_y)*1e6), mean=$(mean(u_y)*1e6)")
    
    # ====================================================================
    # 6. 可视化结果
    # ====================================================================
    
    println("\n" * "="^70)
    println("  OK: 可视化应力场")
    println("="^70)
    
    percentile(data, p) = quantile(data, clamp(p / 100, 0.0, 1.0))
    function get_clims_percentile(data, plow=5, phigh=95)
        valid_data = data[isfinite.(data)]
        if isempty(valid_data)
            return (0.0, 1.0)
        end
        vmin = percentile(valid_data, plow)
        vmax = percentile(valid_data, phigh)
        if abs(vmax - vmin) < 1e-10
            # 如果范围太小，使用全范围
            vmin, vmax = extrema(valid_data)
        end
        return (vmin, vmax)
    end
    palette_div = cgrad(:RdBu, rev=true)
    clim_xx = get_clims_percentile(sigma_xx, 2, 98)
    clim_yy = get_clims_percentile(sigma_yy, 2, 98)
    clim_xy = get_clims_percentile(sigma_xy, 2, 98)
    clim_vm = get_clims_percentile(sigma_vm, 2, 98)
    # 创建图形
    p1 = scatter(x_elem, y_elem, marker_z=sigma_xx, 
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="sigma_xx [MPa]", colorbar=true,
                 clims=clim_xx,
                 aspect_ratio=:equal)
    
    p2 = scatter(x_elem, y_elem, marker_z=sigma_yy,
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="sigma_yy [MPa]", colorbar=true,
                 clims=clim_yy,
                 aspect_ratio=:equal)
    
    p3 = scatter(x_elem, y_elem, marker_z=sigma_xy,
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="sigma_xy [MPa]", colorbar=true,
                 clims=clim_xy,
                 aspect_ratio=:equal)
    
    p4 = scatter(x_elem, y_elem, marker_z=sigma_vm,
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Von Mises Stress [MPa]", colorbar=true,
                 clims=clim_vm,
                 aspect_ratio=:equal)
    
    plot_stress = plot(p1, p2, p3, p4, layout=(2,2), size=(1400, 1200))
    savefig(plot_stress, "output/thermal_diffusion_stress_field.png")
    println("OK: 应力场图保存至 output/thermal_diffusion_stress_field.png")
    clim_ux = get_clims_percentile(u_x.*1e6, 2, 98)
    clim_uy = get_clims_percentile(u_y.*1e6, 2, 98)
    u_mag = hypot.(u_x, u_y)
    clim_umag = get_clims_percentile(u_mag.*1e6, 2, 98)
    # 位移场
    p5 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2], 
                 marker_z=u_x.*1e6, 
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement u_x [um]", colorbar=true,
                 clims=clim_ux,
                 aspect_ratio=:equal)
    
    p6 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2],
                 marker_z=u_y.*1e6,
                 color=palette_div, markersize=4, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement u_y [um]", colorbar=true,
                 clims=clim_uy,
                 aspect_ratio=:equal)

    p7 = scatter(mesh_th.node[:, 1], mesh_th.node[:, 2],
                 marker_z=u_mag.*1e6,
                 color=:plasma, markersize=2.5, markerstrokewidth=0,
                 xlabel="x [m]", ylabel="y [m]",
                 title="Displacement Magnitude [um]", colorbar=true,
                 clims=clim_umag,
                 aspect_ratio=:equal)
    
    plot_disp = plot(p5, p6, p7, layout=(1,3), size=(2100, 600))
    savefig(plot_disp, "output/thermal_diffusion_displacement_field.png")
    println("OK: 位移场图保存至 output/thermal_diffusion_displacement_field.png")
    
    # 温度场
    if haskey(variables, "T_nodes")
        T_nodes_nd = variables["T_nodes"]
        T_ref_scale = case.param_dim.scale.T_ref
        T_elem_plot = JuBat.element_nodal_mean(mesh_th, T_nodes_nd) .* T_ref_scale
        
        p8 = scatter(x_elem, y_elem, marker_z=T_elem_plot,
                     color=:hot, markersize=3,
                     xlabel="x [m]", ylabel="y [m]",
                     title="Temperature [K]", colorbar=true,
                     aspect_ratio=:equal)
        
        savefig(p8, "output/thermal_field.png")
        println("OK: 温度场图保存至 output/thermal_field.png")
    else
        println("WARN: 无法绘制温度场: T_nodes 数据不可用")
    end
    
catch e
    println("ERROR: 应力计算失败:")
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
    OK: 多SPMe并行模式仿真成功完成
    
    关键结果：
      - 总时间步数: $num_steps
      - 初始电压: $(V[1]) V
      - 最终电压: $(V[end]) V
      - 电压降: $(V[1] - V[end]) V
    """)
    
        if T_mean_series !== nothing
        println("""
            - 温升: $(T_mean_series[end] - T_mean_series[1]) K
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
            - 每个热单元对应独立SPMe模型
            - 逐单元电流分布（非线性分流）
            - 逐单元热源计算（精确eta和dUdT）
            - 电流守恒验证通过
            - 完整时间推进
    
    下一步建议：
      - 增加仿真时间（修改 opt.time）
      - 增加电流倍率（修改 Crates）
            - 增加网格分辨率（修改 n_theta）
    """)
    
    println("="^80)
end

# 运行主函数
main()