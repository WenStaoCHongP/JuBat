# 诊断脚本：分析温度反向增长的原因
# 目标：提取过电位、热源、散热功率等关键数据，验证理论分析

using Plots, CSV, DataFrames
include("../src/JuBat.jl")

println("=" ^70)
println("温度反向增长诊断分析")
println("=" ^70)

# 设置参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_l = 2.5
opt = JuBat.Option()
opt.thermalmodel = "lumped"
opt.dtType = "auto"
opt.jacobi = "update"
opt.model = "SPMe"

Crates = 1
I1C = 5
opt.Current = x-> I1C * Crates
opt.dt = [1 10] / Crates
opt.time = [0 3600/Crates]

# 创建case并初始化
case = JuBat.SetCase(param_dim, opt)
y0 = JuBat.ModelInitialisation(case)
yt = copy(y0)

# 数据记录数组
times = Float64[]
voltages = Float64[]
temperatures = Float64[]
eta_n_hist = Float64[]  # 负极过电位
eta_p_hist = Float64[]  # 正极过电位
eta_diff_hist = Float64[]  # 过电位差 (eta_p - eta_n)
Q_rxn_hist = Float64[]  # 反应热
Q_rev_hist = Float64[]  # 可逆热
Q_total_hist = Float64[]  # 总产热
Q_conv_hist = Float64[]  # 对流散热
dT_dt_hist = Float64[]  # 温度变化率
current_hist = Float64[]  # 实际电流

# 手动时间推进（模拟Solve函数的逻辑）
t = opt.time[1]
t_end = opt.time[2]
dt = opt.dt[1]
dt_max = opt.dt[2]
T = case.param.cell.T0  # 初始温度（无量纲）

println("\n开始模拟：")
println("  初始温度: $(T * case.param_dim.scale.T_ref) K")
println("  电流: $(I1C * Crates) A")
println("  截止电压: $(param_dim.cell.v_l) V")
println("  时间步: $(dt) - $(dt_max) s (nd)\n")

step = 0
while t < t_end
    step += 1
    
    # 求解电化学
    M, K, F, variables = JuBat.SPMe(case, yt, t; jacobi=opt.jacobi)
    
    # 更新温度到variables
    variables["temperature"] = [T]
    
    # 求解热方程
    MT, Q_net = JuBat.ThermalLumped(case, variables)
    
    # 提取关键变量
    V_cell = variables["cell voltage"] * case.param.scale.phi
    I_app = variables["cell current"]
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    
    # 计算热源分量
    Q_rxn = abs(I_app * (eta_p - eta_n))
    Q_rev = abs(I_app) * T * (case.param.PE.dUdT(csp_surf) - case.param.NE.dUdT(csn_surf))
    Q_total = Q_rxn + Q_rev  # 欧姆热在LGM50-SPMe中忽略
    
    # 计算散热
    h = case.param.cell.h
    A_cool = case.param.cell.cooling_surface
    T_amb = case.param.cell.T_amb
    T_K = T * case.param_dim.scale.T_ref
    Q_conv = h * A_cool * (T_K - T_amb)  # 单位：W
    
    # 转换为W
    Q_total_W = Q_total * case.param.scale.I_typ * case.param.scale.phi
    
    # 计算温度变化率
    dT_dt = Q_net / MT[1,1]  # 无量纲温度变化率
    
    # 记录数据
    push!(times, t * case.param_dim.scale.t0)
    push!(voltages, V_cell)
    push!(temperatures, T_K)
    push!(eta_n_hist, eta_n)
    push!(eta_p_hist, eta_p)
    push!(eta_diff_hist, eta_p - eta_n)
    push!(Q_rxn_hist, Q_rxn)
    push!(Q_rev_hist, Q_rev)
    push!(Q_total_hist, Q_total_W)
    push!(Q_conv_hist, Q_conv)
    push!(dT_dt_hist, dT_dt * case.param_dim.scale.T_ref / case.param_dim.scale.t0)  # K/s
    push!(current_hist, I_app * case.param_dim.cell.I1C)  # A
    
    # 诊断输出（仅前10步和关键时刻）
    if step <= 10 || step % 20 == 0 || V_cell < 2.6
        println("  步 $step | t=$(round(t*case.param_dim.scale.t0,digits=1))s | " *
                "V=$(round(V_cell,digits=3))V | T=$(round(T_K,digits=2))K | " *
                "η_diff=$(round(eta_p-eta_n,digits=4)) | " *
                "Q=$(round(Q_total_W,digits=2))W | Q_conv=$(round(Q_conv,digits=2))W | " *
                "dT/dt=$(round(dT_dt_hist[end],digits=4))K/s")
    end
    
    # 检查截止条件
    if V_cell < param_dim.cell.v_l
        println("\n达到截止电压 $(param_dim.cell.v_l) V，停止模拟")
        break
    end
    
    # 后向欧拉时间推进
    dt_actual = dt  # 简化：使用固定时间步
    A_thermal = (1.0 / dt_actual) * MT[1,1] + 0  # 无刚度项（集总模型）
    rhs_thermal = (1.0 / dt_actual) * MT[1,1] * T + Q_net
    T = rhs_thermal / A_thermal
    
    # 推进电化学状态
    A_ec = (1.0 / dt_actual) * M + K
    rhs_ec = (1.0 / dt_actual) * M * yt + F
    yt = A_ec \ rhs_ec
    
    # 推进时间
    t += dt_actual
end

println("\n" * "=" ^70)
println("分析完成，共 $step 步")
println("=" ^70)

# ====================
# 数据分析与可视化
# ====================

println("\n关键统计：")
println("  最高温度: $(round(maximum(temperatures),digits=2)) K")
println("  最低温度: $(round(minimum(temperatures),digits=2)) K")
println("  温度范围: $(round(maximum(temperatures)-minimum(temperatures),digits=2)) K")
println("  最大过电位差: $(round(maximum(eta_diff_hist),digits=4))")
println("  最小过电位差: $(round(minimum(eta_diff_hist),digits=4))")
println("  最大产热: $(round(maximum(Q_total_hist),digits=2)) W")
println("  最小产热: $(round(minimum(Q_total_hist),digits=2)) W")

# 找到温度下降的时间段
temp_decrease_indices = findall(x -> x < 0, diff(temperatures))
if !isempty(temp_decrease_indices)
    println("\n发现 $(length(temp_decrease_indices)) 个温度下降时间步：")
    for (i, idx) in enumerate(temp_decrease_indices[1:min(5, end)])
        println("  步 $idx: T=$(round(temperatures[idx],digits=2))K → $(round(temperatures[idx+1],digits=2))K, " *
                "η_diff=$(round(eta_diff_hist[idx],digits=4)) → $(round(eta_diff_hist[idx+1],digits=4))")
    end
else
    println("\n未发现温度下降时间段")
end

# 绘图
println("\n生成诊断图表...")

# 图1：温度与产热/散热对比
p1 = plot(layout=(3,1), size=(1000, 900), margin=5Plots.mm)
plot!(p1[1], times, temperatures, label="Temperature", lw=2, ylabel="T [K]", legend=:topleft)
plot!(p1[2], times, Q_total_hist, label="Q_gen", lw=2, ylabel="Power [W]", legend=:topleft)
plot!(p1[2], times, Q_conv_hist, label="Q_conv", lw=2, linestyle=:dash)
plot!(p1[2], times, Q_total_hist .- Q_conv_hist, label="Q_net", lw=2, linestyle=:dot)
plot!(p1[3], times, dT_dt_hist .* 1000, label="dT/dt", lw=2, ylabel="dT/dt [mK/s]", xlabel="Time [s]", legend=:topleft)
hline!(p1[3], [0], label="", linestyle=:dash, color=:black, alpha=0.3)
savefig(p1, "diagnose_temperature_heat.png")

# 图2：过电位分析
p2 = plot(layout=(2,1), size=(1000, 600), margin=5Plots.mm)
plot!(p2[1], times, eta_p_hist, label="η_p", lw=2, ylabel="Overpotential", legend=:topleft)
plot!(p2[1], times, eta_n_hist, label="η_n", lw=2)
plot!(p2[2], times, eta_diff_hist, label="η_p - η_n", lw=2, ylabel="Overpotential Diff", xlabel="Time [s]", legend=:topleft, color=:red)
savefig(p2, "diagnose_overpotential.png")

# 图3：电压-温度相关性
p3 = scatter(voltages, temperatures, xlabel="Voltage [V]", ylabel="Temperature [K]", 
             marker=:circle, markersize=3, label="", title="T vs V trajectory")
savefig(p3, "diagnose_voltage_temperature.png")

# 图4：过电位-产热相关性
p4 = scatter(eta_diff_hist, Q_total_hist, xlabel="η_p - η_n", ylabel="Q_gen [W]",
             marker=:circle, markersize=3, label="", title="Heat vs Overpotential")
savefig(p4, "diagnose_overpotential_heat.png")

println("✓ 图表已保存:")
println("  - diagnose_temperature_heat.png")
println("  - diagnose_overpotential.png")
println("  - diagnose_voltage_temperature.png")
println("  - diagnose_overpotential_heat.png")

# 导出CSV数据
df = DataFrame(
    time_s = times,
    voltage_V = voltages,
    temperature_K = temperatures,
    eta_n = eta_n_hist,
    eta_p = eta_p_hist,
    eta_diff = eta_diff_hist,
    Q_rxn = Q_rxn_hist,
    Q_rev = Q_rev_hist,
    Q_total_W = Q_total_hist,
    Q_conv_W = Q_conv_hist,
    dT_dt_Ks = dT_dt_hist,
    current_A = current_hist
)
CSV.write("diagnose_data.csv", df)
println("✓ 数据已导出: diagnose_data.csv")

println("\n" * "=" ^70)
println("诊断完成！")
println("=" ^70)
