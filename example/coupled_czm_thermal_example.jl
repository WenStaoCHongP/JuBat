# coupled_czm_thermal_example.jl
# 
# 电化学-热-内聚力多物理场耦合仿真示例
# 
# 特性:
# - 使用双网格生成（热网格合并节点，CZM网格保留重合）
# - 间隙导热模型（损伤相关的界面热阻）
# - 失效单元电流重分配
# - SOH监控和自动终止
# - Mode I only 内聚力模型
#
# 作者: JuBat Team
# 日期: 2024

using Pkg
Pkg.activate(".")
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Statistics
using Printf
using Plots

println("="^60)
println("电化学-热-内聚力多物理场耦合仿真")
println("="^60)

# ========================================================================
# 1. 加载参数
# ========================================================================
println("\n[1] 加载参数...")

param_dim = JuBat.ChooseCell("Jellyroll")
println("  间隙导热参数:")
println("    h_c0 = $(param_dim.cohesive.h_c0) W/(m^2 K)")
println("    k_air = $(param_dim.cohesive.k_air) W/(m K)")
println("  内聚力参数:")
# println("    sigma_max = $(param_dim.cohesive.σ_max_n/1e6) MPa")  # TODO Chunk 2 Task 2.1
# println("    delta_0 = $(param_dim.cohesive.δ_0_n*1e6) um")      # TODO Chunk 2 Task 2.1
# println("    delta_c = $(param_dim.cohesive.δ_c_n*1e6) um")      # TODO Chunk 2 Task 2.1

# ========================================================================
# 2. 创建仿真选项
# ========================================================================
println("\n[2] 创建仿真选项...")

opt = JuBat.Option()
opt.model = "SPMe"
opt.dtType = "auto"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.cool_method = "tab"
opt.per_element_spme = true

# CZM选项
opt.czm_enabled = true
opt.czm_model = "model1"              # 只考虑法向脱粘
opt.czm_update_interval = 1             # 每步更新损伤
opt.czm_soh_threshold = 0.8             # SOH终止阈值 80%
opt.czm_inner_exit_only = true          # 仅内圈单元在断裂时退出
opt.czm_fix_inner = false               # 力学边界：内圈自由，外圈固定
opt.czm_iter_method = "basic"     # 使用弧长法处理后峰软化（配合粘性正则化）
opt.czm_load_steps = 10                 # 载荷子步数
opt.debug_coupling = false               # 打印 CZM 每步诊断
opt.czm_area_loss_enabled = false        # 启用渐进式有效面积损失
opt.czm_area_loss_threshold = 0.83      # 面积开始缩减的 D 阈值
# 粘性正则化（推荐配合 load_substep 或 arc_length 使用）
opt.czm_viscous_enabled = false         # 关闭粘性正则化（回归测试）
opt.czm_visc_tau = 0                  # 松弛时间

println("  CZM选项:")
println("    czm_enabled = $(opt.czm_enabled)")
println("    czm_model = $(opt.czm_model)")
println("    czm_soh_threshold = $(opt.czm_soh_threshold)")
println("    czm_iter_method = $(opt.czm_iter_method)")
println("    czm_viscous_enabled = $(opt.czm_viscous_enabled)")
println("    czm_fix_inner = $(opt.czm_fix_inner) (边界条件: 内圈$(opt.czm_fix_inner ? "固定" : "自由"), 外圈固定)")
println("    czm_area_loss_enabled = $(opt.czm_area_loss_enabled)")
println("    czm_area_loss_threshold = $(opt.czm_area_loss_threshold)")
if opt.czm_viscous_enabled
    println("    czm_visc_tau = $(opt.czm_visc_tau) s")
end

# ========================================================================
# 3. 生成双网格
# ========================================================================
println("\n[3] 生成统一网格...")

# 先创建Case以获取归一化参数
case_temp = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case_temp.param; nθ=80, czm_enabled=true, gsorder=2)

println("  热网格单元数: $(mesh_data.ne)")
println("  热网格节点数: $(mesh_data.nnode)")
println("  界面节点对数: $(length(mesh_data.interface_pairs))")
println("  内圈单元数: $(sum(mesh_data.is_inner_layer))")
println("  外圈单元数: $(sum(.!mesh_data.is_inner_layer))")

# ========================================================================
# 4. 初始化Case
# ========================================================================
println("\n[4] 初始化Case...")

# 创建Case并安装活动热网格
case = JuBat.SetCase(param_dim, opt)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
mesh_th = case.mesh["thermal2D"]

# mesh_data 保持为局部变量用于后续处理

# ========================================================================
# 5. 创建CZM网格
# ========================================================================
println("\n[5] 创建CZM网格...")

czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_th, case.param)

println("  内聚力单元数: $(czm_mesh.n_cohesive)")
println("  层数: $(czm_mesh.n_layers)")

# ========================================================================
# 6. 设置循环参数
# ========================================================================
println("\n[6] 设置循环参数...")

cycle_opt = JuBat.CycleOption(
    n_cycles = 1,           # 循环次数
    SOC_init = 0.65,         # 初始SOC 90%
    t_discharge = 1800.0,   # 放电时间 1小时
    t_charge = 1800.0,      # 充电时间 1小时
    t_rest1 = 600.0,        # 静置1时间 5分钟
    t_rest2 = 600.0,        # 静置2时间 5分钟
    I_discharge = 5.0,      # 放电电流 1C (5A)
    I_charge = 5.0,         # 充电电流 1C (5A)
    V_lower = 2.5,          # 放电截止电压
    V_upper = 4.2,          # 充电截止电压
    dt_cycle = [0.1, 5.0], # 时间步范围
    reset_T_each_cycle = false
)

println("  循环次数: $(cycle_opt.n_cycles)")
println("  初始SOC: $(cycle_opt.SOC_init * 100)%")
println("  放电电流: $(cycle_opt.I_discharge)A ($(cycle_opt.I_discharge/5.0)C)")
println("  充电电流: $(cycle_opt.I_charge)A")

# ========================================================================
# 7. 运行循环仿真
# ========================================================================
println("\n[7] 运行循环仿真...")

result = JuBat.solve_cycling(case, cycle_opt, czm_mesh; verbose=true, save_detailed=true)

# ========================================================================
# 7.5. CSV 导出
# ========================================================================
println("\n[7.5] 导出CSV文件...")

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)
csv_dir = joinpath(output_dir, "csv", "czm_study_2")

# 配置CSV导出选项：仅输出每个阶段首尾步，指定循环完整输出
csv_opt = JuBat.CsvExportOptions(
    :phase_ends,                # 仅阶段首尾
    1,                          # save_every (仅 :custom 模式使用)
    [1],       # 第1个和最后一个循环完整输出
    String[]                    # 不跳过任何文件
)

JuBat.export_cycling_csv(result, case, czm_mesh;
                          output_dir=csv_dir, overwrite=true,
                          csv_opt=csv_opt)

# ========================================================================
# 8. 结果分析
# ========================================================================
println("\n" * "="^60)
println("结果分析")
println("="^60)

println("\n容量衰减:")
for i in 1:result.n_cycles
    @printf("  循环%d: 放电容量 %.3fAh, SOH %.1f%%\n", 
            i, result.capacity_discharge[i], result.soh[i] * 100)
end

println("\n损伤演化:")
for i in 1:result.n_cycles
    @printf("  循环%d: D_max %.2f%%, D_mean %.2f%%, 断裂单元 %d\n",
            i, result.D_max[i] * 100, result.D_mean[i] * 100, result.n_fractured[i])
end

println("\n温度:")
for i in 1:result.n_cycles
    @printf("  循环%d: T_max %.2fK\n", i, result.T_max[i])
end

if result.soh_terminated
    println("\n[!] 仿真因SOH低于阈值而终止")
end

# ========================================================================
# 9. 绘图
# ========================================================================
println("\n[8] 绘图...")

# 容量和SOH
p1 = plot(result.cycle_idx, result.capacity_discharge,
          xlabel="Cycle", ylabel="Discharge Capacity (Ah)",
          label="Capacity", marker=:circle, lw=2)

p2 = plot(result.cycle_idx, result.soh .* 100,
          xlabel="Cycle", ylabel="SOH (%)",
          label="SOH", marker=:square, lw=2, color=:red)

# 损伤
p3 = plot(result.cycle_idx, result.D_max .* 100,
          xlabel="Cycle", ylabel="Damage (%)",
          label="D_max", marker=:diamond, lw=2)
plot!(p3, result.cycle_idx, result.D_mean .* 100,
      label="D_mean", marker=:circle, lw=2, ls=:dash)

# 温度
p4 = plot(result.cycle_idx, result.T_max,
          xlabel="Cycle", ylabel="Max Temperature (K)",
          label="T_max", marker=:star, lw=2, color=:orange)

# 组合图
p_combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600))

# 保存
savefig(p_combined, joinpath(output_dir, "coupled_czm_thermal_results.png"))
println("  结果图已保存到: $(joinpath(output_dir, "coupled_czm_thermal_results.png"))")

# ========================================================================
# 10. 间隙导热分析
# ========================================================================
println("\n[9] 间隙导热分析...")

# 计算所有内聚力单元的有效换热系数（按界面类型取 per-interface 参数）
czm_param_cache = case.czm_param_cache !== nothing ? case.czm_param_cache :
    JuBat.compute_czm_params_per_interface(case)
h_eff_all = [JuBat.compute_element_gap_conductance(czm_mesh, i,
    czm_param_cache.by_interface[czm_mesh.cohesive_elements[i].interface_type])
    for i in 1:czm_mesh.n_cohesive]

println("  有效换热系数统计:")
println("    最小值: $(minimum(h_eff_all)) W/(m^2 K)")
println("    最大值: $(maximum(h_eff_all)) W/(m^2 K)")
println("    平均值: $(mean(h_eff_all)) W/(m^2 K)")

# 断裂单元
fractured = JuBat.get_fractured_elements(czm_mesh)
println("  断裂单元数: $(length(fractured))")

# 活跃热单元
active = JuBat.get_active_elements(czm_mesh, case.geometry)
println("  活跃热单元数: $(length(active)) / $(length(case.geometry.element_layer))")

println("\n" * "="^60)
println("仿真完成")
println("="^60)

