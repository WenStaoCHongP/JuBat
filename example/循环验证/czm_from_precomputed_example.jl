"""
从预计算数据进行CZM仿真示例

功能：
- 从CSV文件读取预计算的温度和SOC数据
- 基于这些数据计算CZM损伤演化
- 与全耦合仿真结果进行对比验证

优势：
- 大幅加快仿真速度（跳过电化学-热求解）
- 可以多次运行不同的CZM参数而不需要重新计算电化学-热
- 便于参数研究和敏感性分析

使用方法：
1. 先运行 export_cycle_data_example.jl 生成数据
2. 然后运行本文件进行CZM仿真

日期：2025
"""

using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

"""
    compute_czm_from_precomputed(cycle_data, param_dim, czm_mesh; verbose=true, n_cycles=1)

从预计算数据计算CZM损伤演化（可重复多循环）。

# 参数
- `cycle_data`: 从CSV加载的循环数据字典
- `param_dim`: 电池维度参数
- `czm_mesh`: CZM网格对象
- `n_cycles`: 循环次数（重复使用相同温度/SOC曲线）

# 返回
- 损伤历史数据
"""
function compute_czm_from_precomputed(cycle_data::Dict, param_dim, czm_mesh; verbose::Bool=true, n_cycles::Int=1, snapshot_cycles::Vector{Int}=Int[])
    n_steps = cycle_data["n_steps"]
    times = cycle_data["times"]
    phases = cycle_data["phases"]
    T_nodes = cycle_data["T_nodes"]  # n_steps × nT
    soc_n = cycle_data["soc_n"]      # n_steps × ne
    soc_p = cycle_data["soc_p"]      # n_steps × ne
    
    ne = size(soc_n, 2)
    nT = size(T_nodes, 2)
    
    if verbose
        println("="^60)
        println("从预计算数据进行CZM仿真")
        println("="^60)
        @printf("  时间步数: %d\n", n_steps)
        @printf("  单元数: %d, 节点数: %d\n", ne, nT)
        @printf("  CZM单元数: %d\n", czm_mesh.n_cohesive)
        @printf("  循环次数: %d\n", n_cycles)
    end
    
    # 获取CZM参数
    cohesive = param_dim.cohesive
    
    # 计算有效材料参数
    E_eff = (param_dim.NE.E * param_dim.NE.thickness + param_dim.PE.E * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    ν_eff = (param_dim.NE.nu * param_dim.NE.thickness + param_dim.PE.nu * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    α_eff = (param_dim.NE.alphaT * param_dim.NE.thickness + param_dim.PE.alphaT * param_dim.PE.thickness) / 
            (param_dim.NE.thickness + param_dim.PE.thickness)
    β_n = param_dim.NE.Omega / 3.0
    β_p = param_dim.PE.Omega / 3.0
    
    # 参考温度和SOC
    T_ref = param_dim.cell.T0
    soc_ref_n = param_dim.NE.cs0 / param_dim.NE.cs_max
    soc_ref_p = param_dim.PE.cs0 / param_dim.PE.cs_max
    
    # 初始化CZM位移场
    ndof_czm = 2 * czm_mesh.nnode
    u_czm = zeros(Float64, ndof_czm)
    
    # 损伤历史记录（循环维度）
    D_max_cycle = zeros(Float64, n_cycles)
    D_mean_cycle = zeros(Float64, n_cycles)
    n_fractured_cycle = zeros(Int, n_cycles)

    # 记录阶段内（用于诊断）
    D_max_hist = zeros(Float64, n_steps)
    D_mean_hist = zeros(Float64, n_steps)
    n_fractured_hist = zeros(Int, n_steps)

    # 最大平均温度时刻（跨所有循环）
    max_T_mean = -Inf
    T_nodes_at_max_Tmean = nothing
    time_at_max_Tmean = 0.0
    cycle_at_max_Tmean = 1

    # 裂纹（损伤）快照
    damage_snapshots = Dict{Int, Vector{Float64}}()
    
    # 更新间隔
    update_interval = 10
    
    if verbose
        println("\n开始CZM计算...")
    end
    
    # 时间步循环（多循环重复同一曲线）
    for cyc in 1:n_cycles
        if verbose
            println("\n--- 循环 $cyc / $n_cycles ---")
        end

        for step in 1:n_steps
            # 提取当前时间步的温度和SOC
            T_nodes_step = T_nodes[step, :]
            soc_n_step = soc_n[step, :]
            soc_p_step = soc_p[step, :]
        
            # 计算温度变化和SOC变化
            dT_elem = zeros(Float64, ne)
            Δsoc_n_elem = zeros(Float64, ne)
            Δsoc_p_elem = zeros(Float64, ne)
        
            for e in 1:ne
                # 计算单元平均温度变化
                nodes = czm_mesh.bulk_element[e, :]
                T_elem = 0.0
                valid_nodes = 0
                for n in nodes
                    if n <= nT
                        T_elem += T_nodes_step[n]
                        valid_nodes += 1
                    end
                end
                if valid_nodes > 0
                    T_elem /= valid_nodes
                    dT_elem[e] = T_elem - T_ref
                end
                
                # SOC变化
                Δsoc_n_elem[e] = soc_n_step[e] - soc_ref_n
                Δsoc_p_elem[e] = soc_p_step[e] - soc_ref_p
            end
        
            # 每隔一定步数更新CZM
            if step % update_interval == 0 || step == 1 || step == n_steps
                # 外力向量
                F_ext = zeros(Float64, ndof_czm)
                
                # CZM求解
                try
                    result = JuBat.solve_czm_step(
                        czm_mesh, F_ext, E_eff, ν_eff, cohesive, param_dim, u_czm;
                        α_eff=α_eff, β_n=β_n, β_p=β_p,
                        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
                        max_iter=30, tol=1e-6
                    )
                    u_czm = result.displacement
                catch e
                    if verbose
                        @printf("  WARN: 循环 %d 步 %d CZM求解失败: %s\n", cyc, step, string(e))
                    end
                end
            end
        
            # 记录损伤统计
            stats = JuBat.get_damage_statistics(czm_mesh)
            D_max_hist[step] = stats.max_D
            D_mean_hist[step] = stats.mean_D
            n_fractured_hist[step] = stats.n_fractured

            # 记录最大平均温度时刻
            T_mean = mean(T_nodes_step)
            if T_mean > max_T_mean
                max_T_mean = T_mean
                T_nodes_at_max_Tmean = copy(T_nodes_step)
                time_at_max_Tmean = times[step]
                cycle_at_max_Tmean = cyc
            end
            
            # 进度输出
            if verbose && (step % 50 == 0 || step == n_steps)
                @printf("  步 %d/%d: t=%.1fs, D_max=%.4f%%, D_mean=%.4f%%\n",
                        step, n_steps, times[step], D_max_hist[step]*100, D_mean_hist[step]*100)
            end
        end

        # 每个循环结束，记录循环级损伤
        D_max_cycle[cyc] = D_max_hist[end]
        D_mean_cycle[cyc] = D_mean_hist[end]
        n_fractured_cycle[cyc] = n_fractured_hist[end]

        if cyc in snapshot_cycles
            damage_snapshots[cyc] = [s.D for s in czm_mesh.damage_states]
        end
    end
    
    if verbose
        println("\nCZM计算完成")
        @printf("  最终最大损伤: %.4f%%\n", D_max_cycle[end] * 100)
        @printf("  最终平均损伤: %.4f%%\n", D_mean_cycle[end] * 100)
        @printf("  断裂单元数: %d / %d\n", n_fractured_cycle[end], czm_mesh.n_cohesive)
    end
    
    return Dict(
        "times" => times,
        "phases" => phases,
        "D_max" => D_max_hist,
        "D_mean" => D_mean_hist,
        "n_fractured" => n_fractured_hist,
        "D_max_cycle" => D_max_cycle,
        "D_mean_cycle" => D_mean_cycle,
        "n_fractured_cycle" => n_fractured_cycle,
        "n_cycles" => n_cycles,
        "T_nodes_at_max_Tmean" => T_nodes_at_max_Tmean,
        "max_T_mean" => max_T_mean,
        "time_at_max_Tmean" => time_at_max_Tmean,
        "cycle_at_max_Tmean" => cycle_at_max_Tmean,
        "damage_snapshots" => damage_snapshots
    )
end


function main()
    println("="^70)
    println("从预计算数据进行CZM仿真")
    println("="^70)
    
    # ========================================================================
    # 1. 加载预计算数据
    # ========================================================================
    println("\n[1/4] 加载预计算数据...")
    
    data_dir = joinpath(@__DIR__, "..", "output", "cycle_data")
    
    if !isdir(data_dir)
        println("\n  ERROR: 数据目录不存在: $data_dir")
        println("  请先运行 export_cycle_data_example.jl 生成数据")
        return nothing
    end
    
    cycle_data = nothing
    try
        cycle_data = JuBat.load_cycle_data_from_csv(data_dir; prefix="cycle")
        println("  OK: 数据加载成功")
        @printf("  时间步数: %d\n", cycle_data["n_steps"])
        @printf("  节点数: %d\n", cycle_data["nT"])
        @printf("  单元数: %d\n", cycle_data["ne"])
    catch e
        println("\n  ERROR: 数据加载失败: $e")
        println("  请先运行 export_cycle_data_example.jl 生成数据")
        return nothing
    end
    
    # ========================================================================
    # 2. 创建CZM网格
    # ========================================================================
    println("\n[2/4] 创建CZM网格...")
    
    # 电池参数
    param_dim = JuBat.ChooseCell("Jellyroll")
    
    # 打印内聚力参数
    println("  内聚力参数：")
        # @printf("    sigma_max_n = %.1f MPa, delta_c_n = %.4f um\n",                                   # TODO Chunk 2 Task 2.1
        #     param_dim.cohesive.σ_max_n / 1e6, param_dim.cohesive.δ_c_n * 1e6)                           # TODO Chunk 2 Task 2.1
        # @printf("    tau_max_t = %.1f MPa, delta_c_t = %.4f um\n",                                      # TODO Chunk 2 Task 2.1
        #     param_dim.cohesive.τ_max_t / 1e6, param_dim.cohesive.δ_c_t * 1e6)                           # TODO Chunk 2 Task 2.1
    
    # 从加载的数据重建网格
    nT = cycle_data["nT"]
    ne = cycle_data["ne"]
    node_coords = cycle_data["node_coords"]
    element_connectivity = cycle_data["element_connectivity"]
    
    # 创建热网格结构（用于CZM）
    mesh_th = (
        node = node_coords,
        element = element_connectivity,
        nlen = nT
    )
    
    # 注意：由于CZM网格创建需要完整的mesh对象，我们需要创建一个完整的mesh
    # 这里我们重新创建热网格来获取正确的CZM网格
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    # 裂纹模式选择："model1" 仅法向；"mix" 混合模式
    opt.czm_model = "model1"
    @printf("  CZM 模式: %s\n", opt.czm_model)
    
    case = JuBat.SetCase(param_dim, opt)
    n_theta = 60
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, nθ_czm=80, gsorder=2)
    JuBat.setup_thermal2D_mesh!(case, mesh_data)
    mesh_th_full = case.mesh["thermal2D"]

    # 创建CZM网格
    czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, mesh_data.thermal2D, case.param)
    
    println("  OK: CZM网格创建完成")
    @printf("  内聚力单元数: %d\n", czm_mesh.n_cohesive)
    
    # ========================================================================
    # 3. 运行CZM计算
    # ========================================================================
    println("\n[3/4] 运行CZM计算...")
    
    # 多循环设置（重复使用同一温度/SOC曲线）
    n_cycles = 10
    snapshot_cycles = [1, 5, 10]
    damage_result = compute_czm_from_precomputed(
        cycle_data, param_dim, czm_mesh; verbose=true, n_cycles=n_cycles, snapshot_cycles=snapshot_cycles
    )
    
    # ========================================================================
    # 4. 结果可视化
    # ========================================================================
    println("\n[4/4] 生成结果图像...")
    
    # 确保输出目录存在
    output_dir = joinpath(@__DIR__, "..", "output")
    isdir(output_dir) || mkdir(output_dir)
    
    times = damage_result["times"]
    D_max = damage_result["D_max"]
    D_mean = damage_result["D_mean"]
    D_max_cycle = damage_result["D_max_cycle"]
    D_mean_cycle = damage_result["D_mean_cycle"]
    n_cycles = damage_result["n_cycles"]
    
    # 图1: 损伤演化
    p1 = plot(times, D_max .* 100,
              xlabel="Time (s)", ylabel="Damage (%)",
              label="D_max", linewidth=2, color=:red,
              title="CZM Damage Evolution (from precomputed data)")
    plot!(p1, times, D_mean .* 100,
          label="D_mean", linewidth=2, color=:blue, linestyle=:dash)
    
    savefig(p1, joinpath(output_dir, "czm_precomputed_damage.png"))
    println("  OK: 保存 output/czm_precomputed_damage.png")
    
    # 图2: 分阶段损伤
    phases = damage_result["phases"]
    
    # 找出各阶段的时间范围
    discharge_idx = findall(x -> x == "discharge", phases)
    rest1_idx = findall(x -> x == "rest", phases)
    charge_idx = findall(x -> x == "charge", phases)
    
    p2 = plot(size=(900, 600), layout=(1, 1))
    
    if !isempty(discharge_idx)
        plot!(p2, times[discharge_idx], D_max[discharge_idx] .* 100,
              label="Discharge", linewidth=2, color=:red)
    end
    if !isempty(rest1_idx)
        # 区分静置1和静置2
        rest_times = times[rest1_idx]
        rest_D = D_max[rest1_idx]
        if length(rest_times) > 1
            mid_idx = div(length(rest_times), 2)
            plot!(p2, rest_times[1:mid_idx], rest_D[1:mid_idx] .* 100,
                  label="Rest 1", linewidth=2, color=:gray)
            if mid_idx < length(rest_times)
                plot!(p2, rest_times[mid_idx+1:end], rest_D[mid_idx+1:end] .* 100,
                      label="Rest 2", linewidth=2, color=:gray, linestyle=:dash)
            end
        else
            plot!(p2, rest_times, rest_D .* 100,
                  label="Rest", linewidth=2, color=:gray)
        end
    end
    if !isempty(charge_idx)
        plot!(p2, times[charge_idx], D_max[charge_idx] .* 100,
              label="Charge", linewidth=2, color=:blue)
    end
    
    plot!(p2, xlabel="Time (s)", ylabel="D_max (%)",
          title="Damage Evolution by Phase", legend=:topright)
    
    savefig(p2, joinpath(output_dir, "czm_precomputed_damage_by_phase.png"))
    println("  OK: 保存 output/czm_precomputed_damage_by_phase.png")

    # 图3: 多循环损伤演化（每循环一次）
    cycles = 1:n_cycles
    p3 = plot(cycles, D_max_cycle .* 100,
              xlabel="Cycle Number", ylabel="Damage (%)",
              label="D_max", linewidth=2, color=:red,
              title="CZM Damage Evolution (per cycle)")
    plot!(p3, cycles, D_mean_cycle .* 100,
          label="D_mean", linewidth=2, color=:blue, linestyle=:dash)
    savefig(p3, joinpath(output_dir, "czm_precomputed_damage_by_cycle.png"))
    println("  OK: 保存 output/czm_precomputed_damage_by_cycle.png")

    # 图4: 最大平均温度对应的二维温度场
    T_nodes_at_max_Tmean = damage_result["T_nodes_at_max_Tmean"]
    if T_nodes_at_max_Tmean !== nothing
        node_coords = cycle_data["node_coords"]
        pT = scatter(node_coords[:, 1], node_coords[:, 2],
                     marker_z=T_nodes_at_max_Tmean,
                     xlabel="x (m)", ylabel="y (m)",
                     title="Temperature Field at Max Mean T (cycle $(damage_result["cycle_at_max_Tmean"]))",
                     colorbar_title="T (K)", markersize=4)
        savefig(pT, joinpath(output_dir, "czm_max_mean_T_field.png"))
        println("  OK: 保存 output/czm_max_mean_T_field.png")
    end

    # 图5: 裂纹（损伤）演化（循环 1/5/10）
    # 计算内聚力单元中心
    function cohesive_centers(czm_mesh)
        centers = zeros(Float64, czm_mesh.n_cohesive, 2)
        for (i, elem) in enumerate(czm_mesh.cohesive_elements)
            n = elem.nodes
            x = mean(czm_mesh.node[n, 1])
            y = mean(czm_mesh.node[n, 2])
            centers[i, 1] = x
            centers[i, 2] = y
        end
        return centers
    end

    centers = cohesive_centers(czm_mesh)
    damage_snapshots = damage_result["damage_snapshots"]
    for c in snapshot_cycles
        D_snapshot = get(damage_snapshots, c, nothing)
        D_snapshot === nothing && continue
        pD = scatter(centers[:, 1], centers[:, 2],
                     marker_z=D_snapshot,
                     xlabel="x (m)", ylabel="y (m)",
                     title="Damage Map after Cycle $c",
                     colorbar_title="D", markersize=4)
        savefig(pD, joinpath(output_dir, "czm_damage_cycle_$(c).png"))
        println("  OK: 保存 output/czm_damage_cycle_$(c).png")
    end
    
    # ========================================================================
    # 结果汇总
    # ========================================================================
    println("\n" * "="^70)
    println("CZM仿真结果汇总")
    println("="^70)
    
    @printf("\n  总时间步数: %d\n", length(times))
    @printf("  仿真时间: %.1f s → %.1f s\n", times[1], times[end])
    @printf("  最终最大损伤: %.4f%%\n", D_max_cycle[end] * 100)
    @printf("  最终平均损伤: %.4f%%\n", D_mean_cycle[end] * 100)
    @printf("  断裂单元数: %d / %d\n", damage_result["n_fractured_cycle"][end], czm_mesh.n_cohesive)
    
    println("\n  输出文件:")
    println("    1. output/czm_precomputed_damage.png - 损伤演化曲线")
    println("    2. output/czm_precomputed_damage_by_phase.png - 分阶段损伤")
    println("    3. output/czm_precomputed_damage_by_cycle.png - 循环损伤演化")
    println("    4. output/czm_max_mean_T_field.png - 最大平均温度场")
    println("    5. output/czm_damage_cycle_1.png / _5.png / _10.png - 裂纹演化")
    
    println("\n" * "="^70)
    println("验证说明")
    println("="^70)
    println("\n  此结果可与全耦合仿真(czm_cycle_example.jl)进行对比。")
    println("  由于预计算方法跳过了电化学-热求解，计算速度大幅提升。")
    println("  对于多循环仿真或参数研究，建议使用此方法。")
    
    return damage_result, czm_mesh
end

# 运行主函数
if abspath(PROGRAM_FILE) == @__FILE__
    damage_result, czm_mesh = main()
end