"""
Jellyroll SPMe-热耦合 vs 圆环热模型对比验证

验证流程：
1. 运行 Jellyroll SPMe-二维分布式热耦合仿真，提取时间依赖的平均内热源
2. 将平均热源输入到简化的圆环热模型（ring2D）
3. 对比两种模型的温度演化曲线，验证热模型的正确性

日期：2026-03-14
"""

using LinearAlgebra, Statistics, Plots, Printf

include(joinpath(@__DIR__, "..", "..", "src", "JuBat.jl"))
using .JuBat

const SEP = "="^80

# ============================================================================
# 辅助函数
# ============================================================================

"""线性插值"""
function interp_linear(x::Vector{Float64}, y::Vector{Float64}, xi::Vector{Float64})
    yi = similar(xi)
    j = 1
    @inbounds for k in eachindex(xi)
        xk = xi[k]
        while j < length(x) - 1 && x[j + 1] < xk
            j += 1
        end
        x1, x2 = x[j], x[j + 1]
        y1, y2 = y[j], y[j + 1]
        yi[k] = abs(x2 - x1) < 1e-15 ? y1 : y1 + (xk - x1) / (x2 - x1) * (y2 - y1)
    end
    return yi
end

"""计算单元面积"""
function compute_element_areas(mesh)
    A_elem = zeros(Float64, size(mesh.element, 1))
    @inbounds for g in 1:length(mesh.gs.detJ)
        e = mesh.gs.ele[g]
        A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    return A_elem
end

"""计算单元平均值（从节点值）"""
function compute_element_mean(mesh, nodal_values)
    ne, nn_per_elem = size(mesh.element)
    elem_mean = zeros(Float64, ne)
    @inbounds for e in 1:ne
        s = zero(eltype(nodal_values))
        for i in 1:nn_per_elem
            s += nodal_values[mesh.element[e, i]]
        end
        elem_mean[e] = s / nn_per_elem
    end
    return elem_mean
end

"""面积加权平均"""
area_weighted_mean(values, areas) = sum(values .* areas) / sum(areas)

"""单元质心半径 (mesh.node 为 x,y; 返回每单元 r = sqrt(cx^2+cy^2))"""
function element_centroid_radii(mesh)
    ne = size(mesh.element, 1)
    r_c = zeros(Float64, ne)
    for e in 1:ne
        nodes = mesh.element[e, :]
        cx = sum(mesh.node[n, 1] for n in nodes) / length(nodes)
        cy = sum(mesh.node[n, 2] for n in nodes) / length(nodes)
        r_c[e] = sqrt(cx^2 + cy^2)
    end
    return r_c
end

# ============================================================================
# 第一部分：Jellyroll SPMe-热耦合仿真
# ============================================================================

function run_jellyroll_simulation()
    println(SEP)
    println("第一部分：Jellyroll SPMe-二维分布式热耦合仿真")
    println(SEP)

    # 参数与选项配置
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l, param_dim.cell.v_h = 2.5, 4.2

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn, opt.Ns, opt.Np = 10, 5, 10
    opt.Nrn, opt.Nrp = 10, 10
    opt.gsorder, opt.dimension = 2, 1
    opt.mechanicalmodel = "none"
    opt.Current = _ -> 5.0  # 1C 放电
    opt.time, opt.dt = [0.0, 60.0], [0.5, 10.0]
    opt.dtType, opt.jacobi, opt.solveType = "auto", "update", "Crank-Nicolson"
    opt.thermal_enabled = true
    opt.thermalmodel, opt.thermal_dim = "distributed2D", "2D"
    opt.cool_method = "none"
    opt.per_element_spme, opt.debug_coupling, opt.czm_enabled = true, false, false

    # 创建案例和网格
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh = case.mesh["thermal2D"]

    println("Jellyroll网格创建完成: $(size(mesh.element, 1)) 个单元")
    println("开始求解...")
    result = JuBat.Solve(case)

    # 提取结果（量纲：t [s]，T [K]，q_nd 无量纲，A_elem 无量纲面积 L²_ref）
    t = result["time [s]"]
    T_nodes_hist = result["thermal2D T_nodes history [K]"]
    q_nd_hist = result["heat_source_fields"]

    A_elem = compute_element_areas(mesh)

    T_vol = [area_weighted_mean(compute_element_mean(mesh, T_nodes_hist[:, k]), A_elem) for k in eachindex(t)]
    q_mean_nd = [area_weighted_mean(q_nd_hist[:, k], A_elem) for k in eachindex(t)]
    q_mean_phys = q_mean_nd .* case.param_dim.scale.q

    println("Jellyroll仿真完成")
    @printf("  时间步数: %d, 温升: %.2f K\n", length(t), T_vol[end] - T_vol[1])
    @printf("  温度: %.2f → %.2f K\n", T_vol[1], T_vol[end])
    @printf("  热源范围: [%.2e, %.2e] W/m^3\n", minimum(q_mean_phys), maximum(q_mean_phys))

    return (; t, T_vol, q_mean_phys, q_mean_nd, q_nd_hist, A_elem, A_total=sum(A_elem), mesh, case, param_dim, q_ref=case.param_dim.scale.q)
end

# ============================================================================
# 第二部分：圆环热模型仿真
# ============================================================================

function run_ring_simulation(jellyroll_result)
    println("\n$SEP")
    println("第二部分：圆环热模型仿真")
    println(SEP)

    (; t, q_mean_phys, param_dim, mesh, case, A_elem, A_total, q_nd_hist) = jellyroll_result

    # Ring.jl 已与 Jellyroll 归一化一致
    param_ring = JuBat.ChooseCell("Ring")
    for prop in [:Rin, :Rout, :h, :T0, :T_amb, :lambda_r, :lambda_t]
        setproperty!(param_ring.cell, prop, getproperty(param_dim.cell, prop))
    end

    opt = JuBat.Option()
    opt.model, opt.thermal_enabled, opt.thermalmodel = "thermal", true, "ring2D_polar"
    opt.time, opt.dt, opt.solveType = [t[1], t[end]], [0.2, 2.0], "backward"

    case_ring = JuBat.SetCase(param_ring, opt)
    q_ref_ring = param_ring.scale.q

    # 创建圆环网格
    mesh_data = JuBat.ring_mesh(case_ring.param; ntheta=40, nr=20, gsorder=2)
    mesh_ring = mesh_data.mesh
    case_ring.mesh["thermal2D"] = mesh_ring

    # 总功率桥接：使 Ring 总热输入与 Jellyroll 一致
    A_elem_ring = compute_element_areas(mesh_ring)
    A_total_ring = sum(A_elem_ring)
    ne = size(mesh_ring.element, 1)
    q_mean_phys_ring = q_mean_phys .* A_total ./ A_total_ring

    # 热源径向分布：用 Jellyroll 单元质心半径 r 分箱
    r_centroid_jelly = element_centroid_radii(mesh)
    nr, ntheta = mesh_data.nr, mesh_data.ntheta
    r_edges = mesh_data.r
    ne_jelly = size(mesh.element, 1)
    nstep = length(t)
    q_bin_nd = fill(NaN64, nr, nstep)
    for b in 1:nr
        in_bin = [e for e in 1:ne_jelly if (r_edges[b] <= r_centroid_jelly[e] < r_edges[b+1])]
        for k in 1:nstep
            if isempty(in_bin)
                q_bin_nd[b, k] = q_mean_phys_ring[k] / q_ref_ring
            else
                q_bin_nd[b, k] = sum(q_nd_hist[e, k] * A_elem[e] for e in in_bin) / sum(A_elem[e] for e in in_bin)
            end
        end
    end
    # 空 bin 用相邻或全局均值填充
    for k in 1:nstep
        for b in 1:nr
            if isnan(q_bin_nd[b, k])
                q_bin_nd[b, k] = q_mean_phys_ring[k] / q_ref_ring
            end
        end
    end
    q_mean_nd_ring = q_mean_phys_ring ./ q_ref_ring
    q_ring_nd = zeros(Float64, ne, nstep)
    for k in 1:nstep
        mean_ring_k = sum(q_bin_nd[div(e - 1, ntheta) + 1, k] * A_elem_ring[e] for e in 1:ne) / A_total_ring
        for e in 1:ne
            ir = div(e - 1, ntheta) + 1
            q_ring_nd[e, k] = q_mean_nd_ring[k] * (q_bin_nd[ir, k] / max(mean_ring_k, 1e-30))
        end
    end

    # 热惯量校准：用 Jellyroll 面积加权有效体积热容
    jparam = jellyroll_result.case.param
    fks_j = JuBat.jellyroll_element_properties(mesh, jparam)[2]
    rho_c_layers = [getproperty(jparam, L).rho * getproperty(jparam, L).heat_Q for L in [:NE, :SP, :PE, :PCC, :NCC]]
    rho_c_target = sum(fks_j * rho_c_layers .* A_elem) / A_total

    rho_c_ring = case_ring.param.cell.heat_Q / max(case_ring.param.cell.volume, 1e-30)
    inertia_scale = rho_c_target / max(rho_c_ring, 1e-30)
    case_ring.param.cell.heat_Q *= inertia_scale

    scale_ring = param_ring.scale
    T_ref = scale_ring.T_ref

    println("圆环网格创建完成: $ne 个单元")
    @printf("  热惯量校准系数: %.4f\n", inertia_scale)

    # 热源：已按半径分布存入 q_ring_nd(ne × nstep)，时间上插值
    q_ring_nd_min, q_ring_nd_max = minimum(q_ring_nd), maximum(q_ring_nd)
    q_ring_mean_t0 = area_weighted_mean(q_ring_nd[:, 1], A_elem_ring)
    rel_spread = (q_ring_nd_max - q_ring_nd_min) / max(q_ring_mean_t0, 1e-30) * 100
    println("  无量纲热源(径向分布)范围: [$q_ring_nd_min, $q_ring_nd_max], 相对差异: $(round(rel_spread; digits=1))%")

    # 初始化变量（首时刻径向分布）
    variables = Dict{String,Any}(
        "T_nodes" => fill(param_ring.cell.T0 / T_ref, mesh_ring.nlen),
        "thermal2D outer_nodes" => mesh_data.outer_nodes,
        "heat_source_fields" => copy(q_ring_nd[:, 1])
    )

    # 热源更新函数：按无量纲时间插值各单元无量纲热源
    t_vec_nd = collect(t) ./ param_ring.scale.t0
    function update_heat_source(t_nd, vars)
        q_at_t = [interp_linear(t_vec_nd, q_ring_nd[e, :], [t_nd])[1] for e in 1:ne]
        vars["heat_source_fields"] = q_at_t
    end

    println("开始圆环热模型时间推进...")

    ring_sol = JuBat.Solve(case_ring;
        thermal_variables=variables,
        thermal_update_fn=update_heat_source,
        thermal_record=true,
        polar_mesh_data=mesh_data)

    # ring_sol.time 现在是无量纲时间，需要转换回物理时间
    t_ring_phys = ring_sol.time .* param_ring.scale.t0
    A_elem_ring = compute_element_areas(mesh_ring)
    T_vol_ring = [area_weighted_mean(compute_element_mean(mesh_ring, ring_sol.T_hist[:, k] .* T_ref), A_elem_ring)
                  for k in 1:length(t_ring_phys)]

    println("圆环热模型求解完成")
    @printf("  时间步数: %d, 温升: %.2f K\n", length(t_ring_phys), T_vol_ring[end] - T_vol_ring[1])
    @printf("  温度: %.2f → %.2f K\n", T_vol_ring[1], T_vol_ring[end])

    return (; t=t_ring_phys, T_vol=T_vol_ring, T_hist=ring_sol.T_hist, mesh=mesh_ring,
            case=case_ring, param_dim=param_ring, T_ref)
end

# ============================================================================
# 第三部分：结果对比与可视化
# ============================================================================

function compare_and_visualize(jellyroll_result, ring_result)
    println("\n$SEP")
    println("第三部分：结果对比与可视化")
    println(SEP)

    out_dir = joinpath(@__DIR__, "..", "..", "output", "thermal_verify")
    isdir(out_dir) || mkpath(out_dir)

    t = jellyroll_result.t
    T_jelly = jellyroll_result.T_vol
    q_mean = jellyroll_result.q_mean_phys
    t_ring = ring_result.t  # 现在是物理时间 [s]
    T_ring = ring_result.T_vol

    T_ring_interp = interp_linear(collect(t_ring), T_ring, collect(t))
    dT = T_jelly .- T_ring_interp
    err_max, err_mean = maximum(abs.(dT)), mean(abs.(dT))
    err_rel = err_mean / (maximum(T_jelly) - minimum(T_jelly) + 1e-12)

    println("\n统计对比:")
    @printf("  Jellyroll 温升: %.4f K, 圆环模型温升: %.4f K\n", T_jelly[end] - T_jelly[1], T_ring[end] - T_ring[1])
    @printf("  温差: 最大 %.4f K, 平均 %.4f K, 相对 %.4f%%\n", err_max, err_mean, 100 * err_rel)

    # 图1：温度演化
    p1 = plot(t, T_jelly; label="Jellyroll SPMe-Thermal", linewidth=2,
              xlabel="Time [s]", ylabel="Temperature [K]", title="Temperature Evolution", legend=:bottomright)
    plot!(p1, t_ring, T_ring; label="Ring2D Thermal", linewidth=2, linestyle=:dash)
    savefig(p1, joinpath(out_dir, "jellyroll_vs_ring_temperature.png"))

    # 图2：热源
    p2 = plot(t, q_mean; label="Mean Heat Source", linewidth=2, color=:red,
              xlabel="Time [s]", ylabel="Heat Source [W/m³]", title="Input Heat Source")
    savefig(p2, joinpath(out_dir, "jellyroll_vs_ring_heat_source.png"))

    # 图3：温差
    p3 = plot(t, dT; label="ΔT (Jellyroll - Ring)", linewidth=2, color=:purple,
              xlabel="Time [s]", ylabel="ΔT [K]", title="Temperature Error", fillrange=0, fillalpha=0.3)
    hline!(p3, [0]; linestyle=:dot, color=:black, label=false)
    savefig(p3, joinpath(out_dir, "jellyroll_vs_ring_temperature_error.png"))

    # 图4：综合对比
    p4 = plot(layout=(2, 2), size=(1000, 800))
    plot!(p4[1], t, T_jelly; label="Jellyroll", linewidth=2)
    plot!(p4[1], t_ring, T_ring; label="Ring", linewidth=2, linestyle=:dash)
    plot!(p4[1]; xlabel="Time [s]", ylabel="Temperature [K]", title="Temperature Evolution")
    plot!(p4[2], t, q_mean; linewidth=2, color=:red, label=false,
          xlabel="Time [s]", ylabel="Heat Source [W/m³]", title="Mean Heat Source")
    plot!(p4[3], t, dT; linewidth=2, color=:purple, label=false,
          xlabel="Time [s]", ylabel="ΔT [K]", title="Temperature Error", fillrange=0, fillalpha=0.3)
    bar!(p4[4], 1:2, [T_jelly[end] - T_jelly[1], T_ring[end] - T_ring[1]];
         color=[:blue, :orange], xlabel="Model", ylabel="ΔT [K]", title="Temperature Rise",
         xticks=(1:2, ["Jellyroll", "Ring2D"]))

    savefig(p4, joinpath(out_dir, "jellyroll_vs_ring_summary.png"))

    # 保存CSV
    open(joinpath(out_dir, "jellyroll_vs_ring_data.csv"), "w") do io
        println(io, "time_s,T_jelly_K,T_ring_K,T_diff_K,q_mean_Wm3")
        for k in eachindex(t)
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g\n", t[k], T_jelly[k], T_ring_interp[k], dT[k], q_mean[k])
        end
    end

    println("输出文件已保存至: $out_dir")
    return (; err_max, err_mean, err_rel)
end

# ============================================================================
# 主函数
# ============================================================================

function main()
    println("\n$SEP")
    println("Jellyroll SPMe-热耦合 vs 圆环热模型对比验证")
    println(SEP)

    jellyroll_result = run_jellyroll_simulation()
    ring_result = run_ring_simulation(jellyroll_result)
    stats = compare_and_visualize(jellyroll_result, ring_result)

    println("\n$SEP")
    println("验证完成")
    println(SEP)
    @printf("""
    结果摘要:
      - Jellyroll 温升: %.4f K
      - 圆环模型温升:   %.4f K
      - 最大温差:       %.4f K
      - 平均温差:       %.4f K
      - 相对误差:       %.4f%%
    """,
        jellyroll_result.T_vol[end] - jellyroll_result.T_vol[1],
        ring_result.T_vol[end] - ring_result.T_vol[1],
        stats.err_max, stats.err_mean, 100 * stats.err_rel)

    return jellyroll_result, ring_result, stats
end

main()
