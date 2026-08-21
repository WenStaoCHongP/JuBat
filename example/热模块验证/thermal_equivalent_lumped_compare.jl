"""
2D 等效集总量 vs Thermal.jl(lumped) 对比脚本

输出：
- output/thermal_equivalent_lumped_compare.csv
- output/thermal_equivalent_temperature_compare.png
- output/thermal_equivalent_power_compare.png
"""

using LinearAlgebra, Statistics, Printf, CSV, Plots

include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat

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
            a = (xk - x1) / (x2 - x1)
            yi[k] = y1 + a * (y2 - y1)
        end
    end
    return yi
end

function derivative_centered(t::Vector{Float64}, y::Vector{Float64})
    n = length(t)
    dy = zeros(Float64, n)
    if n == 1
        return dy
    end
    dy[1] = (y[2] - y[1]) / max(1e-12, (t[2] - t[1]))
    @inbounds for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / max(1e-12, (t[i + 1] - t[i - 1]))
    end
    dy[end] = (y[end] - y[end - 1]) / max(1e-12, (t[end] - t[end - 1]))
    return dy
end

function get_outer_edges(mesh, param)
    # 使用 JuBat 内置的边界节点识别函数
    is_inner, is_outer = JuBat.identify_boundary_nodes(mesh, param)

    seen = Set{Tuple{Int, Int}}()
    edges = Tuple{Int, Int}[]

    @inbounds for e in 1:size(mesh.element, 1)
        n1, n2, n3, n4 = mesh.element[e, :]
        for (a, b) in ((n1, n2), (n2, n3), (n3, n4), (n4, n1))
            if is_outer[a] && is_outer[b]
                key = a < b ? (a, b) : (b, a)
                if !(key in seen)
                    push!(seen, key)
                    push!(edges, key)
                end
            end
        end
    end
    return edges, is_outer
end

function build_option(; thermalmodel::String="distributed2D")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"
    opt.Current = _ -> 5.0
    opt.time = [0.0, 3600.0]
    opt.dt = [0.5, 10.0]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    opt.thermal_enabled = true
    opt.thermalmodel = thermalmodel  # 使用传入的参数
    opt.cool_method = "surface"
    opt.per_element_spme = true
    opt.debug_coupling = false
    opt.czm_enabled = false
    return opt
end

function run_distributed2d_result()
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = build_option(thermalmodel="distributed2D")
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh = case.mesh["thermal2D"]

    result = JuBat.Solve(case)

    t = result["time [s]"]
    T_nodes_hist = result["thermal2D T_nodes history [K]"]
    q_nd_hist = result["heat_source_fields"]

    A_elem = zeros(Float64, size(mesh.element, 1))
    @inbounds for g in 1:length(mesh.gs.detJ)
        e = mesh.gs.ele[g]
        A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end

    H = case.param_dim.cell.width
    q_ref = case.param_dim.scale.q  # 统一能量尺度热源参考 (P_ref / L^3)
    q_phys_hist = q_nd_hist .* q_ref

    # 几何尺度参数（网格坐标已归一化，需要转换为物理尺度）
    scale_L = case.param_dim.scale.L  # 长度尺度 [m]

    T_vol = zeros(Float64, length(t))
    T_edge_avg = zeros(Float64, length(t))  # 边界平均温度诊断
    P_internal = zeros(Float64, length(t))
    P_boundary_outer = zeros(Float64, length(t))
    P_boundary_surface = zeros(Float64, length(t))
    P_boundary_total = zeros(Float64, length(t))
    P_net = zeros(Float64, length(t))
    P_rxn = zeros(Float64, length(t))
    P_ohmic = zeros(Float64, length(t))
    P_reversible = zeros(Float64, length(t))

    h = case.param_dim.cell.h
    Tamb = case.param_dim.cell.T_amb
    Rout = case.param_dim.cell.Rout
    outer_edges, is_outer = get_outer_edges(mesh, case.param)

    # 调试输出：检查边界边识别
    r_nodes = sqrt.(mesh.node[:, 1].^2 .+ mesh.node[:, 2].^2)
    n_outer_nodes = count(is_outer)
    println("\n[边界边识别调试]")
    println("  Rout (有量纲): $Rout m")
    println("  网格坐标范围: x ∈ [$(minimum(mesh.node[:,1])), $(maximum(mesh.node[:,1]))], y ∈ [$(minimum(mesh.node[:,2])), $(maximum(mesh.node[:,2]))]")
    println("  无量纲半径范围: r ∈ [$(minimum(r_nodes)), $(maximum(r_nodes))]")
    println("  识别到的外边界节点数量: $n_outer_nodes")
    println("  识别到的外边界边数量: $(length(outer_edges))")

    # 详细检查外边界节点
    outer_node_indices = findall(is_outer)
    println("\n[外边界节点详情（前10个）]")
    println("  节点ID |     x     |     y     |    r    |")
    for i in 1:min(10, length(outer_node_indices))
        idx = outer_node_indices[i]
        x, y = mesh.node[idx, 1], mesh.node[idx, 2]
        r = hypot(x, y)
        @printf("  %5d | %9.4f | %9.4f | %7.4f |\n", idx, x, y, r)
    end

    q_rxn_ne = result["thermal2D Q_rxn_NE [W/m3]"]
    q_rxn_pe = result["thermal2D Q_rxn_PE [W/m3]"]
    q_rev_ne = result["thermal2D Q_rev_NE [W/m3]"]
    q_rev_pe = result["thermal2D Q_rev_PE [W/m3]"]
    q_ohm_s_ne = result["thermal2D Q_ohm_s_NE [W/m3]"]
    q_ohm_e_ne = result["thermal2D Q_ohm_e_NE [W/m3]"]
    q_sp = result["thermal2D Q_SP [W/m3]"]
    q_ohm_s_pe = result["thermal2D Q_ohm_s_PE [W/m3]"]
    q_ohm_e_pe = result["thermal2D Q_ohm_e_PE [W/m3]"]
    q_pcc = result["thermal2D Q_PCC [W/m3]"]
    q_ncc = result["thermal2D Q_NCC [W/m3]"]

    @inbounds for k in eachindex(t)
        T_nodes = T_nodes_hist[:, k]
        T_elem = JuBat.element_nodal_mean(mesh, T_nodes)
        T_vol[k] = sum(T_elem .* A_elem) / max(1e-12, sum(A_elem))
        # 注意：A_elem 是无量纲面积 (dA* = dA/L²)，需要乘以 L² 转换为物理面积
        P_internal[k] = sum(q_phys_hist[:, k] .* A_elem) * H * scale_L^2
        P_rxn[k] = sum((q_rxn_ne[:, k] .+ q_rxn_pe[:, k]) .* A_elem) * H * scale_L^2
        P_reversible[k] = sum((q_rev_ne[:, k] .+ q_rev_pe[:, k]) .* A_elem) * H * scale_L^2
        P_ohmic[k] = sum((q_ohm_s_ne[:, k] .+ q_ohm_e_ne[:, k] .+ q_sp[:, k] .+ q_ohm_s_pe[:, k] .+ q_ohm_e_pe[:, k] .+ q_pcc[:, k] .+ q_ncc[:, k]) .* A_elem) * H * scale_L^2

        # 外圆周对流散热：∫ h(T-Tamb)dA_outer
        # mesh.node 是无量纲坐标，边长需乘 scale_L 转换为物理长度
        p_out = 0.0
        total_edge_length = 0.0
        weighted_T_edge = 0.0  # 边长加权边界温度
        for (a, b) in outer_edges
            xa, ya = mesh.node[a, 1], mesh.node[a, 2]
            xb, yb = mesh.node[b, 1], mesh.node[b, 2]
            L_edge = hypot(xb - xa, yb - ya)
            T_edge = 0.5 * (T_nodes[a] + T_nodes[b])
            p_out += h * L_edge * scale_L * H * (T_edge - Tamb)
            weighted_T_edge += T_edge * L_edge
            total_edge_length += L_edge
        end
        P_boundary_outer[k] = p_out
        # 计算边长加权的边界平均温度
        T_edge_avg[k] = total_edge_length > 0 ? weighted_T_edge / total_edge_length : Tamb

        # surface 模式分布式散热：2h/H * ∫(T-Tamb)dV = 2h * ∫(T-Tamb)dA
        # A_elem 是无量纲面积，需乘 L² 转换为物理面积
        P_boundary_surface[k] = 2.0 * h * sum((T_elem .- Tamb) .* A_elem) * scale_L^2

        P_boundary_total[k] = P_boundary_outer[k] + P_boundary_surface[k]
        P_net[k] = P_internal[k] - P_boundary_total[k]
    end

    # 诊断输出：最终时刻的温度场分布
    T_nodes_final = T_nodes_hist[:, end]
    println("\n[最终时刻温度场诊断]")
    println("  总节点数: $(length(T_nodes_final))")
    println("  温度范围: T ∈ [$(minimum(T_nodes_final)), $(maximum(T_nodes_final))] K")
    println("  体积平均温度: $(T_vol[end]) K")

    # 外边界节点温度统计
    T_outer_nodes = T_nodes_final[outer_node_indices]
    println("\n[外边界节点温度统计]")
    println("  外边界节点温度范围: T ∈ [$(minimum(T_outer_nodes)), $(maximum(T_outer_nodes))] K")
    println("  外边界节点平均温度: $(mean(T_outer_nodes)) K")
    println("  外边界节点温度标准差: $(std(T_outer_nodes)) K")

    # 内部节点温度统计
    inner_node_indices = findall(.!is_outer)
    T_inner_nodes = T_nodes_final[inner_node_indices]
    println("\n[内部节点温度统计]")
    println("  内部节点温度范围: T ∈ [$(minimum(T_inner_nodes)), $(maximum(T_inner_nodes))] K")
    println("  内部节点平均温度: $(mean(T_inner_nodes)) K")

    # 温度梯度分析
    println("\n[温度梯度分析]")
    println("  外边界/内部 温度差: $(mean(T_inner_nodes) - mean(T_outer_nodes)) K")

    return t, T_vol, T_edge_avg, P_internal, P_boundary_total, P_net, P_rxn, P_ohmic, P_reversible
end

function run_lumped_result()
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = build_option(thermalmodel="lumped")
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    result = JuBat.Solve(case)

    t = result["time [s]"]
    T = result["temperature [K]"]
    P_internal = result["thermal lumped internal heat [W]"]

    h = case.param_dim.cell.h
    A_cool = case.param_dim.cell.cooling_surface
    Tamb = case.param_dim.cell.T_amb
    P_boundary = h .* A_cool .* (T .- Tamb)
    P_net = P_internal .- P_boundary

    return t, T, P_internal, P_boundary, P_net
end

function main()
    println("="^90)
    println("2D 等效集总量 vs Thermal.jl(lumped) 对比")
    println("="^90)

    t2d, T2d, T_edge, Pin2d, Pout2d, Pnet2d, Prxn2d, Pohm2d, Prev2d = run_distributed2d_result()
    tlp, Tlp, Pinlp, Poutlp, Pnetlp = run_lumped_result()

    # 将 lumped 插值到 2D 时间轴
    Tlp_i = interp_linear(tlp, Tlp, t2d)
    Pinlp_i = interp_linear(tlp, Pinlp, t2d)
    Poutlp_i = interp_linear(tlp, Poutlp, t2d)
    Pnetlp_i = interp_linear(tlp, Pnetlp, t2d)

    # PyBaMM 温度与等效净热源
    ref_tbl = CSV.File(joinpath(@__DIR__, "../../src/data/pybamm_SPMe_LGM50_1.0C.csv"))
    tref = Float64.(getproperty(ref_tbl, Symbol("time [s]")))
    Tref = Float64.(getproperty(ref_tbl, Symbol("temperature [K]")))
    Tref_i = interp_linear(tref, Tref, t2d)

    function get_ref_col(tbl, name::String, n::Int)
        try
            return Float64.(getproperty(tbl, Symbol(name)))
        catch
            return fill(NaN, n)
        end
    end

    Qtotal_ref = get_ref_col(ref_tbl, "Q_total [W]", length(tref))
    Qohm_ref = get_ref_col(ref_tbl, "Q_ohmic [W]", length(tref))
    Qirrev_ref = get_ref_col(ref_tbl, "Q_irreversible [W]", length(tref))
    Qrev_ref = get_ref_col(ref_tbl, "Q_reversible [W]", length(tref))

    Qtotal_ref_i = interp_linear(tref, Qtotal_ref, t2d)
    Qohm_ref_i = interp_linear(tref, Qohm_ref, t2d)
    Qirrev_ref_i = interp_linear(tref, Qirrev_ref, t2d)
    Qrev_ref_i = interp_linear(tref, Qrev_ref, t2d)

    p_dim = JuBat.ChooseCell("Jellyroll")
    m = p_dim.cell.mass
    cp = p_dim.cell.heat_Q
    dTdt_ref = derivative_centered(t2d, Tref_i)
    Pnet_ref = m .* cp .* dTdt_ref

    i2 = min(2, length(t2d))
    println("\n[平均功率(去首点)]")
    @printf("  2D等效   Pin/Pout/Pnet: %.6e / %.6e / %.6e W\n", mean(Pin2d[i2:end]), mean(Pout2d[i2:end]), mean(Pnet2d[i2:end]))
    @printf("  Lumped   Pin/Pout/Pnet: %.6e / %.6e / %.6e W\n", mean(Pinlp_i[i2:end]), mean(Poutlp_i[i2:end]), mean(Pnetlp_i[i2:end]))
    @printf("  PyBaMM等效净热源 Pnet: %.6e W\n", mean(Pnet_ref[i2:end]))
    @printf("  2D分项(反应/欧姆/可逆): %.6e / %.6e / %.6e W\n", mean(Prxn2d[i2:end]), mean(Pohm2d[i2:end]), mean(Prev2d[i2:end]))
    @printf("  PyBaMM分项(不可逆/欧姆/可逆): %.6e / %.6e / %.6e W\n", mean(Qirrev_ref_i[i2:end]), mean(Qohm_ref_i[i2:end]), mean(Qrev_ref_i[i2:end]))

    println("\n[温升]")
    @printf("  2D等效 ΔT: %.6f K\n", T2d[end] - T2d[1])
    @printf("  Lumped ΔT: %.6f K\n", Tlp_i[end] - Tlp_i[1])
    @printf("  PyBaMM ΔT: %.6f K\n", Tref_i[end] - Tref_i[1])

    # ====== 边界温度诊断输出 ======
    println("\n" * "="^90)
    println("边界温度 T_edge 诊断")
    println("="^90)

    p_dim = JuBat.ChooseCell("Jellyroll")
    Tamb = p_dim.cell.T_amb
    h = p_dim.cell.h
    A_cool = p_dim.cell.cooling_surface

    # 选择关键时间点进行诊断
    key_times = [0, div(length(t2d), 4), div(length(t2d), 2), 3*div(length(t2d), 4), length(t2d)]

    println("\n[关键时间点边界温度分析]")
    println("时间(s) | T_vol(K) | T_edge(K) | ΔT=T_vol-T_edge(K) | T_vol-Tamb(K) | T_edge-Tamb(K) | 散热比")
    println("-" * repeat("-", 95))

    for idx in key_times
        if idx < 1; idx = 1; end
        if idx > length(t2d); idx = length(t2d); end

        T_v = T2d[idx]
        T_e = T_edge[idx]
        dT = T_v - T_e
        dT_vol = T_v - Tamb
        dT_edge = T_e - Tamb

        # 散热比 = (T_edge - Tamb) / (T_vol - Tamb)
        dissipation_ratio = dT_vol > 1e-6 ? dT_edge / dT_vol : 0.0

        @printf("%7.1f | %8.2f | %9.2f | %18.2f | %13.2f | %14.2f | %6.2f%%\n",
                t2d[idx], T_v, T_e, dT, dT_vol, dT_edge, dissipation_ratio * 100)
    end

    # 计算等效集总散热和实际散热对比
    println("\n[等效集总散热 vs 实际散热分析]")
    println("时间(s) | Pout_2d(W) | P_equiv(W) | 散热效率 | 理论 h_eff/W/m²K")
    println("-" * repeat("-", 70))

    for idx in key_times
        if idx < 1; idx = 1; end
        if idx > length(t2d); idx = length(t2d); end

        T_v = T2d[idx]
        T_e = T_edge[idx]
        dT_vol = T_v - Tamb
        dT_edge = T_e - Tamb

        # 实际 2D 散热
        P_actual = Pout2d[idx]

        # 等效集总散热 = h * A_cool * (T_vol - Tamb)
        P_equiv = h * A_cool * dT_vol

        # 散热效率
        eff = P_equiv > 1e-12 ? P_actual / P_equiv : 0.0

        # 理论等效换热系数 h_eff = h * (T_vol - Tamb) / (T_edge - Tamb)
        h_eff = dT_edge > 1e-6 ? h * dT_vol / dT_edge : h

        @printf("%7.1f | %10.4f | %10.4f | %8.2f%% | %15.2f\n",
                t2d[idx], P_actual, P_equiv, eff * 100, h_eff)
    end

    println("\n[诊断结论]")
    @printf("  - 边界温度比平均温度低: %.2f - %.2f K (时间范围 0-%.0f s)\n",
            minimum(T2d .- T_edge), maximum(T2d .- T_edge), t2d[end])
    @printf("  - 平均散热效率: %.2f%% (2D实际散热 / 等效集总散热)\n",
            mean(Pout2d ./ (h .* A_cool .* (T2d .- Tamb) .+ 1e-12)) * 100)
    @printf("  - 建议: 使用等效换热系数 h_eff ≈ %.1f × h 可使 2D 散热与集总模型等效\n",
            mean((T2d .- Tamb) ./ (T_edge .- Tamb .+ 1e-6)))

    out_dir = joinpath(@__DIR__, "..", "..", "output", "thermal_equivalent_lumped_compare")
    isdir(out_dir) || mkpath(out_dir)
    out_csv = joinpath(out_dir, "thermal_equivalent_lumped_compare.csv")
    open(out_csv, "w") do io
        println(io, "time_s,T_2d_equiv_K,T_edge_K,T_lumped_K,T_pybamm_K,Pin_2d_W,Pout_2d_W,Pnet_2d_W,Pin_lumped_W,Pout_lumped_W,Pnet_lumped_W,Pnet_pybamm_equiv_W,P_rxn_2d_W,P_ohm_2d_W,P_rev_2d_W,Q_total_pybamm_W,Q_ohmic_pybamm_W,Q_irreversible_pybamm_W,Q_reversible_pybamm_W")
        @inbounds for k in eachindex(t2d)
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                t2d[k], T2d[k], T_edge[k], Tlp_i[k], Tref_i[k], Pin2d[k], Pout2d[k], Pnet2d[k], Pinlp_i[k], Poutlp_i[k], Pnetlp_i[k], Pnet_ref[k],
                Prxn2d[k], Pohm2d[k], Prev2d[k], Qtotal_ref_i[k], Qohm_ref_i[k], Qirrev_ref_i[k], Qrev_ref_i[k])
        end
    end

    pT = plot(t2d, T2d, label="2D equivalent lumped T", linewidth=2, xlabel="Time (s)", ylabel="Temperature (K)", title="Temperature Comparison")
    plot!(pT, t2d, T_edge, label="2D boundary T_edge", linewidth=2, linestyle=:dot)
    plot!(pT, t2d, Tlp_i, label="Thermal.jl lumped T", linewidth=2)
    plot!(pT, t2d, Tref_i, label="PyBaMM T", linewidth=2, linestyle=:dash)
    savefig(pT, joinpath(out_dir, "thermal_equivalent_temperature_compare.png"))

    pP = plot(t2d, Pnet2d, label="2D equivalent net", linewidth=2, xlabel="Time (s)", ylabel="Power (W)", title="Net Heat Source Comparison")
    plot!(pP, t2d, Pnetlp_i, label="Thermal.jl lumped net", linewidth=2)
    plot!(pP, t2d, Pnet_ref, label="PyBaMM equivalent net", linewidth=2, linestyle=:dash)
    savefig(pP, joinpath(out_dir, "thermal_equivalent_power_compare.png"))

    pC = plot(t2d, Prxn2d, label="JuBat reaction", linewidth=2, xlabel="Time (s)", ylabel="Power (W)", title="Heat Components Comparison")
    plot!(pC, t2d, Pohm2d, label="JuBat ohmic", linewidth=2)
    plot!(pC, t2d, Prev2d, label="JuBat reversible", linewidth=2)
    plot!(pC, t2d, Qirrev_ref_i, label="PyBaMM irreversible", linewidth=2, linestyle=:dash)
    plot!(pC, t2d, Qohm_ref_i, label="PyBaMM ohmic", linewidth=2, linestyle=:dash)
    plot!(pC, t2d, Qrev_ref_i, label="PyBaMM reversible", linewidth=2, linestyle=:dash)
    savefig(pC, joinpath(out_dir, "thermal_equivalent_component_compare.png"))

    println("\n已写出: output/thermal_equivalent_lumped_compare.csv")
    println("已写出: output/thermal_equivalent_temperature_compare.png")
    println("已写出: output/thermal_equivalent_power_compare.png")
    println("已写出: output/thermal_equivalent_component_compare.png")
end

main()
