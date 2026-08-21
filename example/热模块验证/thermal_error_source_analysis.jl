"""
热误差来源分析脚本（忽略模型口径差异）

目标：
1) 在与 testexample.jl 一致的 60 s 工况下运行 JuBat
2) 计算并对比：
   - 内热源功率 P_internal [W]
   - 边界散热功率 P_boundary [W]
   - 净热源功率 P_net = P_internal - P_boundary [W]
3) 基于 PyBaMM 温度曲线反推等效净热源：P_net_ref = m*cp*dT/dt [W]
4) 输出 CSV 与图像，辅助诊断温升偏小原因
"""

using LinearAlgebra, Statistics, Printf, CSV, Plots

include(joinpath(@__DIR__, "../src/JuBat.jl"))
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

function get_outer_edges(mesh, Rout::Float64)
    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    r = sqrt.(x .^ 2 .+ y .^ 2)
    tol = max(1e-8, 1e-4 * Rout)
    is_outer = abs.(r .- Rout) .<= tol

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
    return edges
end

function main()
    println("="^88)
    println("JuBat vs PyBaMM 热误差来源分析")
    println("="^88)

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

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

    crates = 1.0
    iapp = 5.0 * crates
    opt.Current = _ -> iapp

    opt.time = [0.0, 60.0]
    opt.dt = [0.5, 10.0]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true
    opt.debug_coupling = false
    opt.czm_enabled = false

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh = case.mesh["thermal2D"]

    result = JuBat.Solve(case)

    t = result["time [s]"]
    nstep = length(t)
    T_nodes_hist = result["thermal2D T_nodes history [K]"]
    q_nd_hist = result["heat_source_fields"]
    component_keys = [
        "thermal2D Q_rxn_NE [W/m3]",
        "thermal2D Q_rev_NE [W/m3]",
        "thermal2D Q_ohm_s_NE [W/m3]",
        "thermal2D Q_ohm_e_NE [W/m3]",
        "thermal2D Q_SP [W/m3]",
        "thermal2D Q_rxn_PE [W/m3]",
        "thermal2D Q_rev_PE [W/m3]",
        "thermal2D Q_ohm_s_PE [W/m3]",
        "thermal2D Q_ohm_e_PE [W/m3]",
        "thermal2D Q_PCC [W/m3]",
        "thermal2D Q_NCC [W/m3]",
    ]

    # 体积与面积
    A_elem = zeros(Float64, size(mesh.element, 1))
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    A_total = sum(A_elem)
    H = case.param_dim.cell.width
    V_total = case.param_dim.cell.volume

    # JuBat 物理热源标度（统一能量尺度）
    q_ref = case.param_dim.scale.q  # P_ref / L^3
    q_phys_hist = q_nd_hist .* q_ref

    # 几何尺度参数（网格坐标已归一化，需要转换为物理尺度）
    scale_L = case.param_dim.scale.L  # 长度尺度 [m]

    # 外边界对流项
    Rout = case.param_dim.cell.Rout
    outer_edges = get_outer_edges(mesh, Rout)

    h = case.param_dim.cell.h
    Tamb = case.param_dim.cell.T_amb

    P_internal = zeros(Float64, nstep)
    P_boundary_outer = zeros(Float64, nstep)
    P_boundary_surface = zeros(Float64, nstep)
    P_boundary_total = zeros(Float64, nstep)
    P_net_jubat = zeros(Float64, nstep)
    T_vol = zeros(Float64, nstep)
    P_components = Dict{String, Vector{Float64}}()
    for key in component_keys
        P_components[key] = zeros(Float64, nstep)
    end

    @inbounds for k in 1:nstep
        T_nodes = T_nodes_hist[:, k]
        T_elem = JuBat.element_nodal_mean(mesh, T_nodes)

        # 体积平均温度
        T_vol[k] = sum(T_elem .* A_elem) / max(1e-12, A_total)

        # 内热源总功率：∫ q dV
        # 注意：A_elem 是无量纲面积 (dA* = dA/L²)，需要乘以 L² 转换为物理面积
        q_phys = q_phys_hist[:, k]
        P_internal[k] = sum(q_phys .* A_elem) * H * scale_L^2

        for key in component_keys
            q_comp = result[key][:, k]
            P_components[key][k] = sum(q_comp .* A_elem) * H * scale_L^2
        end

        # 外圆周对流散热：∫ h(T-Tamb)dA
        # 注意：mesh.node 是无量纲坐标，边长 L 是无量纲的，需要乘以 scale_L 转换为物理长度
        p_out = 0.0
        for (a, b) in outer_edges
            xa, ya = mesh.node[a, 1], mesh.node[a, 2]
            xb, yb = mesh.node[b, 1], mesh.node[b, 2]
            L = hypot(xb - xa, yb - ya)
            T_edge = 0.5 * (T_nodes[a] + T_nodes[b])
            p_out += h * L * scale_L * H * (T_edge - Tamb)
        end
        P_boundary_outer[k] = p_out

        # surface 模式额外分布式散热（对应 ThermalDistributed.jl 里的 conv_factor）
        # 等效物理形式：2h/H * ∫(T-Tamb)dV = 2h * ∫(T-Tamb)dA
        # 注意：A_elem 是无量纲面积，需要乘以 L² 转换为物理面积
        P_boundary_surface[k] = 2.0 * h * sum((T_elem .- Tamb) .* A_elem) * scale_L^2

        P_boundary_total[k] = P_boundary_outer[k] + P_boundary_surface[k]
        P_net_jubat[k] = P_internal[k] - P_boundary_total[k]
    end

    # PyBaMM 温度曲线 -> 等效净热源（m*cp*dT/dt）
    ref_path = joinpath(@__DIR__, "../src/data/pybamm_SPMe_LGM50_1.0C.csv")
    ref_tbl = CSV.File(ref_path)
    t_ref = Float64.(getproperty(ref_tbl, Symbol("time [s]")))
    T_ref = Float64.(getproperty(ref_tbl, Symbol("temperature [K]")))

    T_ref_on_t = interp_linear(t_ref, T_ref, t)
    dTdt_ref = derivative_centered(t, T_ref_on_t)

    m_cell = case.param_dim.cell.mass
    cp_cell = case.param_dim.cell.heat_Q
    P_net_ref = m_cell .* cp_cell .* dTdt_ref

    # 统计
    i2 = min(2, nstep)
    mean_P_internal = mean(P_internal[i2:end])
    mean_P_bnd = mean(P_boundary_total[i2:end])
    mean_P_net_j = mean(P_net_jubat[i2:end])
    mean_P_net_ref = mean(P_net_ref[i2:end])

    ratio_net = mean_P_net_j / max(1e-12, mean_P_net_ref)
    gap_net = mean_P_net_j - mean_P_net_ref

    println("\n[统计: 60 s 区间内平均功率，去掉首点]")
    @printf("  JuBat 内热源 P_internal:      %.6e W\n", mean_P_internal)
    @printf("  JuBat 边界散热 P_boundary:     %.6e W\n", mean_P_bnd)
    @printf("  JuBat 净热源 P_net:           %.6e W\n", mean_P_net_j)
    @printf("  PyBaMM 等效净热源 P_net_ref:  %.6e W\n", mean_P_net_ref)
    @printf("  净热源比值 JuBat/PyBaMM:      %.4f\n", ratio_net)
    @printf("  净热源差值 JuBat-PyBaMM:      %.6e W\n", gap_net)

    println("\n[JuBat 内热源分项占比（平均功率）]")
    comp_mean = Dict{String, Float64}()
    for key in component_keys
        comp_mean[key] = mean(P_components[key][i2:end])
    end
    for key in component_keys
        frac = comp_mean[key] / max(1e-12, mean_P_internal)
        @printf("  %-28s : % .6e W  (%6.2f%%)\n", key, comp_mean[key], 100 * frac)
    end
    gain_required = mean_P_net_ref / max(1e-12, mean_P_net_j)
    @printf("\n  对齐 PyBaMM 净热源所需总增益: %.3f x\n", gain_required)

    @printf("\n[温升对比]\n")
    @printf("  JuBat:  %.6f K -> %.6f K, ΔT = %.6f K\n", T_vol[1], T_vol[end], T_vol[end] - T_vol[1])
    @printf("  PyBaMM: %.6f K -> %.6f K, ΔT = %.6f K\n", T_ref_on_t[1], T_ref_on_t[end], T_ref_on_t[end] - T_ref_on_t[1])

    # 写 CSV
    out_dir = joinpath(@__DIR__, "..", "..", "output", "thermal_error_source_analysis")
    isdir(out_dir) || mkpath(out_dir)
    out_csv = joinpath(out_dir, "thermal_error_breakdown.csv")
    open(out_csv, "w") do io
        println(io, "time_s,T_jubat_vol_K,T_pybamm_interp_K,P_internal_W,P_boundary_outer_W,P_boundary_surface_W,P_boundary_total_W,P_net_jubat_W,P_net_pybamm_equiv_W,P_Q_rxn_NE_W,P_Q_rev_NE_W,P_Q_ohm_s_NE_W,P_Q_ohm_e_NE_W,P_Q_SP_W,P_Q_rxn_PE_W,P_Q_rev_PE_W,P_Q_ohm_s_PE_W,P_Q_ohm_e_PE_W,P_Q_PCC_W,P_Q_NCC_W")
        @inbounds for k in 1:nstep
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                t[k], T_vol[k], T_ref_on_t[k], P_internal[k], P_boundary_outer[k], P_boundary_surface[k],
                P_boundary_total[k], P_net_jubat[k], P_net_ref[k],
                P_components["thermal2D Q_rxn_NE [W/m3]"][k],
                P_components["thermal2D Q_rev_NE [W/m3]"][k],
                P_components["thermal2D Q_ohm_s_NE [W/m3]"][k],
                P_components["thermal2D Q_ohm_e_NE [W/m3]"][k],
                P_components["thermal2D Q_SP [W/m3]"][k],
                P_components["thermal2D Q_rxn_PE [W/m3]"][k],
                P_components["thermal2D Q_rev_PE [W/m3]"][k],
                P_components["thermal2D Q_ohm_s_PE [W/m3]"][k],
                P_components["thermal2D Q_ohm_e_PE [W/m3]"][k],
                P_components["thermal2D Q_PCC [W/m3]"][k],
                P_components["thermal2D Q_NCC [W/m3]"][k])
        end
    end
    println("\n已写出: " * out_csv)

    # 图1：功率分解
    p1 = plot(t, P_internal, label="JuBat Internal", linewidth=2, xlabel="Time (s)", ylabel="Power (W)", title="Heat Power Breakdown")
    plot!(p1, t, P_boundary_total, label="JuBat Boundary Loss", linewidth=2)
    plot!(p1, t, P_net_jubat, label="JuBat Net", linewidth=2, color=:black)
    plot!(p1, t, P_net_ref, label="PyBaMM Net (equiv)", linewidth=2, linestyle=:dash, color=:red)
    savefig(p1, joinpath(out_dir, "thermal_error_power_breakdown.png"))

    # 图2：温度对比
    p2 = plot(t, T_vol, label="JuBat Volume-Avg T", linewidth=2, xlabel="Time (s)", ylabel="Temperature (K)", title="Temperature Comparison")
    plot!(p2, t, T_ref_on_t, label="PyBaMM Interp T", linewidth=2, linestyle=:dash)
    savefig(p2, joinpath(out_dir, "thermal_error_temperature_compare.png"))

    p3 = plot(t, P_components["thermal2D Q_rxn_NE [W/m3]"], label="Q_rxn_NE", linewidth=2,
              xlabel="Time (s)", ylabel="Power (W)", title="JuBat Internal Heat Components")
    plot!(p3, t, P_components["thermal2D Q_rxn_PE [W/m3]"], label="Q_rxn_PE", linewidth=2)
    plot!(p3, t, P_components["thermal2D Q_ohm_e_NE [W/m3]"], label="Q_ohm_e_NE", linewidth=2)
    plot!(p3, t, P_components["thermal2D Q_ohm_e_PE [W/m3]"], label="Q_ohm_e_PE", linewidth=2)
    plot!(p3, t, P_components["thermal2D Q_SP [W/m3]"], label="Q_SP", linewidth=2)
    savefig(p3, joinpath(out_dir, "thermal_error_component_breakdown.png"))

    println("已写出: output/thermal_error_power_breakdown.png")
    println("已写出: output/thermal_error_temperature_compare.png")
    println("已写出: output/thermal_error_component_breakdown.png")
end

main()
