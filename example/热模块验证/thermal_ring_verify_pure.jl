#=
圆环纯热模型验证脚本
====================
验证二维圆环热模型（均匀内热源 + 外边界对流）与解析解的对比。

解析解：
    T(r) = T_f + (q/4k)(r_o² - r²) + (q/2k)r_i²ln(r/r_o) + (q/2h)(r_o - r_i²/r_o)

输出：
    - 误差分析文本报告
    - 径向误差曲线对比图
    - 收敛性分析图
=#

using LinearAlgebra, Statistics, Plots, Printf

if !isdefined(Main, :JuBat)
    include(joinpath(@__DIR__, "..", "..", "src", "JuBat.jl"))
end

# ============================================================================
# 辅助函数
# ============================================================================

"""
圆环热传导解析解
"""
function analytical_solution_ring(r, q, k, h, r_i, r_o, T_f)
    term1 = (q / (4.0 * k)) * (r_o^2 - r^2)
    term2 = (q / (2.0 * k)) * r_i^2 * log(r / r_o)
    term3 = (q / (2.0 * h)) * (r_o - (r_i^2 / r_o))
    return T_f + term1 + term2 + term3
end

"""
计算单元平均热源
"""
function compute_q_elem(mesh, q_func, t)
    ne = size(mesh.element, 1)
    q_elem = zeros(Float64, ne)
    w_elem = zeros(Float64, ne)

    wJ = mesh.gs.weight .* mesh.gs.detJ
    xg = mesh.gs.x[:, 1]
    yg = mesh.gs.x[:, 2]
    rg = hypot.(xg, yg)
    thetag = atan.(yg, xg)
    ele = mesh.gs.ele

    @inbounds for g in 1:length(wJ)
        e = ele[g]
        qg = q_func(rg[g], thetag[g], t)
        q_elem[e] += qg * wJ[g]
        w_elem[e] += wJ[g]
    end

    @inbounds for e in 1:ne
        w_elem[e] > 0.0 && (q_elem[e] /= w_elem[e])
    end

    return q_elem
end

"""
提取径向温度分布
"""
function radial_profile(mesh, T)
    r = hypot.(mesh.node[:, 1], mesh.node[:, 2])
    r_unique = unique(sort(r))
    T_avg = zeros(Float64, length(r_unique))
    for (i, rv) in enumerate(r_unique)
        idx = findall(abs.(r .- rv) .< 1e-10)
        T_avg[i] = mean(T[idx])
    end
    return r_unique, T_avg
end

"""
计算角度变化范围（用于验证轴对称性）
"""
function angular_variation(mesh, T)
    r = hypot.(mesh.node[:, 1], mesh.node[:, 2])
    r_unique = unique(sort(r))
    max_range = 0.0
    max_std = 0.0
    for rv in r_unique
        idx = findall(abs.(r .- rv) .< 1e-10)
        Tv = T[idx]
        max_range = max(max_range, maximum(Tv) - minimum(Tv))
        max_std = max(max_std, std(Tv))
    end
    return max_range, max_std
end

# ============================================================================
# 验证运行函数
# ============================================================================

"""
运行单个验证案例
"""
function run_verify_case(param_dim, ntheta, dr, model)
    opt = JuBat.Option()
    opt.model = "thermal"
    opt.thermal_enabled = true
    opt.thermalmodel = model
    opt.time = [0.0, 3600]
    opt.dt = [1.0, 10]

    case = JuBat.SetCase(param_dim, opt)

    scale = param_dim.scale
    T_ref = scale.T_ref
    q_ref = scale.q

    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    nr = round(Int, (Rout - Rin) / dr)
    mesh_data = JuBat.ring_mesh(case.param, ntheta=ntheta, nr=nr, gsorder=2)
    mesh = mesh_data.mesh
    case.mesh["thermal2D"] = mesh

    q0 = 2.0e5  # W/m³ 均匀热源
    q0_nd = q0 / q_ref
    q_func = (r, theta, t) -> q0

    variables = Dict{String,Any}()
    variables["T_nodes"] = fill(param_dim.cell.T0 / T_ref, mesh.nlen)
    variables["thermal2D outer_nodes"] = mesh_data.outer_nodes

    if model == "ring2D_polar"
        ne = size(mesh.element, 1)
        variables["heat_source_fields"] = fill(q0_nd, ne)
    end

    update_fn = (t, vars) -> begin
        if model == "ring2D_polar"
            vars["heat_source_fields"] = fill(q0_nd, ne)
        else
            vars["heat_source_fields"] = compute_q_elem(mesh, q_func, t) ./ q_ref
        end
    end

    # 稳态验证：直接求解 K*T + F = 0，避免瞬态截断误差
    update_fn(0.0, variables)
    if model == "ring2D_polar"
        MT, KT, FT = JuBat.ThermalPolar2D_Ring(case, variables, mesh_data)
    else
        MT, KT, FT = JuBat.ThermalDistributed2D_Ring(case, variables)
        KT, FT = JuBat.ThermalRing2D_BC(KT, FT, case, variables["thermal2D outer_nodes"], 0.0)
    end
    T = (-KT \ FT) .* T_ref

    # 计算误差
    r_nodes_dim = hypot.(mesh.node[:, 1], mesh.node[:, 2]) .* scale.L
    k_r = param_dim.cell.lambda_r
    h = param_dim.cell.h
    T_f = param_dim.cell.T_amb
    T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)

    err = T .- T_exact
    err_l2 = sqrt(mean(err .^ 2))
    err_linf = maximum(abs.(err))
    err_rel = err_l2 / max(1e-12, maximum(abs.(T_exact)))

    # 径向分布
    r_prof, T_r = radial_profile(mesh, T)
    r_prof_dim = r_prof .* scale.L
    T_r_exact = analytical_solution_ring.(r_prof_dim, q0, k_r, h, Rin, Rout, T_f)

    # 角度对称性
    ang_range, ang_std = angular_variation(mesh, T)

    return (
        mesh = mesh,
        T = T,
        r_prof = r_prof_dim,
        T_r = T_r,
        T_r_exact = T_r_exact,
        err = err,
        err_l2 = err_l2,
        err_linf = err_linf,
        err_rel = err_rel,
        ang_range = ang_range,
        ang_std = ang_std,
        Tmin = minimum(T),
        Tmax = maximum(T)
    )
end

# ============================================================================
# 输出函数
# ============================================================================

"""
打印误差分析报告
"""
function print_error_report(results::NamedTuple)
    println("\n" * "="^70)
    println("圆环热模型验证报告 - 误差分析")
    println("="^70)

    for (name, data) in pairs(results)
        println("\n[$name]")
        @printf("  温度范围 [K]       : [%.4f, %.4f]\n", data.Tmin, data.Tmax)
        @printf("  L2 误差 [K]        : %.6e\n", data.err_l2)
        @printf("  L∞ 误差 [K]        : %.6e\n", data.err_linf)
        @printf("  相对误差 [-]       : %.6e (%.4f%%)\n", data.err_rel, data.err_rel * 100)
        @printf("  角度不对称范围 [K] : %.6e\n", data.ang_range)
        @printf("  角度不对称标准差[K]: %.6e\n", data.ang_std)
    end

    println("\n" * "-"^70)
    println("跨方法一致性检查")
    println("-"^70)

    names = collect(keys(results))
    for i in 1:length(names)
        for j in (i+1):length(names)
            d1, d2 = results[names[i]], results[names[j]]
            l2_diff = sqrt(mean((d1.T .- d2.T).^2))
            linf_diff = maximum(abs.(d1.T .- d2.T))
            @printf("  %s vs %s: L2=%.6e K, L∞=%.6e K\n",
                    names[i], names[j], l2_diff, linf_diff)
        end
    end

    # 验收标准
    println("\n" * "-"^70)
    println("验收标准检查")
    println("-"^70)
    all_pass = true
    for (name, data) in pairs(results)
        pass_l2 = data.err_l2 < 1.0  # L2误差 < 1K
        pass_linf = data.err_linf < 2.0  # L∞误差 < 2K
        pass_sym = data.ang_range < 0.1  # 角度不对称 < 0.1K
        status = pass_l2 && pass_linf && pass_sym ? "PASS" : "FAIL"
        all_pass = all_pass && (status == "PASS")
        @printf("  %-15s: %s (L2:%.2e L∞:%.2e Sym:%.2e)\n",
                name, status, data.err_l2, data.err_linf, data.ang_range)
    end
    println("\n总体状态: ", all_pass ? "✓ PASS" : "✗ FAIL")
    println("="^70)
end

"""
绘制误差曲线
"""
function plot_error_curves(results, out_dir)
    isdir(out_dir) || mkpath(out_dir)

    # 1. 径向温度分布对比
    p_temp = Plots.plot(size=(800, 500))
    colors = [:blue, :red, :green]
    linestyles = [:solid, :solid, :solid]
    for (i, (name, data)) in enumerate(pairs(results))
        Plots.plot!(p_temp, data.r_prof, data.T_r, lw=2, color=colors[i],
                    label=string(name), linestyle=linestyles[i])
    end

    # 添加解析解
    first_data = first(values(results))
    Plots.plot!(p_temp, first_data.r_prof, first_data.T_r_exact, lw=2,
                color=:black, linestyle=:dash, label="Exact Solution")

    Plots.plot!(p_temp, xlabel="r [m]", ylabel="T [K]",
                title="Radial Temperature Profile Comparison")
    savefig(p_temp, joinpath(out_dir, "ring_temperature_profile.png"))

    # 2. 径向误差曲线
    p_err = Plots.plot(size=(800, 500))
    for (i, (name, data)) in enumerate(pairs(results))
        err_r = abs.(data.T_r .- data.T_r_exact)
        Plots.plot!(p_err, data.r_prof, err_r, lw=2, color=colors[i],
                    marker=:o, markersize=3, label=string(name))
    end
    Plots.plot!(p_err, xlabel="r [m]", ylabel="|T - T_exact| [K]",
                title="Radial Error Distribution", yscale=:log10)
    savefig(p_err, joinpath(out_dir, "ring_radial_error.png"))

    # 3. 全场误差分布（仅第一个方法）
    first_name = first(keys(results))
    first_data = first(values(results))

    p_field = Plots.plot(size=(600, 500), aspect_ratio=1)
    x = first_data.mesh.node[:, 1]
    y = first_data.mesh.node[:, 2]
    Plots.scatter!(p_field, x, y, marker_z=first_data.err,
                   markersize=2, color=:coolwarm, label="")
    Plots.plot!(p_field, xlabel="x [m]", ylabel="y [m]",
                title="Error Field ($first_name)", colorbar_title="Error [K]")
    savefig(p_field, joinpath(out_dir, "ring_error_field.png"))

    println("\n图表已保存至: $out_dir")
    println("  - ring_temperature_profile.png (温度分布对比)")
    println("  - ring_radial_error.png (径向误差曲线)")
    println("  - ring_error_field.png (全场误差分布)")
end

"""
绘制收敛性分析
"""
function plot_convergence(param_dim, dr, out_dir)
    isdir(out_dir) || mkpath(out_dir)

    ntheta_list = [20, 40, 80, 120, 160]

    err_l2_fem = Float64[]
    err_l2_polar = Float64[]
    err_linf_fem = Float64[]
    err_linf_polar = Float64[]

    println("\n收敛性分析中...")
    for nθ in ntheta_list
        print("  nθ = $nθ ... ")
        fem = run_verify_case(param_dim, nθ, dr, "ring2D")
        polar = run_verify_case(param_dim, nθ, dr, "ring2D_polar")
        push!(err_l2_fem, fem.err_l2)
        push!(err_l2_polar, polar.err_l2)
        push!(err_linf_fem, fem.err_linf)
        push!(err_linf_polar, polar.err_linf)
        println("L2: FEM=$(fem.err_l2) K, Polar=$(polar.err_l2) K")
    end

    # 绘制收敛曲线
    p_conv = Plots.plot(size=(800, 500))
    Plots.plot!(p_conv, ntheta_list, err_l2_fem, lw=2, marker=:o, label="FEM L2")
    Plots.plot!(p_conv, ntheta_list, err_l2_polar, lw=2, marker=:s, label="Polar L2")
    Plots.plot!(p_conv, ntheta_list, err_linf_fem, lw=2, marker=:o, linestyle=:dash, label="FEM L∞")
    Plots.plot!(p_conv, ntheta_list, err_linf_polar, lw=2, marker=:s, linestyle=:dash, label="Polar L∞")
    Plots.plot!(p_conv, xlabel="nθ (周向网格数)", ylabel="Error [K]",
                title="Convergence with Angular Resolution", yscale=:log10)
    savefig(p_conv, joinpath(out_dir, "ring_convergence.png"))

    println("\n收敛曲线已保存: $(joinpath(out_dir, "ring_convergence.png"))")

    return (ntheta=ntheta_list, err_l2_fem=err_l2_fem, err_l2_polar=err_l2_polar)
end

# ============================================================================
# 主程序
# ============================================================================

function main()
    println("\n" * "="^70)
    println("圆环纯热模型验证")
    println("="^70)
    println("模型: 二维圆环，均匀内热源 + 外边界对流冷却")
    println("解析解: 稳态径向温度分布")

    # 参数设置
    param_dim = JuBat.ChooseCell("Ring")
    scale = param_dim.scale
    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout

    @printf("\n几何参数:\n")
    @printf("  内半径 Rin = %.4f m\n", Rin)
    @printf("  外半径 Rout = %.4f m\n", Rout)
    @printf("  厚度 Δr = %.4f m\n", Rout - Rin)

    @printf("\n热物性参数:\n")
    @printf("  径向导热率 k_r = %.4f W/(m·K)\n", param_dim.cell.lambda_r)
    @printf("  对流换热系数 h = %.4f W/(m²·K)\n", param_dim.cell.h)
    @printf("  环境温度 T_f = %.2f K\n", param_dim.cell.T_amb)
    @printf("  热源强度 q = 2.0×10⁵ W/m³\n")

    # 网格参数
    ntheta = 40
    dr = (Rout - Rin) / 20

    @printf("\n网格参数:\n")
    @printf("  周向单元数 nθ = %d\n", ntheta)
    @printf("  径向单元数 nr = %d\n", round(Int, (Rout - Rin) / dr))

    # 运行验证
    println("\n运行验证...")
    fem = run_verify_case(param_dim, ntheta, dr, "ring2D")
    polar = run_verify_case(param_dim, ntheta, dr, "ring2D_polar")

    results = (FEM_Q4=fem, Polar_FVM=polar)

    # 输出目录
    out_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "thermal_ring_verify"))
    isdir(out_dir) || mkpath(out_dir)

    # 打印误差报告
    print_error_report(results)

    # 绘制误差曲线
    plot_error_curves(results, out_dir)

    # 收敛性分析（按当前调试需求先跳过）
    # conv_data = plot_convergence(param_dim, dr, out_dir)

    println("\n" * "="^70)
    println("验证完成")
    println("="^70)
end

main()
