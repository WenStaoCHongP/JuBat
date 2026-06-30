"""
测试案例：Jellyroll 电池单次放电二维应力场/位移场后处理

功能：
- 多 SPMe 并行 + 二维分布式热模型（与 testexample 同工况，但关闭 CZM）
- 在用户指定的时间节点调用 thermal_diffusion_stress_2D 计算宏观应力/位移场
- 输出一张大图：行 = 场分量（σ_xx/σ_yy/σ_xy/σ_vm/|U|），列 = 时间节点

⚠️ 尺度说明：本脚本输出的是 **极片/电极尺度（coating-scale, 二维平面应力）** 的场，
   由全叠合厚度加权有效模量 E_eff（~5e8 Pa 量级）计算；
   与颗粒尺度 Calstressdisp（颗粒 E ~1e10 Pa）不同，二者不可混用。
   注意：thermal_diffusion_stress_2D 返回 σ 为归一化值，本脚本乘 case.param_dim.scale.E_coat 还原为 Pa。

日期：2026-06-30
"""

using Printf
ENV["PLOTS_DEFAULT_BACKEND"] = get(ENV, "JUBAT_PLOTS_BACKEND", "gr")
using Plots
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 用户可配置参数 ============================================================
# 指定要绘图的物理时间节点（秒）；脚本会自动在 result 时间序列里找最接近的索引
const PLOT_TIMES_S = [600, 1800, 3600]
# ===========================================================================

# --------------------------------------------------------------------------
# Helper: 提取每个 Q4 单元的四角坐标（首尾闭合，5 个点）
# --------------------------------------------------------------------------
function q4_element_polygons(mesh)
    ne = size(mesh.element, 1)
    xs = Vector{Vector{Float64}}(undef, ne)
    ys = Vector{Vector{Float64}}(undef, ne)
    @inbounds for e in 1:ne
        n1, n2, n3, n4 = mesh.element[e, 1], mesh.element[e, 2], mesh.element[e, 3], mesh.element[e, 4]
        xs[e] = [mesh.node[n1,1], mesh.node[n2,1], mesh.node[n3,1], mesh.node[n4,1], mesh.node[n1,1]]
        ys[e] = [mesh.node[n1,2], mesh.node[n2,2], mesh.node[n3,2], mesh.node[n4,2], mesh.node[n1,2]]
    end
    return xs, ys
end

# --------------------------------------------------------------------------
# Helper: 在子图 plt 上绘制单元常数场 field_elem（长度 ne）的填色云图。
# 采用"每个单元单独 plot!(Shape; color=...)"——GR 后端最稳健的方案。
# deform_xy = (U_x, U_y)：若提供，按 def_scale 放大叠加变形网格。
# 注：xs[e]/ys[e] 含 5 个点（第 5 = 第 1 闭合），所以 deform 循环用 mod1(k,4)
#     把第 5 个点也映射到 node 1，保持闭合。
# --------------------------------------------------------------------------
function plot_q4_field!(plt, mesh, field_elem, xs, ys; title="", cmap=:RdBu,
                        clim=nothing, deform_xy=nothing, def_scale=0.0)
    fmin, fmax = (clim === nothing) ? (minimum(field_elem), maximum(field_elem)) : clim
    if fmax - fmin < 1e-30
        fmax = fmin + 1.0  # 避免除零
    end
    cmap_grad = cgrad(cmap)
    for e in 1:length(field_elem)
        xe = copy(xs[e]); ye = copy(ys[e])
        if deform_xy !== nothing
            U_x, U_y = deform_xy
            for k in 1:5
                n = mesh.element[e, mod1(k, 4)]
                xe[k] += def_scale * U_x[n]
                ye[k] += def_scale * U_y[n]
            end
        end
        idx = (field_elem[e] - fmin) / (fmax - fmin)
        col = cmap_grad[clamp(idx, 0.0, 1.0)]
        plot!(plt, Plots.Shape(xe, ye); linecolor=:match, linewidth=0.1,
              color=col, label=false, colorbar_entry=false)
    end
    title!(plt, title); xlabel!(plt, "x [m]"); ylabel!(plt, "y [m]")
end

function main()
    println("="^80)
    println("Jellyroll 单次放电 二维应力/位移场后处理")
    println("="^80)

    # ====================================================================
    # 1. 参数设置（同 testexample，唯一差异：czm_enabled = false）
    # ====================================================================
    println("\n[1/5] 参数设置...")

    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = JuBat.Option()
    Crates = 1.0
    i = 5 * Crates
    opt.Current = x -> i
    opt.model = "SPMe"
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"   # 不启用颗粒尺度应力耦合

    opt.time = [0.0, 3600.0]
    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"

    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true

    opt.czm_enabled = false        # ← 关键差异：关闭 CZM，仅传统固体力学

    println("OK: 参数设置完成")
    @printf("  电流: %.2f A (%.2f C)\n", i, Crates)
    @printf("  仿真时间: %.1f s\n", opt.time[end])
    println("  CZM: 关闭（仅传统固体力学）")

    # ====================================================================
    # 2. 创建案例与 Jellyroll 网格
    # ====================================================================
    println("\n[2/5] 创建案例与 Jellyroll 网格...")

    case = JuBat.SetCase(param_dim, opt)

    n_theta = 360
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=n_theta, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    mesh_th = case.mesh["thermal2D"]

    # czm_enabled=false：不创建 case.czm_mesh

    ne   = size(mesh_th.element, 1)
    nlen = mesh_th.nlen
    println("OK: 网格创建完成")
    @printf("  n_theta=%d, ne=%d, nlen=%d\n", n_theta, ne, nlen)

    # 前置字段核查（避免后续 helper 误用字段名）
    @assert hasproperty(mesh_th, :node)    "thermal2D mesh 缺少 .node 字段，实际字段: $(propertynames(mesh_th))"
    @assert hasproperty(mesh_th, :element) "thermal2D mesh 缺少 .element 字段"
    @assert size(mesh_th.node, 2) >= 2     "mesh_th.node 列数 < 2"
    @assert size(mesh_th.element, 2) == 4  "Q4 单元应有 4 列"

    # ====================================================================
    # 3. 求解
    # ====================================================================
    println("\n[3/5] 运行 SPMe+热 求解器（无 CZM）...")

    t_wall_start = time_ns()
    result = JuBat.Solve(case)
    t_wall_s = (time_ns() - t_wall_start) * 1e-9
    println("OK: 求解完成")
    @printf("  wall-clock: %.2f s\n", t_wall_s)

    # ====================================================================
    # 4. 时间节点选择
    # ====================================================================
    println("\n[4/5] 选择时间节点...")

    t_s = result["time [s]"]
    @assert haskey(result, "thermal2D temperature at nodes [K]") "result 缺少 'thermal2D temperature at nodes [K]'，请确认 opt.thermalmodel == \"distributed2D\""
    @assert haskey(result, "thermal2D element soc_n") "result 缺少 'thermal2D element soc_n'"
    @assert haskey(result, "thermal2D element soc_p") "result 缺少 'thermal2D element soc_p'"

    ti_list = Int[]
    for t_target in PLOT_TIMES_S
        if t_target < t_s[1] || t_target > t_s[end]
            @warn "PLOT_TIMES_S 中 $(t_target) s 超出 [$(t_s[1]), $(t_s[end])] 范围，钳到端点" maxlog=1
        end
        ti = argmin(abs.(t_s .- t_target))
        push!(ti_list, ti)
    end
    unique!(ti_list)
    sort!(ti_list)

    println("OK: 选定时间节点")
    for ti in ti_list
        @printf("  ti=%d  t=%.1f s\n", ti, t_s[ti])
    end

    # ====================================================================
    # 5. 在所有选定时间节点重建场 + 数量级 sanity check
    # ====================================================================
    println("\n[5/5] 重建二维场（$(length(ti_list)) 个时间节点）...")

    L_ref = case.param.scale.L

    # NamedTuple 向量（保留类型信息）
    FieldPack = NamedTuple{(:t, :σ_xx, :σ_yy, :σ_xy, :σ_vm, :U_x, :U_y),
                           Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64},
                                 Vector{Float64}, Vector{Float64}, Vector{Float64}}}
    fields_per_ti = FieldPack[]

    σ_vm_global_max = 0.0
    Umag_global_max = 0.0
    for ti in ti_list
        variables_ti = Dict{String, Union{Array{Float64},Float64}}(
            "T_nodes"                 => result["thermal2D temperature at nodes [K]"][:, ti],
            "thermal2D element soc_n" => result["thermal2D element soc_n"][:, ti],
            "thermal2D element soc_p" => result["thermal2D element soc_p"][:, ti],
        )
        @assert length(variables_ti["T_nodes"]) == nlen
        @assert length(variables_ti["thermal2D element soc_n"]) == ne

        vars_out = JuBat.thermal_diffusion_stress_2D(case, variables_ti)
        # thermal_diffusion_stress_2D 返回的 σ 为归一化值（E_eff 通过 scale.E_coat 归一化，
        # 见 src/CouplingState.jl:220-223 compute_effective_coating_modulus docstring）；
        # 乘 scale.E_coat (~5e8 Pa) 还原为物理 Pa。位移字段函数内已 × L_ref，无需再缩放。
        E_coat_scale = case.param_dim.scale.E_coat
        σ_xx = vars_out["diffusion stress xx"]      .* E_coat_scale
        σ_yy = vars_out["diffusion stress yy"]      .* E_coat_scale
        σ_xy = vars_out["diffusion stress xy"]      .* E_coat_scale
        σ_vm = vars_out["diffusion stress vonMises"] .* E_coat_scale
        U_x  = vars_out["displacement x"]
        U_y  = vars_out["displacement y"]
        @assert length(σ_vm) == ne
        @assert length(U_x) == nlen

        Umag = sqrt.(U_x.^2 .+ U_y.^2)
        σ_vm_global_max = max(σ_vm_global_max, maximum(σ_vm))
        Umag_global_max = max(Umag_global_max, maximum(Umag))

        push!(fields_per_ti, (t=t_s[ti], σ_xx=σ_xx, σ_yy=σ_yy, σ_xy=σ_xy,
                              σ_vm=σ_vm, U_x=U_x, U_y=U_y))

        @printf("  t=%7.1f s: σ_vm_max=%.3e Pa (%.2f MPa), |U|_max=%.3e m (%.3f µm)\n",
                t_s[ti], maximum(σ_vm), maximum(σ_vm)*1e-6,
                maximum(Umag), maximum(Umag)*1e6)
    end

    # === sanity check：量级 + 量纲 ===
    ti_last = ti_list[end]
    soc_n_last = result["thermal2D element soc_n"][:, ti_last]
    @printf("  sanity: param.NE.cs0(normalized)=%.4f  (期望 [0,1], 典型~0.9)\n", case.param.NE.cs0)
    @printf("  sanity: param_dim.NE.cs0=%.1f [mol/m3]  (物理浓度, 仅对比)\n", case.param_dim.NE.cs0)
    @printf("  sanity: soc_n range=[%.4f, %.4f]  (期望 [0,1])\n",
            minimum(soc_n_last), maximum(soc_n_last))
    @printf("  sanity: σ_vm_global_max=%.3e Pa（电极尺度；典型 1e9–1e11，被集流体主导）\n",
            σ_vm_global_max)
    @printf("  sanity: scale.E_coat=%.3e Pa（厚度加权有效模量；PE/NE.E_coat=%.2e/%.2e，PCC/NCC.E=%.2e/%.2e）\n",
            case.param_dim.scale.E_coat,
            case.param_dim.PE.E_coat, case.param_dim.NE.E_coat,
            case.param_dim.PCC.E, case.param_dim.NCC.E)

    # 若所有节点位移为零，疑似力学求解失败
    @assert Umag_global_max > 0 "所有节点位移均为零，疑似力学求解失败（K_mech\\F_mech 异常），请检查 src/Mechanical.jl:289 catch 分支是否触发"

    # 自适应变形放大系数：使最大变形 ≈ 5% L_ref
    DEF_SCALE = 0.05 * L_ref / Umag_global_max
    @printf("  DEF_SCALE=%.3e（自适应，使 |U|_max×DEF_SCALE ≈ 5%% L_ref）\n", DEF_SCALE)

    # ====================================================================
    # 6. 绘图：一张大图，行=场分量，列=时间节点
    # ====================================================================
    println("\n[绘图]")

    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)

    nT = length(fields_per_ti)
    row_labels = ["σ_xx [MPa]", "σ_yy [MPa]", "σ_xy [MPa]", "σ_vm [MPa]", "|U| [μm]"]
    nrow = length(row_labels)
    subplots = Plots.Plot[]

    xs0, ys0 = q4_element_polygons(mesh_th)

    for col in 1:nT
        pack = fields_per_ti[col]
        Umag = sqrt.(pack.U_x.^2 .+ pack.U_y.^2)
        # 节点量 |U| 映射到单元常数场（4 角点平均）
        Umag_elem = zeros(ne)
        for e in 1:ne
            n1, n2, n3, n4 = mesh_th.element[e, 1], mesh_th.element[e, 2],
                             mesh_th.element[e, 3], mesh_th.element[e, 4]
            Umag_elem[e] = 0.25 * (Umag[n1] + Umag[n2] + Umag[n3] + Umag[n4])
        end
        # 单位换算到 MPa / μm
        col_fields = [pack.σ_xx*1e-6, pack.σ_yy*1e-6, pack.σ_xy*1e-6,
                      pack.σ_vm*1e-6, Umag_elem*1e6]
        for row in 1:nrow
            p = plot(legend=false)
            if row < 5
                clim = (minimum(col_fields[row]), maximum(col_fields[row]))
                plot_q4_field!(p, mesh_th, col_fields[row], xs0, ys0;
                    title=(@sprintf("t=%.0f s, %s", pack.t, row_labels[row])),
                    cmap=:RdBu, clim=clim)
            else
                clim = (minimum(col_fields[row]), maximum(col_fields[row]))
                plot_q4_field!(p, mesh_th, col_fields[row], xs0, ys0;
                    title=(@sprintf("t=%.0f s, %s (def ×%.0f)", pack.t, row_labels[row], DEF_SCALE)),
                    cmap=:viridis, clim=clim,
                    deform_xy=(pack.U_x, pack.U_y), def_scale=DEF_SCALE)
            end
            push!(subplots, p)
        end
    end

    p_combined = plot(subplots..., layout=(nrow, nT),
                      size=(450*nT, 350*nrow), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
    out_path = joinpath(output_dir, "jellyroll_stress_displacement.png")
    savefig(p_combined, out_path)
    @printf("  大图已保存: %s\n", out_path)

    println("\n" * "="^80)
    @printf("完成：σ_vm 全局峰值 = %.2f MPa（极片/电极尺度）\n", σ_vm_global_max*1e-6)
    println("="^80)

    # === ANCHOR_END_PARAMS ===
end

main()
