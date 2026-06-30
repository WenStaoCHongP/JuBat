"""
测试案例：Jellyroll 电池单次放电二维应力场/位移场后处理

功能：
- 多 SPMe 并行 + 二维分布式热模型（与 testexample 同工况，但关闭 CZM）
- 在用户指定的时间节点调用 thermal_diffusion_stress_2D 计算宏观应力/位移场
- 输出一张大图：行 = 场分量（σ_xx/σ_yy/σ_xy/σ_vm/|U|），列 = 时间节点

⚠️ 尺度说明：本脚本输出的是 **极片/电极尺度（coating-scale, 二维平面应力）** 的场，
   由全叠合厚度加权有效模量 E_eff（~5e8 Pa 量级）计算；
   与颗粒尺度 Calstressdisp（颗粒 E ~1e10 Pa）不同，二者不可混用。

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
function plot_q4_field!(plt, mesh, field_elem, xs, ys; title="", cmap=:RdYlBu_r,
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

    # === ANCHOR_END_PARAMS ===
end

main()
