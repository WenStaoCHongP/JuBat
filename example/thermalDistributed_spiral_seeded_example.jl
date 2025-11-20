using LinearAlgebra, SparseArrays, Plots, Statistics, DelimitedFiles
include("../src/JuBat.jl")

"""
验证示例：使用“螺旋布种(collector-seeded)”条带热网格对分层热模型进行瞬态验证。
要点：
- 网格：`jellyroll_Q4_mesh(...; crop_mode=:collector_seeded)`，每个 Q4 元素跨越完整层序，自动绑定元素层权重 f_k。
- 热模型：`ThermalDistributed2D`，利用 f_k 在装配中聚合 `ρc` 与各向异性 k；并在示例中自定义分层体热源合成元素平均 q_e。
- 时间推进：后退欧拉 (Backward Euler)。
- 输出：时间历程、最终温度场、热应力分布（元素中心）及用于外部绘图的 CSV。
运行：在仓库根目录执行
    julia --project example/thermalDistributed_spiral_seeded_example.jl
"""
function main()
    # 1) 参数与选项
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPM"             # 仅为构造参数使用；本示例只跑热模块
    opt.thermalmodel = "distributed2D"
    opt.thermal_enabled = true
    opt.thermal_dim = "2D"
    # 启用 collector-seeded 逻辑（对 Solve/装配无强制要求，此处用于一致性）
    opt.collector_seeded = true

    # 2) 构造案例与条带热网格（collector-seeded）
    case = JuBat.SetCase(param_dim, opt)
    # 这里 nx 被解释为每圈分段数（见 Jellyrollmodel.jl 实现注释）
    seg_per_turn = 64
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=seg_per_turn, gsorder=2)
    case.mesh["thermal2D"] = mesh_th

    # 3) 变量与热源：利用网格内置的元素层权重 f_k 合成 q_e（SI → 无量纲由装配处理）
    vars = JuBat.StandardVariables(case, 1)
    # 获取与该 mesh 关联的 f_k（ne×5，顺序 [NE,SP,PE,PCC,NCC]）
    fks = JuBat.jellyroll_get_layer_weights(mesh_th)
    if fks === nothing
        error("网格不是通过 jellyroll_collector_seed_mesh 生成，无法获取层权重")
    end
    ne = size(mesh_th.element, 1)
    # 设定分层体热源（W/m^3）：给不同层不同强度，便于验证 f_k 聚合效果
    # 数值可自由调整；此处故意让 PE > NE > SP，集流体较小
    q_NE = 8.0e1; q_SP = 2.0e0; q_PE = 1.0e1; q_PCC = 5.0e-1; q_NCC = 5.0e-1
    q_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        f_NE, f_SP, f_PE, f_PCC, f_NCC = fks[e,1], fks[e,2], fks[e,3], fks[e,4], fks[e,5]
        q_elem[e] = f_NE*q_NE + f_SP*q_SP + f_PE*q_PE + f_PCC*q_PCC + f_NCC*q_NCC
    end
    vars["heat_source_fields"] = q_elem
    vars["heat_source_units_code"] = 1.0  # SI → 由 ThermalDistributed2D 内部按 q_ref 无量纲化
    # 同步给装配使用的层权重（使 ρc/k 聚合也走 f_k 路径）
    vars["thermal2D layer_weights"] = fks

    # 4) 装配与边界条件（常量）
    MT, KT, FT = JuBat.ThermalDistributed2D(case, vars)
    JuBat.ThermalDistributed2D_BC(KT, FT, case, 0.0)

    # 5) 时间推进（后退欧拉）
    nnode = mesh_th.nlen
    T0 = fill(case.param.cell.T0, nnode)     # 无量纲初温
    dt = 0.2                                 # 热学无量纲时间步（参照 scale.t_th）
    Tend = 10.0
    nsteps = Int(ceil(Tend/dt))

    A = (1.0/dt) .* MT + KT
    fac = lu(A)

    Ts_mean = zeros(Float64, nsteps+1)
    Ts_max  = zeros(Float64, nsteps+1)
    Times   = collect(0:nsteps) .* dt

    T = copy(T0)
    Ts_mean[1] = mean(T)
    Ts_max[1]  = maximum(T)
    for k in 1:nsteps
        rhs = (1.0/dt) .* (MT * T) + FT
        T = fac \ rhs
        Ts_mean[k+1] = mean(T)
        Ts_max[k+1]  = maximum(T)
    end

    # 6) 绘图：温度历程与最终场（K）
    Tref = case.param_dim.scale.T_ref
    Ts_mean_K = Ts_mean .* Tref
    Ts_max_K  = Ts_max  .* Tref
    p1 = plot(Times, Ts_mean_K, label="mean(T)", xlabel="time (nd)", ylabel="Temperature [K]")
    plot!(p1, Times, Ts_max_K, label="max(T)")
    savefig(p1, "thermalDistributed_spiral_seeded_history.png")

    # 按元素中心的平均温度绘色（与应力图一致的风格）
    nelem = size(mesh_th.element, 1)
    Te = zeros(Float64, nelem)
    @inbounds for e in 1:nelem
        nodes = mesh_th.element[e, :]
        Te[e] = mean(T[nodes])
    end
    centers = JuBat.jellyroll_element_centers(mesh_th)
    cx, cy = centers[:,1], centers[:,2]
    Te_K = Te .* Tref
    p2 = scatter(cx, cy, marker_z=Te_K, ms=4, c=:turbo, xlabel="x [m]", ylabel="y [m]", aspect_ratio=1,
                 colorbar_title="T [K]", title="Final T field (collector-seeded, element-averaged)")
    savefig(p2, "thermalDistributed_spiral_seeded_topview.png")

    println("Saved: thermalDistributed_spiral_seeded_history.png")
    println("Saved: thermalDistributed_spiral_seeded_topview.png")

    # 7) 按元素中心绘制“热应力”对比（调用库内方法；如缺失则回退等效计算）
    vars["T_nodes"] = T
    vars["temperature"] = [mean(T)]
    σe = nothing
    try
        vars2 = JuBat.thermal_stress(case, vars)
        σe = vars2["thermal2D element thermal stress"]
    catch err
        @warn "thermal_stress not available; using fallback" err
        param = case.param
        Tref = param.scale.T_ref
        T0d = param.cell.T0
        α_n = hasproperty(param.NE, :alphaT) ? getfield(param.NE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        α_p = hasproperty(param.PE, :alphaT) ? getfield(param.PE, :alphaT) : (hasproperty(param.cell, :alphaT) ? param.cell.alphaT : 0.0)
        E_n = hasproperty(param.NE, :E) ? getfield(param.NE, :E) : 0.0
        E_p = hasproperty(param.PE, :E) ? getfield(param.PE, :E) : 0.0
        ν_n = hasproperty(param.NE, :nu) ? getfield(param.NE, :nu) : 0.0
        ν_p = hasproperty(param.PE, :nu) ? getfield(param.PE, :nu) : 0.0
        wt_den = max(1e-12, (param.NE.thickness + param.PE.thickness))
        E_eff = (E_n * param.NE.thickness + E_p * param.PE.thickness) / wt_den
        ν_eff = (ν_n * param.NE.thickness + ν_p * param.PE.thickness) / wt_den
        α_eff = (α_n * param.NE.thickness + α_p * param.PE.thickness) / wt_den
        σe = zeros(Float64, size(mesh_th.element,1))
        centers = JuBat.jellyroll_element_centers(mesh_th)
        for e in 1:size(mesh_th.element,1)
            nodes = mesh_th.element[e, :]
            Te = mean(T[nodes])
            dT_K = (Te - T0d) * Tref
            σe[e] = E_eff * α_eff * dT_K / max(1e-12, (1.0 - ν_eff))
        end
    end
    centers = JuBat.jellyroll_element_centers(mesh_th)
    cx, cy = centers[:,1], centers[:,2]
    p3 = scatter(cx, cy, marker_z=σe, ms=4, c=:balance, xlabel="x [m]", ylabel="y [m]", aspect_ratio=1,
                 colorbar_title="σ_th [Pa]", title="Thermal stress (collector-seeded)")
    savefig(p3, "thermalDistributed_spiral_seeded_stress_topview.png")
    println("Saved: thermalDistributed_spiral_seeded_stress_topview.png")

    # 8) 导出 CSV
    try
        outdir = @__DIR__
    # nodes: id,x,y
    x = mesh_th.node[:,1]; y = mesh_th.node[:,2]
    nodes_tbl = hcat(collect(1:nnode), x, y)
        writedlm(joinpath(outdir, "spiral_seed_nodes.csv"), nodes_tbl, ',')
        # elements: id,n1,n2,n3,n4
        nelem = size(mesh_th.element,1)
        elems_tbl = hcat(collect(1:nelem), mesh_th.element)
        writedlm(joinpath(outdir, "spiral_seed_elements.csv"), elems_tbl, ',')
    # T_nodes in Kelvin: id, T_K（保留节点温度的 CSV 以兼容外部绘图）
    T_K_nodes = T .* Tref
    T_tbl = hcat(collect(1:nnode), T_K_nodes)
        writedlm(joinpath(outdir, "spiral_seed_T_nodes_K.csv"), T_tbl, ',')
        # element thermal stress: id, sigma
        σ_tbl = hcat(collect(1:nelem), σe)
        writedlm(joinpath(outdir, "spiral_seed_stress_elem.csv"), σ_tbl, ',')
        println("Saved CSVs: spiral_seed_nodes.csv, spiral_seed_elements.csv, spiral_seed_T_nodes_K.csv, spiral_seed_stress_elem.csv")
    catch err
        @warn "Failed to export CSV for plotting" err
    end
end

main()
