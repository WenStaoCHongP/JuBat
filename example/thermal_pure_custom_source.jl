using Plots
gr()  # 确保无GUI环境也可保存图片
using Dates
include("../src/JuBat.jl")

# 纯热分布式2D示例：自定义平均体积热源，不接电化学
# 目标：生成温度场分布图和热应力（元素等效）分布图

# ---------------- 用户可调参数（亦可用环境变量覆盖） ----------------
cell_name     = get(ENV, "JUBAT_CELL", "Jellyroll")
total_time_s  = try parse(Float64, get(ENV, "JUBAT_TIME", "60")) catch; 60.0 end
dt_s          = try parse(Float64, get(ENV, "JUBAT_DT", "1.0")) catch; 1.0 end
n_theta       = try parse(Int, get(ENV, "JUBAT_NTHETA", "180")) catch; 180 end
q_source_Wm3  = try parse(Float64, get(ENV, "JUBAT_Q", "5e4")) catch; 5e4 end  # 常量体热源 [W/m^3]
units_thermal = get(ENV, "JUBAT_UNITS", "SI")  # "SI" 或 "nd"
hc_override   = try parse(Float64, get(ENV, "JUBAT_H", "0")) catch; 0.0 end  # 若>0，覆盖对流换热系数 h
function main()
    # ---------------- 初始化案例与热网格 ----------------
    param_dim = JuBat.ChooseCell(cell_name)
    # 覆盖对流换热（在归一化前设置）
    if hc_override > 0
        try
            param_dim.cell.h = hc_override
        catch
        end
    end
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.collector_seeded = true
    opt.units_thermal = (uppercase(units_thermal) == "SI") ? "SI" : "nd"

    # 构建热网格（collector-seeded 条带），并挂到 case.mesh
    case = JuBat.SetCase(param_dim, opt)
    begin
        nθ = max(12, n_theta)
        mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
        case.mesh["thermal2D"] = mesh_th
    end

    # ---------------- 初始温度场与自定义热源 ----------------
    variables = Dict{String, Union{Array{Float64}, Float64}}()
    mesh = case.mesh["thermal2D"]
    # 无量纲初温：T0* = T0_SI / T_ref
    T0_nd = param_dim.cell.T0 / param_dim.scale.T_ref
    variables["T_nodes"] = fill(T0_nd, mesh.nlen)  # 以环境初温为初始（nd）

    # 为 collector-seeded 网格自动提供层权重，以便各向异性/ρc聚合
    try
        fks = JuBat.jellyroll_get_layer_weights(mesh)
        if fks !== nothing
            variables["thermal2D layer_weights"] = fks
        else
            error("网格不是通过 jellyroll_collector_seed_mesh 生成，无法获取层权重")
        end
    catch
    end

    # 构造用户自定义体热源（元素常量）
    ne = size(mesh.element, 1)
    q_elem = fill(q_source_Wm3, ne)
    if opt.units_thermal == "SI"
        variables["heat_source_fields"] = q_elem
        variables["heat_source_units_code"] = 1.0  # 标记为SI单位
    else
        # 无量纲：q* = q / q_ref
        q_ref = param_dim.scale.q_th
        variables["heat_source_fields"] = q_elem ./ max(q_ref, 1e-16)
        variables["heat_source_units_code"] = 0.0
    end

    # ---------------- 时间推进（仅热） ----------------
    T_nodes = variables["T_nodes"]::Vector{Float64}
    times = collect(0.0:dt_s:total_time_s)

    # Scheme B 缩放：dt_th* = dt_elec / t_th；此示例无电化学，直接用 t_ratio = t0/t_th
    scale = param_dim.scale
    t_ratio = scale.t0 / max(scale.t_th, 1e-16)
    dt_nd = dt_s / max(scale.t0, 1e-16)  # 与主框架一致：电化学时间归一化
    dt_th = dt_nd * t_ratio

    for (_k, _t) in enumerate(times)
        MT, KT, FT = JuBat.ThermalDistributed2D(case, variables)
        JuBat.ThermalDistributed2D_BC(KT, FT, case, 0.0)
        A = (1.0/dt_th) .* MT + KT
        rhs = (1.0/dt_th) .* (MT * T_nodes) + FT
        # 轻微对角正则，锚定到 T_amb
        alpha = 1e-12
        if alpha > 0
            for i in 1:size(A,1)
                A[i,i] += alpha
            end
            T_amb_nd = param_dim.cell.T_amb / scale.T_ref
            rhs .+= alpha .* T_amb_nd
        end
    T_new = A \ rhs
        variables["T_nodes"] = T_new
        T_nodes = T_new
    end

    # 为兼容 SPM/SPMe/P2D 的分支，提供一个标量“平均温度”（无量纲）键
    # 简易平均温度（无量纲）
    let Tn = variables["T_nodes"]
        if isa(Tn, AbstractVector) && length(Tn) > 0
            variables["temperature"] = sum(Tn) / length(Tn)
        end
    end
    # 计算二维单元热应力（由 mechanical.jl 的 thermal_stress 提供）
    variables = JuBat.thermal_stress(case, variables)

    # ---------------- 绘图与输出 ----------------
    outdir = joinpath(@__DIR__, "..", "output")
    isdir(outdir) || mkpath(outdir)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")

    # 温度场（节点散点）
    x = mesh.node[:,1]; y = mesh.node[:,2]
    T_K = variables["T_nodes"] .* scale.T_ref
    pT = scatter(x, y; marker_z=T_K, ms=6, mc=:inferno, markerstrokewidth=0, colorbar=true,
                 title="Temperature field [K]", aspect_ratio=1)
    savefig(pT, joinpath(outdir, "thermal_pure_T_" * timestamp * ".png"))

    # 热应力（元素中心散点）
    centers = JuBat.jellyroll_element_centers(mesh)
    σe = haskey(variables, "thermal2D element thermal stress") ? variables["thermal2D element thermal stress"] : zeros(size(centers,1))
    pS = scatter(centers[:,1], centers[:,2]; marker_z=σe, ms=7, mc=:viridis, markerstrokewidth=0, colorbar=true,
                 title="Thermal stress (element) [Pa]", aspect_ratio=1)
    savefig(pS, joinpath(outdir, "thermal_pure_stress_" * timestamp * ".png"))

    println("Saved outputs to: ", outdir)
    println("Files:")
    println(" - ", joinpath(outdir, "thermal_pure_T_" * timestamp * ".png"))
    println(" - ", joinpath(outdir, "thermal_pure_stress_" * timestamp * ".png"))
end

main()
