"""
CZM 热点定位探针（只读诊断，不改求解路径）

用 testexample 配置跑一次带采样 profile 的完整求解，按源码行聚合采样数，
回答一个问题：CZM 的时间到底花在"逐单元循环装配"上，还是花在"每次 Newton
迭代重新做稀疏 LU 分解"上。

背景：docs/planning-with-files/29_仿真提速-SPMe热CZM 在 2026-08-19 测得
84% 采样落在 CzmSolve.jl 的 `K_bc \\ R_bc`，且装配向量化批次（Task 7）实测
仅 -1.7%、随后被回滚。本探针在 2026-08-30 四层重构后的代码上复核该结论。

日期：2026-08-31
"""

using Printf
using Profile
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function build_case()
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2

    opt = JuBat.Option()
    opt.Current = x -> 5.0
    opt.model = "SPMe"
    opt.Nn = 10; opt.Ns = 5; opt.Np = 10
    opt.Nrn = 10; opt.Nrp = 10
    opt.gsorder = 2
    opt.dimension = 1
    opt.mechanicalmodel = "none"
    opt.time = [0.0, 60]
    opt.dt = [0.5, 10]
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.solveType = "Crank-Nicolson"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.thermal_dim = "2D"
    opt.cool_method = "surface"
    opt.per_element_spme = true
    opt.czm.enabled = true
    opt.czm.model = "mix"
    opt.czm.fix_inner = false
    opt.czm.iter_method = "basic"
    opt.czm.load_steps = 10
    opt.czm.tol = 1e-3

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    return case
end

function main()
    println("="^78)
    println("CZM 热点定位（当前代码，四层重构后）")
    println("="^78)

    Profile.init(n = 20_000_000, delay = 0.005)
    case = build_case()
    Profile.clear()
    @profile JuBat.Solve(case)

    data, lidict = Profile.retrieve()
    total = 0
    counts = Dict{String,Int}()
    for ip in data
        ip == 0 && continue
        frames = get(lidict, ip, nothing)
        frames === nothing && continue
        for fr in frames
            f = string(fr.file)
            (occursin("JuBat", f) && endswith(f, ".jl")) || continue
            key = @sprintf("%s:%d", basename(f), fr.line)
            counts[key] = get(counts, key, 0) + 1
            total += 1
        end
    end

    @printf("\n有效 JuBat 采样帧: %d\n\n", total)
    println("按源码行聚合的前 25 位：")
    @printf("  %-38s %8s %8s\n", "位置", "采样", "占比")
    println("  " * "-"^56)
    for (k, v) in first(sort(collect(counts), by = x -> -x[2]), 25)
        @printf("  %-38s %8d %7.2f%%\n", k, v, 100v / total)
    end

    # 按文件聚合
    byfile = Dict{String,Int}()
    for (k, v) in counts
        byfile[split(k, ':')[1]] = get(byfile, split(k, ':')[1], 0) + v
    end
    println("\n按文件聚合：")
    for (k, v) in first(sort(collect(byfile), by = x -> -x[2]), 12)
        @printf("  %-28s %8d %7.2f%%\n", k, v, 100v / total)
    end
    println("\n" * "="^78)
end

main()
