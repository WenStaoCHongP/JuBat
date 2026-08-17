# 耦合工况下三种 CZM 迭代方法收敛性对比（与 testexample.jl 同配置）
# 运行: julia tools/check_czm_methods_coupled.jl

using Printf
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

const METHODS = ["basic", "load_substep", "arc_length"]

function run_one(method::String)
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
    opt.mechanicalmodel = "none"
    opt.time = [0.0, 3600.0]
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
    opt.czm_enabled = true
    opt.czm_iter_method = method
    opt.czm_load_steps = 10
    opt.czm_tol = 1e-3
    opt.czm_max_iter = 100

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=16, nθ_czm=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    # 逐步求解并统计 CZM 收敛
    n_ok = 0
    n_fail = 0
    max_res = 0.0
    worst_step = 0

    # 复用 Solve 初始化逻辑：直接调用 Solve 并在内部无法拦截时，
    # 用 IO 重定向捕获 @warn —— 改为手动时间步循环太重，直接 Solve + 钩子
    result = JuBat.Solve(case)

    # Solve 不返回逐步 CZM 收敛统计；通过 result 时间步与 CZM 输出间接判断
    has_czm = haskey(result, "czm D_max")
    D_max_end = has_czm ? result["czm D_max"][end] : NaN
    n_steps = length(result["time [s]"])

    return (method=method, n_steps=n_steps, D_max_end=D_max_end, solved=true)
end

function main()
    println("=" ^ 70)
    println("耦合工况 CZM 三种方法对比 (testexample 配置, 3600s)")
    println("=" ^ 70)

    # 更可靠：在 update_czm_damage! 层统计 —— 用 eval 包装
    stats = Dict{String, NamedTuple}()

    for method in METHODS
        println("\n>>> 运行方法: $method")
        n_ok = Ref(0)
        n_fail = Ref(0)
        max_res = Ref(0.0)

        # 包装 update_czm_damage! 计数
        orig = JuBat.update_czm_damage!
        function counting_update!(args...; kwargs...)
            u, conv = orig(args...; kwargs...)
            if conv
                n_ok[] += 1
            else
                n_fail[] += 1
            end
            return u, conv
        end

        try
            # 临时替换（模块内函数）
            Core.eval(JuBat, :(update_czm_damage! = $counting_update!))

            param_dim = JuBat.ChooseCell("Jellyroll")
            param_dim.cell.v_l = 2.5
            param_dim.cell.v_h = 4.2
            opt = JuBat.Option()
            opt.Current = x -> 5.0
            opt.model = "SPMe"
            opt.Nn = 10; opt.Ns = 5; opt.Np = 10
            opt.Nrn = 10; opt.Nrp = 10
            opt.gsorder = 2
            opt.mechanicalmodel = "none"
            opt.time = [0.0, 3600.0]
            opt.dt = [0.5, 10.0]
            opt.dtType = "auto"
            opt.jacobi = "update"
            opt.solveType = "Crank-Nicolson"
            opt.thermal_enabled = true
            opt.thermalmodel = "distributed2D"
            opt.per_element_spme = true
            opt.debug_coupling = false
            opt.czm_enabled = true
            opt.czm_iter_method = method
            opt.czm_load_steps = 10
            opt.czm_tol = 1e-3
            opt.czm_max_iter = 100

            case = JuBat.SetCase(param_dim, opt)
            mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=16, nθ_czm=40, gsorder=2)
            case = JuBat.setup_thermal2D_mesh(case, mesh_data)
            case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

            t0 = time_ns()
            result = JuBat.Solve(case)
            wall_s = (time_ns() - t0) * 1e-9

            D_end = haskey(result, "czm D_max") ? result["czm D_max"][end] : NaN
            stats[method] = (
                n_ok=n_ok[],
                n_fail=n_fail[],
                n_total=n_ok[] + n_fail[],
                D_max_end=D_end,
                wall_s=wall_s,
                n_time_steps=length(result["time [s]"]),
            )
            @printf("  CZM 收敛: %d/%d 步, D_max(end)=%.4f%%, wall=%.1f s\n",
                n_ok[], n_ok[] + n_fail[], D_end * 100, wall_s)
        finally
            Core.eval(JuBat, :(update_czm_damage! = $orig))
        end
    end

    println("\n" * "=" ^ 70)
    println("汇总")
    println("=" ^ 70)
    @printf("%-15s | %-12s | %-10s | %s\n", "方法", "CZM收敛", "D_max终值", "墙钟[s]")
    println("-" ^ 55)
    all_ok = true
    for method in METHODS
        s = stats[method]
        conv_str = @sprintf("%d/%d", s.n_ok, s.n_total)
        ok = s.n_fail == 0
        all_ok &= ok
        mark = ok ? "✓" : "✗"
        @printf("%-15s | %-12s | %9.2f%% | %8.1f %s\n",
            method, conv_str, s.D_max_end * 100, s.wall_s, mark)
    end
    println()
    if all_ok
        println("结论: 三种方法在耦合工况下所有 CZM 步均收敛。")
    else
        println("结论: 存在未收敛的 CZM 步，详见上表。")
    end
end

main()
