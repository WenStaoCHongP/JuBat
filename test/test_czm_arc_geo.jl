using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# Batch 5 Task 2：Crisfield 柱面弧长 geo 路径（Theory §6.10，λ 缩放本征应变增量）

function arc_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    pc = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, pc)
    return case, pc, cache
end

@testset "弧长 geo 弹性区与 basic 一致且 λ 达 1" begin
    case, pc, cache = arc_fixture()
    cm = case.czm_mesh
    ne = size(cm.bulk_element, 1)
    ndof = 2 * cm.nnode
    eig = (dT = fill(1e-6, ne), Δsn = zeros(ne), Δsp = zeros(ne))
    kw = (dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp)
    r_basic, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
        kw..., max_iter = 100, tol = 1e-10, iter_method = "basic", cache = cache,
        geo_nl = true, eigenstrain = eig)
    r_arc, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
        kw..., max_iter = 100, tol = 1e-8, n_load_steps = 10, iter_method = "arc_length",
        cache = cache, geo_nl = true, eigenstrain = eig)
    @test r_basic.converged && r_arc.converged
    @test isfinite(r_arc.residual_norm) && r_arc.residual_norm < 1e-6
    @test isapprox(r_arc.displacement, r_basic.displacement; rtol = 1e-6, atol = 1e-12)
end

@testset "geo 弧长系数进入约束且非法值显式失败" begin
    case, pc, cache = arc_fixture()
    cm = case.czm_mesh
    ne = size(cm.bulk_element, 1)
    ndof = 2 * cm.nnode
    eig = (dT = fill(1e-6, ne), Δsn = zeros(ne), Δsp = zeros(ne))
    kw = (dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp,
          max_iter = 100, tol = 1e-8, n_load_steps = 10, iter_method = "arc_length",
          cache = cache, geo_nl = true, eigenstrain = eig)
    r_half, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
        kw..., arc_length_alpha = 0.5)
    @test r_half.converged
    @test_throws ArgumentError JuBat.solve_czm_step(
        cm, zeros(ndof), pc, case.param, zeros(ndof);
        kw..., arc_length_alpha = 0.0)
end

@testset "geo 弧长先迭代平衡自由芯部卷绕预应力参考态" begin
    case, pc, _ = arc_fixture()
    cm = case.czm_mesh
    ne = size(cm.bulk_element, 1)
    ndof = 2 * cm.nnode
    cache = JuBat.ensure_czm_cache(case, cm, pc; fix_inner=false)
    prestress_full = JuBat.winding_prestress_field(cm, case.param)
    prestress = [(0.2a, 0.2b, 0.2c) for (a, b, c) in prestress_full]
    eig = (dT = fill(1e-6, ne), Δsn = zeros(ne), Δsp = zeros(ne))

    result, _ = JuBat.solve_czm_step(
        cm, zeros(ndof), pc, case.param, zeros(ndof);
        dT_elem = eig.dT, Δsoc_n_elem = eig.Δsn, Δsoc_p_elem = eig.Δsp,
        max_iter = 100, tol = 1e-8, n_load_steps = 10, iter_method = "arc_length",
        cache = cache, geo_nl = true, eigenstrain = eig, prestress = prestress)
    @test result.converged
    @test result.residual_norm <= 1e-8
end

@testset "geo_nl=false + arc_length 行为不变（回归锚）" begin
    case, pc, cache = arc_fixture()
    cm = case.czm_mesh
    ne = size(cm.bulk_element, 1)
    ndof = 2 * cm.nnode
    dT = fill(1e-4, ne)
    r1, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
        dT_elem = dT,
        Δsoc_n_elem = zeros(ne), Δsoc_p_elem = zeros(ne),
        max_iter = 100, tol = 1e-10, n_load_steps = 10, iter_method = "arc_length",
        cache = cache)
    r2, _ = JuBat.solve_czm_step(cm, zeros(ndof), pc, case.param, zeros(ndof);
        dT_elem = dT,
        Δsoc_n_elem = zeros(ne), Δsoc_p_elem = zeros(ne),
        max_iter = 100, tol = 1e-10, n_load_steps = 10, iter_method = "arc_length",
        cache = cache, geo_nl = false)
    @test r1.converged == r2.converged
    @test r1.displacement == r2.displacement   # 逐位（geo_nl=false 路径零漂移）
end

@testset "球面弧长修正跨越合成极限点" begin
    # 一自由度软化平衡路径 λ = u - u³，在 u=1/√3 处有极限点。
    alpha = 0.2
    radius = 0.05
    u = 0.0
    lambda = 0.0
    previous_tangent = Float64[]
    u_hist = [u]
    lambda_hist = [lambda]

    for _ in 1:20
        K = 1.0 - 3.0 * u^2
        tangent = 1.0 / K
        direction = isempty(previous_tangent) || dot([tangent, alpha], previous_tangent) >= 0 ? 1.0 : -1.0
        dlambda = direction * radius / hypot(tangent, alpha)
        u0, lambda0 = u, lambda
        u += tangent * dlambda
        lambda += dlambda
        for _ in 1:30
            R = lambda - (u - u^3)
            du_bar = [u - u0]
            dlambda_bar = lambda - lambda0
            g = dot(du_bar, du_bar) + alpha^2 * dlambda_bar^2 - radius^2
            abs(R) < 1e-12 && abs(g) < 1e-12 && break
            K = 1.0 - 3.0 * u^2
            du_R = [R / K]
            du_F = [1.0 / K]
            du, dlambda_corr, _ = JuBat.spherical_arc_length_correction(
                du_bar, dlambda_bar, du_R, du_F, alpha, radius)
            u += du[1]
            lambda += dlambda_corr
        end
        @test abs(lambda - (u - u^3)) < 1e-10
        step_tangent = [u - u0, alpha * (lambda - lambda0)]
        previous_tangent = step_tangent ./ norm(step_tangent)
        push!(u_hist, u)
        push!(lambda_hist, lambda)
    end

    imax = argmax(lambda_hist)
    @test 1 < imax < length(lambda_hist)
    @test maximum(u_hist) > 0.7
    @test any(diff(lambda_hist) .< 0.0)
end
