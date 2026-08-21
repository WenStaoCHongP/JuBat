using Test
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec 2026-08-20-core-collapse-mechanics-design.md §4.2 / §7 Batch 1
# 开关全关时，新入口必须与既有 K_bulk*u 路径逐位等价；未实现的槽位必须报错。

function build_mech_core_fixture(; nθ::Int=40)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    param_cache = JuBat.compute_czm_params_per_interface(case)
    return case, param_cache
end

# 确定性的非零位移场：正弦扰动，避免全零掩盖 f_int 差异
function make_test_u(nnode::Int)
    u = zeros(Float64, 2 * nnode)
    for n in 1:nnode
        u[2*n-1] = 1e-6 * sinpi(2 * n / nnode)
        u[2*n]   = 1e-6 * cospi(2 * n / nnode)
    end
    return u
end

@testset "新入口与 K_bulk*u 逐位等价（线弹性槽位）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    K_ref = JuBat.assemble_bulk_stiffness(czm_mesh, param_cache)
    f_ref = K_ref * u

    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache)

    # 同一装配同一乘法，必须是逐位相等而非近似相等
    @test f_int == f_ref
    @test K_tan == K_ref
    @test length(f_int) == 2 * czm_mesh.nnode
    @test size(K_tan) == (2 * czm_mesh.nnode, 2 * czm_mesh.nnode)
    @test !any(isnan, f_int)
    @test !any(isnan, K_tan)
end

@testset "零位移给零内力（线弹性无预应力）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u0 = zeros(Float64, 2 * czm_mesh.nnode)

    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u0, param_cache)

    @test all(iszero, f_int)
    @test nnz(K_tan) > 0
end

@testset "缓存不变量：传入 K_bulk_cached 时直接复用同一对象" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    cache = JuBat.ensure_czm_cache(case, czm_mesh, param_cache)
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; K_bulk_cached=cache.K_bulk)

    # === 不是"数值相等"，而是同一对象：证明没有重复装配
    @test K_tan === cache.K_bulk
    @test f_int == cache.K_bulk * u
end

@testset "切线对称性（线弹性 + 平面应力各向同性）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    _, K_tan = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, param_cache)

    @test norm(K_tan - transpose(K_tan), Inf) ≤ 1e-12 * norm(K_tan, Inf)
end

@testset "未实现槽位必须报错，不得静默降级（AGENTS 9.7 / spec §6）" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh
    u = make_test_u(czm_mesh.nnode)

    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; geo_nl=true)
    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache; plasticity=true)
    # mech_state 的消费者在 Batch 3 引入；此前传入非 nothing 即报错
    @test_throws ErrorException JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache, :dummy_state)
end

@testset "位移向量长度不符必须报错，不得截断或补零" begin
    case, param_cache = build_mech_core_fixture()
    czm_mesh = case.czm_mesh

    u_short = zeros(Float64, 2 * czm_mesh.nnode - 1)
    @test_throws DimensionMismatch JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u_short, param_cache)

    u_long = zeros(Float64, 2 * czm_mesh.nnode + 3)
    @test_throws DimensionMismatch JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u_long, param_cache)
end
