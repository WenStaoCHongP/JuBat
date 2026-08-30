using Test
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §7 Batch 2：切线 FD、大转角刚体运动零内力（D9 精确性质）、
# 线性退化（u=0 切线 ≡ 线性刚度）、自由膨胀零应力（D-B2-1）、K_G 压缩方向性。
# 夹具传 case.param（归一化网格，对齐生产）。

function build_geo_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    return case
end

# 刚体旋转位移场：u = (R(φ) − I)·X（相对参考构型的大转角）
function rigid_rotation_u(node::Matrix{Float64}, φ::Float64)
    c, s = cos(φ), sin(φ)
    u = zeros(Float64, 2 * size(node, 1))
    for n in 1:size(node, 1)
        x, y = node[n, 1], node[n, 2]
        u[2*n-1] = (c - 1) * x - s * y
        u[2*n]   = s * x + (c - 1) * y
    end
    return u
end

# 均匀膨胀/压缩位移场：u = κ·X（κ 取 √(1+2ε₀)−1 时为精确自由膨胀幅值）
function uniform_scale_u(node::Matrix{Float64}, κ::Float64)
    u = zeros(Float64, 2 * size(node, 1))
    for n in 1:size(node, 1)
        u[2*n-1] = κ * node[n, 1]
        u[2*n]   = κ * node[n, 2]
    end
    return u
end

@testset "大转角刚体运动零内力（完全 GL 精确性质，D9）" begin
    case = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u = rigid_rotation_u(czm_mesh.node, 30.0 * pi / 180)
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, case.param; geo_nl=true)
    @test norm(f_int) ≤ 1e-10 * norm(K_tan, Inf) * norm(u, Inf)
    @test !any(isnan, f_int) && !any(isnan, K_tan)
end

@testset "自由膨胀零应力（D-B2-1：ε₀ 均匀 + 精确膨胀位移）" begin
    case = build_geo_fixture()
    czm_mesh = case.czm_mesh
    ne = size(czm_mesh.bulk_element, 1)
    param = case.param
    ε₀ = 1e-3
    κ = sqrt(1.0 + 2 * ε₀) - 1.0
    u = uniform_scale_u(czm_mesh.node, κ)
    # α/β 分层化后（eigenstrain_of），全层均匀 ε₀ 只能由统一 alphaT 构造：
    # 测试内临时统一五层 alphaT（夹具每次重建，恢复只为本 testset 内卫生）
    α_backup = (param.PE.alphaT, param.NE.alphaT, param.SP.alphaT,
                param.PCC.alphaT, param.NCC.alphaT)
    α_uniform = param.PE.alphaT
    for layer in (param.PE, param.NE, param.SP, param.PCC, param.NCC)
        layer.alphaT = α_uniform
    end
    eigenstrain = (dT=fill(ε₀ / α_uniform, ne), Δsn=zeros(ne), Δsp=zeros(ne))
    f_int, K_tan = JuBat.assemble_bulk_residual_tangent(
        czm_mesh, u, case.param; geo_nl=true, eigenstrain=eigenstrain)
    param.PE.alphaT = α_backup[1]; param.NE.alphaT = α_backup[2]; param.SP.alphaT = α_backup[3]
    param.PCC.alphaT = α_backup[4]; param.NCC.alphaT = α_backup[5]
    @test norm(f_int) ≤ 1e-10 * norm(K_tan, Inf) * norm(u, Inf)
end

@testset "零位移 GL 切线退化为线性刚度（patch/线性极限）" begin
    case = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u0 = zeros(Float64, 2 * czm_mesh.nnode)
    _, K_geo = JuBat.assemble_bulk_residual_tangent(czm_mesh, u0, case.param; geo_nl=true)
    K_lin = JuBat.assemble_bulk_stiffness(czm_mesh, case.param)
    @test isapprox(Array(K_geo), Array(K_lin); rtol=1e-12)
end

@testset "切线有限差分（中等位移，含 K_G 与大梯度项）" begin
    case = build_geo_fixture()
    czm_mesh = case.czm_mesh
    u = zeros(Float64, 2 * czm_mesh.nnode)
    for n in 1:czm_mesh.nnode
        u[2*n-1] = 5e-3 * sin(3.0 * czm_mesh.node[n, 1] + 1.0)
        u[2*n]   = 5e-3 * cos(2.0 * czm_mesh.node[n, 2])
    end
    f0, K = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, case.param; geo_nl=true)
    h = 1e-7
    ndof = 2 * czm_mesh.nnode
    for dof in [1, 2, 101, div(ndof, 2), ndof - 1]
        up = copy(u); up[dof] += h
        um = copy(u); um[dof] -= h
        fp, _ = JuBat.assemble_bulk_residual_tangent(czm_mesh, up, case.param; geo_nl=true)
        fm, _ = JuBat.assemble_bulk_residual_tangent(czm_mesh, um, case.param; geo_nl=true)
        fd = (fp .- fm) ./ (2 * h)
        @test isapprox(collect(K[:, dof]), fd; rtol=1e-6, atol=1e-8 * max(1.0, norm(fd)))
    end
end

@testset "K_G 方向性：均匀压缩降低压缩向刚度" begin
    case = build_geo_fixture()
    czm_mesh = case.czm_mesh
    κ = -2e-3   # 压缩预载
    u = uniform_scale_u(czm_mesh.node, κ)
    d = uniform_scale_u(czm_mesh.node, 1.0)  # 压缩方向试探场
    K_lin = JuBat.assemble_bulk_stiffness(czm_mesh, case.param)
    _, K_geo = JuBat.assemble_bulk_residual_tangent(czm_mesh, u, case.param; geo_nl=true)
    @test dot(d, K_geo * d) < dot(d, Array(K_lin) * d)
end
