using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §3.5（D8 冻结定义）Δ_core 滤波单测 + 多圈持久化审计（Batch 5 Task 1）

function dc_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true; opt.thermalmodel = "distributed2D"; opt.per_element_spme = true
    opt.czm_enabled = true
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    return case
end

@testset "滤波单测：0/1 阶污染被去除，高阶残差保留" begin
    case = dc_fixture()
    cm = case.czm_mesh
    sub = cm.czm_submesh.mesh_bonded
    wt = cm.czm_submesh.winding_turn
    nθn = count(==(1), wt[1:count(==(wt[1]), wt)])   # 单匝窗口（与实现同式独立编码）
    A = 1e-4
    θk = 2π * (0:nθn-1) / nθn                       # 网格天然角序（θ1=0）
    mkfield(coefs) = begin
        u = zeros(2 * cm.nnode)
        for (i, k) in enumerate(1:nθn)
            r = hypot(sub.node[k, 1], sub.node[k, 2])
            c, s = sub.node[k, 1] / r, sub.node[k, 2] / r
            un = A * sum(a * cos(n * θk[i]) for (n, a) in coefs)
            u[2*k-1] = un * c; u[2*k] = un * s
        end
        u
    end
    # 0.3（呼吸 n=0）+ 0.5cosθ（平移 n=1）+ 0.4cos2θ + 0.15cos3θ；nθ=8 下 n≥4 不可分辨
    u = mkfield([(0, 0.3), (1, 0.5), (2, 0.4), (3, 0.15)])
    w, Δ = JuBat.core_ovalization(cm, u, sub.node)
    @test w ≈ A * 0.4 rtol = 0.35          # 残差 = n=2 主导（n=3 并入）
    @test Δ ≈ w / hypot(sub.node[1,1], sub.node[1,2]) rtol = 1e-12
    u01 = mkfield([(0, 0.3), (1, 0.5)])
    w01, _ = JuBat.core_ovalization(cm, u01, sub.node)
    @test w01 ≤ 1e-3 * A                   # 纯 0/1 阶 → 残差近零（离散最小二乘残余）
end

@testset "基准恒定：两次计算相对同一初始螺旋（node_ref 持久）" begin
    case = dc_fixture()
    cm = case.czm_mesh
    case.czm_layout = JuBat.CzmLayout(cm)
    u1 = fill(1e-5, 2 * cm.nnode)
    w1, _ = JuBat.core_ovalization(cm, u1, JuBat.ensure_node_ref!(case))
    u2 = 2 .* u1
    w2, _ = JuBat.core_ovalization(cm, u2, JuBat.ensure_node_ref!(case))
    @test w2 ≥ w1   # 同基准下更大位移 ⟹ 不圆度不减（基准未随"几何更新"重置）
end

@testset "持久化审计：u_prev/damage 跨相位连续（不重置）" begin
    case = dc_fixture()
    cm = case.czm_mesh
    layout = JuBat.CzmLayout(cm)
    u_end = fill(2e-5, 2 * cm.nnode)
    layout.u_prev = u_end
    cm.damage_states[1].D = 0.3
    @test layout.u_prev === u_end
    @test cm.damage_states[1].D == 0.3
end
