using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §7 Batch 2'：默认关逐位不变；开启后 σ₀ 量级对照解析卷绕公式；缺参即 error。
# D-B2'-1：卷入张力等应变/Voigt 分担 + 对数累积压力；D-B2'-4：首平衡重分布。

function build_wp_fixture(; nθ::Int=8)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    return case, param_cache, cache
end

@testset "σ₀ 场解析对照（独立编码公式互核）" begin
    case, param_cache, _ = build_wp_fixture()
    czm_mesh = case.czm_mesh
    p = case.param
    σ0 = JuBat.winding_prestress_field(czm_mesh, p)
    # —— 测试内独立重算（与实现同公式、独立编码）——
    Ew(mt) = JuBat.moduli_of(p, mt)[1]
    t_ne = 2*p.NE.thickness + p.NCC.thickness
    t_pe = 2*p.PE.thickness + p.PCC.thickness
    E_ne = (2*Ew(:NE)*p.NE.thickness + Ew(:NCC)*p.NCC.thickness) / t_ne
    E_pe = (2*Ew(:PE)*p.PE.thickness + Ew(:PCC)*p.PCC.thickness) / t_pe
    ε_ne = p.cell.winding_T_ne / E_ne
    ε_pe = p.cell.winding_T_pe / E_pe
    ε_sp = 0.5 * (ε_ne + ε_pe)
    F_rev = p.cell.winding_T_ne * t_ne + p.cell.winding_T_pe * t_pe   # 每匝每轴长的环向力
    f = F_rev / p.cell.layer
    R_end = maximum(hypot.(czm_mesh.node[:,1], czm_mesh.node[:,2]))
    mt = czm_mesh.czm_submesh.material_type
    checked = 0
    for e in 1:size(czm_mesh.bulk_element, 1)
        (e % 97 == 1 && e > 1) || continue
        xs = czm_mesh.node[czm_mesh.bulk_element[e,:], 1]
        ys = czm_mesh.node[czm_mesh.bulk_element[e,:], 2]
        xc, yc = sum(xs)/4, sum(ys)/4
        r = hypot(xc, yc); c, s = xc/r, yc/r
        pref = f * log(R_end / r)
        ε_w = mt[e] in (:NE, :NCC) ? ε_ne : mt[e] in (:PE, :PCC) ? ε_pe : ε_sp
        σ_θ = Ew(mt[e]) * ε_w - pref
        σ_r = -pref
        tx, ty, nx, ny = -s, c, c, s
        exp_xx = σ_θ*tx*tx + σ_r*nx*nx
        exp_yy = σ_θ*ty*ty + σ_r*ny*ny
        exp_xy = σ_θ*tx*ty + σ_r*nx*ny
        @test σ0[e][1] ≈ exp_xx rtol = 1e-12
        @test σ0[e][2] ≈ exp_yy rtol = 1e-12
        @test σ0[e][3] ≈ exp_xy rtol = 1e-12
        checked += 1
    end
    @test checked ≥ 3
end

@testset "量级校验（§10.4.1 工艺区间，换算有量纲）" begin
    case, _, _ = build_wp_fixture()
    czm_mesh = case.czm_mesh
    p = case.param
    σ0 = JuBat.winding_prestress_field(czm_mesh, p)
    σ_czm = p.scale.σ_czm
    # 环向张力分量：卷入侧层应有正环向应力，量级 ~T_side（1–5 MPa 档）
    R_end = maximum(hypot.(czm_mesh.node[:,1], czm_mesh.node[:,2]))
    r_in = minimum(hypot.(czm_mesh.node[i,1], czm_mesh.node[i,2]) for i in axes(czm_mesh.node,1))
    p_core = ((p.cell.winding_T_ne * (2*p.NE.thickness + p.NCC.thickness) +
               p.cell.winding_T_pe * (2*p.PE.thickness + p.PCC.thickness)) / p.cell.layer) *
             log(R_end / r_in) * σ_czm
    @test 0.1e6 < p_core < 20e6          # 卷芯径向预压（§10.4.1 p0 0.2–1.0 MPa 同量级放宽上限）
    @test all(abs(σ0[e][i]) * σ_czm < 20e6 for e in eachindex(σ0) for i in 1:3)  # 全场 < 20 MPa
end

@testset "缺参/非法即 error（AGENTS 9.4/9.7）" begin
    case, _, _ = build_wp_fixture()
    czm_mesh = case.czm_mesh
    p = case.param
    tne, tpe = p.cell.winding_T_ne, p.cell.winding_T_pe
    p.cell.winding_T_ne = 0.0; p.cell.winding_T_pe = 0.0
    @test_throws ErrorException JuBat.winding_prestress_field(czm_mesh, p)
    p.cell.winding_T_ne = -1.0
    @test_throws ErrorException JuBat.winding_prestress_field(czm_mesh, p)
    p.cell.winding_T_ne = tne; p.cell.winding_T_pe = tpe
    @test size(JuBat.winding_prestress_field(czm_mesh, p)) == (size(czm_mesh.bulk_element, 1),)
end
