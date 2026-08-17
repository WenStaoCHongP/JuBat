using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ==========================================================================
# 宏观力学模块无量纲化重设计 v2 验证
# 见 docs/planning-with-files/力学模块修改/宏观力学模块无量纲化重设计.md §7
# ==========================================================================

@testset "Scale 锚点定义（重设计 v2 §2.1）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    scale = param_dim.scale
    coh = param_dim.cohesive

    # δ_czm = 2G_c/σ_max（断裂能定义的临界分离，替代旧 scale.L）
    @test scale.δ_czm ≈ 2 * coh.G_c_pe_pcc / coh.σ_max_pe_pcc rtol=1e-12
    @test scale.δ_czm ≉ scale.L   # 不再是叠层厚度

    # 联动定义
    @test scale.σ_czm ≈ coh.σ_max_pe_pcc rtol=1e-12
    @test scale.G_czm ≈ scale.σ_czm * scale.δ_czm rtol=1e-12
    @test scale.K_czm ≈ scale.σ_czm / scale.δ_czm rtol=1e-12
end

@testset "锚定界面归一化恒等式（重设计 v2 §4.2 / §7-4）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    coh = case.param.cohesive

    # σ_max* ≡ 1（自锚点）
    @test coh.σ_max_pe_pcc ≈ 1.0 rtol=1e-12
    # δ_c* ≡ 1（δ 锚点）
    @test coh.δ_c_pe_pcc ≈ 1.0 rtol=1e-12
    # G_c* ≡ 1/2（双线性三角形面积恒等式）
    @test coh.G_c_pe_pcc ≈ 0.5 rtol=1e-12
    # K_n*·δ_0* ≡ 1（刚度-起始分离互逆）
    @test coh.K_n_pe_pcc * coh.δ_0_pe_pcc ≈ 1.0 rtol=1e-9
    # 本构可分辨性下界（§6）：δ_0* ≤ 0.1
    @test coh.δ_0_pe_pcc ≤ 0.1
end

@testset "派生量 Λ / E* / L_ch（重设计 v2 §2.2）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    scale = param_dim.scale

    cache = JuBat.compute_czm_params_per_interface(case)
    pe = cache.by_interface[:PE_PCC]
    ne = cache.by_interface[:NE_NCC]

    # Λ = L/δ_czm，两界面共享
    @test pe.Λ ≈ scale.L / scale.δ_czm rtol=1e-12
    @test ne.Λ ≈ pe.Λ

    # E*（双材料调和平均，σ_czm 归一）
    E_star_pe_dim = 2 * param_dim.PE.E_coat * param_dim.PCC.E /
                    (param_dim.PE.E_coat + param_dim.PCC.E)
    @test pe.E_star ≈ E_star_pe_dim / scale.σ_czm rtol=1e-9

    # L_ch = E*·G_c/σ_max²（L 归一）
    L_ch_pe_dim = E_star_pe_dim * param_dim.cohesive.G_c_pe_pcc /
                  param_dim.cohesive.σ_max_pe_pcc^2
    @test pe.L_ch ≈ L_ch_pe_dim / scale.L rtol=1e-9

    # 位移空间有效刚度不变性（§5 一致性校验）：
    # K_n*·Λ == K_n_dim·L/σ_czm（与旧方案 δ_czm=L 的 K* 逐值相同）
    @test pe.K_n * pe.Λ ≈ param_dim.cohesive.K_n_pe_pcc * scale.L / scale.σ_czm rtol=1e-9
end

@testset "moduli_of 双重再缩放到 σ_czm 空间（重设计 v2 §3）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)

    cache = JuBat.compute_czm_params_per_interface(case)
    # 体刚度模量必须与 CzmInterfaceParams.E_eff 同参考（σ_czm 空间）
    E_pe, ν_pe = JuBat.moduli_of(case.param, :PE)
    @test E_pe ≈ cache.by_interface[:PE_PCC].E_eff rtol=1e-12
    @test ν_pe ≈ case.param.PE.nu_coat
    E_ne, _ = JuBat.moduli_of(case.param, :NE)
    @test E_ne ≈ cache.by_interface[:NE_NCC].E_eff rtol=1e-12
    # 连续层同样转到 σ_czm 空间：还原到物理值应等于原始参数
    E_sp, _ = JuBat.moduli_of(case.param, :SP)
    @test E_sp * param_dim.scale.σ_czm ≈ param_dim.SP.E rtol=1e-9
end

@testset "装配层 Λ 物理等价性（重设计 v2 §5）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = mesh_data.czm_submesh
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)
    czm_mesh = case.czm_mesh
    scale = param_dim.scale

    # 选第一个 cohesive 单元，沿其法向对顶面节点施加微小张开位移（弹性段内）
    elem = czm_mesh.cohesive_elements[1]
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top
    x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
    x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
    Lel = hypot(x2 - x1, y2 - y1)
    tx, ty = (x2 - x1) / Lel, (y2 - y1) / Lel
    nx, ny = -ty, tx

    pe = param_cache.by_interface[elem.interface_type]
    # 位移空间张开量（L 归一），映射到分离空间后远小于 δ_0*（保持弹性）
    δ_u = 0.01 * pe.δ_0_n / pe.Λ

    u = zeros(2 * czm_mesh.nnode)
    for n in (n3, n4)
        u[2n - 1] = δ_u * nx
        u[2n]     = δ_u * ny
    end

    K_coh, f_coh, seps, tracts = JuBat.assemble_czm_system(czm_mesh, u, param_cache)

    δ_n_tilde, δ_t_tilde = seps[1]
    T_n_tilde, _ = tracts[1]

    # (1) 分离换算：δ̃ = Λ·(B·ũ)
    @test δ_n_tilde ≈ pe.Λ * δ_u rtol=1e-9
    @test abs(δ_t_tilde) < 1e-12 * max(1.0, abs(δ_n_tilde))

    # (2) 弹性段本构：T̃ = K_n*·δ̃
    @test T_n_tilde ≈ pe.K_n * δ_n_tilde rtol=1e-9

    # (3) 物理等价：T_phys = T̃·σ_czm == K_n_dim·δ_phys
    δ_phys = δ_u * scale.L
    T_phys = T_n_tilde * scale.σ_czm
    @test T_phys ≈ param_dim.cohesive.K_n_pe_pcc * δ_phys rtol=1e-9
    # 确认确实处于弹性段
    @test δ_phys < param_dim.cohesive.δ_0_pe_pcc

    # (4) 切线刚度含一次 Λ：K_coh 在该单元法向 DOF 上的量级 = Λ·K_n*·O(几何权重)
    dof_n3x = 2 * n3 - 1
    @test isfinite(K_coh[dof_n3x, dof_n3x])
    @test !any(isnan, f_coh)
end

@testset "gap conductance 单位契约（重设计 v2，Λ 转换）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    scale = param_dim.scale

    cache = JuBat.compute_czm_params_per_interface(case)
    pe = cache.by_interface[:PE_PCC]
    coh_dim = param_dim.cohesive

    # 物理间隙 δ_phys 落入分支 2（δ_0 < δ < threshold）
    δ_phys = 1e-8   # [m]
    @test coh_dim.δ_0_pe_pcc < δ_phys < coh_dim.threshold
    D = 0.5

    # 归一化输入：分离空间（δ_czm 归一）
    δ_tilde = δ_phys / scale.δ_czm
    h_nd = JuBat.compute_gap_conductance(D, δ_tilde, pe)

    # 物理还原：h_phys = h_nd · λ_scale / L
    h_phys = h_nd * scale.lambda / scale.L
    two_beta_lambda = 2.0 * coh_dim.beta * coh_dim.lambda_m
    h_expected = coh_dim.h_c0 * (1 - D) + coh_dim.k_air / (δ_phys + two_beta_lambda)
    @test h_phys ≈ h_expected rtol=1e-6
end
