using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ==========================================================================
# 宏观力学模块无量纲化重设计 v2 验证
# 见 docs/planning-with-files/14_力学模块修改/宏观力学模块无量纲化重设计.md §7
# ==========================================================================

@testset "Scale 锚点定义（重设计 v2 §2.1）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    scale = param_dim.scale

    # δ_czm = 2G_c/σ_max（断裂能定义的临界分离，替代旧 scale.L）；锚点 = PCC
    @test scale.δ_czm ≈ 2 * param_dim.PCC.G_c / param_dim.PCC.σ_max rtol=1e-12
    @test scale.δ_czm ≉ scale.L   # 不再是叠层厚度

    # 联动定义
    @test scale.σ_czm ≈ param_dim.PCC.σ_max rtol=1e-12
    @test scale.G_czm ≈ scale.σ_czm * scale.δ_czm rtol=1e-12
    @test scale.K_czm ≈ scale.σ_czm / scale.δ_czm rtol=1e-12
end

@testset "锚定界面归一化恒等式（重设计 v2 §4.2 / §7-4）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    pcc = case.param.PCC

    # σ_max* ≡ 1（自锚点）
    @test pcc.σ_max ≈ 1.0 rtol=1e-12
    # δ_c* ≡ 1（δ 锚点）
    @test pcc.δ_c ≈ 1.0 rtol=1e-12
    # G_c* ≡ 1/2（双线性三角形面积恒等式）
    @test pcc.G_c ≈ 0.5 rtol=1e-12
    # K_n*·δ_0* ≡ 1（刚度-起始分离互逆）
    @test pcc.K_n * pcc.δ_0 ≈ 1.0 rtol=1e-9
    # 本构可分辨性下界（§6）：δ_0* ≤ 0.1
    @test pcc.δ_0 ≤ 0.1
end

@testset "锚定量与 Λ 内联（重设计 v2 §2.2，2026-08-30 重构后）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    scale = param_dim.scale
    pcc = case.param.PCC

    # 锚定界面（PE-PCC = PCC）：σ_czm = PCC.σ_max_dim，δ_czm = 2·G_c/σ_max，δ_c* ≡ 1
    @test scale.σ_czm ≈ param_dim.PCC.σ_max rtol=1e-12
    @test scale.δ_czm ≈ 2 * param_dim.PCC.G_c / param_dim.PCC.σ_max rtol=1e-12
    @test pcc.σ_max ≈ 1.0 rtol=1e-12
    @test pcc.δ_c ≈ 1.0 rtol=1e-12

    # Λ 使用点内联（不再存字段）：位移空间有效刚度不变性（§5 一致性校验）
    Λ = scale.L / scale.δ_czm
    @test pcc.K_n * Λ ≈ param_dim.PCC.K_n * scale.L / scale.σ_czm rtol=1e-9
end

@testset "moduli_of 双重再缩放到 σ_czm 空间（重设计 v2 §3）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)

    # 体刚度模量与涂层模量同参考（σ_czm 空间）
    E_pe, ν_pe = JuBat.moduli_of(case.param, :PE)
    @test E_pe ≈ case.param.PE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm rtol=1e-12
    @test ν_pe ≈ case.param.PE.nu_coat
    E_ne, _ = JuBat.moduli_of(case.param, :NE)
    @test E_ne ≈ case.param.NE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm rtol=1e-12
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
    case.mech = JuBat.MechState(case.czm_mesh)
    czm_mesh = case.czm_mesh
    scale = param_dim.scale

    # 选第一个 cohesive 单元，沿其法向对顶面节点施加微小张开位移（弹性段内）
    elem = czm_mesh.cohesive_elements[1]
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top
    nx, ny = JuBat.cohesive_geometry(czm_mesh)[1].n_vec

    ip = JuBat.collector_params(case.param, elem.interface_type)
    Λ = scale.L / scale.δ_czm
    # 位移空间张开量（L 归一），映射到分离空间后远小于 δ_0*（保持弹性）
    δ_u = 0.01 * ip.δ_0 / Λ

    u = zeros(2 * czm_mesh.nnode)
    for n in (n3, n4)
        u[2n - 1] = δ_u * nx
        u[2n]     = δ_u * ny
    end

    K_coh, f_coh, seps, tracts = JuBat.assemble_czm_system(
        czm_mesh, u, case.param; damage_states=case.mech.damage_states,
        geom_cache=JuBat.cohesive_geometry(czm_mesh),
        ws=JuBat.assembly_workspace(czm_mesh))

    δ_n_tilde, δ_t_tilde = seps[1]
    T_n_tilde, _ = tracts[1]

    # (1) 分离换算：δ̃ = Λ·(B·ũ)
    @test δ_n_tilde ≈ Λ * δ_u rtol=1e-9
    @test abs(δ_t_tilde) < 1e-12 * max(1.0, abs(δ_n_tilde))

    # (2) 弹性段本构：T̃ = K_n*·δ̃
    @test T_n_tilde ≈ ip.K_n * δ_n_tilde rtol=1e-9

    # (3) 物理等价：T_phys = T̃·σ_czm == K_n_dim·δ_phys
    δ_phys = δ_u * scale.L
    T_phys = T_n_tilde * scale.σ_czm
    @test T_phys ≈ param_dim.PCC.K_n * δ_phys rtol=1e-9
    # 确认确实处于弹性段
    @test δ_phys < param_dim.PCC.δ_0

    # (4) 切线刚度含一次 Λ：K_coh 在该单元法向 DOF 上的量级 = Λ·K_n*·O(几何权重)
    dof_n3x = 2 * n3 - 1
    @test isfinite(K_coh[dof_n3x, dof_n3x])
    @test !any(isnan, f_coh)
end

@testset "gap conductance 单位契约（重设计 v2，Λ 内联转换）" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    scale = param_dim.scale
    ip = case.param.PCC

    # 物理间隙 δ_phys 落入分支 2（δ_0 < δ < threshold）
    δ_phys = 1e-8   # [m]
    @test param_dim.PCC.δ_0 < δ_phys < param_dim.PCC.threshold
    D = 0.5

    # 归一化输入：分离空间（δ_czm 归一）
    δ_tilde = δ_phys / scale.δ_czm
    h_nd = JuBat.compute_gap_conductance(D, δ_tilde, ip, case.param)

    # 物理还原：h_phys = h_nd · λ_scale / L
    h_phys = h_nd * scale.lambda / scale.L
    two_beta_lambda = 2.0 * param_dim.PCC.beta * param_dim.PCC.lambda_m
    h_expected = param_dim.PCC.h_c0 * (1 - D) + param_dim.PCC.k_air / (δ_phys + two_beta_lambda)
    @test h_phys ≈ h_expected rtol=1e-6
end
