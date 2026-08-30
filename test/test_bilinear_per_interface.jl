using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ---------------------------------------------------------------------------
# 共享参数夹具：PE-PCC 接面（与原始测试一致）
# ---------------------------------------------------------------------------
const PE_PCC_PARAMS = JuBat.CurrentCollector(
    σ_max = 82e6 / 1e10,
    K_n = 2.4e17 / (1e10 / 1e-6),
    δ_0 = 82e6 / 2.4e17,
    δ_c = 2 * 25.3 / 82e6,
    G_c = 25.3,
    τ_max = 82e6 / 1e10,
    K_t = 2.4e17 / (1e10 / 1e-6),
    δ_0_t = 82e6 / 2.4e17,
    δ_c_t = 2 * 25.3 / 82e6,
    G_c_t = 25.3,
    eta = 1.45,
    h_c0 = 1e7, k_air = 0.026, lambda_m = 70e-9,
    beta = 1.0, threshold = 70e-9,
)

# NE-NCC 接面参数（不同 K_n）——用于交叉验证
const NE_NCC_PARAMS = JuBat.CurrentCollector(;
    σ_max = 2 * PE_PCC_PARAMS.σ_max,
    K_n = 2 * PE_PCC_PARAMS.K_n,
    δ_0 = PE_PCC_PARAMS.δ_0,
    δ_c = PE_PCC_PARAMS.δ_c / 2,
    G_c = PE_PCC_PARAMS.G_c,
    τ_max = PE_PCC_PARAMS.τ_max,
    K_t = PE_PCC_PARAMS.K_t,
    δ_0_t = PE_PCC_PARAMS.δ_0_t,
    δ_c_t = PE_PCC_PARAMS.δ_c_t,
    G_c_t = PE_PCC_PARAMS.G_c_t,
    eta = PE_PCC_PARAMS.eta,
    h_c0 = PE_PCC_PARAMS.h_c0, k_air = PE_PCC_PARAMS.k_air, lambda_m = PE_PCC_PARAMS.lambda_m,
    beta = PE_PCC_PARAMS.beta, threshold = PE_PCC_PARAMS.threshold,
)

# 混合模式专用参数：Mode II 强度不同于 Mode I，以真正激活 BK eta 指数项
const MIX_MODE_PARAMS = JuBat.CurrentCollector(;
    σ_max = PE_PCC_PARAMS.σ_max,
    K_n = PE_PCC_PARAMS.K_n,
    δ_0 = PE_PCC_PARAMS.δ_0,
    δ_c = PE_PCC_PARAMS.δ_c,
    G_c = PE_PCC_PARAMS.G_c,
    τ_max = PE_PCC_PARAMS.τ_max,
    K_t = PE_PCC_PARAMS.K_t,
    # δ_0_t / δ_c_t 与 Mode I 不同 → β^eta 项不为零
    δ_0_t = 3 * PE_PCC_PARAMS.δ_0,
    δ_c_t = 3 * PE_PCC_PARAMS.δ_c,
    G_c_t = PE_PCC_PARAMS.G_c_t,
    eta = PE_PCC_PARAMS.eta,
    h_c0 = PE_PCC_PARAMS.h_c0, k_air = PE_PCC_PARAMS.k_air, lambda_m = PE_PCC_PARAMS.lambda_m,
    beta = PE_PCC_PARAMS.beta, threshold = PE_PCC_PARAMS.threshold,
)

@testset "bilinear_* with CzmInterfaceParams" begin
    params = PE_PCC_PARAMS
    D = JuBat.DamageState()

    # 弹性段：δ_n < δ_0_n
    δ_n_small = params.δ_0 / 2
    T_n, T_t, _, D_new = JuBat.bilinear_traction_state(δ_n_small, 0.0, D, params, "model1")
    @test T_n ≈ params.K_n * δ_n_small
    @test T_t ≈ 0.0
    @test D_new.D ≈ 0.0

    # 软化段：δ_0_n < δ_n < δ_c_n
    δ_n_mid = 0.5 * (params.δ_0 + params.δ_c)
    T_n2, _, _, D_new2 = JuBat.bilinear_traction_state(δ_n_mid, 0.0, D, params, "model1")
    @test 0 < T_n2 < params.σ_max
    @test 0 < D_new2.D < 1

    # 完全失效：δ_n > δ_c_n
    δ_n_big = 2 * params.δ_c
    T_n3, _, _, D_new3 = JuBat.bilinear_traction_state(δ_n_big, 0.0, D, params, "model1")
    @test T_n3 ≈ 0.0
    @test D_new3.D ≈ 1.0

    # 切向（Mode II）
    _, T_t2, _, _ = JuBat.bilinear_traction_state(0.0, params.δ_0_t / 2, D, params, "model1")
    @test T_t2 ≈ params.K_t * params.δ_0_t / 2

    # NE-NCC 接面参数（不同 K_n）应给出不同结果
    params_ne = NE_NCC_PARAMS
    T_n_ne, _, _, _ = JuBat.bilinear_traction_state(params_ne.δ_0 / 2, 0.0, D, params_ne, "model1")
    @test T_n_ne ≈ params_ne.K_n * params_ne.δ_0 / 2
    @test T_n_ne ≠ T_n
end

@testset "softening branch closed-form: T_n = (1-D)·K_n·δ_n" begin
    # Issue 2：验证软化段封闭公式（visc_beta=1.0 → D_visc = D_eq）
    params = PE_PCC_PARAMS
    D = JuBat.DamageState()
    δ_n_mid = 0.5 * (params.δ_0 + params.δ_c)
    T_n2, _, _, D_new2 = JuBat.bilinear_traction_state(δ_n_mid, 0.0, D, params, "model1")

    # 软化段封闭公式：T_n = (1 - D) * K_n * δ_n
    expected_T_n2 = (1 - D_new2.D) * params.K_n * δ_n_mid
    @test T_n2 ≈ expected_T_n2 rtol=1e-6
end

@testset "mixed-mode BK eta exponent (model2)" begin
    # Issue 3：覆盖 Materialmatrix.jl:110-111 的 β^eta 混合模式逻辑
    params = MIX_MODE_PARAMS
    D = JuBat.DamageState()

    # 选取同时超过 Mode I 与 Mode II 弹性极限的分离，落入软化段
    δ_n_mix = 5 * params.δ_0
    δ_t_mix = 5 * params.δ_0  # 同一物理尺度，混合

    T_n_mix, T_t_mix, _, D_mix = JuBat.bilinear_traction_state(δ_n_mix, δ_t_mix, D, params, "model2")

    # 混合模式有效分离：δ_eff = sqrt(δ_n^2 + δ_t^2)
    δ_eff = hypot(δ_n_mix, δ_t_mix)

    # 由于 δ_0_t ≠ δ_0_n，β^eta 项实际生效：δ_0_eff = sqrt(δ_0_n^2 + (δ_0_t^2 - δ_0_n^2) * β^eta)
    β = abs(δ_t_mix) / δ_eff
    δ_0_eff_expected = sqrt(params.δ_0^2 + (params.δ_0_t^2 - params.δ_0^2) * β^params.eta)
    @test δ_eff > δ_0_eff_expected  # 确认进入软化段

    # 弹性段 D=0，软化段 D>0
    @test D_mix.D > 0.0
    @test D_mix.D < 1.0

    # T_n 与 T_t 都应非负（张力/剪力不为压缩）
    @test T_n_mix >= 0
    @test T_t_mix >= 0

    # 交叉验证：δ_n=0 的纯 Mode II 也应触发 D>0
    _, _, _, D_pure_ii = JuBat.bilinear_traction_state(0.0, 5 * params.δ_0_t, D, params, "model2")
    @test D_pure_ii.D > 0.0
end

# ---------------------------------------------------------------------------
# Issue 1：以下 5 个函数此前零覆盖，逐一补齐烟雾测试
# ---------------------------------------------------------------------------

@testset "bilinear_traction wrapper" begin
    params = PE_PCC_PARAMS

    # update=true：应就地更新 damage_state
    D_mut = JuBat.DamageState()
    result = JuBat.bilinear_traction(params.δ_c * 2, 0.0, D_mut, params, "model1"; update=true)
    @test result isa Tuple
    @test length(result) == 3
    T_n, T_t, D_val = result
    @test T_n ≈ 0.0           # 完全失效
    @test D_val ≈ 1.0
    @test D_mut.D ≈ 1.0       # 被就地更新
    @test D_mut.fractured == true

    # update=false：不应修改 damage_state
    D_keep = JuBat.DamageState()
    result2 = JuBat.bilinear_traction(params.δ_0 / 2, 0.0, D_keep, params, "model1"; update=false)
    @test result2 isa Tuple
    @test D_keep.D ≈ 0.0      # 未被修改
    @test D_keep.fractured == false

    # 与 bilinear_traction_state 数值一致性
    D_ref = JuBat.DamageState()
    T_n_s, T_t_s, D_eq_s, _ = JuBat.bilinear_traction_state(
        params.δ_0 / 2, 0.0, D_ref, params, "model1")
    T_n_w, T_t_w, D_eq_w = JuBat.bilinear_traction(
        params.δ_0 / 2, 0.0, JuBat.DamageState(), params, "model1"; update=false)
    @test T_n_w ≈ T_n_s
    @test T_t_w ≈ T_t_s
    @test D_eq_w ≈ D_eq_s
end

@testset "bilinear_tangent returns 2×2 matrix" begin
    params = PE_PCC_PARAMS

    # 弹性段：切线刚度对角元 ≈ K_n（model1）
    D = JuBat.DamageState()
    K = JuBat.bilinear_tangent(params.δ_0 / 2, 0.0, D, params, "model1")
    @test size(K) == (2, 2)
    @test all(isfinite, K)
    @test K[1, 1] ≈ params.K_n      # 弹性段法向刚度
    @test K[2, 2] ≈ params.K_t      # model1 切向不损伤
    @test K[1, 2] ≈ 0.0
    @test K[2, 1] ≈ 0.0

    # 完全失效：切线刚度退化到接近零
    D_failed = JuBat.DamageState()
    D_failed.fractured = true
    K_f = JuBat.bilinear_tangent(params.δ_c * 2, 0.0, D_failed, params, "model1")
    @test size(K_f) == (2, 2)
    @test all(isfinite, K_f)
    @test K_f[1, 1] < params.K_n    # 法向已退化
    @test K_f[1, 1] > 0.0           # 仍有小刚度避免奇异

    # 混合模式（model2）：非对偶耦合，off-diagonal 应非零
    K_mix = JuBat.bilinear_tangent(
        5 * MIX_MODE_PARAMS.δ_0, 5 * MIX_MODE_PARAMS.δ_0,
        JuBat.DamageState(), MIX_MODE_PARAMS, "model2")
    @test size(K_mix) == (2, 2)
    @test all(isfinite, K_mix)
    # 软化段非对角项一般非零
    @test abs(K_mix[1, 2]) > 0 || abs(K_mix[2, 1]) > 0
end

@testset "update_damage batch" begin
    params = PE_PCC_PARAMS

    damages = [JuBat.DamageState() for _ in 1:3]
    seps = [
        (params.δ_0 / 2, 0.0),         # 弹性
        (params.δ_c * 2, 0.0),         # 完全失效
        (params.δ_0, params.δ_0_t),    # 边界
    ]

    new_states = JuBat.update_damage(damages, seps, params, "model1")
    @test length(new_states) == 3
    @test all(s isa JuBat.DamageState for s in new_states)

    @test new_states[1].D ≈ 0.0                    # 弹性段无损伤
    @test new_states[2].D ≈ 1.0                    # 完全失效
    @test new_states[2].fractured == true
    @test 0 ≤ new_states[3].D ≤ 1.0

    # 原输入不应被修改（update_damage 返回新状态）
    @test all(d.D ≈ 0.0 for d in damages)

    # 长度不匹配应抛出 AssertionError
    @test_throws AssertionError JuBat.update_damage(
        damages, [(0.0, 0.0)], params, "model1")
end

@testset "compute_gap_conductance 3-branch piecewise" begin
    # 2026-08-30 重构适配：compute_gap_conductance(D, δ_n, ip, param)——
    # δ_n 在 δ_czm 归一空间，入口 ×inv_Λ 转到 L 空间后走三分支。
    param_dim = JuBat.ChooseCell("Jellyroll")
    case = JuBat.SetCase(param_dim, JuBat.Option())
    param = case.param
    ip = param.PCC
    inv_Λ = param.scale.δ_czm / param.scale.L
    delta0 = ip.δ_0 * inv_Λ          # L 空间损伤起始间隙
    threshold = ip.threshold          # 本就在 L 空间
    two_beta_lambda = 2 * ip.beta * ip.lambda_m

    # 分支 1：delta < delta0 → h = h_c0 + k_air/(delta + 2βλ_m)，与 D 无关
    D_arbitrary = 0.7
    δ_b1 = 0.0
    k1 = JuBat.compute_gap_conductance(D_arbitrary, δ_b1, ip, param)
    expected1 = ip.h_c0 + ip.k_air / (δ_b1 * inv_Λ + two_beta_lambda)
    @test k1 ≈ expected1 rtol=1e-9
    @test k1 > 0
    # 与 D 无关
    @test JuBat.compute_gap_conductance(0.0, δ_b1, ip, param) ≈ k1

    # 分支 2：delta0 ≤ delta < threshold → h = h_c0*(1-D_clamped) + k_air/(delta + 2βλ_m)
    D_b2 = 0.5
    δ_b2 = (delta0 + threshold) / 2 / inv_Λ   # 转回 δ_czm 空间输入
    delta_b2_L = δ_b2 * inv_Λ
    @test delta0 < delta_b2_L < threshold
    k2 = JuBat.compute_gap_conductance(D_b2, δ_b2, ip, param)
    D_clamped_b2 = clamp(D_b2, 0.0, 0.9999)
    expected2 = ip.h_c0 * (1 - D_clamped_b2) + ip.k_air / (delta_b2_L + two_beta_lambda)
    @test k2 ≈ expected2 rtol=1e-9
    @test k2 > 0
    # D 增大 → h 减小（损伤降低接触导热）
    @test JuBat.compute_gap_conductance(0.9, δ_b2, ip, param) < k2

    # 分支 3：delta ≥ threshold → h = h_c0*(1-D_clamped) + k_air/(delta + delta0)
    D_b3 = 1.0
    δ_b3 = threshold * 10 / inv_Λ
    delta_b3_L = δ_b3 * inv_Λ
    @test delta_b3_L >= threshold
    k3 = JuBat.compute_gap_conductance(D_b3, δ_b3, ip, param)
    D_clamped_b3 = clamp(D_b3, 0.0, 0.9999)
    expected3 = ip.h_c0 * (1 - D_clamped_b3) + ip.k_air / (delta_b3_L + delta0)
    @test k3 ≈ expected3 rtol=1e-9
    @test k3 > 0

    # 边界处行为：分支 1↔2 / 2↔3 处均有限
    continuity_D = 0.5
    k_below = JuBat.compute_gap_conductance(continuity_D, delta0 * 0.999 / inv_Λ, ip, param)
    k_above = JuBat.compute_gap_conductance(continuity_D, delta0 * 1.001 / inv_Λ, ip, param)
    @test isfinite(k_below) && isfinite(k_above)

    k_thr_below = JuBat.compute_gap_conductance(continuity_D, threshold * 0.999 / inv_Λ, ip, param)
    k_thr_above = JuBat.compute_gap_conductance(continuity_D, threshold * 1.001 / inv_Λ, ip, param)
    @test isfinite(k_thr_below) && isfinite(k_thr_above)
    # 显式验证：两分支公式分母不同 → 不连续（除非 2βλ_m == delta0）
    if two_beta_lambda ≠ delta0
        @test k_thr_below ≠ k_thr_above
    end
end
@testset "compute_element_gap_conductance / compute_all_gap_conductances" begin
    # 2026-08-30 重构适配：批量接口接收 damage_states 向量（不再挂网格）
    param_dim = JuBat.ChooseCell("Jellyroll")
    case = JuBat.SetCase(param_dim, JuBat.Option())
    param = case.param
    ip = param.PCC

    # 手动填充 3 个损伤状态（对应 3 个 cohesive 单元）
    states = [
        JuBat.DamageState(),                            # D=0, δ_max_n=0
        JuBat.DamageState(),                            # D=0, δ_max_n=0
        JuBat.DamageState(),                            # D=0, δ_max_n=0
    ]
    # 第二个单元：人为设置已损伤状态（L 空间落入分支 2）
    states[2].D = 0.5
    states[2].δ_max_n = ((ip.δ_0 + ip.threshold / (param.scale.δ_czm / param.scale.L)) / 2)
    # 第三个单元：完全失效 + 大 gap（分支 3）
    states[3].D = 1.0
    states[3].δ_max_n = ip.threshold * 10 * param.scale.L / param.scale.δ_czm

    # 单元级
    h1 = JuBat.compute_element_gap_conductance(states, 1, ip, param)
    h2 = JuBat.compute_element_gap_conductance(states, 2, ip, param)
    h3 = JuBat.compute_element_gap_conductance(states, 3, ip, param)

    @test all(isfinite, (h1, h2, h3))
    @test all((h1, h2, h3) .> 0)

    # 与纯函数 compute_gap_conductance 一致
    @test h1 ≈ JuBat.compute_gap_conductance(states[1].D, states[1].δ_max_n, ip, param)
    @test h2 ≈ JuBat.compute_gap_conductance(states[2].D, states[2].δ_max_n, ip, param)
    @test h3 ≈ JuBat.compute_gap_conductance(states[3].D, states[3].δ_max_n, ip, param)

    # 全网格批量：长度匹配、值一致
    h_all = JuBat.compute_all_gap_conductances(states, ip, param)
    @test length(h_all) == 3
    @test h_all ≈ [h1, h2, h3]

    # D=0 的单元导热应大于 D=1 的单元（接触路径未被破坏）
    @test h1 > h3
end
