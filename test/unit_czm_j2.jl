using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec §7 Batch 3：单轴屈服/卸载、KKT、耗散非负、一致切线 FD（平面应力一致返回映射）。
# 本构级测试：直接用物理量纲构造（E=70 GPa Al 箔），不经归一化链。

const E_f = 70e9
const ν_f = 0.33
const σ_y = 60e6
const H_f = 0.0
const C_f = JuBat.plane_stress_C(E_f, ν_f)

@testset "弹性步不变（f≤0 时 σ=Ce，C_ep=C）" begin
    e = 2e-4 .* [1.0, -ν_f, 0.0]   # σ11 = E·2e-4 = 14 MPa < 60 MPa
    σ, C_ep, Δeps_p, Δκ = JuBat.return_mapping_plane_stress(
        e, C_f, σ_y, H_f, (0.0, 0.0, 0.0), 0.0)
    @test σ ≈ C_f * e
    @test C_ep ≈ C_f
    @test Δeps_p == (0.0, 0.0, 0.0) && Δκ == 0.0
end

@testset "屈服后停留屈服面（f=0 精确）+ 耗散非负 + Δγ>0" begin
    # 单轴应力路径：e = [e1, -ν e1, 0] 弹性段保持 σ22=0；加载至 σ11 超屈服
    e1 = 2.0 * σ_y / E_f           # 名义两倍屈服应变
    e = [e1, -ν_f * e1, 0.0]
    σ, C_ep, Δeps_p, Δκ = JuBat.return_mapping_plane_stress(
        e, C_f, σ_y, H_f, (0.0, 0.0, 0.0), 0.0)
    q = JuBat.qbar(σ)
    @test q ≈ σ_y + H_f * (0.0 + Δκ) rtol = 1e-10
    @test σ[1] > 0 && Δκ > 0
    # 耗散非负：σ:Δεᵖ = Δγ·σ̄ ≥ 0
    @test σ[1]*Δeps_p[1] + σ[2]*Δeps_p[2] + σ[3]*Δeps_p[3] ≥ 0
    # 塑性松弛：σ11 介于 0 与试算值 2σ_y 之间（平面应力关联流动 n=[1,−0.5,0]，
    # 返回态非单轴——σ22≠0，故不断言 σ11=σ_y，那仅 ν=0.5 时成立）
    @test 0.0 < σ[1] < 2.0 * σ_y
end

@testset "卸载弹性返回 + 永久应变 + 再加载折返" begin
    e1 = 2.0 * σ_y / E_f
    e_load = [e1, -ν_f * e1, 0.0]
    _, _, Δeps_p, Δκ = JuBat.return_mapping_plane_stress(
        e_load, C_f, σ_y, H_f, (0.0, 0.0, 0.0), 0.0)
    p = (Δeps_p[1], Δeps_p[2], Δeps_p[3])
    # 卸载至一半应变：弹性响应，应力低于屈服面
    e_un = 0.5 .* e_load
    σ, C_ep, Δ2, Δκ2 = JuBat.return_mapping_plane_stress(
        e_un, C_f, σ_y, H_f, p, Δκ)
    @test Δ2 == (0.0, 0.0, 0.0) && Δκ2 == 0.0          # 弹性步
    @test JuBat.qbar(σ) < σ_y
    # 卸载弹性响应：σ = C(e_un − p)（精确，含 c12 耦合）
    @test σ ≈ C_f * (e_un .- collect(p)) rtol = 1e-12
    # 再加载回 e_load：立即回屈服面（κ 不再增长的理想塑性再屈服）
    σ3, _, Δ3, Δκ3 = JuBat.return_mapping_plane_stress(
        e_load, C_f, σ_y, H_f, p, Δκ)
    @test JuBat.qbar(σ3) ≈ σ_y rtol = 1e-10
    @test Δκ3 ≥ 0
end

@testset "一致切线有限差分（塑性流动区）" begin
    # 流动状态点：先加载到屈服并提交塑性状态，再在流动点做 FD
    e1 = 2.0 * σ_y / E_f
    e_load = [e1, -ν_f * e1, 0.0]
    _, _, Δeps_p, Δκ = JuBat.return_mapping_plane_stress(
        e_load, C_f, σ_y, H_f, (0.0, 0.0, 0.0), 0.0)
    p = (Δeps_p[1], Δeps_p[2], Δeps_p[3])
    e_flow = [e1 * 1.05, -ν_f * e1 * 1.02, 3e-4]   # 含剪切分量的流动点
    σ0, C_ep, _, _ = JuBat.return_mapping_plane_stress(
        e_flow, C_f, σ_y, H_f, p, Δκ)
    h = 1e-8
    for j in 1:3
        ep = copy(e_flow); ep[j] += h
        em = copy(e_flow); em[j] -= h
        σp, _, _, _ = JuBat.return_mapping_plane_stress(ep, C_f, σ_y, H_f, p, Δκ)
        σm, _, _, _ = JuBat.return_mapping_plane_stress(em, C_f, σ_y, H_f, p, Δκ)
        fd = (σp .- σm) ./ (2h)
        @test isapprox(C_ep[:, j], fd; rtol = 1e-5, atol = 1e3)
    end
end

@testset "双轴+剪切组合步应力点回面（H>0 变体）" begin
    E2, ν2, σy2, H2 = 110e9, 0.34, 200e6, 1e9   # Cu 箔带硬化
    C2 = JuBat.plane_stress_C(E2, ν2)
    p0 = (1e-3, 1e-3, 5e-4)                      # 已有塑性
    κ0 = 0.02
    e = 3e-2 .* [1.0, 0.6, 0.35]
    σ, C_ep, Δp, Δκ = JuBat.return_mapping_plane_stress(e, C2, σy2, H2, p0, κ0)
    @test JuBat.qbar(σ) ≈ σy2 + H2 * (κ0 + Δκ) rtol = 1e-10
    @test all(isfinite, C_ep)
end
