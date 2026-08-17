include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test
using LinearAlgebra
using Printf
include(joinpath(@__DIR__, "unit_czm_newton.jl"))

function _setup_strip()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    czm_mesh, meta = JuBat.create_unit_czm_strip(case.param; y0=1.0)
    cache = JuBat.compute_czm_params_per_interface(case)
    return case, czm_mesh, meta, cache
end

"""固定全部节点，仅释放并驱动指定 cohesive 的顶面（副本）节点。"""
function _bc_drive_cohesive(czm_mesh, coh_idx; u_n=0.0, u_t=0.0)
    coh = czm_mesh.cohesive_elements[coh_idx]
    n1, n2 = coh.nodes_bottom
    n4, n3 = coh.nodes_top   # nodes_top=[n_lo_c, n_hi_c] → n4,n3
    # 法向/切向（与 assemble 一致）
    x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
    x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
    L = hypot(x2 - x1, y2 - y1)
    tx, ty = (x2 - x1) / L, (y2 - y1) / L
    nx, ny = -ty, tx

    top_set = Set([n3, n4])
    bot_set = Set([n1, n2])
    bc_dofs = Int[]; bc_vals = Float64[]
    for n in 1:czm_mesh.nnode
        if n in top_set
            push!(bc_dofs, 2n - 1); push!(bc_vals, u_n * nx + u_t * tx)
            push!(bc_dofs, 2n);     push!(bc_vals, u_n * ny + u_t * ty)
        else
            # 所有非顶面节点固定（含 bottom 与其余 bulk）
            push!(bc_dofs, 2n - 1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2n);     push!(bc_vals, 0.0)
        end
    end
    return bc_dofs, bc_vals, nx, ny
end

@testset "Mode I monotonic bilinear" begin
    case, czm_mesh, meta, cache = _setup_strip()
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    n_steps = 40
    # 位移空间张开：扫过 δ_c
    u_max = 1.2 * pe.δ_c_n / Λ
    u = zeros(2 * czm_mesh.nnode)
    δ_hist = Float64[]; D_hist = Float64[]
    for s in 1:n_steps
        u_n = u_max * s / n_steps
        bc_dofs, bc_vals, _, _ = _bc_drive_cohesive(czm_mesh, i; u_n=u_n)
        u, seps, tracts, ok, Rn = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        δn, _ = seps[i]; Tn, _ = tracts[i]
        push!(δ_hist, δn)
        push!(D_hist, czm_mesh.damage_states[i].D)
        T_ana = analytic_bilinear_T(δn, pe.K_n, pe.δ_0_n, pe.δ_c_n)
        if δn < pe.δ_0_n * 0.98
            @test Tn ≈ T_ana rtol=1e-2 atol=1e-8
        elseif δn < pe.δ_c_n
            @test Tn ≈ T_ana rtol=5e-2 atol=1e-6
        else
            @test abs(Tn) < 1e-3 * max(pe.K_n * pe.δ_0_n, 1.0)
            @test czm_mesh.damage_states[i].D > 0.99
        end
    end
    δ_max = maximum(δ_hist)
    @printf("[Mode I mono] δ_n: min=%.6e  max=%.6e  (δ_0=%.6e  δ_c=%.6e)  D_max=%.4f\n",
            minimum(δ_hist), δ_max, pe.δ_0_n, pe.δ_c_n, maximum(D_hist))
    @test δ_max > pe.δ_0_n
    @test δ_max > 0.5 * pe.δ_c_n   # 必须进入软化/近断裂，避免近零误通过
    @test maximum(D_hist) > 0.99
end

@testset "Mode I unload reload" begin
    case, czm_mesh, meta, cache = _setup_strip()
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    u = zeros(2 * czm_mesh.nnode)
    u_peak = 2.0 * pe.δ_0_n / Λ   # D ≈ 0.5（δ_0 ≪ δ_c 时）
    local seps
    for s in 1:20
        u_n = u_peak * s / 20
        bc_dofs, bc_vals, _, _ = _bc_drive_cohesive(czm_mesh, i; u_n=u_n)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
    end
    D_peak = czm_mesh.damage_states[i].D
    δ_peak = seps[i][1]
    @printf("[Mode I unload] at peak: δ_n=%.6e  D=%.4f  (δ_0=%.6e)\n",
            δ_peak, D_peak, pe.δ_0_n)
    @test δ_peak > pe.δ_0_n          # 峰值必须越过损伤起始
    @test 0.05 < D_peak < 0.99
    for s in 1:20
        u_n = u_peak * (1 - s / 20)
        bc_dofs, bc_vals, _, _ = _bc_drive_cohesive(czm_mesh, i; u_n=u_n)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        @test czm_mesh.damage_states[i].D >= D_peak - 1e-12
    end
    D_unload = czm_mesh.damage_states[i].D
    @printf("[Mode I unload] after unload/reload: D=%.4f  (unchanged from peak)\n", D_unload)
    @test D_unload ≈ D_peak rtol=1e-8
    for s in 1:10
        u_n = 0.5 * u_peak * s / 10
        bc_dofs, bc_vals, _, _ = _bc_drive_cohesive(czm_mesh, i; u_n=u_n)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        @test czm_mesh.damage_states[i].D ≈ D_unload rtol=1e-8
    end
end

@testset "Mode II tangential" begin
    case, czm_mesh, meta, cache = _setup_strip()
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:czm_mesh.n_cohesive]
    i = findfirst(e -> e.interface_type == :PE_PCC, czm_mesh.cohesive_elements)
    pe = cache.by_interface[:PE_PCC]
    Λ = pe.Λ
    u = zeros(2 * czm_mesh.nnode)
    u_max = 1.2 * pe.δ_c_t / Λ
    δt_hist = Float64[]
    D_hist = Float64[]
    for s in 1:30
        u_t = u_max * s / 30
        bc_dofs, bc_vals, _, _ = _bc_drive_cohesive(czm_mesh, i; u_n=0.0, u_t=u_t)
        u, seps, tracts, ok, _ = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals)
        @test ok
        δn, δt = seps[i]; Tn, Tt = tracts[i]
        push!(δt_hist, abs(δt))
        push!(D_hist, czm_mesh.damage_states[i].D)
        @test abs(δn) < 0.2 * max(abs(δt), pe.δ_0_t)
        if abs(δt) < pe.δ_0_t * 0.98
            @test Tt ≈ pe.K_t * δt rtol=5e-2 atol=1e-8
        end
    end
    δt_max = maximum(δt_hist)
    @printf("[Mode II] |δ_t|: max=%.6e  (δ_0_t=%.6e  δ_c_t=%.6e)  D_max=%.4f\n",
            δt_max, pe.δ_0_t, pe.δ_c_t, maximum(D_hist))
    @test δt_max > pe.δ_0_t          # 必须越过切向损伤起始
end

println("unit_czm_bilinear: ALL PASS")
