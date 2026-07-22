include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "unit_czm_newton.jl"))

"""
1D 串联（单位出平面厚度、宽度 W 约掉）：
T = -Σ(h ε0) / (Σ h/E + Σ 1/(K_n Λ))
δ̃_j = T / K_n,j   （分离空间，与 FEM seps 同参考）
冷却 ε0<0 → T>0 → 张开。
"""
function stack_1d_elastic_openings(
        heights::Vector{Float64},
        E::Vector{Float64},
        ε0::Vector{Float64},
        K_n::Vector{Float64},
        Λ::Vector{Float64})
    @assert length(heights) == length(E) == length(ε0) == 8
    @assert length(K_n) == length(Λ) == 4
    compliance = sum(heights[i] / E[i] for i in 1:8) +
                 sum(1.0 / (K_n[j] * Λ[j]) for j in 1:4)
    T = -sum(heights[i] * ε0[i] for i in 1:8) / compliance
    δ_tilde = [T / K_n[j] for j in 1:4]
    return T, δ_tilde
end

@testset "eigenstrain openings vs 1D analytic (elastic)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    param = case.param
    czm_mesh, meta = JuBat.create_unit_czm_strip(param; y0=1.0)
    cache = JuBat.compute_czm_params_per_interface(case)
    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:4]

    α_eff = cache.by_interface[:PE_PCC].α
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0

    ΔT_end = -0.02
    n_steps = 10
    u = zeros(2 * czm_mesh.nnode)

    # 固定底+顶全部 DOF，并约束所有节点 u_x=0，逼近 1D
    function bc_fixed_ends_1d(meta, czm_mesh)
        bc_dofs = Int[]; bc_vals = Float64[]
        fixed_xy = Set(vcat(meta.bottom_nodes, meta.top_nodes_after_czm))
        for n in 1:czm_mesh.nnode
            push!(bc_dofs, 2n - 1); push!(bc_vals, 0.0)  # u_x = 0
            if n in fixed_xy
                push!(bc_dofs, 2n); push!(bc_vals, 0.0)  # u_y = 0
            end
        end
        return bc_dofs, bc_vals
    end

    for s in 1:n_steps
        ΔT = ΔT_end * s / n_steps
        dT = fill(ΔT, 8)
        Δsoc_n = zeros(8)
        Δsoc_p = zeros(8)

        F_tc = JuBat.assemble_thermal_chemical_load(
            czm_mesh, cache, α_eff, β_n, β_p, dT, Δsoc_n, Δsoc_p)
        bc_dofs, bc_vals = bc_fixed_ends_1d(meta, czm_mesh)
        u, seps, tracts, ok, Rn = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals,
            F_thermo_chem=F_tc)
        @test ok

        ε0 = [α_eff * dT[e] + β_n * Δsoc_n[e] + β_p * Δsoc_p[e] for e in 1:8]
        Eν = [JuBat.moduli_of(param, meta.layer_materials[e]) for e in 1:8]
        # ux=0 + 热载荷因子 E(1+ν)/(1-ν²)=E/(1-ν) → 1D 用 E_eff=E/(1-ν)
        E = Float64[Ev / (1.0 - ν) for (Ev, ν) in Eν]
        K_n = Float64[
            cache.by_interface[czm_mesh.cohesive_elements[i].interface_type].K_n
            for i in 1:4]
        Λs = Float64[
            cache.by_interface[czm_mesh.cohesive_elements[i].interface_type].Λ
            for i in 1:4]
        T_ana, δ_ana = stack_1d_elastic_openings(meta.heights, E, ε0, K_n, Λs)
        @test T_ana > 0

        for i in 1:4
            δn_fem, _ = seps[i]
            Tn_fem, _ = tracts[i]
            δn_ana = δ_ana[i]
            @test abs(δn_fem) < cache.by_interface[:PE_PCC].δ_0_n
            @test czm_mesh.damage_states[i].D ≈ 0.0 atol=1e-12
            # FEM 法向约定可能与 1D 拉伸正号相反；比绝对值与 |T|
            @test abs(δn_fem) ≈ abs(δn_ana) rtol=0.15 atol=1e-10
            @test abs(Tn_fem) ≈ abs(T_ana) rtol=0.15 atol=1e-10
        end
    end
end

println("unit_czm_eigenstrain: ALL PASS")
