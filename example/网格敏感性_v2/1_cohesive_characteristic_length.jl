"""
Script 1: 计算内聚力特征长度 l_c 及确定 nθ 范围

输出：
  - E_eff (体材料等效杨氏模量)
  - l_c  (cohesive 特征长度)
  - 内圈/外圈 nθ 约束
  - 4 个推荐的 CZM nθ 水平
  - 热网格 nθ 候选值 (经验)
"""

using Printf

root_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(root_dir, "src", "JuBat.jl"))
using .JuBat

function main()
    println("=" ^ 70)
    println("Script 1: 内聚力特征长度计算")
    println("=" ^ 70)

    param_dim = JuBat.ChooseCell("Jellyroll")

    # ── 提取参数 ──
    E_NE  = param_dim.NE.E              # Pa
    E_PE  = param_dim.PE.E              # Pa
    t_NE  = param_dim.NE.thickness      # m
    t_PE  = param_dim.PE.thickness      # m
    # G_c   = param_dim.cohesive.G_c_n    # J/m²                              # TODO Chunk 2 Task 2.1
    # σ_max = param_dim.cohesive.σ_max_n  # Pa                                # TODO Chunk 2 Task 2.1
    R_in  = param_dim.cell.Rin          # m
    R_out = param_dim.cell.Rout         # m
    nu_NE = param_dim.NE.nu
    nu_PE = param_dim.PE.nu

    # ── E_eff (复现 Mechanical.jl:181) ──
    E_eff  = (E_NE * t_NE + E_PE * t_PE) / (t_NE + t_PE)
    nu_eff = (nu_NE * t_NE + nu_PE * t_PE) / (t_NE + t_PE)

    println("\n[1] 材料参数")
    @printf("  NE.E = %.2f GPa,  NE.thickness = %.1f μm\n", E_NE*1e-9, t_NE*1e6)
    @printf("  PE.E = %.2f GPa,  PE.thickness = %.1f μm\n", E_PE*1e-9, t_PE*1e6)
    @printf("  E_eff = %.2f GPa  (厚度加权平均)\n", E_eff*1e-9)
    @printf("  ν_eff = %.3f\n", nu_eff)

    # ── l_c ──
    l_c = G_c * E_eff / σ_max^2

    println("\n[2] 内聚力参数")
    @printf("  G_c   = %.2f J/m²\n", G_c)
    @printf("  σ_max = %.2f MPa\n", σ_max*1e-6)
    # @printf("  K_n   = %.2e Pa/m\n", param_dim.cohesive.K_n)              # TODO Chunk 2 Task 2.1
    # @printf("  δ_c_n = %.4e m\n", param_dim.cohesive.δ_c_n)               # TODO Chunk 2 Task 2.1
    @printf("  l_c   = %.2e m = %.1f μm\n", l_c, l_c*1e6)

    # ── nθ 区间 ──
    nθ_outer = ceil(Int, 2π * R_out / l_c)
    nθ_inner = ceil(Int, 2π * R_in  / l_c)

    println("\n[3] nθ 约束 (单元弧长 ≤ l_c)")
    @printf("  R_in  = %.4f mm,  R_out = %.4f mm\n", R_in*1e3, R_out*1e3)
    @printf("  外圈约束: nθ ≥ %d\n", nθ_outer)
    @printf("  内圈约束: nθ ≥ %d\n", nθ_inner)

    # ── 4 等分 ──
    step = max(1, round(Int, (nθ_outer - nθ_inner) / 3))
    nθ_czm = [nθ_inner + i * step for i in 0:3]
    nθ_czm[4] = nθ_outer  # 确保覆盖上界

    println("\n[4] CZM nθ 候选值 (每周)")
    for nθ in nθ_czm
        arc_out = R_out * 2π / nθ
        arc_in  = R_in  * 2π / nθ
        @printf("  nθ = %4d:  外弧 %.4f mm,  内弧 %.4f mm\n", nθ, arc_out*1e3, arc_in*1e3)
    end

    println("\n[5] 热网格 nθ 候选值 (经验)")
    nθ_thermal = [20, 40, 80, 160]
    for nθ in nθ_thermal
        arc_out = R_out * 2π / nθ
        @printf("  nθ = %4d:  外弧 %.4f mm\n", nθ, arc_out*1e3)
    end

    println("\n" * "=" ^ 70)
    println("完成")
    println("=" ^ 70)

    return (E_eff=E_eff, l_c=l_c, nθ_czm=nθ_czm, nθ_thermal=nθ_thermal)
end

main()
