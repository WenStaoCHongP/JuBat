include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
using Test
using LinearAlgebra
using Printf
using Plots
include(joinpath(@__DIR__, "unit_czm_newton.jl"))

"""
绘制条带变形云图：每个 Q4 按单元平均 |u| 填色，几何为放大后的变形构型；
灰虚线为未变形轮廓。位移相对几何通常极小，自动放大使最大位移 ≈ 0.15 H。
"""
function plot_unit_strip_deformed(czm_mesh, meta, u; out_path, title="")
    ne = size(czm_mesh.bulk_element, 1)
    nnode = czm_mesh.nnode
    Ux = [u[2n - 1] for n in 1:nnode]
    Uy = [u[2n] for n in 1:nnode]
    Umag = [hypot(Ux[n], Uy[n]) for n in 1:nnode]
    umax = maximum(Umag)
    H = meta.y_interfaces[end] - meta.y_interfaces[1]
    def_scale = umax > 0 ? (0.15 * H / umax) : 1.0

    field_elem = zeros(ne)
    for e in 1:ne
        ns = czm_mesh.bulk_element[e, :]
        field_elem[e] = sum(Umag[n] for n in ns) / 4
    end
    fmin, fmax = minimum(field_elem), maximum(field_elem)
    if fmax - fmin < 1e-30
        fmax = fmin + 1.0
    end
    cmap_grad = cgrad(:viridis)

    plt = plot(aspect_ratio=:equal, legend=:outertopright, size=(560, 720),
               xlabel="x*", ylabel="y*", title=title,
               colorbar_title="|u|*")

    for e in 1:ne
        ns = czm_mesh.bulk_element[e, :]
        xs = [czm_mesh.node[n, 1] for n in ns]; push!(xs, xs[1])
        ys = [czm_mesh.node[n, 2] for n in ns]; push!(ys, ys[1])
        plot!(plt, xs, ys; linestyle=:dash, linecolor=:gray60, linewidth=0.8,
              label=(e == 1 ? "undeformed" : false))
    end

    for e in 1:ne
        ns = czm_mesh.bulk_element[e, :]
        xs = Float64[]; ys = Float64[]
        for n in ns
            push!(xs, czm_mesh.node[n, 1] + def_scale * Ux[n])
            push!(ys, czm_mesh.node[n, 2] + def_scale * Uy[n])
        end
        push!(xs, xs[1]); push!(ys, ys[1])
        idx = (field_elem[e] - fmin) / (fmax - fmin)
        col = cmap_grad[clamp(idx, 0.0, 1.0)]
        plot!(plt, Plots.Shape(xs, ys); fillcolor=col, linecolor=:black,
              linewidth=0.6, label=false, colorbar_entry=false)
    end

    scatter!(plt, [NaN], [NaN]; zcolor=[fmin, fmax], clims=(fmin, fmax),
             colormap=:viridis, markersize=0, label=false,
             colorbar=true, colorbar_title="|u|*")

    for (k, coh) in enumerate(czm_mesh.cohesive_elements)
        n1, n2 = coh.nodes_bottom
        n4, n3 = coh.nodes_top
        xb = [czm_mesh.node[n1, 1] + def_scale * Ux[n1],
              czm_mesh.node[n2, 1] + def_scale * Ux[n2]]
        yb = [czm_mesh.node[n1, 2] + def_scale * Uy[n1],
              czm_mesh.node[n2, 2] + def_scale * Uy[n2]]
        xt = [czm_mesh.node[n4, 1] + def_scale * Ux[n4],
              czm_mesh.node[n3, 1] + def_scale * Ux[n3]]
        yt = [czm_mesh.node[n4, 2] + def_scale * Uy[n4],
              czm_mesh.node[n3, 2] + def_scale * Uy[n3]]
        plot!(plt, xb, yb; linecolor=:white, linewidth=1.2,
              label=(k == 1 ? "coh bottom" : false))
        plot!(plt, xt, yt; linecolor=:red, linewidth=1.2,
              label=(k == 1 ? "coh top" : false))
    end

    annotate!(plt, meta.width * 0.52, meta.y_interfaces[1] + 0.02 * H,
              text(@sprintf("magnification = %.3g×", def_scale), 8, :left, :white))

    mkpath(dirname(out_path))
    savefig(plt, out_path)
    return out_path, def_scale, umax
end

"""
外侧 PE（层1）顶贴固定 PCC，且全体 ε_x=0、底端 σ_y=0（平面应力）：
u_bot ≈ -ε0 (1+ν) h。ε0<0 收缩 → 底边上移为正。
"""
function analytic_free_coat_on_fixed_collector(h::Float64, ε0::Float64, ν::Float64)
    return -ε0 * (1.0 + ν) * h
end

@testset "eigenstrain (PCC/NCC + top fixed, all ux=0)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    param = case.param
    czm_mesh, meta = JuBat.create_unit_czm_strip(param; y0=1.0)
    cache = JuBat.compute_czm_params_per_interface(case)
    # Δsoc=0.8 + E_coat=500MPa：体太软，原 σ_max 下 δ/δ0≪1。
    # 调整：降低 σ_max（更易起裂），并设定 δc/δ0（控制软化陡峭度，利于收敛）。
    # 物理对应：NE 界面 σ_max ≈ 0.1×92 MPa ≈ 9 MPa；δc 不变则 G_c 同比下降。
    function retune_czm_for_damage!(cache; σ_scale::Float64, δc_over_δ0::Float64)
        for key in (:PE_PCC, :NE_NCC)
            ip = cache.by_interface[key]
            σ = ip.σ_max * σ_scale
            τ = ip.τ_max * σ_scale
            δc = ip.δ_c_n
            δ0 = δc / δc_over_δ0
            Kn = σ / δ0
            δc_t = ip.δ_c_t
            δ0_t = δc_t / δc_over_δ0
            Kt = τ / δ0_t
            cache.by_interface[key] = JuBat.CzmInterfaceParams(
                E_eff=ip.E_eff, ν=ip.ν, α=ip.α, Λ=ip.Λ,
                E_star=ip.E_star, L_ch=ip.L_ch,
                σ_max=σ, K_n=Kn, δ_0_n=δ0, δ_c_n=δc, G_c=0.5 * σ * δc,
                τ_max=τ, K_t=Kt, δ_0_t=δ0_t, δ_c_t=δc_t, G_c_t=0.5 * τ * δc_t,
                η=ip.η, czm_model=ip.czm_model,
                h_c0=ip.h_c0, k_air=ip.k_air, lambda_m=ip.lambda_m,
                beta=ip.beta, threshold=ip.threshold,
            )
        end
    end
    σ_scale = 0.10
    δc_over_δ0 = 50.0
    retune_czm_for_damage!(cache; σ_scale=σ_scale, δc_over_δ0=δc_over_δ0)

    czm_mesh.damage_states = [JuBat.DamageState() for _ in 1:4]

    α_pe = param.PE.alphaT
    α_ne = param.NE.alphaT
    β_n = param.NE.Omega / 3.0
    β_p = param.PE.Omega / 3.0
    _, ν_pe = JuBat.moduli_of(param, :PE)

    T_ref = case.param_dim.scale.T_ref
    ΔT_phys = 10.0
    # 放电方向（Δsoc<0）：PE 膨胀（β_p<0）、NE 收缩（β_n>0），NE_NCC 界面进入受拉
    # （Δsoc=+0.8 时全部界面受压、D≡0；线性弹性下反号给出 coh4 δ_n≈+7.7e-3 > δ0_NE）
    Δsoc_frac = -0.8
    ΔT_end = ΔT_phys / T_ref
    n_steps = 40
    n_hold = 10
    u = zeros(2 * czm_mesh.nnode)
    δ0_pe = cache.by_interface[:PE_PCC].δ_0_n
    δ0_ne = cache.by_interface[:NE_NCC].δ_0_n
    δ_abs_max_over_steps = 0.0
    D_max_last = 0.0
    u_pe1_bot_last = 0.0
    u_pe1_ana_last = 0.0
    seps_last = Vector{Tuple{Float64,Float64}}(undef, 4)
    frac_hold = nothing

    @printf("[eigenstrain] BC: fix PCC+NCC+top (uy=0); all nodes ux=0\n")
    @printf("[eigenstrain] E_coat PE/NE = 500 MPa; Δsoc=%.2f\n", Δsoc_frac)
    @printf("[eigenstrain] CZM retune: σ_max×%.2f, δc/δ0=%.0f  (σ*_NE=%.3e, δ0_NE=%.3e)\n",
            σ_scale, δc_over_δ0, cache.by_interface[:NE_NCC].σ_max, δ0_ne)
    @printf("[eigenstrain] load: ΔT=%.1f K (ΔT*=%.6f)\n", ΔT_phys, ΔT_end)

    function bc_fix_collectors_top_all_ux(meta, nnode)
        bc_dofs = Int[]; bc_vals = Float64[]
        for n in 1:nnode
            push!(bc_dofs, 2n - 1); push!(bc_vals, 0.0)
        end
        for n in Iterators.flatten((meta.pcc_nodes, meta.ncc_nodes,
                                    meta.top_nodes_after_czm))
            push!(bc_dofs, 2n); push!(bc_vals, 0.0)
        end
        return bc_dofs, bc_vals
    end

    total_steps = n_steps + n_hold
    for s in 1:total_steps
        frac = frac_hold === nothing ? min(s / n_steps, 1.0) : frac_hold
        ΔT = ΔT_end * frac
        Δsoc = Δsoc_frac * frac
        dT = fill(ΔT, 8)
        Δsoc_n = zeros(8)
        Δsoc_p = zeros(8)
        for e in 1:8
            mt = meta.layer_materials[e]
            if mt === :PE
                Δsoc_p[e] = Δsoc
            elseif mt === :NE
                Δsoc_n[e] = Δsoc
            end
        end

        F_tc = JuBat.assemble_thermal_chemical_load(
            czm_mesh, cache, dT, Δsoc_n, Δsoc_p)
        bc_dofs, bc_vals = bc_fix_collectors_top_all_ux(meta, czm_mesh.nnode)
        u, seps, tracts, ok, Rn = unit_czm_newton_step!(
            czm_mesh, u, cache; bc_dofs=bc_dofs, bc_vals=bc_vals,
            F_thermo_chem=F_tc, max_iter=200, tol=1e-7, visc_beta=0.3)
        if !ok
            @printf("[eigenstrain] FAIL step %d/%d frac=%.3f Rn=%.3e Dmax=%.3e\n",
                    s, total_steps, frac, Rn,
                    maximum(st.D for st in czm_mesh.damage_states))
        end
        @test ok

        ε0 = [JuBat.eigenstrain_of(param, meta.layer_materials[e], dT[e], Δsoc_n[e], Δsoc_p[e])
              for e in 1:8]
        u_pe1_ana = analytic_free_coat_on_fixed_collector(
            meta.heights[1], ε0[1], ν_pe)
        u_pe1_bot = sum(u[2n] for n in meta.bottom_nodes) / length(meta.bottom_nodes)

        for i in 1:4
            δ_abs_max_over_steps = max(δ_abs_max_over_steps, abs(seps[i][1]))
        end
        D_max_last = maximum(st.D for st in czm_mesh.damage_states)
        u_pe1_bot_last = u_pe1_bot
        u_pe1_ana_last = u_pe1_ana
        seps_last = seps

        # 一旦进入损伤：冻结载荷，避免深软化不收敛
        if frac_hold === nothing && D_max_last > 0.02
            frac_hold = frac
            @printf("[eigenstrain] damage onset at step %d  frac=%.3f  Dmax=%.4e → hold load\n",
                    s, frac, D_max_last)
        end

        if s == total_steps
            u_pcc = maximum(hypot(u[2n - 1], u[2n]) for n in meta.pcc_nodes)
            u_ncc = maximum(hypot(u[2n - 1], u[2n]) for n in meta.ncc_nodes)
            u_top = maximum(abs(u[2n]) for n in meta.top_nodes_after_czm)
            ux_max = maximum(abs(u[2n - 1]) for n in 1:czm_mesh.nnode)
            pe3_bot = czm_mesh.cohesive_elements[2].nodes_top
            ne5_top = czm_mesh.cohesive_elements[3].nodes_bottom
            u_mid_lo = sum(u[2n] for n in pe3_bot) / length(pe3_bot)
            u_mid_hi = sum(u[2n] for n in ne5_top) / length(ne5_top)
            @printf("[eigenstrain] end: ΔT_phys=%.1f K  Δsoc=%.2f (frac=%.3f)\n",
                    ΔT_phys, Δsoc, frac)
            εT_pe = α_pe * ΔT
            εT_ne = α_ne * ΔT
            εchem_p = β_p * Δsoc
            εchem_n = β_n * Δsoc
            @printf("  strain PE: ε_T*=%.6e  ε_chem*=%.6e  ε0*=%.6e\n",
                    εT_pe, εchem_p, εT_pe + εchem_p)
            @printf("  strain NE: ε_T*=%.6e  ε_chem*=%.6e  ε0*=%.6e\n",
                    εT_ne, εchem_n, εT_ne + εchem_n)
            @printf("  coeffs: α_pe*=%.6e  α_ne*=%.6e  β_p*=%.6e  β_n*=%.6e\n",
                    α_pe, α_ne, β_p, β_n)
            @printf("  max|u|_PCC*=%.3e  max|u|_NCC*=%.3e  max|uy|_top*=%.3e  max|ux|*=%.3e\n",
                    u_pcc, u_ncc, u_top, ux_max)
            @printf("  PE1 bottom uy* FEM=%.6e  ana(-ε0(1+ν)h)=%.6e  (ratio=%.3f)\n",
                    u_pe1_bot, u_pe1_ana, u_pe1_bot / u_pe1_ana)
            @printf("  mid faces uy* PE3bot=%.3e  NE5top=%.3e\n", u_mid_lo, u_mid_hi)
            for i in 1:4
                iface = czm_mesh.cohesive_elements[i].interface_type
                δ0_i = iface === :PE_PCC ? δ0_pe : δ0_ne
                Tn, _ = tracts[i]
                @printf("  coh%d %-7s  δ_n=%.6e  δ/δ0=%.4e  T_n=%.6e  D=%.6e\n",
                        i, iface, seps[i][1], seps[i][1] / δ0_i, Tn,
                        czm_mesh.damage_states[i].D)
            end
            @printf("[eigenstrain] D: max=%.6e\n", D_max_last)
            @test u_pcc < 1e-10
            @test u_ncc < 1e-10
            @test u_top < 1e-10
            @test ux_max < 1e-10
            @test abs(seps[1][1]) / δ0_pe < 1e-5
        end
    end

    @printf("[eigenstrain] max_|δ_n|=%.6e  δ0_PE=%.6e  δ0_NE=%.6e  max_|δ|/δ0_PE=%.4e\n",
            δ_abs_max_over_steps, δ0_pe, δ0_ne, δ_abs_max_over_steps / δ0_pe)
    @test abs(u_pe1_ana_last) > 1e-10
    @test u_pe1_bot_last ≈ u_pe1_ana_last rtol=0.35 atol=1e-10
    @test D_max_last > 0.0
    @test any(czm_mesh.damage_states[i].D > 0 && seps_last[i][1] > 0 for i in 1:4)

    out_dir = joinpath(@__DIR__, "../output", "unit_czm_eigenstrain")
    out_png = joinpath(out_dir, "unit_czm_eigenstrain_deformed.png")
    path, mag, umax = plot_unit_strip_deformed(
        czm_mesh, meta, u;
        out_path=out_png,
        title=@sprintf("damage; ΔT=%.0fK, Δsoc=%.2f, Dmax=%.3f",
                       ΔT_phys, Δsoc_frac * (frac_hold === nothing ? 1.0 : frac_hold),
                       D_max_last))
    @printf("[eigenstrain] deformed plot → %s  (umax*=%.3e, mag=%.3g×)\n", path, umax, mag)
    @test isfile(path)
end

println("unit_czm_eigenstrain: ALL PASS")
