using LinearAlgebra, Statistics, Plots
include("../src/JuBat.jl")

# tools/check_thermal_kernels.jl
# Compute and export per-element effective radial/tangential thermal conductivities
# Usage: julia --project=. tools/check_thermal_kernels.jl

function main()
    println("Loading param and mesh...")
    param_dim = JuBat.ChooseCell("Jellyroll")
    # collector-seeded mesh gives consistent layer weights per element
    mesh = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=360, gsorder=2)

    centers = JuBat.jellyroll_element_centers(mesh)
    ne = size(mesh.element,1)

    # get element layer weights if available (collector-seeded cached) or compute by sampling
    fks = JuBat.jellyroll_get_layer_weights(mesh)
    if fks === nothing
        println("No cached layer weights found; sampling element areas to compute layer weights (may be slow)")
        fks = JuBat.jellyroll_element_layer_weights(mesh, param_dim; nsamples_per_dim=6, logic=:spiral)
    end

    # material lambdas (W/mK)
    λ_NE  = getfield(param_dim.NE,  :lambda)
    λ_SP  = getfield(param_dim.SP,  :lambda)
    λ_PE  = getfield(param_dim.PE,  :lambda)
    λ_PCC = getfield(param_dim.PCC,:lambda)
    λ_NCC = getfield(param_dim.NCC,:lambda)

    lambdas = [λ_NE, λ_SP, λ_PE, λ_PCC, λ_NCC]

    lam_r = zeros(Float64, ne)
    lam_t = zeros(Float64, ne)

    for e in 1:ne
        fk = fks[e, :]
        # harmonic-style radial effective: 1 / sum(f_k / λ_k)
        denom = 0.0
        for i in 1:5
            λi = max(lambdas[i], 1e-16)
            denom += fk[i] / λi
        end
        lam_r[e] = denom > 0 ? (1.0 / denom) : 0.0
        # tangential (parallel) arithmetic mean
        lam_t[e] = sum(fk .* lambdas)
    end

    # global spiral effective for comparison
    pgeo = JuBat.jellyroll_spiral_params(param_dim)
    λr_eff_geo = pgeo.λ_r_eff
    λt_eff_geo = pgeo.λ_t_eff

    println("Summary (W/mK):")
    println("  lam_r per-element: mean=$(mean(lam_r)), min=$(minimum(lam_r)), max=$(maximum(lam_r))")
    println("  lam_t per-element: mean=$(mean(lam_t)), min=$(minimum(lam_t)), max=$(maximum(lam_t))")
    println("  spiral params λ_r_eff=$(λr_eff_geo), λ_t_eff=$(λt_eff_geo)")

    # write CSV
    outcsv = "tools/thermal_kernels_per_element.csv"
    open(outcsv, "w") do io
        println(io, "elem,x,y,r,lam_r,lam_t,f_NE,f_SP,f_PE,f_PCC,f_NCC")
        for e in 1:ne
            x = centers[e,1]; y = centers[e,2]
            r = hypot(x,y)
            fk = fks[e, :]
            println(io, join([e, x, y, r, lam_r[e], lam_t[e], fk[1], fk[2], fk[3], fk[4], fk[5]], ","))
        end
    end
    println("Wrote: $outcsv")

    # scatter plots of lam_r and lam_t
    xs = centers[:,1]; ys = centers[:,2]
    plt1 = scatter(xs, ys, zcolor=lam_r, markerstrokewidth=0, ms=6, color=cgrad([:skyblue, :red]), title="Element λ_r (W/mK)", aspect_ratio=1)
    savefig(plt1, "tools/thermal_lambda_r.png")
    plt2 = scatter(xs, ys, zcolor=lam_t, markerstrokewidth=0, ms=6, color=cgrad([:skyblue, :red]), title="Element λ_t (W/mK)", aspect_ratio=1)
    savefig(plt2, "tools/thermal_lambda_t.png")
    println("Saved: tools/thermal_lambda_r.png, tools/thermal_lambda_t.png")

    return nothing
end

main()
