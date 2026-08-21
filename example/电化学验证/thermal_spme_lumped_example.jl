using Plots, CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "JuBat.jl"))

const OUTPUT_DIR = joinpath(PROJECT_ROOT, "output", "thermal_spme_lumped_example")
mkpath(OUTPUT_DIR)

function build_case(cell_name::String)
    param_dim = JuBat.ChooseCell(cell_name)
    param_dim.cell.v_l = 2.5

    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermalmodel = "lumped"
    opt.dtType = "auto"
    opt.jacobi = "update"
    opt.time = [0.0, 3600.0]
    opt.dt = [1.0, 10.0]

    i_app = param_dim.cell.I1C
    opt.Current = _ -> i_app

    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    capacity_ah = result["time [s]"] ./ 3600.0 .* abs(i_app)

    return capacity_ah, result
end

function main()
    ref_path = joinpath(PROJECT_ROOT, "src", "data", "pybamm_SPMe_LGM50_1.0C.csv")
    ref_tbl = CSV.File(ref_path)
    ref_capacity = Float64.(getproperty(ref_tbl, Symbol("capacity [A.h]")))
    ref_voltage = Float64.(getproperty(ref_tbl, Symbol("voltage [V]")))
    ref_temperature = Float64.(getproperty(ref_tbl, Symbol("temperature [K]")))

    cell_specs = [
        (name="LG M50", label="LGM50", color=:black),
        (name="Jellyroll", label="Jellyroll", color=:blue),
    ]

    pV = plot(xlabel="Output capacity [Ah]", ylabel="Cell voltage [V]", legend=:bottomleft)
    ylims!(pV, 2.5, 4.3)
    pT = plot(xlabel="Output capacity [Ah]", ylabel="Temperature [K]", legend=:bottomright)

    for spec in cell_specs
        capacity_ah, result = build_case(spec.name)
        plot!(pV, capacity_ah, result["cell voltage [V]"], label="$(spec.label) (JuBat)", linecolor=spec.color, linewidth=2)
        plot!(pT, capacity_ah, result["temperature [K]"], label="$(spec.label) (JuBat)", linecolor=spec.color, linewidth=2)
    end

    plot!(pV, ref_capacity, ref_voltage, label="PyBaMM SPMe LGM50 1C", linestyle=:dash, linecolor=:red, linewidth=2)
    plot!(pT, ref_capacity, ref_temperature, label="PyBaMM SPMe LGM50 1C", linestyle=:dash, linecolor=:red, linewidth=2)

    savefig(pV, joinpath(OUTPUT_DIR, "thermal_spme_lumped_example-V.pdf"))
    savefig(pT, joinpath(OUTPUT_DIR, "thermal_spme_lumped_example-T.pdf"))
end

main()