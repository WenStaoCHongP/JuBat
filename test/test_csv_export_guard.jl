using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))

@testset "CSV export guarded write" begin
    files_written = String[]
    files_skipped = String[]
    calls = Symbol[]

    result = JuBat.write_csv_guarded!("success.csv", files_written, files_skipped) do
        push!(calls, :success)
    end
    @test result === nothing
    @test calls == [:success]
    @test files_written == ["success.csv"]
    @test isempty(files_skipped)

    @test_logs (:warn, r"Failed to write failure.csv") JuBat.write_csv_guarded!(
        "failure.csv", files_written, files_skipped) do
        error("expected test failure")
    end
    @test files_written == ["success.csv"]
    @test files_skipped == ["failure.csv"]
end

@testset "cohesive damage CSV uses cumulative winding angle" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(
        case.param; nθ=8, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    czm_mesh = JuBat.create_czm_mesh(
        mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    n_coh = czm_mesh.n_cohesive
    snapshot = JuBat.CZMSnapshot(
        0.0, 1, "charge", zeros(2 * czm_mesh.nnode),
        zeros(n_coh), zeros(n_coh), zeros(n_coh), zeros(n_coh), zeros(n_coh),
        true, 1, 0.0, "basic")
    result = (czm_snapshots=[snapshot],)

    mktempdir() do output_dir
        JuBat.write_cohesive_damage_csv(
            result, case, czm_mesh, output_dir, true, JuBat.CsvExportOptions())
        lines = readlines(joinpath(output_dir, "cohesive_damage.csv"))
        @test endswith(first(lines), ",theta_cum_deg")

        actual = [parse(Float64, split(line, ',')[end]) for line in lines[2:end]]
        a = case.param.cell.Rin
        b = case.param.cell.layer / (2pi)
        expected = [
            begin
                n1, n2 = elem.nodes_bottom
                mx = 0.5 * (czm_mesh.node[n1, 1] + czm_mesh.node[n2, 1])
                my = 0.5 * (czm_mesh.node[n1, 2] + czm_mesh.node[n2, 2])
                (hypot(mx, my) - a) / b * 180 / pi
            end
            for elem in czm_mesh.cohesive_elements
        ]
        @test actual ≈ expected rtol=1e-12 atol=1e-10
        @test maximum(actual) > 360.0
    end
end

@testset "CSV export minimal public path" begin
    result = (cycle_results=Any[], soh=Float64[], czm_snapshots=Any[])
    sim_case = (opt=(thermal_enabled=false,), geometry=nothing)

    mktempdir() do output_dir
        files_written = JuBat.export_cycling_csv(
            result, sim_case, nothing; output_dir=output_dir, overwrite=true)
        @test files_written == ["cycle_summary.csv"]
        @test isfile(joinpath(output_dir, "cycle_summary.csv"))
        @test startswith(read(joinpath(output_dir, "cycle_summary.csv"), String), "cycle,phase")
    end
end
