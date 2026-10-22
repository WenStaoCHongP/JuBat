using Test, CSV, DelimitedFiles, Serialization

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "output", "test_cycling_csv_timeline")
mkpath(OUT)
mode = ARGS[1]
tag = length(ARGS) > 1 ? ARGS[2] : "candidate"
include(joinpath(ROOT, "example", mode * ".jl"))

# One mechanical snapshot per phase: equal local timestamps cannot identify phases.
times = (discharge=0.5, rest1=0.5, charge=0.5, rest2=0.5)
if length(ARGS) > 2 && ARGS[3] == "unequal"
    times = (discharge=1.0, rest1=0.5, charge=1.5, rest2=1.0)
end
ENV["JUBAT_SNAPSHOT_CYCLES"] = "1"
if "--reuse" in ARGS
    saved = deserialize(joinpath(OUT, mode * "_" * tag * "_source.bin"))
    result, case, dir = saved.result, saved.case, saved.dir
elseif mode == "czm_cycling"
    ENV["JUBAT_CYCLES"] = "1"
    ENV["JUBAT_EXPORT_FIELD_DATA"] = "1"
    ENV["JUBAT_RUN_TAG"] = "csv_timeline_" * tag
    run = main(; n_theta=20, phase_times=times, dt_cycle=[0.5, 0.5])
    result, case, dir = run.result, run.case, run.outdir
elseif mode == "czm_cycling_n20"
    p = JuBat.ChooseCell("Jellyroll")
    tc = TAB_CASES[1]
    thetas = theta_pos_for(tc.theta_kind, spiral_theta_range(p))
    run = run_tab_case(tc, thetas, 1000.0, joinpath(OUT, mode, tag), nothing;
        n_theta=20, n_cycles=1, strength_scale=0.085, stiffness_scale=0.1,
        gc_scale=2.0, iter_method="gp_basic", tau_seconds="5",
        phase_times=times, dt_cycle=[0.5, 0.5])
    @test run.summary.status == "OK"
    result, case, dir = run.result, run.case, run.casedir
else
    error("unknown example: $mode")
end

if !("--reuse" in ARGS)
    # Persist numeric source data so checker-only changes never require another solve.
    cm=case.mesh["czm"]; th=case.mesh["thermal2D"]; scale=case.param.scale
    source_case=(param=(scale=(L=scale.L,δ_czm=scale.δ_czm,σ_czm=scale.σ_czm),),
        mesh=Dict("czm"=>(node=cm.node,bulk_element=cm.bulk_element,cohesive_elements=cm.cohesive_elements),
                  "thermal2D"=>(nlen=th.nlen,node=th.node,element=th.element)))
    serialize(joinpath(OUT, mode * "_" * tag * "_source.bin"), (;result,case=source_case,dir))
end

rows(name) = collect(CSV.File(joinpath(dir, name); stringtype=String))
near(a, b) = isapprox(a, b; rtol=2e-14, atol=1e-14)
phases = ("discharge", "rest1", "charge", "rest2")
cr = only(result.cycle_results)
phase_results = [getproperty(cr, Symbol(p)) for p in phases]
snaps = result.czm_snapshots
# Independent expectation: every mechanical step matches the real accepted
# output column (update_interval=1), including a tiny final solver step.
expected_phase = String[]
expected_time = Float64[]
for (label, ph) in zip(phases, phase_results)
    @test ph.duration == getproperty(times, Symbol(label))
    for t in vec(ph.solve_result["time [s]"])[2:end]
        push!(expected_phase, label)
        push!(expected_time, ph.t_start + t)
    end
end

@testset failfast=true "Real four-phase CSV timeline and values: $mode / $tag" begin
    @test length(snaps) == length(expected_time)
    @test [s.phase for s in snaps] == expected_phase
    @test all(diff(expected_time) .>= 0)
    ps = rows("phase_summary.csv")
    @test [p.phase for p in ps] == collect(phases)
    @test all(near(p.t_start_s, q.t_start) && near(p.t_end_s, q.t_end) for (p,q) in zip(ps,phase_results))
    @test [p.t_start_s for p in ps][2:end] == [p.t_end_s for p in ps][1:end-1]
    cs = only(rows("cycle_summary.csv"))
    @test near(cs.capacity_discharge_Ah, only(result.capacity_discharge))
    @test near(cs.D_max, only(result.D_max))
    @test near(cs.T_max_K, only(result.T_max))
    for (filename, values) in (
        ("damage_history_cyc1.csv", s->s.damage),
        ("sep_n_history_cyc1.csv", s->s.separation_n .* case.param.scale.δ_czm),
        ("sep_t_history_cyc1.csv", s->s.separation_t .* case.param.scale.δ_czm),
        ("traction_n_history_cyc1.csv", s->s.traction_n .* case.param.scale.σ_czm),
        ("traction_t_history_cyc1.csv", s->s.traction_t .* case.param.scale.σ_czm))
        data = readdlm(joinpath(dir,filename), ',', Float64; skipstart=1)
        @test size(data) == (length(snaps), length(snaps[1].damage)+1)
        @test all(near.(data[:,1], expected_time))
        @test all(isfinite,data)
        @test all(near.(data[:,2:end], permutedims(reduce(hcat, values.(snaps)))))
    end
    @test length(rows("element_map.csv")) == length(case.mesh["czm"].cohesive_elements)
    if mode == "czm_cycling"
        for (filename,key,id,value) in (("node_temperature.csv","thermal2D temperature at nodes [K]",:node_id,:T_K),
                                        ("element_currents.csv","thermal2D element current",:elem_id,:I_e))
            actual=rows(filename); offset=0
            for (label,ph) in zip(phases,phase_results)
                raw=ph.solve_result[key]; ts=vec(ph.solve_result["time [s]"])
                for k in axes(raw,2), j in axes(raw,1)
                    offset+=1; row=actual[offset]
                    @test row.cycle==1 && row.phase==label && getproperty(row,id)==j
                    @test near(row.time_s,ph.t_start+ts[k])
                    @test near(getproperty(row,value),raw[j,k])
                end
            end
            @test offset==length(actual)
        end
        for (filename,n,check) in (
            ("node_displacement.csv",size(case.mesh["czm"].node,1),(r,s,j)->near(r.ux,s.displacement[2j-1]*case.param.scale.L) && near(r.uy,s.displacement[2j]*case.param.scale.L)),
            ("cohesive_damage.csv",length(snaps[1].damage),(r,s,j)->near(r.D,s.damage[j]) && near(r.sep_n_m,s.separation_n[j]*case.param.scale.δ_czm) && near(r.sep_t_m,s.separation_t[j]*case.param.scale.δ_czm) && near(r.traction_n_Pa,s.traction_n[j]*case.param.scale.σ_czm) && near(r.traction_t_Pa,s.traction_t[j]*case.param.scale.σ_czm)))
            actual=rows(filename)
            @test length(actual)==length(snaps)*n
            for (i,s) in enumerate(snaps), j in 1:n
                r=actual[(i-1)*n+j]
                @test r.phase==expected_phase[i] && near(r.time_s,expected_time[i])
                @test check(r,s,j)
            end
        end
        actual=rows("element_strain.csv"); cm=case.mesh["czm"]; n=size(cm.bulk_element,1)
        @test length(actual)==length(snaps)*n
        for (i,s) in enumerate(snaps)
            xx,yy,xy=JuBat.recover_bulk_strain(cm.node,cm.bulk_element,s.displacement)
            for j in 1:n
                r=actual[(i-1)*n+j]
                @test r.phase==expected_phase[i] && near(r.time_s,expected_time[i])
                @test near(r.eps_xx,xx[j]) && near(r.eps_yy,yy[j]) && near(r.gam_xy,xy[j])
            end
        end
        stress=rows("element_stress.csv"); offset=0
        for (label,ph) in zip(phases,phase_results)
            sr=ph.solve_result
            for j in axes(sr["diffusion stress xx [Pa]"],1)
                offset+=1; r=stress[offset]
                @test r.phase==label && near(r.time_s,ph.t_end)
                @test all(near(getproperty(r,col),sr[key][j,end]) for (col,key) in
                    ((:sxx_Pa,"diffusion stress xx [Pa]"),(:syy_Pa,"diffusion stress yy [Pa]"),(:sxy_Pa,"diffusion stress xy [Pa]"),(:vm_Pa,"diffusion stress vonMises [Pa]")))
            end
        end
        @test offset==length(stress)
        @test length(filter(f->endswith(f,".csv"),readdir(dir)))==14
    else
        manifest=Dict(r.key=>string(r.value) for r in rows("case_manifest.csv"))
        maxima=[maximum(s.damage) for s in snaps]
        sel=findfirst(>=(0.999*maximum(maxima)),maxima)
        @test near(parse(Float64,manifest["t_sel_s"]),expected_time[sel])
        @test manifest["phase_sel"]==expected_phase[sel]
        @test all(near(r.D,v) for (r,v) in zip(rows("cloud_damage.csv"),snaps[sel].damage))
        ph=phase_results[findfirst(==(expected_phase[sel]),collect(phases))]
        ts=ph.t_start .+ vec(ph.solve_result["time [s]"])
        k=argmin(abs.(ts.-expected_time[sel]))
        @test all(near(r.T_K,v) for (r,v) in zip(rows("cloud_temperature.csv"),ph.solve_result["thermal2D temperature at nodes [K]"][:,k]))
        @test all(near(r.I_e,v) for (r,v) in zip(rows("cloud_current.csv"),ph.solve_result["thermal2D element current"][:,k]))
        @test length(rows("mesh_nodes.csv"))==case.mesh["thermal2D"].nlen
        @test length(rows("mesh_elements.csv"))==size(case.mesh["thermal2D"].element,1)
        @test !isempty(rows("tab_nodes.csv"))
        @test length(filter(f->endswith(f,".csv"),readdir(dir)))==15
    end
end
println("CSV_TIMELINE_VALIDATED ", dir)
